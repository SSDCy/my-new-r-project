# utils/esm_search.py
"""
ESM2 蛋白质功能嵌入搜索模块
- 支持加载预构建参考库
- 支持在线构建自定义参考库（带缓存）
- 支持根据蛋白ID自动获取序列
- 提供单条和批量相似蛋白检索
- 提供 UniProt 注释批量获取（带 ID 格式检查）
- 缓存损坏自动修复
- 新增：分步批量检索，避免长时间阻塞 UI
"""

import streamlit as st
import numpy as np
import pandas as pd
import torch
import os
import re
import time
import pickle
import hashlib
import requests
from typing import Dict, List, Tuple, Optional

DEBUG = True
def debug(msg: str):
    if DEBUG:
        print(f"[ESM2 DEBUG] {msg}")

# ---------- 依赖检查 ----------
try:
    import esm
    ESM_AVAILABLE = True
    debug("fair-esm 已安装。")
except ImportError:
    ESM_AVAILABLE = False
    debug("fair-esm 未安装，ESM2 功能将不可用。")

def show_install_instructions():
    st.error("**ESM2 功能不可用**，请安装 PyTorch 和 fair-esm：")
    st.code("pip install torch fair-esm", language="bash")
    st.info("安装完成后请重启 Streamlit 应用。")

# ---------- 全局缓存 ----------
@st.cache_resource
def load_pretrained_ref_db():
    db_path = "data/ref_db.npz"
    ann_path = "data/ref_annotations.pkl"
    if os.path.exists(db_path) and os.path.exists(ann_path):
        debug("加载预构建参考库...")
        try:
            data = np.load(db_path, allow_pickle=True)
            ids = data["ids"].tolist()
            embeddings_norm = data["embeddings"]
            with open(ann_path, "rb") as f:
                annotations = pickle.load(f)
            debug(f"参考库加载完成：{len(ids)} 条蛋白")
            return ids, embeddings_norm, annotations
        except Exception as e:
            debug(f"预构建参考库加载失败: {e}")
            return None, None, None
    else:
        debug("预构建参考库文件不存在")
        return None, None, None

def is_pretrained_available():
    return os.path.exists("data/ref_db.npz") and os.path.exists("data/ref_annotations.pkl")

# ---------- ESM2 模型 ----------
@st.cache_resource
def load_esm_model(model_name: str = "esm2_t12_35M_UR50D"):
    if not ESM_AVAILABLE:
        show_install_instructions()
        return None, None, None
    debug(f"加载 ESM2 模型：{model_name} ...")
    try:
        model, alphabet = esm.pretrained.load_model_and_alphabet(model_name)
        model.eval()
        if torch.cuda.is_available():
            model = model.cuda()
            debug("使用 GPU 加速")
        else:
            debug("使用 CPU 运行")
        batch_converter = alphabet.get_batch_converter()
        debug("ESM2 模型加载成功")
        return model, alphabet, batch_converter
    except Exception as e:
        debug(f"模型加载失败: {str(e)}")
        st.error(f"模型加载失败: {str(e)}")
        return None, None, None

def get_esm_embedding(sequence: str, model, batch_converter) -> np.ndarray:
    if not ESM_AVAILABLE:
        raise RuntimeError("ESM2 不可用")
    if len(sequence) < 10:
        raise ValueError(f"序列太短（{len(sequence)} aa），至少需要 10 个氨基酸")
    data = [("query", sequence)]
    _, _, batch_tokens = batch_converter(data)
    if torch.cuda.is_available():
        batch_tokens = batch_tokens.cuda()
    with torch.no_grad():
        results = model(batch_tokens, repr_layers=[12])
    token_embeddings = results["representations"][12]
    embedding = token_embeddings.mean(dim=1).squeeze(0).cpu().numpy()
    return embedding

# ---------- 缓存管理 ----------
def compute_fasta_hash(fasta_text: str) -> str:
    return hashlib.sha256(fasta_text.encode('utf-8')).hexdigest()

def get_cache_path(hash_val: str) -> str:
    os.makedirs("data", exist_ok=True)
    return f"data/custom_cache_{hash_val}.npz"

def load_custom_cache(hash_val: str) -> Optional[Dict[str, np.ndarray]]:
    cache_file = get_cache_path(hash_val)
    if os.path.exists(cache_file):
        debug(f"发现缓存文件: {cache_file}")
        try:
            data = np.load(cache_file, allow_pickle=True)
            ids = data["ids"].tolist()
            embeddings = data["embeddings"]
            lib = {id_: embeddings[i] for i, id_ in enumerate(ids)}
            debug(f"缓存加载成功，包含 {len(lib)} 个蛋白")
            return lib
        except Exception as e:
            debug(f"缓存文件损坏，自动删除: {cache_file}")
            debug(f"错误详情: {e}")
            try:
                os.remove(cache_file)
            except Exception as rm_err:
                debug(f"无法删除损坏缓存: {rm_err}")
            return None
    return None

def save_custom_cache(hash_val: str, embeddings_dict: Dict[str, np.ndarray], ids: list):
    cache_file = get_cache_path(hash_val)
    temp_file = cache_file + ".tmp"
    try:
        os.makedirs(os.path.dirname(cache_file), exist_ok=True)
        embeddings = np.stack([embeddings_dict[id_] for id_ in ids], axis=0)
        np.savez(temp_file, ids=np.array(ids), embeddings=embeddings)
        test_data = np.load(temp_file, allow_pickle=True)
        if "ids" not in test_data or "embeddings" not in test_data:
            raise ValueError("缓存文件验证失败：缺少必要字段")
        if os.path.exists(cache_file):
            os.remove(cache_file)
        os.rename(temp_file, cache_file)
        debug(f"缓存保存成功: {cache_file}")
    except Exception as e:
        debug(f"缓存保存失败: {e}")
        if os.path.exists(temp_file):
            try:
                os.remove(temp_file)
            except:
                pass

# ---------- 构建自定义库 ----------
def build_reference_library(fasta_text: str, model, batch_converter,
                            use_cache: bool = True) -> Tuple[Dict[str, np.ndarray], List[str]]:
    debug("开始从 FASTA 构建参考库...")
    if use_cache:
        hash_val = compute_fasta_hash(fasta_text)
        debug(f"FASTA 哈希值: {hash_val}")
        cached = load_custom_cache(hash_val)
        if cached is not None:
            st.success("Cache loaded successfully, no need to rebuild.")  # 修改为英文
            return cached, list(cached.keys())
        else:
            debug("未找到有效缓存，将重新生成嵌入。")

    sequences, ids = [], []
    current_id, current_seq = None, []
    for line in fasta_text.splitlines():
        line = line.strip()
        if line.startswith(">"):
            if current_id:
                sequences.append("".join(current_seq))
                ids.append(current_id)
            current_id = line[1:].split()[0]
            current_seq = []
        else:
            if current_id is not None:
                current_seq.append(line)
    if current_id:
        sequences.append("".join(current_seq))
        ids.append(current_id)

    debug(f"解析到 {len(ids)} 条序列。")
    if len(ids) == 0:
        st.error("FASTA 文件中未找到任何序列")
        return {}, []

    embeddings = {}
    failed = []
    progress_bar = st.progress(0)
    status_text = st.empty()

    for i, (seq_id, seq) in enumerate(zip(ids, sequences)):
        progress = (i + 1) / len(ids)
        progress_bar.progress(progress)
        status_text.text(f"正在处理 {i+1}/{len(ids)}: {seq_id} (长度={len(seq)})")
        debug(f"处理 {i+1}/{len(ids)}: {seq_id} (长度={len(seq)})")
        try:
            embeddings[seq_id] = get_esm_embedding(seq, model, batch_converter)
        except Exception as e:
            debug(f"生成嵌入失败 ID={seq_id}: {str(e)}")
            failed.append(seq_id)
            continue

    progress_bar.empty()
    status_text.empty()

    if failed:
        st.warning(f"以下蛋白嵌入生成失败 ({len(failed)} 个): {', '.join(failed[:10])}{'...' if len(failed) > 10 else ''}")

    debug(f"自定义库构建完成：成功 {len(embeddings)} 条，失败 {len(failed)} 条")

    if use_cache and len(embeddings) > 0:
        save_custom_cache(hash_val, embeddings, ids)
        st.success("Cache saved for future use.")  # 修改为英文

    return embeddings, ids

# ---------- ID 处理 ----------
def extract_accession(full_id: str) -> str:
    if '|' in full_id:
        parts = full_id.split('|')
        if len(parts) >= 2:
            return parts[1]
    return full_id

def is_uniprot_accession(protein_id: str) -> bool:
    return bool(re.fullmatch(r'[A-Z][A-Z0-9]{5,9}', protein_id))

# ---------- 序列获取 ----------
def get_sequence_for_id(protein_id: str, fasta_text: str = None) -> Optional[str]:
    debug(f"尝试获取蛋白 {protein_id} 的序列...")
    if fasta_text:
        ids_seqs = parse_fasta(fasta_text)
        if protein_id in ids_seqs:
            debug(f"从FASTA文本中找到序列，长度={len(ids_seqs[protein_id])}")
            return ids_seqs[protein_id]
        for fid, fseq in ids_seqs.items():
            if protein_id in fid or fid in protein_id:
                debug(f"通过部分匹配找到序列: {fid}")
                return fseq

    try:
        url = f"https://rest.uniprot.org/uniprotkb/{protein_id}.fasta"
        resp = requests.get(url, timeout=10)
        if resp.status_code == 200:
            lines = resp.text.splitlines()
            seq = "".join(line.strip() for line in lines if not line.startswith(">"))
            if seq:
                debug(f"从UniProt API获取序列成功，长度={len(seq)}")
                return seq
    except Exception as e:
        debug(f"UniProt 请求失败: {str(e)}")
    return None

def parse_fasta(text: str) -> Dict[str, str]:
    ids_seqs = {}
    current_id, current_seq = None, []
    for line in text.splitlines():
        line = line.strip()
        if line.startswith(">"):
            if current_id:
                ids_seqs[current_id] = "".join(current_seq)
            current_id = line[1:].split()[0]
            current_seq = []
        else:
            if current_id:
                current_seq.append(line)
    if current_id:
        ids_seqs[current_id] = "".join(current_seq)
    return ids_seqs

# ---------- 核心检索 ----------
def find_similar_proteins(query_seq: str, top_n: int = 10,
                          model=None, batch_converter=None,
                          use_pretrained: bool = True,
                          custom_library: Dict[str, np.ndarray] = None,
                          custom_annotations: Dict[str, Dict] = None) -> pd.DataFrame:
    debug(f"find_similar_proteins 被调用: use_pretrained={use_pretrained}, "
          f"custom_library size: {len(custom_library) if custom_library else 0}")

    if query_seq is None or len(query_seq) < 10:
        return pd.DataFrame({"Error": ["Sequence too short (< 10 aa)"]})

    try:
        query_emb = get_esm_embedding(query_seq, model, batch_converter)
    except Exception as e:
        return pd.DataFrame({"Error": [f"Query embedding failed: {str(e)}"]})

    query_norm = np.linalg.norm(query_emb)
    if query_norm == 0:
        return pd.DataFrame({"Error": ["Zero embedding"]})
    query_vec = query_emb / query_norm

    ids, embeddings_norm, annotations = None, None, {}
    if use_pretrained:
        ids, embeddings_norm, annotations = load_pretrained_ref_db()
        if ids is None:
            return pd.DataFrame({"Error": ["预构建参考库未找到。"]})
    elif custom_library is not None and len(custom_library) > 0:
        ids = list(custom_library.keys())
        embeddings = np.stack([custom_library[k] for k in ids], axis=0)
        norms = np.linalg.norm(embeddings, axis=1, keepdims=True)
        norms[norms == 0] = 1
        embeddings_norm = embeddings / norms
        if custom_annotations:
            annotations = custom_annotations
    else:
        return pd.DataFrame({"Error": ["No reference library available."]})

    if embeddings_norm is None or len(embeddings_norm) == 0:
        return pd.DataFrame({"Error": ["Empty reference library."]})

    similarities = np.dot(embeddings_norm, query_vec)
    top_idx = np.argsort(similarities)[::-1][:top_n]

    result = {
        "Rank": list(range(1, len(top_idx)+1)),
        "Protein_ID": [ids[i] for i in top_idx],
        "Similarity": similarities[top_idx].tolist()
    }
    genes = [annotations.get(ids[i], {}).get("gene", "") for i in top_idx]
    names = [annotations.get(ids[i], {}).get("protein_name", "") for i in top_idx]
    gos = [annotations.get(ids[i], {}).get("go_terms", "") for i in top_idx]
    ecs = [annotations.get(ids[i], {}).get("ec_numbers", "") for i in top_idx]
    result["Gene"] = genes
    result["Protein_Name"] = names
    result["GO"] = gos
    result["EC"] = ecs
    df = pd.DataFrame(result)
    df["Similarity"] = df["Similarity"].round(4)
    debug(f"检索完成，返回 {len(df)} 条结果")
    return df

# ---------- 批量检索（原有同步版本） ----------
def batch_search_top1(query_seqs: List[str], model, batch_converter,
                      custom_library: Dict[str, np.ndarray]) -> List[Dict]:
    debug(f"[批量检索] 开始处理 {len(query_seqs)} 条查询序列。")
    if not custom_library:
        debug("[批量检索] 参考库为空")
        return []

    ids = list(custom_library.keys())
    embeddings = np.stack([custom_library[k] for k in ids], axis=0)
    norms = np.linalg.norm(embeddings, axis=1, keepdims=True)
    norms[norms == 0] = 1
    embeddings_norm = embeddings / norms

    results = []
    for i, seq in enumerate(query_seqs):
        if i % 10 == 0:
            debug(f"[批量检索] 处理 {i+1}/{len(query_seqs)}")
        try:
            q_emb = get_esm_embedding(seq, model, batch_converter)
        except Exception as e:
            debug(f"[批量检索] 序列 {i+1} 嵌入失败: {e}")
            results.append({'query_index': i, 'best_id': None, 'similarity': 0.0})
            continue

        q_norm = np.linalg.norm(q_emb)
        if q_norm == 0:
            results.append({'query_index': i, 'best_id': None, 'similarity': 0.0})
            continue

        sims = np.dot(embeddings_norm, q_emb / q_norm)
        best_idx = np.argmax(sims)
        results.append({
            'query_index': i,
            'best_id': ids[best_idx],
            'similarity': float(sims[best_idx])
        })

    debug(f"[批量检索] 完成，返回 {len(results)} 条结果。")
    return results

# ---------- 分步批量检索（新增，避免 UI 阻塞） ----------
def init_batch_search_state(session_state, total_sequences: int):
    session_state.batch_search_in_progress = True
    session_state.batch_search_index = 0
    session_state.batch_search_results = []
    session_state.batch_search_total = total_sequences
    debug(f"[分步检索] 初始化：总计 {total_sequences} 条序列")

def run_batch_search_step(query_seqs: List[str], model, batch_converter,
                          custom_library: Dict[str, np.ndarray],
                          session_state, step_size: int = 100) -> bool:
    total = session_state.batch_search_total
    start = session_state.batch_search_index
    end = min(start + step_size, total)
    debug(f"[分步检索] 处理批次 {start+1} 到 {end} / {total}")

    ids = list(custom_library.keys())
    embeddings = np.stack([custom_library[k] for k in ids], axis=0)
    norms = np.linalg.norm(embeddings, axis=1, keepdims=True)
    norms[norms == 0] = 1
    embeddings_norm = embeddings / norms

    for i in range(start, end):
        seq = query_seqs[i]
        try:
            q_emb = get_esm_embedding(seq, model, batch_converter)
        except Exception as e:
            debug(f"[分步检索] 序列 {i+1} 嵌入失败: {e}")
            session_state.batch_search_results.append({
                'query_index': i,
                'best_id': None,
                'similarity': 0.0
            })
            continue

        q_norm = np.linalg.norm(q_emb)
        if q_norm == 0:
            session_state.batch_search_results.append({
                'query_index': i,
                'best_id': None,
                'similarity': 0.0
            })
            continue

        sims = np.dot(embeddings_norm, q_emb / q_norm)
        best_idx = np.argmax(sims)
        session_state.batch_search_results.append({
            'query_index': i,
            'best_id': ids[best_idx],
            'similarity': float(sims[best_idx])
        })

    session_state.batch_search_index = end
    debug(f"[分步检索] 已完成 {end}/{total} 条")

    if end >= total:
        session_state.batch_search_in_progress = False
        debug(f"[分步检索] 全部完成，共 {len(session_state.batch_search_results)} 条结果")
        return True
    else:
        return False

# ---------- UniProt 注释（带 ID 检查） ----------
def fetch_uniprot_annotations_batch(protein_ids: List[str], batch_size: int = 50) -> Dict[str, Dict]:
    clean_ids = [extract_accession(pid) for pid in protein_ids]
    valid_count = sum(is_uniprot_accession(acc) for acc in clean_ids)
    debug(f"[UniProt 注释] 总 ID 数: {len(clean_ids)}，其中标准 UniProt accession 数量: {valid_count}")

    if valid_count < len(clean_ids) * 0.1:
        st.warning(
            "检测到参考库中绝大部分蛋白 ID 不是标准 UniProt 格式（如 `P68871`），"
            "无法从 UniProt 获取功能注释。富集分析将不可用。\n\n"
            "👉 若要体验富集分析，请上传标准 UniProt 蛋白序列（如之前提供的 5 条测试 FASTA）。"
        )
        debug("[UniProt 注释] ID 格式不符，跳过 API 请求。")
        return {}

    results = {}
    for i in range(0, len(clean_ids), batch_size):
        batch = clean_ids[i:i+batch_size]
        query = " OR ".join([f"accession:{acc}" for acc in batch])
        url = f"https://rest.uniprot.org/uniprotkb/search?query={query}&format=json&fields=accession,gene_names,protein_name,go_id,ec"
        debug(f"[UniProt 注释] 请求第 {i//batch_size+1} 批，共 {len(batch)} 个 ID。")

        try:
            resp = requests.get(url, timeout=30)
            if resp.status_code == 200:
                data = resp.json()
                for item in data.get('results', []):
                    acc = item['primaryAccession']
                    go_list = [go['id'] for go in item.get('goTerms', [])]
                    ec_list = item.get('ecNumbers', [])
                    protein_name = ""
                    try:
                        protein_name = item.get('proteinDescription', {}).get('recommendedName', {}).get('fullName', {}).get('value', '')
                    except:
                        pass
                    gene_name = ""
                    try:
                        genes = item.get('genes', [])
                        if genes:
                            gene_name = genes[0].get('geneName', {}).get('value', '')
                    except:
                        pass
                    results[acc] = {
                        'gene': gene_name,
                        'protein_name': protein_name,
                        'go_terms': "; ".join(go_list) if go_list else "",
                        'ec_numbers': "; ".join(ec_list) if ec_list else ""
                    }
                debug(f"[UniProt 注释] 第 {i//batch_size+1} 批成功获取 {len(data.get('results', []))} 条记录。")
            else:
                debug(f"[UniProt 注释] 请求失败，状态码：{resp.status_code}")
        except Exception as e:
            debug(f"[UniProt 注释] 请求异常：{str(e)}")
        time.sleep(0.5)

    debug(f"[UniProt 注释] 总共获取 {len(results)} 个蛋白的注释。")
    return results
