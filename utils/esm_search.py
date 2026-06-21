# utils/esm_search.py
"""
ESM2 蛋白质功能嵌入搜索模块
- 支持加载预构建参考库
- 支持在线构建自定义参考库（带缓存）
- 支持根据蛋白ID自动获取序列
- 提供单条和批量相似蛋白检索
- 提供 UniProt 注释批量获取（带 ID 格式检查）
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
        data = np.load(db_path, allow_pickle=True)
        ids = data["ids"].tolist()
        embeddings_norm = data["embeddings"]
        with open(ann_path, "rb") as f:
            annotations = pickle.load(f)
        debug(f"参考库加载完成：{len(ids)} 条蛋白")
        return ids, embeddings_norm, annotations
    else:
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
        batch_converter = alphabet.get_batch_converter()
        return model, alphabet, batch_converter
    except Exception as e:
        debug(f"模型加载失败: {str(e)}")
        st.error(f"模型加载失败: {str(e)}")
        return None, None, None

def get_esm_embedding(sequence: str, model, batch_converter) -> np.ndarray:
    if not ESM_AVAILABLE:
        raise RuntimeError("ESM2 不可用")
    data = [("query", sequence)]
    _, _, batch_tokens = batch_converter(data)
    if torch.cuda.is_available():
        batch_tokens = batch_tokens.cuda()
    with torch.no_grad():
        results = model(batch_tokens, repr_layers=[12])
    token_embeddings = results["representations"][12]
    embedding = token_embeddings.mean(dim=1).squeeze(0).cpu().numpy()
    return embedding

# ---------- 缓存 ----------
def compute_fasta_hash(fasta_text: str) -> str:
    return hashlib.sha256(fasta_text.encode('utf-8')).hexdigest()

def get_cache_path(hash_val: str) -> str:
    os.makedirs("data", exist_ok=True)
    return f"data/custom_cache_{hash_val}.npz"

def load_custom_cache(hash_val: str) -> Optional[Dict[str, np.ndarray]]:
    cache_file = get_cache_path(hash_val)
    if os.path.exists(cache_file):
        debug(f"发现缓存文件: {cache_file}")
        data = np.load(cache_file, allow_pickle=True)
        ids = data["ids"].tolist()
        embeddings = data["embeddings"]
        lib = {id_: embeddings[i] for i, id_ in enumerate(ids)}
        return lib
    return None

def save_custom_cache(hash_val: str, embeddings_dict: Dict[str, np.ndarray], ids: List[str]):
    cache_file = get_cache_path(hash_val)
    embeddings = np.stack([embeddings_dict[id_] for id_ in ids], axis=0)
    np.savez(cache_file, ids=np.array(ids), embeddings=embeddings)

# ---------- 构建自定义库 ----------
def build_reference_library(fasta_text: str, model, batch_converter,
                            use_cache: bool = True) -> Tuple[Dict[str, np.ndarray], List[str]]:
    debug("开始从 FASTA 构建参考库...")
    if use_cache:
        hash_val = compute_fasta_hash(fasta_text)
        debug(f"FASTA 哈希值: {hash_val}")
        cached = load_custom_cache(hash_val)
        if cached is not None:
            st.success("✅ 已从缓存加载嵌入库，无需重新生成。")
            return cached, list(cached.keys())
        else:
            debug("未找到缓存，将重新生成嵌入。")

    # 解析
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
        return {}, []

    embeddings = {}
    failed = []
    for i, (seq_id, seq) in enumerate(zip(ids, sequences)):
        debug(f"处理 {i+1}/{len(ids)}: {seq_id} (长度={len(seq)})")
        try:
            embeddings[seq_id] = get_esm_embedding(seq, model, batch_converter)
        except Exception as e:
            debug(f"生成嵌入失败 ID={seq_id}: {str(e)}")
            failed.append(seq_id)

    if failed:
        st.warning(f"以下蛋白嵌入生成失败: {', '.join(failed)}")
    debug(f"自定义库构建完成：成功 {len(embeddings)} 条，失败 {len(failed)} 条")

    if use_cache and len(embeddings) > 0:
        save_custom_cache(hash_val, embeddings, ids)
        st.success("💾 嵌入库已保存到本地缓存，下次上传相同文件时将自动加载。")
    return embeddings, ids

# ---------- ID 处理 ----------
def extract_accession(full_id: str) -> str:
    if '|' in full_id:
        parts = full_id.split('|')
        if len(parts) >= 2:
            return parts[1]
    return full_id

def is_uniprot_accession(protein_id: str) -> bool:
    """判断ID是否为UniProt Accession格式（如P68871, Q8NHL6, A0A0B4J2F0）"""
    # 模式：通常6-10个字符，包含数字和字母，以字母开头
    return bool(re.fullmatch(r'[A-Z][A-Z0-9]{5,9}', protein_id))

# ---------- 序列获取 ----------
def get_sequence_for_id(protein_id: str, fasta_text: str = None) -> Optional[str]:
    debug(f"尝试获取蛋白 {protein_id} 的序列...")
    if fasta_text:
        ids_seqs = parse_fasta(fasta_text)
        if protein_id in ids_seqs:
            return ids_seqs[protein_id]
    try:
        url = f"https://rest.uniprot.org/uniprotkb/{protein_id}.fasta"
        resp = requests.get(url, timeout=10)
        if resp.status_code == 200:
            lines = resp.text.splitlines()
            seq = "".join(line.strip() for line in lines if not line.startswith(">"))
            if seq:
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
    # ...（保持不变，与之前相同）...
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
        embeddings_norm = embeddings / norms
        if custom_annotations:
            annotations = custom_annotations
    else:
        return pd.DataFrame({"Error": ["No reference library available."]})

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
    return df

# ---------- 批量检索 ----------
def batch_search_top1(query_seqs: List[str], model, batch_converter,
                      custom_library: Dict[str, np.ndarray]) -> List[Dict]:
    debug(f"[批量检索] 开始处理 {len(query_seqs)} 条查询序列。")
    if not custom_library:
        return []
    ids = list(custom_library.keys())
    embeddings = np.stack([custom_library[k] for k in ids], axis=0)
    norms = np.linalg.norm(embeddings, axis=1, keepdims=True)
    embeddings_norm = embeddings / norms
    results = []
    for i, seq in enumerate(query_seqs):
        debug(f"[批量检索] 处理 {i+1}/{len(query_seqs)}，序列长度={len(seq)}")
        try:
            q_emb = get_esm_embedding(seq, model, batch_converter)
        except:
            results.append({'query_index': i, 'best_id': None, 'similarity': 0.0})
            continue
        q_norm = np.linalg.norm(q_emb)
        if q_norm == 0:
            results.append({'query_index': i, 'best_id': None, 'similarity': 0.0})
            continue
        sims = np.dot(embeddings_norm, q_emb / q_norm)
        best_idx = np.argmax(sims)
        results.append({'query_index': i, 'best_id': ids[best_idx], 'similarity': float(sims[best_idx])})
    debug(f"[批量检索] 完成，返回 {len(results)} 条结果。")
    return results

# ---------- UniProt 注释（带 ID 检查） ----------
def fetch_uniprot_annotations_batch(protein_ids: List[str], batch_size: int = 50) -> Dict[str, Dict]:
    """
    批量获取 UniProt 注释。
    自动检测 ID 格式：如果超过 90% 不是标准 UniProt Accession，则直接返回空，
    避免大量无效 API 请求。
    """
    # 提取纯净 accession
    clean_ids = [extract_accession(pid) for pid in protein_ids]
    # 检查格式
    valid_count = sum(is_uniprot_accession(acc) for acc in clean_ids)
    debug(f"[UniProt 注释] 总 ID 数: {len(clean_ids)}，其中标准 UniProt accession 数量: {valid_count}")
    if valid_count < len(clean_ids) * 0.1:  # 少于10%
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
                    protein_name = item.get('proteinDescription', {}).get('recommendedName', {}).get('fullName', {}).get('value', '')
                    results[acc] = {
                        'go': go_list,
                        'ec': ec_list,
                        'protein_name': protein_name
                    }
                debug(f"[UniProt 注释] 第 {i//batch_size+1} 批成功获取 {len(data.get('results', []))} 条记录。")
            else:
                debug(f"[UniProt 注释] 请求失败，状态码：{resp.status_code}")
        except Exception as e:
            debug(f"[UniProt 注释] 请求异常：{str(e)}")
        time.sleep(0.5)
    debug(f"[UniProt 注释] 总共获取 {len(results)} 个蛋白的注释。")
    return results
