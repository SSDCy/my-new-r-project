# utils/build_ref_db.py
"""
离线构建 Swiss-Prot 参考蛋白嵌入库（修复维度与分页问题）
"""

import numpy as np
import pickle
import os
import requests
import time
from tqdm import tqdm
from Bio import SeqIO
import torch
import esm

DEBUG = True
def debug(msg):
    if DEBUG:
        print(f"[BUILD REF DB] {msg}")

def download_swissprot_fasta():
    url = "https://rest.uniprot.org/uniprotkb/stream?format=fasta&query=%28reviewed:true%29"
    debug("开始下载 Swiss-Prot FASTA...")
    response = requests.get(url, stream=True)
    fasta_path = "swissprot.fasta"
    with open(fasta_path, "wb") as f:
        for chunk in response.iter_content(chunk_size=8192):
            f.write(chunk)
    debug(f"下载完成，文件大小：{os.path.getsize(fasta_path)} bytes")
    return fasta_path

def download_swissprot_annotations():
    """分页获取全部注释，返回 {accession: {...}}"""
    debug("开始下载注释信息...")
    base_url = "https://rest.uniprot.org/uniprotkb/search"
    params = {
        "query": "reviewed:true",
        "fields": "accession,gene_names,protein_name,go_id,ec",
        "size": 500,
        "format": "json"
    }
    annotations = {}
    next_url = base_url + "?query=reviewed:true&fields=accession,gene_names,protein_name,go_id,ec&size=500&format=json"
    count = 0
    while next_url:
        debug(f"请求：{next_url}")
        resp = requests.get(next_url)
        if resp.status_code != 200:
            debug(f"请求失败，状态码：{resp.status_code}")
            break
        data = resp.json()
        for item in data.get("results", []):
            acc = item["primaryAccession"]
            genes = item.get("genes", [])
            gene = genes[0]["geneName"]["value"] if genes else None
            protein_name = None
            if "proteinDescription" in item:
                recommended = item["proteinDescription"].get("recommendedName")
                if recommended:
                    protein_name = recommended.get("fullName", {}).get("value")
            go_terms = []
            if "goTerms" in item:
                go_terms = list(set(term["id"] for term in item["goTerms"]))
            ec_numbers = item.get("ecNumbers", [])
            annotations[acc] = {
                "gene": gene,
                "protein_name": protein_name,
                "go_terms": "; ".join(go_terms),
                "ec_numbers": "; ".join(ec_numbers)
            }
        count += len(data["results"])
        debug(f"已获取 {count} 条注释")
        next_link = data.get("links", {}).get("next")
        next_url = next_link if next_link else None
        time.sleep(0.5)
    debug(f"注释下载完成，共 {len(annotations)} 条")
    return annotations

def build_reference_db():
    debug("===== 开始构建参考库 =====")
    os.makedirs("data", exist_ok=True)

    if not os.path.exists("swissprot.fasta"):
        fasta_file = download_swissprot_fasta()
    else:
        fasta_file = "swissprot.fasta"
        debug("使用已有 FASTA 文件。")

    debug("解析 FASTA...")
    records = list(SeqIO.parse(fasta_file, "fasta"))
    ids = [rec.id.split('|')[1] if '|' in rec.id else rec.id for rec in records]
    sequences = [str(rec.seq) for rec in records]
    debug(f"共 {len(ids)} 条序列")

    if not os.path.exists("data/ref_annotations.pkl"):
        annotations = download_swissprot_annotations()
        with open("data/ref_annotations.pkl", "wb") as f:
            pickle.dump(annotations, f)
    else:
        debug("注释文件已存在，跳过下载。")
        with open("data/ref_annotations.pkl", "rb") as f:
            annotations = pickle.load(f)

    debug("加载 ESM2 模型...")
    model, alphabet = esm.pretrained.esm2_t12_35M_UR50D()
    batch_converter = alphabet.get_batch_converter()
    model.eval()

    # 动态获取嵌入维度
    test_seq = sequences[0]
    data = [("test", test_seq)]
    _, _, batch_tokens = batch_converter(data)
    with torch.no_grad():
        results = model(batch_tokens, repr_layers=[12])
    embedding_dim = results["representations"][12].shape[-1]
    debug(f"嵌入维度：{embedding_dim}")

    debug("生成嵌入向量...")
    embeddings = np.zeros((len(ids), embedding_dim), dtype=np.float32)
    for i in tqdm(range(len(sequences)), desc="生成嵌入"):
        seq = sequences[i]
        data = [("protein", seq)]
        _, _, batch_tokens = batch_converter(data)
        with torch.no_grad():
            results = model(batch_tokens, repr_layers=[12])
        token_embeddings = results["representations"][12]
        embedding = token_embeddings.mean(0).cpu().numpy()
        embeddings[i] = embedding

    norms = np.linalg.norm(embeddings, axis=1, keepdims=True)
    embeddings_norm = embeddings / norms

    np.savez("data/ref_db.npz", ids=ids, embeddings=embeddings_norm)
    debug("嵌入库已保存到 data/ref_db.npz")
    debug("===== 构建完成 =====")

if __name__ == "__main__":
    build_reference_db()
