# utils/string_ppi.py

import os
import re
import time
import json
import tempfile
from typing import Dict, List, Optional

import requests
import pandas as pd
import networkx as nx
from pyvis.network import Network

DEBUG = True


def debug(msg: str):
    if DEBUG:
        print(f"[STRING PPI] {msg}")


UNIPROT_ACC_RE = re.compile(
    r"\b(?:[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9][A-Z][A-Z0-9]{2}[0-9](?:[A-Z][A-Z0-9]{2}[0-9])?)\b"
)

STRING_BASE = "https://version-12-0.string-db.org/api/json"
EBI_BLAST_BASE = "https://www.ebi.ac.uk/Tools/services/rest/ncbiblast"
CALLER_IDENTITY = "ppi_annotation_project"


def _extract_accessions_from_text(text: str) -> List[str]:
    if not text:
        return []
    return list(dict.fromkeys(UNIPROT_ACC_RE.findall(str(text))))


# ---------- Extract UniProt IDs from eggNOG annotation ----------
def extract_uniprot_from_eggnog(eggnog_df: pd.DataFrame, target_ids: List[str]) -> Dict[str, str]:
    """
    Extract UniProt accessions from eggNOG annotation DataFrame.

    Returns:
        {Original_ID: UniProt_Accession}
    """
    debug("开始从 eggNOG 提取 UniProt ID...")
    mapping: Dict[str, str] = {}

    if eggnog_df is None or eggnog_df.empty:
        debug("eggNOG DataFrame 为空")
        return mapping

    if not target_ids:
        debug("target_ids 为空")
        return mapping

    target_set = set(str(x).strip() for x in target_ids)

    id_col = eggnog_df.columns[0]
    candidate_cols = []

    preferred_names = [
        "seed_ortholog",
        "Preferred_name",
        "Description",
        "GOs",
        "KEGG_ko",
        "Master protein IDs",
        "query",
    ]

    for name in preferred_names:
        for col in eggnog_df.columns:
            if str(col).lower() == name.lower():
                candidate_cols.append(col)

    for col in eggnog_df.columns:
        if col not in candidate_cols:
            candidate_cols.append(col)

    debug(f"使用 ID 列: {id_col}")
    debug("候选提取列: " + ", ".join(map(str, candidate_cols[:8])))

    for _, row in eggnog_df.iterrows():
        query = str(row[id_col]).strip()
        if query not in target_set:
            continue

        for col in candidate_cols:
            value = row.get(col, "")
            matches = _extract_accessions_from_text(value)
            if matches:
                mapping[query] = matches[0]
                debug(f"  {query} -> {matches[0]} (from {col})")
                break

    debug(f"提取完成，共 {len(mapping)} 个映射")
    return mapping


# ---------- BLAST mapping via EMBL-EBI Job Dispatcher ----------
def _submit_ebi_blast(seq_id: str, sequence: str) -> Optional[str]:
    fasta = f">{seq_id}\n{sequence}\n"
    data = {
        "email": "anonymous@example.com",
        "title": f"ppi_annotation_{seq_id}",
        "program": "blastp",
        "stype": "protein",
        "database": "uniprotkb",
        "sequence": fasta,
    }

    try:
        resp = requests.post(
            f"{EBI_BLAST_BASE}/run",
            data=data,
            headers={"User-Agent": "ppi_annotation_project/1.0"},
            timeout=60,
        )
        if resp.status_code not in (200, 202):
            debug(f"    BLAST 提交失败 {resp.status_code}: {resp.text[:300]}")
            return None

        job_id = resp.text.strip()
        if not job_id or "<html" in job_id.lower():
            debug(f"    BLAST 未返回有效 job id: {job_id[:200]}")
            return None

        debug(f"    BLAST job id: {job_id}")
        return job_id

    except Exception as e:
        debug(f"    BLAST 提交异常: {e}")
        return None


def _wait_ebi_blast(job_id: str, max_wait: int = 180) -> bool:
    waited = 0
    while waited < max_wait:
        time.sleep(5)
        waited += 5

        try:
            resp = requests.get(f"{EBI_BLAST_BASE}/status/{job_id}", timeout=30)
            status = resp.text.strip().upper()
            debug(f"    BLAST status {waited}s: {status}")

            if status == "FINISHED":
                return True
            if status in {"ERROR", "FAILURE", "NOT_FOUND"}:
                return False

        except Exception as e:
            debug(f"    BLAST 状态查询异常: {e}")

    debug("    BLAST 等待超时")
    return False


def _get_ebi_blast_result(job_id: str) -> Optional[str]:
    # EBI supports several result types. out is pairwise text; tsv is easier if enabled.
    for result_type in ["tsv", "out"]:
        try:
            resp = requests.get(
                f"{EBI_BLAST_BASE}/result/{job_id}/{result_type}",
                timeout=60,
            )
            if resp.status_code == 200 and resp.text.strip():
                return resp.text
            debug(f"    获取 BLAST {result_type} 失败 {resp.status_code}: {resp.text[:200]}")
        except Exception as e:
            debug(f"    获取 BLAST {result_type} 异常: {e}")
    return None


def _parse_blast_result_for_accession(text: str) -> Optional[str]:
    if not text:
        return None

    matches = _extract_accessions_from_text(text)
    if matches:
        return matches[0]

    # Some EBI outputs contain sp|P12345|NAME or tr|A0A...|NAME
    m = re.search(r"\b(?:sp|tr)\|([A-Z0-9]+)\|", text)
    if m:
        return m.group(1)

    return None


def run_blast_mapping(sequence_dict: Dict[str, str], max_retries: int = 2) -> Dict[str, str]:
    """
    Map protein sequences to UniProt accessions using EMBL-EBI NCBI BLAST REST API.
    Returns:
        {Original_ID: UniProt_Accession}
    """
    debug(f"开始异步 BLAST 映射，序列数 {len(sequence_dict)}")
    mapping: Dict[str, str] = {}

    for seq_id, seq in sequence_dict.items():
        seq = re.sub(r"\s+", "", str(seq))
        if not seq:
            debug(f"  {seq_id} 序列为空，跳过")
            continue

        debug(f"  BLAST 查询: {seq_id}, 长度 {len(seq)}")

        success = False
        for attempt in range(1, max_retries + 1):
            debug(f"    尝试 {attempt}/{max_retries}")
            job_id = _submit_ebi_blast(seq_id, seq)
            if not job_id:
                time.sleep(3)
                continue

            if not _wait_ebi_blast(job_id):
                time.sleep(3)
                continue

            result_text = _get_ebi_blast_result(job_id)
            accession = _parse_blast_result_for_accession(result_text or "")
            if accession:
                mapping[seq_id] = accession
                debug(f"    {seq_id} -> {accession}")
                success = True
                break

            debug(f"    {seq_id} 未从 BLAST 结果中提取到 accession")
            success = True
            break

        if not success:
            debug(f"    {seq_id} BLAST 完全失败")

        time.sleep(1)

    debug(f"异步 BLAST 完成，得到 {len(mapping)} 个映射")
    return mapping


# ---------- STRING API ----------
def _map_to_string_ids(identifiers: List[str], species: int) -> List[str]:
    unique_ids = list(dict.fromkeys(str(x).strip() for x in identifiers if str(x).strip()))
    if not unique_ids:
        return []

    url = f"{STRING_BASE}/get_string_ids"
    data = {
        "identifiers": "\n".join(unique_ids),
        "species": species,
        "caller_identity": CALLER_IDENTITY,
    }

    debug(f"发送 STRING ID 映射请求: {url}")
    debug(f"species={species}, input IDs={len(unique_ids)}")

    try:
        resp = requests.post(url, data=data, timeout=60)
        if resp.status_code != 200:
            debug(f"STRING ID 映射失败 {resp.status_code}: {resp.text[:300]}")
            return []

        rows = resp.json()
        string_ids = []
        for row in rows:
            sid = row.get("stringId")
            if sid:
                string_ids.append(sid)

        string_ids = list(dict.fromkeys(string_ids))
        debug(f"STRING ID 映射成功: {len(string_ids)}")
        return string_ids

    except Exception as e:
        debug(f"STRING ID 映射异常: {e}")
        return []


def call_string_api(
    identifiers: List[str],
    species: int = 9606,
    required_score: int = 400,
    network_type: str = "functional",
    add_nodes: int = 0,
) -> Optional[List]:
    """
    Call STRING API and return network JSON.
    Input identifiers may be UniProt accessions; they are mapped to STRING IDs first.
    """
    if not identifiers:
        debug("STRING API: 没有输入 ID")
        return None

    string_ids = _map_to_string_ids(identifiers, species)
    if not string_ids:
        debug("STRING API: 未映射到 STRING IDs")
        return None

    url = f"{STRING_BASE}/network"
    data = {
        "identifiers": "\n".join(string_ids),
        "species": species,
        "required_score": required_score,
        "network_type": network_type,
        "add_nodes": add_nodes,
        "caller_identity": CALLER_IDENTITY,
    }

    debug(f"发送 STRING network 请求: {url}")
    debug(f"species={species}, score={required_score}, STRING IDs={len(string_ids)}")

    try:
        resp = requests.post(url, data=data, timeout=60)
        if resp.status_code != 200:
            debug(f"STRING network 请求失败 {resp.status_code}: {resp.text[:300]}")
            return None

        rows = resp.json()
        debug(f"收到 STRING network 数据: {len(rows)} 条")
        return rows

    except Exception as e:
        debug(f"STRING network 请求异常: {e}")
        return None


# ---------- Build pyvis HTML ----------
def build_ppi_network_html(network_data: List, id_mapping: Optional[Dict[str, str]] = None) -> str:
    """
    Build interactive pyvis network HTML from STRING API JSON.
    """
    debug("开始构建网络图")
    G = nx.Graph()

    for edge in network_data:
        node_a = edge.get("preferredName_A") or edge.get("stringId_A")
        node_b = edge.get("preferredName_B") or edge.get("stringId_B")
        score = float(edge.get("score", 0))

        if node_a and node_b:
            G.add_edge(node_a, node_b, weight=score, title=f"Score: {score:.3f}")

    debug(f"图节点数: {G.number_of_nodes()}, 边数: {G.number_of_edges()}")

    if G.number_of_nodes() == 0:
        return "<h3>No STRING interactions found for the selected inputs.</h3>"

    net = Network(height="700px", width="100%", notebook=False, directed=False)
    net.set_options("""
    var options = {
      "nodes": {
        "color": {
          "border": "rgba(0,0,0,0.5)",
          "background": "rgba(72,199,142,0.8)",
          "highlight": {
            "border": "rgba(0,0,0,1)",
            "background": "rgba(72,199,142,1)"
          }
        },
        "font": {
          "size": 14,
          "face": "Arial"
        }
      },
      "edges": {
        "color": {
          "color": "rgba(150,150,150,0.6)",
          "highlight": "rgba(150,150,150,1)"
        },
        "smooth": {
          "type": "continuous"
        }
      },
      "physics": {
        "barnesHut": {
          "gravitationalConstant": -8000,
          "springLength": 150
        }
      }
    }
    """)

    reverse_map = {v: k for k, v in (id_mapping or {}).items()}

    for node in G.nodes():
        label = node
        if node in reverse_map:
            label = f"{node} ({reverse_map[node]})"
        net.add_node(node, label=label, title=label)

    for u, v, data in G.edges(data=True):
        net.add_edge(u, v, title=data.get("title", ""), value=data.get("weight", 1))

    try:
        tmpfile = tempfile.NamedTemporaryFile(delete=False, suffix=".html")
        net.save_graph(tmpfile.name)
        tmpfile.close()

        with open(tmpfile.name, "r", encoding="utf-8") as f:
            html = f.read()

        os.unlink(tmpfile.name)
        debug("网络图 HTML 生成成功")
        return html

    except Exception as e:
        debug(f"生成网络图失败: {e}")
        return f"<h3>网络图生成失败: {e}</h3>"