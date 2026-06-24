# utils/string_ppi.py
import re
import time
import json
import requests
import pandas as pd
from typing import Dict, List, Optional
import networkx as nx
from pyvis.network import Network
import tempfile
import os

DEBUG = True
def debug(msg: str):
    if DEBUG:
        print(f"[STRING PPI] {msg}")

def extract_uniprot_from_eggnog(eggnog_df: pd.DataFrame, target_ids: List[str], fallback_to_seed: bool = True) -> Dict[str, str]:
    """
    从 eggNOG 结果中提取 STRING 可用的标识符。
    优先级：
    1. sp|ACCESSION| 格式中的 UniProt accession
    2. 任何匹配 UniProt accession 模式的独立单词
    3. 如果 fallback_to_seed=True，尝试从 seed_ortholog 中进一步提取：
       - 如 "4565.Traes_7DL_44F6042FE.1"，尝试取点号后的部分（可能不是 UniProt，但作为最后手段）
    若以上均无法获得有效 ID，则不返回该蛋白的映射。
    """
    debug("Extracting identifiers from eggNOG...")
    mapping = {}
    if eggnog_df is None or eggnog_df.empty:
        debug("eggNOG DataFrame is empty.")
        return mapping

    id_col = eggnog_df.columns[0]
    seed_col = None
    for col in eggnog_df.columns:
        if 'seed_ortholog' in col.lower():
            seed_col = col
            break
    if seed_col is None:
        debug("No seed_ortholog column found.")
        return mapping

    # UniProt accession 模式
    uniprot_pattern = r'\b([OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9][A-Z][A-Z0-9]{2}[0-9](?:[A-Z][A-Z0-9]{2}[0-9])?)\b'

    for _, row in eggnog_df.iterrows():
        query = str(row[id_col]).strip()
        if query not in target_ids:
            continue
        seed = str(row[seed_col]).strip()
        debug(f"Processing {query}: seed_ortholog = {seed}")

        # 1) sp|ACCESSION| 格式
        sp_matches = re.findall(r'sp\|([A-Z0-9]+)\|', seed)
        if sp_matches:
            acc = sp_matches[0]
            if re.fullmatch(uniprot_pattern, acc):
                mapping[query] = acc
                debug(f"  -> {acc} (from sp|...)")
                continue

        # 2) 匹配 UniProt accession 模式
        candidates = re.findall(uniprot_pattern, seed)
        good = [c for c in candidates if len(c) >= 6]
        if good:
            best = max(good, key=len)
            mapping[query] = best
            debug(f"  -> {best} (pattern match)")
            continue

        # 3) 回退：尝试从 seed_ortholog 中提取可能的部分（非标准，但有时可用）
        if fallback_to_seed:
            # 例如 "4565.Traes_7DL_44F6042FE.1"，取点后面的部分，通常是一个基因 ID
            parts = seed.split('.')
            if len(parts) >= 2:
                # 取第二部分，但需要排除纯数字物种 ID
                # 如果第一部分是数字（物种 ID），则第二部分可能是基因 ID
                if parts[0].isdigit():
                    potential = '.'.join(parts[1:])  # 取物种 ID 后面的全部
                else:
                    potential = seed
                # 尝试匹配 UniProt accession 模式
                if re.search(uniprot_pattern, potential):
                    mapping[query] = potential
                    debug(f"  -> {potential} (seed_ortholog fallback)")
                else:
                    debug(f"  -> no valid UniProt pattern in seed_ortholog, skipping")
                    # 不添加映射，让 STRING 使用其他方法
            else:
                debug(f"  -> seed_ortholog has only one part, skipping")
        else:
            debug(f"  -> no valid UniProt found and fallback disabled, skipping")

    debug(f"Extraction complete: {len(mapping)} mappings.")
    return mapping

def run_blast_mapping(sequence_dict: Dict[str, str], max_retries: int = 3, max_wait_per_job: int = 180) -> Dict[str, str]:
    debug(f"Starting async BLAST mapping for {len(sequence_dict)} sequences")
    mapping = {}
    base_url = "https://rest.uniprot.org"
    submit_url = f"{base_url}/blast/run"
    for sid, seq in sequence_dict.items():
        success = False
        for attempt in range(max_retries):
            try:
                form_data = {
                    "query": seq,
                    "format": "json",
                    "matrix": "BLOSUM62",
                    "gapopen": "11",
                    "gapext": "1",
                    "threshold": "1e-5",
                    "alignments": "1"
                }
                resp_submit = requests.post(submit_url, data=form_data,
                                            headers={"User-Agent": "ProteomicsApp/1.0"}, timeout=30)
                if resp_submit.status_code != 200:
                    time.sleep(10)
                    continue
                job_id = resp_submit.json().get("jobId")
                if not job_id:
                    time.sleep(10)
                    continue
                status_url = f"{base_url}/blast/status/{job_id}"
                waited = 0
                while waited < max_wait_per_job:
                    time.sleep(10)
                    waited += 10
                    resp_status = requests.get(status_url, timeout=10)
                    if resp_status.status_code != 200:
                        continue
                    status_json = resp_status.json()
                    if status_json.get("jobStatus") == "FINISHED":
                        break
                else:
                    continue
                result_url = f"{base_url}/blast/result/{job_id}"
                resp_result = requests.get(result_url, params={"format": "json"}, timeout=30)
                if resp_result.status_code == 200:
                    result = resp_result.json()
                    hits = result.get("hits", [])
                    if hits:
                        accession = hits[0].get("acc", "")
                        if accession:
                            mapping[sid] = accession
                            success = True
                            break
            except Exception as e:
                debug(f"  BLAST error: {e}")
            time.sleep(5)
        if not success:
            debug(f"  {sid} BLAST failed completely")
        time.sleep(1)
    return mapping

def call_string_api(identifiers: List[str], species: int = 9606, required_score: int = 400,
                    network_type: str = "functional", add_nodes: int = 0) -> Optional[List]:
    if not identifiers:
        return None
    api_url = "https://string-db.org/api/json/network"
    params = {
        "identifiers": "\r".join(identifiers),
        "species": species,
        "required_score": required_score,
        "network_type": network_type,
        "add_nodes": add_nodes,
        "show_query_node_labels": 1
    }
    debug(f"Querying STRING with species={species}, score={required_score}, {len(identifiers)} IDs")
    try:
        resp = requests.post(api_url, data=params, timeout=30)
        if resp.status_code == 200:
            data = resp.json()
            return data
        else:
            debug(f"STRING returned {resp.status_code}: {resp.text[:200]}")
            return None
    except Exception as e:
        debug(f"STRING request failed: {e}")
        return None

def build_ppi_network_html(network_data: List, id_mapping: Optional[Dict[str, str]] = None) -> str:
    debug(f"Building network with {len(network_data)} edges")
    G = nx.Graph()
    for edge in network_data:
        nodeA = edge.get("preferredName_A") or edge.get("stringId_A")
        nodeB = edge.get("preferredName_B") or edge.get("stringId_B")
        score = float(edge.get("score", 0))
        if nodeA and nodeB:
            G.add_edge(nodeA, nodeB, weight=score, title=f"Score: {score:.3f}")
    net = Network(height="700px", width="100%", notebook=False, directed=False)
    net.set_options("""
    var options = {
      "nodes": {"color": {"border": "rgba(0,0,0,0.5)", "background": "rgba(72,199,142,0.8)", "highlight": {"border": "rgba(0,0,0,1)", "background": "rgba(72,199,142,1)"}}, "font": {"size": 14, "face": "Arial"}},
      "edges": {"color": {"color": "rgba(150,150,150,0.6)", "highlight": "rgba(150,150,150,1)"}, "smooth": {"type": "continuous"}},
      "physics": {"barnesHut": {"gravitationalConstant": -8000, "springLength": 150}}
    }
    """)
    for node in G.nodes():
        label = node
        if id_mapping:
            reverse_map = {v: k for k, v in id_mapping.items()}
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
        return html
    except Exception as e:
        debug(f"Network HTML generation failed: {e}")
        return "<h3>Network generation failed</h3>"
