# utils/eggnog_annotation.py
import os
import tempfile
import shutil
import subprocess
from io import StringIO
import pandas as pd
import streamlit as st
import re

DEBUG = True

def debug(msg: str):
    if DEBUG:
        print(f"[EGGNOG] {msg}")

EMAPPER_PATH = r"D:\my_egg\eggnog-mapper-main\emapper.py"
EMAPPER_PYTHON = r"D:\my_egg\emapper2_env\Scripts\python.exe"
EMAPPER_SCRIPT = r"D:\my_egg\emapper2_env\Scripts\emapper.py"
EGGNOG_DATA_DIR = r"D:\my_egg\eggnog-mapper-main\data"
EMAPPER_BIN_DIR = r"D:\my_egg\emapper2_env\Lib\site-packages\eggnogmapper\bin"
SHORT_TEMP = r"C:\eggtmp"

def check_emapper_installed():
    missing = []
    for path in [EMAPPER_PYTHON, EMAPPER_SCRIPT, EGGNOG_DATA_DIR, EMAPPER_BIN_DIR]:
        if not os.path.exists(path):
            missing.append(path)
    required_db_files = [
        os.path.join(EGGNOG_DATA_DIR, "eggnog.db"),
        os.path.join(EGGNOG_DATA_DIR, "eggnog_proteins.dmnd"),
    ]
    for path in required_db_files:
        if not os.path.exists(path):
            missing.append(path)
    if missing:
        st.error("eggNOG environment incomplete, please check paths:")
        for path in missing:
            st.code(path)
        debug("missing paths: " + " | ".join(missing))
        return False
    return True

def run_eggnog_annotation_local(fasta_text, tax_scope="auto", target_orthologs="all"):
    debug("Starting local eggNOG-mapper v2.1.13 ...")
    if not check_emapper_installed():
        return None

    os.makedirs(SHORT_TEMP, exist_ok=True)
    fasta_path = None
    output_dir = None

    try:
        with tempfile.NamedTemporaryFile(
            mode="w", suffix=".fasta", delete=False, encoding="utf-8", dir=SHORT_TEMP
        ) as f:
            f.write(fasta_text.strip() + "\n")
            fasta_path = f.name

        output_dir = tempfile.mkdtemp(prefix="eggnog_out_", dir=SHORT_TEMP)
        project_name = os.path.splitext(os.path.basename(fasta_path))[0]

        cmd = [
            EMAPPER_PYTHON, EMAPPER_SCRIPT,
            "-i", fasta_path, "--itype", "proteins", "-m", "diamond",
            "--data_dir", EGGNOG_DATA_DIR, "-o", project_name,
            "--output_dir", output_dir, "--temp_dir", SHORT_TEMP,
            "--cpu", "4", "--dmnd_iterate", "no", "--sensmode", "fast",
            "--tax_scope", tax_scope, "--target_orthologs", target_orthologs,
            "--pfam_realign", "none", "--override",
        ]

        env = os.environ.copy()
        env["PATH"] = EMAPPER_BIN_DIR + os.pathsep + env.get("PATH", "")
        env["EGGNOG_DATA_DIR"] = EGGNOG_DATA_DIR
        env["PYTHONIOENCODING"] = "utf-8"

        progress = st.progress(0)
        status = st.empty()

        status.text("Running eggNOG-mapper (DIAMOND search in progress)...")
        progress.progress(0.1)

        proc = subprocess.Popen(
            cmd,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            text=True, encoding="utf-8", errors="replace",
            env=env, bufsize=1, universal_newlines=True,
        )

        for line in proc.stdout:
            line = line.rstrip()
            if line.startswith("CITATION:") or line.startswith("===="):
                continue
            debug(line)
            lower = line.lower()
            if "diamond" in lower or "blastp" in lower:
                status.text("Running DIAMOND alignment...")
                progress.progress(0.35)
            elif "functional annotation" in lower:
                status.text("Generating eggNOG functional annotations...")
                progress.progress(0.75)
            elif "done" in lower or "finished" in lower:
                status.text("Finishing up...")
                progress.progress(0.95)

        proc.wait()

        if proc.returncode != 0:
            st.error(f"eggNOG-mapper failed with return code {proc.returncode}")
            return None

        annotation_file = os.path.join(output_dir, f"{project_name}.emapper.annotations")
        if not os.path.exists(annotation_file):
            found = None
            for root, _, files in os.walk(output_dir):
                for file in files:
                    if file.endswith(".emapper.annotations"):
                        found = os.path.join(root, file)
                        break
                if found:
                    break
            annotation_file = found

        if not annotation_file or not os.path.exists(annotation_file):
            st.error("Annotation result file not found.")
            return None

        df = pd.read_csv(annotation_file, sep="\t", comment="#", header=None)
        possible_columns = [
            "query", "seed_ortholog", "evalue", "score", "eggNOG_OGs",
            "max_annot_lvl", "COG_category", "Description", "Preferred_name",
            "GOs", "EC", "KEGG_ko", "KEGG_Pathway", "KEGG_Module",
            "KEGG_Reaction", "KEGG_rclass", "BRITE", "KEGG_TC", "CAZy",
            "BiGG_Reaction", "PFAMs",
        ]
        if df.shape[1] <= len(possible_columns):
            df.columns = possible_columns[:df.shape[1]]
        if df.shape[1] > 0:
            df.rename(columns={df.columns[0]: "Master protein IDs"}, inplace=True)

        progress.progress(1.0)
        status.text("Annotation completed.")
        return df

    except Exception as e:
        st.error(f"eggNOG annotation error: {e}")
        return None

    finally:
        if fasta_path and os.path.exists(fasta_path):
            try: os.unlink(fasta_path)
            except: pass
        if output_dir and os.path.exists(output_dir):
            try: shutil.rmtree(output_dir, ignore_errors=True)
            except: pass

def parse_eggnog_manual_file(file_content):
    debug("Parsing manually uploaded eggNOG file...")
    lines = file_content.split("\n")
    data_lines = [line for line in lines if not line.startswith("#") and line.strip()]
    if not data_lines:
        st.error("No valid data rows found.")
        return None
    try:
        df = pd.read_csv(StringIO("\n".join(data_lines)), sep="\t")
        if df.columns[0].startswith("#"):
            df.rename(columns={df.columns[0]: "Master protein IDs"}, inplace=True)
        return df
    except Exception as e:
        st.error(f"Failed to parse TSV: {e}")
        return None

# ---------- GO 名称映射 ----------
def load_go_names(obo_path: str) -> dict:
    """解析 GO basic obo 文件，返回 {GO_ID: name} 字典。"""
    debug(f"Loading GO obo file: {obo_path}")
    go_map = {}
    if not os.path.exists(obo_path):
        st.warning(f"GO obo file not found at {obo_path}. Falling back to GO IDs.")
        debug(f"File not found: {obo_path}")
        return go_map

    current_id = None
    current_name = None
    try:
        with open(obo_path, 'r', encoding='utf-8') as f:
            for line in f:
                line = line.strip()
                if line == "[Term]":
                    if current_id and current_name:
                        go_map[current_id] = current_name
                    current_id = None
                    current_name = None
                elif line.startswith("id: "):
                    current_id = line[4:].strip()
                elif line.startswith("name: "):
                    current_name = line[6:].strip()
            if current_id and current_name:
                go_map[current_id] = current_name
        debug(f"GO obo parsing completed. Total GO terms: {len(go_map)}")
    except Exception as e:
        st.error(f"Error parsing GO obo file: {e}")
        debug(f"Error: {e}")
        return {}
    return go_map

# ---------- KO 名称映射（增强解析）----------
def load_ko_names(kegg_file: str) -> dict:
    """
    解析 TBtools.KeggBackEnd 文件，返回 {ko_id: name} 字典。
    支持两种常见格式：
      - ko:K03781\t名称
      - K03781\t名称
    增加调试输出前几行以便排查。
    """
    debug(f"Loading KEGG backend file: {kegg_file}")
    ko_map = {}
    if not os.path.exists(kegg_file):
        st.warning(f"KEGG backend file not found at {kegg_file}. Falling back to KO IDs.")
        debug(f"File not found: {kegg_file}")
        return ko_map

    try:
        with open(kegg_file, 'r', encoding='utf-8') as f:
            lines = f.readlines()
        # 调试：打印前5行
        debug("First 5 lines of KEGG file:")
        for idx, line in enumerate(lines[:5]):
            debug(f"  Line {idx+1}: {line.rstrip()}")

        for i, line in enumerate(lines):
            line = line.strip()
            if not line:
                continue
            parts = line.split('\t')
            if len(parts) < 2:
                continue
            ko_raw = parts[0].strip()
            ko_name = parts[1].strip()
            # 处理 KO ID：如果有 "ko:" 前缀，直接使用；如果没有，尝试匹配 K\d{5} 模式，自动添加 "ko:"
            if ko_raw.startswith("ko:"):
                ko_id = ko_raw
            else:
                match = re.match(r'(K\d{5})', ko_raw)
                if match:
                    ko_id = "ko:" + match.group(1)
                else:
                    # 跳过无法识别的行
                    continue
            ko_map[ko_id] = ko_name
            if i % 2000 == 0:
                debug(f"Parsed KO: {ko_id} -> {ko_name}")
        debug(f"KEGG file parsing completed. Total KO terms: {len(ko_map)}")
    except Exception as e:
        st.error(f"Error parsing KEGG backend file: {e}")
        debug(f"Error: {e}")
        return {}
    return ko_map

def get_go_name(go_id: str, go_map: dict) -> str:
    """获取 GO term 名称，如果不存在则返回原 ID。"""
    if go_id in go_map:
        return f"{go_id} ({go_map[go_id]})"
    else:
        debug(f"GO name not found for {go_id}")
        return go_id

def get_ko_name(ko_id: str, ko_map: dict) -> str:
    """获取 KO 名称，如果不存在则返回原 ID。"""
    if ko_id in ko_map:
        return f"{ko_id} ({ko_map[ko_id]})"
    else:
        debug(f"KO name not found for {ko_id}")
        return ko_id
