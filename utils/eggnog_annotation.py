# utils/eggnog_annotation.py

import os
import tempfile
import shutil
import subprocess
from io import StringIO

import pandas as pd
import streamlit as st

DEBUG = True


def debug(msg: str):
    if DEBUG:
        print(f"[EGGNOG] {msg}")


# Keep this for app.py database status check:
# app.py uses os.path.dirname(EMAPPER_PATH) + "\\data"
EMAPPER_PATH = r"D:\my_egg\eggnog-mapper-main\emapper.py"

# Real working eggNOG-mapper v2 environment
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
        st.error("eggNOG 本地环境不完整，请检查以下路径：")
        for path in missing:
            st.code(path)
        debug("missing paths: " + " | ".join(missing))
        return False

    debug(f"emapper python: {EMAPPER_PYTHON}")
    debug(f"emapper script: {EMAPPER_SCRIPT}")
    debug(f"data dir: {EGGNOG_DATA_DIR}")
    debug(f"bin dir: {EMAPPER_BIN_DIR}")
    return True


def run_eggnog_annotation_local(fasta_text, tax_scope="auto", target_orthologs="all"):
    """Run local eggNOG-mapper v2.1.13 and return annotation DataFrame."""
    debug("开始本地注释流程：使用 eggNOG-mapper v2.1.13 + fast/low-memory mode")

    if not check_emapper_installed():
        return None

    os.makedirs(SHORT_TEMP, exist_ok=True)

    fasta_path = None
    output_dir = None

    try:
        with tempfile.NamedTemporaryFile(
            mode="w",
            suffix=".fasta",
            delete=False,
            encoding="utf-8",
            dir=SHORT_TEMP,
        ) as f:
            f.write(fasta_text.strip() + "\n")
            fasta_path = f.name

        output_dir = tempfile.mkdtemp(prefix="eggnog_out_", dir=SHORT_TEMP)
        project_name = os.path.splitext(os.path.basename(fasta_path))[0]

        debug(f"query FASTA: {fasta_path}")
        debug(f"output dir: {output_dir}")

        cmd = [
            EMAPPER_PYTHON,
            EMAPPER_SCRIPT,
            "-i", fasta_path,
            "--itype", "proteins",
            "-m", "diamond",
            "--data_dir", EGGNOG_DATA_DIR,
            "-o", project_name,
            "--output_dir", output_dir,
            "--temp_dir", SHORT_TEMP,
            "--cpu", "4",
            "--dmnd_iterate", "no",
            "--sensmode", "fast",
            "--tax_scope", tax_scope,
            "--target_orthologs", target_orthologs,
            "--pfam_realign", "none",
            "--override",
        ]

        debug("emapper command: " + " ".join(cmd))

        env = os.environ.copy()
        env["PATH"] = EMAPPER_BIN_DIR + os.pathsep + env.get("PATH", "")
        env["EGGNOG_DATA_DIR"] = EGGNOG_DATA_DIR
        env["PYTHONIOENCODING"] = "utf-8"

        progress = st.progress(0)
        status = st.empty()
        log_box = st.empty()

        status.text("正在运行 eggNOG-mapper 本地注释...")
        progress.progress(0.1)

        proc = subprocess.Popen(
            cmd,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            text=True,
            encoding="utf-8",
            errors="replace",
            env=env,
            bufsize=1,
            universal_newlines=True,
        )

        logs = []
        for line in proc.stdout:
            line = line.rstrip()
            if not line:
                continue

            logs.append(line)
            debug(line)

            lower = line.lower()
            if "diamond" in lower or "blastp" in lower:
                status.text("正在运行 DIAMOND 搜索...")
                progress.progress(0.35)
            elif "functional annotation" in lower:
                status.text("正在生成 eggNOG 功能注释...")
                progress.progress(0.75)
            elif "done" in lower or "finished" in lower:
                status.text("注释流程即将完成...")
                progress.progress(0.95)

            if len(logs) % 5 == 0:
                log_box.code("\n".join(logs[-12:]))

        proc.wait()

        if proc.returncode != 0:
            st.error(f"eggNOG-mapper 运行失败，返回码：{proc.returncode}")
            log_box.code("\n".join(logs[-80:]) if logs else "无输出")
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
            st.error("未找到 eggNOG 注释结果文件。")
            log_box.code("\n".join(logs[-80:]) if logs else "无输出")
            return None

        debug(f"reading annotation file: {annotation_file}")

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
        status.text("注释完成。")
        debug("本地注释流程完成")
        return df

    except Exception as e:
        st.error(f"eggNOG 本地注释异常：{e}")
        debug(f"exception: {e}")
        return None

    finally:
        if fasta_path and os.path.exists(fasta_path):
            try:
                os.unlink(fasta_path)
            except Exception:
                pass

        if output_dir and os.path.exists(output_dir):
            try:
                shutil.rmtree(output_dir, ignore_errors=True)
            except Exception:
                pass


def parse_eggnog_manual_file(file_content):
    """Parse manually uploaded eggNOG annotation TSV file."""
    debug("开始解析手动上传的 eggNOG 文件...")

    lines = file_content.split("\n")
    data_lines = [line for line in lines if not line.startswith("#") and line.strip()]

    if not data_lines:
        st.error("文件中未找到有效数据行。")
        return None

    try:
        df = pd.read_csv(StringIO("\n".join(data_lines)), sep="\t")
        if df.columns[0].startswith("#"):
            df.rename(columns={df.columns[0]: "Master protein IDs"}, inplace=True)

        debug(f"解析成功：{df.shape[0]} 行，{df.shape[1]} 列")
        return df

    except Exception as e:
        st.error(f"解析 TSV 文件失败：{e}")
        debug(f"parse error: {e}")
        return None