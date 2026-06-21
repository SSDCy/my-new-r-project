import requests
import time
import tempfile
import os
import re
from io import StringIO
import streamlit as st
import pandas as pd

CD_ENDPOINT = "https://www.ncbi.nlm.nih.gov/Structure/bwrpsb/bwrpsb.cgi"

def batch_cd_search(fasta_text):
    """自动发送 FASTA 到 NCBI CD-Search，返回注释 DataFrame"""
    if not fasta_text.strip():
        st.warning("FASTA 内容为空")
        return None

    with tempfile.NamedTemporaryFile(mode='w', suffix='.fasta', delete=False) as f:
        f.write(fasta_text)
        temp_path = f.name

    try:
        with open(temp_path, 'rb') as fp:
            params = {
                'db': 'cdd',
                'smode': 'auto',
                'tdata': 'hits',
                'dmode': 'rep',
                'cddefl': 'false',
                'qdefl': 'false',
                'clonly': 'false',
                'useid1': 'true',
            }
            files = {'queries': fp}
            resp = requests.post(CD_ENDPOINT, data=params, files=files, timeout=120)
        if resp.status_code != 200:
            st.error(f"CD-Search 提交失败，HTTP {resp.status_code}")
            return None

        match = re.search(r'QM3-qcdsearch-[A-Za-z0-9-]+', resp.text)
        if not match:
            st.error("未在响应中找到 CD-Search ID")
            return None
        cdsid = match.group(0)

        max_wait = 900
        elapsed = 0
        interval = 20
        while elapsed < max_wait:
            time.sleep(interval)
            elapsed += interval
            status_resp = requests.post(CD_ENDPOINT, data={'cdsid': cdsid, 'tdata': 'hits'}, timeout=120)
            status_lines = [line for line in status_resp.text.split('\n') if '#status' in line]
            if status_lines:
                status = int(re.search(r'(\d+)', status_lines[0]).group(1))
                if status == 0:
                    break
                elif status == 3:
                    continue
            # 其他情况继续等待

        lines = status_resp.text.split('\n')
        data_lines = [l for l in lines if not l.startswith('#') and l.strip()]
        if not data_lines:
            st.warning("CD-Search 返回空结果")
            return None

        df = pd.read_csv(StringIO('\n'.join(data_lines)), sep='\t')
        id_col = df.columns[0]
        df[id_col] = df[id_col].str.replace(r'Q#\d+ -\s*>?', '', regex=True)
        return df
    except Exception as e:
        st.error(f"CD-Search 出错: {e}")
        return None
    finally:
        try:
            if os.path.exists(temp_path):
                os.unlink(temp_path)
        except Exception as e:
            print(f"[WARNING] 无法删除临时文件 {temp_path}: {e}")


def parse_cd_tsv(content):
    """
    解析手动上传的 CD-Search TSV 文本。
    格式说明：
        - 文件可能包含以 '#' 开头的注释行。
        - 数据部分为 Tab 分隔，第一列通常是 Query ID (可能带有 'Q#1 - >' 前缀)。
    返回清洗后的 DataFrame，若失败则返回 None。
    """
    print("[DEBUG] parse_cd_tsv: 开始解析手动上传的 CD-Search 文件")
    lines = content.split('\n')
    # 过滤掉注释行和空行
    data_lines = [l for l in lines if not l.startswith('#') and l.strip()]
    if not data_lines:
        st.warning("文件中未找到有效数据行")
        return None

    try:
        # 检测分隔符（可能是 Tab 或逗号）
        sep = '\t' if '\t' in data_lines[0] else ','
        df = pd.read_csv(StringIO('\n'.join(data_lines)), sep=sep, dtype=str)
        print(f"[DEBUG] 解析到 {df.shape[0]} 行, {df.shape[1]} 列")
        # 清洗第一列（Query ID）：移除 "Q#1 - >" 等前缀
        id_col = df.columns[0]
        df[id_col] = df[id_col].str.replace(r'Q#\d+\s*-\s*>?', '', regex=True).str.strip()
        return df
    except Exception as e:
        st.error(f"解析 TSV 时出错: {e}")
        return None
