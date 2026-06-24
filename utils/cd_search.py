# utils/cd_search.py
# 完整版：自动分批提交 + 本地缓存 + 断点续跑 + 自动合并
# 支持从缓存加载已完成的结果
# 支持在新电脑上重新分批运行

import requests
import time
import tempfile
import os
import re
import glob
import json
from io import StringIO
import streamlit as st
import pandas as pd

DEBUG = True

def debug(msg: str):
    if DEBUG:
        timestamp = time.strftime("%H:%M:%S")
        print(f"[{timestamp}][CD-SEARCH] {msg}", flush=True)

# ==================== 配置 ====================
CD_ENDPOINT = "https://www.ncbi.nlm.nih.gov/Structure/bwrpsb/bwrpsb.cgi"
CD_CHUNK_SIZE = 950  # NCBI 限制每次不超过 1000 条
CD_MAX_WAIT = 900    # 每个批次最多等待 900 秒
CD_POLL_INTERVAL = 20  # 轮询间隔 20 秒

# ==================== 辅助函数 ====================

def get_cache_root():
    """获取 cd_cache 根目录的绝对路径"""
    current_dir = os.path.dirname(os.path.abspath(__file__))
    cache_root = os.path.join(current_dir, "..", "cd_cache")
    os.makedirs(cache_root, exist_ok=True)
    return cache_root


def make_task_id(fasta_text: str) -> str:
    """
    根据 FASTA 内容生成任务 ID。
    使用 SHA256 哈希，确保相同内容生成相同 ID。
    """
    import hashlib
    hash_val = hashlib.sha256(fasta_text.encode('utf-8')).hexdigest()[:16]
    return f"task_{hash_val}"


def parse_fasta_for_cd(fasta_text: str) -> tuple:
    """
    解析 FASTA 文本，返回 (ids, sequences) 列表。
    """
    debug("解析 FASTA 文本...")
    ids = []
    sequences = []
    current_id = None
    current_seq = []
    
    for line in fasta_text.splitlines():
        line = line.strip()
        if not line:
            continue
        if line.startswith(">"):
            if current_id is not None:
                ids.append(current_id)
                sequences.append("".join(current_seq))
            current_id = line[1:].split()[0]
            current_seq = []
        else:
            if current_id is not None:
                current_seq.append(line)
    
    if current_id is not None:
        ids.append(current_id)
        sequences.append("".join(current_seq))
    
    debug(f"解析完成：共 {len(ids)} 条序列")
    return ids, sequences


def split_into_chunks(ids: list, sequences: list, chunk_size: int = CD_CHUNK_SIZE) -> list:
    """
    将序列列表分割成多个批次。
    返回: [(chunk_ids, chunk_sequences, chunk_index), ...]
    """
    n_total = len(ids)
    n_chunks = (n_total + chunk_size - 1) // chunk_size
    chunks = []
    
    for i in range(n_chunks):
        start = i * chunk_size
        end = min(start + chunk_size, n_total)
        chunk_ids = ids[start:end]
        chunk_seqs = sequences[start:end]
        chunks.append((chunk_ids, chunk_seqs, i + 1))
    
    debug(f"分割完成：{n_total} 条序列 → {n_chunks} 个批次（每批最多 {chunk_size} 条）")
    return chunks


def write_fasta(ids: list, sequences: list, filepath: str):
    """写入 FASTA 文件"""
    with open(filepath, 'w', encoding='utf-8') as f:
        for seq_id, seq in zip(ids, sequences):
            f.write(f">{seq_id}\n{seq}\n")
    debug(f"写入 FASTA: {filepath} ({len(ids)} 条序列)")


def submit_cd_search(fasta_path: str) -> str:
    """
    提交 FASTA 文件到 NCBI CD-Search，返回 cdsid。
    失败返回 None。
    """
    debug(f"提交 CD-Search: {fasta_path}")
    
    try:
        with open(fasta_path, 'rb') as fp:
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
            debug(f"提交失败，HTTP {resp.status_code}")
            return None
        
        match = re.search(r'QM3-qcdsearch-[A-Za-z0-9-]+', resp.text)
        if not match:
            debug(f"未找到 CD-Search ID，响应前200字符: {resp.text[:200]}")
            return None
        
        cdsid = match.group(0)
        debug(f"提交成功，cdsid = {cdsid}")
        return cdsid
        
    except Exception as e:
        debug(f"提交异常: {e}")
        return None


def poll_and_fetch(cdsid: str, max_wait: int = CD_MAX_WAIT, interval: int = CD_POLL_INTERVAL) -> str:
    """
    轮询 CD-Search 结果，返回 TSV 文本内容。
    失败返回 None。
    """
    debug(f"开始轮询 {cdsid}，最多等待 {max_wait}s...")
    elapsed = 0
    
    while elapsed < max_wait:
        time.sleep(interval)
        elapsed += interval
        
        try:
            resp = requests.post(CD_ENDPOINT, data={'cdsid': cdsid, 'tdata': 'hits'}, timeout=120)
        except Exception as e:
            debug(f"轮询请求异常: {e}")
            continue
        
        status_lines = [line for line in resp.text.split('\n') if '#status' in line]
        if status_lines:
            try:
                status = int(re.search(r'(\d+)', status_lines[0]).group(1))
            except:
                status = -1
            
            if status == 0:
                debug(f"搜索完成（{elapsed}s）")
                return resp.text
            elif status == 3:
                debug(f"搜索进行中...（{elapsed}s）")
                continue
            else:
                debug(f"未知状态码: {status}（{elapsed}s）")
        
        # 检查是否有错误
        if 'error' in resp.text.lower() or 'too many' in resp.text.lower():
            debug(f"检测到错误响应: {resp.text[:300]}")
            return None
    
    debug(f"轮询超时（{max_wait}s）")
    return None


def parse_cd_response(content: str) -> pd.DataFrame:
    """
    解析 CD-Search 返回的 TSV 内容为 DataFrame。
    """
    lines = content.split('\n')
    data_lines = [l for l in lines if not l.startswith('#') and l.strip()]
    
    if not data_lines:
        debug("响应中无有效数据行")
        return None
    
    try:
        df = pd.read_csv(StringIO('\n'.join(data_lines)), sep='\t', dtype=str)
        # 清洗第一列（Query ID）
        id_col = df.columns[0]
        df[id_col] = df[id_col].str.replace(r'Q#\d+\s*-\s*>?', '', regex=True).str.strip()
        debug(f"解析成功：{len(df)} 条记录")
        return df
    except Exception as e:
        debug(f"解析 TSV 失败: {e}")
        return None


# ==================== 核心功能 ====================

def batch_cd_search_with_cache(fasta_text, progress_callback=None):
    """
    批量 CD-Search：自动分批 → 提交 → 缓存 → 合并。
    
    参数:
        fasta_text: FASTA 格式的序列文本
        progress_callback: 可选的回调函数，用于更新 Streamlit 进度
    
    返回:
        (DataFrame, task_dir): 合并后的结果和任务目录路径
    """
    debug("=" * 60)
    debug("开始批量 CD-Search（分批 + 缓存模式）")
    debug("=" * 60)
    
    # 解析 FASTA
    ids, sequences = parse_fasta_for_cd(fasta_text)
    n_total = len(ids)
    
    if n_total == 0:
        st.error("FASTA 中没有找到任何序列")
        return None, None
    
    # 创建任务目录
    task_id = make_task_id(fasta_text)
    cache_root = get_cache_root()
    task_dir = os.path.join(cache_root, task_id)
    os.makedirs(task_dir, exist_ok=True)
    debug(f"任务目录: {task_dir}")
    
    # 分割批次
    chunks = split_into_chunks(ids, sequences, CD_CHUNK_SIZE)
    n_chunks = len(chunks)
    
    # 进度文件
    progress_file = os.path.join(task_dir, "progress.txt")
    
    def log_progress(msg: str):
        """写入进度文件并打印到终端"""
        timestamp = time.strftime("%Y-%m-%d %H:%M:%S")
        line = f"[{timestamp}] {msg}"
        with open(progress_file, 'a', encoding='utf-8') as f:
            f.write(line + '\n')
        debug(msg)
    
    log_progress(f"===== START: {n_total} 条序列 / {n_chunks} 个批次 =====")
    
    # 处理每个批次
    collected = []
    success_count = 0
    skip_count = 0
    fail_count = 0
    failed_chunks = []
    
    for chunk_idx, (chunk_ids, chunk_seqs, chunk_num) in enumerate(chunks):
        chunk_index = chunk_num  # 从1开始
        debug(f"\n--- 处理批次 {chunk_index}/{n_chunks} ({len(chunk_ids)} 条序列) ---")
        
        # 更新进度
        if progress_callback:
            progress_callback(chunk_index / n_chunks)
        
        # 检查缓存
        cache_file = os.path.join(task_dir, f"chunk_{chunk_index:04d}.tsv")
        
        if os.path.exists(cache_file):
            debug(f"批次 {chunk_index}: 发现缓存文件，直接加载")
            try:
                df = pd.read_csv(cache_file, sep='\t', dtype=str)
                collected.append(df)
                skip_count += 1
                log_progress(f"Chunk {chunk_index}/{n_chunks}: SKIP (cached, {len(df)} rows)")
                continue
            except Exception as e:
                debug(f"缓存文件损坏: {e}，将重新搜索")
                os.remove(cache_file)
        
        # 写入临时 FASTA
        temp_fasta = os.path.join(task_dir, f"chunk_{chunk_index:04d}.fasta")
        write_fasta(chunk_ids, chunk_seqs, temp_fasta)
        
        # 提交搜索
        cdsid = submit_cd_search(temp_fasta)
        
        # 删除临时 FASTA（节省空间）
        try:
            os.remove(temp_fasta)
        except:
            pass
        
        if cdsid is None:
            debug(f"批次 {chunk_index}: 提交失败")
            fail_count += 1
            failed_chunks.append(chunk_index)
            log_progress(f"Chunk {chunk_index}/{n_chunks}: FAIL (submit)")
            time.sleep(30)  # 失败后等待30秒再继续
            continue
        
        # 轮询结果
        content = poll_and_fetch(cdsid)
        
        if content is None:
            debug(f"批次 {chunk_index}: 获取结果失败")
            fail_count += 1
            failed_chunks.append(chunk_index)
            log_progress(f"Chunk {chunk_index}/{n_chunks}: FAIL (fetch)")
            time.sleep(30)
            continue
        
        # 解析结果
        df = parse_cd_response(content)
        
        if df is None or df.empty:
            debug(f"批次 {chunk_index}: 解析结果为空")
            fail_count += 1
            failed_chunks.append(chunk_index)
            log_progress(f"Chunk {chunk_index}/{n_chunks}: FAIL (parse)")
            time.sleep(30)
            continue
        
        # 保存缓存
        df.to_csv(cache_file, sep='\t', index=False)
        collected.append(df)
        success_count += 1
        log_progress(f"Chunk {chunk_index}/{n_chunks}: DONE ({len(df)} rows)")
        debug(f"缓存已保存: {cache_file}")
        
        # 批次之间休息 10 秒，避免请求过快
        time.sleep(10)
    
    # 最终汇总
    debug(f"\n===== 汇总 =====")
    debug(f"成功: {success_count}/{n_chunks}")
    debug(f"跳过(缓存): {skip_count}/{n_chunks}")
    debug(f"失败: {fail_count}/{n_chunks}")
    
    if not collected:
        log_progress("===== END: 所有批次均失败 =====")
        return None, task_dir
    
    # 合并所有结果
    debug("合并所有批次结果...")
    combined = pd.concat(collected, ignore_index=True)
    
    # 保存合并结果
    combined_file = os.path.join(task_dir, "cd_search_all.tsv")
    combined.to_csv(combined_file, sep='\t', index=False)
    debug(f"合并结果已保存: {combined_file}")
    debug(f"合并后共 {len(combined)} 条记录，{combined[combined.columns[0]].nunique()} 个唯一蛋白")
    
    log_progress(f"===== END: {len(combined)} rows, {success_count} success, {skip_count} skip, {fail_count} fail =====")
    
    if progress_callback:
        progress_callback(1.0)
    
    return combined, task_dir


def batch_cd_search_simple(fasta_text):
    """
    简化版：提交单个 FASTA 到 NCBI CD-Search（用于少量序列）。
    返回 DataFrame 或 None。
    """
    debug("开始自动 CD-Search（单次提交）...")
    
    if not fasta_text.strip():
        st.warning("FASTA 内容为空")
        return None
    
    with tempfile.NamedTemporaryFile(mode='w', suffix='.fasta', delete=False, encoding='utf-8') as f:
        f.write(fasta_text)
        temp_path = f.name
    
    debug(f"临时 FASTA: {temp_path}")
    
    try:
        cdsid = submit_cd_search(temp_path)
        if cdsid is None:
            return None
        
        content = poll_and_fetch(cdsid)
        if content is None:
            return None
        
        df = parse_cd_response(content)
        return df
        
    except Exception as e:
        debug(f"异常: {e}")
        return None
    finally:
        try:
            if os.path.exists(temp_path):
                os.unlink(temp_path)
        except:
            pass


# ==================== 手动文件解析 ====================

def parse_cd_tsv(content):
    """
    解析手动上传的 CD-Search TSV 文本。
    """
    debug("开始解析手动上传的 CD-Search TSV 文件")
    lines = content.split('\n')
    data_lines = [l for l in lines if not l.startswith('#') and l.strip()]
    if not data_lines:
        st.warning("文件中未找到有效数据行")
        return None

    try:
        sep = '\t' if '\t' in data_lines[0] else ','
        df = pd.read_csv(StringIO('\n'.join(data_lines)), sep=sep, dtype=str)
        debug(f"解析到 {df.shape[0]} 行, {df.shape[1]} 列")
        id_col = df.columns[0]
        df[id_col] = df[id_col].str.replace(r'Q#\d+\s*-\s*>?', '', regex=True).str.strip()
        debug(f"清洗后共 {len(df)} 条记录")
        return df
    except Exception as e:
        st.error(f"解析 TSV 时出错: {e}")
        return None


# ==================== 缓存加载 ====================

def load_cd_cache(task_name: str = None):
    """
    从本地 cd_cache/ 目录加载已完成的 CD-Search 结果。
    """
    cache_root = get_cache_root()
    debug(f"扫描缓存目录: {cache_root}")
    
    if not os.path.exists(cache_root):
        debug("缓存目录不存在")
        return None, None
    
    subdirs = [d for d in os.listdir(cache_root) if os.path.isdir(os.path.join(cache_root, d))]
    debug(f"发现 {len(subdirs)} 个子目录: {subdirs}")
    
    if not subdirs:
        return None, None
    
    if task_name is not None:
        if task_name not in subdirs:
            debug(f"任务 '{task_name}' 不在缓存中")
            return None, None
        subdirs = [task_name]
    
    all_chunks = []
    loaded = 0
    failed = 0
    
    for subdir in subdirs:
        subdir_path = os.path.join(cache_root, subdir)
        
        # 优先使用合并文件
        combined_file = os.path.join(subdir_path, "cd_search_all.tsv")
        if os.path.exists(combined_file):
            try:
                df = pd.read_csv(combined_file, sep='\t', dtype=str)
                debug(f"加载合并文件: {combined_file} ({len(df)} rows)")
                all_chunks.append(df)
                loaded += 1
                continue
            except Exception as e:
                debug(f"读取合并文件失败: {e}")
        
        # 回退：加载各个分块文件
        chunk_files = sorted(glob.glob(os.path.join(subdir_path, "chunk_*.tsv")))
        debug(f"目录 {subdir}: 找到 {len(chunk_files)} 个分块文件")
        
        for cf in chunk_files:
            try:
                df = pd.read_csv(cf, sep='\t', dtype=str)
                all_chunks.append(df)
                loaded += 1
            except Exception as e:
                debug(f"读取 {cf} 失败: {e}")
                failed += 1
    
    debug(f"加载完成: 成功 {loaded}, 失败 {failed}")
    
    if not all_chunks:
        return None, None
    
    combined = pd.concat(all_chunks, ignore_index=True)
    debug(f"合并后共 {len(combined)} 条记录")
    
    task_used = task_name if task_name else subdirs[0]
    return combined, task_used


def get_cached_tasks():
    """获取所有可用的缓存任务名称"""
    cache_root = get_cache_root()
    if not os.path.exists(cache_root):
        return []
    tasks = [d for d in os.listdir(cache_root) if os.path.isdir(os.path.join(cache_root, d))]
    debug(f"可用缓存任务: {tasks}")
    return tasks
