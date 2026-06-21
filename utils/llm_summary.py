# utils/llm_summary.py
"""
AI自然语言功能摘要模块
- 支持调用 DeepSeek API 生成摘要
- 支持本地规则模板（无需 API Key）
"""

import os
import json
import requests
from typing import Optional, List, Dict

DEBUG = True
def debug(msg: str):
    if DEBUG:
        print(f"[LLM SUMMARY DEBUG] {msg}")

# ---------- 读取 API Key ----------
DEEPSEEK_API_KEY = os.environ.get("DEEPSEEK_API_KEY") or os.getenv("DEEPSEEK_API_KEY")
if not DEEPSEEK_API_KEY:
    debug("⚠️ DEEPSEEK_API_KEY 未设置，将使用本地规则模板")
else:
    debug(f"✅ DEEPSEEK_API_KEY 已加载，前缀: {DEEPSEEK_API_KEY[:8]}...")

# ---------- API 调用函数 ----------
def call_deepseek(prompt: str, max_tokens: int = 300, temperature: float = 0.3) -> Optional[str]:
    """
    调用 DeepSeek API 生成文本
    """
    if not DEEPSEEK_API_KEY:
        debug("API Key 为空，跳过 API 调用")
        return None

    url = "https://api.deepseek.com/v1/chat/completions"
    headers = {
        "Authorization": f"Bearer {DEEPSEEK_API_KEY}",
        "Content-Type": "application/json"
    }
    payload = {
        "model": "deepseek-chat",  # 或 deepseek-reasoner
        "messages": [
            {"role": "system", "content": "你是一位蛋白质功能注释专家，擅长用通俗易懂的中文解释蛋白质的功能。"},
            {"role": "user", "content": prompt}
        ],
        "max_tokens": max_tokens,
        "temperature": temperature,
        "stream": False
    }

    debug(f"调用 DeepSeek API，提示词长度: {len(prompt)} 字符")
    try:
        response = requests.post(url, headers=headers, json=payload, timeout=30)
        if response.status_code == 200:
            result = response.json()
            content = result['choices'][0]['message']['content'].strip()
            debug(f"API 返回成功，内容长度: {len(content)} 字符")
            return content
        else:
            debug(f"API 错误，状态码: {response.status_code}，响应: {response.text[:200]}")
            return None
    except Exception as e:
        debug(f"API 调用异常: {str(e)}")
        return None

# ---------- 构建 Prompt ----------
def build_summary_prompt(
    protein_id: str,
    protein_name: str = "",
    similarity_score: float = None,
    similar_protein: str = "",
    go_terms: List[str] = None,
    ec_numbers: List[str] = None
) -> str:
    """
    根据蛋白信息构建 Prompt
    """
    go_terms = go_terms or []
    ec_numbers = ec_numbers or []

    # 构建信息块
    info_lines = [f"蛋白ID: {protein_id}"]
    if protein_name:
        info_lines.append(f"蛋白名称: {protein_name}")
    if similarity_score is not None and similar_protein:
        info_lines.append(f"与 {similar_protein} 的相似度: {similarity_score:.3f}")
    if go_terms:
        info_lines.append(f"预测的 GO 术语: {', '.join(go_terms)}")
    if ec_numbers:
        info_lines.append(f"预测的 EC 编号: {', '.join(ec_numbers)}")

    info_text = "\n".join(info_lines)

    prompt = f"""
请根据以下蛋白质信息，用 **1-2 句中文** 概括其可能的功能：

{info_text}

要求：
1. 只基于上面提供的信息进行推断，不要编造信息
2. 使用流畅、通俗的生物学术语
3. 如果信息不足，诚实说明“该蛋白功能尚不明确，但与……相似”

请直接输出描述，不要加额外格式。
"""
    return prompt

# ---------- 主函数：生成摘要 ----------
def generate_function_summary(
    protein_id: str,
    protein_name: str = "",
    similarity_score: float = None,
    similar_protein: str = "",
    go_terms: List[str] = None,
    ec_numbers: List[str] = None,
    use_api: bool = True
) -> str:
    """
    生成蛋白质功能摘要。
    - 如果 use_api=True 且 API Key 有效，调用 DeepSeek。
    - 否则使用本地规则模板。
    """
    debug(f"生成摘要: protein_id={protein_id}, use_api={use_api}")

    # 如果 go_terms 是字符串形式，转换为列表
    if isinstance(go_terms, str):
        go_terms = [g.strip() for g in go_terms.split(';') if g.strip()]
    if isinstance(ec_numbers, str):
        ec_numbers = [e.strip() for e in ec_numbers.split(';') if e.strip()]

    # 尝试用 API
    if use_api and DEEPSEEK_API_KEY:
        prompt = build_summary_prompt(
            protein_id, protein_name, similarity_score,
            similar_protein, go_terms, ec_numbers
        )
        result = call_deepseek(prompt)
        if result:
            return result
        else:
            debug("API 调用失败，回退到本地模板")

    # 备用方案：本地规则模板
    return _generate_template_summary(protein_id, protein_name, go_terms, ec_numbers)

# ---------- 备用方案：本地规则模板（无需 API）----------
# 一个简易的 GO → 中文描述映射表
GO_TEMPLATES = {
    "GO:0003674": "分子功能",
    "GO:0003824": "催化活性",
    "GO:0005215": "转运蛋白活性",
    "GO:0005488": "结合",
    "GO:0005515": "蛋白结合",
    "GO:0005524": "ATP 结合",
    "GO:0004672": "蛋白激酶活性",
    "GO:0004713": "蛋白酪氨酸激酶活性",
    "GO:0004683": "钙调蛋白依赖性蛋白激酶活性",
    "GO:0000166": "核苷酸结合",
    "GO:0005525": "GTP 结合",
    "GO:0016740": "转移酶活性",
    "GO:0016787": "水解酶活性",
    "GO:0016491": "氧化还原酶活性",
    "GO:0016301": "激酶活性",
    "GO:0016772": "转移酶活性，转移含磷基团",
    "GO:0032553": "核糖核苷酸结合",
    "GO:0046872": "金属离子结合",
    "GO:0043167": "离子结合",
    "GO:0003774": "运动活性",
    "GO:0005198": "结构分子活性",
    "GO:0003676": "核酸结合",
    "GO:0003700": "转录因子活性",
}

def _generate_template_summary(protein_id, protein_name, go_terms, ec_numbers):
    """
    基于规则的摘要生成（不联网、不花钱的备用方案）
    """
    debug("使用本地模板生成摘要")
    parts = []

    # 添加蛋白名称
    if protein_name:
        parts.append(f"{protein_name}")
    else:
        parts.append(f"蛋白 {protein_id}")

    # 提取功能描述
    func_descs = []
    for go in go_terms:
        if go in GO_TEMPLATES:
            func_descs.append(GO_TEMPLATES[go])
        elif go.startswith("GO:"):
            func_descs.append(go)  # 保留原始ID

    if ec_numbers:
        func_descs.append(f"EC 编号: {', '.join(ec_numbers)}")

    if func_descs:
        parts.append(f"可能具有以下功能: {', '.join(set(func_descs))}")
    else:
        parts.append("尚未有明确的功能注释信息")

    return "；".join(parts) + "。"
