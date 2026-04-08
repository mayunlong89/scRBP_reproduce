"""
ENCODE Bulk RNA-seq Download Script — Fixed Version
====================================================
Fixes for IMR90 / PC3 / MCF7 returning 0 files:

Problem 1: biosample_ontology.term_name 对部分细胞系名称匹配失败
           (ENCODE内部用的是标准化名称，如 "PC-3" 而非 "PC3")
Problem 2: assay_title 精确匹配有时找不到，需要用 assay_slims
Problem 3: 部分实验的文件在 /files/ endpoint 而不嵌在实验里
Problem 4: 需要同时搜索 "polyA plus RNA-seq" 和 "total RNA-seq"
"""

import requests
import os
import json
import time
import sys
from pathlib import Path

BASE_URL = "https://www.encodeproject.org"
HEADERS  = {"Accept": "application/json"}
OUT_DIR  = Path("encode_data")

# ──────────────────────────────────────────────────
# 关键修复 1: ENCODE内部标准名称映射
#   运行 diagnose_cell_line() 可自动发现正确名称
# ──────────────────────────────────────────────────
CELL_LINE_ALIASES = {
    # 你用的名字      →   ENCODE term_name 候选列表（按优先级）
    "GM12878": ["GM12878"],
    "K562":    ["K562"],
    "HCT116":  ["HCT116"],
    "HepG2":   ["HepG2"],
    "MCF7":    ["MCF7"],                    # 有时写作 "MCF-7"
    "IMR90":   ["IMR-90", "IMR90"],         # ENCODE用 "IMR-90"
    "PC3":     ["PC-3", "PC3"],             # ENCODE用 "PC-3"
    "Panc1":   ["Panc1", "PANC-1", "Panc-1"],
}

# ──────────────────────────────────────────────────
# 关键修复 2: 多种assay_title组合
# ──────────────────────────────────────────────────
ASSAY_TITLES = [
    "polyA plus RNA-seq",
    "polyA minus RNA-seq",   # 备用
    "total RNA-seq",          # 备用（某些细胞系只有total）
]


def api_get(url, params=None, retries=3, delay=1.0):
    """带重试的GET请求"""
    for attempt in range(retries):
        try:
            r = requests.get(url, params=params, headers=HEADERS, timeout=30)
            r.raise_for_status()
            return r.json()
        except Exception as e:
            print(f"  [retry {attempt+1}/{retries}] {e}")
            time.sleep(delay * (attempt + 1))
    return None


# ══════════════════════════════════════════════════
# STEP 0: 诊断工具 — 找出ENCODE内部的正确名称
# ══════════════════════════════════════════════════
def diagnose_cell_line(query_name):
    """
    搜索ENCODE中包含该名称的所有biosample term，
    帮助找到正确的内部名称。
    """
    print(f"\n🔍 Diagnosing '{query_name}'...")
    url = f"{BASE_URL}/search/"
    params = {
        "type": "Biosample",
        "term_name": query_name,
        "format": "json",
        "limit": 10,
    }
    data = api_get(url, params)
    if not data:
        print("  API error")
        return []

    found = set()
    for item in data.get("@graph", []):
        name = item.get("term_name", "")
        if name:
            found.add(name)

    # 也搜索RNA-seq实验
    params2 = {
        "type": "Experiment",
        "assay_title": "polyA plus RNA-seq",
        "biosample_ontology.term_name": query_name,
        "format": "json",
        "limit": 5,
    }
    data2 = api_get(url, params2)
    if data2:
        for exp in data2.get("@graph", []):
            bs = exp.get("biosample_summary", "")
            term = exp.get("biosample_ontology", {})
            if isinstance(term, list):
                for t in term:
                    found.add(t.get("term_name", ""))
            elif isinstance(term, dict):
                found.add(term.get("term_name", ""))

    found.discard("")
    print(f"  Found ENCODE names: {sorted(found)}")
    return sorted(found)


def diagnose_all():
    """运行所有细胞系的诊断，打印正确名称映射"""
    print("=" * 60)
    print("ENCODE Cell Line Name Diagnostic")
    print("=" * 60)
    mapping = {}
    for cl in CELL_LINE_ALIASES:
        names = diagnose_cell_line(cl)
        mapping[cl] = names
    print("\n\n📋 Recommended CELL_LINE_ALIASES update:")
    for cl, names in mapping.items():
        print(f'  "{cl}": {names},')
    return mapping


# ══════════════════════════════════════════════════
# STEP 1: 搜索实验
# ══════════════════════════════════════════════════
def search_experiments(term_name, assay_title="polyA plus RNA-seq"):
    """搜索特定细胞系+assay的实验列表"""
    url = f"{BASE_URL}/search/"
    params = {
        "type": "Experiment",
        "assay_title": assay_title,
        "biosample_ontology.term_name": term_name,
        "status": "released",
        "format": "json",
        "limit": "all",
    }
    data = api_get(url, params)
    if not data:
        return []
    experiments = data.get("@graph", [])
    return experiments


def find_experiments_for_cell_line(cell_line_key):
    """
    尝试所有别名 × 所有assay_title，找到实验列表。
    返回 (experiments, matched_term_name, matched_assay)
    """
    aliases = CELL_LINE_ALIASES.get(cell_line_key, [cell_line_key])

    for assay in ASSAY_TITLES:
        for alias in aliases:
            exps = search_experiments(alias, assay)
            if exps:
                print(f"  ✓ Found {len(exps)} experiment(s) for "
                      f"'{alias}' / '{assay}'")
                return exps, alias, assay
            else:
                print(f"  ✗ No results: '{alias}' / '{assay}'")

    return [], None, None


# ══════════════════════════════════════════════════
# STEP 2: 从实验中提取文件
# ══════════════════════════════════════════════════
def get_gene_quant_files_from_experiment(exp_accession):
    """
    获取一个实验的gene quantification TSV文件列表。
    直接查询 /experiments/{accession}/ endpoint。
    """
    url = f"{BASE_URL}/experiments/{exp_accession}/?format=json"
    data = api_get(url)
    if not data:
        return []

    tsv_files = []
    for f in data.get("files", []):
        if (f.get("output_type") == "gene quantifications"
                and f.get("file_format") == "tsv"
                and f.get("status") == "released"
                and f.get("assembly") in ["GRCh38", None, ""]
                and "hg19" not in f.get("assembly", "")):
            tsv_files.append({
                "accession": f["accession"],
                "href": f.get("href", f"/files/{f['accession']}/@@download/{f['accession']}.tsv"),
                "replicate": f.get("biological_replicates", []),
                "assembly": f.get("assembly", "unknown"),
                "genome_annotation": f.get("genome_annotation", ""),
            })
    return tsv_files


def get_all_files_for_cell_line(cell_line_key):
    """获取某细胞系所有gene quantification文件的下载信息"""
    print(f"\n{'='*50}")
    print(f"Processing: {cell_line_key}")
    print(f"{'='*50}")

    experiments, term_name, assay = find_experiments_for_cell_line(cell_line_key)

    if not experiments:
        print(f"  ⚠️  No experiments found for {cell_line_key}")
        print(f"      → Try running diagnose_cell_line('{cell_line_key}')")
        return []

    all_files = []
    for exp in experiments:
        acc = exp.get("accession", "")
        if not acc:
            continue
        files = get_gene_quant_files_from_experiment(acc)
        print(f"  Exp {acc}: {len(files)} gene_quant TSV file(s)")
        for f in files:
            f["experiment"] = acc
            f["cell_line"]  = cell_line_key
            f["term_name"]  = term_name
            f["assay"]      = assay
        all_files.extend(files)

    return all_files


# ══════════════════════════════════════════════════
# STEP 3: 下载文件
# ══════════════════════════════════════════════════
def download_file(file_info, out_dir=OUT_DIR, skip_existing=True):
    """下载单个TSV文件"""
    cell_line = file_info["cell_line"]
    accession = file_info["accession"]
    href      = file_info["href"]

    save_dir  = out_dir / cell_line
    save_dir.mkdir(parents=True, exist_ok=True)
    save_path = save_dir / f"{accession}.tsv"

    if skip_existing and save_path.exists() and save_path.stat().st_size > 1000:
        print(f"  [skip] {accession}.tsv already exists")
        return True

    # 构建下载URL
    if href.startswith("http"):
        url = href
    else:
        url = f"{BASE_URL}{href}"

    try:
        r = requests.get(url, headers=HEADERS, timeout=60, stream=True)
        r.raise_for_status()
        with open(save_path, "wb") as fh:
            for chunk in r.iter_content(chunk_size=8192):
                fh.write(chunk)
        size_kb = save_path.stat().st_size // 1024
        print(f"  [ok] {accession}.tsv  ({size_kb} KB)"
              f"  rep={file_info['replicate']}"
              f"  asm={file_info['assembly']}")
        return True
    except Exception as e:
        print(f"  [error] {accession}: {e}")
        return False


# ══════════════════════════════════════════════════
# STEP 4: 备用方案 — 直接用文件搜索endpoint
# ══════════════════════════════════════════════════
def search_files_directly(term_name, assay_title="polyA plus RNA-seq"):
    """
    直接搜索 /files/ endpoint。
    某些细胞系的文件通过实验endpoint找不到，但能直接搜到。
    """
    url = f"{BASE_URL}/search/"
    params = {
        "type": "File",
        "output_type": "gene quantifications",
        "file_format": "tsv",
        "assembly": "GRCh38",
        "status": "released",
        "dataset.assay_title": assay_title,
        "dataset.biosample_ontology.term_name": term_name,
        "format": "json",
        "limit": "all",
    }
    data = api_get(url, params)
    if not data:
        return []
    return data.get("@graph", [])


def get_all_files_fallback(cell_line_key):
    """
    备用方案：直接搜索File对象。
    当实验搜索找不到文件时使用。
    """
    print(f"\n  [fallback] Searching files directly for {cell_line_key}...")
    aliases = CELL_LINE_ALIASES.get(cell_line_key, [cell_line_key])

    for assay in ASSAY_TITLES:
        for alias in aliases:
            files = search_files_directly(alias, assay)
            if files:
                print(f"  ✓ Fallback found {len(files)} file(s) for '{alias}'")
                result = []
                for f in files:
                    acc  = f.get("accession", "")
                    href = f.get("href", f"/files/{acc}/@@download/{acc}.tsv")
                    result.append({
                        "accession": acc,
                        "href":      href,
                        "replicate": f.get("biological_replicates", []),
                        "assembly":  f.get("assembly", "GRCh38"),
                        "genome_annotation": f.get("genome_annotation", ""),
                        "cell_line": cell_line_key,
                        "term_name": alias,
                        "assay":     assay,
                        "experiment": f.get("dataset", "").split("/")[-2],
                    })
                return result

    return []


# ══════════════════════════════════════════════════
# MAIN
# ══════════════════════════════════════════════════
def download_all_cell_lines(cell_lines=None, out_dir="encode_data",
                             run_diagnosis=False):
    """
    主入口函数。

    Parameters
    ----------
    cell_lines  : list of str, 默认下载所有8个
    out_dir     : 输出目录
    run_diagnosis : True = 先诊断名称再下载（推荐第一次运行时开启）
    """
    if cell_lines is None:
        cell_lines = list(CELL_LINE_ALIASES.keys())

    global OUT_DIR
    OUT_DIR = Path(out_dir)
    OUT_DIR.mkdir(exist_ok=True)

    # 可选：先运行诊断
    if run_diagnosis:
        diagnose_all()
        print("\n请根据上面的输出更新 CELL_LINE_ALIASES，然后重新运行（run_diagnosis=False）")
        return

    summary = {}

    for cl in cell_lines:
        # 方法1：通过实验搜索
        files = get_all_files_for_cell_line(cl)

        # 方法2：如果方法1失败，用文件直搜
        if not files:
            files = get_all_files_fallback(cl)

        if not files:
            print(f"  ❌ {cl}: 未找到任何文件，请手动检查ENCODE portal")
            summary[cl] = 0
            continue

        # 去重（同一个accession可能出现多次）
        seen = set()
        unique_files = []
        for f in files:
            if f["accession"] not in seen:
                seen.add(f["accession"])
                unique_files.append(f)

        print(f"\n  → Downloading {len(unique_files)} unique files for {cl}")
        ok = sum(download_file(f) for f in unique_files)
        summary[cl] = ok
        time.sleep(0.3)  # 避免请求过快

    # 打印汇总
    print("\n" + "=" * 50)
    print("Download Summary")
    print("=" * 50)
    total = 0
    for cl, n in summary.items():
        status = "✓" if n > 0 else "❌"
        print(f"  {status} {cl:12s}: {n} file(s)")
        total += n
    print(f"\n  Total: {total} files → {OUT_DIR}/")

    # 保存元数据
    meta_path = OUT_DIR / "download_summary.json"
    with open(meta_path, "w") as fh:
        json.dump(summary, fh, indent=2)
    print(f"  Metadata saved: {meta_path}")

    return summary


# ══════════════════════════════════════════════════
# 手动下载辅助：生成wget命令
# ══════════════════════════════════════════════════
def generate_wget_commands(cell_line_key, out_dir="encode_data"):
    """
    生成可直接粘贴到终端的wget命令。
    适合API下载失败时手动操作。
    """
    files = get_all_files_for_cell_line(cell_line_key)
    if not files:
        files = get_all_files_fallback(cell_line_key)

    if not files:
        print(f"No files found for {cell_line_key}")
        return

    save_dir = Path(out_dir) / cell_line_key
    print(f"\n# wget commands for {cell_line_key}")
    print(f"mkdir -p {save_dir}")
    for f in files:
        acc  = f["accession"]
        href = f["href"]
        url  = href if href.startswith("http") else f"{BASE_URL}{href}"
        print(f"wget -q -O {save_dir}/{acc}.tsv '{url}'")


# ══════════════════════════════════════════════════
# 验证已下载文件
# ══════════════════════════════════════════════════
def verify_downloads(data_dir="encode_data"):
    """检查每个细胞系下载了多少有效文件"""
    import glob
    data_dir = Path(data_dir)
    print("\n📁 Download Verification")
    print("=" * 50)
    for cl in CELL_LINE_ALIASES:
        tsvs = list((data_dir / cl).glob("*.tsv"))
        valid = [f for f in tsvs if f.stat().st_size > 10_000]
        empty = [f for f in tsvs if f.stat().st_size <= 10_000]
        status = "✓" if valid else "❌"
        print(f"  {status} {cl:12s}: {len(valid)} valid, {len(empty)} empty/small")
        if empty:
            for f in empty:
                print(f"       ⚠️  {f.name} ({f.stat().st_size} bytes) — re-download needed")


# ══════════════════════════════════════════════════
# 使用示例
# ══════════════════════════════════════════════════
if __name__ == "__main__":

    # ── 模式1：诊断名称（推荐第一次运行）──────────────
    # download_all_cell_lines(run_diagnosis=True)

    # ── 模式2：下载全部 ───────────────────────────────
    # download_all_cell_lines()

    # ── 模式3：只下载失败的细胞系 ────────────────────
    # download_all_cell_lines(cell_lines=["IMR90", "PC3", "MCF7"])

    # ── 模式4：验证已有下载 ───────────────────────────
    # verify_downloads("encode_data")

    # ── 模式5：生成wget命令（手动下载） ─────────────
    # generate_wget_commands("IMR90")
    # generate_wget_commands("PC3")
    # generate_wget_commands("MCF7")

    # 默认：诊断 + 下载失败的三个细胞系
    mode = sys.argv[1] if len(sys.argv) > 1 else "fix"

    if mode == "diagnose":
        diagnose_all()

    elif mode == "fix":
        # 只修复空的细胞系
        print("Fixing IMR90, PC3, MCF7 (and HCT116)...")
        download_all_cell_lines(cell_lines=["PC3","HCT116"])

    elif mode == "all":
        download_all_cell_lines()

    elif mode == "verify":
        verify_downloads()

    elif mode == "wget":
        for cl in ["IMR90", "PC3", "MCF7"]:
            generate_wget_commands(cl)
