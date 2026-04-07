#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse, numpy as np
import matplotlib.pyplot as plt

def main():
    ap = argparse.ArgumentParser("Plot TRS null (per RBP, per cell type)")
    ap.add_argument("--npz", required=True, help="*.ct.nulls_perRBP.npz")
    ap.add_argument("--rbp", required=True, help="RBP name, e.g. CELF1")
    ap.add_argument("--cell-type", required=True, help="cell type exactly as in results, e.g. Monocytes")
    ap.add_argument("--bins", type=int, default=40)
    ap.add_argument("--out", required=True, help="output PDF")  # ← 仅修改了这行的描述
    args = ap.parse_args()

    data = np.load(args.npz, allow_pickle=True)
    rbps       = data["rbps"].tolist()
    cell_types = data["cell_types"].tolist()
    K_vec      = data["K_null"]

    trs_real_mm = data["trs_real_mm"]            # (CT × R)
    trs_null_list = data["trs_null_list"]        # len(R), each (CT × K_j)

    try:
        j = rbps.index(args.rbp)
    except ValueError:
        raise SystemExit(f"RBP '{args.rbp}' not found. Choices: {rbps[:5]}... total={len(rbps)}")

    try:
        i = cell_types.index(args.cell_type)
    except ValueError:
        raise SystemExit(f"cell_type '{args.cell_type}' not found. Choices: {cell_types}")

    K = int(K_vec[j])
    null_ct = trs_null_list[j][i, :]            # (K,)
    real_ct = float(trs_real_mm[i, j])

    # empirical p（与主程序一致的右尾）
    p_emp = ( (null_ct > real_ct).sum() + 1.0 ) / (K + 1.0)

    plt.figure(figsize=(6.4,4.8))
    if K > 0:
        plt.hist(null_ct, bins=args.bins, edgecolor="black",color="#BDD2FF")
        y0,y1 = plt.ylim()
        plt.vlines(real_ct, y0, y1, linestyles="dashed")
        plt.ylim(y0,y1)
        title = f"{args.rbp} | {args.cell_type}  (K={K})  REAL={real_ct:.3g}  p={p_emp:.3g}"
    else:
        plt.axvline(real_ct, linestyle="dashed")
        title = f"{args.rbp} | {args.cell_type}  (K=0)  REAL={real_ct:.3g}"

    plt.title(title)
    plt.xlabel("TRS (min–max by per-RBP null)")
    plt.ylabel("Frequency")
    plt.tight_layout()
    plt.savefig(args.out, format="pdf")  # ← 仅修改了这行：保存为PDF矢量图
    print(f"[OK] wrote {args.out}")

if __name__ == "__main__":
    main()

