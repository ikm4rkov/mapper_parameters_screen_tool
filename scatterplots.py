#!/usr/bin/env python3
"""
Revised flexible report generator for UU/UM and uniquely-mapped RNA/DNA counts.

New layout assumptions:
- Raw contacts stored in subdirectories:
    <rawdir>/rna_<rna_mapper>_<rna_run>__dna_<dna_mapper>_<dna_run>/

- Raw files inside each subdir:
    raw_contacts_Unique_RNA.tab.rc
    raw_contacts_Other.tab.rc

- Filtered contacts stored in corresponding subdirectories:
    <filtereddir>/rna_<rna_mapper>_<rna_run>__dna_<dna_mapper>_<dna_run>/

- Filtered file name (fixed):
    filtered_raw_contacts_Unique_RNA.tab.rc

Usage:
python3 scatterplots.py \
    -r raw_contacts \
    -i filtered_contacts \
    -o outprefix \
    -d /path/to/outdir \
    --rna star \
    --dna bwa

RNA/DNA mapper values are lowercased automatically and used for directory filtering.
"""

import os
import argparse
import traceback
import re
from collections import Counter
import pandas as pd
import matplotlib.pyplot as plt

# ---------------- constants ----------------
KNOWN_PAIRTYPES = {"UU", "UM", "MU", "MM", "UN", "NU", "MN", "NM", "NN"}

RAW_UNIQUE_NAME = "raw_contacts_Unique_RNA.tab.rc"
RAW_OTHER_NAME = "raw_contacts_Other.tab.rc"
FILT_UNIQUE_NAME = "filtered_raw_contacts_Unique_RNA.tab.rc"

# ---------------- utils ----------------
def dbg(msg, quiet=False):
    if not quiet:
        print(msg)

# ---------------- pairtype parsing/counting ----------------
def read_pairtype_series(path, quiet=False):
    if not os.path.exists(path):
        dbg(f"[read_pairtype_series] file not found: {path}", quiet)
        return pd.Series([], dtype=str)

    try:
        with open(path, "r", encoding="utf-8") as fh:
            first_line = fh.readline().rstrip("\n")
    except Exception as e:
        dbg(f"[read_pairtype_series] cannot read file header {path}: {e}", quiet)
        return pd.Series([], dtype=str)

    header_tokens = first_line.split("\t")

    # direct header match
    for candidate in ("ATA_pairtype", "pairtype"):
        if candidate in header_tokens:
            try:
                df = pd.read_csv(path, sep="\t", dtype=str, usecols=[candidate])
                s = df[candidate].astype(str).str.strip()
                dbg(f"[read_pairtype_series] used header column '{candidate}' in {path} (n={len(s)})", quiet)
                return s
            except Exception as e:
                dbg(f"[read_pairtype_series] failed reading header column '{candidate}' in {path}: {e}", quiet)

    # heuristic header scan
    try:
        df0 = pd.read_csv(path, sep="\t", dtype=str, nrows=1)
        for c in df0.columns:
            cl = c.lower()
            if "pairtype" in cl or cl.startswith("ata_pairtype"):
                try:
                    df = pd.read_csv(path, sep="\t", dtype=str, usecols=[c])
                    s = df[c].astype(str).str.strip()
                    dbg(f"[read_pairtype_series] used column '{c}' (heuristic) in {path} (n={len(s)})", quiet)
                    return s
                except Exception:
                    pass
    except Exception:
        pass

    # fallback: 2nd column
    try:
        df_nohdr = pd.read_csv(path, sep="\t", dtype=str, header=None, usecols=[1])
        s = df_nohdr.iloc[:, 0].astype(str).str.strip()
        dbg(f"[read_pairtype_series] fallback to column index 1 in {path} (n={len(s)})", quiet)
        return s
    except Exception as e:
        dbg(f"[read_pairtype_series] final fallback failed for {path}: {e}", quiet)
        return pd.Series([], dtype=str)


def count_pairtypes_from_series(s):
    cnt = Counter()
    if s is None:
        return cnt
    try:
        s2 = s.dropna().astype(str).str.strip()
        s2 = s2.str.split(";").str[0]
        s2 = s2[s2.str.len() == 2]
        for v in s2:
            if v in KNOWN_PAIRTYPES:
                cnt[v] += 1
    except Exception:
        pass
    return cnt


def get_pairtype_counts(path, quiet=False):
    s = read_pairtype_series(path, quiet=quiet)
    cnt = count_pairtypes_from_series(s)
    dbg(f"[get_pairtype_counts] {os.path.basename(path)} -> {dict(cnt)}", quiet)
    return cnt

# ---------------- computing metrics ----------------
def compute_unique_R_D(unique_cnt, other_cnt):
    # Updated formula
    # U_m_R = UU + UM + UN
    # U_m_D = UU + MU + NU
    uu = unique_cnt.get("UU", 0)
    um = unique_cnt.get("UM", 0)
    un = other_cnt.get("UN", 0)
    mu = other_cnt.get("MU", 0)
    nu = other_cnt.get("NU", 0)
    U_m_R = int(uu + um + un)
    U_m_D = int(uu + mu + nu)
    return U_m_R, U_m_D

# ---------------- main ----------------
def main():
    p = argparse.ArgumentParser(description="UU/UM and uniquely-mapped RNA/DNA report (mapper-based directory layout)")

    # directories
    p.add_argument("-r", "--rawdir", default="raw_contacts", help="Raw contacts root directory (default: raw_contacts)")
    p.add_argument("-i", "--filtereddir", default="filtered_contacts", help="Filtered contacts root directory (default: filtered_contacts)")
    p.add_argument("-d", "--outdir", required=True, help="Output directory")
    p.add_argument("-o", "--outprefix", required=True, help="Output prefix")

    # mappers
    p.add_argument("--rna", required=True, help="RNA mapper name (e.g. star)")
    p.add_argument("--dna", required=True, help="DNA mapper name (e.g. bwa)")

    p.add_argument("--quiet", action="store_true", help="Suppress debug output")

    args = p.parse_args()
    quiet = args.quiet

    raw_root = args.rawdir
    filt_root = args.filtereddir
    outdir = args.outdir

    rna_mapper = args.rna.lower()
    dna_mapper = args.dna.lower()

    if not os.path.isdir(raw_root):
        print(f"[ERROR] rawdir not found: {raw_root}")
        return
    if not os.path.isdir(filt_root):
        print(f"[ERROR] filtereddir not found: {filt_root}")
        return

    os.makedirs(outdir, exist_ok=True)

    # regex for subdir matching
    # rna_<rna>_<num>__dna_<dna>_<num>
    pat = re.compile(rf"^rna_{rna_mapper}_[^_]+__dna_{dna_mapper}_[^_]+$", re.IGNORECASE)

    raw_subdirs = [d for d in os.listdir(raw_root) if os.path.isdir(os.path.join(raw_root, d)) and pat.match(d)]

    dbg(f"Matched subdirectories (count={len(raw_subdirs)}): {raw_subdirs}", quiet)

    summary = []
    raw_points = []
    filt_points = []
    umapped_points = []

    for sub in sorted(raw_subdirs):
        raw_dir = os.path.join(raw_root, sub)
        filt_dir = os.path.join(filt_root, sub)

        raw_unique = os.path.join(raw_dir, RAW_UNIQUE_NAME)
        raw_other = os.path.join(raw_dir, RAW_OTHER_NAME)
        filt_unique = os.path.join(filt_dir, FILT_UNIQUE_NAME)

        dbg(f"\n--- processing {sub}", quiet)
        dbg(f" raw_unique: {raw_unique}", quiet)
        dbg(f" raw_other : {raw_other}", quiet)
        dbg(f" filt_unique: {filt_unique}", quiet)

        missing = []
        if not os.path.exists(raw_unique): missing.append(raw_unique)
        if not os.path.exists(raw_other): missing.append(raw_other)
        if not os.path.exists(filt_unique): missing.append(filt_unique)

        if missing:
            print(f"[WARN] Skipping {sub} (missing files):")
            for m in missing:
                print(" -", m)
            continue

        try:
            unique_cnt = get_pairtype_counts(raw_unique, quiet=quiet)
            other_cnt = get_pairtype_counts(raw_other, quiet=quiet)
            filt_cnt = get_pairtype_counts(filt_unique, quiet=quiet)
        except Exception as e:
            print(f"[ERROR] counting pairtypes in {sub}: {e}")
            print(traceback.format_exc())
            continue

        uu_raw = int(unique_cnt.get("UU", 0)); um_raw = int(unique_cnt.get("UM", 0))
        uu_filt = int(filt_cnt.get("UU", 0)); um_filt = int(filt_cnt.get("UM", 0))
        U_m_R, U_m_D = compute_unique_R_D(unique_cnt, other_cnt)

        raw_points.append((um_raw, uu_raw))
        filt_points.append((um_filt, uu_filt))
        umapped_points.append((U_m_D, U_m_R))

        summary.append({
            "dataset": sub,
            "UU_raw": uu_raw,
            "UM_raw": um_raw,
            "UU_filtered": uu_filt,
            "UM_filtered": um_filt,
            "U_m_R": U_m_R,
            "U_m_D": U_m_D
        })

    # -------- outputs --------
    out_tsv = os.path.join(outdir, f"{args.outprefix}.tsv")
    try:
        df = pd.DataFrame(summary)
        df.to_csv(out_tsv, sep="\t", index=False)
        print(f"[OK] Saved summary: {out_tsv} (rows={len(df)})")
    except Exception as e:
        print(f"[ERROR] writing summary TSV: {e}")
        print(traceback.format_exc())

    # UU/UM plots
    out_png_uuum = os.path.join(outdir, f"{args.outprefix}_UU_UM.png")
    try:
        fig, axes = plt.subplots(1, 2, figsize=(12, 6))
        fig.suptitle("UU vs UM Counts (Raw and Filtered)")

        ax = axes[0]
        for um, uu in raw_points:
            ax.scatter(um, uu)
        ax.set_xlabel("UM count")
        ax.set_ylabel("UU count")
        ax.set_title("Raw Unique_RNA")

        ax = axes[1]
        for um, uu in filt_points:
            ax.scatter(um, uu)
        ax.set_xlabel("UM count")
        ax.set_ylabel("UU count")
        ax.set_title("Filtered Unique_RNA")

        plt.tight_layout(rect=[0,0,1,0.96])
        plt.savefig(out_png_uuum, dpi=300)
        plt.close(fig)
        print(f"[OK] Saved plot: {out_png_uuum}")
    except Exception as e:
        print(f"[ERROR] making UU/UM plot: {e}")
        print(traceback.format_exc())

    # Unique RNA vs DNA
    out_png_um = os.path.join(outdir, f"{args.outprefix}_UniqueRNA_DNA.png")
    try:
        fig2, ax2 = plt.subplots(figsize=(6,6))
        for dna, rna in umapped_points:
            ax2.scatter(dna, rna)
        ax2.set_xlabel("U_m_D (uniquely mapped DNA)")
        ax2.set_ylabel("U_m_R (uniquely mapped RNA)")
        ax2.set_title("Uniquely mapped RNA vs DNA")
        plt.tight_layout()
        plt.savefig(out_png_um, dpi=300)
        plt.close(fig2)
        print(f"[OK] Saved plot: {out_png_um}")
    except Exception as e:
        print(f"[ERROR] making UniqueRNA_DNA plot: {e}")
        print(traceback.format_exc())

    dbg(f"Finished. Processed {len(summary)} datasets.", quiet)


if __name__ == "__main__":
    main()

