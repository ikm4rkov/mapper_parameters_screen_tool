#!/usr/bin/env python3
"""
Flexible report generator for UU/UM and uniquely-mapped RNA/DNA counts.

Usage:
  python3 make_report_flexible.py \
      -r /path/to/rawdir \
      -i /path/to/filtereddir \
      -d /path/to/outdir \
      -o outprefix

Use --quiet to suppress debug output.
"""
import os
import argparse
import traceback
from collections import Counter
import pandas as pd
import matplotlib.pyplot as plt

# Known pairtype tokens
KNOWN_PAIRTYPES = {"UU", "UM", "MU", "MM", "UN", "NU", "MN", "NM", "NN"}

def dbg(msg, quiet=False):
    if not quiet:
        print(msg)

# ---------------- file discovery helpers ----------------
def possible_names_for_filtered(base):
    """
    Return list of plausible filtered filenames given base (without suffix)
    The pipeline tends to use:
     <base>_Unique_RNAfiltered_<base>_Unique_RNA.tab.rc
    but we are tolerant about missing extensions etc.
    """
    names = []
    # canonical
    names.append(f"{base}_Unique_RNAfiltered_{base}_Unique_RNA.tab.rc")
    names.append(f"{base}_Unique_RNAfiltered_{base}_Unique_RNA")  # maybe without ext
    # sometimes repeated without the second "_Unique_RNA"
    names.append(f"{base}_Unique_RNAfiltered_{base}.tab.rc")
    names.append(f"{base}_Unique_RNAfiltered_{base}")
    # older style: maybe just base + "filtered" variations
    names.append(f"{base}_Unique_RNAfiltered.tab.rc")
    names.append(f"{base}_Unique_RNAfiltered")
    return names

def possible_names_for_other(base):
    names = [
        f"{base}_Other.tab.rc",
        f"{base}_Other",
        f"{base}Other.tab.rc",
        f"{base}Other"
    ]
    return names

def find_existing(pathdir, names):
    """Return first existing full path for one of names (list), or None."""
    for n in names:
        fp = os.path.join(pathdir, n)
        if os.path.exists(fp):
            return fp
    return None

# ---------------- pairtype parsing/counting ----------------
def read_pairtype_series(path, quiet=False):
    """
    Try to read the pairtype series from a TSV file.
    Strategy:
      1) read header (first line) and if contains 'ATA_pairtype' or 'pairtype' use that column
      2) else read with pandas (header=0) and look for those column names
      3) else read without header and take column index 1 (second column)
    Returns a pandas Series of stripped strings (may be empty).
    """
    if not os.path.exists(path):
        dbg(f"[read_pairtype_series] file not found: {path}", quiet)
        return pd.Series([], dtype=str)

    # Try to read header line
    try:
        with open(path, "r", encoding="utf-8") as fh:
            first_line = fh.readline().rstrip("\n")
    except Exception as e:
        dbg(f"[read_pairtype_series] cannot read file header {path}: {e}", quiet)
        return pd.Series([], dtype=str)

    header_tokens = first_line.split("\t")
    # If header explicitly contains column name
    for candidate in ("ATA_pairtype", "pairtype"):
        if candidate in header_tokens:
            try:
                df = pd.read_csv(path, sep="\t", dtype=str, usecols=[candidate])
                s = df[candidate].astype(str).str.strip()
                dbg(f"[read_pairtype_series] used header column '{candidate}' in {path} (n={len(s)})", quiet)
                return s
            except Exception as e:
                dbg(f"[read_pairtype_series] failed reading header column '{candidate}' in {path}: {e}", quiet)
                break

    # Try reading as header=0 and see if any column name contains pairtype substring
    try:
        df0 = pd.read_csv(path, sep="\t", dtype=str, nrows=1)
        cols = list(df0.columns)
        for c in cols:
            if "pairtype" in c.lower() or c.lower().startswith("ata_pairtype"):
                try:
                    df = pd.read_csv(path, sep="\t", dtype=str, usecols=[c])
                    s = df[c].astype(str).str.strip()
                    dbg(f"[read_pairtype_series] used column '{c}' (header heuristics) in {path} (n={len(s)})", quiet)
                    return s
                except Exception:
                    pass
    except Exception:
        # may not have header or pandas fails - fall through
        pass

    # Last fallback: read without header and take the second column (index 1)
    try:
        df_nohdr = pd.read_csv(path, sep="\t", dtype=str, header=None, usecols=[1])
        s = df_nohdr.iloc[:, 0].astype(str).str.strip()
        dbg(f"[read_pairtype_series] fallback to column index 1 in {path} (n={len(s)})", quiet)
        return s
    except Exception as e:
        dbg(f"[read_pairtype_series] final fallback failed for {path}: {e}", quiet)
        return pd.Series([], dtype=str)

def count_pairtypes_from_series(s):
    """Return Counter of pairtype occurrences from pandas Series s"""
    cnt = Counter()
    if s is None:
        return cnt
    try:
        # Filter only known two-letter codes
        s2 = s.dropna().astype(str).str.strip()
        # If there are verbose tags (e.g., 'UU;something') strip to first token
        s2 = s2.str.split(";").str[0]
        # Keep only tokens of length 2
        s2 = s2[s2.str.len() == 2]
        for v in s2:
            if v in KNOWN_PAIRTYPES:
                cnt[v] += 1
    except Exception:
        pass
    return cnt

def get_pairtype_counts(path, quiet=False):
    """Return Counter for known pairtypes in file path (quiet controls debug)."""
    s = read_pairtype_series(path, quiet=quiet)
    cnt = count_pairtypes_from_series(s)
    dbg(f"[get_pairtype_counts] {os.path.basename(path)} -> counts: {dict(cnt)}", quiet)
    return cnt

# ---------------- computing UU/UM and U_m_R / U_m_D ----------------
def get_uu_um_from_counts(cnt):
    return int(cnt.get("UU", 0)), int(cnt.get("UM", 0))

def compute_unique_R_D(unique_cnt, other_cnt):
    """
    Now using updated formula:
      U_m_R = UU + UM + UN
      U_m_D = UU + MU + NU
    """
    uu = unique_cnt.get("UU", 0)
    um = unique_cnt.get("UM", 0)
    un = other_cnt.get("UN", 0)
    mu = other_cnt.get("MU", 0)
    nu = other_cnt.get("NU", 0)

    U_m_R = int(uu + um + un)
    U_m_D = int(uu + mu + nu)
    return U_m_R, U_m_D

# ---------------- main --------------------------------------------
def main():
    p = argparse.ArgumentParser(description="Flexible UU/UM and unique RNA/DNA report (tolerant discovery)")
    p.add_argument("-r", "--rawdir", required=True, help="Raw dir (contains Unique_RNA and Other files)")
    p.add_argument("-i", "--filtereddir", required=True, help="Filtered dir (contains filtered Unique_RNA files)")
    p.add_argument("-d", "--outdir", required=True, help="Output directory")
    p.add_argument("-o", "--outprefix", required=True, help="Output prefix")
    p.add_argument("--quiet", action="store_true", help="Suppress debug output")
    args = p.parse_args()
    quiet = args.quiet

    dbg("=== ARGUMENTS ===", quiet)
    dbg(f"rawdir    : {args.rawdir}", quiet)
    dbg(f"filtereddir: {args.filtereddir}", quiet)
    dbg(f"outdir    : {args.outdir}", quiet)
    dbg(f"outprefix : {args.outprefix}", quiet)
    dbg("=================", quiet)

    if not os.path.isdir(args.rawdir):
        print(f"[ERROR] rawdir not found: {args.rawdir}")
        return
    if not os.path.isdir(args.filtereddir):
        print(f"[ERROR] filtereddir not found: {args.filtereddir}")
        return
    os.makedirs(args.outdir, exist_ok=True)

    raw_entries = sorted(os.listdir(args.rawdir))
    filt_entries = sorted(os.listdir(args.filtereddir))
    dbg(f"Found {len(raw_entries)} entries in rawdir (sample): {raw_entries[:40]}", quiet)
    #dbg(f"Found {len(filt_entries)} entries in filtereddir (sample): {filt_entries[:40]}", quiet)

    summary = []
    raw_points = []
    filt_points = []
    umapped_points = []

    # Candidate raw unique files: contains "Unique_RNA" and not 'filtered'/'cigar_stat'/'out' etc.
    candidates = []
    for fn in raw_entries:
        name_l = fn.lower()
        if "unique_rna" in name_l and ("filtered" not in name_l) and ("cigar_stat" not in name_l) and fn.endswith(".tab.rc"):
            candidates.append(fn)
    # if none found with .tab.rc, relax to any filename containing Unique_RNA
    if not candidates:
        for fn in raw_entries:
            name_l = fn.lower()
            if "unique_rna" in name_l and ("filtered" not in name_l) and ("cigar_stat" not in name_l):
                candidates.append(fn)
    dbg(f"Candidate raw unique files (count={len(candidates)}): {candidates[:40]}", quiet)

    for fname in candidates:
        dbg(f"\n--- candidate: {fname}", quiet)
        # Attempt to get base: try to remove trailing suffix "_Unique_RNA.tab.rc" or "Unique_RNA"
        name_noext = fname
        if name_noext.endswith(".tab.rc"):
            name_noext = name_noext[:-7]
        # find "Unique_RNA" occurrence
        idx = name_noext.find("Unique_RNA")
        if idx == -1:
            # try lowercase
            idx = name_noext.lower().find("unique_rna")
        if idx == -1:
            dbg(f" can't determine base for {fname}, skipping", quiet)
            continue
        base = name_noext[:idx].rstrip("_")
        dbg(f" inferred base='{base}'", quiet)

        # Paths to try
        raw_unique_path = os.path.join(args.rawdir, fname)
        other_candidates = possible_names_for_other(base)
        raw_other_path = find_existing(args.rawdir, other_candidates)
        # If not found in rawdir try some variants (maybe 'Other' in filtered dir)
        if raw_other_path is None:
            raw_other_path = find_existing(args.filtereddir, other_candidates)

        filtered_candidates = possible_names_for_filtered(base)
        filtered_path = find_existing(args.filtereddir, filtered_candidates)
        # sometimes filtered file is in raw dir (rare), try there too
        if filtered_path is None:
            filtered_path = find_existing(args.rawdir, filtered_candidates)

        dbg(f" expected raw_unique: {raw_unique_path} -> exists={os.path.exists(raw_unique_path)}", quiet)
        dbg(f" expected raw_other  (tried): {other_candidates} -> found={raw_other_path}", quiet)
        dbg(f" expected filtered   (tried): {filtered_candidates} -> found={filtered_path}", quiet)

        missing = []
        if not os.path.exists(raw_unique_path):
            missing.append(raw_unique_path)
        if raw_other_path is None or not os.path.exists(raw_other_path):
            missing.append("(other) " + (raw_other_path or "NOT_FOUND"))
        if filtered_path is None or not os.path.exists(filtered_path):
            missing.append("(filtered) " + (filtered_path or "NOT_FOUND"))

        if missing:
            print(f"[WARN] Skipping base='{base}' because missing files:")
            for m in missing:
                print("  -", m)
            continue

        # Count pairtypes
        try:
            unique_cnt = get_pairtype_counts(raw_unique_path, quiet=quiet)
            other_cnt = get_pairtype_counts(raw_other_path, quiet=quiet)
            filt_cnt = get_pairtype_counts(filtered_path, quiet=quiet)
        except Exception as e:
            print(f"[ERROR] counting pairtypes for base='{base}': {e}")
            print(traceback.format_exc())
            continue

        uu_raw = int(unique_cnt.get("UU", 0)); um_raw = int(unique_cnt.get("UM", 0))
        uu_filt = int(filt_cnt.get("UU", 0)); um_filt = int(filt_cnt.get("UM", 0))
        U_m_R, U_m_D = compute_unique_R_D(unique_cnt, other_cnt)

        dbg(f"[RESULT] base={base}: UU_raw={uu_raw}, UM_raw={um_raw}, UU_filt={uu_filt}, UM_filt={um_filt}, U_m_R={U_m_R}, U_m_D={U_m_D}", quiet)

        raw_points.append((um_raw, uu_raw))
        filt_points.append((um_filt, uu_filt))
        umapped_points.append((U_m_D, U_m_R))

        summary.append({
            "base": base,
            "raw_unique_file": os.path.basename(raw_unique_path),
            "raw_other_file": os.path.basename(raw_other_path),
            "filtered_file": os.path.basename(filtered_path),
            "UU_raw": uu_raw,
            "UM_raw": um_raw,
            "UU_filtered": uu_filt,
            "UM_filtered": um_filt,
            "U_m_R": U_m_R,
            "U_m_D": U_m_D
        })

    # write summary
    out_tsv = os.path.join(args.outdir, f"{args.outprefix}.tsv")
    try:
        df = pd.DataFrame(summary)
        df.to_csv(out_tsv, sep="\t", index=False)
        print(f"[OK] Saved summary: {out_tsv} (rows={len(df)})")
    except Exception as e:
        print(f"[ERROR] writing summary TSV: {e}")
        print(traceback.format_exc())

    # Plot UU/UM raw vs filtered
    out_png_uuum = os.path.join(args.outdir, f"{args.outprefix}_UU_UM.png")
    try:
        fig, axes = plt.subplots(1, 2, figsize=(12, 6))
        fig.suptitle("UU vs UM Counts (Raw and Filtered)")
        ax = axes[0]
        for um, uu in raw_points:
            ax.scatter(um, uu)
        ax.set_xlabel("UM count"); ax.set_ylabel("UU count"); ax.set_title("Raw Unique_RNA")
        ax = axes[1]
        for um, uu in filt_points:
            ax.scatter(um, uu)
        ax.set_xlabel("UM count"); ax.set_ylabel("UU count"); ax.set_title("Filtered Unique_RNA")
        plt.tight_layout(rect=[0,0,1,0.96])
        plt.savefig(out_png_uuum, dpi=300)
        plt.close(fig)
        print(f"[OK] Saved plot: {out_png_uuum}")
    except Exception as e:
        print(f"[ERROR] making UU/UM plot: {e}")
        print(traceback.format_exc())

    # Plot U_m_R vs U_m_D
    out_png_um = os.path.join(args.outdir, f"{args.outprefix}_UniqueRNA_DNA.png")
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

    dbg(f"Finished. Processed {len(summary)} bases.", quiet)

if __name__ == "__main__":
    main()

