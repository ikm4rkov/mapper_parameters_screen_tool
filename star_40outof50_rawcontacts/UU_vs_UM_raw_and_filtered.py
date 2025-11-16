#!/usr/bin/env python3
import os
import pandas as pd
import matplotlib.pyplot as plt

def classify_file(filepath, root_dir):
    """
    Classify files as 'raw', 'filtered', or 'skip' based on location and name.
    """
    name = os.path.basename(filepath).lower()
    parent = os.path.basename(os.path.dirname(filepath)).lower()

    # Skip irrelevant files
    if any(s in name for s in ["filteredcigar_stat", "filtered_out", "other"]):
        return "skip"

    # Filtered data: only in 'filtered/*Unique*/' subdirectories
    if parent.startswith("run_") and "filtered" in os.path.dirname(filepath):
        if "unique" in parent:
            return "filtered"
        else:
            return "skip"

    # Raw data: files directly in root directory
    if os.path.dirname(filepath) == root_dir:
        if "filtered" not in name and name.endswith(".tab.rc"):
            return "raw"

    return "skip"

def get_color(filename):
    """Return bright color based on 'Unique' or 'Other' in filename."""
    name = filename.lower()
    if "unique" in name:
        return "blue"
    return "orange"

def count_pairtypes(filepath):
    """Count UU and UM in ATA_pairtype or pairtype column."""
    try:
        with open(filepath, "r", encoding="utf-8") as f:
            header_line = f.readline().strip()
        headers = header_line.split("\t")

        col_name = next((c for c in ["ATA_pairtype", "pairtype"] if c in headers), None)
        if not col_name:
            print(f"⚠️  Warning: no valid pairtype column found in {filepath}")
            return 0, 0

        df = pd.read_csv(filepath, sep="\t", usecols=[col_name], dtype=str)
        counts = df[col_name].value_counts()
        uu_count = counts.get("UU", 0)
        um_count = counts.get("UM", 0)
        return uu_count, um_count
    except Exception as e:
        print(f"⚠️  Skipping {filepath}: {e}")
        return 0, 0

def main():
    root_dir = os.getcwd()
    filtered_dir = os.path.join(root_dir, "filtered")

    # Collect all relevant files
    files_to_process = []

    # Raw files: directly in root directory
    for f in os.listdir(root_dir):
        if f.endswith(".tab.rc"):
            files_to_process.append(os.path.join(root_dir, f))

    # Filtered files: search in filtered/ subdirectories
    if os.path.isdir(filtered_dir):
        for subdir, _, files in os.walk(filtered_dir):
            for f in files:
                if f.endswith(".tab.rc"):
                    files_to_process.append(os.path.join(subdir, f))

    data = {"raw": [], "filtered": []}
    summary_records = []

    for path in files_to_process:
        file_type = classify_file(path, root_dir)
        if file_type == "skip":
            continue

        uu_count, um_count = count_pairtypes(path)
        color = get_color(path)
        data[file_type].append((path, uu_count, um_count, color))
        summary_records.append({
            "filename": os.path.basename(path),
            "file_type": file_type,
            "UU_count": uu_count,
            "UM_count": um_count,
            "color": color
        })

    summary_df = pd.DataFrame(summary_records)
    if not summary_df.empty:
        summary_df.sort_values(by="UU_count", ascending=False, inplace=True)
        out_tsv = os.path.join(root_dir, "UU_vs_UM_summary.tsv")
        summary_df.to_csv(out_tsv, sep="\t", index=False)
        print(f"✅ Saved summary table to: {out_tsv}")
    else:
        print("⚠️  No valid data found to summarize.")
        return

    # Create plot
    fig, axes = plt.subplots(1, 2, figsize=(12, 6), sharex=True, sharey=True)
    fig.suptitle("UU vs UM Counts (Raw and Filtered Contacts)")

    # RAW contacts
    ax = axes[0]
    for _, uu, um, color in data["raw"]:
        ax.scatter(um, uu, color=color)
    ax.set_xlabel("UM count")
    ax.set_ylabel("UU count")
    ax.set_title("Raw Contacts")

    # FILTERED contacts
    ax = axes[1]
    for _, uu, um, color in data["filtered"]:
        ax.scatter(um, uu, color=color)
    ax.set_xlabel("UM count")
    ax.set_ylabel("UU count")
    ax.set_title("Filtered Contacts")

    plt.tight_layout(rect=[0, 0, 1, 0.96])
    out_png = os.path.join(root_dir, "UU_vs_UM_raw_and_filtered.png")
    plt.savefig(out_png, dpi=300)
    print(f"✅ Saved plot as: {out_png}")

if __name__ == "__main__":
    main()

