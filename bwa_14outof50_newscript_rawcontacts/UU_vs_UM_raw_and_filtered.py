#!/usr/bin/env python3
import os
import pandas as pd
import matplotlib.pyplot as plt

def classify_file(filename):
    """
    Classify file as 'raw', 'filtered', or 'skip' based on naming rules.
    """
    name = filename.lower()

    # Skip metadata, 'Other', and filtered_out files
    if any(s in name for s in ["filteredcigar_stat", "filtered_out", "other"]):
        return "skip"

    # Filtered data
    if "filteredfiltered" in name:
        return "filtered"

    # Raw data
    if "dna_output" in name and "filtered" not in name:
        return "raw"

    return "skip"  # default safety

def get_color(filename):
    """Return bright color based on 'Other' or 'Unique' in filename."""
    name = filename.lower()
    if "unique" in name:
        return "blue"
    else:
        return "orange"  # bright fallback

def count_pairtypes(filepath):
    """Count UU and UM in ATA_pairtype or pairtype column."""
    try:
        # Read header to find column name
        with open(filepath, "r", encoding="utf-8") as f:
            header_line = f.readline().strip()
        headers = header_line.split("\t")

        col_name = None
        for candidate in ["ATA_pairtype", "pairtype"]:
            if candidate in headers:
                col_name = candidate
                break

        if col_name is None:
            print(f"Warning: No valid pairtype column found in {filepath}")
            return 0, 0

        df = pd.read_csv(filepath, sep="\t", usecols=[col_name], dtype=str)
        counts = df[col_name].value_counts()
        uu_count = counts.get("UU", 0)
        um_count = counts.get("UM", 0)
        return uu_count, um_count

    except Exception as e:
        print(f"Skipping {filepath}: {e}")
        return 0, 0

def main():
    files = [f for f in os.listdir(".") if f.endswith(".tab.rc")]
    if not files:
        print("No .tab.rc files found in the current directory.")
        return

    data = {"raw": [], "filtered": []}
    summary_records = []

    for file in files:
        file_type = classify_file(file)
        if file_type == "skip":
            continue

        uu_count, um_count = count_pairtypes(file)
        color = get_color(file)
        data[file_type].append((file, uu_count, um_count, color))
        summary_records.append({
            "filename": file,
            "file_type": file_type,
            "UU_count": uu_count,
            "UM_count": um_count,
            "color": color
        })

    # Sort summary by UU count (descending)
    summary_df = pd.DataFrame(summary_records)
    if not summary_df.empty:
        summary_df.sort_values(by="UU_count", ascending=False, inplace=True)
        summary_df.to_csv("UU_vs_UM_summary.tsv", sep="\t", index=False)
        print("Saved summary to: UU_vs_UM_summary.tsv")
    else:
        print("No valid data found to summarize.")

    # Create a figure with two subplots (shared axes)
    fig, axes = plt.subplots(1, 2, figsize=(12, 6), sharex=True, sharey=True)
    fig.suptitle("UU vs UM Counts (Raw and Filtered Contacts)")

    # Plot RAW contacts
    ax = axes[0]
    for _, uu, um, color in data["raw"]:
        ax.scatter(um, uu, color=color)
    ax.set_xlabel("UM count")
    ax.set_ylabel("UU count")
    ax.set_title("Raw Contacts")

    # Plot FILTERED contacts
    ax = axes[1]
    for _, uu, um, color in data["filtered"]:
        ax.scatter(um, uu, color=color)
    ax.set_xlabel("UM count")
    ax.set_ylabel("UU count")
    ax.set_title("Filtered Contacts")

    plt.tight_layout(rect=[0, 0, 1, 0.96])
    output_file = "UU_vs_UM_raw_and_filtered.png"
    plt.savefig(output_file, dpi=300)
    print(f"Saved plot as: {output_file}")

if __name__ == "__main__":
    main()

