import os
import json
import re
import shutil

def main():
    files = [f for f in os.listdir(".") if os.path.isfile(f)]

    archive_dir = "original_filenames"
    os.makedirs(archive_dir, exist_ok=True)

    # Revised regex:
    # - dna_output followed by ANYTHING until __rna_output (1+ underscores)
    # - rna_output followed by ANYTHING until __output_
    pattern = re.compile(
        r"^(dna_output.*?__rna_output_+.*?)__output_(Unique_RNA|Other)\.tab\.rc$"
    )

    groups = {}
    mapping = {}

    # ---- Grouping ----
    for f in files:
        m = pattern.match(f)
        if m:
            signature = m.group(1)  # the dna→rna parameter block
            groups.setdefault(signature, []).append(f)

    # ---- Renaming + Moving ----
    for i, (signature, fnames) in enumerate(groups.items(), start=1):
        mapping[str(i)] = fnames

        for f in fnames:
            suffix = "Unique_RNA.tab.rc" if "Unique_RNA" in f else "Other.tab.rc"
            new_name = f"dna_output_{i}_rna_output_{i}_{suffix}"

            print(f"\nProcessing: {f}")
            print(f"  → New name: {new_name}")

            # move original to archive
            original_path = f
            archived_path = os.path.join(archive_dir, f)
            shutil.move(original_path, archived_path)

            # create renamed copy
            shutil.copy(archived_path, new_name)

    # ---- Save JSON ----
    with open("name_map.json", "w") as fw:
        json.dump(mapping, fw, indent=4)

    print("\nDONE.")
    print("Original long filenames moved to:", archive_dir)
    print("Mapping stored in name_map.json")

if __name__ == "__main__":
    main()

