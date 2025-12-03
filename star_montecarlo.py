import subprocess
from pathlib import Path
import os
from concurrent.futures import ThreadPoolExecutor, as_completed
import json
import argparse
import random

# ============================================================
# Constants
# ============================================================

STAR_EXEC = "STAR"
THREADS = 5
STAR_THREADS = 1
GENOME_DIR = "index"

# Parameter space
PARAMS = {
    "--scoreGap": [-10, -5, 0, 5, 10],
    "--scoreGapNoncan": [-16, -12, -8, -4, 0],
    "--scoreGapGCAG": [-8, -6, -4, -2, 0],
    "--scoreGapATAC": [-16, -12, -8, -4, 0],
    "--scoreDelOpen": [-4, -3, -2, -1, 0],
    "--scoreInsOpen": [-4, -3, -2, -1, 0],
    "--scoreStitchSJshift": [-4, -2, 1, 3, 5],
    "--seedSearchStartLmax": [25, 37.5, 50, 62.5, 75],
    "--seedSearchStartLmaxOverLread": [0.5, 0.75, 1.0, 1.25, 1.5],
    "--seedPerReadNmax": [500, 750, 1000, 1250, 1500],
    "--seedPerWindowNmax": [25, 37.5, 50, 62.5, 75],
    "--seedSplitMin": [6, 9, 12, 15, 18],
    "--winAnchorDistNbins": [4, 6, 9, 11, 13],
    "--outFilterScoreMin": [-5, -2.5, 0, 2.5, 5]
}

# ============================================================
# Helpers
# ============================================================

def sample_random_config(seed=None):
    if seed is not None:
        random.seed(seed)
    return {param: random.choice(vals) for param, vals in PARAMS.items()}


def sanitize(v: str) -> str:
    return str(v).replace(",", "-").replace("=", "-").replace(".", "_")


def make_dirname(cfg: dict) -> str:
    parts = []
    for k, v in cfg.items():
        parts.append(f"{k.strip('-')}{sanitize(v)}")
    return "_".join(parts)


def run_star(run_id, config, args):
    # Output directory name composed of parameters
    dir_name = make_dirname(config)
    out_dir = Path(args.output_root) / dir_name
    out_dir.mkdir(parents=True, exist_ok=True)

    # temp prefix required by STAR
    tmp_prefix = str(out_dir / "tmp_")

    # final output paths
    sam_out = out_dir / "output.sam"
    log_file = out_dir / "log.txt"
    config_file = out_dir / "config.txt"

    # base STAR command
    cmd = [
        STAR_EXEC,
        "--runThreadN", str(args.star_threads),
        "--genomeDir", GENOME_DIR,
        "--readFilesIn", args.input,
        "--outFileNamePrefix", tmp_prefix,
        "--outSAMtype", "SAM",
        "--outSAMunmapped", "Within",
        "--outSAMattributes", "NH", "HI", "AS", "nM", "NM", "MD", "jM", "jI", "XS", "MC",
    ]

    # add sampled parameters
    for p, v in config.items():
        cmd.extend([p, str(v)])

    # optional MARGI mode
    if args.margi == 1:
        cmd.extend([
            "--chimOutType", "Junctions",
            "--chimSegmentMin", "5"
        ])

    # write config file
    with open(config_file, "w") as cf:
        for p, v in config.items():
            cf.write(f"{p} {v}\n")

    # run STAR
    with open(log_file, "w") as log:
        subprocess.run(cmd, stdout=log, stderr=log)

    # identify STAR SAM (it ends with "Aligned.out.sam")
    generated_sam = Path(tmp_prefix + "Aligned.out.sam")
    if generated_sam.exists():
        generated_sam.rename(sam_out)


# ============================================================
# Main
# ============================================================

def main():
    global GENOME_DIR, THREADS, STAR_THREADS

    parser = argparse.ArgumentParser(
        description="Run Monte Carlo STAR alignments with random parameters."
    )

    parser.add_argument("--input", required=True, help="Input FASTQ")
    parser.add_argument("--part", default="rna")
    parser.add_argument("--seed", type=int)
    parser.add_argument("--config-log", default="generated_parameters.json")
    parser.add_argument("--index", help="STAR genome index directory")
    parser.add_argument("--threads", type=int)
    parser.add_argument("--star-threads", type=int)
    parser.add_argument("--runs", type=int, default=100)
    parser.add_argument("--margi", type=int, default=0)
    parser.add_argument("--output-root", type=str, default=".",
                        help="Directory to place per-run subdirectories")

    args = parser.parse_args()

    # apply options
    if args.index:
        GENOME_DIR = args.index
    if args.threads:
        THREADS = args.threads
    if args.star_threads:
        if args.star_threads < 1:
            parser.error("--star-threads must be >=1")
        STAR_THREADS = args.star_threads
    args.star_threads = STAR_THREADS

    # random seed
    if args.seed:
        random.seed(args.seed)

    # generate configs
    configs = [sample_random_config((args.seed or 1) + i)
               for i in range(args.runs)]

    # save config log
    with open(args.config_log, "w") as jf:
        json.dump(configs, jf, indent=4)

    # parallel execution
    with ThreadPoolExecutor(max_workers=THREADS) as exe:
        futures = [
            exe.submit(run_star, i, cfg, args)
            for i, cfg in enumerate(configs)
        ]

        for i, fut in enumerate(as_completed(futures), 1):
            try:
                fut.result()
                print(f"[{i}/{args.runs}] run completed")
            except Exception as e:
                print(f"[{i}/{args.runs}] run failed: {e}")


if __name__ == "__main__":
    main()

