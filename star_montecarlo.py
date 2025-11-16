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
THREADS = 5          # Number of concurrent STAR jobs
STAR_THREADS = 1      # Threads per STAR process (--runThreadN)
GENOME_DIR = "index"
OUT_ROOT = "star_montecarlo_outputs"

# Parameter space derived from common STAR presets + extensions
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
# Helper functions
# ============================================================

def sample_random_config(seed=None):
    """Sample a random STAR parameter configuration."""
    if seed is not None:
        random.seed(seed)
    return {param: random.choice(values) for param, values in PARAMS.items()}


def build_output_dir(run_id, part):
    """Create grouped output directory (10 runs per folder)."""
    group = run_id // 10 + 1
    dir_name = os.path.join(OUT_ROOT, f"{group}{part}")
    Path(dir_name).mkdir(parents=True, exist_ok=True)
    return dir_name


def run_star(run_id, config, args):
    """Run a single STAR alignment with given configuration."""
    out_dir = build_output_dir(run_id, args.part)
    prefix = f"run_{run_id:04d}_"
    out_prefix = os.path.join(out_dir, prefix)

    cmd = [
        STAR_EXEC,
        "--runThreadN", str(args.star_threads),
        "--genomeDir", GENOME_DIR,
        "--readFilesIn", args.input,
        "--outSAMunmapped", "Within",
        "--outFileNamePrefix", out_prefix,
    ]

    # Add sampled parameters
    for param, value in config.items():
        cmd.extend([param, str(value)])

    # --------------------------------------------------------
    # MARGI mode: add STAR-specific flags
    # --------------------------------------------------------
    if args.margi == 1:
        cmd.extend([
            "--chimOutType", "Junctions",
            "--chimSegmentMin", "5"
        ])

    log_file = os.path.join(out_dir, f"{prefix}log.txt")

    # Run STAR and log output
    with open(log_file, "w") as log:
        subprocess.run(cmd, stdout=log, stderr=log)


# ============================================================
# Main logic
# ============================================================

def main():
    global GENOME_DIR, THREADS, STAR_THREADS

    parser = argparse.ArgumentParser(
        description="Run Monte Carlo STAR alignments with random parameters."
    )

    parser.add_argument("--input", type=str, required=True,
                        help="Input reads file path")
    parser.add_argument("--part", type=str, default="rna",
                        help="Suffix for output directories (e.g. rna, dna)")
    parser.add_argument("--seed", type=int,
                        help="Random seed for reproducibility")
    parser.add_argument("--config-log", type=str,
                        default="generated_parameters.json",
                        help="Output file to log all sampled configs")
    parser.add_argument("--index", type=str, default=None,
                        help="Override genome directory (default: index)")
    parser.add_argument("--threads", type=int, default=None,
                        help="Number of parallel STAR runs")
    parser.add_argument("--star-threads", type=int, default=None,
                        help="Threads per STAR process (default=1)")
    parser.add_argument("--runs", type=int, default=100,
                        help="Number of Monte Carlo runs (default=100)")

    # --------------------------
    # Add MARGI mode option
    # --------------------------
    parser.add_argument(
        "-m", "--margi",
        type=int,
        default=0,
        help="MARGI mode: 1 adds STAR chimera flags, 0 disables (default=0)"
    )

    args = parser.parse_args()

    # Validate margi mode
    if args.margi not in (0, 1):
        parser.error("--margi must be 0 or 1")

    if args.margi == 1:
        print("MARGI mode ON → adding STAR chimera flags")
    else:
        print("MARGI mode OFF")

    # ---- Apply overrides ----
    if args.index is not None:
        GENOME_DIR = args.index
        print(f"Using custom genome directory: {GENOME_DIR}")

    if args.threads is not None:
        THREADS = args.threads
        print(f"Using {THREADS} parallel executor threads.")

    if args.star_threads is not None:
        if args.star_threads < 1:
            parser.error("--star-threads must be at least 1")
        STAR_THREADS = args.star_threads
    args.star_threads = STAR_THREADS

    if args.seed is not None:
        random.seed(args.seed)
        print(f"Using random seed: {args.seed}")

    num_runs = args.runs
    print(f"Launching {num_runs} STAR runs with {THREADS} parallel workers "
          f"and {STAR_THREADS} threads per STAR process.")

    # ---- Generate random configurations ----
    configs = [sample_random_config(i + (args.seed or 1)) for i in range(num_runs)]

    # ---- Save configs ----
    with open(args.config_log, "w") as f:
        json.dump(configs, f, indent=4)
    print(f"Saved all sampled configurations to {args.config_log}")

    # ---- Parallel execution ----
    with ThreadPoolExecutor(max_workers=THREADS) as executor:
        futures = [
            executor.submit(run_star, run_id, config, args)
            for run_id, config in enumerate(configs)
        ]

        for i, future in enumerate(as_completed(futures), 1):
            try:
                future.result()
                print(f"[{i}/{num_runs}] run completed")
            except Exception as exc:
                print(f"[{i}/{num_runs}] run failed: {exc}")


# ============================================================
# Entrypoint
# ============================================================

if __name__ == "__main__":
    main()

