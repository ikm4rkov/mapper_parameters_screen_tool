import random
import subprocess
import os
import argparse
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor, as_completed

# ============================================================
# Constants and parameter space
# ============================================================

BOWTIE2_EXEC = "bowtie2"
THREADS = 5              # Number of concurrent bowtie2 jobs
BOWTIE_THREADS = 1        # Threads per single Bowtie2 process (-p)
INDEX = "Sscrofa11_1"
OUT_ROOT = "bowtie2_wrapper_test_bams1"
CONFIG_LOG_FILE = "all_sampled_configs.tsv"

# Parameter space
PARAMS = {
    "-D": [5, 10, 15, 20, 25],
    "-R": [1, 2, 3, 4],
    "-N": [0, 1],
    "-L": [10, 20, 24, 25],
    "-i": ["S,0,2.50", "S,1,1.15", "S,1,0.50", "S,1,0.30", "S,0.5,1.50"],
    "--ma": [1, 2, 3],
    "--mp": ["4,2", "6,2", "6,3", "8,2"],
    "--score-min": [
        "L,-1.0,-0.6", "L,0,-0.6", "L,1.0,-0.6", "L,2.0,-0.6",
        "L,3.0,-0.6", "L,4.0,-0.5", "L,1.0,-0.8", "L,2.0,-0.8",
    ],
}

# ============================================================
# Helper functions
# ============================================================

def sample_random_config():
    """Sample a random parameter configuration."""
    return {param: random.choice(values) for param, values in PARAMS.items()}


def build_output_dir(run_id, part):
    """Create output directory group (10 runs per folder)."""
    group = run_id // 10 + 1
    dir_name = os.path.join(OUT_ROOT, f"{group}{part}")
    Path(dir_name).mkdir(parents=True, exist_ok=True)
    return dir_name


def run_bowtie2(run_id, config, args):
    """Run a single Bowtie2 mapping with given configuration."""
    out_dir = build_output_dir(run_id, args.part)

    prefix = f"run_{run_id:04d}_"
    sam_out = os.path.join(out_dir, f"{prefix}alignments.sam")
    log_file = os.path.join(out_dir, f"{prefix}log.txt")
    params_file = os.path.join(out_dir, f"{prefix}params.txt")

    # Base Bowtie2 command
    cmd = [
        BOWTIE2_EXEC,
        "--end-to-end",
        "-x", INDEX,
        "-U", args.input,
        "-k", "5",
        "-p", str(args.bowtie_threads),
        "-S", sam_out,
    ]

    # Add MARGI mode if enabled
    if args.margi == 1:
        cmd.append("-SP5M")

    # Write sampled parameters and extend command
    with open(params_file, "w") as param_log:
        for param, value in config.items():
            cmd.extend([param, str(value)])
            param_log.write(f"{param} {value}\n")

    # Execute Bowtie2 and log output
    with open(log_file, "w") as log:
        subprocess.run(cmd, stdout=log, stderr=log)


def save_all_configs(configs, output_path):
    """Save all sampled configurations to a TSV file."""
    if not configs:
        return
    with open(output_path, "w") as f:
        headers = list(configs[0].keys())
        f.write("run_id\t" + "\t".join(headers) + "\n")
        for i, config in enumerate(configs):
            values = [str(config[param]) for param in headers]
            f.write(f"{i}\t" + "\t".join(values) + "\n")


# ============================================================
# Main logic
# ============================================================

def main():
    global INDEX, THREADS, BOWTIE_THREADS

    parser = argparse.ArgumentParser(
        description="Run Monte Carlo Bowtie2 alignments with random parameters."
    )
    parser.add_argument("--input", type=str, required=True,
                        help="Input reads file path")
    parser.add_argument("--part", type=str, default="dna",
                        help="Suffix for output directories (e.g. dna, rna)")
    parser.add_argument("--seed", type=int,
                        help="Random seed for reproducibility")
    parser.add_argument("--config-log", type=str, default=CONFIG_LOG_FILE,
                        help="Output file to log all sampled configs")
    parser.add_argument("--index", type=str, default=None,
                        help="Override genome index (default: Sscrofa11_1)")
    parser.add_argument("--threads", type=int, default=None,
                        help="Number of parallel Bowtie2 runs")
    parser.add_argument("--bowtie-threads", type=int, default=None,
                        help="Threads per Bowtie2 process (-p)")
    parser.add_argument("--margi", type=int, choices=[0, 1], default=0,
                        help="Enable MARGI mode (-SP5M). 1=on, 0=off (default: 0)")

    args = parser.parse_args()

    # ---- Apply arguments ----
    if args.index is not None:
        INDEX = args.index
        print(f"Using custom genome index: {INDEX}")

    if args.threads is not None:
        THREADS = args.threads
        print(f"Using {THREADS} parallel executor threads.")

    if args.bowtie_threads is not None:
        if args.bowtie_threads < 1:
            parser.error("--bowtie-threads must be at least 1")
        BOWTIE_THREADS = args.bowtie_threads
    args.bowtie_threads = BOWTIE_THREADS

    if args.seed is not None:
        random.seed(args.seed)
        print(f"Using seed: {args.seed}")

    if args.margi == 1:
        print("MARGI mode enabled: adding -SP5M to Bowtie2 command.")

    print(f"Launching 100 Bowtie2 runs with {THREADS} workers "
          f"and {BOWTIE_THREADS} threads per process (-p).")

    # ---- Generate configs and save ----
    configs = [sample_random_config() for _ in range(100)]
    save_all_configs(configs, args.config_log)

    # ---- Parallel execution ----
    with ThreadPoolExecutor(max_workers=THREADS) as executor:
        futures = [
            executor.submit(run_bowtie2, run_id, config, args)
            for run_id, config in enumerate(configs)
        ]

        for i, future in enumerate(as_completed(futures), 1):
            try:
                future.result()
                print(f"[{i}/100] run completed")
            except Exception as exc:
                print(f"[{i}/100] run failed: {exc}")


# ============================================================
# Entrypoint
# ============================================================

if __name__ == "__main__":
    main()

