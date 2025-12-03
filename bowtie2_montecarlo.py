import random
import subprocess
import os
import argparse
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor, as_completed

# ============================================
# Constants and parameter space
# ============================================

BOWTIE2_EXEC = "bowtie2"
THREADS = 5
BOWTIE_THREADS = 1
INDEX = "Sscrofa11_1"
CONFIG_LOG_FILE = "all_sampled_configs.tsv"

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

# ============================================
# Helpers
# ============================================

def sample_random_config():
    return {param: random.choice(values) for param, values in PARAMS.items()}


def sanitize_value(v: str) -> str:
    return v.replace(",", "-").replace("=", "-").replace(".", "_")


def make_config_dirname(config):
    parts = []
    for k, v in config.items():
        sv = sanitize_value(str(v))
        parts.append(f"{k.strip('-')}{sv}")
    return "_".join(parts)


def run_bowtie2(run_id, config, args):
    # Directory name based on params
    dir_name = make_config_dirname(config)
    out_dir = Path(args.output_root) / dir_name
    out_dir.mkdir(parents=True, exist_ok=True)

    sam_out = out_dir / "output.sam"
    log_file = out_dir / "log.txt"
    config_file = out_dir / "config.txt"

    # Base bowtie2 command
    cmd = [
        BOWTIE2_EXEC,
        "--end-to-end",
        "-x", INDEX,
        "-U", args.fastq_inputs,
        "-k", "5",
        "-p", str(args.bowtie_threads),
        "-S", str(sam_out),
    ]

    if args.margi == 1:
        cmd.append("-SP5M")

    # Write config file and extend command
    with open(config_file, "w") as cf:
        for param, value in config.items():
            cf.write(f"{param} {value}\n")
            cmd.extend([param, str(value)])

    # Run bowtie2
    with open(log_file, "w") as log:
        subprocess.run(cmd, stdout=log, stderr=log)


def save_all_configs(configs, output_path):
    if not configs:
        return
    with open(output_path, "w") as f:
        headers = list(configs[0].keys())
        f.write("run_id\t" + "\t".join(headers) + "\n")
        for i, cfg in enumerate(configs):
            row = [str(cfg[h]) for h in headers]
            f.write(f"{i}\t" + "\t".join(row) + "\n")


# ============================================
# Main
# ============================================

def main():
    global INDEX, THREADS, BOWTIE_THREADS

    parser = argparse.ArgumentParser(
        description="Run Monte Carlo Bowtie2 alignments with random parameters."
    )

    parser.add_argument("--fastq-inputs", required=True,
                        help="FASTQ input file")
    parser.add_argument("--part", type=str, default="dna")
    parser.add_argument("--seed", type=int)
    parser.add_argument("--config-log", type=str, default=CONFIG_LOG_FILE)
    parser.add_argument("--index", type=str)
    parser.add_argument("--threads", type=int)
    parser.add_argument("--bowtie-threads", type=int)
    parser.add_argument("--margi", type=int, choices=[0, 1], default=0)
    parser.add_argument("--output-root", type=str, default=".",
                        help="Directory where subdirectories will be created")

    args = parser.parse_args()

    # Apply arguments
    if args.index:
        INDEX = args.index

    if args.threads:
        THREADS = args.threads

    if args.bowtie_threads:
        if args.bowtie_threads < 1:
            parser.error("--bowtie-threads must be >= 1")
        BOWTIE_THREADS = args.bowtie_threads
    args.bowtie_threads = BOWTIE_THREADS

    if args.seed is not None:
        random.seed(args.seed)

    # Sample configurations
    configs = [sample_random_config() for _ in range(100)]
    save_all_configs(configs, args.config_log)

    # Run in parallel
    with ThreadPoolExecutor(max_workers=THREADS) as executor:
        futures = [
            executor.submit(run_bowtie2, run_id, config, args)
            for run_id, config in enumerate(configs)
        ]

        for i, fut in enumerate(as_completed(futures), 1):
            try:
                fut.result()
                print(f"[{i}/100] run completed")
            except Exception as exc:
                print(f"[{i}/100] run failed: {exc}")


if __name__ == "__main__":
    main()

