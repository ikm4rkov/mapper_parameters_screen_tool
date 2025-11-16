#!/usr/bin/env python3
import argparse
import subprocess
import sys
import os

def main():
    parser = argparse.ArgumentParser(
        description="Master script for selecting and running mapper worker scripts."
    )

    # Common control options
    parser.add_argument("--imargi", required=True, type=str,
                        help="Set to true or false. If true, exits with 'not implemented yet'.")
    parser.add_argument("--mode", required=True, type=str,
                        choices=["bowtie2", "star", "hisat2", "bwa"],
                        help="Select mapper: bowtie2, star, hisat2, or bwa.")

    # Bowtie2 / STAR options
    parser.add_argument("--input", type=str, help="Path to input FASTQ file.")
    parser.add_argument("--part", type=str, choices=["dna", "rna"], help="Type of sequencing part.")
    parser.add_argument("--seed", type=str, help="Random seed for generation.")
    parser.add_argument("--config-log", type=str, help="Optional rename of log file.")
    parser.add_argument("--index", type=str, help="Path to genome index.")
    parser.add_argument("--threads", type=str, help="Number of threads.")

    # BWA / HISAT2 options
    parser.add_argument("-f", type=str, help="Comma-separated list of FASTQ input files.")
    parser.add_argument("-r", type=str, help="Path to reference genome index.")
    parser.add_argument("-t", type=str, help="Threads number.")
    parser.add_argument("-o", type=str, help="Output directory.")

    args = parser.parse_args()

    # Handle --imargi
    if args.imargi.lower() == "true":
        print("not implemented yet")
        sys.exit(0)
    elif args.imargi.lower() != "false":
        print("Error: --imargi must be either true or false.")
        sys.exit(1)

    # Get directory of master script
    script_dir = os.path.dirname(os.path.abspath(__file__))

    # Dispatch based on mode
    if args.mode == "bowtie2":
        worker_script = os.path.join(script_dir, "bowtie2_montecarlo.py")
        cmd = [
            "python3", worker_script,
            "--input", args.input,
            "--part", args.part,
            "--seed", args.seed,
            "--index", args.index,
            "--threads", args.threads
        ]
        if args.config_log:
            cmd += ["--config-log", args.config_log]

    elif args.mode == "star":
        worker_script = os.path.join(script_dir, "star_montecarlo.py")
        cmd = [
            "python3", worker_script,
            "--input", args.input,
            "--part", args.part,
            "--seed", args.seed,
            "--index", args.index,
            "--threads", args.threads
        ]
        if args.config_log:
            cmd += ["--config-log", args.config_log]

    elif args.mode == "bwa":
        worker_script = os.path.join(script_dir, "bwa_mapping.sh")
        cmd = [
            "bash", worker_script,
            "-f", args.f,
            "-r", args.r,
            "-t", args.t,
            "-o", args.o
        ]

    elif args.mode == "hisat2":
        worker_script = os.path.join(script_dir, "hisat2_mapping.sh")
        cmd = [
            "bash", worker_script,
            "-f", args.f,
            "-r", args.r,
            "-t", args.t,
            "-o", args.o
        ]

    else:
        print(f"Unknown mode: {args.mode}")
        sys.exit(1)

    # Check if worker script exists
    if not os.path.exists(worker_script):
        print(f"Error: Worker script not found: {worker_script}")
        sys.exit(1)

    # Print and execute
    print("Executing:", " ".join(cmd))
    subprocess.run(cmd, check=True)

if __name__ == "__main__":
    main()

