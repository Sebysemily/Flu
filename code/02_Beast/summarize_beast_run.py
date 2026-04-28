#!/usr/bin/env python3
import argparse
import os


def main() -> None:
    parser = argparse.ArgumentParser(description="Summarize duplicated BEAST runs into one scenario-level sentinel.")
    parser.add_argument("--scenario", required=True)
    parser.add_argument("--run-done", dest="run_dones", action="append", required=True)
    parser.add_argument("--seed", dest="seeds", action="append", required=True, type=int)
    parser.add_argument("--out", required=True)
    args = parser.parse_args()

    out_dir = os.path.dirname(args.out)
    if out_dir:
        os.makedirs(out_dir, exist_ok=True)

    with open(args.out, "w", encoding="utf-8") as handle:
        handle.write(f"scenario\t{args.scenario}\n")
        handle.write(f"replicate_count\t{len(args.run_dones)}\n")
        handle.write(f"run_done_files\t{';'.join(os.path.abspath(path) for path in args.run_dones)}\n")
        handle.write(f"seeds\t{';'.join(str(seed) for seed in args.seeds)}\n")


if __name__ == "__main__":
    main()
