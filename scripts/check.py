#!/usr/bin/env python3
# check.py
# Check presence of 06_{alg}_clustered.rds for a given sample directory (or multiple).
# Usage:
#   python check.py -s /path/to/Human_Breast_Cancer_BlockASection1
#   python check.py -s sample_dir1 sample_dir2 -c config.yaml --only-missing --missing-file miss.txt

import argparse
import os
import sys
import yaml
from pathlib import Path

# Candidate normalization folders to auto-detect if user not specifying --norms
CANDIDATE_NORMS = ["LogNormalize", "SCTransform", "RC", "CLR", "SpaNorm"]

def load_config(cfg_path):
    if not cfg_path:
        return {}
    p = Path(cfg_path)
    if not p.exists():
        print(f"Warning: config file {cfg_path} not found. Using command-line arguments only.", file=sys.stderr)
        return {}
    with p.open() as fh:
        return yaml.safe_load(fh)

def detect_norms(sample_dir, explicit_norms=None):
    if explicit_norms:
        return explicit_norms
    found = []
    for n in CANDIDATE_NORMS:
        if (Path(sample_dir) / n).is_dir():
            found.append(n)
    # if nothing found, fallback to candidates (so we check all possible)
    return found if found else CANDIDATE_NORMS

def build_expected_paths_for_sample(sample_dir, norms, feats, algs):
    expected = []
    sample_dir = Path(sample_dir)
    sample_name = sample_dir.name
    for norm in norms:
        for feat in feats:
            for alg in algs:
                # path: SAMPLE_DIR/{norm}/{feat}/{alg}/06_{alg}_clustered.rds
                p = sample_dir / norm / feat / alg / f"06_{alg}_clustered.rds"
                expected.append({
                    "sample_dir": str(sample_dir),
                    "sample": sample_name,
                    "norm": norm,
                    "feat": feat,
                    "alg": alg,
                    "path": str(p)
                })
    return expected

def main():
    p = argparse.ArgumentParser(description="Check 06_{alg}_clustered.rds presence for sample dir(s)")
    p.add_argument("-s", "--samples", nargs="+", required=True,
                   help="One or more sample directories (full path). Example: -s data/Human_Breast_Cancer_BlockASection1")
    p.add_argument("-c", "--config", default="config.yaml", help="Path to config.yaml (default: config.yaml)")
    p.add_argument("--feats", nargs="+", help="Optional override for features list (space-separated)")
    p.add_argument("--algs", nargs="+", help="Optional override for algorithm keys (space-separated)")
    p.add_argument("--norms", nargs="+", help="Optional override for normalization methods to check (space-separated)")
    p.add_argument("--only-missing", action="store_true", help="Only print missing paths")
    p.add_argument("--missing-file", help="Write missing paths to this file (one per line)")
    p.add_argument("--exist-file", help="Write existing paths to this file (one per line)")
    args = p.parse_args()

    cfg = load_config(args.config)

    # features: priority: --feats > config.features.methods
    if args.feats:
        feats = args.feats
    else:
        feats = cfg.get("features", {}).get("methods", [])
        if not feats:
            print("ERROR: No features specified (use --feats or provide features.methods in config.yaml)", file=sys.stderr)
            sys.exit(1)

    # algorithms: priority --algs > config.algorithms keys
    if args.algs:
        algs = args.algs
    else:
        alg_map = cfg.get("algorithms", {})
        if not alg_map:
            print("ERROR: No algorithms found in config.yaml (or use --algs)", file=sys.stderr)
            sys.exit(1)
        algs = list(alg_map.keys())

    # For each sample dir passed, build expected paths (auto-detect norms if not provided)
    all_expected = []
    for sd in args.samples:
        if not os.path.isdir(sd):
            print(f"WARNING: sample path does not exist or not a directory: {sd}", file=sys.stderr)
            continue
        norms = detect_norms(sd, explicit_norms=args.norms)
        expected = build_expected_paths_for_sample(sd, norms, feats, algs)
        all_expected.extend(expected)

    total = len(all_expected)
    existing = []
    missing = []

    for e in all_expected:
        if os.path.exists(e["path"]):
            existing.append(e["path"])
        else:
            missing.append(e["path"])

    # Summary per sample (grouped)
    from collections import defaultdict
    per_sample = defaultdict(lambda: {"total":0, "existing":0, "missing":0})
    for e in all_expected:
        s = e["sample_dir"]
        per_sample[s]["total"] += 1
        if os.path.exists(e["path"]):
            per_sample[s]["existing"] += 1
        else:
            per_sample[s]["missing"] += 1

    # Print results
    print("="*60)
    print(f"Checked {total} expected '06_*_clustered.rds' files for {len(per_sample)} sample dir(s).")
    print(f"  Found: {len(existing)}    Missing: {len(missing)}")
    print("-"*60)
    for sd, stats in per_sample.items():
        print(f"Sample dir: {sd}")
        print(f"  expected: {stats['total']}, found: {stats['existing']}, missing: {stats['missing']}")
    print("="*60)

    if args.only_missing:
        for pth in missing:
            print(pth)
    else:
        if missing:
            print("\nMissing files (first 200 shown):")
            for pth in missing[:200]:
                print("  " + pth)
        else:
            print("\nNo missing files.")

    if args.exist_file:
        with open(args.exist_file, "w") as fh:
            fh.write("\n".join(existing))
        print(f"Wrote existing-file list to {args.exist_file}")

    if args.missing_file:
        with open(args.missing_file, "w") as fh:
            fh.write("\n".join(missing))
        print(f"Wrote missing-file list to {args.missing_file}")

    # exit non-zero when missing (useful for CI)
    if missing:
        sys.exit(2)
    else:
        sys.exit(0)

if __name__ == "__main__":
    main()

