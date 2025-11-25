#!/usr/bin/env bash
set -euo pipefail

# make_top_files.sh
# Usage:
#   ./make_top_files.sh [BASE_DIR]            # e.g. ./make_top_files.sh results/Mouse_embryo_stereoseq
#   ./make_top_files.sh --force [BASE_DIR]    # force overwrite existing top_*.csv
#   ./make_top_files.sh [BASE_DIR] 500 1000   # custom N list
#   ./make_top_files.sh --force               # use default BASE_DIR and default Ns, overwrite
#
# If BASE_DIR omitted, default is "results/Mouse_embryo_stereoseq"
# Default Ns: 500 1000 2000 5000 10000

# defaults
DEFAULT_BASE="results/Mouse_embryo_stereoseq"
FORCE=0

# parse optional --force
if [[ "${1:-}" == "--force" ]]; then
  FORCE=1
  shift
fi

# next optional arg may be BASE_DIR (a directory) or a number (N)
BASE_DIR="${DEFAULT_BASE}"
if [[ "${1:-}" != "" ]]; then
  # if first non-flag argument is a directory, treat as BASE_DIR
  if [[ -d "$1" && ! "$1" =~ ^[0-9]+$ ]]; then
    BASE_DIR="$1"
    shift
  fi
fi

# remaining args (if any) are Ns; else default list
if [ $# -ge 1 ]; then
  Ns=("$@")
else
  Ns=(500 1000 2000 5000 10000)
fi

# validate BASE_DIR
if [ ! -d "${BASE_DIR}" ]; then
  echo "ERROR: BASE_DIR not found: ${BASE_DIR}" >&2
  exit 1
fi

echo "Base directory: ${BASE_DIR}"
echo "Will create top files for N = ${Ns[*]}"
echo "Force overwrite = ${FORCE}"

# find all 04_*.csv under BASE_DIR and process
find "${BASE_DIR}" -type f -name "04_*.csv" | while IFS= read -r src; do
  dir=$(dirname "$src")
  base=$(basename "$src")
  echo "=== Processing: $src ==="

  # count total lines (including header)
  total_lines=$(wc -l < "$src" | tr -d '[:space:]' || true)
  if [ -z "$total_lines" ] || [ "$total_lines" -eq 0 ]; then
    echo "  WARN: empty file: $src, skipping."
    continue
  fi

  for N in "${Ns[@]}"; do
    # validate N is positive integer
    if ! [[ "${N}" =~ ^[0-9]+$ ]]; then
      echo "  SKIP: invalid N value: ${N}" >&2
      continue
    fi

    out="${dir}/top_${N}.csv"
    if [ -f "$out" ] && [ "${FORCE}" -ne 1 ]; then
      echo "  SKIP: $out exists (use --force to overwrite)"
      continue
    fi

    # data lines (excluding header)
    data_lines=$(( total_lines > 0 ? total_lines - 1 : 0 ))
    if [ "$data_lines" -le 0 ]; then
      echo "  NOTE: no data rows in $src; copying header only to $out"
      head -n 1 "$src" > "$out"
    elif [ "$data_lines" -lt "$N" ]; then
      echo "  NOTE: data rows ($data_lines) < $N; copying entire file to $out"
      cp -f "$src" "$out"
    else
      # header + first N data rows
      head -n 1 "$src" > "$out"
      sed -n "2,$((N+1))p" "$src" >> "$out"
      echo "  WROTE top $N -> $out"
    fi
  done
done

echo "Done. To inspect generated files:"
echo "  find \"${BASE_DIR}\" -type f -name 'top_*.csv' | sort | head"

