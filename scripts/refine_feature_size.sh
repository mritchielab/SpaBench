#!/usr/bin/env bash
set -euo pipefail
module load snakemake
module load python
module load R

# run_refine.sh
# Usage:
#   ./run_refine.sh                 # process all samples under ./results/
#   ./run_refine.sh -s SAMPLE_NAME  # process single sample by name
#   ./run_refine.sh -P SAMPLE_PATH  # process single sample by full path
#   ./run_refine.sh -s NAME -j 24 -c 4  # adjust snakemake jobs/cores
#   ./run_refine.sh -n "500 1000"    # only process specified top_n values (no default)
#   ./run_refine.sh -m "BANKSY,leiden" # only run specified clustering methods
#
# Options:
#   -s SAMPLE_NAME    Sample directory name under results/ (e.g. Mouse_brain_squidpy)
#   -P SAMPLE_PATH    Full path to a sample directory (e.g. /vast/.../results/Mouse_brain_squidpy)
#   -j JOBS           snakemake -j (default 48)
#   -c CORES          cores for metrics_summary (default 8)  # kept for compatibility though metrics_summary not used
#   -n "N1 N2 ..."    space- or comma-separated list of top_n values to run (no default)
#   -m "M1 M2 ..."    space- or comma-separated list of clustering methods to run (no default)
#   -h                show this help
# --------------------------------------------------

print_help() {
  sed -n '1,200p' "$0" | sed -n '1,200p' >/dev/stderr
  cat <<'EOF' >/dev/stderr

Examples:
  ./run_refine.sh
  ./run_refine.sh -s Mouse_brain_squidpy
  ./run_refine.sh -P /vast/.../results/Mouse_brain_squidpy -j 32 -c 8
  ./run_refine.sh -n "500 1000 2000"
  ./run_refine.sh -m "BANKSY,leiden,louvain_refine" -n 500 -s Mouse_brain_squidpy -j 80 -c 12
EOF
}

# default settings (can be overridden by options)
ROOT_BASE_DEFAULT="/vast/projects/SpatialTranscriptomicsBenchmarking/Snakemake_workflow/results"
SMK_JOBS=48
SMK_SUMMARY_CORES=8

# arrays (edit if needed)
NORMS=(LogNormalize CLR RC SCTransform SpaNorm)
METHODS=(moransi nnSVG markvariogram SPARK "SPARK-X" vst disp mvp)
# default ALGS; may be overridden by -m option
ALGS=(BANKSY BayesSpace louvain leiden SLM louvain_refine SpaGCN SpatialPCA stlearn SpatialLeiden)
SIZES=(500 1000 2000 5000 10000)

# allowed clustering method names for validation
ALLOWED_METHODS=(BANKSY leiden louvain_refine SpaGCN SpatialPCA BayesSpace louvain SLM SpatialLeiden stlearn)

# parse options
SAMPLE=""
SAMPLE_PATH=""
N_SIZES_RAW=""  # raw string from -n option (no default)
M_ALGS_RAW=""   # raw string from -m option (no default)
while getopts ":s:P:j:c:n:m:h" opt; do
  case "${opt}" in
    s) SAMPLE="${OPTARG}" ;;
    P) SAMPLE_PATH="${OPTARG}" ;;
    j) SMK_JOBS="${OPTARG}" ;;
    c) SMK_SUMMARY_CORES="${OPTARG}" ;;
    n) N_SIZES_RAW="${OPTARG}" ;;
    m) M_ALGS_RAW="${OPTARG}" ;;
    h) print_help; exit 0 ;;
    \?) echo "Invalid option: -${OPTARG}" >&2; print_help; exit 1 ;;
    :) echo "Option -${OPTARG} requires an argument." >&2; print_help; exit 1 ;;
  esac
done
shift $((OPTIND-1))

# If -n provided, parse and override SIZES (no default per your request)
if [ -n "${N_SIZES_RAW:-}" ]; then
    # allow comma-separated or space-separated input; normalize commas to spaces
    N_SIZES_RAW="${N_SIZES_RAW//,/ }"
    # read into array
    read -r -a PARSED_SIZES <<< "${N_SIZES_RAW}"
    # validate each entry is a positive integer
    SIZES=()
    for s in "${PARSED_SIZES[@]}"; do
        if ! [[ "${s}" =~ ^[0-9]+$ ]]; then
            echo "Invalid -n size value: '${s}'. Must be positive integer." >&2
            exit 1
        fi
        SIZES+=("${s}")
    done
fi

# If -m provided, parse and override ALGS (no default)
if [ -n "${M_ALGS_RAW:-}" ]; then
    M_ALGS_RAW="${M_ALGS_RAW//,/ }"
    read -r -a PARSED_ALGS <<< "${M_ALGS_RAW}"
    # validate each provided alg is allowed
    NEW_ALGS=()
    for a in "${PARSED_ALGS[@]}"; do
        # trim whitespace
        a_trimmed="$(echo "${a}" | xargs)"
        if [ -z "${a_trimmed}" ]; then
            continue
        fi
        found=0
        for allowed in "${ALLOWED_METHODS[@]}"; do
            if [ "${a_trimmed}" = "${allowed}" ]; then
                found=1
                break
            fi
        done
        if [ "${found}" -eq 0 ]; then
            echo "Invalid clustering method in -m: '${a_trimmed}'. Allowed: ${ALLOWED_METHODS[*]}" >&2
            exit 1
        fi
        NEW_ALGS+=("${a_trimmed}")
    done
    if [ "${#NEW_ALGS[@]}" -eq 0 ]; then
        echo "No valid clustering methods provided to -m." >&2
        exit 1
    fi
    ALGS=("${NEW_ALGS[@]}")
fi

# logging helper
log() { printf '%s [%s] %s\n' "$(date '+%Y-%m-%d %H:%M:%S')" "${SAMPLE:-ALL}" "$*"; }

# root base default; will be overridden if SAMPLE_PATH provided
ROOT_BASE="${ROOT_BASE_DEFAULT}"

# helper: find top csv for given sample/norm/method/size
find_top_csv() {
    local sample="$1"; local norm="$2"; local method="$3"; local size="$4"
    local d="${ROOT_BASE}/${sample}/${norm}/${method}"
    local f1="${d}/top_${size}.csv"
    local f2="${d}/${size}.csv"
    if [ -f "${f1}" ]; then echo "${f1}"; return 0; fi
    if [ -f "${f2}" ]; then echo "${f2}"; return 0; fi
    if [ -d "${d}" ]; then
        shopt -s nullglob
        for cand in "${d}"/*.csv "${d}"/*.tsv "${d}"/*.txt; do
            echo "${cand}"; shopt -u nullglob; return 0
        done
        shopt -u nullglob
    fi
    return 1
}

install_feature_csv() {
    local src="$1"; local sample="$2"; local norm="$3"; local method="$4"
    local dest_dir="results/${sample}/${norm}/${method}"
    local dest="${dest_dir}/04_${method}.csv"
    mkdir -p "${dest_dir}"
    cp -f "${src}" "${dest}"
    log "copied: ${src} -> ${dest}"
}

# touch helper: only touch if file exists (do not create new files)
touch_if_exists() {
    for f in "$@"; do
        if [ -e "${f}" ]; then
            touch "${f}"
            log "touched: ${f}"
        fi
    done
}

# touch only the cluster & annotate inputs for a given alg (based on Snakefile)
# Arguments: sample norm method alg
touch_inputs_for_alg() {
    local sample="$1"; local norm="$2"; local method="$3"; local alg="$4"

    # common derived paths
    local feat_csv="results/${sample}/${norm}/${method}/04_${method}.csv"
    local reduced_rds="data/${sample}/${norm}/${method}/04_reduced.rds"
    local select_rds="data/${sample}/${norm}/${method}/03_${method}.rds"

    # sample-level inputs
    local raw_seurat="data/${sample}/01_raw_seurat.rds"
    local raw_anndata="data/${sample}/01_raw_anndata.h5ad"
    local sc_anndata="data/${sample}/01_sc_raw_anndata.h5ad"
    local positions="data/${sample}/01_positions.csv"

    # annotation inputs
    local clustered_rds="data/${sample}/${norm}/${method}/${alg}/06_${alg}_clustered.rds"
    local clusters_csv="results/${sample}/${norm}/${method}/${alg}/07_${alg}_clusters.csv"
    local truth_file="data/Ground_truth/${sample}/labels.tsv"

    # Determine which cluster rule inputs to touch (match your Snakefile)
    case "${alg}" in
        SpaGCN)
            # SpaGCN inputs: raw_h5ad, positions, features_csv
            touch_if_exists "${raw_anndata}" "${positions}" "${feat_csv}"
            ;;
        stlearn)
            # stlearn inputs: raw_sc_anndata, features_csv
            touch_if_exists "${sc_anndata}" "${feat_csv}"
            ;;
        SpatialLeiden)
            # SpatialLeiden inputs: raw_sc_anndata, features_csv
            touch_if_exists "${sc_anndata}" "${feat_csv}"
            ;;
        BANKSY|BayesSpace|louvain|leiden|SLM|louvain_refine|SpatialPCA)
            # These use reduced rds as cluster input (data/.../04_reduced.rds)
            touch_if_exists "${reduced_rds}"
            ;;
        *)
            # default fallback: touch reduced rds and feature csv if present
            touch_if_exists "${reduced_rds}" "${feat_csv}"
            ;;
    esac

    # Touch annotation inputs (if exist) - annotate_clusters inputs in Snakefile:
    #   seurat = data/{sample}/{norm}/{feat}/{alg}/06_{alg}_clustered.rds
    #   clusters = results/{sample}/{norm}/{feat}/{alg}/07_{alg}_clusters.csv
    #   truth = data/Ground_truth/{sample}/labels.tsv
    touch_if_exists "${clustered_rds}" "${clusters_csv}" "${truth_file}"
}

# Determine samples to process
if [ -n "${SAMPLE_PATH:-}" ]; then
    if [ ! -d "${SAMPLE_PATH}" ]; then
        echo "SAMPLE_PATH set but not a directory: ${SAMPLE_PATH}" >&2
        exit 1
    fi
    ROOT_BASE="$(dirname "${SAMPLE_PATH}")"
    sample_name="$(basename "${SAMPLE_PATH}")"
    SAMPLES=("${sample_name}")
    SAMPLE="${sample_name}"
else
    if [ -n "${SAMPLE:-}" ]; then
        SAMPLES=("${SAMPLE}")
    else
        SAMPLES=()
        for d in results/*; do
            [ -d "$d" ] || continue
            sample_name=$(basename "$d")
            SAMPLES+=("${sample_name}")
        done
    fi
fi

if [ "${#SAMPLES[@]}" -eq 0 ]; then
    echo "No samples found to process. Exiting." >&2
    exit 1
fi

log "ROOT_BASE = ${ROOT_BASE}"
log "Samples to process: ${SAMPLES[*]}"
log "snakemake -j set to ${SMK_JOBS}, metrics_summary cores ${SMK_SUMMARY_CORES}"
log "SIZES to run: ${SIZES[*]}"
log "ALGS to run: ${ALGS[*]}"

# Main loop
for sample in "${SAMPLES[@]}"; do
    SAMPLE="${sample}"
    log "=== START sample: ${sample} ==="

    for size in "${SIZES[@]}"; do
        log "----- size = ${size} -----"

        # 1) copy feature files for all norms/methods for this sample
        for norm in "${NORMS[@]}"; do
            for method in "${METHODS[@]}"; do
                src="$(find_top_csv "${sample}" "${norm}" "${method}" "${size}" || true)"
                if [ -z "${src}" ]; then
                    log "WARN: no top file for sample=${sample} norm=${norm} method=${method} size=${size}; skipping."
                    continue
                fi
                install_feature_csv "${src}" "${sample}" "${norm}" "${method}"
            done
        done

        # 2) Run snakemake for each clustering algorithm, for concrete norm+method combos
        for alg in "${ALGS[@]}"; do
            for norm in "${NORMS[@]}"; do
                for method in "${METHODS[@]}"; do
                    target="results/${sample}/${norm}/${method}/${alg}/09_${alg}_metrics.csv"
                    feat_csv="results/${sample}/${norm}/${method}/04_${method}.csv"
                    if [ ! -f "${feat_csv}" ]; then
                        log "SKIP snakemake for ${target} because feature csv missing: ${feat_csv}"
                        continue
                    fi

                    # ---------------- touch inputs for cluster & annotation (per Snakefile) ----------
                    touch_inputs_for_alg "${sample}" "${norm}" "${method}" "${alg}"
                    # ----------------------------------------------------------------------------------

                    log "snakemake -> target=${target}"
                    if snakemake -j "${SMK_JOBS}" -p "${target}"; then
                        log "snakemake succeeded for target=${target}"
                        # If snakemake generated the metrics file, append it to aggregated per-size per-alg file
                        if [ -f "${target}" ]; then
                            dest_dir="results/${sample}/metrics/${size}"
                            mkdir -p "${dest_dir}"
                            dest_file="${dest_dir}/${alg}_metrics.csv"

                            # Read header of the produced metrics file
                            header_line="$(head -n 1 "${target}")"

                            if [ ! -f "${dest_file}" ]; then
                                # create aggregated file with new header: sample,norm,method,clustering + original header
                                printf '%s\n' "\"sample\",\"norm\",\"method\",\"clustering\",${header_line}" > "${dest_file}"
                                # append data lines (skip header) and prepend metadata columns
                                tail -n +2 "${target}" | awk -v s="${sample}" -v n="${norm}" -v m="${method}" -v a="${alg}" '{print "\""s"\",\""n"\",\""m"\",\""a"\"," $0}' >> "${dest_file}"
                                log "Created aggregated metrics: ${dest_file} (from ${target})"
                            else
                                # Append new rows (skip header) with metadata columns
                                tail -n +2 "${target}" | awk -v s="${sample}" -v n="${norm}" -v m="${method}" -v a="${alg}" '{print "\""s"\",\""n"\",\""m"\",\""a"\"," $0}' >> "${dest_file}"
                                log "Appended metrics to: ${dest_file} (from ${target})"
                            fi
                        else
                            log "WARN: expected metrics file not found after snakemake: ${target}"
                        fi
                    else
                        log "ERROR: snakemake failed for target=${target}. Continuing."
                    fi
                done
            done
        done

        # Note: metrics_summary step removed per request; aggregation handled incrementally above.

    done # sizes

    log "=== DONE sample: ${sample} ==="
done

log "All samples processed."
