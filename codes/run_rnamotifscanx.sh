#!/bin/bash
#SBATCH -p intel
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 16
#SBATCH --mem=16GB
#SBATCH -t 24:00:00
#SBATCH -J rnamotif_prepare_scan
#SBATCH -o slurm-%j.out

set -euo pipefail

echo "Job ${SLURM_JOB_ID:-N/A} on ${HOSTNAME} (cwd: $(pwd))"

# ========= USER PATHS =========
PDB_ROOT="/home/s081p868/scratch/RNA_Structure_Evaluation/predictions/farfar2/str"
DATA_DIR="/home/s081p868/scratch/RNA_Structure_Evaluation/data"
OUT_ROOT="/home/s081p868/scratch/RNA_Structure_Evaluation/RNAMotifScanX_out/farfar_pdb"
FASTA_ALL="${DATA_DIR}/human_seqs_non3d_rfams.fa"

# ========= RNAMotifScanX ENV =========
export RNAMOTIFSCANX_PATH="/home/s081p868/scratch/RNAMotifScanX-release"
export RNAVIEW="$RNAMOTIFSCANX_PATH/thirdparty/RNAVIEW"

SCAN_BIN="${RNAMOTIFSCANX_PATH}/bin/scan"
MODELS_DIR="${RNAMOTIFSCANX_PATH}/models"

module load anaconda3
source activate py27

# Use Python 2.7 explicitly
if [[ -n "${CONDA_PREFIX:-}" && -x "${CONDA_PREFIX}/bin/python" ]]; then
  PY_EXE="${CONDA_PREFIX}/bin/python"
else
  PY_EXE="$(command -v python2 || true)"
fi
if [[ -z "${PY_EXE:-}" ]]; then
  echo "ERROR: Could not find Python 2.7 interpreter. Is env 'py27' available?" >&2
  exit 2
fi
echo "[debug] Using Python: ${PY_EXE} ($(${PY_EXE} -V 2>&1))"

# Optional RNAMotifScanX env
if [[ -f "${RNAMOTIFSCANX_PATH}/set_env.sh" ]]; then
  source "${RNAMOTIFSCANX_PATH}/set_env.sh" "${RNAMOTIFSCANX_PATH}"
fi

mkdir -p "${OUT_ROOT}"

# ---- Constant 4-char tag + chain A ----
TEMP_TAG="ABCD"   # FASTA header will be >ABCD_A

# Extract one record and write FASTA with header ">${TEMP_TAG}_A" into OUT_ROOT
extract_to_temp_fa() {
  local header="$1" outfa="$2"
  awk -v ID="$header" -v TAG="${TEMP_TAG}" '
    BEGIN{want=0}
    /^>/{
      hdr=$0; sub(/^>/,"",hdr); sub(/ .*/, "", hdr)
      if (hdr==ID) { print ">" TAG "_A"; want=1; next } else { want=0 }
    }
    want{ print }
  ' "${FASTA_ALL}" > "${outfa}"
}

# Snapshot helpers (list files/dirs at depth 1 of OUT_ROOT)
snap_files() { find "${OUT_ROOT}" -maxdepth 1 -mindepth 1 -type f -printf '%f\n' | sort; }
snap_dirs()  { find "${OUT_ROOT}" -maxdepth 1 -mindepth 1 -type d -printf '%f\n' | sort; }

# Run RNAMotifScanX scan for a given workdir (per-case folder)
run_scan() {
  local workdir="$1"
  local models_dir="${MODELS_DIR}"
  local pdb_dir="${workdir}"
  local out_dir="${workdir}/Res_motifs_orig"

  mkdir -p "${out_dir}"

  # Choose inputs produced by PrepareInput.py
  # Use the first *.rmsx.in_ and the first *.rmsx.nch in workdir
  # (adjust patterns if your PrepareInput names differ)
  local in_file nch_file
  in_file="$(ls -1 "${pdb_dir}"/*.rmsx.in 2>/dev/null | head -n1 || true)"
  nch_file="$(ls -1 "${pdb_dir}"/*.rmsx.nch 2>/dev/null | head -n1 || true)"

  if [[ -z "${in_file}" || -z "${nch_file}" ]]; then
    echo "WARNING: Could not find .rmsx.in_ and/or .rmsx.nch in ${pdb_dir}. Skipping scan."
    return 0
  fi

  # Loop all model .struct files
  shopt -s nullglob
  local struct_file base_name
  for struct_file in "${models_dir}"/*.struct; do
    base_name="$(basename "${struct_file}" .struct)"
    mkdir -p "${out_dir}/${base_name}"
    echo "[scan] ${base_name}"
    # Prefer ${SCAN_BIN}; fallback to ./bin/scan if not found
    if [[ -x "${SCAN_BIN}" ]]; then
      "${SCAN_BIN}" "${struct_file}" "${in_file}" --map_pdb="${nch_file}" --num_threads "${SLURM_CPUS_PER_TASK:-1}" > "${out_dir}/${base_name}/result.log"
    else
      ./bin/scan "${struct_file}" "${in_file}" --map_pdb="${nch_file}" --num_threads "${SLURM_CPUS_PER_TASK:-1}" > "${out_dir}/${base_name}/result.log"
    fi
  done
}

# Collect PDBs
mapfile -t all_pdbs < <(find "${PDB_ROOT}" -type f -name "*.pdb" | sort)
if [[ ${#all_pdbs[@]} -eq 0 ]]; then
  echo "No .pdb files found under ${PDB_ROOT}"
  exit 1
fi

# Support array runs: if SLURM_ARRAY_TASK_ID is set, process only that index; else process all
if [[ -n "${SLURM_ARRAY_TASK_ID:-}" ]]; then
  start_idx="${SLURM_ARRAY_TASK_ID}"
  end_idx="${SLURM_ARRAY_TASK_ID}"
else
  start_idx=0
  end_idx=$((${#all_pdbs[@]} - 1))
fi

for (( idx=${start_idx}; idx<=${end_idx}; idx++ )); do
  pdb="${all_pdbs[$idx]}"
  base="$(basename "${pdb}" .pdb)"
  echo "=== [${idx}] Processing: ${base} ==="

  workdir="${OUT_ROOT}/${base}"
  mkdir -p "${workdir}"

  # Clean leftovers with our tag in OUT_ROOT
  rm -f "${OUT_ROOT}/${TEMP_TAG}.fa" \
        "${OUT_ROOT}/${TEMP_TAG}.pdb" \
        "${OUT_ROOT}/${TEMP_TAG}.pdb.mca" \
        "${OUT_ROOT}/${TEMP_TAG}.pdb.out" || true

  # Create fresh ABCD.fa / ABCD.pdb in OUT_ROOT
  tmp_fa="${OUT_ROOT}/${TEMP_TAG}.fa"
  tmp_pdb="${OUT_ROOT}/${TEMP_TAG}.pdb"

  extract_to_temp_fa "${base}" "${tmp_fa}"
  if ! grep -q "^>${TEMP_TAG}_A$" "${tmp_fa}"; then
    echo "WARNING: No FASTA record found for '${base}' in ${FASTA_ALL} — skipping."
    rm -f "${tmp_fa}"
    continue
  fi
  cp -f "${pdb}" "${tmp_pdb}"

  # ---- Snapshot OUT_ROOT before run ----
  before_files="$(mktemp)"; before_dirs="$(mktemp)"
  snap_files > "${before_files}"
  snap_dirs  > "${before_dirs}"

  # ---- Run PrepareInput inside OUT_ROOT ----
  (
    cd "${OUT_ROOT}"
    echo "[debug] PrepareInput in: $(pwd)"
    ls -l "${TEMP_TAG}.pdb" "${TEMP_TAG}.fa" || true
    rm -f "${TEMP_TAG}.pdb.mca" "${TEMP_TAG}.pdb.out" || true
    "${PY_EXE}" "${RNAMOTIFSCANX_PATH}/scripts/PrepareInput.py" "${TEMP_TAG}.pdb" "${TEMP_TAG}.fa"
  )

  # ---- Snapshot OUT_ROOT after run ----
  after_files="$(mktemp)"; after_dirs="$(mktemp)"
  snap_files > "${after_files}"
  snap_dirs  > "${after_dirs}"

  # ---- Compute new items (files & dirs) ----
  new_files="$(mktemp)"; new_dirs="$(mktemp)"
  comm -13 "${before_files}" "${after_files}" > "${new_files}" || true
  comm -13 "${before_dirs}"  "${after_dirs}"  > "${new_dirs}"  || true

  # ---- Copy NEW FILES into workdir, then delete from OUT_ROOT ----
  while IFS= read -r f; do
    [[ -z "${f}" ]] && continue
    cp -f "${OUT_ROOT}/${f}" "${workdir}/"
    rm -f "${OUT_ROOT}/${f}"
  done < "${new_files}"

  # ---- Copy NEW DIRS into workdir, then delete from OUT_ROOT ----
  while IFS= read -r d; do
    [[ -z "${d}" ]] && continue
    if [[ -d "${workdir}/${d}" ]]; then
      cp -rf "${OUT_ROOT}/${d}/." "${workdir}/${d}/"
      rm -rf "${OUT_ROOT}/${d}"
    else
      cp -r "${OUT_ROOT}/${d}" "${workdir}/"
      rm -rf "${OUT_ROOT}/${d}"
    fi
  done < "${new_dirs}"

  # Ensure inputs ABCD.* also land in workdir (for auditing), then remove from OUT_ROOT
  cp -f "${OUT_ROOT}/${TEMP_TAG}.fa"  "${workdir}/${TEMP_TAG}.fa"  || true
  cp -f "${OUT_ROOT}/${TEMP_TAG}.pdb" "${workdir}/${TEMP_TAG}.pdb" || true
  rm -f "${OUT_ROOT}/${TEMP_TAG}.fa" "${OUT_ROOT}/${TEMP_TAG}.pdb" || true

  # Cleanup temp lists
  rm -f "${before_files}" "${before_dirs}" "${after_files}" "${after_dirs}" "${new_files}" "${new_dirs}"

  # =============== RUN SCAN for this case ===============
  echo "[scan] models_dir=${MODELS_DIR}"
  echo "[scan] pdb_dir=${workdir}"
  echo "[scan] output_dir=${workdir}/Res_motifs"
  run_scan "${workdir}"

  echo "Done: ${base} → results in ${workdir}"
done

echo "All done. Outputs under: ${OUT_ROOT}"
