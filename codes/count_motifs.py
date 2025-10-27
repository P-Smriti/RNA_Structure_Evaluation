#!/usr/bin/env python3
"""
Update fasta_mapping_with_length.csv with num_motifs_farfar2 counts.

Steps:
  1) Walk FARFAR2 output under farfar_pdb/URS*/Res_motifs_orig/*/result.log
     and count motif hits as (# of non-header, non-blank lines per result.log).
  2) Sum per-URS across all motif subfolders => num_motifs_farfar2.
  3) Overwrite the 'num_motifs_farfar2' column for matching file_name rows in
     the provided mapping CSV, preserving order and other columns.

Usage:
  python3 update_farfar2_counts_in_mapping.py \
    --farfar-root /home/s081p868/scratch/RNA_Structure_Evaluation/RNAMotifScanX_out/farfar_pdb \
    --mapping-csv /home/s081p868/scratch/RNA_Structure_Evaluation/data/fasta_mapping_with_length.csv \
    --out-csv /home/s081p868/scratch/RNA_Structure_Evaluation/data/fasta_mapping_with_length.updated.csv

To update in-place:
  python3 update_farfar2_counts_in_mapping.py --farfar-root ... \
    --mapping-csv ... --out-csv ...
  mv ...updated.csv ...with_length.csv
"""

import argparse
import csv
from pathlib import Path
from collections import defaultdict

def count_results(log_path: Path) -> int:
    """Count non-header, non-blank lines in a result.log."""
    if not log_path.is_file():
        return 0
    n = 0
    with log_path.open("r", encoding="utf-8", errors="ignore") as f:
        for line in f:
            s = line.strip()
            if not s or s.startswith("#"):
                continue
            n += 1
    return n

def collect_farfar_counts(farfar_root: Path) -> dict:
    """
    Return { 'URS...': total_count_across_all_motifs } by scanning:
      farfar_root/URS*/Res_motifs_orig/*/result.log
    """
    counts = defaultdict(int)
    urs_dirs = sorted([p for p in farfar_root.glob("URS*") if p.is_dir()])
    for urs in urs_dirs:
        motifs_root = urs / "Res_motifs_orig"
        if not motifs_root.is_dir():
            continue
        total = 0
        for motif_dir in motifs_root.iterdir():
            if not motif_dir.is_dir():
                continue
            total += count_results(motif_dir / "result.log")
        counts[urs.name] = total
    return dict(counts)

def main():
    ap = argparse.ArgumentParser(description="Update mapping CSV with FARFAR2 motif counts")
    ap.add_argument("--farfar-root", required=True, help="Path to farfar_pdb directory")
    ap.add_argument("--mapping-csv", required=True, help="Existing mapping CSV to update")
    ap.add_argument("--out-csv", required=True, help="Where to write the updated CSV")
    args = ap.parse_args()

    farfar_root = Path(args.farfar_root).expanduser().resolve()
    mapping_csv = Path(args.mapping_csv).expanduser().resolve()
    out_csv = Path(args.out_csv).expanduser().resolve()

    if not farfar_root.is_dir():
        raise SystemExit(f"[ERROR] Not a directory: {farfar_root}")
    if not mapping_csv.is_file():
        raise SystemExit(f"[ERROR] Mapping CSV not found: {mapping_csv}")

    # 1) Collect counts
    farfar_counts = collect_farfar_counts(farfar_root)
    print(f"[INFO] Collected FARFAR2 counts for {len(farfar_counts)} URS folders.")

    # 2) Read mapping CSV, overwrite num_motifs_farfar2 by file_name
    with mapping_csv.open("r", newline="", encoding="utf-8") as f:
        reader = csv.DictReader(f)
        fieldnames = list(reader.fieldnames or [])
        if "file_name" not in fieldnames:
            raise SystemExit("[ERROR] mapping CSV missing 'file_name' column.")
        if "num_motifs_farfar2" not in fieldnames:
            fieldnames.append("num_motifs_farfar2")  # add if missing

        rows = []
        updated, missing = 0, 0
        for row in reader:
            fname = row.get("file_name", "").strip()
            val = farfar_counts.get(fname, 0)
            row["num_motifs_farfar2"] = str(val)
            rows.append(row)
            if fname in farfar_counts:
                updated += 1
            else:
                missing += 1

    out_csv.parent.mkdir(parents=True, exist_ok=True)
    with out_csv.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)

    print(f"[OK] Wrote updated CSV to: {out_csv}")
    print(f"[STATS] Updated rows: {updated}; No counts found (set 0): {missing}")

if __name__ == "__main__":
    main()
