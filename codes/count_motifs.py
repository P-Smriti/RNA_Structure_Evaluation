#!/usr/bin/env python3
"""
Update mapping CSV with motif hit totals from one or more tools AND
write per-motif breakdown CSVs.

Scans:
  <ROOT>/URS*/Res_motifs/{c-loop_consensus,e-loop_consensus,k-turn_consensus,
                         reverse-kturn_consensus,sarcin-ricin_consensus}/result.log

Counts rows as non-empty, non-header (header may begin with '#' or not).

Adds/updates these columns when the corresponding root is provided:
  - num_motifs_farfar2
  - num_motifs_rhofold
  - num_motifs_alphafold3

Also writes per-tool CSVs with columns:
  file_name,
  num_motifs_c-loop_consensus,
  num_motifs_e-loop_consensus,
  num_motifs_k-turn_consensus,
  num_motifs_reverse-kturn_consensus,
  num_motifs_sarcin-ricin_consensus,
  <tool-total-column>
"""

import argparse
import csv
from pathlib import Path
from collections import defaultdict

EXPECTED_MOTIFS = [
    "c-loop_consensus",
    "e-loop_consensus",
    "k-turn_consensus",
    "reverse-kturn_consensus",
    "sarcin-ricin_consensus",
]

def count_results(log_path: Path) -> int:
    """Count non-header, non-blank lines in result.log (header can be '#fragment_ID...' or 'fragment_ID...')."""
    if not log_path.is_file():
        return 0
    n = 0
    first_nonempty_seen = False
    with log_path.open("r", encoding="utf-8", errors="ignore") as f:
        for line in f:
            s = line.strip()
            if not s:
                continue
            if not first_nonempty_seen:
                first_nonempty_seen = True
                head = s.lstrip("#").strip().lower()
                if head.startswith("fragment_id"):
                    # header row
                    continue
            if s.startswith("#"):
                continue
            n += 1
    return n

def collect_counts_and_breakdown(root: Path):
    """
    Return:
      totals: {URS -> total_count}
      per_motif: {URS -> {motif -> count}}  (for EXPECTED_MOTIFS only)
    """
    totals = defaultdict(int)
    per_motif = defaultdict(lambda: {m: 0 for m in EXPECTED_MOTIFS})
    if not root or not root.is_dir():
        return {}, {}
    for urs in sorted(p for p in root.glob("URS*") if p.is_dir()):
        mroot = urs / "Res_motifs"
        if not mroot.is_dir():
            continue
        total = 0
        # ensure we only count the five specific motifs; missing dirs -> 0
        for motif in EXPECTED_MOTIFS:
            d = mroot / motif
            c = count_results(d / "result.log") if d.is_dir() else 0
            per_motif[urs.name][motif] = c
            total += c
        totals[urs.name] = total
    return dict(totals), dict(per_motif)

def write_per_motif_csv(out_path: Path, per_motif: dict, total_col: str):
    """
    Write CSV in the exact format requested:
    file_name,num_motifs_c-loop_consensus,...,num_motifs_sarcin-ricin_consensus,<total_col>
    """
    if not per_motif:
        return
    header = ["file_name"] + [f"num_motifs_{m}" for m in EXPECTED_MOTIFS] + [total_col]
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with out_path.open("w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=header)
        w.writeheader()
        for urs in sorted(per_motif.keys()):
            row = {"file_name": urs}
            total = 0
            for m in EXPECTED_MOTIFS:
                v = int(per_motif[urs].get(m, 0))
                row[f"num_motifs_{m}"] = v
                total += v
            row[total_col] = total
            w.writerow(row)

def main():
    ap = argparse.ArgumentParser(description="Update mapping CSV with motif totals + write per-motif CSVs")
    ap.add_argument("--mapping-csv", required=True, help="Existing mapping CSV to update")
    ap.add_argument("--out-csv",     required=True, help="Where to write the updated CSV")

    ap.add_argument("--farfar-root",    help="Path to farfar_pdb directory")
    ap.add_argument("--rhofold-root",   help="Path to rhofold_pdb directory")
    ap.add_argument("--alphafold-root", help="Path to alphafold_pdb directory")
    args = ap.parse_args()

    mapping_csv = Path(args.mapping_csv).expanduser().resolve()
    out_csv     = Path(args.out_csv).expanduser().resolve()
    if not mapping_csv.is_file():
        raise SystemExit(f"[ERROR] Mapping CSV not found: {mapping_csv}")

    # Collect totals + per-motif breakdowns (any args may be None)
    farfar_root    = Path(args.farfar_root).expanduser().resolve()    if args.farfar_root    else None
    rhofold_root   = Path(args.rhofold_root).expanduser().resolve()   if args.rhofold_root   else None
    alphafold_root = Path(args.alphafold_root).expanduser().resolve() if args.alphafold_root else None

    f_tot, f_break = collect_counts_and_breakdown(farfar_root)    if farfar_root    else ({}, {})
    r_tot, r_break = collect_counts_and_breakdown(rhofold_root)   if rhofold_root   else ({}, {})
    a_tot, a_break = collect_counts_and_breakdown(alphafold_root) if alphafold_root else ({}, {})

    print(f"[INFO] FARFAR2 URS counted:   {len(f_tot)}")
    print(f"[INFO] RhoFold  URS counted:  {len(r_tot)}")
    print(f"[INFO] AlphaFold3 URS counted:{len(a_tot)}")

    # Write per-motif CSVs at the requested default locations
    if farfar_root:
        write_per_motif_csv(farfar_root.parent / "farfar2_motif_counts_1.csv", f_break, "num_motifs_farfar2")
    if rhofold_root:
        write_per_motif_csv(rhofold_root.parent / "rhofold_motif_counts_1.csv", r_break, "num_motifs_rhofold")
    if alphafold_root:
        write_per_motif_csv(alphafold_root.parent / "alphafold3_motif_counts_1.csv", a_break, "num_motifs_alphafold3")

    # Update totals in mapping CSV
    with mapping_csv.open("r", newline="", encoding="utf-8") as f:
        reader = csv.DictReader(f)
        fieldnames = list(reader.fieldnames or [])

        if f_tot and "num_motifs_farfar2" not in fieldnames:
            fieldnames.append("num_motifs_farfar2")
        if r_tot and "num_motifs_rhofold" not in fieldnames:
            fieldnames.append("num_motifs_rhofold")
        if a_tot and "num_motifs_alphafold3" not in fieldnames:
            fieldnames.append("num_motifs_alphafold3")

        rows = []
        for row in reader:
            name = (row.get("file_name") or "").strip()
            if f_tot:
                row["num_motifs_farfar2"]  = str(f_tot.get(name, 0))
            if r_tot:
                row["num_motifs_rhofold"]  = str(r_tot.get(name, 0))
            if a_tot:
                row["num_motifs_alphafold3"] = str(a_tot.get(name, 0))
            rows.append(row)

    out_csv.parent.mkdir(parents=True, exist_ok=True)
    with out_csv.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)

    print(f"[OK] Wrote updated mapping to: {out_csv}")
    if farfar_root:
        print(f"[OK] Wrote FARFAR2 per-motif CSV to: {farfar_root.parent / 'farfar2_motif_counts_1.csv'}")
    if rhofold_root:
        print(f"[OK] Wrote RhoFold per-motif CSV to: {rhofold_root.parent / 'rhofold_motif_counts_1.csv'}")
    if alphafold_root:
        print(f"[OK] Wrote AlphaFold3 per-motif CSV to: {alphafold_root.parent / 'alphafold3_motif_counts_1.csv'}")

if __name__ == "__main__":
    main()
