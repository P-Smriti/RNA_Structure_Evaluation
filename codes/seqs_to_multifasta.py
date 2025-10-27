#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import re
from textwrap import wrap
from pathlib import Path

# ---- Hardcoded paths (as requested) ----
nfile   = "/home/s081p868/scratch/RNA_Structure_Evaluation/data/human_seqs_non3d_rfams.txt"
outfile = "/home/s081p868/scratch/RNA_Structure_Evaluation/data/human_seqs_non3d_rfams.fa"

WRAP_LEN = None  # set to an int (e.g., 80) if you want FASTA line-wrapping

def sanitize_id(s: str) -> str:
    """Replace any non [A-Za-z0-9_.-] with underscore."""
    return re.sub(r"[^A-Za-z0-9_.-]", "_", s)

def _bump_suffix(base: str, taken: set) -> str:
    i = 1
    while True:
        cand = f"{base}_{i}"
        if cand not in taken:
            return cand
        i += 1

def parse_to_fasta(infile: str, outfile: str, wrap_len=None):
    infile = Path(infile)
    outfile = Path(outfile)

    taken = set()  # headers already used (to avoid accidental duplicates)

    with infile.open() as fin, outfile.open("w") as fout:
        current_family = None

        for raw_line in fin:
            line = raw_line.strip()
            if not line:
                continue

            # Family line, e.g., "Family: RF00246"
            if line.startswith("Family:"):
                m = re.match(r"Family:\s*(\S+)", line)
                if m:
                    current_family = m.group(1).rstrip(",")  # e.g., RF00246
                continue

            # Sequence line starts with URS... then the sequence
            if line.startswith("URS"):
                parts = line.split()
                if len(parts) < 2:
                    continue

                original_id = parts[0]               # e.g., URS00000A7BEB_9606/9-93
                seq = parts[1].replace("-", "")      # remove gaps
                family = current_family or ""         # e.g., RF00246

                # Build header to match PDB basename pattern:
                # URS00000A7BEB_9606_9-93_RF00246
                core = sanitize_id(original_id)       # converts "/" -> "_", spaces -> "_", etc.
                header_id = f"{core}_{family}" if family else core

                # Ensure uniqueness inside this run
                if header_id in taken:
                    header_id = _bump_suffix(header_id, taken)
                taken.add(header_id)

                # Write FASTA
                fout.write(f">{header_id}\n")
                if wrap_len:
                    for chunk in wrap(seq, wrap_len):
                        fout.write(chunk + "\n")
                else:
                    fout.write(seq + "\n")

if __name__ == "__main__":
    parse_to_fasta(nfile, outfile, wrap_len=WRAP_LEN)
    print(f"MultiFASTA written to {outfile}")
