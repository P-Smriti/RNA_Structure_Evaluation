#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import re
from textwrap import wrap

infile = "/home/s081p868/scratch/RNA_Structure_Evaluation/data/human_seqs_non3d_rfams.txt"
outfile = "/home/s081p868/scratch/RNA_Structure_Evaluation/data/human_seqs_non3d_rfams.fasta"

def parse_to_fasta(infile, outfile, wrap_len=None):
    with open(infile) as f, open(outfile, "w") as out:
        current_family = None
        for line in f:
            line = line.strip()
            if not line:
                continue
            # detect family line
            if line.startswith("Family:"):
                m = re.match(r"Family:\s*(\S+)", line)
                if m:
                    current_family = m.group(1).rstrip(",")  # remove trailing comma if present
            # detect sequence line
            elif line.startswith("URS"):
                parts = line.split()
                if len(parts) >= 2 and current_family:
                    seq_id = parts[0]
                    seq = parts[1].replace("-", "")  # remove gaps
                    header = f">{seq_id} {current_family}"
                    out.write(header + "\n")
                    if wrap_len:
                        for chunk in wrap(seq, wrap_len):
                            out.write(chunk + "\n")
                    else:
                        out.write(seq + "\n")

if __name__ == "__main__":
    parse_to_fasta(infile, outfile, wrap_len=None)
    print(f"MultiFASTA written to {outfile}")