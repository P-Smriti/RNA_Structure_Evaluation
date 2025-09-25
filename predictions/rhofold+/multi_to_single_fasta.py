from pathlib import Path
import re

in_path = Path("human_sequences.fasta")  # change to your multi-FASTA
outdir = Path("split_fa")
outdir.mkdir(parents=True, exist_ok=True)

with in_path.open() as fh:
    out = None
    for line in fh:
        if line.startswith(">"):
            if out and not out.closed:
                out.close()
            header = line[1:].strip()             # drop ">"
            safe = re.sub(r"[\/|]", "_", header)  # replace / and |
            out = (outdir / f"{safe}.fasta").open("w")
            out.write(line)                       # write header as-is
        else:
            if out:
                out.write(line)
    if out and not out.closed:
        out.close()

