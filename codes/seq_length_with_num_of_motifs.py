

import os
import pandas as pd

input_dir = r"../data/split_fa"
output_excel = r"results.csv"

data = []

for filename in sorted(os.listdir(input_dir)):
    if filename.endswith(".fasta"):
        filepath = os.path.join(input_dir, filename)
        with open(filepath, "r") as f:
            lines = f.read().splitlines()
            if not lines:
                continue

            header = lines[0].strip()
            if header.startswith(">"):
                header = header[1:]
            clean_name = header.replace("/", "_").replace("|", "_")

            seq = "".join(lines[1:]).replace(" ", "").replace("\n", "")
            seq_length = len(seq)

            data.append({
                #"fasta_file": filename,
                "file_name": clean_name,
                "sequence_length": seq_length,
                "num_motifs_farfar2": 1,     # placeholder
                "num_motifs_rhofold": 1,     # placeholder
                "num_motifs_alphafold3": 1   # placeholder
            })

df = pd.DataFrame(data)
df.to_csv(output_excel, index=False)

print(f"Excel saved successfully to: {output_excel}")
