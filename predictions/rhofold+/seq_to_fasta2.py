input_file = "human_remaining_final.txt"
output_file = "human_multi_fasta.fasta"

fasta_entries = []
current_family = None

with open(input_file, "r") as f:
    for line in f:
        line = line.strip()
        if not line:
            continue

        if line.startswith("Family:"):
            current_family = line.split()[1].strip(",")
            continue

        parts = line.split()
        if len(parts) == 2 and current_family:
            seq_id, seq = parts
            seq = seq.replace("-", "")   # remove gaps
            fasta_entries.append(f">{seq_id}|{current_family}\n{seq}")

with open(output_file, "w") as f:
    f.write("\n".join(fasta_entries))

print(f"Multi-FASTA written to {output_file}")
