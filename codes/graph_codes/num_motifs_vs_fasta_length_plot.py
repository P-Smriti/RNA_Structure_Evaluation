import pandas as pd
import matplotlib.pyplot as plt

def plot_motif_histogram(file_path, motif_column):
    """
    Reads a CSV file, plots a histogram of sequence length vs. number of motifs.

    Args:
        file_path (str): The path to the CSV file.
        motif_column (str): The name of the column containing motif counts
                            ("num_motifs_farfar2", "num_motifs_rhofold", or "num_motifs_alphafold3").
    """
    try:
        df = pd.read_csv(file_path)
    except FileNotFoundError:
        print(f"Error: File not found at {file_path}")
        return
    except Exception as e:
        print(f"Error reading CSV file: {e}")
        return

    if motif_column not in ["num_motifs_farfar2", "num_motifs_rhofold", "num_motifs_alphafold3"]:
        print(f"Error: Invalid motif column '{motif_column}'. Please choose from 'num_motifs_farfar2', 'num_motifs_rhofold', or 'num_motifs_alphafold3'.")
        return

    if "sequence_length" not in df.columns:
        print("Error: 'sequence_length' column not found in the CSV file.")
        return

    if motif_column not in df.columns:
        print(f"Error: '{motif_column}' column not found in the CSV file.")
        return

    # Create bins of 10 for sequence length
    bins = range(0, df['sequence_length'].max() + 10, 10)

    plt.figure(figsize=(12, 6))
    plt.hist(df['sequence_length'], weights=df[motif_column], bins=bins, edgecolor='black')
    plt.xlabel("Sequence Length")
    plt.ylabel(f"Number of Motifs ({motif_column})")
    plt.title(f"Histogram of Sequence Length vs. {motif_column}")
    plt.grid(axis='y', alpha=0.75)
    plt.show()
    plt.savefig("/home/s081p868/scratch/RNA_Structure_Evaluation/data/sequence_length_vs_motif.png", dpi=300, bbox_inches='tight')

#plt.show()
