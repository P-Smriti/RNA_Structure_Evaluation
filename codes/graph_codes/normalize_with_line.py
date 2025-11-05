import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

def plot_motif_histogram(csv_path, y_column, bin_width=10, save_path=None):
    """
    Plot a normalized histogram (bars) with a line overlay.
    X-axis: sequence length (binned)
    Y-axis: average motifs per structure in each bin

    Parameters
    ----------
    csv_path : str
        Path to CSV file.
    y_column : str
        One of: 'num_motifs_farfar2', 'num_motifs_rhofold', 'num_motifs_alphafold3'
    bin_width : int
        Width of bins (default = 10)
    save_path : str or None
        Optional path to save the figure.
    """

    allowed = ["num_motifs_farfar2", "num_motifs_rhofold", "num_motifs_alphafold3"]
    if y_column not in allowed:
        raise ValueError("y_column must be one of: " + ", ".join(allowed))

    # Read CSV
    df = pd.read_csv(csv_path)
    df["sequence_length"] = pd.to_numeric(df["sequence_length"], errors="coerce")
    df[y_column] = pd.to_numeric(df[y_column], errors="coerce")
    df = df.dropna(subset=["sequence_length", y_column])

    if df.empty:
        raise ValueError("No valid data in file.")

    # Define bins
    min_len = int(np.floor(df["sequence_length"].min()))
    max_len = int(np.ceil(df["sequence_length"].max()))
    start = (min_len // bin_width) * bin_width
    end_edge = ((max_len // bin_width) + 1) * bin_width
    bins = np.arange(start, end_edge + bin_width, bin_width)

    # Bin the data
    binned = pd.cut(df["sequence_length"], bins=bins, right=False)

    # Compute per-bin values
    counts = df.groupby(binned)[y_column].count()
    sums = df.groupby(binned)[y_column].sum()
    normalized = (sums / counts).fillna(0)

    # Ensure all bins appear
    intervals = pd.IntervalIndex.from_breaks(bins, closed="left")
    normalized = normalized.reindex(intervals, fill_value=0)

    # Labels and midpoints
    labels = [f"{int(iv.left)}-{int(iv.right - 1)}" for iv in normalized.index]
    midpoints = [iv.mid for iv in normalized.index]

    # Plot
    plt.figure(figsize=(10, 5))
    plt.bar(labels, normalized.values, color="skyblue", alpha=0.7, label="Avg motifs/bin")
    plt.plot(range(len(normalized)), normalized.values, color="darkblue", marker="o", linewidth=2, label="Trend line")

    plt.xlabel("Sequence length (bins of %d)" % bin_width)
    plt.ylabel("Average %s per structure" % y_column.replace("_", " "))
    plt.title("Normalized %s vs sequence length" % y_column.replace("_", " "))
    plt.xticks(rotation=45, ha="right")
    plt.grid(axis="y", linestyle="--", alpha=0.4)
    plt.legend()
    plt.tight_layout()

    if save_path:
        plt.savefig(save_path, dpi=150, bbox_inches="tight")

    plt.show()

    # Return data used for plotting
    return pd.DataFrame({
        "length_bin": labels,
        "midpoint": midpoints,
        "avg_motifs": normalized.values
    })

plot_column_name = "num_motifs_farfar2"
df = plot_motif_histogram(
    "../../data/fasta_mapping_with_length_updated.csv",
    plot_column_name,
    bin_width=10,
    save_path=f"normalized_{plot_column_name[4:]}_hist_line.png"
)
print(df)

