import os
import glob
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

BS_DIR  = "out/bs_plots"
OUT_DIR = "out/bs_plots/png"
os.makedirs(OUT_DIR, exist_ok=True)

csv_files = sorted(glob.glob(os.path.join(BS_DIR, "*.csv")))
if not csv_files:
    print(f"No CSV files found in {BS_DIR}")
    exit(1)

print(f"Found {len(csv_files)} band structure files")

for csv_path in csv_files:
    df = pd.read_csv(csv_path)

    # --- Extract high-symmetry tick positions and labels ---
    # Use BOTH conditions: label is non-empty AND path_coord is non-zero
    # (Gamma at coord=0 is the exception — keep it via label check alone)
    mask = df["label"].notna() & (df["label"].astype(str).str.strip() != "")
    ticks  = df.loc[mask, "path_coord"].tolist()
    labels = df.loc[mask, "label"].tolist()

    # Replace "G" or "Gamma" with proper Gamma symbol
    labels = [r"$\Gamma$" if l in ("G", "Gamma") else l for l in labels]

    # --- Extract U value from filename for title ---
    basename = os.path.splitext(os.path.basename(csv_path))[0]
    u_str    = basename.split("_U")[-1]
    title    = f"Band Structure   U = {u_str}"

    # --- Plot ---
    fig, ax = plt.subplots(figsize=(7, 5))

    band_cols = [c for c in df.columns if c.startswith("band_")]
    for col in band_cols:
        ax.plot(df["path_coord"], df[col], color="steelblue", linewidth=0.8)

    # High-symmetry vertical lines
    for t in ticks:
        ax.axvline(t, color="black", linewidth=0.6, linestyle="-")

    # Tick marks only at high-symmetry points, no minor ticks elsewhere
    ax.set_xticks(ticks)
    ax.set_xticklabels(labels, fontsize=12)
    ax.tick_params(axis="x", which="major", direction="in", length=6, width=0.8)

    ax.set_xlim(df["path_coord"].iloc[0], df["path_coord"].iloc[-1])
    ax.set_ylabel("Energy", fontsize=11)
    ax.set_title(title, fontsize=11)

    # Clean up top/right spines (optional, common for band structures)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    plt.tight_layout()

    out_path = os.path.join(OUT_DIR, basename + ".png")
    fig.savefig(out_path, dpi=150)
    plt.close(fig)
    print(f"  Saved: {out_path}")

print("Done.")