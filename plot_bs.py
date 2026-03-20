import os
import glob
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

BS_DIR  = "out/bs_plots"
DOS_DIR = "out/dos"
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
    mask   = df["label"].notna() & (df["label"].astype(str).str.strip() != "")
    ticks  = df.loc[mask, "path_coord"].tolist()
    labels = df.loc[mask, "label"].tolist()
    labels = [r"$\Gamma$" if l in ("G", "Gamma") else l for l in labels]

    # --- Chemical potential ---
    mu = df["mu"].iloc[0]

    # --- Extract U value and match DOS file ---
    basename = os.path.splitext(os.path.basename(csv_path))[0]
    u_str    = basename.split("_U")[-1]
    title    = f"U = {u_str}"

    dos_filename = basename.replace("bs_", "dos_") + ".csv"
    dos_path     = os.path.join(DOS_DIR, dos_filename)

    # --- Plot ---
    has_dos = os.path.exists(dos_path)
    fig, axes = plt.subplots(1, 2 if has_dos else 1,
                             figsize=(10, 5) if has_dos else (7, 5),
                             sharey=True)
    ax_bs  = axes[0] if has_dos else axes
    ax_dos = axes[1] if has_dos else None

    # --- Band structure ---
    band_cols = [c for c in df.columns if c.startswith("band_")]
    for col in band_cols:
        ax_bs.plot(df["path_coord"], df[col], color="steelblue", linewidth=0.8)

    ax_bs.axhline(mu, color="red", linewidth=0.8, linestyle="--", label=r"$\mu$")
    ax_bs.legend(fontsize=10, frameon=False)

    for t in ticks:
        ax_bs.axvline(t, color="black", linewidth=0.6, linestyle="-")

    ax_bs.set_xticks(ticks)
    ax_bs.set_xticklabels(labels, fontsize=12)
    ax_bs.tick_params(axis="x", which="major", direction="in", length=6, width=0.8)
    ax_bs.set_xlim(df["path_coord"].iloc[0], df["path_coord"].iloc[-1])
    ax_bs.set_ylabel("Energy", fontsize=11)
    ax_bs.set_title(f"Band Structure   {title}", fontsize=11)
    ax_bs.spines["top"].set_visible(False)
    ax_bs.spines["right"].set_visible(False)

    # --- DOS ---
    if has_dos:
        dos_df = pd.read_csv(dos_path)
        ax_dos.plot(dos_df["dos"], dos_df["energy"], color="steelblue", linewidth=0.8)
        ax_dos.axhline(mu, color="red", linewidth=0.8, linestyle="--")
        ax_dos.set_xlabel("DOS", fontsize=11)
        ax_dos.set_title(f"DOS   {title}", fontsize=11)
        ax_dos.spines["top"].set_visible(False)
        ax_dos.spines["right"].set_visible(False)
        ax_dos.tick_params(axis="x", which="major", direction="in", length=6, width=0.8)

    plt.tight_layout()

    out_path = os.path.join(OUT_DIR, basename + ".png")
    fig.savefig(out_path, dpi=150)
    plt.close(fig)
    print(f"  Saved: {out_path}")

print("Done.")