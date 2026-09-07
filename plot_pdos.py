import pandas as pd
import matplotlib.pyplot as plt

CSV_PATH = "out/projected_dos.csv"
OUT_PATH = "out/projected_dos.png"

# ---- PRB / APS style ----
plt.rcParams.update({
    "font.family": "serif",
    "font.serif": ["Times New Roman", "Times", "DejaVu Serif"],
    "mathtext.fontset": "stix",
    "font.size": 9,
    "axes.labelsize": 9,
    "axes.titlesize": 9,
    "xtick.labelsize": 8,
    "ytick.labelsize": 8,
    "legend.fontsize": 8,
    "axes.linewidth": 0.8,
    "lines.linewidth": 1.0,
    "xtick.direction": "in",
    "ytick.direction": "in",
    "xtick.top": True,
    "ytick.right": True,
    "xtick.minor.visible": True,
    "ytick.minor.visible": True,
    "xtick.major.size": 4,
    "ytick.major.size": 4,
    "xtick.minor.size": 2,
    "ytick.minor.size": 2,
    "savefig.dpi": 300,
})

df = pd.read_csv(CSV_PATH)
E  = df["energy"]
mu = df["mu"].iloc[0]

# Basis ordering: layer x spin x (yz, xz, xy)
orbitals = ["xy", "xz", "yz"]
orb_colors = {"yz": "tab:green", "xz": "tab:purple", "xy": "tab:blue"}
lstyle = {"yz": "dashed", "xz": "solid", "xy": "solid"}

# PRB double-column width is ~7.0 in; use a modest height for a 2-panel row
fig, axes = plt.subplots(1, 2, figsize=(7.0, 3.0), sharey=True, sharex=True)

for ax, layer in zip(axes, ("L1", "L2")):
    for orb in orbitals:
        up = df[f"{layer}_up_{orb}"]
        dn = df[f"{layer}_dn_{orb}"]
        ax.plot(E,  up, color=orb_colors[orb], linewidth=1.0, label=orb, linestyle=lstyle[orb])
        ax.plot(E, -dn, color=orb_colors[orb], linewidth=1.0, linestyle=lstyle[orb])

    ax.axvline(mu, color="black", linewidth=0.7, linestyle="--", label=r"$\mu$")
    ax.axhline(0.0, color="black", linewidth=0.6)
    ax.set_xlabel("Energy (eV)")
    ax.set_title(f"Layer {layer[-1]}", fontsize=9, loc="left")
    # APS style: full box, ticks on all sides (set above via rcParams);
    # do NOT hide spines here.

axes[0].set_ylabel("DOS (arb. units)")
axes[0].legend(loc="upper right", frameon=False, ncol=2, handlelength=1.5,
               columnspacing=1.0, borderaxespad=0.3)

fig.text(0.02, 0.93, "(a)", fontsize=9, fontweight="bold")
fig.text(0.52, 0.93, "(b)", fontsize=9, fontweight="bold")

plt.tight_layout()
plt.savefig(OUT_PATH)
plt.savefig(OUT_PATH.replace(".png", ".pdf"))
print(f"Saved: {OUT_PATH}")