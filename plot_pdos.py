import pandas as pd
import matplotlib.pyplot as plt

CSV_PATH = "out/projected_dos.csv"
OUT_PATH = "out/projected_dos.png"

df = pd.read_csv(CSV_PATH)
E  = df["energy"]
mu = df["mu"].iloc[0]

# Basis ordering: layer x spin x (yz, xz, xy)
orbitals = ["xy", "xz", "yz"]
orb_colors = {"yz": "tab:green", "xz": "tab:purple", "xy": "tab:blue"}
lstyle = {"yz": "dashed", "xz": "solid", "xy": "solid"}

# Two layers, each with spin-up (plotted positive) / spin-down (plotted negative)
fig, axes = plt.subplots(1, 2, figsize=(11, 5), sharey=True, sharex=True)

for ax, layer in zip(axes, ("L1", "L2")):
    for orb in orbitals:
        up = df[f"{layer}_up_{orb}"]
        dn = df[f"{layer}_dn_{orb}"]
        ax.plot(E,  up, color=orb_colors[orb], linewidth=1.0, label=orb, linestyle=lstyle[orb])
        ax.plot(E, -dn, color=orb_colors[orb], linewidth=1.0, linestyle=lstyle[orb])

    ax.axvline(mu, color="black", linewidth=0.8, linestyle="--", label=r"$\mu$")
    ax.axhline(0.0, color="black", linewidth=0.6)
    ax.set_xlabel("Energy", fontsize=11)
    ax.set_title(f"Projected DOS - Layer {layer[-1]}  (up > 0 / down < 0)", fontsize=11)
    ax.legend(fontsize=9, frameon=False, ncol=2)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

axes[0].set_ylabel("DOS  (spin up > 0,  spin down < 0)", fontsize=11)

plt.tight_layout()
plt.savefig(OUT_PATH, dpi=150)
plt.show()
print(f"Saved: {OUT_PATH}")
