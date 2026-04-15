import pandas as pd
import matplotlib.pyplot as plt

# --- Load parameters from comment lines ---
params = {}
with open("out/mca_lam_sweep.csv") as f:
    for line in f:
        if not line.startswith("#"):
            break
        key, val = line[1:].split("=")
        params[key.strip()] = float(val.strip())

param_text = "\n".join(f"{k} = {v:g}" for k, v in params.items())

# --- Load data ---
df = pd.read_csv("build/out/mca_lam_sweep.csv",
                 comment="#")

fig, axes = plt.subplots(1, 2, figsize=(12, 5))

# --- Left: E_MCA vs lambda ---
ax = axes[0]
ax.plot(df["lam"], df["E_MCA"].abs(),
        color="steelblue", linewidth=1.5, zorder=1)
ax.scatter(df["lam"], df["E_MCA"].abs(),
           color="steelblue", edgecolors="white",
           s=60, linewidths=1.2, zorder=2)
ax.set_xlabel(r"SOC strength $\lambda$ / $t_1$", fontsize=13)
ax.set_ylabel(r"$|E_{[110]} - E_{[001]}|$ / $t_1$", fontsize=13)
ax.set_title("MCA Energy vs SOC", fontsize=14)
ax.grid(True, alpha=0.3, linestyle="--")
ax.tick_params(labelsize=11)

# --- Right: S vs lambda for both directions ---
ax = axes[1]
ax.plot(df["lam"], df["S_110"],
        color="tomato", linewidth=1.5, label=r"$\hat{n} = [110]$", zorder=1)
ax.scatter(df["lam"], df["S_110"],
           color="tomato", edgecolors="white",
           s=60, linewidths=1.2, zorder=2)
ax.plot(df["lam"], df["S_001"],
        color="steelblue", linewidth=1.5, label=r"$\hat{n} = [001]$", zorder=1)
ax.scatter(df["lam"], df["S_001"],
           color="steelblue", edgecolors="white",
           s=60, linewidths=1.2, zorder=2)
ax.set_xlabel(r"SOC strength $\lambda$ / $t_1$", fontsize=13)
ax.set_ylabel(r"Magnetisation $S$", fontsize=13)
ax.set_title("Magnetisation vs SOC", fontsize=14)
ax.legend(fontsize=11)
ax.grid(True, alpha=0.3, linestyle="--")
ax.tick_params(labelsize=11)

fig.text(0.5, -0.02, param_text, ha="center", va="top", fontsize=10,
         family="monospace",
         bbox=dict(boxstyle="round", facecolor="lightyellow", alpha=0.6))

plt.tight_layout()
plt.savefig("out/mca_lam_sweep.png",
            dpi=150, bbox_inches="tight")
plt.show()

print("Plot saved to out/mca_lam_sweep.png")
