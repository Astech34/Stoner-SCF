import pandas as pd
import matplotlib.pyplot as plt

# --- Load data ---
df = pd.read_csv("build/out/stoner_U_sweep.csv")

# --- Plot ---
fig, ax = plt.subplots(figsize=(8, 5))

ax.plot(
    df["U"], df["S_final"],
    color="steelblue",
    linewidth=1.5,
    zorder=1
)
ax.scatter(
    df["U"], df["S_final"],
    color="steelblue",
    edgecolors="white",
    s=60,
    linewidths=1.2,
    zorder=2
)

# --- Labels ---
ax.set_xlabel("Hubbard $U$ / $t_1$", fontsize=13)
ax.set_ylabel("Magnetisation $S$", fontsize=13)
ax.set_title("Stoner SCF: Magnetisation vs Hubbard $U$", fontsize=14)

# --- Formatting ---
ax.set_xlim(df["U"].min() - 0.3, df["U"].max() + 0.3)
ax.set_ylim(-0.02, df["S_final"].max() * 1.1 + 0.01)
ax.axhline(0, color="black", linewidth=0.8, linestyle="--", alpha=0.4)
ax.grid(True, alpha=0.3, linestyle="--")
ax.tick_params(labelsize=11)

plt.tight_layout()
plt.savefig("build/out/stoner_U_sweep.png", dpi=150)
plt.show()

print("Plot saved to out/stoner_U_sweep.png")