import pandas as pd
import matplotlib.pyplot as plt

CSV_PATH = "out/band_structure.csv"

df = pd.read_csv(CSV_PATH)

mask   = df["label"].notna() & (df["label"].astype(str).str.strip() != "")
ticks  = df.loc[mask, "path_coord"].tolist()
labels = df.loc[mask, "label"].tolist()
labels = [r"$\Gamma$" if l in ("G", "Gamma") else l for l in labels]

mu = df["mu"].iloc[0]

fig, ax = plt.subplots(figsize=(7, 5))

band_cols = [c for c in df.columns if c.startswith("band_")]
for col in band_cols:
    ax.plot(df["path_coord"], df[col], color="steelblue", linewidth=0.8)

ax.axhline(mu, color="red", linewidth=0.8, linestyle="--", label=r"$\mu$")
ax.legend(fontsize=10, frameon=False)

for t in ticks:
    ax.axvline(t, color="black", linewidth=0.6)

ax.set_xticks(ticks)
ax.set_xticklabels(labels, fontsize=12)
ax.tick_params(axis="x", which="major", direction="in", length=6, width=0.8)
ax.set_xlim(df["path_coord"].iloc[0], df["path_coord"].iloc[-1])
ax.set_ylabel("Energy", fontsize=11)
ax.set_title("Band Structure", fontsize=11)
ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)

plt.tight_layout()
plt.savefig("out/band_structure.png", dpi=150)
plt.show()
print("Saved: out/band_structure.png")
