import sys
import re

import numpy as np
import subprocess
import shutil
import os
from pathlib import Path

# Matches: <indent><key><spaces>=<spaces><value><trailing comment>
LINE_RE = re.compile(
    r'^(?P<indent>\s*)'
    r'(?P<key>[A-Za-z_][A-Za-z0-9_]*)'
    r'(?P<eq>\s*=\s*)'
    r'(?P<value>[^#\n]*?)'
    r'(?P<trail>\s*(#.*)?)\s*$'
)


def parse_updates(args):
    """Parse ['U=2.5', 'J=0.4'] -> {'U': '2.5', 'J': '0.4'}"""
    updates = {}
    for arg in args:
        if '=' not in arg:
            print(f"Skipping malformed arg (no '='): {arg}")
            continue
        key, value = arg.split('=', 1)
        updates[key.strip()] = value.strip()
    return updates


def update_file(path: Path, updates: dict):
    lines = path.read_text().splitlines(keepends=True)
    remaining = set(updates.keys())
    out_lines = []

    for line in lines:
        # Skip full-line comments and blank lines untouched
        stripped = line.strip()
        if stripped.startswith('#') or not stripped:
            out_lines.append(line)
            continue

        m = LINE_RE.match(line)
        if not m:
            # Doesn't look like a key=value line; leave as-is
            out_lines.append(line)
            continue

        key = m.group('key')
        if key in updates:
            new_value = updates[key]
            newline = '\n' if line.endswith('\n') else ''
            out_lines.append(
                f"{m.group('indent')}{key}{m.group('eq')}{new_value}{m.group('trail')}{newline}"
            )
            remaining.discard(key)
        else:
            out_lines.append(line)

    path.write_text(''.join(out_lines))

    if remaining:
        print(f"Warning: these keys were not found in {path} and were NOT added: {sorted(remaining)}")

    applied = set(updates.keys()) - remaining
    if applied:
        print(f"Updated: {', '.join(f'{k}={updates[k]}' for k in sorted(applied))}")

def archive_csv(src_csv: Path, dest_path: Path):
    """
    Copy src_csv to dest_dir, renaming it to include the parameter name/value.
    e.g. archive_csv(Path("out/rho.csv"), Path("out/archive"), "lam", 0.05)
         -> out/archive/rho_lam0.05.csv
    """
    #dest_dir.mkdir(parents=True, exist_ok=True)

    stem = src_csv.stem       # "rho"
    suffix = src_csv.suffix   # ".csv"

    if not src_csv.exists():
        raise FileNotFoundError(f"Expected output CSV not found: {src_csv}")

    shutil.move(src_csv, dest_path)
    print(f"Archived: {dest_path}")
    return dest_path

env = os.environ.copy()
env["OMP_NUM_THREADS"] = "32"

command = ["./build/Stoner-SCF"]

if __name__ == "__main__":
    for j in np.linspace(0.1, 1, 10):
        for cf in np.linspace(-0.2, 0.2, 10):
            for ne in [8.0, 9.0, 10.0]:
                update_file(Path("params.in"), {"J": str(j), "delta_cf1": str(cf), "delta_cf2": str(cf), "N_target": str(ne)})

                try:
                    result = subprocess.run(command, env=env, check=True, text=True)
                except subprocess.CalledProcessError as e:
                    print(f"Run failed for j={j}, cf={cf}, ne={ne}: {e}")
                    continue  # or break, depending on what you want

                print("Success!")
                try:
                    archive_csv(Path("out/n_electron_sweep.csv"), Path(f"out/n_sweep_csv/nsweep_j{j}_{cf}_{ne}.csv"))
                except FileNotFoundError as e:
                    print(e)