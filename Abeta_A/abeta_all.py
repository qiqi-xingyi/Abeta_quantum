# --*-- conding:utf-8 --*--
# @time:5/5/25 16:58
# @Author : Yuqi Zhang
# @Email : yzhan135@kent.edu
# @File:abeta_all.py

#!/usr/bin/env python3
import os, re, glob

# ---- user-defined: full Abeta-42 one-letter sequence ----
sequence = "DAEFRHDSGYEVHHQKLVFFAEDVGSNKGAIIGLMVGGVVIA"
n_residues = len(sequence)  # should be 42

# data structures for accumulating coords
# coords[i] will be a list of [x,y,z] for residue index i (1-based)
coords = {i: [] for i in range(1, n_residues+1)}

if __name__ == '__main__':


    # find all windowed XYZ files
    pattern = "Abeta_A_window_*.xyz"
    for filepath in glob.glob(pattern):
        # extract start index from filename
        m = re.search(r"window_(\d+)\.xyz$", filepath)
        if not m:
            continue
        start_idx = int(m.group(1))  # 1-based

        with open(filepath) as f:
            lines = [l.strip() for l in f if l.strip()]
        # skip the header line (atom count) and possible comment line
        data_lines = []
        # guess: first line is count (integer), rest are data
        if re.fullmatch(r"\d+", lines[0]):
            data_lines = lines[1:]
        else:
            data_lines = lines

        # parse each atom/residue line
        for local_idx, line in enumerate(data_lines, start=1):
            parts = line.split()
            if len(parts) < 4:
                continue
            res_one_letter = parts[0]
            x, y, z = map(float, parts[1:4])
            global_idx = start_idx + local_idx - 1
            if not (1 <= global_idx <= n_residues):
                # out of bounds: skip or warn
                continue
            # sanity check: sequence match
            expected = sequence[global_idx-1]
            if res_one_letter != expected:
                print(f"Warning: file {filepath} line {local_idx+1} residue code "
                      f"{res_one_letter} ≠ expected {expected} at position {global_idx}")
            coords[global_idx].append((x, y, z))

    # compute averaged coordinates
    averaged = {}
    for i in range(1, n_residues+1):
        pts = coords[i]
        if not pts:
            raise ValueError(f"No coordinates found for residue #{i} ({sequence[i-1]})")
        # average each dimension
        xs, ys, zs = zip(*pts)
        averaged[i] = (sum(xs)/len(xs), sum(ys)/len(ys), sum(zs)/len(zs))

    # write out merged XYZ
    outname = "merged_abeta42.xyz"
    with open(outname, "w") as out:
        out.write(f"{n_residues}\n")
        out.write("Merged Abeta-42 from sliding windows (averaged overlaps)\n")
        for i in range(1, n_residues+1):
            aa = sequence[i-1]
            x, y, z = averaged[i]
            out.write(f"{aa} {x:.6f} {y:.6f} {z:.6f}\n")

    print(f"Done → {outname}")
