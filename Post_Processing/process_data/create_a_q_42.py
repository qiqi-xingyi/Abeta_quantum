# --*-- conding:utf-8 --*--
# @time:5/5/25 17:32
# @Author : Yuqi Zhang
# @Email : yzhan135@kent.edu
# @File:create_a_q_42.py

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
stitch_af_quantum_rigid.py

方案：
• Quantum PDB:   1‑7 + 26‑end
• AlphaFold CIF: 8‑25
在 7↔8、25↔26 处用 (N,CA,C) 三个主链原子做刚性超定位，保证连续性。
"""

import numpy as np
from Bio.PDB import (
    MMCIFParser, PDBParser, Superimposer,
    Structure, Model, Chain, Residue, PDBIO
)

# ————— 文件路径 —————
AF_CIF = "fold_42/fold_42_model_0.cif"
QT_PDB = "Abeta_quantum.pdb"
OUT_PDB = "Abeta42.pdb"

# ————— 链与区段 —————
CHAIN = "A"
HELIX_START, HELIX_END = 8, 28          # AF 螺旋 8‑25
BACKBONE = ("N", "CA", "C")             # 用这 3 个原子做对齐

# ------------------------------------------------------------------
def backbone_atoms(res, names=BACKBONE):
    """返回 Residue 中指定名称的 Atom 对象列表（按给定顺序）"""
    return [res[a] for a in names if a in res]

def copy_residue(res):
    """深拷贝 Residue（保留 id/resname/segid）"""
    nr = Residue.Residue(res.id, res.resname, res.segid)
    for atom in res:
        nr.add(atom.copy())
    return nr

def rigid_align_segment(moving_residues, fixed_anchor_res, moving_anchor_res):
    """
    对 moving_residues 做刚性变换，使 moving_anchor_res 的骨架
    精确贴到 fixed_anchor_res。
    """
    fixed_atoms   = backbone_atoms(fixed_anchor_res)
    moving_atoms  = backbone_atoms(moving_anchor_res)

    sup = Superimposer()
    sup.set_atoms(fixed_atoms, moving_atoms)

    # 收集需要整体变换的原子
    all_atoms = [a for res in moving_residues for a in res]
    sup.apply(all_atoms)

    print(f"  align RMSD={sup.rms:.3f} Å  (res {moving_anchor_res.id[1]} ↔ {fixed_anchor_res.id[1]})")

def main():
    # ---------- 读取 ----------
    af = MMCIFParser(QUIET=True).get_structure("AF", AF_CIF)
    qt = PDBParser(QUIET=True).get_structure("QT", QT_PDB)

    af_chain = af[0][CHAIN]
    qt_chain = qt[0][CHAIN]

    # ---------- 分段 ----------
    res_N = [copy_residue(qt_chain[(" ", i, " ")])        # 1–7
             for i in range(1, HELIX_START)]
    res_AF = [copy_residue(af_chain[(" ", i, " ")])       # 8–25
              for i in range(HELIX_START, HELIX_END+1)]
    max_qt = max(r.id[1] for r in qt_chain)
    res_C = [copy_residue(qt_chain[(" ", i, " ")])        # 26–end
             for i in range(HELIX_END+1, max_qt+1)]

    # ---------- 刚性对齐 ----------
    # • N 段整体旋转+平移，使 QT res7 → AF res8
    rigid_align_segment(
        moving_residues = res_N,
        fixed_anchor_res = af_chain[(" ", HELIX_START, " ")],      # AF 8
        moving_anchor_res= res_N[-1]                               # QT 7
    )
    # • C 段整体旋转+平移，使 QT res26 → AF res25
    rigid_align_segment(
        moving_residues = res_C,
        fixed_anchor_res = af_chain[(" ", HELIX_END, " ")],        # AF 25
        moving_anchor_res= res_C[0]                                # QT 26
    )

    # ---------- 组装 ----------
    merged = Structure.Structure("STITCHED")
    mdl = Model.Model(0); merged.add(mdl)
    chain_new = Chain.Chain(CHAIN); mdl.add(chain_new)
    for seg in (res_N, res_AF, res_C):
        for r in seg:
            chain_new.add(r)

    # ---------- 写出 ----------
    io = PDBIO(); io.set_structure(merged); io.save(OUT_PDB)
    print(f"\n✅ 拼接并对齐完成 → {OUT_PDB}")

if __name__ == "__main__":
    main()
