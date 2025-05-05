# --*-- conding:utf-8 --*--
# @time:5/5/25 17:14
# @Author : Yuqi Zhang
# @Email : yzhan135@kent.edu
# @File:create_42.py


import numpy as np
from Bio.PDB import (
    MMCIFParser, PDBParser,
    Structure, Model, Chain, Residue, PDBIO,
    Superimposer
)

# —— 固定路径 ——
AI_CIF_FILE   = "fold_42/fold_42_model_0.cif"
QT_PDB_FILE   = "Abeta_42.pdb"
OUTPUT_PDB    = "AI_Q_abeta42.pdb"

# —— 螺旋区段（按实验/模拟常见 Helix‐1） ——
CHAIN_ID     = "A"
HELIX_START  = 8
HELIX_END    = 25

# —— 碰撞阈值（Å） ——
CLASH_THRESHOLD = 1.5

def load_structures():
    """读取 AI CIF 和 Quantum PDB"""
    cif_p = MMCIFParser(QUIET=True)
    ai_struc = cif_p.get_structure("AI", AI_CIF_FILE)
    pdb_p = PDBParser(QUIET=True)
    qt_struc = pdb_p.get_structure("QT", QT_PDB_FILE)
    return ai_struc, qt_struc

def align_helix(ai_struc, qt_struc):
    """刚性对齐 AI 螺旋段到 Quantum 结构中对应的 Helix 残基位置"""
    ai_chain = ai_struc[0][CHAIN_ID]
    qt_chain = qt_struc[0][CHAIN_ID]

    def ca_atom(chain, resseq):
        return chain[(" ", resseq, " ")]["CA"]

    # 取两端 Cα 原子进行对齐
    moving = [ ca_atom(ai_chain, HELIX_START),
               ca_atom(ai_chain, HELIX_END) ]
    fixed  = [ ca_atom(qt_chain, HELIX_START),
               ca_atom(qt_chain, HELIX_END) ]

    sup = Superimposer()
    sup.set_atoms(fixed, moving)

    # 对 AI 螺旋段所有原子应用同样变换
    helix_atoms = []
    for resseq in range(HELIX_START, HELIX_END+1):
        res = ai_chain[(" ", resseq, " ")]
        helix_atoms.extend(list(res.get_atoms()))
    sup.apply(helix_atoms)
    print(f"螺旋段对齐完成，RMSD = {sup.rms:.3f} Å")

def merge_structures(ai_struc, qt_struc):
    """将对齐后的 AI 螺旋插入 Quantum 骨架"""
    merged = Structure.Structure("MERGED")
    new_model = Model.Model(0)
    merged.add(new_model)

    qt_model = qt_struc[0]
    # 复制 Quantum 结构非螺旋区段
    for ch in qt_model:
        new_ch = Chain.Chain(ch.id)
        new_model.add(new_ch)
        for res in ch:
            het, idx, icode = res.id
            if ch.id == CHAIN_ID and HELIX_START <= idx <= HELIX_END:
                continue
            nr = Residue.Residue(res.id, res.resname, res.segid)
            for atom in res:
                nr.add(atom.copy())
            new_ch.add(nr)

    # 插入已对齐的 AI 螺旋段
    ai_chain = ai_struc[0][CHAIN_ID]
    tgt_chain = merged[0][CHAIN_ID]
    for resseq in range(HELIX_START, HELIX_END+1):
        res = ai_chain[(" ", resseq, " ")]
        nr = Residue.Residue(res.id, res.resname, res.segid)
        for atom in res:
            nr.add(atom.copy())
        tgt_chain.add(nr)

    return merged

def detect_clashes(structure):
    """检测原子碰撞，小于阈值则报告"""
    atoms = [a for m in structure for c in m for r in c for a in r]
    coords = np.array([a.get_coord() for a in atoms])
    n = len(atoms)
    clashes = []
    for i in range(n):
        res_i = atoms[i].get_parent().id[1]
        diffs = coords[i+1:] - coords[i]
        d2 = np.sum(diffs*diffs, axis=1)
        for j_off, dist2 in enumerate(d2):
            if dist2 < CLASH_THRESHOLD**2:
                j = i + 1 + j_off
                res_j = atoms[j].get_parent().id[1]
                if res_i == res_j:
                    continue
                clashes.append((i, j, np.sqrt(dist2)))
    if not clashes:
        print("✅ 未检测到原子碰撞。")
    else:
        print(f"⚠️ 检测到 {len(clashes)} 处可能碰撞（< {CLASH_THRESHOLD} Å）：")
        for idx, (i, j, d) in enumerate(clashes[:10], 1):
            a1, a2 = atoms[i], atoms[j]
            print(f"  {idx}. Residue {a1.get_parent().id[1]}:{a1.name}"
                  f" – {a2.get_parent().id[1]}:{a2.name}, 距离={d:.2f} Å")
        if len(clashes) > 10:
            print(f"  ... 还有 {len(clashes)-10} 处冲突未列出。")

def write_pdb(structure):
    io = PDBIO()
    io.set_structure(structure)
    io.save(OUTPUT_PDB)
    print(f"合并并检查完成 → {OUTPUT_PDB}")

if __name__ == "__main__":
    ai_struc, qt_struc = load_structures()
    align_helix(ai_struc, qt_struc)
    merged = merge_structures(ai_struc, qt_struc)
    detect_clashes(merged)
    write_pdb(merged)
