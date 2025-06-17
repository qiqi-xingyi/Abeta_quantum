###############################################
# render5_surface_sticks.pml
#  stitched_aligned_abeta42.pdb  → 5 张 4K / 600 ppi
###############################################

# ---------- 0. 加载与全局设置 ----------
load Abeta42.pdb, prot
bg_color white
set ray_opaque_background, off
set antialias, 2
set ray_shadows, on
set cartoon_fancy_helices, 1
set stick_radius, 0.22
set sphere_scale, 0.28
orient prot

# ---------- 1–3. Cartoon（默认配色，三视角） ----------
hide everything
show cartoon, prot
zoom prot

# (1) 正视角
ray 3840,2160
png cartoon_A.png, dpi=600

# (2) Y 轴 +120°
rotate y, 120, prot
ray 3840,2160
png cartoon_B.png, dpi=600

# (3) 在 (2) 基础上 Z 轴 +90°
rotate z, 90, prot
ray 3840,2160
png cartoon_C.png, dpi=600

# ---------- 4. 半透明表面 + 内部 sticks（残基彩虹色） ----------
orient prot
hide everything
# 内部分子结构：sticks（无球！）
show sticks, prot
util.cbay prot          # PyMOL 内置按残基彩虹上色
# 半透明外表面
show surface, prot
set transparency, 0.4
# 确保内部绝无球体
hide spheres, prot
set stick_ball, off
ray 3840,2160
png surface_innersticks.png, dpi=600

# ---------- 5. Ball‑and‑Stick（元素色，含球体） ----------
orient prot
hide everything
show sticks,   prot
show spheres,  prot          # 这一步才带球体
util.cpk  prot               # CPK 元素色：C/N/O/S/H…
ray 3840,2160
png ballstick_element.png, dpi=600

# ---------- 完成 ----------
print "\n已生成 5 个 4K / 600 ppi PNG："
print "  cartoon_A.png"
print "  cartoon_B.png"
print "  cartoon_C.png"
print "  surface_innersticks.png"
print "  ballstick_element.png"

