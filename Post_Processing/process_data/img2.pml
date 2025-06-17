###############################################
# render_custom_ballstick.pml
# Abeta42.pdb → ballstick_custom.png  (4K / 600 ppi)
###############################################

# ---------- 0. 加载与全局设置 ----------
load Abeta42.pdb, prot
bg_color white
set ray_opaque_background, off
set antialias, 2
set ray_shadows, on
set stick_radius, 0.22
set sphere_scale, 0.28
orient prot

# ---------- 1. Ball-and-Stick 模型 ----------
hide everything
show sticks, prot
show spheres, prot

# 按 CPK 配色，其余原子使用默认配色
util.cpk prot

# 特殊原子着色：
# 氢原子染蓝色
color blue, elem H
# 碳原子染黄色
color yellow, elem C

# ---------- 2. 渲染与保存 ----------
ray 3840,2160
png ballstick_custom.png, dpi=600

# ---------- 完成 ----------
print "\n已生成球棍模型： ballstick_custom.png"