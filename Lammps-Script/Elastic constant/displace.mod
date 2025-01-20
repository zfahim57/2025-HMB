# NOTE: This script should not need to be
# modified. See in.elastic for more info.
#
# Find which reference length to use

if "${dir} == 1" then &
   "variable len0 equal ${lx0}" 
if "${dir} == 2" then &
   "variable len0 equal ${ly0}" 
if "${dir} == 3" then &
   "variable len0 equal ${lz0}" 
if "${dir} == 4" then &
   "variable len0 equal ${lz0}" 
if "${dir} == 5" then &
   "variable len0 equal ${lz0}" 
if "${dir} == 6" then &
   "variable len0 equal ${ly0}" 

# Reset box and simulation parameters

#clear
#box tilt large
read_dump dump.npt_${T}K ${N} x y z vx vy vz
include potential.mod
reset_timestep 0

# Negative deformation

variable delta equal -${up}*${len0}
variable deltaxy equal -${up}*xy
variable deltaxz equal -${up}*xz
variable deltayz equal -${up}*yz
if "${dir} == 1" then &
   "change_box all x delta 0 ${delta} xy delta ${deltaxy} xz delta ${deltaxz} remap units box"
if "${dir} == 2" then &
   "change_box all y delta 0 ${delta} yz delta ${deltayz} remap units box"
if "${dir} == 3" then &
   "change_box all z delta 0 ${delta} remap units box"
if "${dir} == 4" then &
   "change_box all yz delta ${delta} remap units box"
if "${dir} == 5" then &
   "change_box all xz delta ${delta} remap units box"
if "${dir} == 6" then &
   "change_box all xy delta ${delta} remap units box"

# Relax atoms positions

#minimize ${etol} ${ftol} ${maxiter} ${maxeval}

# NVT equilibration

# Obtain new stress tensor

variable step equal step
variable temp equal temp
variable press equal press
variable pxx equal pxx
variable pyy equal pyy
variable pzz equal pzz
variable pxy equal pxy
variable pxz equal pxz
variable pyz equal pyz

#dump equidump all custom ${tdumpeq} dump.nvt_${dir}_${T}K id type x y z #vx vy vz fx fy fz 
fix nvt all nvt temp ${T} ${T} $(100.0*dt)
fix aves1 all ave/time 1 ${tequihalf} ${tequi} v_pxx v_pyy v_pzz v_pxy v_pxz v_pyz file aves1neg_${dir}_${T}K.press

fix output all print 10           &
                      "${step} ${temp} ${press} ${pxx} ${pyy} ${pzz} ${pxy} ${pxz} ${pyz}"       &
                      title "# step temp press pxx pyy pzz pxy pxz pyz" screen no &
                      file equi_nvtneg_${T}K_${dir}.dat

run ${tequi}
unfix nvt
#undump equidump
unfix output

# Obtain new stress tensor

variable tmp equal pe
variable e1 equal ${tmp}
variable tmp equal press
variable p1 equal ${tmp}
variable tmp equal f_aves1[1]
variable pxx1 equal ${tmp}
variable tmp equal f_aves1[2]
variable pyy1 equal ${tmp}
variable tmp equal f_aves1[3]
variable pzz1 equal ${tmp}
variable tmp equal f_aves1[4]
variable pxy1 equal ${tmp}
variable tmp equal f_aves1[5]
variable pxz1 equal ${tmp}
variable tmp equal f_aves1[6]
variable pyz1 equal ${tmp}

# Compute elastic constant from pressure tensor

variable C1neg equal ${d1}
variable C2neg equal ${d2}
variable C3neg equal ${d3}
variable C4neg equal ${d4}
variable C5neg equal ${d5}
variable C6neg equal ${d6}

unfix aves1

# Reset box and simulation parameters

#clear
#box tilt large
#read_restart ../equi_${T}K/restart.npt_${T}K
read_dump dump.npt_${T}K ${N} x y z vx vy vz
include potential.mod

# Obtain new stress tensor

variable step equal step
variable temp equal temp
variable press equal press
variable pxx equal pxx
variable pyy equal pyy
variable pzz equal pzz
variable pxy equal pxy
variable pxz equal pxz
variable pyz equal pyz

# Positive deformation

variable delta equal ${up}*${len0}
variable deltaxy equal ${up}*xy
variable deltaxz equal ${up}*xz
variable deltayz equal ${up}*yz
if "${dir} == 1" then &
   "change_box all x delta 0 ${delta} xy delta ${deltaxy} xz delta ${deltaxz} remap units box"
if "${dir} == 2" then &
   "change_box all y delta 0 ${delta} yz delta ${deltayz} remap units box"
if "${dir} == 3" then &
   "change_box all z delta 0 ${delta} remap units box"
if "${dir} == 4" then &
   "change_box all yz delta ${delta} remap units box"
if "${dir} == 5" then &
   "change_box all xz delta ${delta} remap units box"
if "${dir} == 6" then &
   "change_box all xy delta ${delta} remap units box"

# Relax atoms positions
# NVT equilibration

#dump equidump all custom ${tdumpeq} dump.nvtpos_${dir}_${T}K id type x y z #vx vy vz fx fy fz
fix nvt all nvt temp ${T} ${T} $(100.0*dt)
fix aves1 all ave/time 1 ${tequihalf} ${tequi} v_pxx v_pyy v_pzz v_pxy v_pxz v_pyz file aves1pos_${dir}_${T}K.press

fix output all print 10           &
                      "${step} ${temp} ${press} ${pxx} ${pyy} ${pzz} ${pxy} ${pxz} ${pyz}"       &
                      title "# step temp press pxx pyy pzz pxy pxz pyz" screen no &
                      file equi_nvtpos_${T}K_${dir}.dat

run ${tequi}
unfix nvt
#unfix aves1
#undump equidump
unfix output

# Obtain new stress tensor

variable tmp equal f_aves1[1]
variable pxx1 equal ${tmp}
variable tmp equal f_aves1[2]
variable pyy1 equal ${tmp}
variable tmp equal f_aves1[3]
variable pzz1 equal ${tmp}
variable tmp equal f_aves1[4]
variable pxy1 equal ${tmp}
variable tmp equal f_aves1[5]
variable pxz1 equal ${tmp}
variable tmp equal f_aves1[6]
variable pyz1 equal ${tmp}

# Compute elastic constant from pressure tensor

variable C1pos equal ${d1}
variable C2pos equal ${d2}
variable C3pos equal ${d3}
variable C4pos equal ${d4}
variable C5pos equal ${d5}
variable C6pos equal ${d6}

# Combine positive and negative

variable C1${dir} equal 0.5*(${C1neg}+${C1pos})
variable C2${dir} equal 0.5*(${C2neg}+${C2pos})
variable C3${dir} equal 0.5*(${C3neg}+${C3pos})
variable C4${dir} equal 0.5*(${C4neg}+${C4pos})
variable C5${dir} equal 0.5*(${C5neg}+${C5pos})
variable C6${dir} equal 0.5*(${C6neg}+${C6pos})

unfix aves1 

# Delete dir to make sure it is not reused

variable dir delete
                                
