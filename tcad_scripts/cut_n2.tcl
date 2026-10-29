#/home/TESLA/Synopsys1/STDB/hemt/GaN_Long_Channel/cut_n2.tcl
#S-Visual: Output from Tcl command history.
load_file /home/TESLA/Synopsys1/STDB/hemt/GaN_Long_Channel/n2_Equi_des.tdr -fod
create_plot -dataset n2_Equi_des
select_plots {Plot_n2_Equi_des}
# Plot_n2_Equi_des
# Plot_n2_Equi_des
# n2_Equi_des
create_cutline -plot Plot_n2_Equi_des -type y -at 3.5
# C1(n2_Equi_des)
create_plot -dataset C1(n2_Equi_des) -1d
select_plots {Plot_C1(n2_Equi_des)}
# Plot_C1(n2_Equi_des)
# Plot_C1(n2_Equi_des)
set_axis_prop -plot Plot_C1(n2_Equi_des) -axis y -type log
# 0
create_curve -axisX X -axisY DopingConcentration -dataset {C1(n2_Equi_des)} -plot Plot_C1(n2_Equi_des)
# Curve_1
set_axis_prop -plot Plot_C1(n2_Equi_des) -axis x -range {-5.41147 4.04147}
# 0
set_best_look {Plot_C1(n2_Equi_des)}
# 0
zoom_plot -plot Plot_C1(n2_Equi_des) -reset
zoom_plot -plot Plot_C1(n2_Equi_des) -reset
select_plots {Plot_n2_Equi_des}
# Plot_n2_Equi_des
export_variables {X eDensity ConductionBandEnergy eQuasiFermiPotential xMoleFraction} -dataset $cut -filename /home/TESLA/Synopsys1/STDB/hemt/GaN_Long_Channel/device1_n2_equilibrium_cut_y3.5um.csv -overwrite
# can't read "cut": no such variable
set cut [create_cutline -plot Plot_n2_Equi_des -type y -at 3.5]
# C2(n2_Equi_des)
export_variables {X eDensity ConductionBandEnergy eQuasiFermiPotential xMoleFraction} -dataset $cut -filename /home/TESLA/Synopsys1/STDB/hemt/GaN_Long_Channel/device1_n2_equilibrium_cut_y3.5um.csv -overwrite
# /home/TESLA/Synopsys1/STDB/hemt/GaN_Long_Channel/device1_n2_equilibrium_cut_y3.5um.csv
set p1 [create_plot -1d]
create_curve -plot $p1 -dataset $cut -axisX X -axisY eDensity
set_axis_prop -plot $p1 -axis y -type log
set_axis_prop -plot $p1 -axis x -range {-1.13 -1.05}
# 0
set_best_look {Plot_1}
# 0
undo
set_axis_prop -plot Plot_1 -axis y -type linear
# 0
set_axis_prop -plot Plot_1 -axis y -type log
# 0
set_axis_prop -plot Plot_1 -axis x -type log
# 0
set_axis_prop -plot Plot_1 -axis x -type linear
# 0
set_axis_prop -plot Plot_1 -axis y -type linear
# 0
set_axis_prop -plot Plot_1 -axis x -type log
# 0
set_axis_prop -plot Plot_1 -axis x -type linear
# 0
set_axis_prop -plot Plot_1 -axis y -type log
# 0
set_axis_prop -plot Plot_1 -axis x -range {-1.13 -1.05}
# 0
