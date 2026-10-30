# cut_n2.tcl - Sentaurus Visual: vertical cut through the Device 1 equilibrium solution (n2, VGS = VDS = 0)
# Cut at Y = 3.5 um (gate-drain gap), running along X (depth). Makes the CSV used in step1c_fig1b.ipynb.
load_file /home/TESLA/Synopsys1/STDB/hemt/GaN_Long_Channel/n2_Equi_des.tdr
create_plot -dataset n2_Equi_des
set cut [create_cutline -plot Plot_n2_Equi_des -type y -at 3.5]
export_variables {X eDensity ConductionBandEnergy eQuasiFermiPotential xMoleFraction} -dataset $cut -filename /home/TESLA/Synopsys1/STDB/hemt/GaN_Long_Channel/device1_n2_equilibrium_cut_y3.5um.csv -overwrite