* n16_full.cmd - Device 1 (Fig. 5(a), ref. [11]): gate sweep at V_DS = 0 for the compact-model check (Step 3e)
* Base: n13_full.cmd. Changes: isothermal (no current flows at V_DS = 0, so no self-heating),
* gate work function 6.415 eV [CALIBRATE, n13-n15], drain held at 0 V,
* electrons under the gate written to the .plt at every bias point (CurrentPlot).

File {
  Grid=      "n1_msh.tdr"
  Parameter= "sdevice_n13.par"
  Plot=      "n16_des.tdr"
  Current=   "n16_des.plt"
  Output=    "n16_des.log"
}

* [ASSUME] surface donors at the graded-AlGaN/SiN surface (same as n13)
Physics (MaterialInterface= "AlGaN/Nitride") {
  Traps ( ( Donor Conc=3.3e13 Level FromMidBandGap EnergyMid=1.2
            eXsection=1e-30 hXsection=1e-14 Add2TotalDoping(ChargedTraps) ) )
}

* [ASSUME] acceptor traps in the GaN (same as n13)
Physics (Region="channel") {
  Traps ( (Acceptor Level Conc=1e15 EnergyMid=0.45 FromMidBandGap eXSection=1e-15 hXSection=1e-15) )
}

* no polarization sheet at the GaN bottom (same as n13)
Physics (MaterialInterface="GaN/Oxide") {
  PiezoElectric_Polarization(activation= 0.0)
}

Physics {
  Temperature= 300.00
  AreaFactor=  1e+03
  Fermi
  Thermionic
  Piezoelectric_Polarization(strain)
  DefaultParametersFromFile
  EffectiveIntrinsicDensity(NoBandgapNarrowing)
  Recombination ( SRH Radiative )
  Mobility ( DopingDependence(Masetti) HighFieldSaturation( CaugheyThomas ) )
  Aniso ( Poisson direction(SimulationSystem)= (1 0 0) )
  eBarrierTunneling "GateNLM"
  eBarrierTunneling "SourceNLM"
  eBarrierTunneling "DrainNLM"
  Traps( (Acceptor SFactor = "DeepLevels" Level EnergyMid=+0.45 FromMidBandGap eXSection=1e-15 hXSection=1e-15) )
}

Plot {
  Potential eQuasiFermi eDensity hDensity SpaceCharge
  ConductionBandEnergy ValenceBandEnergy xMoleFraction
  ElectricField eCurrent
}

* Electrons under the gate: window X -1.12 to -1.02 um (grade top to 80 nm into the GaN),
* Y 1.1 to 1.6 um (inside the 1.0-1.7 um gate, away from its edges).
* Sheet density (cm^-2) = Average (cm^-3) x window height (0.1 um = 1e-5 cm). Integrate is a cross-check.
CurrentPlot {
  eDensity( Average( Window[ (-1.12 1.1) (-1.02 1.6) ] ) )
  eDensity( Integrate( Window[ (-1.12 1.1) (-1.02 1.6) ] ) )
}

Math {
  NumberOfThreads= 8
  ParallelLicense (Wait)
  Wallclock
  Digits= 5
  ErrRef(electron)= 1.000e+04
  ErrRef(hole)= 1.000e+04
  RHSMin= 1.000e-10
  Extrapolate
  ExtendedPrecision(80)
  TensorGridAniso(aniso)
  ComputeDopingConcentration
  RefDens_eGradQuasiFermi_EParallelToInterface_HFS= 1.000e+08
  RefDens_hGradQuasiFermi_EParallelToInterface_HFS= 1.000e+08
  Method= Blocked
  SubMethod= ILS(set=12 )
  ILSrc= "set(12) {
          iterative (gmres(100), tolrel=1e-10, tolunprec=1e-4, maxit=200);
          preconditioning (ilut(1e-08,-1),left);
          ordering (symmetric=nd, nonsymmetric=mpsilst);
          options(verbose=0, refineresidual=10);
          };"
  NonLocal "GateNLM"   ( Electrode="gate"   Digits= 4 Length= 12e-7 EnergyResolution= 1e-3 )
  NonLocal "SourceNLM" ( Electrode="source" Digits= 4 Length= 12e-7 EnergyResolution= 1e-3 )
  NonLocal "DrainNLM"  ( Electrode="drain"  Digits= 4 Length= 12e-7 EnergyResolution= 1e-3 )
}

Electrode {
  { Name="gate"   Voltage= 0 Schottky Workfunction= 6.415 Resist= 10. }
  { Name="drain"  Voltage= 0 Schottky Workfunction= 4.3   Resist= 400. }
  { Name="source" Voltage= 0 Schottky Workfunction= 4.3   Resist= 400. }
}

Solve {
  Coupled (Iterations= 10000 LinesearchDamping= 1e-5) { Poisson }
  Coupled (Iterations= 100) { Poisson Electron Hole }
  * ramp 1: gate 0 -> -3 V (device off)
  Quasistationary ( InitialStep= 0.01 MaxStep= 0.05 MinStep= 1e-7 Increment= 1.3
                    Goal { Name="gate" Voltage= -3.0 } ) { Coupled { Poisson Electron Hole } }
  * ramp 2: gate -3 -> +1 V, the data we use (steps of at most 0.05 V)
  Quasistationary ( InitialStep= 0.002 MaxStep= 0.0125 MinStep= 1e-7 Increment= 1.3
                    Goal { Name="gate" Voltage= 1.0 } ) { Coupled { Poisson Electron Hole } }
  Plot (FilePrefix= "n16_end")
}