* Step 13d 2026-11-24: band/density snapshots for run @NAME@ at the second turn-on (user choice A+B, 24 Nov).
*   Built from lg_ext_des.cmd (Step 13c), unchanged except:
*   1. output names @NAME@s_* (Step 13b/13c files untouched); same mesh @NAME@_msh.tdr;
*   2. Plot snapshots in the +0.8 -> +2.0 V sweep at Time t, V_GS = 0.8 + 1.2 t:
*      1.000 1.600 1.700 1.750 1.775 1.785 1.800 1.900 V (a200: steep step between +1.766 and +1.785 V, Step 13c).
*   Read with lg_snap.py (per-layer n_s and E_C - E_Fn under the gate). Snapshot syntax as n8_full.cmd (Step 6b-5).
*   Check: IdVg_@NAME@s_des.plt must reproduce IdVg_@NAME@x_des.plt (same deck, same steps apart from the plot times).
*   Hypothesis under test [ASSUME]: the second turn-on is electrons filling the x = 0.27 AlGaN under the gate.
* Physics, electrodes, traps, Math: identical to lg_des.cmd / lg_ext_des.cmd.

File {
  Grid=    "@NAME@_msh.tdr"
  Parameter="sdevice_n5.par"

  Plot=   "@NAME@s_des.tdr"
  Current="@NAME@s_des.plt"
  Output= "@NAME@s_des.log"
}

* Gate: Schottky, work function 5.45 eV, the value calibrated on Structure A (Sec. II, run n4)
* [CALIBRATE there; start value here, to be set from Device 2's turn-on].
* Source/drain: ohmic on n+ GaN [P, Fig. 5(b)]; Resist as in the Sec. II deck [ASSUME].
Electrode {
  { Name="gate"   Voltage= 0 Schottky Workfunction= 5.45 Resist= 10. }
  { Name="drain"  Voltage= 0 Resist= 10. }
  { Name="source" Voltage= 0 Resist= 10. }
}

* [ASSUME] surface donors (from the Synopsys example), at the barrier/SiN surface
Physics (MaterialInterface= "AlGaN/Nitride") {
  Traps (
    (
     Donor
     Conc=3e13
     Level
     FromMidBandGap
     EnergyMid=0.36
     eXsection=1e-30
     hXsection=1e-14
     Add2TotalDoping(ChargedTraps)
    )
  )
}

* [ASSUME] acceptor traps in the 40 nm GaN channel (not given in [14] or the paper)
Physics (Region="channel") {
  Traps(
    (Acceptor Level Conc= 1e15 EnergyMid= 0.45 FromMidBandGap
     eXSection=1e-15 hXSection=1e-15 )
  )
}

* [ASSUME] the Oxide block is only a stand-in for the unspecified substrate;
* switch off the polarization sheet at its interface with the Al0.04GaN buffer.
Physics (MaterialInterface="AlGaN/Oxide") {
  PiezoElectric_Polarization(activation= 0.0)
}

Plot {
  Potential eQuasiFermi hQuasiFermi
  eDensity hDensity SpaceCharge
  Current eCurrent hCurrent CurrentPotential
  ElectricField SemiconductorElectricField InsulatorElectricField
  eMobility hMobility
  eVelocity hVelocity
  Doping DonorConcentration AcceptorConcentration
  PE_Polarization/vector
  PE_Charge PiezoCharge
  ConductionBandEnergy ValenceBandEnergy
  ElectronAffinity
  BandGap EffectiveBandGap BandGapNarrowing
  QCEffectiveBandGap
  xMolefraction yMolefraction
  LatticeTemperature eTemperature hTemperature
  lHeatFlux TotalHeat
  eJouleHeat hJouleHeat JouleHeat
  PeltierHeat ThomsonHeat RecombinationHeat
  DissipationRateDensity
  OpticalAbsorptionHeat
  ThermalConductivity LatticeHeatCapacity
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
  Transient= BE
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
  NonLocal "GateNLM"   (
    Electrode="gate"
    Digits= 4
    Length= 12e-7
    EnergyResolution= 1e-3
  )
}

Physics {
  Temperature=  300.00
  AreaFactor=    1e+03

  Fermi
  Thermionic
  Piezoelectric_Polarization(strain)
  DefaultParametersFromFile
  EffectiveIntrinsicDensity(NoBandgapNarrowing)

  Recombination (
    SRH
    Radiative
  )

  Mobility (
    DopingDependence(Masetti)
    HighFieldSaturation(
      CaugheyThomas
    )
  )

  Aniso (
    Poisson
    direction(SimulationSystem)= (1 0 0)
  )

  eBarrierTunneling "GateNLM"

  Traps(
    (Acceptor SFactor = "DeepLevels" Level  EnergyMid=+0.45 FromMidBandGap
     eXSection=1e-15 hXSection=1e-15 )
  )
}

Solve {
  Coupled (Iterations= 10000 LinesearchDamping= 1e-5) { Poisson }
  Coupled (Iterations= 100) { Poisson Electron Hole }
  * drain to 5 V at V_GS = 0 (Fig. 6(b) bias [P]), as Step 13b and n7
  Quasistationary (
    InitialStep= 0.01 MaxStep= 0.02 MinStep= 1e-6
    Goal { Name= "drain" Voltage= 5 }
  ) { Coupled { Poisson Electron Hole } }
  * gate to +0.8 V, coarse (as the Step 13b ramp to -6 V)
  Quasistationary (
    InitialStep= 0.01 MaxStep= 0.05 MinStep= 1e-5
    Goal { Name= "gate" Voltage= 0.8 }
  ) { Coupled { Poisson Electron Hole } }
  * the sweep: +0.8 V to +2.0 V, at most 1.2 V x 0.0175 = 21 mV per step (first step 10.5 mV, as Step 13b)
  NewCurrentFile= "IdVg_"
  Quasistationary (
    InitialStep= 0.00875 MaxStep= 0.0175 MinStep= 1e-6
    Goal { Name= "gate" Voltage= 2 }
  ) { Coupled { Poisson Electron Hole }
      Plot ( FilePrefix= "@NAME@s_vg" Time= (0.16667; 0.66667; 0.75; 0.79167; 0.8125; 0.82083; 0.83333; 0.91667) NoOverwrite )
  }
}
