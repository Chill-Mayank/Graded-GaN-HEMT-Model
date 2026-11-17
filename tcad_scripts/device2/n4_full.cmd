* n4 2026-11-10: n3's classical equilibrium + strain fix. Step 6b-2.
*   Classical: n3's density-gradient step failed (|RHS| x 1e10) even from the converged classical start.
*   Density gradient is not in the paper's TCAD list (Fermi, polarization, SRH [P]) and Device 1 ran without
*   it, so Device 2 runs classical; eQuantumPotential is removed from Physics and Solve (same result as
*   n3_Classical, where it was on but not solved).
*   The one physics change: the relaxed layer is the 0.8 um Al0.04GaN buffer, not the 40 nm GaN
*   (sdevice_n4.par) [Ambacher et al., JAP 85, 3222, 1999].
* n3 2026-11-10: n2 with a staged equilibrium solve [NUM]; physics unchanged. Step 6b-2.
*   n2 diverged in { Poisson eQuantumPotential } (|step| grew to 88 V in 100 iterations) and then in the
*   full coupled solve, so its plot is not a solution. Standard continuation: solve the classical
*   drift-diffusion equilibrium first, then add the density-gradient equation on top of it.
* n2 2026-11-10: Device 2 (paper Fig. 5(b)), equilibrium only (V_GS = V_DS = 0), Step 6b-2.
* Built from GaN_SecII/n4_full.cmd (Sec. II, 2026-10-15). Same physics, traps and solver; the one
* change is forced by the new structure: source and drain sit on regrown n+ GaN [P, Fig. 5(b)], so
* they are ohmic contacts [standard; Sentaurus Device User Guide T-2022.03, Electrode] instead of
* 4.3 eV Schottky contacts with barrier tunnelling. The source/drain tunnelling lines are removed.
* Tags: [P] paper  [R14] ref. [14]  [EX] Synopsys GaN example  [CALIBRATE]  [ASSUME]  [NUM]

File {
  Grid=    "n1_msh.tdr"
  Parameter="sdevice_n4.par"

  Plot=   "n4_des.tdr"
  Current="n4_des.plt"
  Output= "n4_des.log"
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
  Plot(FilePrefix="n4_Equi")
}
