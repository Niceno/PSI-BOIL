LIST OF Tests

ALGORITHM
  NavierStokes:       Navier Stokes equations
  GhostFluidMethod:   Ghost fluid method
  Simple:             SIMPLE algorithm
  TwoDomains:         Two independent domains
  TwoDomainsVelocity: Two domains linked with insert boundary

UTILITIES
  LookUpTable:        Table for material properties (e.g. temperature dependent)
  Pathline:           Draw particle path lines
  Randomflow-inlet:   Randomflow generated for entier inlet plane
  Random-plane:       Randomflow generated in each decomposed domain
  Restart:            How to restart simulation by Yohei
  VariableProperties: Variable properties

IMMERSED BOUNDARY METHOD
  Body:
  Car:

APPLICATIONS
  Bubble3D_ReferenceFrame  Rising bubble in stagnant liquid with reference frame (PID controller)
  MHD:                Magnetro hydro-dynamics. Rising air bubble in liquid metal under magnetic field 
  PhaseChange:        Phase change simulation
    |--NucleateBoil   Nucleate boiling simulation
  RisingBubble:       rising bubble using ColorFunction by Bojan
  Zeppelin:           zeppelin

SCRIPT (Automatically run simulations)
  CheckVOF:           Verifications and validations for VOF
