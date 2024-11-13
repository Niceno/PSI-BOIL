This directory includes the simulation of falling liquid film on wall with phase change
- Air (mon condensable gas) is included in the gas phase
- The variable eps: the volume fraction of air in the gas phase [-]

Known bug:
- Restart from *.bck has problem. (It can restart, but instable)

EquilibriumCondition
   Initial and inlet: T = 60 [C] and eps = 0.8 [-] uniform
   Phase change must not take place because Tsat(eps=0.8) = 60 [C]

Superheat
   Inlet temperature: 100 [C] for both phases
   Inlet eps: 0.3; Tint(eps=0.3) = 90 [C]
   Degree of superheat = 100 - 90 = 10 [K]

ExpYan
   Similar to the experimental condition of Yan et al. IJHMT 1991, 1105-1111
   Domain size: experiment: x (channel half width)= 0.0075 ~ 0.025
                            z (height)            = 1.0
                simulation: x (channel half width)= 0.0024
                            z (height)            = 0.0128

