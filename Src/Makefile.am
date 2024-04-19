bin_PROGRAMS = Boil

Boil_SOURCES = main.cpp    \
               Ravioli/*.h

SUBDIRS = Equation/Staggered/Momentum      \
          Equation/Staggered               \
	  Solver/Additive                  \
          Solver/Linear/Gauss              \
          Solver/Linear/Iterative          \
          Solver/Linear/Krylov             \
          Solver/Linear                    \
          Solver/Preconditioner            \
          Equation/Centered/Enthalpy       \
          Equation/Centered/EnthalpyFD     \
          Equation/Centered/EnthalpyTif    \
          Equation/Centered/Concentration  \
          Equation/Centered/Pressure       \
          Equation/Centered/PhaseChange    \
          Equation/Centered/PhaseChangeVOF \
          Equation/Centered/CIPCSL2        \
          Equation/Centered/VOF            \
          Equation/Centered/Distance       \
          Equation/Centered                \
          Equation/Dispersed               \
          Equation/Floodfill               \
          Equation/Heaviside               \
          Equation/Lagrangian              \
          Equation/Nucleation              \
          Equation/Pathline                \
          Equation/Tifmodel                \
          Equation/Topology                \
          Monitor/Rack                     \
          Monitor/Location                 \
          Matter                           \
          LookUpTable                      \
          Grid                             \
          Plot                             \
          Plot/TEC                         \
          Plot/VTK                         \
          Domain                           \
          Body                             \
          Model                            \
          Field/Vector                     \
          Field/Scalar                     \
          Field/ScalarBool                 \
          Field/VectorBool                 \
          Field/ScalarInt                  \
          Timer                            \
          Formula                          \
          Profile                          \
          RandomFlow                       \
          Global                           \
          Boundary                         \
          Ravioli                          \
          SimulationTime                   \
          Parallel/Out                     \
          Parallel                         \
          Template                         \
          Custom/IF97

LDADD = Equation/Staggered/Momentum/libMomentum.a            \
        Equation/Staggered/libStaggered.a                    \
	Solver/Additive/libAdditive.a                        \
        Solver/Linear/Gauss/libGauss.a                       \
        Solver/Linear/Iterative/libIterative.a               \
        Solver/Linear/Krylov/libKrylov.a                     \
        Solver/Linear/libLinear.a                            \
        Solver/Preconditioner/libPreconditioner.a            \
        Equation/Centered/Enthalpy/libEnthalpy.a             \
        Equation/Centered/EnthalpyFD/libEnthalpyFD.a         \
        Equation/Centered/EnthalpyTif/libEnthalpyTif.a       \
        Equation/Centered/Concentration/libConcentration.a   \
        Equation/Centered/Pressure/libPressure.a             \
        Equation/Centered/PhaseChange/libPhaseChange.a       \
        Equation/Centered/PhaseChangeVOF/libPhaseChangeVOF.a \
        Equation/Centered/CIPCSL2/libCipcsl2.a               \
        Equation/Centered/VOF/libVof.a                       \
        Equation/Centered/Distance/libDistance.a             \
        Equation/Centered/libCentered.a                      \
        Equation/Dispersed/libDispersed.a                    \
        Equation/Floodfill/libFloodfill.a                    \
	Equation/Heaviside/libHeaviside.a                    \
	Equation/Lagrangian/libLagrangian.a                  \
	Equation/Nucleation/libNucleation.a                  \
	Equation/Pathline/libPathline.a                      \
        Equation/Tifmodel/libTif.a                           \
        Equation/Topology/libTopology.a                      \
        Monitor/Rack/libRack.a                               \
        Monitor/Location/libLocation.a                       \
        Matter/libMatter.a                                   \
        LookUpTable/libLookUpTable.a                         \
        Grid/libGrid.a                                       \
        Plot/libPlot.a                                       \
        Plot/TEC/libPlotTEC.a                                \
        Plot/VTK/libPlotVTK.a                                \
        Domain/libDomain.a                                   \
        Body/libBody.a                                       \
        Model/libModel.a                                     \
        Field/Vector/libVector.a                             \
        Field/Scalar/libScalar.a                             \
        Field/VectorBool/libVectorBool.a                     \
        Field/ScalarBool/libScalarBool.a                     \
        Field/ScalarInt/libScalarInt.a                       \
        Timer/libTimer.a                                     \
        Formula/libFormula.a                                 \
        Profile/libProfile.a                                 \
        RandomFlow/libRandomFlow.a                           \
        Global/libGlobal.a                                   \
        Boundary/libBoundary.a                               \
        Ravioli/libRavioli.a                                 \
        SimulationTime/libSimulationTime.a                   \
        Parallel/Out/libOut.a                                \
        Parallel/libParallel.a                               \
        Template/libTemplate.a                               \
        Domain/libDomain.a                                   \
        Custom/IF97/libIF97.a
### Domain/libDomain.a is define twice because of compile error in the case that
### vector_bnd_extract_global is used.
