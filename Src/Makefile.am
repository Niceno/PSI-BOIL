bin_PROGRAMS = Boil

Boil_SOURCES = main.cpp    \
               Ravioli/*.h

SUBDIRS = Equation/Staggered/Momentum     \
          Equation/Staggered              \
          Solver                          \
          Equation/Centered/Enthalpy      \
          Equation/Centered/EnthalpyFD    \
          Equation/Centered/Concentration \
          Equation/Centered/Pressure      \
          Equation/Centered/ColorFunction \
          Equation/Centered/LevelSet      \
          Equation/Centered/PhaseChange   \
          Equation/Centered/PhaseField    \
          Equation/Centered/ColorCIP      \
          Equation/Centered/CIPCSL2       \
          Equation/Centered/Distance      \
          Equation/Centered/VOF           \
          Equation/Centered               \
          Equation/Dispersed              \
          Equation/Floodfill              \
          Monitor/Rack                    \
          Monitor/Location                \
          Matter                          \
          LookUpTable                     \
          Grid                            \
          Plot                            \
          Plot/GMV                        \
          Plot/TEC                        \
          Plot/VTK                        \
          Domain                          \
          Board                           \
          Body                            \
          Model                           \
          Field/Vector                    \
          Field/Scalar                    \
          Field/ScalarBool                \
          Field/VectorBool                \
          Field/ScalarInt                 \
          Timer                           \
          Formula                         \
          Profile                         \
          RandomFlow                      \
          Global                          \
          Boundary                        \
          Ravioli                         \
          SimulationTime                  \
          Parallel/Out                    \
          Parallel                        \
          Template

LDADD = Equation/Staggered/Momentum/libMomentum.a          \
        Equation/Staggered/libStaggered.a                  \
        Solver/libSolver.a                                 \
        Equation/Centered/Enthalpy/libEnthalpy.a           \
        Equation/Centered/EnthalpyFD/libEnthalpyFD.a       \
        Equation/Centered/Concentration/libConcentration.a \
        Equation/Centered/Pressure/libPressure.a           \
        Equation/Centered/ColorFunction/libColorFunction.a \
        Equation/Centered/LevelSet/libLevelSet.a           \
        Equation/Centered/PhaseChange/libPhaseChange.a     \
        Equation/Centered/PhaseField/libPhaseField.a       \
        Equation/Centered/ColorCIP/libColorCIP.a           \
        Equation/Centered/CIPCSL2/libCipcsl2.a             \
        Equation/Centered/Distance/libDistance.a           \
        Equation/Centered/VOF/libVof.a                \
        Equation/Centered/libCentered.a                    \
        Equation/Dispersed/libDispersed.a                  \
        Equation/Floodfill/libFloodfill.a                  \
        Monitor/Rack/libRack.a                             \
        Monitor/Location/libLocation.a                     \
        Matter/libMatter.a                                 \
        LookUpTable/libLookUpTable.a                       \
        Grid/libGrid.a                                     \
        Plot/libPlot.a                                     \
        Plot/GMV/libPlotGMV.a                              \
        Plot/TEC/libPlotTEC.a                              \
        Plot/VTK/libPlotVTK.a                              \
        Domain/libDomain.a                                 \
        Board/libBoard.a                                   \
        Body/libBody.a                                     \
        Model/libModel.a                                   \
        Field/Vector/libVector.a                           \
        Field/Scalar/libScalar.a                           \
        Field/VectorBool/libVectorBool.a                   \
        Field/ScalarBool/libScalarBool.a                   \
        Field/ScalarInt/libScalarInt.a                     \
        Timer/libTimer.a                                   \
        Formula/libFormula.a                               \
        Profile/libProfile.a                               \
        RandomFlow/libRandomFlow.a                         \
        Global/libGlobal.a                                 \
        Boundary/libBoundary.a                             \
        Ravioli/libRavioli.a                               \
        SimulationTime/libSimulationTime.a                 \
        Parallel/Out/libOut.a                              \
        Parallel/libParallel.a                             \
        Template/libTemplate.a
