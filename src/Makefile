# Disable all of make's built-in rules (similar to Fortran's implicit none)
MAKEFLAGS += --no-builtin-rules --no-builtin-variables

# Compiler settings
fortran_compiler := ifort 
# compiler_flags := /debug /Od /Qdiag-disable:10448 /check:all /traceback
# compiler_flags := /Qdiag-disable:10448 /nologo /debug:full /O2 /assume:buffered_io /fpp /DDebug /reentrancy:threaded /free /warn:all /real-size:64 /Qauto /fpe:0 /fp:strict /fpconstant /module:"Debug\\" /object:"Debug\\" /Fd"Debug\vc160.pdb" /check:pointer /check:bounds /check:uninit /check:format /check:output_conversion /check:stack /libs:static /threads /Qmkl:sequential /c

# Select the operating system
ifeq ($(OS),Windows_NT)
    # Windows
    compiler_flags := /Qdiag-disable:10448 /nologo /debug:full /O2 /assume:buffered_io /fpp /DDebug /reentrancy:threaded /free /warn:all /real-size:64 /Qauto /fpe:0 /fp:strict /fpconstant /check:pointer /check:bounds /check:uninit /check:format /check:output_conversion /check:stack /libs:static /threads #/Qmkl:sequential /c
    EXE = .exe
    remove_files_command := del /f
else
    # Assume Linux
    compiler_flags := -O0 -g -check all -traceback
    EXE = .out
    remove_files_command := rm -f
endif 

# set the files to remove

# list of all source files
source_files := GlobalConstants.f90\
                A3DLinearElasticity.f90\
                MCSS_Soil_Model.f90\
                GeoMath.f90\
                MatrixMath.f90\
                FileIO.f90\
                String.f90\
                Counters.f90\
                Feedback.f90\
                ReadCalculationData.f90\
                ElemCalcTETRA.f90\
                ElemCalcQUAD.f90\
                ReadMaterialData.f90\
                ElemCalcTRI.f90\
                getversion.f90\
                InitialiseKernel.f90\
                ElemCalc.f90\
                MeshInfo.f90\
                ReadGeometryData.f90\
                ISORT.f90\
                ElemConnections.f90\
                Particle.f90\
                MPMData.f90\
                RotBoundCond.f90\
                MPMStresses.f90\
                MPMDynViscousBoundary.f90\
                TwoLayerFormulation.f90\
                MPMDYN2PhaseSP.f90\
                MPMDYN3PhaseSP.f90\
                MPMDynContact.f90\
                ReadMPMData.f90\
                WriteMPMData.f90\
                MPMMeshAdjustment.f90\
                MPMConvPhase.f90\
                WriteTestData.f90\
                MPMEmptyElements.f90\
                AdjustParticleDiscretisation.f90\
                MPMInit.f90\
                MPMDYNBTSig.f90\
                RigidBody.f90\
                LagrangianPhase.f90\
                Liquid.f90\
                ExternalSoilModel.f90\
                MPMStrainSmoothing.f90\
                MPMDYNStresses.f90\
                MPMDYNConvPhase.f90\
                MPMExcavation.f90\
                BuildBJacDet.f90\
                BuildDElastic.f90\
                BuildLoad.f90\
                GetPrinStress.f90\
                GetStrain.f90\
                InitialiseElementType.f90\
                timing.f90\
                WriteVTKASCII.f90\
                WriteVTKBinary.f90\
                WriteVTKOutput.f90\
                WriteVTK2Layer.f90\
                WriteNodalData.f90\
                WriteResultData.f90\
                ErrorHandler.f90\
                MPMDynamicExplicit.f90\
                Kernel.f90\
                Anura3D.f90

# Object files
object_files := $(patsubst %.f90,%.obj,$(source_files))

# Target executable
TARGET := Anura3D$(EXE)

# Compilation rule for Fortran files
%.obj: %.f90
	$(fortran_compiler) $(compiler_flags) -c $< -o $@

# Linking rule
$(TARGET): $(object_files)
	$(fortran_compiler) $(compiler_flags) $^ -o $@

# Phony target for cleaning
.PHONY: clean
clean:
	$(remove_files_command) *.obj *.mod *genmod*.f90 *.pdb *.ilk $(TARGET)