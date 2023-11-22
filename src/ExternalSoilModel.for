	!*****************************************************************************
    !                                       ____  _____  
    !           /\                         |___ \|  __ \ 
    !          /  \   _ __  _   _ _ __ __ _  __) | |  | |
    !         / /\ \ | '_ \| | | | '__/ _` ||__ <| |  | |
    !        / ____ \| | | | |_| | | | (_| |___) | |__| |
    !       /_/    \_\_| |_|\__,_|_|  \__,_|____/|_____/ 
    !
    !
	!	Anura3D - Numerical modelling and simulation of large deformations 
    !   and soil–water–structure interaction using the material point method (MPM)
    !
    !	Copyright (C) 2023  Members of the Anura3D MPM Research Community 
    !   (See Contributors file "Contributors.txt")
    !
    !	This program is free software: you can redistribute it and/or modify
    !	it under the terms of the GNU Lesser General Public License as published by
    !	the Free Software Foundation, either version 3 of the License, or
    !	(at your option) any later version.
    !
    !	This program is distributed in the hope that it will be useful,
    !	but WITHOUT ANY WARRANTY; without even the implied warranty of
    !	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    !	GNU Lesser General Public License for more details.
    !
    !	You should have received a copy of the GNU Lesser General Public License
    !	along with this program.  If not, see <https://www.gnu.org/licenses/>.
	!
    !*****************************************************************************


module ModExternalSoilModel
!**********************************************************************
!
!  Function: Contains the routines related to calling user-defined soil models in external DLLs.
!     $Revision: 9064 $
!     $Date: 2023-08-07 13:55 -0400 (Moore, 08 Aug 2023) $
!
!**********************************************************************
      
use ModMPMData
use ModGlobalConstants
use ModReadCalculationData
use ModReadMaterialData
use ModMPMInit
use user32
use kernel32
use ModMeshInfo

contains


subroutine StressSolid(IDpt, IDel, BMatrix,IEntityID)
    !**********************************************************************
    !
    !    Function:  calculate stresses at material point using external soil models
    !
    !*********************************************************************        
    
    implicit none
        
    integer(INTEGER_TYPE), intent(in) :: IDpt ! global integration/material point number
    integer(INTEGER_TYPE), intent(in) :: IDel ! global element number
    ! B-matrix at the considered integration point (here only used if ApplyObjectiveStress=TRUE)
    real(REAL_TYPE), dimension(NVECTOR, ELEMENTNODES), intent(in) :: BMatrix    
    integer(INTEGER_TYPE), intent(in) :: IEntityID ! entity ID (here only used if ApplyObjectiveStress=TRUE)
    
    ! local variables
    character(len=80) :: cmname
    integer(INTEGER_TYPE) :: I ! counter
    integer(INTEGER_TYPE) :: IDset ! ID of material parameter set
    integer(INTEGER_TYPE) :: ntens ! Dimension of stress vector to pass to ESM 
    integer(INTEGER_TYPE), parameter :: nAddVar = 12
    real(REAL_TYPE), dimension(NPROPERTIES) :: props ! array of material properties
    real(REAL_TYPE), dimension(nAddVar) :: AdditionalVar
    real(REAL_TYPE), dimension(MatParams(MaterialIDArray(IDpt))%UMATDimension) :: Stress, StrainIncr ! stress and strain increment in integration/material point
    real(REAL_TYPE), dimension(NTENSOR) :: Sig0, StressIncr, StressPrinc, TempStrainIncr, TempStrainIncrPrevious
    real(REAL_TYPE), dimension(NSTATEVAR) :: StateVar ! state parameters in integration/material
    real(REAL_TYPE) :: Eunloading, PlasticMultiplier
    character(len=64) :: NameModel ! name of the constitutive model
    logical :: IsUndrEffectiveStress
    real(REAL_TYPE) :: DSigWP ! Change of water pressure at integration point 
    real(REAL_TYPE) :: DSigGP ! Change of gas pressure at integration point 
    real(REAL_TYPE) :: Bulk ! Bulk modulus
    real(REAL_TYPE) :: DEpsVol ! Incremental volumetric strain (water)

    pointer (p, ESM)             
          
    ! get constitutive model in integration/material point
    IDset = MaterialIDArray(IDpt) ! is the material number stored in $$MATERIAL_INDEX in the GOM-file
    NameModel = MatParams(IDset)%MaterialModel ! name of constitutive model as specified in GOM-file
    ntens = MatParams(IDset)%UMATDimension  ! 2D or 3D formulation of the External soil model   
          
    ! get strain increments in integration/material point
    TempStrainIncr = GetEpsStep(Particles(IDpt)) ! incremental strain vector assigned to point
    
    if (CalParams%ApplyImplicitQuasiStatic) then
        if (CalParams%ImplicitIntegration%Iteration > 1) then
            do I = 1, NTENSOR
                TempStrainIncrPrevious(I) = GetEpsStepPreviousI(Particles(IDpt), I)
            end do
            
            TempStrainIncr = TempStrainIncr - TempStrainIncrPrevious
            
        end if
    end if
        
    StrainIncr = 0.0

    do I=1, NTENSOR
    StrainIncr(I) = StrainIncr(I) + TempStrainIncr(I)
    enddo 
        
    DEpsVol = StrainIncr(1) + StrainIncr(2) + StrainIncr(3) ! volumetric strain, valid for 2D and 3D
          
    IsUndrEffectiveStress = &
    !code version 2016 and previous
    ((CalParams%ApplyEffectiveStressAnalysis.and.(trim(MatParams(IDSet)%MaterialType)=='2-phase')) .or. &
    !code version 2017.1 and following
    (trim(MatParams(IDSet)%MaterialType)==SATURATED_SOIL_UNDRAINED_EFFECTIVE))
          
    ! initalise water pressure (only needed for undrained analyses)
    DSigWP = 0.0d0
    DSigGP = 0.0d0
    ! for effective stress analysis
    if (IsUndrEffectiveStress) then
        if (Particles(IDpt)%Porosity > 0.0) then
        Bulk = Particles(IDpt)%BulkWater / Particles(IDpt)%Porosity ! kN/m2
        DSigWP = Bulk * DEpsVol
        else
        DSigWP = 0.0
        end if
        call AssignWatandGasPressureToGlobalArray(IDpt, DSigWP, DSigGP)
    end if ! effective stress analysis
          
    ! get stresses in integration/material point      
    do I = 1, NTENSOR
    Sig0(I) = SigmaEff0Array(IDpt, I) ! get initial stress of current step assigned to point 
    end do
    Stress=0.0
    do I=1, NTENSOR
        Stress(I) = Stress(I) + Sig0(I) !To use 3D UMAT also for 2D formulations
    enddo 
          
    ! initialise state variables (only for very first time and load step)

    !!! Commented out on 08/07/23 - reason: Update from the 2023 release

    ! if ((CalParams%IStep == 1).and.(CalParams%TimeStep == 1)) then
    ! StateVar = MatParams(IDset)%ESM_Statvar_in ! array of NSTATEVAR initial value of the State Variables for the external soil model in UMAT/VUMAT format
    ! else 
    StateVar = ESMstatevArray(IDpt,:)
    ! end if 
    
    ! Add the emailed code here    
    if (IsUndrEffectiveStress) then
        Particles(IDPt)%WaterPressure = Particles(IDPt)%WaterPressure + DSigWP
    end if 
    
    
    !call AssignWatandGasPressureToGlobalArray(IDpt, DSigWP, DSigGP) !Note that the subroutine checks Cavitation Threshold & Gas Pressure
          
    !get values of variables of interest for UMAT model
    AdditionalVar(1) = Particles(IDPt)%Porosity
    AdditionalVar(2) = Particles(IDPt)%WaterPressure
    AdditionalVar(3) = Particles(IDPt)%WaterPressure0 
    AdditionalVar(4) = Particles(IDPt)%GasPressure
    AdditionalVar(5) = Particles(IDPt)%GasPressure0
    AdditionalVar(6) = Particles(IDPt)%DegreeSaturation
    AdditionalVar(7) = CalParams%TotalRealTime
    AdditionalVar(8) = CalParams%OverallRealTime
    AdditionalVar(9) = CalParams%TimeIncrement
    AdditionalVar(10) = CalParams%IStep
    AdditionalVar(11) = CalParams%TimeStep
    AdditionalVar(12) = Particles(IDpt)%BulkWater
          
    ! get name of DLL
    cmname = MatParams(IDSet)%SoilModelDLL
    ! get material properties  
    props = MatParams(IDSet)%ESM_Solid ! This sets props to be the material parameters for an external soil model
         
    if (trim(NameModel)//char(0) == trim('linear_elasticity')//char(0)) then
    props(1) = Particles(IDpt)%ShearModulus ! shear modulus, G
    cmname = UMAT_LINEAR_ELASTICITY
    elseif (trim(NameModel)//char(0) == trim(ESM_MOHR_COULOMB_STANDARD)//char(0)) then
    props(1) = Particles(IDpt)%ShearModulus ! shear modulus, G
    props(2) = MatParams(IDSet)%PoissonRatio 
    props(3) = SIN(MatParams(IDSet)%FrictionAngle*(Pi/180.0)) 
    props(4) = Particles(IDpt)%CohesionCosPhi 
    props(5) = SIN(MatParams(IDSet)%DilatancyAngle*(Pi/180.0))
    props(6) = MatParams(IDSet)%TensileStrength
    cmname = UMAT_MOHR_COULOMB_STANDARD
    endif    
    
    ! initialise UMAT
    
    
    if (NameModel == ESM_ARB_Model_MohrCoulombStrainSoftening) then
        ! Abdel's (Unreleased) Mohr-Coulomb Strain Softening Formulation
        call ESM_MohrCoulombStrainSoftening(IDpt, IDel, IDset, Stress, Eunloading, PlasticMultiplier, StrainIncr, NSTATEVAR, StateVar,&
                                            nAddVar, AdditionalVar,cmname, NPROPERTIES, props, CalParams%NumberOfPhases, ntens)

    else if (NameModel == ESM_NON_ASSOC_MOHR_COULOMB) then 
        ! Luis's (2022) Viscoplastic Mohr-Coulomb
        call ESM_VPSS_MC(IDpt, IDel, IDset, Stress, Eunloading, PlasticMultiplier, StrainIncr, NSTATEVAR, StateVar,&
                        nAddVar, AdditionalVar,cmname, NPROPERTIES, props, CalParams%NumberOfPhases, ntens)
    
    else 
        p = GetProcAddress(MatParams(IDSet)%SoilModelDLLHandle, "ESM"C) ! Pointing to the ESM .dll 

        call ESM(IDpt, IDel, IDset, Stress, Eunloading, PlasticMultiplier, StrainIncr, NSTATEVAR, StateVar,&
                    nAddVar, AdditionalVar,cmname, NPROPERTIES, props, CalParams%NumberOfPhases, ntens)
    end if
    
    ! save unloading stiffness in Particles array  
    Particles(IDpt)%ESM_UnloadingStiffness = Eunloading
                 
    if (IsUndrEffectiveStress) then
        Particles(IDpt)%BulkWater = AdditionalVar(12)
    end if
    
    call SetIPL(IDpt, IDel, int(PlasticMultiplier))


    
    ! to use objective stress definition
    if (CalParams%ApplyObjectiveStress) then ! Consider large deformation terms
    call Hill(IdEl, ELEMENTNODES, IncrementalDisplacementSoil(1:Counters%N, IEntityID),  &
                     ReducedDof, ElementConnectivities, BMatrix, Sig0(1:NTENSOR), Stress(1:NTENSOR), DEpsVol)
    end if ! objective stress            
            
    ! write new stresses to global array
    do I=1, NTENSOR
        StressIncr(I) = Stress(I) - Sig0(I)
    enddo             
                               
    ! save updated state variables and in Particles array
    ESMstatevArray(IDpt,:) = StateVar
          
    call CalculatePrincipalStresses(IDpt, Stress(1:NTENSOR), StressPrinc)
    call AssignStressStrainToGlobalArrayESM(IDpt, NTENSOR, StressIncr, StressPrinc, StrainIncr)

    ! write plasticity state to global array
    !  call SetIPL(IDpt, IDel, int(StateVar(50)))
    if (CalParams%ApplyBulkViscosityDamping) then
    RateVolStrain(IDEl) = DEpsVol / CalParams%TimeIncrement
    call CalculateViscousDamping_interface(IDpt, IDEl)
    end if  

end subroutine StressSolid

subroutine CalculateViscousDamping_interface(ParticleID, IEl)
    !**********************************************************************
    !
    !> Computes a pressure term introducing bulk viscosity damping to the equation of motion.
    !>
    !! \param[in] ParticleID ID of considered material point.
    !! \param[in] IEl ID of element of the considered material point.
    !! \param[in] DilationalWaveSpeed Current wave speed computed for the considered material point.
    !
    !*********************************************************************

    implicit none
    integer(INTEGER_TYPE), intent(in) :: ParticleID, IEl
    real(REAL_TYPE) :: ViscousDampingPressure = 0.0
    real(REAL_TYPE) :: Density = 0.0
    real(REAL_TYPE) :: ElementLMinLocal = 0.0
    real(REAL_TYPE) :: RateVolStrainLocal = 0.0
    real(REAL_TYPE) :: MaterialIndex = 0.0
    real(REAL_TYPE) :: DilationalWaveSpeed = 0.0
    logical :: IsUndrEffectiveStress
    if (.not.CalParams%ApplyBulkViscosityDamping) return
    MaterialIndex = MaterialIDArray(ParticleID)

     IsUndrEffectiveStress = &
        !code version 2016 and previous
        ((CalParams%ApplyEffectiveStressAnalysis.and.(trim(MatParams(MaterialIndex)%MaterialType)=='2-phase')) .or. &
        !code version 2017.1 and following
        (trim(MatParams(MaterialIndex)%MaterialType)==SATURATED_SOIL_UNDRAINED_EFFECTIVE))
    
    !if (CalParams%ApplyEffectiveStressAnalysis
     if (IsUndrEffectiveStress &
    .or.((CalParams%NumberOfPhases==2).or.(CalParams%NumberOfPhases==3))) then
      Density = MatParams(MaterialIndex)%DensityMixture / 1000.0
    else
      Density = (1 - MatParams(MaterialIndex)%InitialPorosity) * MatParams(MaterialIndex)%DensitySolid / 1000.0
    end if
    ElementLMinLocal = ElementLMin(IEl)
    RateVolStrainLocal = RateVolStrain(IEl)

    call GetWaveSpeed(ParticleID, DilationalWaveSpeed)
    ViscousDampingPressure = CalParams%BulkViscosityDamping1 *  &
      Density * DilationalWaveSpeed * ElementLMinLocal * RateVolStrainLocal
    if ((RateVolStrainLocal < 0.0).and.(CalParams%BulkViscosityDamping2 > 0.0)) then
      ViscousDampingPressure = ViscousDampingPressure + &
        Density * (CalParams%BulkViscosityDamping2 * ElementLMinLocal * RateVolStrainLocal)**2
    end if
    Particles(ParticleID)%DBulkViscousPressure = ViscousDampingPressure
end subroutine CalculateViscousDamping_interface
        
 !*************************************************************************************************************************************
 !*******************************Everything below this line has been added when compared to the original (Jonathan Moore)**************
 !*************************************************************************************************************************************

Subroutine ESM_VPSS_MC(NPT,NOEL,IDSET,STRESS,EUNLOADING,PLASTICMULTIPLIER, DSTRAN,NSTATEV,STATEV,NADDVAR,ADDITIONALVAR,CMNAME,NPROPS,PROPS,NUMBEROFPHASES,NTENS)
        
    !DEC$ ATTRIBUTES DLLEXPORT, ALIAS:"ESM" :: ESM
    !---Updating these variables to follow the updated fortran format 
    implicit none 
    !!!!! These intents may need to be changed to intent(in) !!!!!
    integer, intent(in) :: nstatev, ntens, nprops
    integer :: naddvar
    integer, intent(in) :: npt
    integer, intent(in) :: noel
    integer, intent(in) :: idset
    double precision, dimension(ntens), intent(inout) :: stress
    double precision, intent(inout) :: eunloading, plasticmultiplier
    double precision, dimension(ntens), intent(inout) :: dstran
    double precision, dimension(nstatev), intent(inout) :: statev
    double precision, dimension(naddvar) :: additionalvar 
    character(len=80) :: cmname !This is a parameter do not set intent
    double precision, dimension(nprops), intent(inout) :: props
    integer:: numberofphases    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    !---Local variables required in standard UMAT
    integer :: IStep, TimeStep
    double precision :: sse, spd, scd ! specific elastic strain energy, plastic dissipation, creep dissipation
    double precision :: rpl ! only for fully coupled thermal analysis: volumetric heat generation
    double precision, dimension(:), allocatable :: ddsddt ! only for fully coupled thermal analysis: variation of stress increment due to temperature
    double precision, dimension(:), allocatable :: drplde ! only for fully coupled thermal analysis: variation of volumetric heat generation due to strain increment
    double precision                            :: drpldt ! only for fully coupled thermal analysis: variation of volumetric heat generation due to temperature
    double precision, dimension(:), allocatable :: stran
    double precision, dimension(:), allocatable :: time
    double precision, dimension(:), allocatable :: predef
    double precision, dimension(:), allocatable :: dpred    
    double precision, dimension(:), allocatable :: coords
    double precision, dimension(:,:), allocatable :: ddsdde ! Jacobian matrix of the constitutive model (tangent stiffness matrix in case of MC)
    double precision, dimension(:,:), allocatable :: drot
    double precision, dimension(:,:), allocatable :: dfgrd0
    double precision, dimension(:,:), allocatable :: dfgrd1
    double precision :: pnewdt, dtime, temp, dtemp, celent
    double precision :: Value ! auxiliary variable holding any real valued number
    double precision :: Porosity, WaterPressure, WaterPressure0, GasPressure, GasPressure0, DegreeSaturation  
    integer          :: ndi, nshr, layer, kspt, kstep, kinc     
    
    !---Local variables defned by the user
    ! e.g. integer :: var_local	  
    !---User can define here additional variables     
    
    allocate( ddsddt(ntens), drplde(ntens), stran(ntens), time(2), predef(1), dpred(1), &
          coords(3), ddsdde(ntens,ntens), drot(3,3), dfgrd0(3,3), dfgrd1(3,3) )
    
    !Initialization
    Eunloading = 0.0
    PlasticMultiplier = 0.0

    !Rename additional variables
    Porosity = AdditionalVar(1)
    WaterPressure = AdditionalVar(2)
    WaterPressure0 = AdditionalVar(3)
    GasPressure = AdditionalVar(4)
    GasPressure0 = AdditionalVar(5)
    DegreeSaturation = AdditionalVar(6)
    time(1) = AdditionalVar(7)   !TotalRealTime
    time(2) = AdditionalVar(8)   !OverallTotalTime
    dtime = AdditionalVar(9)     !TimeIncrement
    IStep = AdditionalVar(10)    
    TimeStep = AdditionalVar(11)   !Note: Very first time and load step: Istep=1 and TimeStep=1   
    
    !Call the UMAT
    call UMAT_VPSS_MC(stress, statev, ddsdde, sse, spd, scd, rpl, ddsddt, drplde, drpldt, stran, dstran, time, dtime, temp, &
            dtemp, predef, dpred, cmname, ndi, nshr, ntens, nstatev, props, nprops, coords, pnewdt, celent, dfgrd0, &      
            dfgrd1, noel, npt, layer, kspt, kstep, kinc)
    
    !---Definition of Eunloading -> required to define the max time step
          Eunloading = max(ddsdde(1,1),ddsdde(2,2),ddsdde(3,3))
    !---Always define this value to run the simulation

    ! PlasticMultiplier can be given as an output because plastic points can be plotted as a result
    return
end subroutine ESM_VPSS_MC
  
  
Subroutine UMAT_VPSS_MC(STRESS, STATEV, DDSDDE, SSE, SPD, SCD, RPL, DDSDDT, DRPLDE, &
			                DRPLDT, STRAN, DSTRAN, TIME, DTIME, TEMP, DTEMP, PREDEF, DPRED,&
			                CMNAME, NDI, NSHR, NTENS, &	NSTATEV, PROPS, NPROPS, COORDS, &
			                PNEWDT, CELENT, DFGRD0, DFGRD1, NOEL, NPT, &
			                LAYER, KSPT, KSTEP, KINC)
    !DEC$ ATTRIBUTES DLLEXPORT, ALIAS:"UMAT" :: UMAT
    implicit none

	!Defining inputs
    ! Arguments:
    !          I/O  Type
    !  PROPS    I   R()  : List with model parameters
    !  DSTRAN   I   R()  : Strain increment
    !  DTIME    I   R()  : Time increment
    !  DDSDDE   O   R(,) : Material stiffness matrix
    !  STRESS  I/O  R()  : stresses
    !  STATEV  I/O  R()  : state variables
    !
    
    integer, intent(in) :: NSTATEV, NPROPS, NPT
    integer :: NTENS
    integer, intent(in) :: NOEL 
	real (Real_Type), dimension(NTENS), intent(inout) :: STRESS
    real(Real_Type), dimension(NSTATEV), intent(inout) :: STATEV
    real(Real_Type), dimension(NTENS, NTENS), intent(inout) :: DDSDDE
	real(Real_Type), intent(in) :: sse, spd, scd ! specific elastic strain energy, plastic dissipation, creep dissipation
    real(Real_Type), intent(in) :: rpl ! only for fully coupled thermal analysis: volumetric heat generation
	real(Real_Type), dimension(NTENS), intent(in) :: DDSDDT
    real(Real_Type), dimension(NTENS), intent(in) :: DRPLDE
	real(Real_Type)                               :: DRPLDT
    real(Real_Type), dimension(NTENS), intent(in) :: STRAN
    real(Real_Type), dimension(NTENS), intent(in) :: DSTRAN
    real(Real_Type), dimension(2), intent(in) :: TIME
    real(Real_Type), dimension(1), intent(in) :: PREDEF
    real(Real_Type), dimension(1), intent(in) :: DPRED
    real(Real_Type), dimension(NPROPS), intent(in) :: PROPS
    real(Real_Type), dimension(3), intent(in) :: COORDS
    real(Real_Type), dimension(3,3), intent(in) :: DFGRD0
    real(Real_Type), dimension(3,3), intent(in) :: DFGRD1
    REAL(Real_Type), intent(in) :: PNEWDT,  TEMP, DTEMP, CELENT
	double precision, intent(in) :: DTIME
    character(len = 80), intent(in):: CMNAME
	integer, intent(in) :: NDI, NSHR, LAYER, KSPT, KSTEP, KINC  
  
    ! Local variables:
    !
    !  DE        : Linear Elastic constitutive matrix
    !  dSig	     : Stress increment vector
    !  Sig	     : Stress vector
    !  dEpsE     : Elastic strain increment vector
    !  dEpsP     : Plastic strain increment vector
    !  dEps      : Total strain increment vector
    !  EpsE      : Elastic strain vector
    !  EpsP      : Plastic strain vector
    !  Eps       : Total strain vector
    !  EpsRate	 : Total strain rate tensor
    !
    integer     			:: N_S, N_i
    real(Real_Type), dimension(6,6) :: DE
    real(Real_Type), dimension(6)   :: Sig
    real(Real_Type), dimension(6)   :: dEpsP
    real(Real_Type), dimension(6)   :: EpsP, ERate
    real(Real_Type), dimension(6) :: dEpsE, EpsE, dEps, Eps
    logical :: switch_smooth, switch_original, switch_yield
    real(Real_Type) :: G_0, enu, eM_tc, eN, D_min, eh, alpha_G, alpha_K, alpha_D, D_part, G_s, RefERate !SSMC props local variables, (props)
    real(Real_Type) :: G, bk, eta_y, DP, eI_coeff, Sum_rate  !SSMC state variables (statv)
    double precision :: Error_yield_1, Error_yield_2, Error_Euler_max, Error_Yield_last, Error_Yield_max, num_OS_Iterations
    double precision :: F1, F2, bK_0
    integer :: i
        
        
    !
	! Viscoplastic NAMC model with smoothed outer surface in pi plane
	!
	! Contents of PROPS(10) NAMC with HSR
	!  1 : G_0				shear modulus
	!  2 : enu				Poisson's ratio
	!  3 : eM_tc            Critical stress ratio for triaxial compression
    !  4 : eN               Nova's vol coupling coefficient
	!  5 : D_min			Minimum dilation
	!  6 : eh				hardening parameter
	!  7 : alpha_G			Shear modulus viscosity factor
	!  8 : alpha_K			Bulk modulus viscosity factor
	!  9 : alpha_D			dilation viscosity
	! 10 : D_part			Particle diameter
	! 11 : G_s				Specific gravity      
	! 12 : RefERate			Reference strain rate
	! 13 : Switch_smooth	Boolean switch for activating strain rate smoothing
	! 14 : N_S				Degree of smoothening
    ! 15 : switch_original	Changes from Wang's to Olzak&Perzyna consistency 
	G_0			= PROPS(1)         ! shear modulus
    enu			= PROPS(2)         ! Poisson's ratio
    eM_tc		= PROPS(3)         ! Critical stress ratio
    eN          = PROPS(4)         ! Nova's vol coupling coefficient
    D_min       = PROPS(5)         ! Minimum dilation
    eh			= PROPS(6)         ! hardening parameter
	alpha_G	    = PROPS(7)	       ! Shear modulus viscosity factor
	alpha_K	    = PROPS(8)	       ! Bulk modulus viscosity factor
	if ((alpha_K==0.0d0).and.(alpha_G>0.0d0)) then ! in case user uses 0 it will be computed from nu
		alpha_K= 2.5*alpha_G
	endif
	alpha_D		= PROPS(9)    ! dilation angle viscosity
	if ((alpha_D==0.0d0).and.(alpha_G>0.0d0)) then ! in case user uses 0 it will be computed from nu
		alpha_D= alpha_G
    endif
    D_part		= PROPS(10)			! Associated particle diameter, value in mm
    G_s  		= PROPS(11)			! Specific gravity
	RefERate	= PROPS(12)		    ! reference strain rate
	if (RefERate==0.0d0) then
		RefERate=2.5e-5
	endif
	call dbltobool(PROPS(13), switch_smooth)  ! switch for activating strain rate smoothening
	N_S=PROPS(14)							  ! Degree of smoothening
    call dbltobool(PROPS(15), switch_original)! (1 for Wang, 0 for Olzak&Perzyna)
    
    G			    = STATEV(1)			         ! shear modulus
	bK			    = STATEV(2)			         ! bulk modulus
	eta_y		    = STATEV(3)                  ! friction ratio
    Dp		        = STATEV(4)                  ! Dilation
    eI_coeff		= STATEV(5)                  ! Inertial coefficient
    call dbltobool(STATEV(6),switch_yield)       ! Point is yielding
	do i=1,6
		EpsP(i)     =STATEV(6+i)				 !Plastic strain               
	end do
	N_i             =STATEV(13)			         ! Current number of strain rate sums
	SUM_rate        =STATEV(14)			         ! Current sum of strain rates

    num_OS_Iterations = 0
    
	!_____Error control state parameters__________________________________________________
	Error_yield_1=0.0d0                                                                   !
    Error_yield_2=0.0d0                                                                   !
    Error_Euler_max=0.0d0															      !
	Error_Yield_last=0.0d0															      !
	Error_Yield_max=0.0d0																  !
	!							                                                          !
	!_____________________________________________________________________________________!		
	if (DTIME==0.0d0) then
		ERate= 0.0d0    ! Current strain rate
	else
		ERate= (1/DTIME)*DSTRAN ! Current strain rate
    end if
    
	bK_0= 2*G_0*(1+ENU)/(3*(1-2*ENU))

	!***********************************************************************************
	!Call the refined modified Euler algorithm
	call NAMC_HSR(NOEL, G_0, enu, eM_tc, eN, D_min, eh, alpha_G, alpha_K, alpha_D, D_part, G_s,&
                 switch_smooth, RefERate, N_S, switch_original,&
				 G, bK, eta_y, Dp, EpsP, eI_coeff, switch_yield, N_i, SUM_rate,&
				 DSTRAN, STRESS, Sig, Erate, DTIME,&
				 Error_yield_1, Error_yield_2, Error_Euler_max, Error_Yield_last, &
				 Error_Yield_max)
	!************************************************************************************
	!************************************************************************************

	!Stress and state variables updating
	Do i=1,NTENS
      STRESS(i) = Sig(i)
    End Do
    STATEV(1) = G
    STATEV(2) = bK
    STATEV(3) = eta_y
	STATEV(4) = Dp
    STATEV(5) =	eI_coeff	
	STATEV(6) = logic2dbl(switch_yield)
	do i=1,6
		STATEV(6+i) = EpsP(i)
	end do			
	STATEV(13)=N_i
	STATEV(14)=SUM_rate

	!_____Error control state parameters__________________________________________________
	!     Comment if not wanted                                                           !
	!_____________________________________________________________________________________!	
	STATEV(28)=Error_yield_1
    STATEV(29)=Error_yield_2
    STATEV(30)=Error_Euler_max
	STATEV(31)=Error_Yield_last
	STATEV(32)=Error_Yield_max
    
	!************************************************************************************
	!************************************************************************************
	!Tangent stiffness matrix to be returned done by elastic stiffness
    F1  = bK+(4*G/3)
    F2  = bK-(2*G/3)
    DDSDDE = 0.0
    DDSDDE(1:3,1:3) = F2
    DDSDDE(1,1) = F1
    DDSDDE(2,2) = F1
    DDSDDE(3,3) = F1
    DDSDDE(4,4) = G
    DDSDDE(5,5) = G
    DDSDDE(6,6) = G
	!*************************************************************************************
	!End of UMAT
	Return
end subroutine UMAT_VPSS_MC
    
!_____________________________________________________________________________________________    
!##    ##    ###    ##     ##  ######     ##     ##  ######  ########  
!###   ##   ## ##   ###   ### ##    ##    ##     ## ##    ## ##     ## 
!####  ##  ##   ##  #### #### ##          ##     ## ##       ##     ## 
!## ## ## ##     ## ## ### ## ##          #########  ######  ########  
!##  #### ######### ##     ## ##          ##     ##       ## ##   ##   
!##   ### ##     ## ##     ## ##    ##    ##     ## ##    ## ##    ##  
!##    ## ##     ## ##     ##  ######     ##     ##  ######  ##     ##  

subroutine NAMC_HSR(NOEL, G_0, nu, M_tc, N, D_min, h, alpha_G, alpha_K, alpha_D, D_part, G_s,&
                     switch_smooth, RefERate, N_S, switch_original, &
					 G, K, eta_y, Dp, EpsP, I_coeff, switch_yield, N_i, SUM_rate,&
					 dEps, Sig_0, Sig, Erate, DTIME,&
					 Error_yield_1, Error_yield_2, Error_Euler_max, Error_Yield_last, &
					 Error_Yield_max)
    
    !*********************************************************************************
    !
    ! Elasto-plastic constitutive model with strain softening and strain rate effects,
	! Non associated MOHR-COULOMB 
    ! Explicit MODIFIED EULER INTEGRATION SCHEME with automatic error control.
    ! Final correction of the yield surface drift (END OF STEP CORRECTION).
    !
    !*********************************************************************************
	implicit none
    
	!input variables
	integer, intent(in) :: NOEL !Global ID of Gauss point or particle
    integer, intent(in) :: N_S
    double precision, intent(in)::  G_0, nu, M_tc, N, D_min, alpha_G, alpha_K, &
									alpha_D, D_part, G_s, &
									RefERate, DTIME
    double precision, dimension(6), intent(in)::  dEps
    double precision, dimension(6), intent(inout):: Sig_0
    
    logical, intent(in):: switch_smooth, switch_original           
    
	!output variables
    integer, intent(inout):: N_i
    
    double precision, intent(inout):: h, G, K, eta_y, Dp, I_coeff, &
                                      Error_yield_1, Error_yield_2, Error_Euler_max,&
									  Error_Yield_last, Error_Yield_max, SUM_rate
    double precision, dimension(6), intent(inout):: Sig, Erate, EpsP
    
    logical, intent(inout):: switch_yield

	!local variables
    integer:: counter, MAXITER, SubStepping_MaxIter
    double precision:: p, q, theta, M, I_act, I_0, Gu, Ku, eta_yu, Du, dI, &
                       G1, K1, eta_y1, Dp1, G2, K2, eta_y2, Dp2, &
                       p_t, q_t, dI_t, dI_TT, I_TT, dD1, dD2
    double precision:: epsq_rate, epsq_p, eps_v
    double precision:: F0, FT, alpha
    double precision:: FTOL, STOL, DTmin , LTOL, R_TT, qR, RTOL
    double precision:: dummyvar(3), D1, D2
    double precision:: DE(6,6), dSig_el(6), Sig_t(6), dEps_t(6), dEps_TT(6), &
                       Sig1(6), Sig2(6), dEpsp1(6), dEpsp2(6), dEpsp(6), &
                       dSig1(6), dSig2(6), Epspt(6)
    double precision:: T, DT
    
    logical:: ApplyStrainRateUpdates=.false., &!Check if the strain rate path crosses the reference line
              IsUnloading, & !If true material unloads (not elastic unloading)
              Failed=.false. !If rel residual is too high failed is true
    !print *, NOEL
	!______________________________________________________________________________
    ! Error tolerances
    FTOL=1.0e-8		!Yield surface tolerance
    STOL=1.0e-3		!Relative error tolerance
    DTmin=1.0e-9	!Minimum pseudo-time increment, originally Dtmin = 1.0e-9
	LTOL=0.01d0		!Tolerance for elastic unloading	  
	MAXITER=20		!Max. number of iterations
    RTOL = STOL * 1.0e-1
    !______________________________________________________________________________
	!______________________________________________________________________________
    ! Initialization of error trackers
    Error_yield_1   =0.0	!Drift after first approximation
    Error_yield_2   =0.0	!Drift after second approximation
	Error_Euler_max =0.0    !Max relative residual
	Error_Yield_max =0.0    !max drift after averaging
	Error_Yield_last=0.0    !max abs drift after drift correction		
    !______________________________________________________________________________
    !Initialization of state variables
    call Get_invariants(Sig_0, p, q, theta)
	call Get_strain_invariants(EpsP, eps_v, epsq_p)!plastic deviatoric strain
    I_0=RefERate
    !call Get_I_coeff(D_part, G_s, -100.0, RefERate, I_0)!Reference inertial coefficient
    call Get_M(M_tc, theta, M)!Get M
    if (G==0.0d0) then
        G=G_0
    endif
    if (K==0.0d0) then
        K=2*G_0*(1+nu)/(3*(1-2*nu))
    endif        
    !if (I_coeff==0.0d0) then
    !    call Get_I_coeff(D_part, G_s, p, RefERate, I_coeff)
    !endif

    ! Dp testing: This should be the first time that the dilatancy can be updated
    if (Dp==0.0d0) then ! Since plastic deviatoric strain epsq_p will be zero in elastic Dp will be be zero 
        call Get_Dp(h, D_min, I_coeff, I_0, epsq_p, alpha_D, ApplyStrainRateUpdates, Dp)
    endif
    
    !print *, eta_y
    if (eta_y==0.0d0) then
        call Get_invariants(dEps, dummyvar(1), dummyvar(2), theta)
        call Get_M(M_tc, theta, M)
        eta_y=M-Dp*(1.0-N)
    endif
    !print *, eta_y

    !_____________________________________________________________________________
    !Evaluate yield function at initial stress-state variable
    
    call YieldFunction(q, p, eta_y, F0)
    
    !_____________________________________________________________________________
    !Compute the current deviatoric strain rate
    
    call TwoNormTensor_strain(Erate,6,dummyvar(1))!Norm of the strain rate
    TrackStrainRate = dummyvar(1) !Track the strain rate
    
	 !Compute a smoothed strain rate if switch_smooth is true
	if (switch_smooth) then
		N_i=N_i+1
        !if (Noel ==1) then
        !    print *, N_i
        !endif
		if (N_i<N_s) then !not enough values
			Sum_rate=sum_rate+dummyvar(1) !accumulate the strain rate
			dummyvar(2)=Sum_rate/n_i !takes the average
		else !enough values
			Sum_rate=Sum_rate*(1.0-1.0/N_s) !approximate sum without one term
			Sum_rate=Sum_rate+dummyvar(1)
			dummyvar(2)=Sum_rate/N_s !averaged strain rate
        endif
        if (dummyvar(1)==0.0d0) then !in case in first step this value is zero
            Erate=0.0d0
        else
            Erate=(dummyvar(2)/dummyvar(1))*Erate !corrected strain rate tensor
        endif
        call TwoNormTensor_strain(Erate,6,dummyvar(1)) ! recalculate the two norm so the updated erate value can be tracked
        TrackSmoothStrainRate = dummyvar(1) !track the smoothed strain rate
    endif
            
    !________________________________________________________________________________
    
    !Update state variables due to HSR
    !Store state variables
    Gu=G
    Ku=K
    eta_yu=eta_y
    Du=Dp        
    
    call Get_strain_invariants(Erate, dummyvar(1), epsq_rate)! deviatoric strain rate
    TrackShearStrainRate = epsq_rate ! Track shear strain rate
    call Get_I_coeff(D_part, G_s, p, epsq_rate, I_act)!actual inertial coefficient
    TrackInertialCoefficient = I_act ! Track the intertial coefficient
    !I_act = 0 ! Test to see if this is the reason the models get different results for no updates and D_part = 0
    dI=I_act-I_coeff !change in inertial coefficient
    call check4crossing(I_coeff, I_act, dI, I_0, ApplyStrainRateUpdates) !Check if update needed    
    
    if (ApplyStrainRateUpdates) then
        ! Just using this as a way to catch if apply strain rate updates is flipped to true
        DeformCateg = 10.0
    endif
        
    ! Dp Testing: Second time that dilatancy can be updated.
    if (applystrainrateupdates) then !update
        !h=h*(i_act/i_0)**alpha_g
        call update_gk(g_0, nu, i_act, i_0, alpha_g, alpha_k, gu, ku) !updates the elastic properties
        call get_dp(h, d_min, i_act, i_0, epsq_p, alpha_d, applystrainrateupdates, du) !updates the dilation
        eta_yu=m-du*(1.0-n) !updates eta
    endif
    !_________________________________________________________________________________
    
	! Fill elastic material matrix
	D1  = Ku+(4*Gu/3)
	D2  = Ku-(2*Gu/3)
	DE  = 0.0
	DE(1:3,1:3) = D2
	DE(1,1) = D1
	DE(2,2) = D1
	DE(3,3) = D1
	DE(4,4) = Gu
	DE(5,5) = Gu
	DE(6,6) = Gu
           
    !Get elastic predictor stress
    call MatVec(DE,6,dEps,6,dSig_el)
    call AddVec(Sig_0, dSig_el, 1.0, 1.0, 6, Sig_t)
    !Get new invariant stresses
    call Get_invariants(Sig_t, p_t, q_t, theta)
    !Evaluate yield function
    call YieldFunction(q_t, p_t, eta_yu, FT)
    
    !Track variables
    !Trackq_t = q_t
    !Trackp_t = p_t
    TrackEta_y = eta_y ! Track eta_y for plotting

    !___________________________________________________________________________
    !Now check elastic loading, unloading
	if (FT<-FTOL) then !Elastic behavior
       !Update state parameters
        G=Gu
        K=Ku
        eta_y=eta_yu
        Dp=Du
        I_coeff=I_act

        !Update stress
        Sig=Sig_t
        switch_yield=.false.
        
        !Track Loading Variables
        DeformCateg = 100.0
        !print *, FT
        !TrackFT = FT
    !___________________________________________________________________________
	!*************************  Plastic behavior  ****************************** 
    !***************************************************************************
    !***************************************************************************
        
    !___________________________________________________________________________
    !***********Checking surface bisection and viscoplastic unloading***********
    !===========================================================================
    else !plastic behavior            
        !if (F0<-FTOL) then !Elasto-plastic transition                
        !    call Newton_Raphson(G, K, eta_y, Dp, G_0, nu, M, M_tc, N, D_min, h, D_part, &
								!G_s, epsq_p, I_coeff, I_act, I_0, alpha_K, alpha_G, alpha_D,&
								!Gu, Ku, eta_yu, Du, dEps, MAXITER, FTOL,&
								!F0, Sig_0, alpha)
        !    !Track Deformation category
        !    DeformCateg = 200.0
        !    
        !
        !else!pure plastic deformations or unloading
        !    call Check_Unloading(M_tc, eta_y, eta_yu, dI, Sig_0, dSig_el, LTOL, IsUnloading)  !Determines if is unloading path
        !    if (IsUnloading) then !Find new alpha
        !        call Newton_Raphson(G, K, eta_y, Dp, G_0, nu, M, M_tc,N, D_min, h, D_part, &
								!	G_s, epsq_p, I_coeff, I_act, I_0, alpha_K, alpha_G, alpha_D,&
								!	Gu, Ku, eta_yu, Du, dEps, MAXITER, FTOL,&
								!	F0, Sig_0, alpha)
        !        !Track Deformation category
        !        DeformCateg = 300.0
        !    
        !    else !Pure plasticity
        !        alpha=0.0d0
        !        !Track Deformation category
        !        DeformCateg = 400.0
        !    endif                
        !endif 
    
        ! Now that the updated state parameters have been found plug integrate the 
        ! Need to make dEpsP a state variable (Done)         
        call Ortiz_Simo_Integration(G_0, nu, M_tc, M, N, D_min, h, G, K, eta_y, Dp, &
                                I_0, I_coeff, dI, alpha_G, alpha_K, alpha_D, Sig_0, EpsP, dEps, &
                                FTOL, NOEL)
        Sig = Sig_0 ! Update the stresses
        
        ! Track varaibles for output
        TrackDp = Dp

    endif
    !Track variables
    TrackShearModulus = Gu
    call Get_invariants(Sig, p, q, theta) ! Recalculate invariants so that they can be stored
    Trackq_t = q
    Trackp_t = p
end subroutine NAMC_HSR
!*******************************************************************************************

!___________________________________________________________________________________________
  !     ##    ## ######## ##      ## ########  #######  ##    ## 
  !     ###   ## ##       ##  ##  ##    ##    ##     ## ###   ## 
  !     ####  ## ##       ##  ##  ##    ##    ##     ## ####  ## 
  !     ## ## ## ######   ##  ##  ##    ##    ##     ## ## ## ## 
  !     ##  #### ##       ##  ##  ##    ##    ##     ## ##  #### 
  !     ##   ### ##       ##  ##  ##    ##    ##     ## ##   ### 
  !     ##    ## ########  ###  ###     ##     #######  ##    ## 
  !########     ###    ########  ##     ##  ######   #######  ##    ## 
  !##     ##   ## ##   ##     ## ##     ## ##    ## ##     ## ###   ## 
  !##     ##  ##   ##  ##     ## ##     ## ##       ##     ## ####  ## 
  !########  ##     ## ########  #########  ######  ##     ## ## ## ## 
  !##   ##   ######### ##        ##     ##       ## ##     ## ##  #### 
  !##    ##  ##     ## ##        ##     ## ##    ## ##     ## ##   ### 
  !##     ## ##     ## ##        ##     ##  ######   #######  ##    ##

subroutine Newton_Raphson(G, K, eta_y, Dp, G_0, nu, M, M_tc, No, D_min, h, D_par, &
						  G_s, eps_q, I_coeff, I_act, I_0, k_K, k_G, k_D,&
                          Gu, Ku, eta_yu, Du, dEps, MAXITER, FTOL, &
						  F0, Sig_0, alpha)
    !**************************************************************************************
    !  This subroutine determine the elastic proportion with consistent state and elastic *
    !  parameters                                                                         *
    !**************************************************************************************
    implicit none
    !input
    integer, intent(in):: MAXITER
    double precision, intent(in):: G_0, nu, M, M_tc, No, D_min, h, D_par, G_s, &
                                eps_q, k_G, k_K, k_D, I_act, I_0, FTOL
    double precision, dimension(6), intent(in):: dEps
    !output
    double precision, intent(inout):: G, K, eta_y, Dp, I_coeff, &
                                    Gu, Ku, eta_yu, Du, F0, Sig_0(6)
    double precision, intent(out):: alpha
    !local variables
    logical:: ApplyStrainRateUpdates
    integer:: n, i
    double precision:: FT, dG, dK, deta, dD, dI, &
                    dI_alpha, I_alpha, &
                    p_alpha, q_alpha, dummyvar, F_prime
    double precision:: dEps_alpha(6), dSig_alpha(6), Sig_alpha(6), &
                    n_vec(6), L, dSigdAlpha(6) 
    double precision:: D1, D2, DE(6,6), dDdG(6,6), dDdK(6,6), aux(6,6)
    
    !Initialize parameters
    FT=1000
    alpha=0.5d0
    n=0
    !Determine Changes in state parameters
    dG=Gu-G
    dK=Ku-K
    deta=eta_yu-eta_y
    dD=Du-Dp
    dI=I_act-I_coeff    
    do while ((abs(FT)>=FTOL).and.(n<=MAXITER))
        !Store variables
        Gu=G
        Ku=K
        eta_yu=eta_y
        Du=Dp
        F0=FT
        
        n=n+1
        !___________________________________________________________________
        !Compute trial strains and stresses
        dEps_alpha=alpha*dEps
        dI_alpha=alpha*dI
        I_alpha=I_coeff+dI_alpha
        !___________________________________________________________________
        !Evaluate rate crossing
        call check4crossing(I_coeff, I_alpha, dI_alpha,I_0, ApplyStrainRateUpdates)
        if (ApplyStrainRateUpdates) then !Update parameters
            call Update_GK(G_0, nu, I_alpha, I_0, k_G, k_K, Gu, Ku)
            call Get_Dp(h, D_min, I_alpha, I_0, eps_q, k_D, ApplyStrainRateUpdates, Du)
            eta_yu=M-Du*(1.0-No)
        endif
        !___________________________________________________________________
        !Update trial stress 
        !Ensemble elastic matrix
        D1  = Ku+(4*Gu/3)
        D2  = Ku-(2*Gu/3)
        DE  = 0.0
        DE(1:3,1:3) = D2
        DE(1,1) = D1
        DE(2,2) = D1
        DE(3,3) = D1
        DE(4,4) = Gu
        DE(5,5) = Gu
        DE(6,6) = Gu  
        call MatVec(DE,6,dEps_alpha,6,dSig_alpha)
        call AddVec(Sig_0, dSig_alpha, 1.0, 1.0, 6, Sig_alpha)
        !___________________________________________________________
        !Evaluate yield function
        call Get_invariants(Sig_alpha, p_alpha, q_alpha, dummyvar)
        call YieldFunction(q_alpha, p_alpha, eta_yu, FT)
        
        !Evaluate n=dF/dSig and L=dF/dXs
        call Get_dF_to_dSigma(M_tc, eta_yu, Sig_alpha, n_vec)!dF/dSig
        L=-p_alpha*(1.0-No) !L=dF/dXs        
        !____________________________________________________________
        !Evaluate F'(alpha)
        !Elastic matrix derivatives
        dDdK=0.0
        dDdK(1:3,1:3) = 1.0
        dDdG=0.0
        dDdG(1:3,1:3)=-2./3.
        dDdG(1,1)=4./3.
        dDdG(2,2)=4./3.
        dDdG(3,3)=4./3.
        dDdG(4,4)=1.0
        dDdG(5,5)=1.0
        dDdG(5,5)=1.0
        aux=alpha* (dK*dDdK+dG*dDdG)
        aux=DE+aux
        call MatVec(aux, 6, dEps, 6, dSigdAlpha)
        F_prime=0.0d0
        do i=1,6 !dot product
            F_prime=F_prime+n_vec(i)*dSigdAlpha(i)
        enddo
        F_prime=F_prime+L*dD!F'(alpha)
        !______________________________________________________________
        !Update alpha
        alpha=alpha-FT/F_prime    
        if ((alpha>1.0d0).or.(alpha<0.0d0)) then
            alpha=0.0d0
        endif  
    end do
    !Update variables to return
    G=Gu
    K=Ku
    eta_y=eta_yu
    Dp=Du
    I_coeff=I_alpha
    Sig_0=Sig_alpha    
end subroutine Newton_Raphson
!*********************************************************************************
!*********************************************************************************
!*********************************************************************************
!*********************************************************************************
!*********************************************************************************
!*********************************************************************************
!*********************************************************************************
!*********************************************************************************
!*********************************************************************************
!*********************************************************************************
subroutine Ortiz_Simo_Integration(G_0, nu, M_tc, M, No, D_min, h, G, K, eta_y, Dp, &
                                I_0, I, dI, k_G, k_K, k_D, Sig, EpsP, dEps, &
                                FTOL, NOEL)

    !------------------------Function Information------------------------
    ! Ortiz and Simo (1986) integration scheme (Algorithm on pg. 6)
    ! Returns:
    ! Increment of plastic strain, plastic strain, and change of state parameters
    !----------------------End Function Information----------------------

    implicit none
    ! Variable Info
        !!!------Begin External variables----!!!
        ! G_0: Initial Shear Modulus
        ! nu: Poisson ratio
        ! M_tc: Critical stress ratio for triaxial loading
        ! M: Current critical stress ratio
        ! No: Nova's volumetric coupling coefficient
        ! D_min: Minimum dilatancy
        ! h: Dilatancy hardening parameter
        ! G: Current Shear modulus
        ! K: Current Bulk Modulus
        ! eta_y: Current Stress ratio
        ! Dp: Current dilatancy
        ! I_0: Reference interial coefficient
        ! I: Current inertial coefficient
        ! dI: Increment of the inertial coeff
        ! k_G: Shear modulus viscosity coeff
        ! k_K: Bulk modulus viscosity coeff
        ! k_D: Dilatancy viscosity coeff
        ! Sig: Input stress
        ! EpsP: Plastic strain
        ! dEps: Increment of strain
        ! dEpsp: Increment of Plastic Strain
        ! FTOL: Tolerance for initial predictor
    !!!------End External variables-----!!!

    !!!------Begin new output variables------!!!
    !!!---------End Returned calculated variables-----!!!

    !--------------Input variables--------------!
    
    ! Input scalar values
    double precision, intent(in):: G_0, nu, M_tc, No, D_min, h, &
                                k_G, k_K, k_D, FTOL
    integer, intent(in) :: NOEL
    ! Input vector values
    double precision, dimension (6), intent(in):: dEps

    !-------------End Input Variables-----------!

    !--------------Output Variables-------------!
    ! In/Out scalar values
    double precision, intent(inout):: G, K, eta_y, Dp, I_0, I, dI, M

    ! In/Out Vector values
    double precision, dimension(6), intent(inout):: Sig, EpsP !, dEpsP

    ! Out scalar values
    !double precision, intent(out):: 

    ! Out vector values
    ! double precision, dimension(6), intent(out):: dSig

    !-------------End Output Variables--------!

    !-------------local Variables-------------!
    ! Local scalar values
    double precision:: I_f, F, p, q, epsv_p, epsq_p, eta_yu, Du, Mu, dummyVal, a_Dot_m, L, H_term, &
                       D1, D2, b, dLambda, dD

    logical:: ApplyStrainRateUpdate = .false.
    integer:: counter, MaxIter

    ! Local vector values
    double precision, dimension(6):: dEpsE, dummyVec, dEpsPu, EpsPu, Sigu, &
                                    m_vec, n_vec, DE_m, a, dSig

    ! Local matrix values
    double precision, dimension(6,6):: DE

    !------------End local Variables----------!

    ! Store variables for updating
    Sigu = Sig
    EpsPu = EpsP
    ! dEpsPu = dEpsP
    eta_yu = eta_y
    Du = Dp
    Mu = M
    
    ! Apply strain rate updates
    I_f=I+dI
    call check4crossing(I,  I_f, dI, I_0, ApplyStrainRateUpdate)
    
    call Get_strain_invariants(EpsPu, epsv_p, epsq_p)

    if (ApplyStrainRateUpdate) then !Update parameters
        call Update_GK(G_0, nu, I_f, I_0, k_G, k_K, G, K)
        call Get_Dp(h, D_min, I_f, I_0, epsq_p, k_D, ApplyStrainRateUpdate, Du)
        eta_yu = Mu-du*(1.0 * No)
    endif

    ! Turn off strain rate affects
    ApplyStrainRateUpdate = .false.

    !--------------------Compute elastic predictor---------------------------!
    dEpsE = dEps ! This is Luis's option for predicting the elastic strain

    !-----------Ensemble elastic matrix------------------!
    D1  = K+(4*G/3)
    D2  = K-(2*G/3)
    DE  = 0.0
    DE(1:3,1:3) = D2
    DE(1,1) = D1
    DE(2,2) = D1
    DE(3,3) = D1
    DE(4,4) = G
    DE(5,5) = G
    DE(6,6) = G 
    !----------End Ensemble elastic matrix--------------!

    ! Calc the stress predictor
    call MatVec(DE, 6, dEpsE, 6, dSig)

    !Update the stresses
    Sigu = Sigu + dSig

    !-----------------End Compute elastic predictor------------------------!
    ! Since all strain is assumed to be elastic no plastic strain increment

    !-------------------Begin Yielding Check--------------------------!

    ! Compute stress invariants
    call Get_invariants(Sigu, p, q, dummyVal)
    
    ! M = M_tc*(1 + 0.25(cos(1.5 * theta + 0.25 *pi))
    call Get_M(M_tc, dummyVal, Mu)
    
    eta_yu = Mu - Du * (1.0 -No)

    ! Compute the value of the yield Function
    call YieldFunction(q, p, eta_yu, F)

    if (abs(F) < FTOL) then
        ! Prediction is correct, stress and strain values can be updated and returned
        
        ! Store that it passed on the first iteration
        num_OS_Iterations = 1
        
        ! Update Sig, EpsP, dEpsP, eta_y, Dp
        Sig = Sigu
        EpsP = EpsPu
        ! dEpsP = dEpsPu
        eta_y = eta_yu
        Dp = Du

        ! Exit out of the Subroutine and return values
        return
    endif
    ! else
        ! Predicted stress too large, iterations are needed
    !-------------------End Yielding Check--------------------------!


    !-----------------------Begin Plastic Descent-----------------------!
    
    ! Set num_OS_Iterations to next index
    num_OS_Iterations = 2
    
    ! Set Maximum number of iterations
    MaxIter = 100000
    counter = 0
    
    call Get_dD_to_dEpsP(D_min, h, I_0, k_D, epsq_p, epsv_p, &
                                dEpsPu, I, ApplyStrainRateUpdate, a) !a=dD/dEpsP

    ! Compute stress invariants
    call Get_invariants(Sigu, p, q, dummyVal)
    
    do while (abs(F) >= FTOL .and. counter <= MaxIter) 
        !---------------------Begin Compute derivatives--------------------------!
        call Get_dF_to_dSigma(Mu, eta_yu, Sigu, n_vec) !n=dF/dSig
        call Get_dP_to_dSigma(Du, Sigu, m_vec) !m=dP/dSig
        L = -p * (1-No) !dF/Xs = dF/dDp = Xi in Ortiz & Simo
        
        !-----------------------End Compute derivatives--------------------------!

        !----------------Compute Denominator of dLambda--------------------------!
        ! Denominator = v:D:r - Xi.h == n_vec:DE:m_vec - H_term
        
        ! Compute DE.m_vec
        call MatVec(DE, 6, m_vec, 6,  DE_m)

        ! Compute n_vec.DE.m_vec
        call DotProduct_2(n_vec, DE_m, 6, dummyVal)

        ! compute Xi.h = H_term = L.a.m
        call DotProduct_2(a, m_vec, 6, a_Dot_m)
        H_term = L * a_Dot_m

        !--------------End Compute Denominator of dLambda-----------------------!

        ! Compute the change in the plastic potential (lambda)
        dLambda = F/(dummyVal - H_term) 

        ! Compute the stress update
        Sigu = Sigu - dLambda * DE_m

        ! Compute stress invariants
        call Get_invariants(Sigu, p, q, dummyVal)
                    
        ! Update M
        call Get_M(M_tc, dummyVal, Mu)
        
        ! calc dLambda * m (increment of plastic strain)
        dEpsPu = dLambda * m_vec
        
        ! Accumulate plastic strain
        EpsPu = EpsPu + dEpsPu
    
        ! Calc strain invariants
        call Get_strain_invariants(EpsPu, epsv_p, epsq_p)
        
        ! Calc a = dD/dEpsP
        call Get_dD_to_dEpsP(D_min, h, I_0, k_D, epsq_p, epsv_p, &
                                dEpsPu, I, ApplyStrainRateUpdate, a) !a=dD/dEpsP
        
        ! Update dilatancy
        dD = 0.0
        !call DotProduct_2(a, dEpsPu, 6, dD) !plastic hard/softening
        Du = Du + dD
        
        ! Update eta_y
        eta_yu = Mu - Du * (1.0 -No)

        ! Calc the yield function
        call YieldFunction(q, p, eta_yu, F)
        
        ! Update Counter to end while loop
        Counter = Counter + 1
        
    end do
    
    ! Calc the change in Sig
     dSig = Sigu - Sig ! Don't think I need this because I'm not substepping

    !Calc dEpsP
    dEpsPu = EpsPu - EpsP

    ! Update Sig, dEpsP, EpsP, eta_y, Dp
    Sig = Sigu
    ! dEpsP = dEpsPu
    EpsP = EpsPu
    eta_y = eta_yu
    
    dD = Du -Dp
    Dp = Du
    M = Mu
    ! Store Yield function value at end of iteration
    TrackFVal = F
    
    ! Track 2-Norm of plastic strain EpsP
    call TwoNormTensor(EpsP, 6, NormEpsP)
end subroutine Ortiz_Simo_Integration
!______________________________________________________________________________________    
!######## ##     ## ##       ######## ########      ######  
!##       ##     ## ##       ##       ##     ##    ##    ## 
!##       ##     ## ##       ##       ##     ##    ##       
!######   ##     ## ##       ######   ########      ######  
!##       ##     ## ##       ##       ##   ##            ## 
!##       ##     ## ##       ##       ##    ##     ##    ## 
!########  #######  ######## ######## ##     ##     ######  
!   ###    ##        ######    #######  ########  #### ######## ##     ## ##     ## 
!  ## ##   ##       ##    ##  ##     ## ##     ##  ##     ##    ##     ## ###   ### 
! ##   ##  ##       ##        ##     ## ##     ##  ##     ##    ##     ## #### #### 
!##     ## ##       ##   #### ##     ## ########   ##     ##    ######### ## ### ## 
!######### ##       ##    ##  ##     ## ##   ##    ##     ##    ##     ## ##     ## 
!##     ## ##       ##    ##  ##     ## ##    ##   ##     ##    ##     ## ##     ## 
!##     ## ########  ######    #######  ##     ## ####    ##    ##     ## ##     ## 

subroutine Euler_Algorithm(G_0, nu, M_tc, M, No,  D_min, h, Dpart, Gs,&
                                G, K, eta_y, Dp, &	
                            erate, I_0, I, dI, k_G, k_K, k_D, dtime, DT, &
                            switch_original, &
                            Sig, EpsP, dEps, dD, dEpsp, dSig)
    !************************************************************************
    ! Euler's algorithm (finite difference) integration						*
    ! Returns:																*
    ! Increment of stresses and change of state parameters                  *                           
    !************************************************************************  
    implicit none
    !input
    logical, intent(in):: switch_original
    double precision, intent(in):: G_0, nu, M_tc, M, No, D_min, h, Dpart, Gs,&
                                erate(6), I_0, k_G, k_K, k_D, dtime, DT,&
                                dEps(6)
    double precision, intent(inout):: I, dI
    !output
    double precision, intent(inout):: G, K, eta_y, Dp, EpsP(6) , Sig(6)
    double precision, intent(out):: dD, dSig(6), dEpsp(6)
    !local variables
    logical:: ApplyStrainRateUpdate
    double precision:: I_f, DE(6,6), D1, D2, dSig_el(6), num, lambda, &
                    n_vec(6), p, q, dummyvar, L, b, m_vec(6), a(6), &
                    epsq_p, epsv_p, Hard, den, dummyvec(6), Hvp, erate_q

    !________________________________________________________________________
    !Evaluate for strain rate updating
    I_f=I+dI
    call check4crossing(I,  I_f, dI, I_0, ApplyStrainRateUpdate)
    
    if (ApplyStrainRateUpdate) then
        call Update_GK(G_0, nu, I, I_0, k_G, k_K, G, K)
    endif
    !________________________________________________________________________
    
    !________________________________________________________________________
    !Compute invariants and derivatives
    call Get_invariants(Sig, p, q, dummyvar)
    call Get_strain_invariants(EpsP, epsv_p, epsq_p)
    call Get_dF_to_dSigma(M_tc, eta_y, Sig, n_vec) !n=dFdSig
    call Get_dP_to_dSigma(Dp, Sig, m_vec) !m=dP/dSig
    L = -p*(1.0-No) !L=dF/dXs
    call Get_dD_to_dEpsP(D_min, h, I_0, k_D, epsq_p, epsv_p, &
                        EpsP, I, ApplyStrainRateUpdate, a) !a=dD/dEpsq^p   
    call Get_dD_to_dI(D_min, h, I_0, k_D, epsq_p, I, b) !b=dXs/dI in place of dXs/dErate
    !________________________________________________________________________
    
    !________________________________________________________________________
    !Compute elastic predictor	
    !Ensemble elastic matrix    
    D1  = K+(4*G/3)
    D2  = K-(2*G/3)
    DE  = 0.0
    DE(1:3,1:3) = D2
    DE(1,1) = D1
    DE(2,2) = D1
    DE(3,3) = D1
    DE(4,4) = G
    DE(5,5) = G
    DE(6,6) = G 
    
    !Calculate the predictor
    call MatVec(DE, 6, dEps, 6, dSig_el)
    !________________________________________________________________________
    
    !________________________________________________________________________
    !compute numerator
    call DotProduct_2(dSig_el, n_vec, 6, num)! dSig_el . n_vec  
    
    if ((ApplyStrainRateUpdate).and.(.not.switch_original)) then !Dashpot visc.
        num=num+L*b*dI
    endif
    !_______________________________________________________________________
    
    !_______________________________________________________________________
    !Compute denominator
    
    !Compute hardening modulus H=-L.a.m
    call DotProduct_2(a, m_vec, 6, dummyvar)
    Hard=-L*dummyvar !hardening modulus
    
    !compute Den=n.D.m +H
    call MatVec(DE, 6, m_vec, 6, dummyvec)
    call DotProduct_2(n_vec, dummyvec, 6, den)
    den=den+Hard
    
    if ((ApplyStrainRateUpdate).and.(switch_original)) then !original
    !compute viscous hardening modulus Hvp=-L.bv.m/(dt*DT)
        call Get_strain_invariants(erate,dummyvar,erate_q)
        call Get_dEpsq_to_dEps(erate_q,erate,dummyvec)
        dummyvec=b*Dpart*sqrt(Gs/abs(p))*dummyvec
        call DotProduct_2(dummyvec, m_vec, 6, dummyvar)
        Hvp=-L*dummyvar/(dtime*DT)
        den=den+Hvp
    endif
    !_____________________________________________________________________
    
    !_____________________________________________________________________
    !compute viscoplastic mult. lambda, dEpsp, dSig, and dD
    
    lambda=num/den! viscoelastic multiplier    
    dEpsp=lambda*m_vec! Flow rule
    dummyvec=dEps-dEpsp
    call MatVec(DE, 6, dummyvec, 6, dSig) !stress increment
    dD=0.0
    call DotProduct_2(a, dEpsp, 6, dD)!plastic hard/softening
    dD=dD+b*dI !rate hard/softening
    !_____________________________________________________________________
    
    !_____________________________________________________________________
    !Update EpsP, Sig, D, and eta_y
    EpsP=EpsP+dEpsp
    Sig=Sig+dSig
    Dp=Dp+dD
    eta_y=M-Dp*(1.0-No)
    !____________________________________________________________________    
end subroutine Euler_Algorithm
!*************************************************************************
    
    
!________________________________________________________________________   
    
! ######  ######## ########  ########  ######   ######  
!##    ##    ##    ##     ## ##       ##    ## ##    ## 
!##          ##    ##     ## ##       ##       ##       
! ######     ##    ########  ######    ######   ######  
!      ##    ##    ##   ##   ##             ##       ## 
!##    ##    ##    ##    ##  ##       ##    ## ##    ## 
! ######     ##    ##     ## ########  ######   ######  
!########  ########  #### ######## ######## 
!##     ## ##     ##  ##  ##          ##    
!##     ## ##     ##  ##  ##          ##    
!##     ## ########   ##  ######      ##    
!##     ## ##   ##    ##  ##          ##    
!##     ## ##    ##   ##  ##          ##    
!########  ##     ## #### ##          ## 
    
subroutine Stress_Drift_Correction(G_0, nu, M_tc, M, No, D_min, h, Dpart, Gs, &
                                   G, K, eta_y, Dp, &	
                                   erate, I_0, I_f, I, dI, k_G, k_K, k_D, dtime, DT, &
                                   switch_original, MAXITER, F0, FTOL, &
                                   Sig, EpsP, dEps)
    !************************************************************************
    ! Corrects any drift caused during the Mod. Euler's procedure       	*
    ! Returns corrected stresses and state parameters                       *                 
    !************************************************************************
    implicit none
    !input
    logical, intent(in):: switch_original
    integer, intent(in):: MAXITER
    double precision, intent(in):: G_0, nu, M_tc, M, No, D_min, h, Dpart, Gs,&
                                erate(6), I_0, k_G, k_K, k_D, dtime, DT,&
                                FTOL, dEps(6)
    double precision, intent(inout):: I, I_f, dI
    !output
    double precision, intent(inout):: G, K, eta_y, Dp, EpsP(6) , Sig(6), F0
    !local variables
    logical:: ApplyStrainRateUpdate
    integer:: n
    double precision::epsq_p, epsv_p, n_vec(6), m_vec(6), L, a(6), b,&
                    DE(6,6), D1, D2, Den,  Hard, Hvp, dlambda, &
                    dummyvec(6), dummyvar, FC, &
                    Du, eta_yu, Sigu(6), dD, dSig(6), dEpsp(6),Epspu(6), &
                    p, q, erate_q
    !________________________________________________________________________
    !Evaluate for strain rate updating    
    call check4crossing(I, I_f, dI, I_0, ApplyStrainRateUpdate)
    
    if (ApplyStrainRateUpdate) then
        call Update_GK(G_0, nu, I_f, I_0, k_G, k_K, G, K)
    endif
    !________________________________________________________________________
    
    !________________________________________________________________________
    !Initialize values
    n=0 !counter
    FC=F0    
    !________________________________________________________________________
    
    do while ((abs(FC)>FTOL).and.(n<MAXITER))
        !___________________________________________________________________
        !update values
        F0=FC
        n=n+1
        !___________________________________________________________________
        
        !________________________________________________________________________
        !Compute invariants and derivatives
        call Get_invariants(Sig, p, q, dummyvar)
        call Get_strain_invariants(EpsP, epsv_p, epsq_p)
        call Get_dF_to_dSigma(M_tc, eta_y, Sig, n_vec) !n=dFdSig
        call Get_dP_to_dSigma(Dp, Sig, m_vec) !m=dP/dSig
        L=-p*(1.0-No) !L=dF/dXs
        call Get_dD_to_dEpsP(D_min, h, I_0, k_D, epsq_p, epsv_p, &
                            EpsP, I_f, ApplyStrainRateUpdate, a) !a=dD/dEpsq^p
        call Get_dD_to_dI(D_min, h, I_0, k_D, epsq_p, I_f, b) !b=dXs/dI in place of dXs/dErate
        !________________________________________________________________________
        
        !___________________________________________________________________
        !compute denominator
        
        !Ensemble elastic matrix    
        D1  = K+(4*G/3)
        D2  = K-(2*G/3)
        DE  = 0.0
        DE(1:3,1:3) = D2
        DE(1,1) = D1
        DE(2,2) = D1
        DE(3,3) = D1
        DE(4,4) = G
        DE(5,5) = G
        DE(6,6) = G 
        
        !Compute hardening modulus H=-L.a.m
        call DotProduct_2(a, m_vec, 6, dummyvar)
        Hard=-L*dummyvar !hardening modulus
    
        !compute Den=n.D.m +H
        call MatVec(DE, 6, m_vec, 6, dummyvec)
        call DotProduct_2(n_vec, dummyvec, 6, den)
        den=den+Hard
    
        if ((ApplyStrainRateUpdate).and.(switch_original)) then !original
        !compute viscous hardening modulus Hvp=-L.bdI/dErate.m/(dt*DT)
            call Get_strain_invariants(erate,dummyvar,erate_q)
            call Get_dEpsq_to_dEps(erate_q,erate,dummyvec)
            dummyvec=b*Dpart*sqrt(Gs/abs(p))*dummyvec
            call DotProduct_2(dummyvec, m_vec, 6, dummyvar)
            Hvp=-L*dummyvar/(dtime*DT)
            den=den+Hvp
        endif        
        !___________________________________________________________________
        
        !___________________________________________________________________
        !Compute DeltaLambda=F0/den and updated Sig, D, and eta_y
        dlambda=F0/den      !error of plastic multiplier
        dEpsp=dlambda*m_vec !error in plastic strain tensor
        Epspu=EpsP-dEpsp
        call MatVec(DE, 6, dEpsp, 6, dSig) !stress drift
        dD=0.0
        call DotProduct_2(a, dEpsp, 6, dD)!plastic hard/softening
        !if (switch_original) dD=dD+b*dI
        
        !update
        Sigu=Sig-dSig
        Du=Dp-dD
        eta_yu=M-Du*(1.0-No)
        !__________________________________________________________________
        
        !__________________________________________________________________
        !Evaluate yield function
        call Get_invariants(Sigu, p, q, dummyvar)
        call YieldFunction(q, p, eta_yu, FC)
        !__________________________________________________________________
        !Evaluate change direction
        if (abs(FC)>abs(F0)) then
            call DotProduct_2(n_vec, n_vec, 6, den)
            dlambda=F0/den
            dEpsp=dlambda*m_vec !error in plastic strain tensor
            Epspu=EpsP-dEpsp
            call MatVec(DE, 6, dEpsp, 6, dSig) !stress drift
            Sigu=Sig-dSig
            Du=Dp
            eta_yu=eta_y
        end if
        !__________________________________________________________________
        
        !__________________________________________________________________
        !Update the state parameters
        EpsP=Epspu
        Sig=Sigu
        Dp=Du
        eta_y=eta_yu
        !__________________________________________________________________
    end do
    
    
end subroutine Stress_Drift_Correction
!******************************************************************************
    
    
    
    
    
!______________________________________________________________________________  
!##     ##    ###    ########  ########  
!##     ##   ## ##   ##     ## ##     ## 
!##     ##  ##   ##  ##     ## ##     ## 
!######### ##     ## ########  ##     ## 
!##     ## ######### ##   ##   ##     ## 
!##     ## ##     ## ##    ##  ##     ## 
!##     ## ##     ## ##     ## ########    
! ######   #######  ######## ########       ###    ##    ## ########  
!##    ## ##     ## ##          ##         ## ##   ###   ## ##     ## 
!##       ##     ## ##          ##        ##   ##  ####  ## ##     ## 
! ######  ##     ## ######      ##       ##     ## ## ## ## ##     ## 
!      ## ##     ## ##          ##       ######### ##  #### ##     ## 
!##    ## ##     ## ##          ##       ##     ## ##   ### ##     ## 
! ######   #######  ##          ##       ##     ## ##    ## ########  
!########  ######## ########  #### ##     ##    ###    ######## #### ##     ## 
!##     ## ##       ##     ##  ##  ##     ##   ## ##      ##     ##  ##     ## 
!##     ## ##       ##     ##  ##  ##     ##  ##   ##     ##     ##  ##     ## 
!##     ## ######   ########   ##  ##     ## ##     ##    ##     ##  ##     ## 
!##     ## ##       ##   ##    ##   ##   ##  #########    ##     ##   ##   ##  
!##     ## ##       ##    ##   ##    ## ##   ##     ##    ##     ##    ## ##   
!########  ######## ##     ## ####    ###    ##     ##    ##    ####    ###    
!########  ######  
!##       ##    ## 
!##       ##       
!######    ######  
!##             ## 
!##       ##    ## 
!########  ###### 
!Derivatives inclosed here
	subroutine Get_Dp(h, D_min, I, I_0, eps_q, k, ApplyRateUpdating, D)
	!*********************************************************************
	! Returns the dilation for current inertial coefficient and dev.     *
	! strain															 *
	!*********************************************************************
    implicit none
    logical, intent(in):: ApplyRateUpdating
    double precision, intent(in):: h, D_min, I, I_0, eps_q, k
    !out
    double precision, intent(out):: D
    !local variables
    double precision:: D_mm
    if (ApplyRateUpdating) then
		D_mm=D_min*(I/I_0)**k !strain/rate hardening
    else
        D_mm=D_min
    endif
    
    D=h*D_mm*eps_q*exp(1.0-h*eps_q) !hardening rule
    end subroutine Get_Dp
    
    subroutine Update_GK(G_0, nu, I, I_0, k_G, k_K, G, K)
	!*********************************************************************
	! Returns updated elastic modulus                                    *
	!																	 *
	!*********************************************************************
    implicit none
    !input
    double precision, intent(in):: G_0, nu, I, I_0, k_G, k_K
    !output
    double precision, intent(out):: G, K
    !local variables
    double precision:: K_0
    G=G_0*(I/I_0)**k_G! updated modulus
    
    K_0=2*G_0*(1+nu)/(3*(1-2*nu))!bulk modulus
    K=K_0*(I/I_0)**k_K! updated modulus
    end subroutine Update_GK
    
    
    subroutine Get_dF_to_dSigma(M_tc, eta_y, Sig, n_vec)
	!************************************************************************
	! Returns the derivative of the yield function with respect to the		*
	! stress tensor 														*
    ! n=dF/dSigma =dF/dp*dp/dSigma+ dF/dq*dq/dSigma +dF/dtheta*dtheta/dSigma*
    ! n is a (1X6) vector													*                           
	!************************************************************************  
    implicit none
    !input
    double precision, intent(in):: M_tc, eta_y, Sig(6)
    !output
    double precision, dimension(6):: n_vec
    !local variables
    double precision:: p, q, theta, pi=2.0*acos(0.0d0), &
                       J2, J3, dJ3dsig(6), dfdtheta, &
					   dpdsig(6), dqdsig(6), dev(6), dev2(6), &
					   TrS2, II(6), dthetadSig(6), COS_3THETA
    !Get the invariants
    call Get_invariants(Sig, p, q, theta)
    !Get dF/dp=eta_y and dF/dq=1
    !Get dF/dtheta
    dfdtheta=0.45*p*M_tc*((cos(1.5d0*theta+0.25d0*pi))**0.2)*sin(1.5d0*theta+0.25d0*pi)
    !___________________________________________________________________________
    !1) Get dp/dSig=1/3 Imat
    dpdsig=0.0d0
    dpdsig(1)=1.0/3.0
    dpdsig(2)=1.0/3.0
    dpdsig(3)=1.0/3.0
    
    !2) Get dq/dsig= 2 *dev/3*q
    dev=Sig
    dev(1)=dev(1)-p
    dev(2)=dev(2)-p
    dev(3)=dev(3)-p
    
    dqdSig=(3.0/(2.0*q))*dev
    
    !3) Get dtheta/dSigma= (1/3cos3theta) d/dsigma((J3/2) * (3/J2)^1.5)
    J2=(q**2)/3.0
    J3=dev(1)*dev(2)*dev(3)-dev(1)*dev(6)**2-dev(2)*dev(4)**2-dev(3)*dev(5)**2+2.0*dev(4)*dev(5)*dev(6)
    !Fill S.S
    dev2(1)=dev(1)**2+dev(4)**2+dev(5)**2
    dev2(2)=dev(2)**2+dev(4)**2+dev(6)**2
    dev2(3)=dev(3)**2+dev(5)**2+dev(6)**2
    dev2(4)=dev(4)*(dev(1)+dev(2))+dev(5)*dev(6)
    dev2(5)=dev(5)*(dev(1)+dev(3))+dev(4)*dev(6)
    dev2(6)=dev(6)*(dev(2)+dev(3))+dev(4)*dev(5)
    !Compute dJ3dSig
    TrS2=dev2(1)+dev2(2)+dev2(3)  
    II=0.0d0!Identity tensor
    II(1)=1.0
    II(2)=1.0
    II(3)=1.0
    dJ3dsig=dev2-(TrS2*II/3.0d0)
    !Compute dtheta/dsig
    
    dthetadSig=dJ3dsig-(1.5*J3/J2)*dev
    COS_3THETA=cos(3.0*theta)
    dthetadSig=(sqrt(3.0)/(2.0*COS_3THETA*J2**1.5))*dthetadSig
    !__________________________________________________________________
    !Get n_vec=dF/dSig
	n_vec=(eta_y*dpdsig)+dqdSig+(dfdtheta*dthetadSig) !n_vec=dF/dSig
    end subroutine Get_dF_to_dSigma
    
    
    subroutine Get_dD_to_dI(D_min, h, I_0, kD, eps_q, I, b)
    !************************************************************************
	! Returns the derivative of the Dilation with respect to the inertial	*
	! coefficient 															*
    ! b=dD/dI																*
    ! b is a scalar															*                           
	!************************************************************************ 
    implicit none
    !input
    double precision, intent(in):: D_min, h, I_0, kD, eps_q, I
    !output
    double precision, intent(out)::b
    !local variables
    
    b=h*D_min*eps_q*exp(1-h*eps_q)*kD*((I/I_0)**(kD-1.0))/I_0
       
    end subroutine Get_dD_to_dI
    
    subroutine Get_dP_to_dSigma(D, Sig, m_vec)
	!************************************************************************
	! Returns the derivative of the plastic potential function with respect *
	! to the stress tensor													*
    ! m=dP/dSigma =dP/dp*dp/dSigma+ dP/dq*dq/dSigma							*
    ! m is a (1X6) vector													*                           
	!************************************************************************  
    implicit none
    !input
    double precision, intent(in):: D, Sig(6)
    !output
    double precision, dimension(6):: m_vec
    !local variables
    double precision:: p, q, theta, pi=2.0*acos(0.0d0), &
                       J2, J3, dJ3dsig(6), dfdtheta, &
					   dpdsig(6), dqdsig(6), dev(6), dev2(6), &
					   TrS2, II(6), dthetadSig(6), COS_3THETA
    !Get the invariants
    call Get_invariants(Sig, p, q, theta)
    !Get dP/dp=-D and dF/dq=1
    !___________________________________________________________________________
    !1) Get dp/dSig=1/3 Imat
    dpdsig=0.0d0
    dpdsig(1)=1.0/3.0
    dpdsig(2)=1.0/3.0
    dpdsig(3)=1.0/3.0
    
    !2) Get dq/dsig= 2 *dev/3*q
    dev=Sig
    dev(1)=dev(1)-p
    dev(2)=dev(2)-p
    dev(3)=dev(3)-p
    
    dqdSig=(3.0/(2.0*q))*dev  
    
    !__________________________________________________________________
    !Get m_vec=dP/dSig
	m_vec=(-D*dpdsig)+dqdSig !m_vec=dP/dSig
    end subroutine Get_dP_to_dSigma
    
    
    
    subroutine Get_dD_to_dEpsP(D_min, h, I_0, k_D, epsq_p, epsv_p, &
						        EpsP, I, ApplyStrainRateUpdate, a)
   	!************************************************************************
	! Returns the derivative of the Dilation with respect to the plastic    *
	! strain																*
    ! a=dXs/dEpsp= dD/dEpsq * dEPsq/dEpsp									*
    ! a is a (1X6) vector													*                           
	!************************************************************************
    implicit none
	!input
	logical, intent(in):: ApplyStrainRateUpdate
	double precision, intent(in):: D_min, h, I_0, k_D, epsq_p, epsv_p, &
                                    EpsP(6), I
	!output
	double precision, intent(out):: a(6)
	!local variables
	double precision:: D, dDdEpsq_p, dev(6),dEpsq_pdEpsp(6) 
	
	!________________________________________________________________________
	!Get dD/dEpsq_p
	if (ApplyStrainRateUpdate) then
		D=D_min*(I/I_0)**k_D
	else
		D=D_min
	end if
	
	dDdEpsq_p=h*D*exp(1.0-h*epsq_p)*(1.0-h*epsq_p)
    !_______________________________________________________________________
    
    !_______________________________________________________________________
    !Get dEpsp_Q/dEpsp= 
	dev=EpsP
    dev(1)=dev(1)-(epsv_p/3.0)
    dev(2)=dev(2)-(epsv_p/3.0)
    dev(3)=dev(3)-(epsv_p/3.0) !deviatoric stress tensor
    
    if (epsq_p>0.0d0) then !in case of zero plastic strain
        dEpsq_pdEpsp=(2.0/(3.0*epsq_p))*dev
    else
        dEpsq_pdEpsp=0.0d0
    endif    
    !_______________________________________________________________________
    
    !_______________________________________________________________________
    !Get a=dXs/dEpsp
    
    a=dDdEpsq_p*dEpsq_pdEpsp
    !______________________________________________________________________    
	end subroutine Get_dD_to_dEpsP
	
	subroutine Get_dEpsq_to_dEps(Epsq, Eps, dEqdEpsq)
	!************************************************************************
	! Returns the derivative of the deviatoric strain with respect to the   *
	! deviatoric strain	tensor						     					*
    ! dEqdEpsq is a (1X6) vector											*                           
	!************************************************************************
	implicit none
	!input
	double precision, intent(in):: Epsq, Eps(6)
	!output
	double precision, intent(out):: dEqdEpsq(6)
	!local variables
	double precision:: evol, dev(6)
	
	evol=Eps(1)+Eps(2)+Eps(3)!vol strain
    
	dev=Eps
	dev(1)=dev(1)-evol/3.0
	dev(2)=dev(2)-evol/3.0
	dev(3)=dev(3)-evol/3.0 !deviatoric strain tensor
	
	if (Epsq>0.0d0) then !in case of zero plastic strain
        dEqdEpsq=(2.0/(3.0*Epsq))*dev
    else
        dEqdEpsq=0.0d0
    endif  	
	end subroutine Get_dEpsq_to_dEps
!**********************************************************************************
	
    
!_________________________________________________________________________________   
!######## ##     ## ##    ##  ######  ######## ####  #######  ##    ##  ######  
!##       ##     ## ###   ## ##    ##    ##     ##  ##     ## ###   ## ##    ## 
!##       ##     ## ####  ## ##          ##     ##  ##     ## ####  ## ##       
!######   ##     ## ## ## ## ##          ##     ##  ##     ## ## ## ##  ######  
!##       ##     ## ##  #### ##          ##     ##  ##     ## ##  ####       ## 
!##       ##     ## ##   ### ##    ##    ##     ##  ##     ## ##   ### ##    ## 
!##        #######  ##    ##  ######     ##    ####  #######  ##    ##  ######  
 

    Subroutine Get_strain_invariants(Eps, Eps_v, Eps_q)
	!*********************************************************************
	! Takes the strain tensor and returns deviatoric and vol. strain     *
	!																	 *
	!********************************************************************* 
    implicit none
    !input
    double precision, dimension(6), intent(in):: Eps
    !output
    double precision, intent(out):: Eps_v, Eps_q
    !local variables
	double precision:: dev(6)
    Eps_v=Eps(1)+Eps(2)+Eps(3)! vol strain
    
	dev=Eps 
	dev(1)=dev(1)-(Eps_v/3.0)
	dev(2)=dev(2)-(Eps_v/3.0)
	dev(3)=dev(3)-(Eps_v/3.0)!deviatoric strain tensor
	
    call TwoNormTensor_strain(dev, 6, Eps_q)
    Eps_q=Eps_q*sqrt(2.0/3.0) ! dev strain
    end subroutine  Get_strain_invariants
    
    
    subroutine Get_invariants(Sig, p, q, theta)
	!*********************************************************************
	! Takes the stress tensor Sig and return invariants p, q, and theta  *
	!																	 *
	!*********************************************************************
	implicit none
	!input variables
	double precision, dimension(6), intent(in):: Sig
	!output variables
	double precision, intent(out)::p, q, theta
	!local variables
	double precision:: dev(6), J2, J3, sin3theta    

	p=(Sig(1)+Sig(2)+Sig(3))/3.0 !mean stress
	dev=Sig
	dev(1)=dev(1)-p !computes deviatoric stress tensor
	dev(2)=dev(2)-p
	dev(3)=dev(3)-p
	
    call TwoNormTensor(dev, 6, J2)
    J2=(J2**2)/2.0 !J_2 invariant
    q=sqrt(3*J2) ! deviatoric stress
    
    !J3 stress invariant
    J3=dev(1)*dev(2)*dev(3)-dev(1)*dev(6)**2-dev(2)*dev(4)**2-dev(3)*dev(5)**2+2.0*dev(4)*dev(5)*dev(6)
    
    !sin3theta
    if (J2>0.0d0) then
       sin3theta=0.5*J3*(3.0/J2)**(1.5d0) 
    else !Assume triaxial compression
       sin3theta=-1.0d0 
    endif
	if (sin3theta<-1.0) sin3theta=-1.0d0
    if (sin3theta>1.0) sin3theta=1.0d0
        
    
    theta=-asin(sin3theta)/3.0d0 !Lode's angle

    end subroutine Get_invariants
    
    
    subroutine Check_Unloading(M_tc, eta_y, eta_yu, dI, Sig, dSig,&
								LTOL, IsUnloading)
	!*********************************************************************
	! Returns true if stress path is viscoplastic unloading              *
	!																	 *
	!*********************************************************************  
    !input
    implicit none
    double precision, intent(in):: M_tc, eta_y, eta_yu, dI, &
                                   Sig(6), dSig(6), LTOL
    !output
    logical, intent(out):: IsUnloading
    !local variables
    double precision:: deta, n_vec(6), n_norm, Sig_norm,&
						dSIg_inner_n, beta, phi
    IsUnloading=.false.
    
    deta=eta_yu-eta_y!change in state parameter    
    call Get_dF_to_dSigma(M_tc, eta_yu, Sig, n_vec)!Normal to surface
    call TwoNormTensor(n_vec, 6, n_norm) !norm of n_vec
    call TwoNormTensor(dSig, 6, Sig_norm) !norm of dSig
    call TensorInnerProduct(dSig, n_vec, 6,dSIg_inner_n) !inner product
    
    beta=acos(deta/(n_norm*Sig_norm))!conical aperture is a plane for inviscid mat.
    phi=acos(dSIg_inner_n/(n_norm*Sig_norm))!angle between stress rate and normal
    
    if (phi-beta>LTOL) IsUnloading=.true. !condition for unloading
    
    end subroutine Check_Unloading
    
    
    subroutine Get_I_coeff(D_part, G_s, p, eps_rate, I)
	!*********************************************************************
	! Returns the inertial coefficient                                   *
	!																	 *
	!********************************************************************* 
    implicit none
    !input
    double precision, intent(in):: D_part, G_s, p, eps_rate
    !output
    double precision, intent(out):: I
    !local variables
    I=D_part*eps_rate*sqrt(G_s/abs(p))
    end subroutine Get_I_coeff
    
    subroutine Get_M(M_tc, theta, M)
	!*********************************************************************
	! Returns M															 *
	!																	 *
	!*********************************************************************
    implicit none
    !in
    double precision, intent(in):: M_tc, theta
    !out
    double precision, intent(out):: M
    !local
    double precision:: COS_VAL, pi=2*acos(0.0d0)
    COS_VAL=cos(1.5*theta+0.25*pi)
    M=M_tc*(1+0.25*COS_VAL**1.2)
    end subroutine Get_M
    
    
    subroutine YieldFunction(q, p, eta_y, F)
	!*********************************************************************
	! Returns the value of the yield function evaluated at q, p , eta    *
	!																	 *
	!*********************************************************************
    implicit none
    !in
    double precision, intent(in):: q, p, eta_y
    !out
    double precision, intent(out):: F
    !local variables
    
    F=q+eta_y*p !sign is due to compression being negative in UMAT
	end subroutine YieldFunction
!***********************************************************************************************
	
!_______________________________________________________________________________________________
!##     ##    ###    ######## ##     ## 
!###   ###   ## ##      ##    ##     ## 
!#### ####  ##   ##     ##    ##     ## 
!## ### ## ##     ##    ##    ######### 
!##     ## #########    ##    ##     ## 
!##     ## ##     ##    ##    ##     ## 
!##     ## ##     ##    ##    ##     ## 
    
    Subroutine TensorInnerProduct(TensorA, TensorB, N, Re)
	!***********************************************************************
	!
	!     Calculate 2NormTensor = sqrt(A:A)
	!
	! I   Tensor  : (Square or vector of dimension N)
	! I   N     :   Number of elements
	! O   2Norm : Resulting norm
	!
	!***********************************************************************
	  Implicit Double Precision (A-H,O-Z)
      Dimension TensorA(N), TensorB(N)
	!***********************************************************************
	  X=N/2
      Re=0.0d0
	  Do I=1,X
		  Re=Re+TensorA(I)*TensorB(I)
	  end Do
	  Do I=X+1,N
		  Re=Re+2*(TensorA(I)*TensorB(I))
	  end do
    end subroutine TensorInnerProduct
    
    
    Subroutine TwoNormTensor(Tensor, N, TwoNorm)
	!***********************************************************************
	!
	!     Calculate 2NormTensor = sqrt(A:A)
	!
	! I   Tensor  : (Square or vector of dimension N)
	! I   N     :   Number of elements
	! O   2Norm : Resulting norm
	!
	!***********************************************************************
	  Implicit Double Precision (A-H,O-Z)
      Dimension Tensor(N)
	!***********************************************************************
	  X=N/2
      TwoNorm=0.0d0
	  Do I=1,X
		  TwoNorm=TwoNorm+Tensor(I)*Tensor(I)
	  end Do
	  Do I=X+1,N
		  TwoNorm=TwoNorm+2*(Tensor(I)*Tensor(I))
	  end do
	  TwoNorm=sqrt(TwoNorm)

    end subroutine TwoNormTensor
        
    Subroutine TwoNormTensor_strain(Tensor, N, TwoNorm)
	!***********************************************************************
	!
	!     Calculate 2NormTensor = sqrt(A:A)
	!
	! I   Tensor  : (Square or vector of dimension N)
	! I   N     :   Number of elements
	! O   2Norm : Resulting norm
	!
	!***********************************************************************
	  Implicit Double Precision (A-H,O-Z)
      Dimension Tensor(N)
	!***********************************************************************
	  X=N/2
      TwoNorm=0.0d0
	  Do I=1,X
		  TwoNorm=TwoNorm+Tensor(I)*Tensor(I)
	  end Do
	  Do I=X+1,N
		  TwoNorm=TwoNorm+0.5*(Tensor(I)*Tensor(I))!The convention in UMAT is to use engineering shear strains
	  end do
	  TwoNorm=sqrt(TwoNorm)

	end subroutine TwoNormTensor_strain
    
    
	Subroutine MatVec(xMat,IM,Vec,N,VecR)
	!***********************************************************************
	!
	!     Calculate VecR = xMat*Vec
	!
	! I   xMat  : (Square) Matrix (IM,*)
	! I   Vec   : Vector
	! I   N     : Number of rows/colums
	! O   VecR  : Resulting vector
	!
	!***********************************************************************
      Implicit Double Precision (A-H,O-Z)
      Dimension xMat(IM,*),Vec(*),VecR(*)
	!***********************************************************************
      Do I=1,N
        X=0
        Do J=1,N
          X=X+xMat(I,J)*Vec(J)
        End Do
        VecR(I)=X
      End Do
      Return
    End Subroutine MatVec
    

	 Subroutine AddVec(Vec1,Vec2,R1,R2,N,VecR)
	!***********************************************************************
	!
	!     Calculate VecR() = R1*Vec1()+R2*Vec2()
	!
	! I   Vec1,
	! I   Vec2  : Vectors
	! I   R1,R2 : Multipliers
	! I   N     : Number of rows
	! O   VecR  : Resulting vector
	!
	!***********************************************************************
      Implicit Double Precision (A-H,O-Z)
      Dimension Vec1(*),Vec2(*),VecR(*)
	!***********************************************************************
      Do I=1,N
        X=R1*Vec1(I)+R2*Vec2(I)
        VecR(I)=X
      End Do
      Return
    End Subroutine AddVec
    
    Subroutine DotProduct_2(VecA, VecB,N, Dp)
	!***********************************************************************
	!
	!     Calculate the dot product of A(Nx1) and B(1xN)
	!
	! I   VecA VecB  : Vectors
	! I   N     :   Dimension
	! O   Dp : Dot product
	!
	!***********************************************************************
	  Implicit Double Precision (A-H,O-Z)
      Dimension VecA(N), VecB(N)
	!***********************************************************************
	  Dp=0.0d0
	  Do I=1,N
			  Dp=Dp+VecA(I)*VecB(I)
	  end do

	end subroutine DotProduct_2

	subroutine dbltobool(A,B)
	!******************************************************************
	! Takes a double which values are either 1.0 or 0.0 and returns a *
	! Boolean
	!******************************************************************
	implicit none
	double precision, intent(in):: A
	logical, intent(out):: B
	if (A<1.0) then
		B=.false.
	else
		B=.true.
	endif
	end subroutine dbltobool

	function logic2dbl(a)
	  logical, intent(in) :: a

	  if (a) then
		logic2dbl = 1.d0
	  else
		logic2dbl = 0.d0
	  end if
	end function logic2dbl
	
	Subroutine check4crossing(IErate0I, IErateI, dErate_eff,RateRef, Apply)
    !******************************************************************
	! determines if strain rate updating must occur                   *
	! Boolean                                                         *
	!******************************************************************
    ! IErate0I: The previous (initial reference strain rate. State parameter pulled from the previous time step
    ! IErateI: The new inertial coefficient for this strain rate
    ! dErate_eff: The increment of strain rate change (This is calculated in this subroutine)
    ! RateRef: Reference strain rate other variables are compared to 
    ! Apply: Boolean keeping track to determine if strain rate updates should be applied
		implicit none
		double precision, intent(inout):: IErate0I, IErateI, dErate_eff
        double precision, intent(in)   :: RateRef
		logical:: cond1, cond2, cond3
		logical, intent(out)::Apply
		  Apply=.false.
          ! If the rate from the last time step is less than or equalt ot the reference rate, update the previous time step value to be the reference rate
		  if(IErate0I<=RateRef) IErate0I=RateRef
          
          ! If the current rate is less than the reference rate than update the current rate to be the reference rate
		  if (IErateI<=RateRef) IErateI=RateRef
          
          ! Cond1 - Checks if the rate has moved from slower than reference to faster than reference on this time step (Rate increased)
		  cond1=(IErate0I==RateRef).and.(IErateI>RateRef)
          
          ! Cond2 - Checks if the rate has moved from greater than the reference to slower than the reference (Rate slowed)
		  cond2=(IErate0I>RateRef).and.(IErateI==RateRef)
		  
          ! Calculate the rate increment
          dErate_eff=IErateI-IErate0I
          
          ! Cond3 - Check if the current and previous value is greater than the reference rate. if they are that means that rate affects should be applied
		  cond3=(IErate0I>RateRef).and.(IErateI>RateRef)
          
          ! Check if any of the conditions are true, if so strain rate effects need to be applied
		  if (cond1.or.cond2.or.cond3) Apply=.true.
	end Subroutine check4crossing

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!    
 !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     
 !███╗░░░███╗░█████╗░░██████╗░██████╗
 !████╗░████║██╔══██╗██╔════╝██╔════╝
 !██╔████╔██║██║░░╚═╝╚█████╗░╚█████╗░
 !██║╚██╔╝██║██║░░██╗░╚═══██╗░╚═══██╗
 !██║░╚═╝░██║╚█████╔╝██████╔╝██████╔╝
 !╚═╝░░░░░╚═╝░╚════╝░╚═════╝░╚═════╝░
      
      
    SUBROUTINE ESM_MohrCoulombStrainSoftening(NPT,NOEL,IDSET,STRESS,EUNLOADING,PLASTICMULTIPLIER,&
                                            DSTRAN,NSTATEV,STATEV,NADDVAR,ADDITIONALVAR,CMNAME,NPROPS,PROPS,NUMBEROFPHASES,NTENS)

        !DEC$ ATTRIBUTES DLLEXPORT, ALIAS:"ESM" :: ESM
        implicit double precision (a-h, o-z) 
        CHARACTER*80 CMNAME         
        DIMENSION STRESS(NTENS),&
        DSTRAN(NTENS),STATEV(NSTATEV),ADDITIONALVAR(NADDVAR),PROPS(NPROPS)
        !NPT(1),NOEL(1),IDSET(1),EUNLOADING(1),PLASTICMULTIPLIER(1),NUMBEROFPHASES(1)

        !---Local variables required in standard UMAT
        integer :: IStep, TimeStep
        double precision, dimension(:), allocatable :: ddsddt ! only for fully coupled thermal analysis: variation of stress increment due to temperature
        double precision, dimension(:), allocatable :: drplde ! only for fully coupled thermal analysis: variation of volumetric heat generation due to strain increment
        double precision, dimension(:), allocatable :: stran
        double precision, dimension(:), allocatable :: time
        double precision, dimension(:), allocatable :: predef
        double precision, dimension(:), allocatable :: dpred    
        double precision, dimension(:), allocatable :: coords
        double precision, dimension(:,:), allocatable :: ddsdde ! Jacobian matrix of the constitutive model (tangent stiffness matrix in case of MC)
        double precision, dimension(:,:), allocatable :: drot
        double precision, dimension(:,:), allocatable :: dfgrd0
        double precision, dimension(:,:), allocatable :: dfgrd1
        double precision :: sse, spd, scd ! specific elastic strain energy, plastic dissipation, creep dissipation
        double precision :: rpl ! only for fully coupled thermal analysis: volumetric heat generation
        double precision :: drpldt ! only for fully coupled thermal analysis: variation of volumetric heat generation due to temperature
        double precision :: pnewdt, dtime, temp, dtemp, celent
        double precision :: Value ! auxiliary variable holding any real valued number
        double precision :: Porosity, WaterPressure, WaterPressure0, GasPressure, GasPressure0, DegreeSaturation  

    
        integer :: ndi, nshr, layer, kspt, kstep, kinc     

        
        
        allocate( ddsddt(ntens), drplde(ntens), stran(ntens), time(2), predef(1), dpred(1),  &
              coords(3), ddsdde(ntens,ntens), drot(3,3), dfgrd0(3,3), dfgrd1(3,3) )
    
        ! Initialization
        Eunloading = 0.0
        PlasticMultiplier = 0.0
     
        !Rename additional variables
        Porosity = AdditionalVar(1)
        WaterPressure = AdditionalVar(2)
        WaterPressure0 = AdditionalVar(3)
        GasPressure = AdditionalVar(4)
        GasPressure0 = AdditionalVar(5)
        DegreeSaturation = AdditionalVar(6)
        time(1) = AdditionalVar(7)   !TotalRealTime
        time(2) = AdditionalVar(8)   !OverallTotalTime
        dtime = AdditionalVar(9)     !TimeIncrement
        IStep = AdditionalVar(10)    
        TimeStep = AdditionalVar(11)   !Note: Very first time and load step: Istep=1 and TimeStep=1   
        
        IDTask = 0
        
        IF((IStep==1).and.(TimeStep==1)) IDTask = 1
     
        IF (IDTask == 1) then ! initialisation of state variables
            STATEV(1)=PROPS(3)
            STATEV(2)=PROPS(5)
            STATEV(3)=PROPS(7)
        END IF ! IDTask = 1
      
        !---Call the UMAT      
        call umat_MohrCoulombStrainSoftening(stress, statev, ddsdde, sse, spd, scd, rpl, ddsddt, drplde, drpldt, stran, dstran, time, dtime, temp, &
           dtemp, predef, dpred, cmname, ndi, nshr, ntens, nstatev, props, nprops, coords, drot, pnewdt, celent, dfgrd0, &
           dfgrd1, noel, npt, layer, kspt, kstep, kinc)

        !---Definition of Eunloading -> required to define the max time step
        Eunloading = max(ddsdde(1,1),ddsdde(2,2),ddsdde(3,3))
        !---Always define this value to run the simulation

        ! PlasticMultiplier can be given as an output because plastic points can be plotted as a result
        return
    end subroutine ESM_MohrCoulombStrainSoftening
      
    !*USER SUBROUTINES
    SUBROUTINE UMAT_MohrCoulombStrainSoftening(STRESS,STATEV,DDSDDE,SSE,SPD,SCD,&
                                                RPL,DDSDDT,DRPLDE,DRPLDT,&
                                                STRAN,DSTRAN,TIME,DTIME,TEMP,DTEMP,PREDEF,DPRED,CMNAME,&
                                                NDI,NSHR,NTENS,NSTATEV,PROPS,NPROPS,COORDS,DROT,PNEWDT,&
                                                CELENT,DFGRD0,DFGRD1,NOEL,NPT,LAYER,KSPT,KSTEP,KINC)

        !DEC$ ATTRIBUTES DLLEXPORT, ALIAS:"UMAT" :: UMAT
        !INCLUDE 'ABA_PARAM.INC'
    
        CHARACTER*80 CMNAME
        DIMENSION STRESS(NTENS),STATEV(NSTATEV),&
        DDSDDE(NTENS,NTENS),DDSDDT(NTENS),DRPLDE(NTENS),&
        STRAN(NTENS),DSTRAN(NTENS),TIME(2),PREDEF(1),DPRED(1),&
        PROPS(NPROPS),COORDS(3),DROT(3,3),DFGRD0(3,3),DFGRD1(3,3)


        ! Arguments:
        !          I/O  Type
        !  PROPS    I   R()  : List with model parameters
        !  DSTRAN   I   R()  : Strain increment
        !  DDSDDE   O   R(,) : Material stiffness matrix
        !  STRESS  I/O  R()  : stresses
        !  STATEV  I/O  R()  : state variables
        !
        !
        !---  Local variables
        !
        Dimension :: DE(6,6), dSig(6), Sig(6), dEpsP(6), EpsP(6)

        !
        ! Mohr-Coulomb Strain Softening model 
        !
        ! Contents of PROPS(9) MCSS
        !  1 : G       shear modulus
        !  2 : ENU     Poisson's ratio
        !  3 : cp      peak cohesion
        !  4 : cr      residual cohesion
        !  5 : phip    peak friction angle 
        !  6 : phir    residual friction angle
        !  7 : psip    peak dilation angle
        !  8 : psir    residual dilation angle 
        !  9 : factor  shape factor
        !
        Rad  = 45d0 / datan(1d0)
        !*
        !* ... start correction routine
        !*
        G      = PROPS(1)         ! shear modulus
        ENU    = PROPS(2)         ! Poisson's ratio
        cp     = PROPS(3)         ! peak cohesion
        cr     = PROPS(4)         ! residual cohesion
        phip   = PROPS(5)/Rad     ! peak friction angle (rad)
        phir   = PROPS(6)/Rad     ! residual friction angle (rad)
        psip   = PROPS(7)/Rad     ! peak dilation angle (rad)
        psir   = PROPS(8)/Rad     ! residual dilation angle (rad)
        factor = PROPS(9)         ! shape factor
        
        c    = STATEV(1)          ! cohesion 
        phi  = STATEV(2)          ! friction angle
        psi  = STATEV(3)          ! dilatancy angle
        Do i = 1,NTENS
        EpsP(i) = STATEV(3+i) 
        end do

        ipl     =   0
        !*
        ! Fill elastic material matrix
        F1  = 2*G*(1-ENU)/(1-2*ENU)
        F2  = 2*G*( ENU )/(1-2*ENU)
        DE  = 0.0
        DE(1:3,1:3) = F2
        DE(1,1) = F1
        DE(2,2) = F1
        DE(3,3) = F1
        DE(4,4) = G
        DE(5,5) = G
        DE(6,6) = G
        !*
        ! elastic stress increment
        Call MatVec( DE, 6, DSTRAN, 6, dSig)
        ! elastic stress
        Call AddVec_MohrCoulombStrainSoftening( STRESS, dSig, 1d0, 1d0, 6, Sig )
        
        call MOHRStrainSoftening(IntGlo,F1,F2,G,cp,cr,phip,phir,psip,psir,factor,c,phi,psi,stress,EpsP,DSTRAN,dEpsP,Sig,IPL)

        !*
        !* ... stress state parameters update
        !*
        Do i=1,NTENS
          STRESS(i) = Sig(i)
        End Do
        
        STATEV(1) = c   
        STATEV(2) = phi  
        STATEV(3) = psi   
        Do i = 1,NTENS
        STATEV(3+i) = EpsP(i)
        end do

        !*
        !* ... Tangent stiffness matrix to be returned (done by elastic stiffness)
        !*
        G       =   PROPS(1)       ! G
        ENU     =   PROPS(2)       ! nu
        F1  = 2*G*(1-ENU)/(1-2*ENU)
        F2  = 2*G*( ENU )/(1-2*ENU)
        DDSDDE = 0.0
        DDSDDE(1:3,1:3) = F2
        DDSDDE(1,1) = F1
        DDSDDE(2,2) = F1
        DDSDDE(3,3) = F1
        DDSDDE(4,4) = G
        DDSDDE(5,5) = G
        DDSDDE(6,6) = G
        !*
        !* ... end UMAT routine
        !*
        Return
    end SUBROUTINE UMAT_MohrCoulombStrainSoftening

!***********************************************************************
      Subroutine MOHRStrainSoftening(IntGlo,D1,D2, GG,cp,cr,phip,phir, &
       psip,psir,factor,c,phi,psi,Sig0,EpsP,DEps,DEpsP,SigC,IPL)
        !**********************************************************************
        !
        ! Elastoplastic constitutive model with STRAIN SOFTENING, based on the 
        ! MOHR-COULOMB criterion (considering modifications of Abbo & Sloan (1995))
        ! Following Ortiz and Simo (1986) to determine stress update
        !
        !**********************************************************************

        implicit none

        !Local variables
        integer :: i,n,m,it
        double precision :: F,F0,F2 !Evaluation of the Yield function
        double precision :: alpha !Elastic Strain proportion
        double precision :: SSTOL !Tolerance Relative Error
        double precision :: YTOL !Tolerance Relative Error of the yield function evaluation
        double precision :: SPTOL !Tolerance Softening parameters
        double precision :: Rn !Relative error
        double precision :: T,DT,T1,beta,DTmin !Substepping parameters
        double precision :: c1,phi1,psi1,c2,phi2,psi2
        double precision :: ctol,phitol,psitol !c,phi,psi tolerances
        double precision :: Dcr,Dphir,Dpsir !Diference between current and residial values
        double precision :: moduleEr,moduleSigDSig
        double precision :: EpsPEq,EpsPEq1,EpsPEq2 !Equivalent Plastic Deformation
        double precision :: DEpsPEq !Derivative Strain in function of Equivalent Plastic Deformation
        double precision, dimension(6) :: SigYield, SigYield2
        double precision, dimension(6) :: DSigPP,DSigP1,DSigP2
        double precision, dimension(6) :: DEpsPP,DEpsPP1,DEpsPP2
        double precision, dimension(6) :: DEpsS,DEpsSS
        double precision, dimension(6) :: EpsP1,EpsP2
        double precision, dimension(6) :: DEpsPEqDPS,DEpsPEqDPS1
        double precision, dimension(6) :: sumSg,Er
        double precision, dimension(3) :: DSPDPEq,DSPDPEq1 !Variation of softening parameters (c,phi,psi) in function of plastic strain
        !In variables
        integer, intent(in) :: IntGlo !Global ID of Gauss point or particle
        double precision, intent(in) :: D1,D2,GG !Elastic Parameters
        double precision, intent(in) :: cp,cr,phip,phir,psip,psir,factor !Softening parameter
        !Inout variables
        double precision, intent(inout), dimension(6) :: DEps !Incremental total strain
        double precision, intent(inout):: c,phi,psi !cohesion,friction angle and dilatancy angle
        double precision, intent(inout), dimension(6) :: EpsP !Accumulated Plastic Strain
        double precision, intent(inout), dimension(6) :: Sig0 !Initial Stress
        double precision, intent(inout), dimension(6) :: SigC !Final Stress
        double precision, intent(inout), dimension(6) :: DEpsP !Incremental plastic strain

        !Out variables
        integer, intent(out) :: IPL

        !Initialization
        DEpsPEq = 0.0d0
        EpsPEq = 0.0d0
        SigYield = 0.0d0
        DEpsP = 0.0d0
        F = 0.0d0
        it = 0

        if (c > cp.or.phi > phip.or.psi > psip) then
            c = min(c,cp)
            phi = min(phi,phip)
            psi = min(psi,psip)
        end if
        
        if (c < cr.or.phi < phir.or.psi < psir) then
            c = max(c,cr)
            phi = max(phi,phir)
            psi = max(psi,psir)
        end if

        !Tolerances
        SSTOL = 0.01d0 !Tolerance Relative Error (10-3 to 10-5)
        YTOL = 1e-8 !Tolerance Error on the Yield surface (10-6 to 10-9)
        SPTOL = 0.01d0 !Tolerance Softening Parameters (0.0001d0)
        ctol = abs(cp-cr)*SPTOL
        phitol = abs(phip-phir)*SPTOL
        psitol = abs(psip-psir)*SPTOL
        DTmin = 0.000000001d0
        
        !Check the yield function value
        call DetermineYieldFunctionValue(IntGlo,SigC,c,phi,F)
        
        !If F<0 then the behaviour is elastic --> Return
        if (F <= YTOL) then
            IPL = 0
            return
        end if

        !If F>0, the behaviour is elastoplastic --> Continue
        Dcr = abs(c - cr)
        Dphir = abs(phi - phir)
        Dpsir = abs(psi - psir)
        !Check if we are in residual conditions or in softening conditions
        if (Dcr <= ctol.and.Dphir <= phitol.and.Dpsir <= psitol) then
            IPL = 1 !IPL=1 Residual conditions --> no changes of the strength parameters
            c = cr
            phi = phir
            psi = psir
        else
            IPL = 2 !IPL=2 Softening Conditions --> changes of the strength parameters
        end if
    
        
        call MCSS_Ortiz_Simo_Integration(GG, D1, D2, IntGlo, Sig0, c, phi, psi, factor, DEps, EpsP, dEpsP, cr, phir, psir, cp, phip, &
                                              psip, ctol, phitol, psitol, YTOL)
        
        ! State parameters {phi, psi, c} updated inside ortiz-simo
        ! EpsP updated inside of the integration
        ! dEpsP updated inside of ortiz-Simo
        
        ! Update the Stress
        SigC = Sig0
        
        ! Increment of elastic stress not updated
    end subroutine MOHRStrainSoftening
      
    Subroutine MCSS_Ortiz_Simo_Integration(G, D1, D2, IntGlo, Sig, c, phi, psi, factor, dEps, EpsP, dEpsP, cr, phir, psir, cp, phip, &
                                              psip, ctol, phitol, psitol, FTOL)
        !**********************************************************************
        ! Function:  To update the stress (Sig) and plastic strain (EpsP)
        ! Follows Ortiz and Simo (1986) https://doi.org/10.1002/nme.1620230303
        !
        ! Last Modified: 11/10/2023
        ! Author: Jonathan J. Moore
        !**********************************************************************
    ! Subroutine: Determines the change in stress (Sig), increment of plastic strain dEpsP
    
          ! List the input variables
              ! G: Shear modulus
              ! D1, D2: Components of the elastic stiffness  matrix
              ! IntGlo: Global iD of the Gauss point or particle
              ! Sig: Current stress state
              ! c: Cohesion
              ! phi: Friction angle
              ! psi: Dilatancy angle
              ! factor: Softening parameter
              ! dEps: Total strain increment
              ! EpsP:  Accumulated plastic strain
              ! cr: residual cohesion value 
              ! phir: residual friction angle
              ! psir: residual dilatancy angle
              ! cp: peak cohesion value
              ! phip: peak friction angle
              ! psip: peak dilatancy angle
              ! ctol: cohesion tolerance around the residual value (cr)
              ! phitol: friction angle tolerance around the residual value (phir)
              ! psitol: dilatnacy angle tolerance around the residual value (psir)
      
          ! -------------------------- Variable Definitions --------------------------
          ! ------------- Scalar Values -------------
          ! In 
          double precision, intent(in) :: G, D1, D2, factor, cr, phir, psir, cp, phip, psip, ctol, phitol, psitol, FTOL
          integer :: IntGlo
          
          ! In/Out 
          double precision, intent(inout) :: c, phi, psi
          ! Out 
          ! double precision, intent(out) ::
      
          ! ------------- Vector Values -------------
          ! In 
          !double precision, intent(inout), dimension(6) :: 
          ! In/Out 
          double precision, intent(inout), dimension(6) :: Sig, dEps, EpsP, dEpsP
      
          ! ------------- Local Variables -------------
          ! Variable definitions
              ! F: Yield surface value
              ! cu:  Updated cohesion value
              ! phiu: Updated friction angle value
              ! psiu: updated dilatancy value
              ! p: Mean stress 
              ! J: Deviatoric stress 
              ! Lode: Lode angle
              ! S3TA: ??
              ! dummyVal_1, dummyVal_2, dummyVal_3: Free variables
              ! H: Hardening parameter (dF/dLambda)
              ! D1, D2: Dummy values  to store stfiffness matrix components
              ! epsPEq: Equivalent plastic  strain (constant scaled norm of the plastic strain)
              ! dLambda: Increment of the plastic multiplier
              ! Maxiter: The maximum number of time the gradient descent  method should be used
              ! counter: Track current number of iterations
          
              ! dummyVec_6d: length 6 free  vector
              ! dEpsPu: Updated increment of plastic strain
              ! EpsPu: Updated value of the total plastic strain
              ! Sigu: Updated stress value
              ! dSigu: Updated increment of stress
              ! m_vec: Normal to the plastic potential (dP/dSig)
              ! n_vec: Normal to the Yield surface (dF/dSig)
              ! DE_m: Elastic stiffness matrix times the plastic potential normal
              ! DEpsPEqDPS: Derivative of the Equivalent plastic strain  wrt to the plastic strain (dEpsPEq/dEpsP)
              ! DSPDPEq: Derivative of the state  parameters wrt the equivalent strain (dXs/dEpsEq)
      
              ! dummyVec_3d: length 3 free vector
      
              ! dFdSP: Derivative of the yield surface wrt  to the  state parameters (dF/dXs)
      
              ! DE: Stiffness matricx
      
          ! Local scalar values
      
          double precision :: F, cu, Phiu, Psiu , J, Lode, S3TA, dummyVal_1, dummyVal_2, dummyVal_3, H, epsPEq, dLambda
          integer:: MaxIter, counter
      
          ! Local vector values
          double precision, dimension(6):: dummyVec_6d, dEpsPu, EpsPu, Sigu, dSigu, &
                                          m_vec, n_vec, DE_m, DEpsPEqDPS
      
          double precision, dimension(3):: dummyVec_3d, DSPDPEq
      
          double precision, dimension(2):: dFdSP ! Derivative of the yield function with respect to the softening parameters (phi, c)
      
          ! Local matrix values
          double precision, dimension(6,6):: DE
      
          ! -------------------------- Begin Calculations --------------------------
      
      
          ! Store variables for updating
          Sigu = Sig
          EpsPu = EpsP
      
          cu = c ! Updated cohesion
          Phiu = Phi! Updated friction angle
          Psiu = Psi ! Updated dilatancy
      
          ! Form the stiffness matix
          DE = 0.0
          DE(1:3,1:3) = D2
          DE(1,1) = D1
          DE(2,2) = D1
          DE(3,3) = D1
          DE(4,4) = G
          DE(5,5) = G
          DE(6,6) = G 
      
          ! Keep the State varaibles constant
      
          ! Calc the elastic predictor for the stresses 
          ! (Assumes that all of strain increment is elastic therfore there is no change in the equivalent plastic strain)
          call MatVec_MohrCoulombStrainSoftening(DE, 6, dEps, 6, dSigu)
      
          ! Update the stresses
          Sigu = Sigu + dSigu
      
          ! Evalue the yield surface
          call DetermineYieldFunctionValue(IntGlo, Sigu, cu, phiu, F)
      
          if (abs(F) < FTOL) then
              ! Prediction of the stress and strain values are correct and the values can be updated and returned
      
              ! Update Sig, EpsP, dEpsP
              Sig = Sigu
              EpsP(:) = 0
              dEpsP(:) = 0
      
              ! Update yield surface values
              c = cu
              phi = phiu
              psi = psiu
      
              ! Exit the subroutine
              return
          end if
      
          ! Max number of plastic descent iterations
          MaxIter = 100000
          counter = 0
          
          do while(abs(F) >= FTOL .and. counter <= MaxIter)
              call CalculateInvariants(IntGlo, Sigu, p, J, Lode, S3TA)
      
              ! Calc the equivalent plastic strain
              call CalculateEpsPEq(EpsPu, epsPEq)     
      
              ! Calc n_vec, m_vec
              call CalculateDerivativesYieldFunctAndPlasticPotential(Sigu, p, J, Lode, S3TA, cu, phiu, psiu, n_vec, m_vec)
      
              ! dF/dXs
              call CalculateDerivativesYieldFunctSofteningParameters(p, J, Lode, S3TA, cu, phiu, dFdSP)
      
              ! dXs/dEpsPEq
              call CalculateDerivativesStrSoftParamRespectEquivalentPlasticStrain(factor,cp,cr,&
                                      phip,phir,psip,psir,EpsPEq,DSPDPEq)
      
              ! Calc dEpsPEq/dEpsP
              call CalculateDerivativesEquivalentPlasticStrainRespectPlasticStrain(EpsPu, EpsPEq, DEpsPEqDPS)
      
              ! Calc D * m
              call MatVec_MohrCoulombStrainSoftening(DE, 6, m_vec, 6, DE_m)
      
              ! Compute n_vec.DE.m_vec
              call DotProduct_2(n_vec, DE_m, 6, dummyVal_1)
      
              ! Make a 1x3 vector to store dF/dXs
              dummyVec_3d(:) = 0
              dummyVec_3d(1) = dFdSP(1)
              dummyVec_3d(2) = dFdSP(2)
      
              ! Calc the dot product between dF/dXs.dXs/dEpsPEq
              call DotProduct_2(dummyVec_3d, DSPDPEq, 3, dummyVal_2)
      
              ! Calc the dot product between dEpsPEq/dEpsP.dP/dSig
              call DotProduct_2(DEpsPEqDPS, m_vec, 6, dummyVal_3)
      
              ! Need to calc the hardening/softening parameter (H)
              ! H = dF/dXs.dXs/dEpsEq * dEpsEq/dEpsP.dP/dSig
              H = dummyVal_2 * dummyVal_3
      
              ! calc dLambda (Increment of the plastic multiplier)
              dLambda = F/(dummyVal_1 - H)
              !dLambda = F/(dummyVal_1)
      
              ! Compute the stress update
              Sigu = Sigu - dLambda * DE_m
      
              ! Accumulate plastic strain
              EpsPu = EpsPu + dLambda * m_vec
      
              ! Calc the equivalent plastic strain
              call CalculateEpsPEq(EpsPu, epsPEq)     
      
              ! Update the state parameters (c, phi, psi)
              call CalculateSofteningParameters(epsPEq,factor,cp,cr,phip,phir,psip,psir,cu,phiu,psiu)
      
              ! Calc the yield function value
              call DetermineYieldFunctionValue(IntGlo,Sigu,cu,phiu,F)
      
              ! Update the counter
              counter = counter + 1
          end do
      
          ! Retun the integrated parameters
          Sig = Sigu
          dEpsP = EpsPu-EpsP
          EpsP = EpsPu
      
          c = cu
          Phi = Phiu
          Psi = Psiu
      
    end Subroutine MCSS_Ortiz_Simo_Integration

      Subroutine DetermineElasticProportionPegasusMethod(IntGlo,Sig,DSig,DEps,c,phi,YTOL,alpha)
        !**********************************************************************
        !
        ! The PEGASUS METHOD method is used  
        !
        !**********************************************************************

        implicit none

        !Local variables
        integer :: Its     
        double precision :: alpha0,alpha1,F0,F1,F
        double precision, dimension(6) :: Sig0,Sig1,SigNew
        !In variables
        double precision, intent(in), dimension(6) :: Sig, DSig
        double precision, intent(in), dimension(6) :: DEps
        double precision, intent(in) :: c,phi,YTOL
        integer, intent(in) :: IntGlo       !Global ID of Gauss point or particle
        !Out variables
        double precision, intent(out) :: alpha

        alpha0 = 0.0d0
        alpha1 = 1.0d0

        Sig0 = Sig + alpha0*DSig ! = Sig0
        Sig1 = Sig + alpha1*DSig ! = SigE
            
        call DetermineYieldFunctionValue(IntGlo,Sig0,c,phi,F0)
        call DetermineYieldFunctionValue(IntGlo,Sig1,c,phi,F1)

        F=YTOL+1000
        Its = 0 !Counter

        do while (abs(F) > YTOL.and.Its < 1000)
            alpha = alpha1 - F1*(alpha1-alpha0)/(F1-F0)
            SigNew = Sig + alpha*DSig
        
            call DetermineYieldFunctionValue(IntGlo,SigNew,c,phi,F)

            if ((F*F1) < 0.0d0) then
                alpha0 = alpha1
                F0 = F1
            else
                F0 = F1*F0/(F1+F)
            end if

            alpha1 = alpha
            F1 = F
            Its = Its + 1

        end do  
        if (Its >= 1000) then
                alpha = 0.0d0
        end if
      end subroutine DetermineElasticProportionPegasusMethod

      Subroutine CalculateInvariants(IntGlo,Sig,p,J,Lode,S3TA)
        !**********************************************************************
        !
        ! Calcuation of the invariants (defined as Abbo & Sloan (1995))
        !
        !**********************************************************************

        implicit none

        !Local variables
        double precision :: Sx,Sy,Sz,SqTxy,SqTyz,SqTxz,suma,h1,h2,J2,J3
        double precision, parameter :: C00000 = 0.0D0
        double precision, parameter :: C00001 = 1.0D0
        double precision, parameter :: C00P16 = 0.166666666666666D0
        double precision, parameter :: C00002 = 2.0D0
        double precision, parameter :: C00003 = 3.0D0
        double precision, parameter :: CP3333 = 0.333333333333333D0
        double precision, parameter :: C00IR3 = 0.577350269189626D0
        double precision, parameter :: TINY = 0.000000000000001D0
        !In variables
        double precision, intent(in), dimension(6) :: Sig
        integer, intent(in) :: IntGlo !Global ID of Gauss point or particle
        !Out variables
        double precision, intent(out) :: p,J,Lode,S3TA !Invariants

        p = C00000
        J = C00000
        Lode = C00000

        !Mean stress (p)
        p = CP3333 * (Sig(1) + Sig(2) + Sig(3))

        !Deviatoric stress (J)
        Sx = Sig(1) - p
        Sy = Sig(2) - p
        Sz = Sig(3) - p
        suma = (Sig(1)-Sig(2))*(Sig(1)-Sig(2))+(Sig(1)-Sig(3))*(Sig(1)-Sig(3))+(Sig(2)-Sig(3))*(Sig(2)-Sig(3))
        SqTxy =  Sig(4) * Sig(4)
        SqTyz =  Sig(5) * Sig(5)
        SqTxz =  Sig(6) * Sig(6)

        J2 = C00P16 * suma + SqTxy + SqTyz + SqTxz
        J3 = Sx*Sy*Sz + C00002 * Sig(4)*Sig(5)*Sig(6) - Sx*SqTyz - Sy*SqTxz - Sz*SqTxy
        J = SQRT(J2)

        !Lode's angle (Lode)
        if (J2 > C00000) then

            h1 = -C00003/(C00002*C00IR3)
            h2 = J3/(J*J*J)
            S3TA = h1*h2
            if (S3TA < -C00001) then
                S3TA = -C00001
            else if (S3TA > C00001) then
                S3TA = C00001
        end if
            Lode = CP3333*asin(S3TA)
        else  !Special case of zero deviatoric stress
            Lode = C00000
            S3TA = C00000
        end if

      end subroutine CalculateInvariants

      Subroutine DetermineYieldFunctionValue(IntGlo,Sig,c,phi,F)
        !**********************************************************************
        !
        ! In this subroutine the yield function evaluated is a smooth hyperbolic approximation to the
        ! Mohr-Coulomb yield criterion (Abbo and Sloan, 1995).
        !
        ! The edges of the hexagonal pyramid and the tip have been smoothed.
        ! There are two parameters aSmooth (smoothes the tip) and ATTRAN(smoothes the edges)
        ! In this case aSmooth=0.0005*c*cot(phi) and LodeT=25�.
        ! If aSmooth=0 and LodeT=30� the original Mohr-Coulomb is obtained.
        !
        !**********************************************************************

        implicit none

        !Local variables
        double precision ::  p,J,Lode,S3TA !Invariants
        double precision ::  COH, SPHI, CPHI, COTPHI, STA, CTA, K, aSmooth, ASPHI2, SGN, A, B
        double precision, parameter :: C00001 = 1.0d0 !Parameters
        double precision, parameter :: C00003 = 3.0d0
        double precision, parameter :: C00P50 = 0.0005d0
        double precision, parameter :: C00000 = 0.0d0
        double precision, parameter :: C00IR3 = 0.577350269189626d0
        double precision, parameter :: C000P1 = 0.00000000001D0
        !Constants for rounded K function (for LodeT=25)
        !double precision, parameter :: A1 = 1.432052062044227d0
        !double precision, parameter :: A2 = 0.406941858374615d0
        !double precision, parameter :: B1 = 0.544290524902313d0
        !double precision, parameter :: B2 = 0.673903324498392d0
        !double precision, parameter :: ATTRAN = 0.436332312998582d0 !Smoothing parameter: LodeT in radians
        !Constants for rounded K function (for LodeT=29.5)
        double precision, parameter :: A1 = 7.138654723242414d0
        double precision, parameter :: A2 = 6.112267270920612d0
        double precision, parameter :: B1 = 6.270447753139589d0
        double precision, parameter :: B2 = 6.398760841429403d0
        double precision, parameter :: ATTRAN = 0.514872129338327d0 !Smoothing parameter: LodeT in radians
        !Constants for rounded K function (for LodeT=30)
        !double precision, parameter :: A1 = -138300705.446275
        !double precision, parameter :: A2 = -138300706.472675
        !double precision, parameter :: B1 = -138300706.3123
        !double precision, parameter :: B2 = 0.192450089729875
        !double precision, parameter :: ATTRAN = 0.523598776 !Smoothing parameter: LodeT in radians
        !In variables
        double precision, intent(in), dimension(6) :: Sig
        double precision, intent(in) :: c,phi
        integer, intent(in) :: IntGlo !Global ID of Gauss point or particle

        !Out variables
        double precision, intent(out) :: F

        F = C00000

        !Calculation of the invariants (p',J,Lode)
        call CalculateInvariants(IntGlo,Sig,p,J,Lode,S3TA)

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !!Evaluation of the yield function with Smoothing!!!!!!!!!!!!!!!!!!!!!!!!
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !Material parameters
        COH = c     !Cohesion
        SPHI = sin(phi) 
        CPHI = cos(phi)
        COTPHI = CPHI/SPHI
        aSmooth = C00P50*COH*COTPHI !Smoothing parameter
        ASPHI2 = aSmooth*aSmooth*SPHI*SPHI
        if (abs(phi) == C00000) then
                ASPHI2 = C00P50*C00P50*COH*COH*CPHI*CPHI
        end if

        !Calculate K function
        if (abs(Lode) < ATTRAN) then
            STA = sin(Lode)
            CTA = cos(Lode)
            K = CTA - STA*SPHI*C00IR3
        else
            SGN = SIGN(C00001,Lode)
            A = A1 + A2*SGN*SPHI
            B = B1*SGN + B2*SPHI
            K = A - B*S3TA
        end if

        !Calculate value of Hyperbolic Yield function
        F = p*SPHI + sqrt(J*J*K*K+ASPHI2) - COH*CPHI
      
      end subroutine DetermineYieldFunctionValue


      Subroutine CalculateDerivativesYieldFunctAndPlasticPotential(Sig,p,J,Lode,S3TA,c,phi,psi,DFDSig,DPPDSig)
        !**********************************************************************
        !
        ! Calculation of the derivatives of the yield function (F) and the plastic potencial punction (P).
        ! Based on Abbo & Sloan (1995)
        !
        !**********************************************************************

        implicit none

        !Local variables
        integer :: i
        double precision :: COH, SPHI, CPHI, TPHI, COTPHI, STA, CTA, A, B,&
                            D, aSmooth, ASPHI2, SGN, T3TA, C3TA, J2, psi2
        double precision ::   K, dKdLode
        double precision :: SPSI, CPSI, TPSI, COTPSI, ASPSI2
        double precision :: i1, i2, Sx, Sy, Sz
        double precision :: DFDp,DFDJ,DFDLode !Derivatives F respect Invariants
        double precision :: DPDp,DPDJ,DPDLode !Derivatives P respect Invariants
        double precision :: C1, C2, C3
        double precision, dimension(6):: DpDSig,DJDSig,DJ3DSig !Derivatives Invariants

        double precision, parameter :: C00001 = 1.0D0 !Parameters
        double precision, parameter :: C000P5 = 0.5D0
        double precision, parameter :: C00P50 = 0.0005D0
        double precision, parameter :: C00000 = 0.0D0
        double precision, parameter :: C00003 = 3.0D0
        double precision, parameter :: C00004 = 4.0D0
        double precision, parameter :: C00002 = 2.0D0
        double precision, parameter :: CP3333 = 0.333333333333333D0
        double precision, parameter :: C00IR3 = 0.577350269189626D0
        double precision, parameter :: C0R3I2 = 0.866025403784439D0
        double precision, parameter :: C000P1 = 0.000000000000001D0 
        double precision, parameter :: J0 = 0.001D0 
        !Constants for rounded K function (for LodeT=25)
        !double precision, parameter :: A1 = 1.432052062044227d0
        !double precision, parameter :: A2 = 0.406941858374615d0
        !double precision, parameter :: B1 = 0.544290524902313d0
        !double precision, parameter :: B2 = 0.673903324498392d0
        !double precision, parameter :: ATTRAN = 0.436332312998582d0 !Smoothing parameter: LodeT in radians
        !Constants for rounded K function (for LodeT=29.5)
        double precision, parameter :: A1 = 7.138654723242414d0
        double precision, parameter :: A2 = 6.112267270920612d0
        double precision, parameter :: B1 = 6.270447753139589d0
        double precision, parameter :: B2 = 6.398760841429403d0
        double precision, parameter :: ATTRAN = 0.514872129338327d0 !Smoothing parameter: LodeT in radians
        !Constants for rounded K function (for LodeT=30)
        !double precision, parameter :: A1 = -138300705.446275
        !double precision, parameter :: A2 = -138300706.472675
        !double precision, parameter :: B1 = -138300706.3123
        !double precision, parameter :: B2 = 0.192450089729875
        !double precision, parameter :: ATTRAN = 0.523598776 !Smoothing parameter: LodeT in radians
        !In variables
        double precision, intent(in) ::  c,phi,psi !Soft Parameters
        double precision, intent(in), dimension(6) :: Sig
        !Out variables
        double precision, intent(out), dimension(6) :: DFDSig, DPPDSig !Derivatives respect Sigma
        !Inout variables
        double precision, intent(inout) :: p,J,Lode,S3TA !Invariants

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !!!!! DFDSig = C1*DPDSig + C2*DJDSig + C3*DJ3DSig  !!!!!!!!!!!!!!!!!!!!!
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        !Material parameters
        COH = c !Cohesion
        SPHI = sin(phi) 
        CPHI = cos(phi)
        COTPHI = CPHI/SPHI
        aSmooth = C00P50*COH*COTPHI !Smoothing parameter
        ASPHI2 = aSmooth*aSmooth*SPHI*SPHI
        if (abs(phi) == C00000) then
            ASPHI2 = C00P50*C00P50*COH*COH*CPHI*CPHI
        end if

        if (J == C00000) then
            J2 = C000P1
            J = sqrt(J2)
        else
            J2 = J*J
        end if

        CTA = cos(Lode)
        C3TA = CTA*(C00004*CTA*CTA-C00003)
        T3TA = S3TA/C3TA

        !Calculate K function and its derivative
        if (abs(Lode) < ATTRAN) then
            STA = S3TA/(C00004*CTA*CTA-C00001)
            K = CTA - STA*SPHI*C00IR3
            dKdLode =  - STA - C00IR3*SPHI*CTA
        else
            SGN = SIGN(C00001,Lode) ! It puts the Lode's sign to the number 1
            A = A1 + A2*SGN*SPHI
            B = B1*SGN + B2*SPHI
            K = A - B*S3TA
            dKdLode = - C00003*B*C3TA
        end if
        
        !Derivative Dp/DSig
        DpDSig(1) = CP3333
        DpDSig(2) = CP3333
        DpDSig(3) = CP3333
        DpDSig(4) = C00000
        DpDSig(5) = C00000
        DpDSig(6) = C00000
        
        !Derivative DJ/DSig
        i1 = C000P5/J
        if (J < 0.0001) then
            i1 = 0.0d0
        end if
        Sx = Sig(1)-p
        Sy = Sig(2)-p
        Sz = Sig(3)-p

        DJDSig(1) = i1 * Sx
        DJDSig(2) = i1 * Sy
        DJDSig(3) = i1 * Sz
        DJDSig(4) = i1 * C00002 * Sig(4)
        DJDSig(5) = i1 * C00002 * Sig(5)
        DJDSig(6) = i1 * C00002 * Sig(6)

        !Derivative DJ3/DSig
        i2 = CP3333*J*J
        DJ3DSig(1) = (Sy*Sz - Sig(5)*Sig(5) + i2)
        DJ3DSig(2) = (Sx*Sz - Sig(6)*Sig(6) + i2)
        DJ3DSig(3) = (Sx*Sy - Sig(4)*Sig(4) + i2)
        DJ3DSig(4) = C00002*(Sig(5)*Sig(6) - Sz*Sig(4))
        DJ3DSig(5) = C00002*(Sig(6)*Sig(4) - Sx*Sig(5))
        DJ3DSig(6) = C00002*(Sig(4)*Sig(5) - Sy*Sig(6))

        D = J*K/(sqrt(J2*K*K + ASPHI2))

        !C1F
        C1 = SPHI
        !C2F
        C2 = D*K - T3TA*D*dKdLode
        !C3F
        C3 = -C0R3I2*dKdLode*D/(J2*C3TA)
        
        !DFDSig!
        do i=1,6
            DFDSig(i) = C1*DpDSig(i) + C2*DJDSig(i) + C3*DJ3DSig(i)
        end do

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !!!!! DPPDSig = DFDSig (if associated Flow Rule)  !!!!!!!!!!!!!!!!!!!!!!
        !!!!! or
        !!!!! DPPDSig = C1*DPDSig + C2*DJDSig + C3*DJ3DSig  !!!!!!!!!!!!!!!!!!!!
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        if (abs(J) < J0) then
            psi2 = phi - abs(J)*(phi - psi)/J0
        else
            psi2 = psi
        end if

        if (phi == psi2) then !If Associated Flow Rule, then DPPDSig = DFDSig
            DPPDSig = DFDSig

        else !If Non-Associated Flow Rule, then calculate...
            !Material parameters
            SPSI = sin(psi2) 
            CPSI = cos(psi2)
            if (SPSI<0.0001) then
                COTPSI=0
            else
                COTPSI = CPSI/SPSI
            end if
            aSmooth = C00P50*COH*COTPSI !Smoothing parameter
            ASPSI2 = aSmooth*aSmooth*SPSI*SPSI
            if (abs(psi2) == C00000)then
                ASPSI2 = C00000
            end if

            !Calculate K function and its derivative
            if (abs(Lode) <= ATTRAN) then
                K = CTA - STA*SPSI*C00IR3
                dKdLode = - STA - C00IR3*SPSI*CTA
            else
                A = A1 + A2*SGN*SPSI
                B = B1*SGN + B2*SPSI
                K = A - B*S3TA
                dKdLode = - C00003*B*C3TA
            end if

            D = J*K/(sqrt(J*J*K*K + ASPSI2))

            !C1F
            C1 = SPSI
            !C2F
            C2 = D*K - T3TA*D*dKdLode
            !C3F
            C3 = -C0R3I2*dKdLode*D/(J2*C3TA)

            !DPPDSig
            do i=1,6
                DPPDSig(i) = C1*DpDSig(i) + C2*DJDSig(i) + C3*DJ3DSig(i)
            end do

        end if

      end subroutine CalculateDerivativesYieldFunctAndPlasticPotential


      Subroutine CalculateDerivativesYieldFunctSofteningParameters(p,J,Lode,S3TA,c,phi,DFDSP)
        !**********************************************************************
        !
        ! Calculation of the derivatives of the yield function (F) with respect the strength parameters
        ! The strength parameters are: cohesion (COH) and friction angle (PHI)
        !
        !**********************************************************************

        implicit none

        !Local variables
        double precision :: COH, SPHI, CPHI, TPHI, COTPHI, STA, CTA, A, B,&
                            Denom, Num, aSmooth, ASPHI2, SGN
        double precision :: K, dKdPhi, dadc, dadPhi
        double precision, parameter :: C00001 = 1.0D0 !Parameters
        double precision, parameter :: C00P50 = 0.0005D0
        double precision, parameter :: C00000 = 0.0D0
        double precision, parameter :: C00003 = 3.0D0
        double precision, parameter :: C00002 = 2.0D0
        double precision, parameter :: C00IR3 = 0.577350269189626D0
        double precision, parameter :: C000P1 = 0.00000000001D0
        !Constants for rounded K function (for LodeT=25)
        !double precision, parameter :: A1 = 1.432052062044227d0
        !double precision, parameter :: A2 = 0.406941858374615d0
        !double precision, parameter :: B1 = 0.544290524902313d0
        !double precision, parameter :: B2 = 0.673903324498392d0
        !double precision, parameter :: ATTRAN = 0.436332312998582d0 !Smoothing parameter: LodeT in radians
        !Constants for rounded K function (for LodeT=29.5)
        double precision, parameter :: A1 = 7.138654723242414d0
        double precision, parameter :: A2 = 6.112267270920612d0
        double precision, parameter :: B1 = 6.270447753139589d0
        double precision, parameter :: B2 = 6.398760841429403d0
        double precision, parameter :: ATTRAN = 0.514872129338327d0 !Smoothing parameter: LodeT in radians
        !Constants for rounded K function (for LodeT=30)
        !double precision, parameter :: A1 = -138300705.446275
        !double precision, parameter :: A2 = -138300706.472675
        !double precision, parameter :: B1 = -138300706.3123
        !double precision, parameter :: B2 = 0.192450089729875
        !double precision, parameter :: ATTRAN = 0.523598776 !Smoothing parameter: LodeT in radians

        !In variables
        double precision, intent(in) :: p,J,Lode,S3TA !Invariants
        double precision, intent(in) :: c,phi !Soft Parameters
        !Out variables
        double precision, intent(out), dimension(2) :: DFDSP !Derivatives respect Soft Parameters


        !Material parameters
        COH = c !Cohesion
        SPHI = sin(phi) 
        CPHI = cos(phi)
        COTPHI = CPHI/SPHI

        !Calculate aSmooth and its derivatives
        if (abs(phi) == C00000) then
            COTPHI = C00000
            dadc = C00000
            dadPhi = C00000
        else
            dadc = C00P50*CPHI/SPHI
            dadPhi = - C00P50*COH/(SPHI*SPHI)
        end if
        aSmooth = C00P50*COH*COTPHI !Smoothing parameter
        ASPHI2 = aSmooth*aSmooth*SPHI*SPHI
        if (abs(phi) == C00000) then
        ASPHI2 = C00P50*C00P50*COH*COH*CPHI*CPHI
        end if

        !Calculate K function and its derivatives
        if (abs(Lode) <= ATTRAN) then
            STA = sin(Lode)
            CTA = cos(Lode)
            K = CTA - STA*SPHI*C00IR3
            dKdPhi = - C00IR3*CPHI*STA
        else
            SGN = SIGN(C00001,Lode) !It puts the Lode's sign to the number 1
            A = A1 + A2*SGN*SPHI
            B = B1*SGN + B2*SPHI
            K = A - B*S3TA
            dKdPhi = A2*SGN*CPHI - B2*CPHI*S3TA
        end if

        !Operating..
        Denom = (sqrt(J*J*K*K + ASPHI2))
        Num =  J*J*K*dKdPhi + aSmooth*SPHI*SPHI*dadPhi + aSmooth*aSmooth*SPHI*CPHI

        !Derivative DF/Dc
        DFDSP(1) = aSmooth*SPHI*SPHI*dadc/Denom - CPHI

        !Derivative DF/Dphi
        DFDSP(2) = p*CPHI + Num/Denom + COH*SPHI

        if (J <= C00000) then
            DFDSP(1) = - CPHI
            DFDSP(2) = p*CPHI + COH*SPHI
        end if

        end subroutine CalculateDerivativesYieldFunctSofteningParameters


        subroutine CalculateDerivativesStrSoftParamRespectEquivalentPlasticStrain(factor,cp,cr,&
                                phip,phir,psip,psir,EpsPEq,DSPDPEq)
        !**********************************************************************
        !
        ! Calculation of the derivatives of the strength parameters with respect
        ! the equivalent plastic shear strain
        !
        !**********************************************************************

        implicit none

        !In Variables
        double precision, intent(in) :: EpsPEq
        double precision, intent(in) :: factor,cp,cr,phip,phir,psip,psir
        !Out Variables
        double precision, intent(out), dimension(3):: DSPDPEq
        
        !Derivative Cohesion respect Equivalent Plastic Strain (Dc/DPEq)
        DSPDPEq(1) = -factor * (cp - cr) * (exp(-factor*EpsPEq))
        !Derivative Friction angle respect Equivalent Plastic Strain (Dphi/DPEq)
        DSPDPEq(2) = -factor * (phip - phir) * (exp(-factor*EpsPEq))
        !Derivative Dilatancy angle respect Equivalent Plastic Strain (Dpsi/DPEq)
        DSPDPEq(3) = -factor * (psip - psir) * (exp(-factor*EpsPEq))

      end subroutine CalculateDerivativesStrSoftParamRespectEquivalentPlasticStrain

      Subroutine CalculateDerivativesEquivalentPlasticStrainRespectPlasticStrain(EpsP,EpsPEq,DEpsPEqDPS)
        !**********************************************************************
        !
        ! Calculation of the derivatives of the equivalent plastic shear strain
        ! with respect the plastic strain
        !
        !**********************************************************************

        implicit none

        !Local Variables
        double precision :: k1, k2, k3
        double precision :: EpsPM
        double precision, dimension(3) :: EpsDev
        !In Variables
        double precision, intent(in), dimension(6) :: EpsP
        double precision, intent(in) :: EpsPEq
        !Out Variables
        double precision, intent(out), dimension(6):: DEpsPEqDPS
        
        if (EpsPEq < 0.00000000001d0) then
            k1 = 0.0d0
        else
            k1 = 2.0d0/(3.0d0*EpsPEq)
        end if
        
        k2 = k1 * 1.0d0/3.0d0
        k3 = k1 * 2.0d0

        EpsPM = k2 * (EpsP(1) + EpsP(2) + EpsP(3))
        EpsDev(1) = EpsP(1)-EpsPM
        EpsDev(2) = EpsP(2)-EpsPM
        EpsDev(3) = EpsP(3)-EpsPM

        DEpsPEqDPS(1) = k2 * ( 2.0d0*EpsDev(1) - EpsDev(2) - EpsDev(3))
        DEpsPEqDPS(2) = k2 * (-EpsDev(1) + 2.0d0*EpsDev(2) - EpsDev(3))
        DEpsPEqDPS(3) = k2 * (-EpsDev(1) - EpsDev(2) + 2.0d0*EpsDev(3))
        DEpsPEqDPS(4) = k3 * EpsP(4)
        DEpsPEqDPS(5) = k3 * EpsP(5)
        DEpsPEqDPS(6) = k3 * EpsP(6)

        end subroutine CalculateDerivativesEquivalentPlasticStrainRespectPlasticStrain


        subroutine CalculateEpsPEq(EpsP,EpsPEq)
        !**********************************************************************
        !
        ! Calculation of the equivalent plastic shear strain 
        !
        !**********************************************************************
        
        implicit none

        !Local variables
        double precision:: EpsPM, C1, C2
        double precision, dimension(3) :: EpsDev
        !In variables
        double precision, intent(in), dimension(6) :: EpsP
        !Out variables
        double precision, intent(out) :: EpsPEq
        
        !EpsPEq = ((2/3)ep:ep)^(1/2), ep is the deviatoric plastic strain
        
        EpsPM = (1.0d0/3.0d0) * (EpsP(1) + EpsP(2) + EpsP(3))
        EpsDev(1) = EpsP(1)-EpsPM
        EpsDev(2) = EpsP(2)-EpsPM
        EpsDev(3) = EpsP(3)-EpsPM
        C1 = 2.0d0/3.0d0
        C2 = C1 * 2.0d0
        
        EpsPEq = sqrt(C1*EpsDev(1)*EpsDev(1) + C1*EpsDev(2)*EpsDev(2) +  C1*EpsDev(3)*EpsDev(3) +&
                        C2*EpsP(4)*EpsP(4) + C2*EpsP(5)*EpsP(5) + C2*EpsP(6)*EpsP(6))

      end subroutine CalculateEpsPEq
      

      !Subroutine CalculateIncrementSofteningParameters(DSPDPEq,DEpsPEqDPS,DEpsP,Dh)
      !!**********************************************************************
      !!
      !! Calculation of the increment of the strenght parameters due to the softening
      !!
      !!**********************************************************************
      !
      !implicit none
      !
      !!Local variables
      !double precision :: k
      !!In variables
      !double precision, intent(in), dimension(3) :: DSPDPEq
      !double precision, intent(in), dimension(6) :: DEpsPEqDPS
      !double precision, intent(in), dimension(6) :: DEpsP
      !!Out variables
      !double precision, intent(out), dimension(3) :: Dh
      !
      !
      !k = DEpsPEqDPS(1)*DEpsP(1) + DEpsPEqDPS(2)*DEpsP(2) + DEpsPEqDPS(3)*DEpsP(3) + 
      !*       DEpsPEqDPS(4)*DEpsP(4) + DEpsPEqDPS(5)*DEpsP(5) + DEpsPEqDPS(6)*DEpsP(6)
      
      
      !Dh(1) = DSPDPEq(1)*k
      !Dh(2) = DSPDPEq(2)*k
      !Dh(3) = DSPDPEq(3)*k
      
      !Dh(1) = min (Dh(1) , 0.0d0)
      !Dh(2) = min (Dh(2) , 0.0d0)
      !Dh(3) = min (Dh(3) , 0.0d0)
      
      !end subroutine CalculateIncrementSofteningParameters


      Subroutine CalculateSofteningParameters(EpsPEq,factor,cp,cr,phip,phir,psip,psir,c,phi,psi)
        !**********************************************************************
        !
        ! Calculation of strenght parameters (c, phi, psi)
        !
        !**********************************************************************

        implicit none

        !In variables
        double precision, intent(in) :: EpsPEq,factor,cp,cr,phip,phir,psip,psir
        !Out variables
        double precision, intent(out) :: c,phi,psi  

        c = cr + (cp-cr)*exp(-factor*EpsPEq) 
        phi = phir + (phip-phir)*exp(-factor*EpsPEq) 
        psi = psir + (psip-psir)*exp(-factor*EpsPEq) 

      end subroutine CalculateSofteningParameters

      Subroutine DetermineDSigAndDEpsP(IntGlo,D1,D2,GG,c,phi,psi,Sig,DEpsPEqDPS,DSPDPEq,DEps,DSig,DEpsP)
        !**********************************************************************
        !
        ! Calculation of the stress increment and plastic strain increment
        !
        !         dSig = Dep * dEps
        !         dEpsP = Lambda * DPDSig
        !
        !**********************************************************************

        implicit none

        !Local variables
        integer :: i,k
        double precision :: A,Ai,Denom,Fact,LambdaNum,Lambda
        double precision :: p,J,Lode,S3TA !Invariants
        double precision, dimension(6,6) :: Num,Num1,Prod
        double precision, dimension(6) :: Denom1
        double precision, dimension(6) :: DPPDSig !Derivatives Plastic potential respect net stress
        double precision, dimension(6) :: DFDSig !Derivatives Yield function respect net stress
        double precision, dimension(2) :: DFDSP !Derivatives Yield function respect Soft Parameters
        double precision, dimension(6,6) :: Dep !Elastoplastic Constitutive Matrix 
        !In Variables
        double precision, intent(in) :: c,phi,psi !Softening parameters
        double precision, intent(in) :: D1,D2,GG !Elastic parameters
        double precision, intent(in), dimension(6):: DEpsPEqDPS
        double precision, intent(in), dimension(6) :: Sig
        double precision, intent(in), dimension(3) :: DSPDPEq !Derivatives respect Equivalent Plastic Strain
        double precision, intent(in), dimension(6) :: DEps
        integer, intent(in) :: IntGlo !Global ID of Gauss point or particle
        !Out Variables
        double precision, intent(out), dimension(6) :: DSig
        double precision, intent(out), dimension(6) :: DEpsP

        call CalculateInvariants(IntGlo,Sig,p,J,Lode,S3TA)
        call CalculateDerivativesYieldFunctAndPlasticPotential(Sig,p,J,Lode,S3TA,c,phi,psi,DFDSig,DPPDSig)
        call CalculateDerivativesYieldFunctSofteningParameters(p,J,Lode,S3TA,c,phi,DFDSP)

        !Parameter A (H = -A --> A>0 softening / A<0 hardening)
        A = 0.0d0
        Ai = (DFDSP(1)*DSPDPEq(1) + DFDSP(2)*DSPDPEq(2))
        do i=1,6
        A = A + Ai * DEpsPEqDPS(i) * DPPDSig(i)
        end do

        !Elastoplastic Constitutive Matrix (Dep)
        do i=1,6
            do k=1,6
                Prod(i,k) =  DPPDSig(i) * DFDSig(k)
            end do
        end do

        Num1(1,1) = D1*Prod(1,1) + D2*Prod(2,1) + D2*Prod(3,1)
        Num1(1,2) = D1*Prod(1,2) + D2*Prod(2,2) + D2*Prod(3,2)
        Num1(1,3) = D1*Prod(1,3) + D2*Prod(2,3) + D2*Prod(3,3)
        Num1(1,4) = D1*Prod(1,4) + D2*Prod(2,4) + D2*Prod(3,4)
        Num1(1,5) = D1*Prod(1,5) + D2*Prod(2,5) + D2*Prod(3,5)
        Num1(1,6) = D1*Prod(1,6) + D2*Prod(2,6) + D2*Prod(3,6)

        Num1(2,1) = D2*Prod(1,1) + D1*Prod(2,1) + D2*Prod(3,1)
        Num1(2,2) = D2*Prod(1,2) + D1*Prod(2,2) + D2*Prod(3,2)
        Num1(2,3) = D2*Prod(1,3) + D1*Prod(2,3) + D2*Prod(3,3)
        Num1(2,4) = D2*Prod(1,4) + D1*Prod(2,4) + D2*Prod(3,4)
        Num1(2,5) = D2*Prod(1,5) + D1*Prod(2,5) + D2*Prod(3,5)
        Num1(2,6) = D2*Prod(1,6) + D1*Prod(2,6) + D2*Prod(3,6)

        Num1(3,1) = D2*Prod(1,1) + D2*Prod(2,1) + D1*Prod(3,1)
        Num1(3,2) = D2*Prod(1,2) + D2*Prod(2,2) + D1*Prod(3,2)
        Num1(3,3) = D2*Prod(1,3) + D2*Prod(2,3) + D1*Prod(3,3)
        Num1(3,4) = D2*Prod(1,4) + D2*Prod(2,4) + D1*Prod(3,4)
        Num1(3,5) = D2*Prod(1,5) + D2*Prod(2,5) + D1*Prod(3,5)
        Num1(3,6) = D2*Prod(1,6) + D2*Prod(2,6) + D1*Prod(3,6)

        Num1(4,1) = GG*Prod(4,1)
        Num1(4,2) = GG*Prod(4,2)
        Num1(4,3) = GG*Prod(4,3)
        Num1(4,4) = GG*Prod(4,4)
        Num1(4,5) = GG*Prod(4,5)
        Num1(4,6) = GG*Prod(4,6)

        Num1(5,1) = GG*Prod(5,1)
        Num1(5,2) = GG*Prod(5,2)
        Num1(5,3) = GG*Prod(5,3)
        Num1(5,4) = GG*Prod(5,4)
        Num1(5,5) = GG*Prod(5,5)
        Num1(5,6) = GG*Prod(5,6)

        Num1(6,1) = GG*Prod(6,1)
        Num1(6,2) = GG*Prod(6,2)
        Num1(6,3) = GG*Prod(6,3)
        Num1(6,4) = GG*Prod(6,4)
        Num1(6,5) = GG*Prod(6,5)
        Num1(6,6) = GG*Prod(6,6)



        Num(1,1) = D1*Num1(1,1) + D2*Num1(1,2) + D2*Num1(1,3)
        Num(1,2) = D2*Num1(1,1) + D1*Num1(1,2) + D2*Num1(1,3)
        Num(1,3) = D2*Num1(1,1) + D2*Num1(1,2) + D1*Num1(1,3)
        Num(1,4) = GG*Num1(1,4)
        Num(1,5) = GG*Num1(1,5)
        Num(1,6) = GG*Num1(1,6)

        Num(2,1) = D1*Num1(2,1) + D2*Num1(2,2) + D2*Num1(2,3)
        Num(2,2) = D2*Num1(2,1) + D1*Num1(2,2) + D2*Num1(2,3)
        Num(2,3) = D2*Num1(2,1) + D2*Num1(2,2) + D1*Num1(2,3)
        Num(2,4) = GG*Num1(2,4)
        Num(2,5) = GG*Num1(2,5)
        Num(2,6) = GG*Num1(2,6)

        Num(3,1) = D1*Num1(3,1) + D2*Num1(3,2) + D2*Num1(3,3)
        Num(3,2) = D2*Num1(3,1) + D1*Num1(3,2) + D2*Num1(3,3)
        Num(3,3) = D2*Num1(3,1) + D2*Num1(3,2) + D1*Num1(3,3)
        Num(3,4) = GG*Num1(3,4)
        Num(3,5) = GG*Num1(3,5)
        Num(3,6) = GG*Num1(3,6)

        Num(4,1) = D1*Num1(4,1) + D2*Num1(4,2) + D2*Num1(4,3)
        Num(4,2) = D2*Num1(4,1) + D1*Num1(4,2) + D2*Num1(4,3)
        Num(4,3) = D2*Num1(4,1) + D2*Num1(4,2) + D1*Num1(4,3)
        Num(4,4) = GG*Num1(4,4)
        Num(4,5) = GG*Num1(4,5)
        Num(4,6) = GG*Num1(4,6)

        Num(5,1) = D1*Num1(5,1) + D2*Num1(5,2) + D2*Num1(5,3)
        Num(5,2) = D2*Num1(5,1) + D1*Num1(5,2) + D2*Num1(5,3)
        Num(5,3) = D2*Num1(5,1) + D2*Num1(5,2) + D1*Num1(5,3)
        Num(5,4) = GG*Num1(5,4)
        Num(5,5) = GG*Num1(5,5)
        Num(5,6) = GG*Num1(5,6)

        Num(6,1) = D1*Num1(6,1) + D2*Num1(6,2) + D2*Num1(6,3)
        Num(6,2) = D2*Num1(6,1) + D1*Num1(6,2) + D2*Num1(6,3)
        Num(6,3) = D2*Num1(6,1) + D2*Num1(6,2) + D1*Num1(6,3)
        Num(6,4) = GG*Num1(6,4)
        Num(6,5) = GG*Num1(6,5)
        Num(6,6) = GG*Num1(6,6)



        Denom1(1) = DFDSig(1)*D1 + DFDSig(2)*D2 + DFDSig(3)*D2
        Denom1(2) = DFDSig(1)*D2 + DFDSig(2)*D1 + DFDSig(3)*D2
        Denom1(3) = DFDSig(1)*D2 + DFDSig(2)*D2 + DFDSig(3)*D1
        Denom1(4) = DFDSig(4)*GG
        Denom1(5) = DFDSig(5)*GG
        Denom1(6) = DFDSig(6)*GG

        Denom =   Denom1(1)*DPPDSig(1) + Denom1(2)*DPPDSig(2) + &
                    Denom1(3)*DPPDSig(3) + Denom1(4)*DPPDSig(4) + &
                Denom1(5)*DPPDSig(5) + Denom1(6)*DPPDSig(6) - A

        Fact = 1d0/Denom

        !Dep
        Dep(1,1) = D1 - Fact*Num(1,1)
        Dep(1,2) = D2 - Fact*Num(1,2)
        Dep(1,3) = D2 - Fact*Num(1,3)
        Dep(1,4) = -Fact*Num(1,4)
        Dep(1,5) = -Fact*Num(1,5)
        Dep(1,6) = -Fact*Num(1,6)

        Dep(2,1) = D2 - Fact*Num(2,1)
        Dep(2,2) = D1 - Fact*Num(2,2)
        Dep(2,3) = D2 - Fact*Num(2,3)
        Dep(2,4) = -Fact*Num(2,4)
        Dep(2,5) = -Fact*Num(2,5)
        Dep(2,6) = -Fact*Num(2,6)

        Dep(3,1) = D2 - Fact*Num(3,1)
        Dep(3,2) = D2 - Fact*Num(3,2)
        Dep(3,3) = D1 - Fact*Num(3,3)
        Dep(3,4) = -Fact*Num(3,4)
        Dep(3,5) = -Fact*Num(3,5)
        Dep(3,6) = -Fact*Num(3,6)

        Dep(4,1) = -Fact*Num(4,1)
        Dep(4,2) = -Fact*Num(4,2)
        Dep(4,3) = -Fact*Num(4,3)
        Dep(4,4) = GG - Fact*Num(4,4)
        Dep(4,5) = -Fact*Num(4,5)
        Dep(4,6) = -Fact*Num(4,6)

        Dep(5,1) = -Fact*Num(5,1)
        Dep(5,2) = -Fact*Num(5,2)
        Dep(5,3) = -Fact*Num(5,3)
        Dep(5,4) = -Fact*Num(5,4)
        Dep(5,5) = GG - Fact*Num(5,5)
        Dep(5,6) = -Fact*Num(5,6)

        Dep(6,1) = -Fact*Num(6,1)
        Dep(6,2) = -Fact*Num(6,2)
        Dep(6,3) = -Fact*Num(6,3)
        Dep(6,4) = -Fact*Num(6,4)
        Dep(6,5) = -Fact*Num(6,5)
        Dep(6,6) = GG - Fact*Num(6,6)

        !!!!!!!!! Calculate Plastic multipliler(Lambda)!!!!!!!!!!!!!!!!!
        LambdaNum =   Denom1(1)*DEps(1) + Denom1(2)*DEps(2) + &
                    Denom1(3)*DEps(3) + Denom1(4)*DEps(4) + &
                    Denom1(5)*DEps(5) + Denom1(6)*DEps(6) 
        Lambda =  LambdaNum/Denom

        !!!!!!!!! Determine DSig --> (DSig = Dep*dEps) !!!!!!!!!!!
        do i=1,6
            DSig(i) = 0.0d0
            do k=1,6
                DSig(i) =  DSig(i) + Dep(i,k) * DEps(k)
            end do
        end do

        !!!!!!!!! Determine DEpsP --> (DEpsP = Lambda*DPDSig) !!!!!!!!!!!!
        do i=1,6
            DEpsP(i) = Lambda * DPPDSig(i)
        end do

      end subroutine DetermineDSigAndDEpsP

      Subroutine EndOfStepCorrection(IntGlo,D1,D2,GG,IPL,F,Sig,DSPDPEq,DEpsPEqDPS,EpsP,c,phi,psi)
        !**********************************************************************
        !
        ! Final correction of the yield surface drift (END OF STEP CORRECTION).
        ! The stresses, the plastic strain and the strength parameters are corrected.
        !
        !**********************************************************************

        implicit none

        !Local variables
        integer :: i
        double precision :: p,J,Lode,S3TA !Invariants
        double precision :: Lambda,param,c2,phi2,psi2,F2
        double precision :: Denom,A,Ai
        double precision, dimension(2) :: DFDSP
        double precision, dimension(6) :: DPPDSig,DFDSig,Sig2,DEpsP,EpsP2
        double precision, dimension(6) :: Denom1
        double precision, dimension(3) :: Dh
        !In Variables
        integer, intent(in) :: IntGlo,IPL !Global ID of Gauss point or particle
        double precision, intent(in):: D1,D2,GG
        double precision, intent(in), dimension(3) :: DSPDPEq !Derivatives respect Equivalent Plastic Strain
        double precision, intent(in), dimension(6) :: DEpsPEqDPS !Derivatives respect Equivalent Plastic Strain
        !InOut Variables
        double precision, intent(inout):: c,phi,psi
        double precision, intent(inout), dimension(6) :: Sig
        double precision, intent(inout), dimension(6) :: EpsP
        double precision, intent(inout):: F

        call CalculateInvariants(IntGlo,Sig,p,J,Lode,S3TA)
        call CalculateDerivativesYieldFunctAndPlasticPotential(Sig,p,J,Lode,S3TA,c,phi,psi,DFDSig,DPPDSig)
        call CalculateDerivativesYieldFunctSofteningParameters(p,J,Lode,S3TA,c,phi,DFDSP)

        !Parameter A (hardening/softening parameter)
        A = 0.0d0
        Ai = (DFDSP(1)*DSPDPEq(1) + DFDSP(2)*DSPDPEq(2))
        do i=1,6
            A = A + Ai * DEpsPEqDPS(i) * DPPDSig(i)
        end do

        Denom1(1) = DPPDSig(1)*D1 + DPPDSig(2)*D2 + DPPDSig(3)*D2
        Denom1(2) = DPPDSig(1)*D2 + DPPDSig(2)*D1 + DPPDSig(3)*D2
        Denom1(3) = DPPDSig(1)*D2 + DPPDSig(2)*D2 + DPPDSig(3)*D1
        Denom1(4) = DPPDSig(4)*GG
        Denom1(5) = DPPDSig(5)*GG
        Denom1(6) = DPPDSig(6)*GG

        Denom = Denom1(1)*DFDSig(1) + Denom1(2)*DFDSig(2) + &
                Denom1(3)*DFDSig(3) + Denom1(4)*DFDSig(4) + &
                Denom1(5)*DFDSig(5) + Denom1(6)*DFDSig(6) - A

        Lambda = F/Denom !factor correction

        Sig2 = Sig - Lambda * Denom1 ! Sig2 = Sig + fact * Denom1 Stress corrected
        DEpsP = Lambda * DPPDSig
        EpsP2 = EpsP + DEpsP

        if (IPL == 1)then
            Dh = 0.0d0
        else
            param = DEpsPEqDPS(1) * DEpsP(1) + DEpsPEqDPS(2) * DEpsP(2) + DEpsPEqDPS(3) * DEpsP(3) + &
                DEpsPEqDPS(4) * DEpsP(4) + DEpsPEqDPS(5) * DEpsP(5) + DEpsPEqDPS(6) * DEpsP(6)
            Dh(1) = min (DSPDPEq(1)*param, 0.0d0)
            Dh(2) = min (DSPDPEq(2)*param, 0.0d0)
            Dh(3) = min (DSPDPEq(3)*param, 0.0d0)
        end if

        c2 = c + Dh(1)
        phi2 = phi + Dh(2)
        psi2 = psi + Dh(3)

        call DetermineYieldFunctionValue(IntGlo,Sig2,c2,phi2,F2)
        
        if ((abs(F2) > abs(F)).or.(Denom == 0.0d0)) then !NormalCorrectionScheme
            Denom = 0.0d0
            Denom = DFDSig(1)*DFDSig(1) + DFDSig(2)*DFDSig(2) + &
                    DFDSig(3)*DFDSig(3) + DFDSig(4)*DFDSig(4) + &
                    DFDSig(5)*DFDSig(5) + DFDSig(6)*DFDSig(6)

            Lambda = F/Denom
            Sig = Sig - Lambda * DFDSig
            DEpsP = Lambda * DPPDSig
            EpsP = EpsP + DEpsP
            call DetermineYieldFunctionValue(IntGlo,Sig,c,phi,F)
        else
            Sig = Sig2
            EpsP = EpsP2
            c = c2
            phi = phi2
            psi = psi2
            F = F2
        end if

      end subroutine EndOfStepCorrection

      Subroutine CalculatePrincipalStresses(IntGlo,Sig,SigPrin)
        !**********************************************************************
        !
        ! Implemented in the frame of the MPM project.
        !
        !**********************************************************************

        implicit none

        !Local variables
        double precision, dimension(3) :: xN1,xN2,xN3
        double precision :: Sig1,Sig2,Sig3,p,q
        !In Variables
        integer, intent(in) :: IntGlo ! Global ID of Gauss point or particle
        double precision, intent(in), dimension(6) :: Sig
        !Out Variables
        double precision, intent(out), dimension(6) :: SigPrin

        call PrincipalSig(1,Sig,xN1,xN2,xN3,Sig1,Sig2,Sig3,P,Q)

        If (Sig1 >= Sig2.and.Sig2 >= Sig3) then
            SigPrin(1) = Sig1
            SigPrin(2) = Sig2
            SigPrin(3) = Sig3
        else if (Sig1 >= Sig3.and.Sig3 >= Sig2) then
            SigPrin(1) = Sig1
            SigPrin(2) = Sig3
            SigPrin(3) = Sig2
        else if (Sig3 >= Sig1.and.Sig1 >= Sig2) then
            SigPrin(1) = Sig3
            SigPrin(2) = Sig1
            SigPrin(3) = Sig2
        else if (Sig3 >= Sig2.and.Sig2 >= Sig1) then
            SigPrin(1) = Sig3
            SigPrin(2) = Sig2
            SigPrin(3) = Sig1
        else if (Sig2 >= Sig1.and.Sig1 >= Sig3) then
            SigPrin(1) = Sig2
            SigPrin(2) = Sig1
            SigPrin(3) = Sig3
        else if (Sig2 >= Sig3.and.Sig3 >= Sig1) then
            SigPrin(1) = Sig2
            SigPrin(2) = Sig3
            SigPrin(3) = Sig1
        end if

            SigPrin(4) = 0.0d0
            SigPrin(5) = 0.0d0
            SigPrin(6) = 0.0d0

      end subroutine CalculatePrincipalStresses

      Subroutine PrincipalSig(IOpt,S,xN1,xN2,xN3,S1,S2,S3,P,Q)
        Implicit Double Precision (A-H,O-Z)
        Dimension S(*),xN1(*),xN2(*),xN3(*)

        If (iOpt.Eq.1) Then
            Call Eig_3_MohrCoulombStrainSoftening(0,S,xN1,xN2,xN3,S1,S2,S3,P,Q) ! with Eigenvectors
        Else
            Call Eig_3a_MohrCoulombStrainSoftening(0,S,S1,S2,S3,P,Q) ! no Eigenvectors
        End If
        Return
      end subroutine PrincipalSig
      
      Subroutine Eig_3_MohrCoulombStrainSoftening(iOpt,St,xN1,xN2,xN3,S1,S2,S3,P,Q)
        Implicit Double Precision (A-H,O-Z)
        Dimension St(6),A(3,3),V(3,3),xN1(3),xN2(3),xN3(3)
        !     *          xN1(3),xN2(3),xN3(3)
        !
        ! Get Eigenvalues/Eigenvectors for 3*3 matrix
        ! Wim Bomhof 15/11/'01
        ! PGB : adaption to Principal stress calculation
        !
        ! Applied on principal stresses, directions
        ! Stress vector St(): XX, YY, ZZ, XY, YZ, ZX
        !
        A(1,1) = St(1) ! xx
        A(1,2) = St(4) ! xy = yx
        A(1,3) = St(6) ! zx = xz

        A(2,1) = St(4) ! xy = yx
        A(2,2) = St(2) ! yy
        A(2,3) = St(5) ! zy = yz

        A(3,1) = St(6) ! zx = xz
        A(3,2) = St(5) ! zy = yz
        A(3,3) = St(3) ! zz

        ! Set V to unity matrix
        V(1,1) = 1
        V(2,1) = 0
        V(3,1) = 0

        V(1,2) = 0
        V(2,2) = 1
        V(3,2) = 0

        V(1,3) = 0
        V(2,3) = 0
        V(3,3) = 1


        abs_max_s=0.0
        Do i=1,3
            Do j=1,3
            if (abs(a(i,j)) .Gt. abs_max_s) abs_max_s=abs(a(i,j))
            End Do
        End Do
        Tol = 1d-20 * abs_max_s
        it = 0
        itmax = 50
        Do While ( it.Lt.itMax .And. abs(a(1,2))+abs(a(2,3))+abs(a(1,3)) .Gt. Tol )
        !     *           abs(a(1,2))+abs(a(2,3))+abs(a(1,3)) .Gt. Tol )
            it=it+1
            Do k=1,3
            If (k .Eq. 1) Then
                ip=1
                iq=2
            Else If (k .Eq.2) Then
                ip=2
                iq=3
            Else
                ip=1
                iq=3
            End If
            If (a(ip,iq) .Ne. 0.0) Then
                tau=(a(iq,iq)-a(ip,ip))/(2.0*a(ip,iq))
                If (tau .Ge.0.0) Then
                sign_tau=1.0
                Else
                sign_tau=-1.0
                End If
                t=sign_tau/(abs(tau)+sqrt(1.0+tau*tau))
                c=1.0/sqrt(1.0+t*t)
                s=t*c
                a1p=c*a(1,ip)-s*a(1,iq)
                a2p=c*a(2,ip)-s*a(2,iq)
                a3p=c*a(3,ip)-s*a(3,iq)
                a(1,iq)=s*a(1,ip)+c*a(1,iq)
                a(2,iq)=s*a(2,ip)+c*a(2,iq)
                a(3,iq)=s*a(3,ip)+c*a(3,iq)
                a(1,ip)=a1p
                a(2,ip)=a2p
                a(3,ip)=a3p

                v1p=c*v(1,ip)-s*v(1,iq)
                v2p=c*v(2,ip)-s*v(2,iq)
                v3p=c*v(3,ip)-s*v(3,iq)
                v(1,iq)=s*v(1,ip)+c*v(1,iq)
                v(2,iq)=s*v(2,ip)+c*v(2,iq)
                v(3,iq)=s*v(3,ip)+c*v(3,iq)
                v(1,ip)=v1p
                v(2,ip)=v2p
                v(3,ip)=v3p

                ap1=c*a(ip,1)-s*a(iq,1)
                ap2=c*a(ip,2)-s*a(iq,2)
                ap3=c*a(ip,3)-s*a(iq,3)
                a(iq,1)=s*a(ip,1)+c*a(iq,1)
                a(iq,2)=s*a(ip,2)+c*a(iq,2)
                a(iq,3)=s*a(ip,3)+c*a(iq,3)
                a(ip,1)=ap1
                a(ip,2)=ap2
                a(ip,3)=ap3
            End If ! a(ip,iq)<>0
            End Do ! k
        End Do ! While
        ! principal values on diagonal of a
        S1 = a(1,1)
        S2 = a(2,2)
        S3 = a(3,3)
        ! Derived invariants
        P = (S1+S2+S3)/3
        Q = Sqrt( ( (S1-S2)**2 + (S2-S3)**2 + (S3-S1)**2 ) / 2 )

        ! Sort eigenvalues S1 <= S2 <= S3
        is1 = 1
        is2 = 2
        is3 = 3
        if (s1.Gt.s2) Then
            t   = s2
            s2  = s1
            s1  = t
            it  = is2
            is2 = is1
            is1 = it
        End If
        if (s2.Gt.s3) Then
            t   = s3
            s3  = s2
            s2  = t
            it  = is3
            is3 = is2
            is2 = it
        End If
        if (s1.Gt.s2) Then
            t   = s2
            s2  = s1
            s1  = t
            it  = is2
            is2 = is1
            is1 = it
        End If
        Do i=1,3
            xN1(i) = v(i,is1) ! first  column
            xN2(i) = v(i,is2) ! second column
            xN3(i) = v(i,is3) ! third  column
        End Do
        Return
      end Subroutine Eig_3_MohrCoulombStrainSoftening ! Eig_3

      Subroutine Eig_3a_MohrCoulombStrainSoftening(iOpt,St,S1,S2,S3,P,Q) ! xN1,xN2,xN3,
        Implicit Double Precision (A-H,O-Z)
        Dimension St(6),A(3,3)   !  V(3,3),xN1(3),xN2(3),xN3(3)
        !
        ! Get Eigenvalues ( no Eigenvectors) for 3*3 matrix
        ! Wim Bomhof 15/11/'01
        !
        ! Applied on principal stresses, directions
        ! Stress vector XX, YY, ZZ, XY, YZ, ZX
        !
        A(1,1) = St(1) ! xx
        A(1,2) = St(4) ! xy = yx
        A(1,3) = St(6) ! zx = xz

        A(2,1) = St(4) ! xy = yx
        A(2,2) = St(2) ! yy
        A(2,3) = St(5) ! zy = yz

        A(3,1) = St(6) ! zx = xz
        A(3,2) = St(5) ! zy = yz
        A(3,3) = St(3) ! zz

        abs_max_s=0.0
        Do i=1,3
            Do j=1,3
            if (abs(a(i,j)) .Gt. abs_max_s) abs_max_s=abs(a(i,j))
            End Do
        End Do
        Tol = 1d-20 * abs_max_s
        If (iOpt.Eq.1) Tol = 1d-50*abs_max_s
        it=0
        itmax = 50
        Do While ( it.lt.itmax .And.&
                    abs(a(1,2))+abs(a(2,3))+abs(a(1,3)) .Gt. Tol )

            it=it+1
            Do k=1,3
            If (k .Eq. 1) Then
                ip=1
                iq=2
            Else If (k .Eq.2) Then
                ip=2
                iq=3
            Else
                ip=1
                iq=3
            End If
            If (a(ip,iq) .Ne. 0.0) Then         ! ongelijk nul ?
                tau=(a(iq,iq)-a(ip,ip))/(2.0*a(ip,iq))
                If (tau .Ge.0.0) Then
                sign_tau=1.0
                Else
                sign_tau=-1.0
                End If
                t=sign_tau/(abs(tau)+sqrt(1.0+tau*tau))
                c=1.0/sqrt(1.0+t*t)
                s=t*c
                a1p=c*a(1,ip)-s*a(1,iq)
                a2p=c*a(2,ip)-s*a(2,iq)
                a3p=c*a(3,ip)-s*a(3,iq)
                a(1,iq)=s*a(1,ip)+c*a(1,iq)
                a(2,iq)=s*a(2,ip)+c*a(2,iq)
                a(3,iq)=s*a(3,ip)+c*a(3,iq)
                a(1,ip)=a1p
                a(2,ip)=a2p
                a(3,ip)=a3p

                ap1=c*a(ip,1)-s*a(iq,1)
                ap2=c*a(ip,2)-s*a(iq,2)
                ap3=c*a(ip,3)-s*a(iq,3)
                a(iq,1)=s*a(ip,1)+c*a(iq,1)
                a(iq,2)=s*a(ip,2)+c*a(iq,2)
                a(iq,3)=s*a(ip,3)+c*a(iq,3)
                a(ip,1)=ap1
                a(ip,2)=ap2
                a(ip,3)=ap3
            End If ! a(ip,iq)<>0
            End Do ! k
        End Do ! While
        ! principal values on diagonal of a
        S1 = a(1,1)
        S2 = a(2,2)
        S3 = a(3,3)
        ! Derived invariants
        P = (S1+S2+S3)/3
        Q = Sqrt( ( (S1-S2)**2 + (S2-S3)**2 + (S3-S1)**2 ) / 2 )

        if (s1.Gt.s2) Then
            t   = s2
            s2  = s1
            s1  = t
        End If
        if (s2.Gt.s3) Then
            t   = s3
            s3  = s2
            s2  = t
        End If
        if (s1.Gt.s2) Then
            t   = s2
            s2  = s1
            s1  = t
        End If
        Return
      end Subroutine Eig_3a_MohrCoulombStrainSoftening

      Subroutine MatVec_MohrCoulombStrainSoftening(xMat,IM,Vec,N,VecR)
        !C***********************************************************************
        !C
        !C     Calculate VecR = xMat*Vec
        !C
        !C I   xMat  : (Square) Matrix (IM,*)
        !C I   Vec   : Vector
        !C I   N     : Number of rows/colums
        !C O   VecR  : Resulting vector
        !C
        !C***********************************************************************
            Implicit Double Precision (A-H,O-Z)
            Dimension xMat(IM,*),Vec(*),VecR(*)
        !C***********************************************************************
            Do I=1,N
                X=0
                Do J=1,N
                X=X+xMat(I,J)*Vec(J)
                End Do
                VecR(I)=X
            End Do
            Return
      end Subroutine MatVec_MohrCoulombStrainSoftening    ! Subroutine MatVec

      Subroutine AddVec_MohrCoulombStrainSoftening(Vec1,Vec2,R1,R2,N,VecR)
        !C***********************************************************************
        !C
        !C     Calculate VecR() = R1*Vec1()+R2*Vec2()
        !C
        !C I   Vec1,
        !C I   Vec2  : Vectors
        !C I   R1,R2 : Multipliers
        !C I   N     : Number of rows
        !C O   VecR  : Resulting vector
        !C
        !C***********************************************************************
            Implicit Double Precision (A-H,O-Z)
            Dimension Vec1(*),Vec2(*),VecR(*)
        !C***********************************************************************
            Do I=1,N
                X=R1*Vec1(I)+R2*Vec2(I)
                VecR(I)=X
            End Do
            Return
      end Subroutine AddVec_MohrCoulombStrainSoftening    ! Subroutine AddVec


end module ModExternalSoilModel