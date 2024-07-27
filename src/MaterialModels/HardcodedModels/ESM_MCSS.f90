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
!   and soil�water�structure interaction using the material point method (MPM)
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

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 !███╗░░░███╗░█████╗░░██████╗░██████╗
 !████╗░████║██╔══██╗██╔════╝██╔════╝
 !██╔████╔██║██║░░╚═╝╚█████╗░╚█████╗░
 !██║╚██╔╝██║██║░░██╗░╚═══██╗░╚═══██╗
 !██║░╚═╝░██║╚█████╔╝██████╔╝██████╔╝
 !╚═╝░░░░░╚═╝░╚════╝░╚═════╝░╚═════╝░
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


module MOD_MCSS_ESM
   !**********************************************************************
   !
   ! Module: contains all functions and subroutines required for the Mohr-Coloumb Strain Softening constitutive model
   !
   ! Note: Integer type and real type is not specified in this module. This was done because the ESMs are usually compiled into external .dlls
   !        and the type wouldn't be specfied there
   ! TODO: Add integer and real type. Also determine if doubles are really needed for this computation. As they are cast into reals does it matter?
   !
   !     $Revision: ????? $
   !     $Date: 2023-12-28 11:41 +0500 (WaveHello, 28 Dec 2023) $
   !
   !**********************************************************************
   !TODO: Figure out why only this module created a .mod file
   implicit none
   private ! Makes all function private to this module (No other modules can get access)
   public ESM_MohrCoulombStrainSoftening ! Overides private status for specific subroutine

contains
   SUBROUTINE ESM_MohrCoulombStrainSoftening(NPT,NOEL,IDSET,STRESS,EUNLOADING,PLASTICMULTIPLIER,&
      DSTRAN,NSTATEV,STATEV,NADDVAR,ADDITIONALVAR,CMNAME,NPROPS,PROPS,NUMBEROFPHASES,NTENS)

      CHARACTER*80 CMNAME
      integer :: NPT, NOEL, IDSET, NSTATEV, NADDVAR, NPROPS, NUMBEROFPHASES, NTENS
      integer :: IStep, TimeStep
      double precision :: Eunloading, PLASTICMULTIPLIER
      double precision :: STRESS(NTENS), DSTRAN(NTENS), STATEV(NSTATEV), ADDITIONALVAR(NADDVAR), &
         PROPS(NPROPS)

      !---Local variables required in standard UMAT
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


      integer :: ndi, nshr, layer, kspt, kstep, kinc, IDTask



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

   SUBROUTINE UMAT_MohrCoulombStrainSoftening(STRESS, STATEV, DDSDDE, SSE, SPD, SCD, &
      RPL, DDSDDT, DRPLDE, DRPLDT, STRAN, DSTRAN, TIME, DTIME, TEMP,&
      DTEMP, PREDEF, DPRED, CMNAME, NDI, NSHR, NTENS, NSTATEV, PROPS,&
      NPROPS, COORDS, DROT, PNEWDT, CELENT, DFGRD0, DFGRD1, NOEL, NPT,&
      LAYER, KSPT, KSTEP, KINC)

      CHARACTER*80 CMNAME

      double precision :: SSE, SPD, SCD ! specific elastic strain energy, plastic dissipation, creep dissipation
      double precision :: rpl ! only for fully coupled thermal analysis: volumetric heat generation
      double precision :: drpldt ! only for fully coupled thermal analysis: variation of volumetric heat generation due to temperature
      double precision :: pnewdt, dtime, temp, dtemp, celent
      integer :: ndi, nshr, layer, kspt, kstep, kinc, ntens, nstatev, nprops, noel, &
         npt, ipl, intGlo

      double precision :: STRESS(NTENS), STATEV(NSTATEV), DDSDDE(NTENS,NTENS), &
         DDSDDT(NTENS),DRPLDE(NTENS),STRAN(NTENS),DSTRAN(NTENS),&
         TIME(2),PREDEF(1),DPRED(1), PROPS(NPROPS),COORDS(3),DROT(3,3),&
         DFGRD0(3,3),DFGRD1(3,3)


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

      double precision :: Rad ! Store the value of a radian

      !Soil properties
      double precision :: G, ENU, cp, cr, phip, phir, psip, psir, factor, c, phi, psi, &
         F1, F2
      ! Stress variables
      double precision :: p, YTOL, Euler_DT_min
      
      ! Define elastic matrix (DE), stress increment, stress, plastic strain incremenent, and total plastic strain
      double precision :: DE(6,6), dSig(6), Sig(6), dEpsP(6), EpsP(6)

      integer :: i, integration_flag, num_integration_iters
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
      ! 10 : integration_flag - Controls the stress integration scheme 0 for Euler and 1 for Ortiz-Simo
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
      integration_flag = PROPS(10) ! Stress Integration flag
      YTOL = PROPS(11) ! Yield surface tolerance
      num_integration_iters = PROPS(12) ! Number of stress integration iterations allowed
      Euler_DT_min = PROPS(13) ! Minimum pseudo-time step for Sloan integration
      
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
      Sig = STRESS + dSig

      call MOHRStrainSoftening(IntGlo,F1,F2,G,cp,cr,phip,phir,psip,psir,factor,c,phi,psi,stress,EpsP,&
                                DSTRAN,dEpsP,Sig,IPL, integration_flag, YTOL, num_integration_iters, Euler_DT_min)

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

   Subroutine MOHRStrainSoftening(IntGlo,D1,D2, GG,cp,cr,phip,phir, &
      psip,psir,factor,c,phi,psi,Sig0,EpsP,DEps,DEpsP,SigC,IPL, integration_flag, YTOL, &
       num_integration_iters, Euler_DT_min)
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
      integer, parameter :: ONE =1, ZERO = 0
      double precision :: F,F0,F2 !Evaluation of the Yield function
      double precision :: alpha !Elastic Strain proportion
      double precision :: SSTOL !Tolerance Relative Error
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
      double precision, dimension(6) :: DSigPP,DSigP1,DSigP2, DSIGE
      double precision, dimension(6) :: DEpsPP,DEpsPP1,DEpsPP2
      double precision, dimension(6) :: DEpsS,DEpsSS
      double precision, dimension(6) :: EpsP1,EpsP2
      double precision, dimension(6) :: DEpsPEqDPS,DEpsPEqDPS1
      double precision, dimension(6) :: sumSg,Er
      double precision, dimension(3) :: DSPDPEq,DSPDPEq1 !Variation of softening parameters (c,phi,psi) in function of plastic strain
      !In variables
      integer, intent(in) :: IntGlo !Global ID of Gauss point or particle
      integer, intent(in) :: integration_flag, num_integration_iters
      double precision, intent(in) :: D1,D2,GG !Elastic Parameters
      double precision, intent(in) :: cp,cr,phip,phir,psip,psir,factor !Softening parameter
      double precision, intent(in) :: YTOL !Tolerance Error on the Yield surface (10-6 to 10-9)
      double precision, intent(in) :: Euler_DT_min ! Minimum psuedo-time step size
      
      !Inout variables
      double precision, intent(inout), dimension(6) :: DEps !Incremental total strain
      double precision, intent(inout):: c,phi,psi !cohesion,friction angle and dilatancy angle
      double precision, intent(inout), dimension(6) :: EpsP !Accumulated Plastic Strain
      double precision, intent(inout), dimension(6) :: Sig0 !Initial Stress
      double precision, intent(inout), dimension(6) :: SigC !Final Stress
      double precision, intent(inout), dimension(6) :: DEpsP !Incremental plastic strain
    
      !Out variables
      integer, intent(out) :: IPL

      if (num_integration_iters == 0) then
          print *, "Number of MCSS stress iterations set to zero (Param 12 in .GOM file)"
      end if
      
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
      SPTOL = 0.01d0 !Tolerance Softening Parameters (0.0001d0)
      ctol = abs(cp-cr)*SPTOL
      phitol = abs(phip-phir)*SPTOL
      psitol = abs(psip-psir)*SPTOL
      
      if (Euler_DT_min <= 1e-10) then
          print *, "Input DTmin, MCSS param 13 is less then 1e-10, that's too small, DTmin set to 1e-8"
          DTmin = 1e-8
      else
        DTmin = Euler_DT_min
      endif
        
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

      select case(integration_flag)

      case(ZERO)
         !Determine the proportion (alpha) of the stress increment that lies within the yield function.
         !The PEGASUS ALGORITHM SCHEME FOR CONVENTIONAL ELASTOPLASTIC MODELS has been used
         call DetermineYieldFunctionValue(IntGlo,Sig0,c,phi,F0)

         if (F0 < -YTOL) then !In this Time increment there is part of elastic behavior
            call DetermineElasticProportionPegasusMethod(IntGlo,Sig0,DSigE,DEps,c,phi,YTOL,alpha)
         else
            alpha = 0.0d0 !In this time increment all is plastic behavior
         end if

         !Calculate the direction of the stress path--> missing
         !It is assumed that the direction is always outside the yield surface.

         !Determine the elastic portion of the stress increment
         DSigE = alpha * DSigE !Correct Incremental Elastic Stress

         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         !Determine the plastic portion of the stress increment.
         !The method used is the MODIFIED EULER INTEGRATION SCHEME with error control
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

         !Initialise parameters
         SigYield = Sig0 + DSigE !Sigma on the Yield surface
         DEpsS = (1.0d0-alpha) * DEps !Incremental Plastic Strain

         T = 0.0d0
         DT = 1.0d0

         !Start the plastification
         Do while (T <= 1.0d0)
            m = 0 !Counter
            Rn = 100

            call CalculateEpsPEq(EpsP,EpsPEq) !Determine Equivalent plastic Strain (EpsPEq)

            Do while (Rn > SSTOL.and.m < num_integration_iters)
               !1)Calculation of the portion of the plastic strain increment (DEpsPP)
               DEpsSS = DT * DEpsS !Portion of the plastic strain increment

               !Calculate a first estimate of the associated stress
               !hardening/softening parameter changes
               call CalculateDerivativesStrSoftParamRespectEquivalentPlasticStrain(factor,cp,cr,phip,phir,psip,psir,&
                  EpsPEq,DSPDPEq)
               call CalculateDerivativesEquivalentPlasticStrainRespectPlasticStrain(EpsP,EpsPEq,DEpsPEqDPS)
               call DetermineDSigAndDEpsP(IntGlo,D1,D2,GG,c,phi,psi,SigYield,DEpsPEqDPS,DSPDPEq,DEpsSS,DSigP1,DEpsPP1)
               EpsP1 = EpsP + DEpsPP1

               call CalculateEpsPEq(EpsP1,EpsPEq1) !Determine Equivalent plastic Strain (EpsPEq)

               !if (IPL == 1) then !Residual conditions --> no changes of the strength parameters
               !    c1 = c
               !    phi1 = phi
               !    psi1 = psi
               !else !IPL=2 Softening Conditions --> changes of the strength parameters
               call CalculateSofteningParameters(EpsPEq1,factor,cp,cr,phip,phir,psip,psir,c1,phi1,psi1)
               !end if

               !2)Calculate a second estimate of the associated stress
               !hardening/softening parameter changes
               SigYield2 = SigYield + DSigP1

               call CalculateDerivativesStrSoftParamRespectEquivalentPlasticStrain(factor,cp,cr,phip,phir,psip,psir,&
                  EpsPEq1,DSPDPEq1)
               call CalculateDerivativesEquivalentPlasticStrainRespectPlasticStrain(EpsP1,EpsPEq1,DEpsPEqDPS1)
               call DetermineDSigAndDEpsP(IntGlo,D1,D2,GG,c1,phi1,psi1,SigYield2,DEpsPEqDPS1,DSPDPEq1,DEpsSS,DSigP2,DEpsPP2)
               EpsP2 = EpsP + DEpsPP2

               call CalculateEpsPEq(EpsP2,EpsPEq2) !Determine Equivalent plastic Strain (EpsPEq)

               !if (IPL == 1) then !Residual conditions --> no changes of the strength parameters
               !    c2 = c
               !    phi2 = phi
               !    psi2 = psi
               !else  !IPL=2 Softening Conditions --> changes of the strength parameters
               call CalculateSofteningParameters(EpsPEq2,factor,cp,cr,phip,phir,psip,psir,c2,phi2,psi2)
               !end if

               !3)Obtain a more accurate modified Euler estimate of the changes in stress,
               !plastic strain and hardening/softening parameters
               DSigPP = 0.5d0 * (DSigP1 + DSigP2)

               !Calculation of the relative error
               Er = 0.5d0 * (DSigP1 - DSigP2)
               moduleEr = sqrt(Er(1)*Er(1)+Er(2)*Er(2)+ Er(3)*Er(3)+ Er(4)*Er(4)+Er(5)*Er(5)+Er(6)*Er(6))

               sumSg = SigYield + DSigPP
               moduleSigDSig = sqrt(sumSg(1)*sumSg(1) + sumSg(2)*sumSg(2) + sumSg(3)*sumSg(3)+ &
                  sumSg(4)*sumSg(4) + sumSg(5)*sumSg(5) + sumSg(6)*sumSg(6))

               !Check the relative error (Rn) of the new stresses, with the defined tolerance (SSTOL)
               Rn = (moduleEr/moduleSigDSig)

               ! check whether decreasing of DT is possible, if not exit loop
               if (DT == DTmin) then
                  exit
               end if

               !4)If Rn>SSTOL, the loop is not finished and the substep is recalculated smaller
               if (Rn > SSTOL) then
                  beta = max (0.9d0*(sqrt(SSTOL/Rn)), 0.1d0)
                  beta = min (beta, 1.1d0)
                  DT = max (DT*beta, DTmin)
                  m = m + 1 ! Update counter
               end if

            end do

            !Update the accumulated stresses, plastic strain and softening parameters
            SigYield = SigYield + DSigPP
            DEpsPP = 0.5d0 * (DEpsPP1 + DEpsPP2)
            DEpsP = DEpsP + DEpsPP
            EpsP = EpsP + DEpsPP

            call CalculateEpsPEq(EpsP,EpsPEq) !Determine Equivalent plastic Strain (EpsPEq)

            call CalculateSofteningParameters(EpsPEq,factor,cp,cr,phip,phir,psip,psir,c,phi,psi)

            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            !!!!!!!!!!!!!!!!!!!!!!! END OF STEP CORRECTION!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            !Check if we are on/in the yield surface, otherwise we are still outside (F>0)
            !and a correction is needed.
            call DetermineYieldFunctionValue(IntGlo,SigYield,c,phi,F)
            n=0 !Counter
            do while (abs(F) > YTOL.and.n < 10) !The correction is needed
               n = n + 1
               call CalculateEpsPEq(EpsP,EpsPEq)             !Determine Equivalent plastic Strain (EpsPEq)
               call CalculateDerivativesStrSoftParamRespectEquivalentPlasticStrain(factor,cp,cr,phip,phir,psip,psir,&
                  EpsPEq,DSPDPEq)
               call CalculateDerivativesEquivalentPlasticStrainRespectPlasticStrain(EpsP,EpsPEq,DEpsPEqDPS)
               call EndOfStepCorrection(IntGlo,D1,D2,GG,IPL,F,SigYield,DSPDPEq,DEpsPEqDPS,EpsP,c,phi,psi)
            end do
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

            !The substep is updated
            T1 = T + DT

            !If T1>1 the calculation is finished
            If (T1 >= 1d0) then
               SigC = SigYield   !Determine Final stresses
               return
            end if

            !If T1<1, calculation of the next substep DT
            beta = min (0.9d0*(sqrt(SSTOL/Rn)), 1.1d0)
            if (m > 1) then ! the previous step failed
               beta = min (beta, 1.0d0)
               DT = beta * DT
               it = it+1
            else
               DT = beta * DT
               it = 0
            end if
            DT = max (DT, DTmin)
            DT = min (DT, 1.0d0-T1)
            T = T1

         end do  !If T=1 the loop is finishedv   !Determine the proportion (alpha) of the stress increment that lies within the yield function.
         
       case(ONE)

         call MCSS_Ortiz_Simo_Integration(GG, D1, D2, IntGlo, Sig0, c, phi, psi, factor, DEps, EpsP, dEpsP, cr, phir, psir, cp, phip, &
            psip, ctol, phitol, psitol, YTOL, num_integration_iters)

         ! State parameters {phi, psi, c} updated inside ortiz-simo
         ! EpsP updated inside of the integration
         ! dEpsP updated inside of ortiz-Simo

         ! Update the Stress
         SigC = Sig0

         ! Increment of elastic stress not updated
      end select
   end subroutine MOHRStrainSoftening

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
   
   Subroutine DetermineElasticProportionPegasusMethod(IntGlo,Sig,DSig,DEps,c,phi,YTOL,alpha)
      !**********************************************************************
      !
      ! The PEGASUS METHOD method is used
      !
      !**********************************************************************

      implicit none

      !Local variables
      integer :: Its, max_iterations
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
      max_iterations = 1000 ! Maximum newton iterations
      
      do while (abs(F) > YTOL.and.Its < max_iterations)
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
      if (Its >= max_iterations) then
         alpha = 0.0d0
      end if
   end subroutine DetermineElasticProportionPegasusMethod

   Subroutine MCSS_Ortiz_Simo_Integration(G, D1, D2, IntGlo, Sig, c, phi, psi, factor, dEps, EpsP, dEpsP, cr, phir, psir, cp, phip, &
      psip, ctol, phitol, psitol, FTOL, max_iterations)
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
      ! max_iterations: The maximum number of time the gradient descent  method should be used

      ! -------------------------- Variable Definitions --------------------------
      ! ------------- Scalar Values -------------
      ! In
      double precision, intent(in) :: G, D1, D2, factor, cr, phir, psir, cp, phip, psip, ctol, phitol, psitol, FTOL
      integer, intent(in) :: IntGlo, max_iterations

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
      double precision :: p
      integer:: counter

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
      dSigu = matmul(DE, dEps)

      ! Update the stresses
      Sigu = Sigu + dSigu

      ! Evalue the yield surface
      call DetermineYieldFunctionValue(IntGlo, Sigu, cu, phiu, F)

      if (F <= FTOL) then
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
      counter = 0

      do while(abs(F) >= FTOL .and. counter <= max_iterations)
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
         DE_m = matmul(DE, m_vec)

         ! Compute n_vec.DE.m_vec
         dummyVal_1 = dot_product(n_vec, DE_m)

         ! Make a 1x3 vector to store dF/dXs
         dummyVec_3d(:) = 0
         dummyVec_3d(1) = dFdSP(1)
         dummyVec_3d(2) = dFdSP(2)

         ! Calc the dot product between dF/dXs.dXs/dEpsPEq
         dummyVal_2 = dot_product(dummyVec_3d, DSPDPEq)

         ! Calc the dot product between dEpsPEq/dEpsP.dP/dSig
         dummyVal_3 = dot_product(DEpsPEqDPS, m_vec)

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
      ! TODO: Store rounded params in module array and pass those values to the
      ! associated functions
      ! There's a couple functions that use these values
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
      ! TODO: shorten code length store-> -factor * (exp(-factor*EpsPEq)) as
      !       softening factor
      !**********************************************************************

      implicit none

      ! In Variables
      double precision, intent(in) :: EpsPEq
      double precision, intent(in) :: factor,cp,cr,phip,phir,psip,psir
      ! Out Variables
      double precision, intent(out), dimension(3):: DSPDPEq

      ! Derivative Cohesion respect Equivalent Plastic Strain (Dc/DPEq)
      DSPDPEq(1) = -factor * (cp - cr) * (exp(-factor*EpsPEq))
      ! Derivative Friction angle respect Equivalent Plastic Strain (Dphi/DPEq)
      DSPDPEq(2) = -factor * (phip - phir) * (exp(-factor*EpsPEq))
      ! Derivative Dilatancy angle respect Equivalent Plastic Strain (Dpsi/DPEq)
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
      ! TODO: Turn the long calculation into a matrix vector product
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

   Subroutine PrincipalSig(IOpt, S, xN1, xN2, xN3, S1, S2, S3, P, Q)
      implicit none

      !! TODO: Need to assign intents for each variable

      integer :: iOPt
      ! Pricipal stress values
      double precision :: s1, s2, s3
      ! Stress invariants
      double precision :: P, Q
      double precision ::  S(:), xN1(:), xN2(:), xN3(:)

      If (iOpt.Eq.1) Then
         ! Calculates Eigenvalues and eigenvectors
         Call Eig_3_MohrCoulombStrainSoftening(0, S, xN1, xN2, xN3, S1, S2, S3, P, Q)
      Else

         Call Eig_3a_MohrCoulombStrainSoftening(0, S, S1, S2, S3, P, Q) ! no Eigenvectors
      End If
      Return
   end subroutine PrincipalSig

   Subroutine Eig_3_MohrCoulombStrainSoftening(iOpt,St,xN1,xN2,xN3,S1,S2,S3,P,Q)
      implicit none
      !TODO: Remove iO

      integer :: iOpt, i, j, k, it, itmax, is3, is2, is1, iq, ip, sign_tau
      ! ^ Optional Integer value (Not used)
      double precision :: v3p, v2p, v1p, t, tol, tau, s, s1, s2, s3, c, ap1, ap2, ap3, abs_max_s, a1p, a2p, a3p
      ! Stress invariants p, q
      double precision :: P, Q
      double precision :: St(6),A(3,3),V(3,3),xN1(3),xN2(3),xN3(3)
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
      implicit none

      integer :: iOpt

      ! Principal stress values
      double precision :: s1, s2, s3

      ! Stress invariants
      double precision :: P, Q

      double precision :: St(6),A(3,3)   !  V(3,3),xN1(3),xN2(3),xN3(3)

      ! Local Variables
      integer :: i, j, k, it, itmax, is3, is2, is1, iq, ip, sign_tau
      double precision :: v3p, v2p, v1p, t, tol, tau, s, c, ap1, ap2, ap3, abs_max_s, a1p, a2p, a3p

      ! Get Eigenvalues ( no Eigenvectors) for 3*3 matrix
      ! Wim Bomhof 15/11/'01
      !
      ! Applied on principal stresses, directions
      ! Stress vector XX, YY, ZZ, XY, YZ, ZX

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
end module MOD_MCSS_ESM
