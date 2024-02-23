module MOD_ESM_Bardesie
   implicit none

contains
!    BARODESY FOR CLAY
!    Copyright (C) 2020  M. Bode, B. Schneider-Muntau, W. Fellin
!
!    This program is free software: you can redistribute it and/or modify
!    it under the terms of the GNU General Public License as published by
!    the Free Software Foundation, either version 3 of the License, or
!    (at your option) any later version.
!
!    This program is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!    GNU General Public License for more details.
!
!    You should have received a copy of the GNU General Public License
!    along with this program.  If not, see !<https://www.gnu.org/licenses/>.
!
!------------------------------------------------------------------------------
   subroutine umat(stress,statev,ddsdde,sse,spd,scd,&
      rpl,ddsddt,drplde,drpldt,&
      stran,dstran,time,dtime,temp,dtemp,predef,dpred,cmname,&
      ndi,nshr,ntens,nstatv,props,nprops,coords,drot,pnewdt,&
      celent,dfgrd0,dfgrd1,noel,npt,layer,kspt,jstep,kinc)
!------------------------------------------------------------------------------
! user subroutine for Abaqus (VERSION 2019)
!------------------------------------------------------------------------------
!
! Implemented constitutive model:
! -----------------------------
! Barodesie for clay
!
!
! The constitutive law is based on the article
! Medicus G. and Fellin W. (2017)
! An improved version of barodesy for clay
! Acta Geotechnica
!
! The Implementation of Strength Reduction is based on the article
! Bode et al (2020)
! Application of Barodesy - extended by the intergranular strain concept
! Proceedings of IACMAG conference 2020
!
!
! The Intergranular Strain extension is based on the article
! Bode et al.(2019)
! An intergranular strain concept for material models formulated as rate equations
! International Journal for Numerical and Analytical Methods in Geomechanics
!
!
! Second order Work calculation is based on the article
! Medicus et al (2018)
! Second order Work in barodesy
! Acta Geotechnica
!
!
! Asymptotic states and peak states are based on the article
! Medicus G. (2020)
! Asymptotic states and peak states in barodesy
! Geotechnic letters
!
! The implementation of the tangent opearator is based on the article
! Fellin, W. and Ostermann, A. (2002):
! Consistent tangent operators for constitutive rate equations.
! International Journal for Numerical and Analytical Methods in Geomechanics
!
!
!
! ----------------------------------------------------------------------------
! Material constants: n = index of material constants (props) in ABAQUS input file
!
!        ---------------------------------------------------------------------
!        n
!        ---------------------------------------------------------------------
!        1      phi_c
!        2      N
!        3      lambda
!        4      Kappa
!        5      Sigma*
!        6      m_R (intergranular strain)
!        7      m_T (intergranular strain)
!        8      R   (intergranular strain)
!        9      beta_r (intergranular strain)
!       10      chi  (intergranular strain)

!        ----------------------------------------------------------------------
!
! Solution dependent state variables (statev):
! definition via sdvini or initial conditions in input file
!
!        ---------------------------------------------------------------------
!        n
!        ---------------------------------------------------------------------
!
!        1	actual void ratio (input and output)
!        2	delta_11
!        3	delta_22
!        4	delta_33
!        5	delta_12
!        6	delta_13
!        7	delta_23
!        8  IS Mod (input)
!        9  SR Mod (input)
!  	  10  attempt for dtSub (output and input for next step)
!       11  p_t
!       12  mobilized friction angle (output)
!       13  2nd Invariant of Strain Tensor (output)
!       14  OCR p_e/p (output)
!       15  Dilatanz delta=trD/valD
!       16  Lode Angle
!       17  W2 second order work (output)
!       18  nSub
!       19  rho (magnitude of intergranular strain)
!       20  delta0:D0 (direcional change for IS)
!       21  EEM G
!       22  EEM nu
!       23  FoR factor of strength reduction (output)
!       24  phi_red (output)
!       25  phi_peak_red (output)
!       26  isASBS
!       27  isCorr (stress correction for SR)
!
!	----------------------------------------------------------------------
! Authors:

!     M. Bode, manuel.bode@uibk.ac.at
!     B. Schneider-Muntau, barbara.schneider-muntau@uibk.ac.at
!     W. Fellin
!     Devision of Geotechnical and Tunnel Engineering
!
!     University of Innsbruck
!
! This implementation uses ideas from:
!     H. Huegel (1995): estimation of a practical step size for
!                       forward Euler with constant step size
!                       the subroutine for evaluation of constitutive
!                       law based on his umat
!     D. Roddeman (1997): idea of reducing time substep when constitutive
!                         law is not defined
!     D. Masin: elastic stiffness for low stress states
!
! Last changes:
! 2/2006: Bug fix in subroutine evolut()
!         pertubation in y --> symmetric pertubation in D
! 9/2012 B. Schneider-Muntau
!       based on the subroutine for hypoplasticity the new constituve equation barodesy was implemented
! 12/2012 D. Masin
!       corrections for Barodesy, optimisation for PLAXIS
! 09/2015 B. Schneider-Muntau
!       new barodesy-Version Medicus 2015
! 12/2017 M. Bode
!       implementation of solution for zero stress level
!       restruction for better performance and clearness
!       optimisation for use in PLAXIS and ABAQUS
! 03/2018 M. Bode
!       Implementation of strength reduction and second order work
!       Additional output variables
! 05/2020 M.Bode
!       Bugfix in output variables
!       Implementation of the intergranular strain concept
!       Incoorporation of the ASBS
!       Improvement of the strength reduction process
!----------------------------------------------------------------------------

      implicit none


! declaration of ABAQUS specific variables:
!-----------------------------------------

      character*80 cmname
      integer ntens, ndi, nshr, nstatv, nprops, noel, npt,&
         layer, kspt, jstep(4), kinc
      double precision stress(ntens), statev(nstatv),&
         ddsdde(ntens,ntens), ddsddt(ntens), drplde(ntens),&
         stran(ntens), dstran(ntens), time(2), predef(1), dpred(1),&
         props(nprops), coords(3), drot(3,3), dfgrd0(3,3), dfgrd1(3,3),&
         pnewdt, celent,dtime, temp,dtemp,sse, spd, scd, rpl, drpldt


! declaration of internal variables:
!----------------------------------

      integer tension, isNaN, dtime_tmp
      double precision ameanstress, initstress(6),initvoid,zero, p_t,S(3),stress_pt(ntens)
      double precision phi_a,phi_p
      double precision initprops(nprops), TR(3,3)
      double precision initstatev(nstatv),props_red(nprops),isASBS, isPEAK
! initial OCR
      double precision initOCR


! integration Variables
      integer i, error, maxnint, nasv, nasvdim, &
         nfasv, nydim, nyact, j, nAddVar,naccst
      double precision D(3,3), valD, theta, &
         tolintT, tolabsT, dtsub, hmin , pmin,ddsdde_mod,pmin_voidr


! declaration of internal functions:
!-----------------------------------
      double precision valu33


! nasvdim denotes the maximum number of additional state variables
! if more than 18 are used change nasvdim and dimension of Q,
! asvr, asvrh in subroutine evolut
      parameter (nasvdim = 18)
      parameter (nydim = 42+7*nasvdim)

! additional state variables
      double precision  asv(nasvdim)
! number of first additional state variable in statev field
      parameter (nfasv = 1)
! solution vector (stresses, Jacobian information, additional state variables)
      double precision  y(nydim),y_sav(nydim)
! dummy variables for integrator
      double precision yp(nydim), v(nydim), y2(nydim), yh(nydim)
      external evolut

! AddVar: additional Variables for communication between main UMAT and subroutines
!----------------------------------------------------------------------------
      parameter (nAddVar = 11)
      double precision AddVar(nAddVar),  AddVar_tmp(nAddVar)

!     AddVar(1) ... n_Sub
!     AddVar(2) ... outputwrite (debugging variable)
!     AddVar(3) ... ISmod
!     AddVar(4) ... p_t
!     AddVar(5) ... EEM-G
!     AddVar(6) ... EEM-nu
!     AddVar(7) ... empty
!     AddVar(8) ... isCorrected
!     AddVar(9) ... empty
!     AddVar(10)... temperature -> FoR
!     AddVar(11)... p_min



! activate output for debugging:
!-------------------------------
! switch for printing information
      logical prsw, elprsw

!  parameter (prsw=.true.)
      parameter (prsw=.false.)

! print informations about time integration, useful when problems occur
      elprsw = .false.
      if (prsw) then
! print only in some defined elements and nodes
         if ((noel.eq.101).and.(npt.eq.1)) elprsw = .true.
      endif

! Error-Management:
! ----------------
! error =  0 ... no problem in time integration
! error =  1 ... problems in evaluation of the time rate, (e.g. undefined
!                stress state), reduce time integration substeps
! error =  2 ... getObjTR propagets new substepsize
!
! error =  3 ... problems in time integration, reduce abaqus load increment
!                (cut-back)
! error = 10 ... severe error in evaluating of time rate (e.g. wrong material
!                constants), terminate calculation
! tension = 1... tension stress in output. return to initial values
!
! isNaN = 1 ... NaN values in Input Data


      error = 0
      tension=0
      isNaN = 0
      isASBS = 0.0
      ddsdde_mod = 0.0

! define number of additional state variables
!---------------------------------------------
      call define(nasv)
      nyact = 42 + 7*nasv
      if (nyact.gt.nydim) then
         write(6,*) 'UMAT: nasvdim too small, program terminated'
         call XIT
      endif
      zero=0.0d0


! check material parameters in input
!---------------------------------------
      call check_param(props,nprops,error)



! calculate initial mean stress
!-----------------------
      ameanstress=-(stress(1)+stress(2)+stress(3))/3.0

! minimal Stress used for Barodesy
      pmin = 0.5
! minimal Stress used for Voidratio determination
      pmin_voidr = pmin


! ------------------------------------------------------------------------------------------------
! Initialise Variables
! ------------------------------------------------------------------------------------------------


! initialise AddVar array for internal communication
!-----------------------------------------------------
      do i = 1,nAddVar
         AddVar(i) = 0.0d0
      enddo
      AddVar(3) = statev(8) ! IS mode
      AddVar(4) = statev(11) !p_t
      AddVar(5) = statev(21) !G_EEM
      AddVar(6) = statev(22) !nu_EEM
      AddVar(7) = ddsdde_mod
      AddVar(11) = pmin



! Calculate EEM Parameters:
      if (AddVar(5).le.0.0.or.AddVar(6).le.0.0) then
         call get_EEM_params(AddVar,nAddVar,props,nprops)
      endif

      if (elprsw) AddVar(2) = 1


! remember inital values of stress,props and statev
!---------------------------------------------------
      do i=1,ntens
         initstress(i) = stress(i)
      enddo
      do i=1,nstatv
         initstatev(i) = statev(i)
      enddo
      initvoid=statev(1)

      do i=1,nprops
         initprops(i) = props(i)
      enddo
      do i=1,nAddVar
         AddVar_tmp(i) =AddVar(i)
      enddo


! initialise void ratio for initial OCR:
!---------------------------------------
      call get_voidr(statev(1),ameanstress,props,nprops,pmin_voidr)


! calculate reduced material parameters for strength reduction
!-------------------------------------------------------------
      call strength_red(temp,dtemp,props,nprops,statev,nstatv)
      AddVar(10) = statev(23)
      AddVar_tmp(10) = statev(23)


! initialise intergranular strain tensor
!----------------------------------------
      call ini_IS(statev,nstatv,props,nprops,drot,AddVar,nAddVar)


! Internal variables:
!-----------------------------

! vector of additional state variables:
      do i=1,nasv
         asv(i) = statev(i-1+nfasv)
      enddo

! calculate stretching tensor D
      call getD(D,ndi,nshr,ntens,dstran,dtime)
      valD = valu33(D)



! initialise solution vektor y and save initial values y_sav
      call iniy(y,nydim,asv,nasv,ndi,nshr,ntens,stress)
      do i = 1,nydim
         y_sav(i) = y(i)
      enddo


      if (elprsw) then
         write(6,*)
         write(6,*) '==================================================='
         write(6,*) 'Call of umat:'
         write(6,*) '==================================================='
         write(6,*)
         call wrista(3,y,nydim,D,dtime,coords,statev,nstatv,&
            props,nprops,noel,npt,ndi,nshr,jstep(1),kinc)
      endif

! Integration-Parameters:
!------------------------

! tolerance for stress error
      tolintT =  1.0d-3
! absolute tolerance for stress error
      tolabsT =  1.d-3
! suggested time substep size
      dtsub = statev(10)
      if (dtsub.le.0.0) dtsub = 1.d-4
! maximum number of time substeps
      maxnint = 10000
! minimal time substep size
      hmin = 1.d-10
! parameter of the numerical differentiation: sqrt(macheps)*||D||
! double precision
      theta = 1.0d-7 * max(valD,1.d0)
! initialise integration counter
      naccst = 0



! check for NaN in Input
!-------------------------------

      do i = 1,ntens
         if(stress(i).ne.stress(i)) then
            isNaN = 1
         endif
      enddo

      do i = 1,ntens
         if (dstran(i).ne.dstran(i) ) then
            isNaN = 1
         endif
      enddo

      do i = 1,nprops
         if (props(i).ne.props(i) ) then
            isNaN = 1
         endif
      enddo

      if (isNaN.eq.1) then
         write(6,*) 'UMAT: Found NaN in Input Data'
         call wrista(2,y,nydim,D,dtime,coords,statev,nstatv,&
            props,nprops,noel,npt,ndi,nshr,jstep,kinc)
         call xit()
      endif

! check initial stress state
!----------------------------
      call check_tension(stress,tension,ntens,0.0,S)

! elastic stiffness for zeroD:
!-----------------------------
      if (valD.eq.0.d0) then

! Calculate approximation of jacobian from barodesy
         if(tension.lt.1) then
            dtime_tmp = dtime
            call getDZjac(ddsdde,ntens,ndi,nshr, &
               stress,yp,nydim,asv,nasv,props,nprops, &
               theta,dtime_tmp,dtime_tmp,tolabsT,AddVar_tmp,nAddVar)
            if (error.gt.0) then
               call wrista(2,y,nydim,D,dtime,coords,statev,nstatv,&
                  props,nprops,noel,npt,ndi,nshr,jstep,kinc)
               call xit()
            endif

         else
! Use elastic stiffness matrix in case of initial tension states
            call elaststiff(ddsdde,ameanstress,props,nprops,ntens,tolabsT,AddVar,nAddVar)
            statev = initstatev
         endif
!       return

      else !valD.ne.0


! Correct ouf of ASBS states in case of strength reduction
         call Stress_correction(stress,statev(1),ntens,statev(27),statev(9),valD,AddVar,nAddVar,props,nprops)


! -----------------------------------------------------------------------------
! Start of Time integration
! -----------------------------------------------------------------------------



         if (tension.eq.0) then

! local extrapolation based on forward Euler, variable substeps,
! consistent Jacobian and error estimation
            if ((dtsub.le.hmin).or.(dtsub.gt.dtime)) then
               dtsub = dtime
            endif

            call eulexp(y,nyact,asv,nasv,dtsub,dtime,tolintT,tolabsT,&
               maxnint,hmin,theta,D,valD,props,nprops,error,&
               yp,v,y2,yh,elprsw,AddVar,nAddVar,naccst)


            if (error.eq.3) then
! number of sub steps larger than nmaxint or sub step size to small (h < hmin)
! ABAQUS: reduce abaqus load increment
               pnewdt = 0.25d0
! Plaxis: take solution as crude approximation, maybe near tension states
               write(6,*) 'UMAT: ask ABAQUS for reduced step size (error=3)'

               call wrista(2,y,nydim,D,dtime,coords,statev,nstatv,&
                  props,nprops,noel,npt,ndi,nshr,jstep,kinc)
               dtSub = 0

            elseif (error.eq.10) then
               write(6,*) 'error = ',error
               call wrista(2,y,nydim,D,dtime,coords,statev,nstatv,&
                  props,nprops,noel,npt,ndi,nshr,jstep,kinc)


! Dont use STOP, use abaqus-subroutine XIT to terminate the program correctly
               call XIT
            endif

         endif ! tension.eq.0


! caluclate one elastic step in case of initial tension (lowstress)
         if (tension.eq.1) then
            if (elprsw) write(6,*)'Initial lowstress: use elastic calculation'

            ddsdde_mod = 1.0
            AddVar(7) = ddsdde_mod

            call evolut(nydim,y_sav,yp,asv,nasv,theta,D,valD,props,nprops,&
               dtime,tolabsT,error,dtime_tmp,AddVar,nAddVar)
            do i = 1,nydim
               y(i) = y_sav(i) + dtime*yp(i) ! explicit euler forward integration

            enddo


         endif

! calculate elastic stiffness matrix for jacobian
!-------------------------------------------------
         if (ddsdde_mod.eq.1.0) then
            if (elprsw) write(6,*) 'UMAT: Initially in Tension, calc elastic'
            ameanstress=-(y(1)+y(2)+y(3))/3 !mean stress after integration
            call elaststiff(ddsdde,ameanstress,props,nprops,ntens,tolabsT,AddVar,nAddVar)

         endif !ddsdde_mod

! calculate output from solution vector y
!-----------------------------------------
         call solout(stress,ntens,ndi,nshr,&
            asv,nasv,ddsdde,dtime,y,nydim,ddsdde_mod)

         do i=1,nasv
            statev(i-1+nfasv) = asv(i)
         enddo


! -----------------------------------------------------------------------------
! End of time integration
! -----------------------------------------------------------------------------

      endif !case valD=0

! ------------------------------------------------------------------------------------------------
! Prepare Output
! ------------------------------------------------------------------------------------------------

! generate output stress with consideration of numerical cohesion
!-----------------------------------------------------------------
      ameanstress = -1.0/3.0*(stress(1)+stress(2)+stress(3))
      p_t = AddVar(4)
      do i = 1,ntens
         stress_pt(i) = stress(i)
         if (i.le.3) then
            stress_pt(i) = stress_pt(i)-p_t
         endif
      enddo

! check for tension in solution stresses, consider here numerical cohesion
!---------------------------------------
      tension = 0
      call check_tension(stress_pt,tension,ntens,0.0,S)
      if(tension.eq.1 .or. error.ge.1) then
         if (elprsw) write(6,*) 'UMAT: Reset to initial values, error = ',error
         if (elprsw) write(6,*) 'UMAT: Reset to initial values, tension = ',tension
         if (elprsw) then
            call wrista(2,y,nydim,D,dtime,coords,statev,nstatv,&
               props,nprops,noel,npt,ndi,nshr,jstep,kinc)
         endif

         do i = 1,ntens
            stress(i)=initstress(i)
         enddo

         do i = 1,nstatv
            statev(i)=initstatev(i)
         enddo

         do i=1,nprops
            props(i) = initprops(i)
         enddo

         return ! Return here: Do not calculate any updateted Satev

      end if

! check for NAN in solution
!---------------------------------------
      do i = 1,ntens
         if(stress(i).ne.stress(i)) then
            isNaN = 1

         endif
      enddo

      do i = 1,ntens
         do j = 1,ntens
            if(ddsdde(i,j).ne.ddsdde(i,j)) then
               isNaN = 1
            endif
         enddo
      enddo

      if (isNaN.eq.1) then
         write(6,*) 'UMAT: Found NaN in Output Data'
         call wrista(2,y,nydim,D,dtime,coords,statev,nstatv,&
            props,nprops,noel,npt,ndi,nshr,jstep,kinc)
         call xit()
      endif


! reset initial props in case of strength reduction
! -------------------------------------------------

      do i=1,nprops
         props_red(i) = props(i) ! save reduced values
         props(i) = initprops(i) ! reset to initial props
      enddo


! Correct ouf of ASBS states in case of strength reduction
!----------------------------------------------------------
      call Stress_correction(stress,statev(1),ntens,statev(27),statev(9),valD,AddVar,nAddVar,props_red,nprops)


! ------------------------------------------------------------------------------------------------
! Generate output state variables
! ------------------------------------------------------------------------------------------------

      call phimob(ndi,nshr,ntens,stress,statev(12),zero) ! equivalent matsuoka nakai friction angle
      call eps2Inv(statev(13),dstran,stran,ntens,noel,npt) !gamma_s: deviatoric strain
      call get_delta(D,valD,statev(15))  ! dilatancy
      call check_ASBS(isPeak,isASBS,statev(1),stress,ntens,props_red,nprops,&
         AddVar,nAddVar,phi_p,statev(25)) ! check for ASBS states, statev(25): phi_ASBS

! Set isPeak or isASBS condition for output
      statev(26) = 0.0
      if (isASBS.ge.1.0) then
         statev(26) = 1.0 ! ASBS Points
      endif


      call getOCR(statev(14),props,nprops,stress,ntens,statev(1)) ! overconsildation ratio
      call get_W2(D,stress,initstress,ntens,dtime,statev(17))   ! second order work
      statev(10) = min(dtsub,dtime) !propagate substepsize for next increment
      statev(23) = temp + dtemp  ! Strength reduction factor
      statev(18) = naccst ! number of substeps for actual increment
      call get_Lode(stress,statev(16),ntens) ! Lode angle
      statev(11) = AddVar(4) ! p_t
      call get_IS_out(statev,nstatv,props,nprops,D,valD,AddVar,nAddVar) ! Output from intergranular strain

      statev(21) = AddVar(5) ! G_EEM
      statev(22) = AddVar(6) ! nu_EEM


      return
   end !UMAT


! ============================================================================
! End of main subroutine UMAT
! ============================================================================

!-----------------------------------------------------------------------------
   subroutine getDZjac(ddsdde,ntens,ndi,nshr,&
      stress,yp,nydim,asv,nasv,props,nprops,&
      theta,dtime,h,tolabs,AddVar,nAddVar)
!-----------------------------------------------------------------------------
! generation of a crude approximation of the Jacobian matrix for D=0
!

      implicit none
      integer ntens, ndi, nshr, nasv, nprops, nydim, nAddVar
      double precision ddsdde(ntens,ntens), stress(ntens),&
         asv(nasv), props(nprops), theta,&
         dtime, yp(nydim),hprop, AddVar(nAddVar), &
         asvr(nasv)

      integer i, j, k, error,l,tension
      double precision T(3,3), Dh1(3,3), Dh2(3,3), TR1(3,3),&
         TR2(3,3), B(6,6), valDh1, valDh2, h, tolabs,pmin
      double precision  asv_tmp(nasv),AddVar_tmp(nAddVar), ameanstress
      hprop = 0.0
      h = 1e-10
      error = 0
      pmin = AddVar(11)

      ! Save initial Values of asv and AddVar
      do i = 1,nasv
         asv_tmp(i) = asv(i)
         asvr(i) = 0.0
      enddo
      do i = 1,nAddVar
         AddVar_tmp(i) = AddVar(i)
      enddo

106   format(1X,3(a9,f12.4,2X))
      call getT(T,ndi,nshr,ntens,stress)
!     call check_tension(stress,tension,ntens,S)
      ameanstress = -1.0/3.0*(T(1,1)+T(2,2)+T(3,3))
      do j=1,6
         do i=1,3
            do k=1,3
               Dh1(i,k) = 0.0
               Dh2(i,k) = 0.0
            enddo
         enddo
! perturbation of D --> Dh
         if (j.le.3) then
            Dh1(j,j) = Dh1(j,j) - theta
            Dh2(j,j) = Dh2(j,j) + theta
         else
            i = j/3
            k = j-i-1
            Dh1(i,k) = Dh1(i,k) - theta
            Dh2(i,k) = Dh2(i,k) + theta
         endif
         valDh1 = theta
         valDh2 = theta

! calculate perturbed objective time derivative of stresses
         call getobjtr(TR1,asvr,error,T,Dh1,valDh1,&
            asv,nasv,props,nprops,h,tolabs,hprop,AddVar,nAddVar)
         if (error.eq.1.or.error.gt.2) then
            write(6,*) 'UMAT: Error during first call with D=0.   '
            write(6,*) 'Error = ', error
            write(6,*) 'Stresses:'
            write(6,*)
            write(6,106) 'T(1,1) = ',T(1,1),'T(1,2) = ',T(1,2),'T(1,3) = ',&
               T(1,3)
            write(6,106) 'T(2,1) = ',T(2,1),'T(2,2) = ',T(2,2),'T(2,3) = ',&
               T(2,3)
            write(6,106) 'T(3,1) = ',T(3,1),'T(3,2) = ',T(3,2),'T(3,3) = ',&
               T(3,3)
            do i = 1,nprops
               write(6,*) 'props(',i,') = ',props(i)
            enddo
            do i = 1,nasv
               write(6,*) 'asv(',i,') = ',asv(i)
            enddo
            do i = 1,nAddVar
               write(6,*) 'AddVar(',i,') = ',AddVar(i)
            enddo
            return
         endif

         do i = 1,nasv
            asv(i) = asv_tmp(i)
         enddo
         do i = 1,nAddVar
            AddVar(i) = AddVar_tmp(i)
         enddo


         call getobjtr(TR2,asvr,error,T,Dh2,valDh2,&
            asv,nasv,props,nprops,h,tolabs,hprop,AddVar,nAddVar)
         if (error.eq.1.or.error.gt.2) then
            write(6,*) 'UMAT: Error during second call with D=0.   '
            write(6,*) 'Error = ', error
            write(6,*) 'Stresses:'
            write(6,*)
            write(6,106) 'T(1,1) = ',T(1,1),'T(1,2) = ',T(1,2),'T(1,3) = ',&
               T(1,3)
            write(6,106) 'T(2,1) = ',T(2,1),'T(2,2) = ',T(2,2),'T(2,3) = ',&
               T(2,3)
            write(6,106) 'T(3,1) = ',T(3,1),'T(3,2) = ',T(3,2),'T(3,3) = ',&
               T(3,3)
            do i = 1,nprops
               write(6,106) 'props(',i,') = ',props(i)
            enddo
            do i = 1,nasv
               write(6,106) 'asv(',i,') = ',asv(i)
            enddo
            do i = 1,nAddVar
               write(6,106) 'AddVar(',i,') = ',AddVar(i)
            enddo
            return
         endif

         do i = 1,nasv
            asv(i) = asv_tmp(i)
         enddo
         do i = 1,nAddVar
            AddVar(i) = AddVar_tmp(i)
         enddo
         do i=1,3
            do k=1,3
               TR2(i,k) = ( TR2(i,k) - TR1(i,k) ) / theta / 2
            enddo
         enddo
         do i=1,3
            B(i,j) = TR2(i,i)
         enddo
         B(4,j) = TR2(1,2)
         B(5,j) = TR2(1,3)
         B(6,j) = TR2(2,3)
      enddo
      do i=1,ndi
         do j=1,ndi
            ddsdde(i,j) = B(i,j)
         enddo
      enddo

      do i=ndi+1,ndi+nshr
         do j=1,ndi
            ddsdde(i,j) = B(i+3-ndi,j)
         enddo
      enddo

      do i=1,ndi
         do j=ndi+1,ndi+nshr
            ddsdde(i,j) = B(i,j+3-ndi)
         enddo
      enddo

      do i=ndi+1,ndi+nshr
         do j=ndi+1,ndi+nshr
            ddsdde(i,j) = B(i+3-ndi,j+3-ndi)
         enddo
      enddo

      return
   end


!------------------------------------------------------------------------------
   double precision function valu33(x)
!------------------------------------------------------------------------------
! value of a symmetric 3x3 matrix
!------------------------------------------------------------------------------
      implicit none
      double precision x(3,3)
      valu33 = 0.0d0
      valu33 = sqrt(x(1,1)**2.0+x(2,2)**2.0+x(3,3)**2.0+&
         2.0d0*(x(1,2)**2.0+x(1,3)**2.0+x(2,3)**2.0))
      return
   end !valu33

!-----------------------------------------------------------------------------
   subroutine iniy(y,nydim,asv,nasv,ndi,nshr,ntens,stress)
!-----------------------------------------------------------------------------
! initializes the vector of state variables
!-----------------------------------------------------------------------------
      implicit none
      integer nydim, nasv, ndi, nshr, ntens
      double precision y(nydim), asv(nasv), stress(ntens)

      integer i

      do i=1,nydim
         y(i) = 0.0
      enddo

      do i=1,ndi
         y(i) = stress(i)
      enddo
      if (nshr.ge.1) then
         y(4) = stress(ndi+1)
      endif
      if (nshr.ge.2) then
         y(5) = stress(ndi+2)
      endif
      if (nshr.ge.3) then
         y(6) = stress(ndi+3)
      endif
! additional state variables
      do i=1,nasv
         y(42+i) = asv(i)
      enddo

      return
   end !iniy


!------------------------------------------------------------------------------
   subroutine solout(stress,ntens,ndi,nshr, &
      asv,nasv,ddsdde,dtime,y,nydim,ddsdde_mod)
!------------------------------------------------------------------------------
! transform the vector of state variables to UMAT output
!------------------------------------------------------------------------------
      implicit none
      integer nydim, nasv, ndi, nshr, ntens
      double precision y(nydim), asv(nasv), stress(ntens),&
         ddsdde(ntens,ntens), dtime ,ddsdde_mod

      integer i,j
      do i = 1,ntens
         stress(i) = 0.0
      enddo
! updated stresses
      do i=1,ndi
         stress(i) = y(i)
      enddo
      if (nshr.ge.1) stress(ndi+1) = y(4)
      if (nshr.ge.2) stress(ndi+2) = y(5)
      if (nshr.ge.3) stress(ndi+3) = y(6)

! additional state variables
      do i=1,nasv
         asv(i) = y(42+i)
      enddo

! only if consistent tangent was calculated during integration process
      if (ddsdde_mod.eq.0 )then

! Jacobian
         do i=1,ndi
            do j=1,ndi
               ddsdde(i,j) = y(i+6*j)/dtime
            enddo
         enddo

         do i=ndi+1,ndi+nshr
            do j=1,ndi
               ddsdde(i,j) = y((i+3-ndi)+6*j)/dtime
            enddo
         enddo

         do i=1,ndi
            do j=ndi+1,ndi+nshr
               ddsdde(i,j) = y(i+6*(j+3-ndi))/dtime
            enddo
         enddo

         do i=ndi+1,ndi+nshr
            do j=ndi+1,ndi+nshr
               ddsdde(i,j) = y((i+3-ndi)+6*(j+3-ndi))/dtime
            enddo
         enddo

      endif

      return
   end !solout


!------------------------------------------------------------------------------
   subroutine getT(T,ndi,nshr,ntens,stress)
!------------------------------------------------------------------------------
! transform stresses from vector into matrix form
!------------------------------------------------------------------------------
      implicit none
      integer ndi, nshr, ntens
      double precision stress(ntens), T(3,3)

      integer i,j

      do i=1,3
         do j=1,3
            T(i,j) = 0.0
         enddo
      enddo

      do i=1,ndi
         T(i,i) = stress(i)
      enddo

      if (nshr.ge.1) then
         T(1,2) = stress(ndi+1)
         T(2,1) = T(1,2)
      endif
      if (nshr.ge.2) then
         T(1,3) = stress(ndi+2)
         T(3,1) = T(1,3)
      endif
      if (nshr.ge.3) then
         T(2,3) = stress(ndi+3)
         T(3,2) = T(2,3)
      endif

      return
   end !getT


!------------------------------------------------------------------------------
   subroutine getD(D,ndi,nshr,ntens,dstran,dtime)
!------------------------------------------------------------------------------
! strain increment into strain rate (from vector into matrix form)
!------------------------------------------------------------------------------
      implicit none
      integer ndi, nshr, ntens
      double precision D(3,3), dstran(ntens), dtime

      integer i,j


      D(1,1) = 0.0
      D(2,1) = 0.0
      D(3,1) = 0.0

      D(1,2) = 0.0
      D(2,2) = 0.0
      D(3,2) = 0.0

      D(1,3) = 0.0
      D(2,3) = 0.0
      D(3,3) = 0.0


      do i=1,ndi
         D(i,i) = dstran(i)/dtime
      enddo
      if (nshr.ge.1) then
         D(1,2) = 0.5d0*dstran(ndi+1)/dtime
         D(2,1) = D(1,2)
      endif
      if (nshr.ge.2) then
         D(1,3) = 0.5d0*dstran(ndi+2)/dtime
         D(3,1) = D(1,3)
      endif
      if (nshr.ge.3) then
         D(2,3) = 0.5d0*dstran(ndi+3)/dtime
         D(3,2) = D(2,3)
      endif

      return
   end !getD

!------------------------------------------------------------------------------
   subroutine evolut(n,y,yp,asv,nasv,theta,D,valD,props,nprops,&
      h,tolabs,error,hprop,AddVar,nAddVar)
!------------------------------------------------------------------------------
! evolution equations of the state variables
!------------------------------------------------------------------------------
      implicit none
      integer n, nprops, nasv,ntens, nAddVar
      double precision y(n), yp(n), D(3,3), valD
      double precision props(nprops), theta, asv(nasv),void, AddVar(nAddVar)
      integer error

      integer i,ii,j,k
      double precision T(3,3),Th(3,3), TR(3,3), Tvec(6), TRh(3,3), Dh(3,3), &
         valDh, valu33, h, &
         tolabs, hprop, AddVar_tmp(nAddVar),ddsdde_mod
!      double precision  Q(18), asvr(18), asvrh(18)
      double precision  Q(nasv), asvr(nasv), asvrh(nasv)

      error = 0
      hprop = h
      ddsdde_mod = AddVar(7)
      do i=1,3
         T(i,i)=y(i)
      enddo

      T(1,2) = y(4)
      T(1,3) = y(5)
      T(2,3) = y(6)
      T(2,1) = y(4)
      T(3,1) = y(5)
      T(3,2) = y(6)


      do i=1,nasv
         asv(i) = y(i+42)
         asvr(i) = 0.0
         asvrh(i) = 0.0
         Q(i) = 0.d0
      enddo


      call getobjtr(TR,asvr,error,T,D,valD,asv,nasv,props,nprops,&
         h,tolabs,hprop,AddVar,nAddVar)


      if (error.ge.1) then
         return
      endif


      do i=1,3
         yp(i)=TR(i,i)
      enddo
      yp(4) = TR(1,2)
      yp(5) = TR(1,3)
      yp(6) = TR(2,3)
      do i=1,nasv
         y(i+42) = asv(i)
         yp(42+i) = asvr(i)
      enddo

      if (ddsdde_mod.eq.0.0) then

! time derivative of consistent Jacobian
!----------------------------------------

         hprop = h
!
!
         do j=1,6
!! differentiation of the i'th stress y(i)
!! with respect to the j'th stretching
! B_ij = y(i+6*j) = stress(i) / (dstran(j) / dtime) etc.
! B_11 = y(7) = partial T_11 / partial D_11 etc.
            do i=1,6
               Tvec(i) = y(i)+theta*y(i+6*j)
            enddo !i
            do i=1,3
               Th(i,i) = Tvec(i)
            enddo !i
            Th(1,2) = Tvec(4)
            Th(1,3) = Tvec(5)
            Th(2,3) = Tvec(6)
            Th(2,1) = Tvec(4)
            Th(3,1) = Tvec(5)
            Th(3,2) = Tvec(6)
!  additional state variables
!  G_1 = partial e / partial D_11
            do i=1,nasv
               Q(i) = y(42+i)+theta*y(42+nasv+i+nasv*(j-1))
            enddo !i
            do i=1,3
               do k=1,3
                  Dh(i,k) = D(i,k)
               enddo !k
            enddo !i
! perturbation of D --> Dh
            if (j.le.3) then
               Dh(j,j) = Dh(j,j) + theta
            else
               i = j/3
               k = j-i-1
! pertubation in y --> symmetric pertubation in D
               Dh(i,k) = Dh(i,k) + theta/2.0d0
               Dh(k,i) = Dh(k,i) + theta/2.0d0
            endif
            valDh = valu33(Dh)

! do not change Values in AddVar
            do i = 1,nAddVar
               AddVar_tmp(i) = AddVar(i)
            enddo
            AddVar_tmp(2) = 0.0
! calculate perturbed objective time derivative of stresses


            call getobjtr(TRh,asvrh,error,Th,Dh,valDh,Q,nasv,props,nprops,&
               h,tolabs,hprop,AddVar_tmp,nAddVar)
            if (error.ge.1) then
               return
            endif
! dB = 1/theta ( TR(T+theta B, Dh, Q) - TR(T,D,e) )
            do i=1,3
               yp(i+6*j) = ( TRh(i,i) - TR(i,i) )/theta
            enddo !i
            yp(4+6*j) = ( TRh(1,2) - TR(1,2) )/theta
            yp(5+6*j) = ( TRh(1,3) - TR(1,3) )/theta
            yp(6+6*j) = ( TRh(2,3) - TR(2,3) )/theta
! dG
            do i=1,nasv
               ii = 42+nasv+i
               yp(ii+nasv*(j-1)) = ( asvrh(i)-asvr(i) ) / theta
            enddo !i

         enddo !j

      endif

      return
   end !evolut


!------------------------------------------------------------------------------
   subroutine eps2Inv(epsInv,dstran,stran,ntens,noel,npt)
!------------------------------------------------------------------------------
! calculate second Invariant of updated deviatoric strain
!------------------------------------------------------------------------------
      implicit none
      integer  ntens, i, noel,npt
      double precision epsInv
      double precision stran(ntens),dstran(ntens),temp(6)
      double precision stran_6(6),dstran_6(6)
      do i = 1,6
         temp(i) = 0.0
      enddo

      do i = 1,ntens
         temp(i)=stran(i)+dstran(i)
      end do

      epsInv = (2./9.*((temp(1)-temp(2))**2+(temp(2)-temp(3))**2+(temp(3)-temp(1))**2)+&
         1./3.*((temp(4)*temp(4))+(temp(5)*temp(5))+(temp(6)*temp(6))))**(0.5)
      return
   end !eps2Inv


!------------------------------------------------------------------------------
   subroutine phimob(ndi,nshr,ntens,stress,phim,Tc)
!------------------------------------------------------------------------------
! calculate the mobilised friction angle with the principal stresses
!------------------------------------------------------------------------------
      implicit none
      integer ndi, nshr, ntens
      double precision stress(ntens), phim, Tc

      integer lstr,i
      double precision tmin, tmax, ps(3), PI, zero
      double precision P,Q, stress_6(6)
      parameter(PI=3.141592653589793d0, zero=1d-10)
      double precision k_s, k_m,sin_phi_m, S(3), xN1(3),xN2(3),xN3(3)
      lstr = 1


! eigenvalues of the stress tensor as principal stress
!-----------------------------------------------------
      do i = 1,6
         stress_6(i) = 0.0
      enddo
      do i = 1,ntens
         stress_6(i) = stress(i)
      enddo

! Calculate mobilised friction angle from Barodesy  R-Function
      call Eig_3(1,stress_6,xN1,xN2,xN3,S(1),S(2),S(3),P,Q)
! deviatoric stress factor
      k_s = (log(-S(1)*((-S(1)*S(2)*S(3))**(-1.0/3.0))))**2.0 + (log(-S(2)/((-S(1)*S(2)*S(3))**(1.0/3.0))))**2.0 + (log(-S(3)/((-S(1)*S(2)*S(3))**(1.0/3.0))))**2.0

! Mobilised friction angle for Barodesy
      k_m = exp(sqrt(3.0/2.0*k_s))
      if (1+k_m.le.zero) then
         phim = 0.d0
      else
         sin_phi_m = (1-k_m)/(1+k_m)
         phim = abs(asin(sin_phi_m))*180/PI
      endif

      return
   end !phimob


!------------------------------------------------------------------------------
   subroutine getOCR(OCR,props,nprops,stress,ntens,voidr)
!------------------------------------------------------------------------------
! calculate actual OCR with Hvorslev equivalent pressure p_e as Output variable
!------------------------------------------------------------------------------
      integer nprops, ntens
      double precision OCR, props(nprops),stress(ntens),p,p_e,voidr

      p = (stress(1)+stress(2)+stress(3))
      p = -p/3.0
! Hvorslev equivalent pressure
!-----------------------------
      p_e = exp((props(2)- log(1+voidr))/props(3))

      if (p.gt.0) then
         OCR = p_e/p
      else
         OCR = 0.0
      end if
      if (OCR.ne.OCR) OCR = 0.0

   end ! getOCR

!------------------------------------------------------------------------------
   subroutine get_delta(D,valD,delta)
!------------------------------------------------------------------------------
! Calculate dilatancy delta = trD/valD
!------------------------------------------------------------------------------
      double precision D(3,3), valD, delta, trD
      trD = D(1,1) + D(2,2) + D(3,3)
      delta = 0.0
      if (valD.ne.0.0) then
         delta = trD/valD
      endif
      return
   end !get_delta

!------------------------------------------------------------------------------
   subroutine get_W2(D,stress1,stress0,ntens,dt,W2)
!------------------------------------------------------------------------------
! calculate second order work W2=tr(D * TR)
!
! TR(3,3) ... Stressrate over inkrementsize dt: TR =(T1-T0)/dt
! D(3,3)  ... stretching Tensor
! M.Bode 2018
!------------------------------------------------------------------------------
      integer i,ntens
      double precision W2 ,dt
      double precision D(3,3),TR(3,3),DTR(3,3),stress0(ntens),stress1(ntens),timerate(6)

      W2 = 0.0
      if(dt.ne.0) then

         do i=1,6
            timerate(i) = 0.0
         enddo

         do i=1,ntens
            timerate(i) = (stress1(i)-stress0(i))/dt
         enddo

         do i=1,3
            TR(i,i) = timerate(i)
         enddo
         TR(1,2) = timerate(4)
         TR(2,1) = timerate(4)
         TR(1,3) = timerate(5)
         TR(3,1) = timerate(5)
         TR(2,3) = timerate(6)
         TR(3,2) = timerate(6)

         call aikbkj(D,TR,DTR)

         W2=DTR(1,1)+DTR(2,2)+DTR(3,3)

      endif

      return
   end !get_W2

!------------------------------------------------------------------------------
   subroutine get_Lode(stress,theta,ntens)
!------------------------------------------------------------------------------
! calculate Lode Angle theta, output in degree
! M. Bode, 2018
! Definition of Lode angle:
! Triaxial compression: 30 degree
! pure shear : 0 degree
! Triaxial extension: -30 degree
!------------------------------------------------------------------------------
      integer i,j
      double precision stress(ntens), T(3,3), S(3,3),SS(3,3),stress_6(6)
      double precision theta, detS, P ,PI,sin3theta, J2, J3,Q

      PI=3.141592653589793
      theta = 0.0
      do i = 1,3
         do j = 1,3
            T(i,j) = 0.0d0
         enddo
      enddo
      do i = 1,6
         stress_6(i) = 0.0
      enddo
      do i = 1,ntens
         stress_6(i) = stress(i)
      enddo

      call getT(T,3,3,6,stress_6)

      P = -1.0/3.0*(stress(1)+stress(2)+stress(3))

      q = sqrt(1.0/2.0*((T(1,1)-T(2,2))**2+(T(2,2)-T(3,3))**2 + (T(3,3)-T(1,1))**2+6*((T(1,2)**2+T(2,3)**2+T(3,1)**2))))

      ! for isotropic stress state: theta = 0.0
      if (abs(Q).le.10**-10) then
         theta = 0.0d0
         return
      endif

      do i = 1,3
         do j = 1,3
            S(i,j) = T(i,j)
         enddo !j
      enddo !i

      do i = 1,3
         S(i,i) = S(i,i) + p
      enddo !i

      call aikbkj(S,S,SS)
      J2 = (SS(1,1)+SS(2,2)+SS(3,3))/2
      if (J2.gt.10**-7) then
         call det3(S,detS)
         J3 = detS
         sin3theta = J3*sqrt((3/J2)**(3))/2

         !sometimes numerical errors appear with acos=1, becaus its not defined greater than 1
         if (sin3theta.ge.1-1e-5.and.sin3theta.le.1+1e-5) sin3theta = 1-1e-8
         if (sin3theta.le.-1+1e-5.and.sin3theta.ge.-1-1e-5) sin3theta = -1+1e-8
!	in hydrostatic state lode angle is not defined, so it's set zero
         if (sin3theta.le.1.and.sin3theta.ge.-1) then
            theta = -asin(sin3theta)*180/PI/3.0

         endif


         if (theta.ne.theta) theta = 0.0
      endif


   end !get_Lode

!------------------------------------------------------------------------------
   subroutine get_IS_out(statev,nstatev,props,nprops,D,valD,AddVar,nAddVar)
!------------------------------------------------------------------------------
! calculate Intergranular Strain output variables
! rho = magnitude of IS
! omega = angle between strainrate direction and IS direction
! M.Bode 2019
!------------------------------------------------------------------------------
      integer i,j, nstatev,nprops,nAddVar
      double precision statev(nstatev),props(nprops),AddVar(nAddVar)
      double precision rho, valdelta, delta(3,3), D(3,3), D0(3,3)
      double precision delta0_D0, valD,delta0(3,3),zeros(3,3), isMod

      zeros = reshape( (/ 0,0,0 , 0,0,0 , 0,0,0 /), (/3,3/))
      R_ID = props(8)
      m_R = props(6)
      isMod =AddVar(3)
      rho = 0.0
      delta0_D0 = 0.0
      valdelta = 0.0
      delta = zeros


      do i = 1,3
         delta(i,i) = statev(i+1)
      enddo
      delta(1,2) = 0.5*statev(5)
      delta(2,1) = 0.5*statev(5)
      delta(1,3) = 0.5*statev(6)
      delta(3,1) = 0.5*statev(6)
      delta(2,3) = 0.5*statev(7)
      delta(3,2) = 0.5*statev(7)


      if (isMod.ge.1.0) then
         valdelta = valu33(delta)
         if (valD.gt.0) then
!          D0 = D/valD
            call normalize33(D,valD,D0)
         else
            D0 = zeros
         endif
!          valdelta = valu33(delta)
         valdelta = sqrt(delta(1,1)**2+delta(2,2)**2+delta(3,3)**2+&
            2.0d0*(delta(1,2)**2+delta(1,3)**2+delta(2,3)**2))
         if (valdelta.gt.0) then



            if (valdelta.gt.R_ID) then
               rho = 1.0*valdelta/R_ID


               delta = delta/valdelta*R_ID
               valdelta = R_ID
               do i = 1,3
                  statev(i+1) = delta(i,i)
               enddo
               statev(5) = 2.0*delta(1,2)
               statev(6) = 2.0*delta(1,3)
               statev(7) = 2.0*delta(2,3)

               call normalize33(delta,valdelta,delta0)
               call tensprod2dd2(delta0,D0,delta0_D0)
               rho = 1.0*valdelta/R_ID

            endif !valdelta>R_ID

            call normalize33(delta,valdelta,delta0)
            call tensprod2dd2(delta0,D0,delta0_D0)
            rho = 1.0*valdelta/R_ID

         else
            rho = 0.0
            delta0_D0 = 0.0
         endif
      else !no intergranular strain
         rho = 0.0
         delta0_D0 = 0.0
      endif
      statev(19) = rho
      statev(20) = delta0_D0
   end !get_rho


!------------------------------------------------------------------------------
   subroutine eulexp(y,n,asv,nasv,h,time,tol,tolabs,maxnint,hmin,&
      theta,D,valD,props,nprops,error,yp,v,y2,yh,elprsw,&
      AddVar,nAddVar,naccst)
!------------------------------------------------------------------------------
!  numerical solution of y'=f(y)
!  forward Euler with local extrapolation
! Fellin and Ostermann
!------------------------------------------------------------------------------

      implicit none
      integer n, nasv, nprops, maxnint, error,maxnintsav,nAddVar

      double precision y(n), h, D(3,3), valD, asv(nasv),y1sav(n),&
         yp(n), v(n), y2(n), yh(n),y3(n),y2p(n),&
         props(nprops), theta, time, tol,&
         tolabs , AddVar(nAddVar),ysav(n)
      logical elprsw

      integer i, naccst, nrejst
      double precision sci, errt, hh, actt, fhnew, hsav, hmin,hprop,hminsav
      logical final

      if (h.eq.time) then
         final = .true.
      else
         final = .false.
      endif
      naccst = 0
      fhnew = 0.0d0
      nrejst = 0
      error = 0
      actt = 0.0d0
      hsav = h
      hprop = h
      hminsav = hmin
      maxnintsav =maxnint
      do i = 1,n
         ysav(i) = y(i)
         y1sav(i) = y(i)
         y2(i) = 0.d0
         y3(i) = 0.0d0
         yp(i) = 0.d0
         y2p(i) = 0.d0
         v(i) = 0.d0
      enddo


!  two Euler steps with size h/2

      if (elprsw) write(6,*) 'time integration starts with h = ',h
      if (elprsw) write(6,*) 'end time', time


! reset substep


!	first call of time rate
      call evolut(n,y,yp,asv,nasv,theta,D,valD,props,nprops,h,&
         tolabs,error,hprop,AddVar,nAddVar)
      if (error.ge.3) then
         if (elprsw) then
            write(6,*) 'eulexp: call evolut for h/2-step --> error > 3'
         endif

         return

      elseif (error.eq.1) then
         ! undefined stress state in evolut
         write(6,*) 'UMAT - first call of constitutive law:'
         write(6,*) 'state variables (stress state) given by ABAQUS'
         write(6,*) 'can not be handled by constitutive law!'
         error = 10
         return


      endif
! start new substep
20    continue

      if ((h.lt.hmin).and.(.not. final)) then
         if (elprsw) write(6,*) 'time sub step too small, reject step'
         error = 3
         return
      endif
      hh = 0.5d0*h

      do i=1,n
         y2(i) = y(i) + hh*yp(i)
         y1sav(i) = y(i)
      enddo

! second call of timerate
      call evolut(n,y2,y2p,asv,nasv,theta,D,valD,props,nprops,h,tolabs,&
         error,hprop,AddVar,nAddVar)


      if (error.eq.1) then
         ! undefined stress state in evolut (tr T <= 0), h = h/2
         h = h/2.0
         if (elprsw) then
            write(6,*) 'eulexp: stress rate undefined at half sub step'
            write(6,*) 'stress(i): ',(y2(i),i=1,6)
            write(6,*) 'reduce sub step size'
         endif
         do i = 1,n
            y(i) = y1sav(i)
         enddo
         goto 20
      elseif (error.eq.2) then
         ! stress rate exceeds linear area, h = hprop
         h = hprop
         ! allow more and smaller substeps in case of error 2
         hmin = 1.d-16
         maxnint = maxnintsav*3.0
         if (elprsw) then
            write(6,*) 'eulexp: stepsize to big for material Model'
            write(6,*) 'stress(i): ',(y2(i),i=1,6)
            write(6,*) 'reduce sub step size to ', hprop
         endif
         do i = 1,n
            y(i) = y1sav(i)
         enddo
         goto 20
      elseif (error.ge.3) then
         if (elprsw) then
            write(6,*) 'eulexp: call evolut at halft substep --> error > 3'
         endif
         return
      endif

!

      do i=1,n
         y3(i) =y2(i) + hh*y2p(i)
      enddo
!
!  difference of step with size h to y2
!
      do i=1,n
         v(i) = y3(i) - y(i) - h*yp(i)
      enddo

!  error estimate (of the first 6 components, stresses)
!
      errt = 0.0d0
      do i=1,6
         sci = max( abs(y3(i)), abs(y(i)) ) + tolabs
         errt = errt + ( v(i)/sci )**2.0d0
      enddo
      errt = sqrt(errt)
      errt = max(errt,1.d-16)
!
! calculate factor of step change (times a safety factor)
!
      fhnew = 0.90d0*sqrt(tol/errt)

      if (errt.gt.tol) then
!
! do not accept step size, reject step
!
         final = .false.
         nrejst = nrejst + 1

         if (elprsw) write(6,*) 'reject: h, actt', h,actt

         h = max(0.2d0, fhnew)*h
         do i = 1,n
            y(i) = y1sav(i)
         enddo
         goto 20
      else ! errt.le.tol
!
! update and try to calculate time rate at the end of the substep
! second oder update
!
         do i=1,n
            yh(i) = y3(i) + v(i)
            v(i) = 0.0d0
         enddo
         call evolut(n,yh,v,asv,nasv,theta,D,valD,props,nprops,&
            h,tolabs,error,hprop,AddVar,nAddVar)
         if (error.ge.3) then
            if (elprsw) then
               write(6,*) 'eulexp: test call evolut end of step --> error > 3'
            endif
            return
         elseif (error.eq.1) then
! updated state variables cannot be handled by constitutive law:
! reduce step size, reject step
            h = h/2
            if (elprsw) then
               write(6,*)'eulexp: stress rate undefined at end of sub step'
               write(6,*)'stress(i): ',(yh(i),i=1,6)
               write(6,*)'reduce sub step size'
            endif
            do i = 1,n
               y(i) = y1sav(i)
            enddo
            goto 20
         elseif (error.eq.2) then
            ! stress rate exceeds linear area, h = hprop
            h = hprop
            hmin = 1.d-16
            maxnint = maxnintsav*3.0
            if (elprsw) then
               write(6,*) 'eulexp: stepsize to big for material Model'
               write(6,*) 'stress(i): ',(yh(i),i=1,6)
               write(6,*) 'reduce sub step size to ', hprop
            endif
            do i = 1,n
               y(i) = y1sav(i)
            enddo
            goto 20
         else ! error.eq.0
! updated state variables can be handled by constitutive law: accept
! step size and update state variables, reuse the above objective time
! rate of state variables for new step
            actt = actt + h
            naccst = naccst + 1
            if (naccst.gt.maxnint) then
               if (elprsw) write(6,*) 'number of time substeps ',naccst,&
                  ' too big, reject step'
               error = 3
               return
            endif
            if (elprsw) write(6,*) 'accept: h, actt', h,actt


            do i=1,n
               y(i) = yh(i)
               ysav(i) = yh(i)
               yp(i) = v(i)
               v(i) = 0.0d0
            enddo



         endif !errt.gt.tol


         if (final) then
            h = max(hsav,hminsav)


            if (elprsw) then
               write(6,*) 'end of time integration'
               write(6,*) 'number of accepted substeps: ', naccst
               write(6,*) 'number of rejected substeps: ', nrejst
               write(6,*) 'time integration proposes ',&
                  'new start h = ',h
            endif

            return
         endif
! suggested new step size limited by a factor of 5
         h = min(5.d0,fhnew)*h
         if (actt+h.ge.time) then
            final = .true.
! substep size for next time step
! last suggested substep, reduced by a safety factor
            hsav = 0.5d0*h
            h = time - actt
         endif
         goto 20
      endif

      return
   end !eulexp


!------------------------------------------------------------------------------
   subroutine wrista(mode,y,nydim,D,dtime,coords,statev,nstatv,props,&
      nprops,noel,npt,ndi,nshr,jstep,kinc)
!------------------------------------------------------------------------------
! Generates Output for debugging
!------------------------------------------------------------------------------
      implicit none
      integer mode, nydim, nstatv, nprops, noel, &
         npt, ndi, nshr, jstep, kinc
      double precision y(nydim), D(3,3), coords(3), statev(nstatv),&
         props(nprops), dtime

      integer i


      if (mode.eq.2) then
         write(6,*) '==================================================='
         write(6,*) 'ERROR: abaqus job failed during call of umat'
         write(6,*) '==================================================='
         write(6,*) 'state dumb:'
         write(6,*)
      endif
      write(6,111) 'Step: ',jstep, 'increment: ',kinc,&
         'element: ', noel, 'Integration point: ',npt
      if (mode.eq.2) then
         write(6,*) 'Co-ordinates of material point:'
         write(6,104) 'x1 = ',coords(1),' x2 = ',coords(2),' x3 = ',&
            coords(3)
         write(6,*)
         write(6,*) 'Material parameters:'
         write(6,*)
         do i=1,nprops
            write(6,105) 'prop(',i,') = ',props(i)
         enddo
         write(6,*)
         write(6,102) 'No. of mean components:  ',ndi
         write(6,102) 'No. of shear components: ',nshr
         write(6,*)
      endif
      if ((mode.eq.2).or.(mode.eq.3)) then
         write(6,*) 'Stresses:'
         write(6,*)
         write(6,106) 'T(1,1) = ',y(1),'T(1,2) = ',y(4),'T(1,3) = ',&
            y(5)
         write(6,106) 'T(2,1) = ',y(4),'T(2,2) = ',y(2),'T(2,3) = ',&
            y(6)
         write(6,106) 'T(3,1) = ',y(5),'T(3,2) = ',y(6),'T(3,3) = ',&
            y(3)
         write(6,*)
         write(6,*) 'Stretching rate:'
         write(6,*)
         write(6,101) 'D(1,1) = ',D(1,1),'D(1,2) = ',D(1,2),'D(1,3) = ',&
            D(1,3)
         write(6,101) 'D(2,1) = ',D(2,1),'D(2,2) = ',D(2,2),'D(2,3) = ',&
            D(2,3)
         write(6,101) 'D(3,1) = ',D(3,1),'D(3,2) = ',D(3,2),'D(3,3) = ',&
            D(3,3)
         write(6,*)
         write(6,*) 'Time increment:'
         write(6,*)
         write(6,108) 'dtime = ',dtime
         write(6,*)
         write(6,*) 'Void ratio:'
         write(6,*)
         write(6,109) 'e = ',y(43)
         write(6,*)
         if (props(6).ge.1) then
            write(6,*) 'Intergranular strain:'
            write(6,*)
            write(6,101) 'd(1,1) = ',y(44),'d(1,2) = ',1.0/2.0*y(47),'d(1,3) = ',&
               1.0/2.0*y(48)
            write(6,101) 'd(2,1) = ',1.0/2.0*y(47),'d(2,2) = ',y(45),'d(2,3) = ',&
               1.0/2.0*y(49)
            write(6,101) 'd(3,1) = ',1.0/2.0*y(48),'d(3,2) = ',1.0/2.0*y(49),'d(3,3) = ',&
               y(46)
            write(6,*)
         endif
      endif

101   format(1X,3(a9,e11.4,2X))
102   format(1X,a25,i1)
103   format(1X,a7,i5)
104   format(1X,3(a5,f10.4,2X))
105   format(1X,a5,i2,a4,f20.3)
106   format(1X,3(a9,f12.4,2X))
107   format(1X,3(a10,f12.4,2X))
108   format(1X,a8,f16.8)
109   format(1X,a4,f10.4)
110   format(1X,a5,f10.4)
111   format(1X,a6,i4,2X,a11,i4,2X,a9,i10,2X,a19,i4)

      return
   end ! wrista


!------------------------------------------------------------------------------
   subroutine normalize33(M,valM,M0)
!------------------------------------------------------------------------------
! normalize a 3*3 matrix
!------------------------------------------------------------------------------
      implicit none
      integer i, j
      double precision M(3,3), valM(1), M0(3,3),zero
      zero=0.0

      if(abs(valM(1)).ge.zero) then
         do i=1,3
            do j=1,3
               M0(i,j) = M(i,j)/valM(1)
            enddo
         enddo
      endif
      if(abs(valM(1)).le.zero) then
         do i=1,3
            do j=1,3
               M0(i,j) = 0.0d0
            enddo
         enddo
      endif

      return
   end  !normalize33


!-----------------------------------------------------------------------------
   subroutine check_tension(stress,tension,ntens,tol,S)
!-----------------------------------------------------------------------------
      implicit none
!
      integer tension,testnan,i,ntens
!
      double precision stress(ntens),pmean
      double precision S(3),P,Q,tmin
      double precision minstress ,tol,stress_6(6)
      double precision xN1(3),xN2(3),xN3(3)
      minstress=tol
      !    minstress=0.0
      tension = 0

      pmean=-(stress(1)+stress(2)+stress(3))/3.0d0
      S(1) = 0.0
      S(2) = 0.0
      S(3) = 0.0
!       check for positive mean stress
      if(pmean .lt. minstress) then
         tension=1
         return
      end if
!

      do i = 1,6
         stress_6(i) = 0.0
      enddo
      do i = 1,ntens
         stress_6(i) = stress(i)
      enddo
! Eigenvalues of stress tensor as principal stresses
!----------------------------------------------------

      call Eig_3(1,stress_6,xN1,xN2,xN3,S(1),S(2),S(3),P,Q)

! check for tension in principal stresses
!-----------------------------------

      do i=1,3
         if(S(i) .gt. -minstress) then
            tension=1
         endif
      enddo

      return
   end !check tension

!***********************************************************************
   Subroutine Eig_3(iOpt,St,xN1,xN2,xN3,S1,S2,S3,P,Q)
      Implicit Double Precision (A-H,O-Z)
      Dimension St(6),A(3,3),V(3,3),&
         xN1(3),xN2(3),xN3(3)
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
      V(1,1) = 1.0
      V(2,1) = 0.0
      V(3,1) = 0.0

      V(1,2) = 0.0
      V(2,2) = 1.0
      V(3,2) = 0.0

      V(1,3) = 0.0
      V(2,3) = 0.0
      V(3,3) = 1.0


      abs_max_s=0.0
      Do i=1,3
         Do j=1,3
            if (abs(a(i,j)) .Gt. abs_max_s) abs_max_s=abs(a(i,j))
         End Do
      End Do
      Tol = 1d-20 * abs_max_s
      it = 0
      itmax = 150
      Do While ( it.Lt.itMax .And.&
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
            If (abs(a(ip,iq)) .gt. Tol) Then
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
   End ! Eig_3

!------------------------------------------------------------------------------
   subroutine aikbkj(a,b,c)
!------------------------------------------------------------------------------
! matrix multiplication: a*b=c
!------------------------------------------------------------------------------
      implicit none
      integer i, j, k
      double precision a(3,3), b(3,3), c(3,3)
      c = reshape( (/ 0.0,0.0,0.0 , 0.0,0.0,0.0 , 0.0,0.0,0.0 /), (/3,3/))

      do i=1,3
         do j=1,3
            do k=1,3
               c(i,j) = c(i,j) + a(i,k)*b(k,j)
            enddo
         enddo
      enddo

      return
   end ! aikbkj

!------------------------------------------------------------------------------
   subroutine det3(A,detA)
!------------------------------------------------------------------------------
! Determinant of 3x3 Matrix A
!------------------------------------------------------------------------------
      double precision A(3,3), detA, detA1,detA2

      detA1 = A(1,1)*A(2,2)*A(3,3)+A(1,2)*A(2,3)*A(3,1)+A(1,3)*A(2,1)*A(3,2)
      detA2 = A(3,1)*A(2,2)*A(1,3)+A(3,2)*A(2,3)*A(1,2)+A(3,3)*A(2,1)*A(1,2)
      detA = detA1-detA2
      return
   end ! det3


!------------------------------------------------------------------------------
   subroutine tensprod2dd2(A,B,C)
!------------------------------------------------------------------------------
! inner product for 2nd order tensors A(3,3):B(3,3) = C
!------------------------------------------------------------------------------
      integer i,j
      double precision A(3,3), B(3,3), C
      C = 0.d0
      do i=1,3
         do j=1,3
            C = C + A(i,j)*B(i,j)
         enddo
      enddo
      return


   end ! tensprod2dd2

!------------------------------------------------------------------------------
   subroutine tensprod4dd2(A,B,C)
!------------------------------------------------------------------------------
! inner product for 4th and 2nd order tensors A(3,3,3,3):B(3,3) = C
!------------------------------------------------------------------------------
      double precision A(3,3,3,3), B(3,3), C(3,3)
      integer i, j, k, l

      do i=1,3
         do j=1,3
            C(i,j) = 0.d0
            do k=1,3
               do l=1,3
                  C(i,j) = C(i,j) + A(i,j,k,l)*B(k,l)
               enddo
            enddo
         enddo
      enddo

      return

   end !tensprod4dd2


!-----------------------------------------------------------------------------
   subroutine tensprod22(A,B,C)
!-----------------------------------------------------------------------------
! tensorial product of a 2nd order tensors A and B
! C_ijkl = T_ij D_kl
!-----------------------------------------------------------------------------
      implicit none
      real(8) A(3,3), B(3,3), C(3,3,3,3)
      integer i, j, k, l

      do i=1,3
         do j=1,3
            do k=1,3
               do l=1,3
                  C(i,j,k,l) = 0.0
                  C(i,j,k,l) = A(i,j)*B(k,l)
               enddo
            enddo
         enddo
      enddo

      return
   end ! tensprod22


! =============================================================================
! Following subroutines have to be adapted when changing the constitutive law
! =============================================================================

!------------------------------------------------------------------------------
   subroutine define(nasv)
!------------------------------------------------------------------------------
! number of additional state variables
! must be less than  18 (otherwise change nasvdim in umat and
! dimension of Q, asvr, asvrh in subroutine evolut)
!------------------------------------------------------------------------------

      implicit none
      integer nasv
      nasv = 7
      return
   end ! define


!----------------------------------------------------------------------------------------------
   subroutine check_param(props,nprops,error)
!-----------------------------------------------------------------------------------------------
! checks material parameters for barodesy for clay
!-----------------------------------------------------------------------------------------------
! input variables
      integer error, nprops
      double precision props(nprops)

! material parameters
      double precision phic, eN,lambda,kappa,sigst
      double precision PI
      PI=3.141592653589793
! Get Properties
      if (nprops.ge.5) then
         phic = props(1)*PI/180
         eN   = props(2)
         lambda   = props(3)
         kappa   = props(4)
         Sigst  = props(5)
      else
         write(6,*) 'UMAT: To less material constants.'
         write(6,*) 'nprops = ', nprops
         write(6,*) 'You have to define at least nprops=5 constants:'
         write(6,*) &
            'phi_c [degree], N, Lambda, Kappa, Sigma*'
         error = 10
         return
      endif

! check user input on severe errors
! which would lead to a crash in evaluation of constitutive law
      if (phic.le.zero) then
         write(6,*) 'UMAT: phic = props(1) must be > 0!, phic = ',phic
         error = 10
      endif
      if (phic.ge.90) then
         write(6,*) 'UMAT: phic = props(1) must be < 90'
         error = 10
      endif
      if (eN.le.zero) then
         write(6,*) 'UMAT: N = props(2) must be > 0'
         error = 10
      endif
      if (lambda.le.zero) then
         write(6,*) 'UMAT: lambda = props(3) must be > 0'
         error = 10
      endif
      if (kappa.le.zero) then
         write(6,*) 'UMAT: kappa = props(4) must be > 0'
         error = 10
      endif
      if (Sigst.le.zero) then
         write(6,*) 'UMAT: Sigma*= props(5) must be > 0'
         error = 10
      endif
      if (Sigst.ne.1) then
         write(6,*) 'UMAT: Use Sigma*/= 1 only when NOT calculating in kPa'
      endif

      if (error.eq.10) then
         write(6,*) 'Program terminated.'
         call xit()
      endif

   end !check params


!------------------------------------------------------------------------------
   subroutine get_EEM_params(AddVar,nAddVar,props,nprops)
!------------------------------------------------------------------------------
! calculate parameters for external elastic model (EEM) according to Bode et al 2019
! M.Bode 2019
!------------------------------------------------------------------------------
      implicit none

      integer i,j,nprops,nAddVar
      double precision  PI, Kc, K0,props(nprops)
      double precision c1,c2,c3,c4,c5,c6,lambda,kappa, alpha, G_EEM,nu_EEM,phic,m,realzero
      double precision D_cu(3,3), valD_cu, D0_cu(3,3), R_cu(3,3),m_R,K_e,valu33
      double precision AddVar(nAddVar), eMat(3,3),eHeMat(3,3),valR_cu,R0_cu(3,3),KG

      PI=3.141592653589793
      realzero = 0.0d0


! material parameters
      phic = props(1)*PI/180
      lambda   = props(3)
      kappa   = props(4)

      Kc = (1-sin(phic))/(1+sin(phic))
      K0 = 1-sin(phic)

      c2= (  3*(Kc*(1-Kc)*K0*(1-K0))**0.5  +  3*kc*(1-K0)  )  /  (2*(Kc-K0))
      c1 = (Kc)  /  (c2**2 *(1-Kc))
      c4 = 1
      c5 = 1/Kc
      c3 = (-1*3**0.5/lambda + 3**0.5/kappa)  /  (2**(c5*lambda) + (0.002)**(c5*lambda) -2)
      c6 = (1)/(2*(-1*3**0.5/kappa/c3-(1-2**(lambda*c5))))


! calculate elastic parameters for the EEM concept
!---------------------------------------------

! bulk modulus K
!----------------
      K_e = 0.5*(1.0/kappa-2.0/(3.0**0.5)*c3*(2.0**(lambda*c5)-1.0)+1.0/lambda)

! shear modulus G:
!------------------
      D_cu(:,:) = 0.0
      D_cu(1,1) = -1.0
      D_cu(2,2) = 0.5
      D_cu(3,3) = 0.5
      valD_cu = valu33(D_cu)
      call normalize33(D_cu,valD_cu,D0_cu)


! calculate R-Funktion:
      m = 0.0
      KG = 1.0 - (1.0)/(1.0+c1*(m-c2)**2.0)
      alpha=realzero
      alpha = (log(KG))  /  ((1.5)**0.5)


      do i=1,3
         do j=1,3
            eMat(i,j) = alpha*D0_cu(i,j)
         enddo ! j
      enddo ! i

! e^Matrix
      call expmbaro(eMat,eHeMat)
      do i=1,3
         do j=1,3
            R_cu(i,j) = -1.0 * eHeMat(i,j)
         enddo ! j
      enddo ! i

      valR_cu = valu33(R_cu)
      call normalize33(R_cu,valR_cu,R0_cu)


! shear modulus
      G_EEM = c3/(2*sqrt(2.0))*(R0_cu(1,1)-R0_cu(2,2))

! poissons ratio
!----------------
      nu_EEM = (3.0*K_e-2.0*G_EEM)/(6.0*K_e+2.0*G_EEM)


! save calculated values
!------------------------

      AddVar(5) = G_EEM
      AddVar(6) = nu_EEM

      return
   end ! get_EEM_params


!------------------------------------------------------------------------------
   subroutine getObjTR(TR,asvr,error,T,D,valD,asv,nasv,props,nprops,&
      dt,tolabs,hprop,AddVar,nAddVar)
!------------------------------------------------------------------------------
! calculate objective time rate of stresses and additional state variables
! abaqus requires Green-McInnis-Naghdi stress rate
! M.Bode 2019
!------------------------------------------------------------------------------
      implicit none

! Input Variables
! ----------------
      integer error, nasv,nprops,nAddVar
      double precision T(3,3), D(3,3),TR(3,3), TR_e(3,3)
      double precision asvr(nasv),asv(nasv),props(nprops),AddVar(nAddVar)
      double precision valD,dt,tolabs,hprop

      integer isNAN, tension, TRmod,lowstress, ISmod
      double precision FoS, p_t,voidr,p_mean,  p_t_mod, trT,p_tol,trD
      double precision stress(6), tol_pt,S(3),valu33
      double precision D0(3,3)
      double precision Einh(3,3), zeros(3,3), II(3,3,3,3), I_I(3,3,3,3)

      integer i,j ,k,l

! State boundary surface
!--------------------------
      double precision isPEAK,isASBS,phi_m,isSBScorrect,phi_p,phi_a


! Intergranular Strain

      double precision delta(3,3),dt_delta(3,3), rho, valdelta, delta0(3,3)
      double precision delta0_delta0(3,3,3,3),term_beta_r(3,3,3,3)
      double precision delta0_D0
      double precision R_IS, beta_IS, m_R,m_elast

! Variables for debugging
      double precision T_new(3,3),S_new(3)

      Einh = reshape( (/ 1.0,0.0,0.0 , 0.0,1.0,0.0 , 0.0,0.0,1.0 /), (/3,3/))
      zeros = reshape( (/ 0.0,0.0,0.0 , 0.0,0.0,0.0 , 0.0,0.0,0.0 /), (/3,3/))

! Initialise variables
      p_mean = 0.0
      error = 0
      tension = 0
      lowstress = 0
      TRmod = 0
      TR = zeros
      rho = 0.0
      delta0_D0 = 0.0
      m_elast = 1.0
      p_tol = AddVar(11)
      do i = 1,3
         do j = 1,3
            do k = 1,3
               do l = 1,3
                  II(i,j,k,l) = 0.0
               enddo ! l
            enddo ! k
         enddo ! j
      enddo ! i

      do i = 1,3
         do j = 1,3
            II(i,j,i,j) = 1.0d0
         enddo ! j
      enddo ! i


! Initialise variables from asv
      voidr = asv(1)

      do i = 1,nasv
         asvr(i) = 0.0d0
      enddo

! Initialise variables from AddVar

      ISmod = INT(AddVar(3))
      FoS = AddVar(10)




! check input for NAN
!----------------------

      isNAN = 0
      do i = 1,3
         do j = 1,3
            if (T(i,j).ne.T(i,j)) then
               isNaN = 1
            endif
         enddo
      enddo
      do i = 1,nasv
         if (asv(i).ne.asv(i) ) then
            isNaN = 1
         endif
      enddo
      if (isNaN.eq.1) then
         write(6,*) 'UMAT: NaN found in getObjTR'
         error = 1
      endif

! check actual void ratio
!--------------------------
      if (voidr.le.0) then
         write(6,*) 'UMAT: illegal (negative) actual void ratio detected'
         write(6,*) 'e = ', voidr
         error = 3
         write(6,*) 'Try reduced ABAQUS step size (error = ',error,')'
         write(6,*) ''
         return
      endif

! initialise numerical cohesion
!-----------------------------------

      p_mean = -(T(1,1)+T(2,2)+T(3,3))*1.0/3.0

      p_t = AddVar(4)
      do i=1,3
         T(i,i) = T(i,i) - p_t
      enddo
      p_mean = p_mean + p_t
      trT = T(1,1)+T(2,2)+T(3,3)


! initialise stretching rate abbreviations
!-----------------------------------------
      trD = D(1,1)+D(2,2)+D(3,3)
      call normalize33(D,valD,D0)


!------------------------------------------------
! definition of special cases:
!------------------------------------------------
! TRmod:
! TRmod = 0: normal calculation
! TRmod = 1: lowstress
! TRmod = 2: Intergranular Strain


! undefined stress state (tension)
!---------------------------------
      stress(1)=T(1,1)
      stress(2)=T(2,2)
      stress(3)=T(3,3)
      stress(4)=T(2,1)
      stress(5)=T(1,3)
      stress(6)=T(2,3)
      call check_tension(stress,tension,3,0.0,S)

! check for tension
!---------------------
      if (tension.eq.1) then
         write(6,*) 'getobjTR: Tension stress in Input'
         write(6,*) 'T = ',stress
         error = 1
         return
      endif


! check for lowstress
!---------------------
      call check_tension(stress,lowstress,3,p_tol,S)
      if (lowstress.eq.1) then
         TRmod = 1
      endif


! check for intergranular strain
!-------------------------------
      if (ISmod.ge.1.and.FoS.le.1.0) then
         TRmod = 2
         R_IS = props(8)
         beta_IS = props(9)
         m_R = props(6)
! check for out-of-SBS states
!-------------------------------
         call check_ASBS(isPeak,isASBS,voidr,stress,6,props,nprops,&
            AddVar,nAddVar,phi_p,phi_a)

         if (isASBS.eq.1) then
            AddVar(8) = isASBS
         endif

         if (lowstress.eq.1) then
            TRmod = 1 ! perform elastic calculation
            m_elast = m_R ! stiffnes increase for elastic calculation
         endif

! initialise tensor of intergranular strain (delta)
         delta = zeros
         delta0 = zeros

         do i = 1,3
            delta(i,i) = asv(i+1)
         enddo

         delta(1,2) = 0.5*asv(5)
         delta(2,1) = 0.5*asv(5)
         delta(1,3) = 0.5*asv(6)
         delta(3,1) = 0.5*asv(6)
         delta(2,3) = 0.5*asv(7)
         delta(3,2) = 0.5*asv(7)

! calculate magnitude and normalized intergranular strain
         valdelta = 0.0d0
         valdelta = valu33(delta)

         rho = 0.0d0
         rho = valdelta/R_IS

         if (rho.gt.1.0) then

            rho = 1.0d0

            do i = 1,3
               do j = 1,3
                  delta(i,j) = delta(i,j)/valdelta*R_IS
               enddo !j
            enddo !i

!update asv in case of valdelta > R_ID
            do i = 1,3
               asv(i+1) = delta(i,i)
            enddo !i
            asv(5) = 2.0*delta(1,2)
            asv(6) = 2.0*delta(1,3)
            asv(7) = 2.0*delta(2,3)

            valdelta = R_IS


         endif !check rho

         if (isASBS.eq.1.0) then
            rho = 1.0d0
            do i = 1,3
               do j = 1,3
                  delta(i,j) = D0(i,j)*R_IS
               enddo !j
            enddo !i

!update asv in case of out-of ASBS
            do i = 1,3
               asv(i+1) = delta(i,i)
            enddo !i
            asv(5) = 2.0*delta(1,2)
            asv(6) = 2.0*delta(1,3)
            asv(7) = 2.0*delta(2,3)

            valdelta = R_IS


         endif !check ASBS

         call normalize33(delta,valdelta,delta0)
! Calculate directional interpolation delta0:D and d:D
         call tensprod2dd2(delta0,D0,delta0_D0)

      endif !check ISmod


      if (error.gt.0) then
         do i=1,3
            T(i,i) = T(i,i) + p_t
         enddo !i
         return
      endif


!------------------------------------------------
! Calculate Stress rate:
! Interception of TR-mode
!------------------------------------------------

      select case (TRmod)
       case (0) !basic calculation
         call getObjTR_BaroCL(TR,T, D,voidr,valD,D0, props, nprops,&
            AddVar, nAddVar)

       case (1) ! lowstress

         call getObjTR_elast(TR_e,T, D, voidr, props, nprops,&
            AddVar, nAddVar,lowstress)
         TR = TR_e * m_elast
       case (2) ! Intergranular strain

         call getObjTR_IS(TR,T, D, voidr,valD,D0, delta,rho, props, nprops, dt,hprop, error,&
            ISmod,AddVar, nAddVar,valdelta,delta0,delta0_D0)

      end select


      if (error.gt.0) then
         write(6,*)'getObjTR: Error after calculation of stress rate'
         write(6,*)'error =',error
         write(6,*)'TRmod =',TRmod
         return
      endif


! evolution equations for state variables
!------------------------------------------------

      asvr(1) = ( 1.0d0 + voidr)*trD

      if (ISmod.ge.1.0) then
! timerate for intergranular strain
         dt_delta =  zeros

         if (delta0_D0.gt.0) then

            call tensprod22(delta0,delta0,delta0_delta0)
            term_beta_r = II-delta0_delta0*rho**beta_IS
            call tensprod4dd2(term_beta_r,D,dt_delta)

         else
            dt_delta = D

         endif

! intergranular strain tensor
         do i = 1,3
            asvr(i+1) = dt_delta(i,i)
         enddo
         asvr(5) = 2.0*dt_delta(1,2)
         asvr(6) = 2.0*dt_delta(1,3)
         asvr(7) = 2.0*dt_delta(2,3)
      endif


! reset numerical cohesion
!------------------------------------------------
      do i=1,3
         T(i,i) = T(i,i) + p_t
      enddo


   end ! getObjTR


!----------------------------------------------------------------------------------------------
   subroutine getObjTR_elast(TR,T, D, voidr, props, nprops, &
      AddVar, nAddVar,lowstress)
!----------------------------------------------------------------------------------------------
! calculate elastic stress rate using the EEM approach
!
! M. Bode 2019
!-----------------------------------------------------------------------------------------------

! Input Variables
! ----------------
      integer nprops,nAddVar, lowstress
      double precision T(3,3), D(3,3),TR(3,3),D0(3,3)
      double precision props(nprops),AddVar(nAddVar)
      double precision voidr,ISmod,p

!
      double precision valu33,p_min,m_R
      double precision Einh(3,3), zeros(3,3), II(3,3,3,3), I_I(3,3,3,3)

      integer i,j ,k ,l

! EEM approach
      double precision G_e, nu_e, trD,trT
      double precision  M_e(3,3,3,3), TR_e(3,3)

! initialise variables
!----------------------
      Einh = reshape( (/ 1.0,0.0,0.0 , 0.0,1.0,0.0 , 0.0,0.0,1.0 /), (/3,3/))
      zeros = reshape( (/ 0.0,0.0,0.0 , 0.0,0.0,0.0 , 0.0,0.0,0.0 /), (/3,3/))
      ISmod = AddVar(3)
! build unity 4th order tensors
      do i = 1,3
         do j = 1,3
            do k = 1,3
               do l = 1,3
                  II(i,j,k,l) = 0.0d0
               enddo ! l
            enddo ! k
         enddo ! j
      enddo ! i
      do i = 1,3
         do j = 1,3
            II(i,j,i,j) = 1.0d0
         enddo ! j
      enddo ! i
      call tensprod22(Einh,Einh,I_I)

! calculate abbreviations
!------------------------------
      trD = D(1,1)+D(2,2)+D(3,3)
      trT = T(1,1)+T(2,2)+T(3,3)
      p = -1.0/3.0*trT

      p_min = AddVar(11)
      if (lowstress.eq.1) then
         p = p_min
      endif


! Check elastic parameters are already calculated:
!----------------------------------------------------
      G_e = AddVar(5)*p
      nu_e = AddVar(6)

! calculate elastic stiffeness matrix
!-------------------------------------

      do i= 1,3
         do j= 1,3
            do k = 1,3
               do l = 1,3
                  M_e(i,j,k,l) = 2.0*G_e*(II(i,j,k,l)+nu_e/(1.0-2.0*nu_e)*I_I(i,j,k,l))
               enddo ! l
            enddo ! k
         enddo ! j
      enddo ! i


! calculate elastic stress rate
      call tensprod4dd2(M_e,D,TR_e)

      TR = TR_e

   end !getObjTR_elast


!----------------------------------------------------------------------------------------------
   subroutine getObjTR_BaroCL(TR,T, D,voidr,valD,D0, props, nprops,&
      AddVar, nAddVar)
!----------------------------------------------------------------------------------------------
! calculate stress rate for Barodesy for Clay (Medicus G., and Fellin W., 2017)
!
! M.Bode,B.Schneider-Muntau, 2019
!----------------------------------------------------------------------------------------------

! input variables
      integer nprops,nAddVar
      double precision TR(3,3), T(3,3), D(3,3), D0(3,3), props(nprops), AddVar(nAddVar)
      double precision voidr, valD

! material parameters
      double precision c1,c2,c3,c4,c5, phic,lambda,kappa,sigst,eN,Kc,K0
! abbreviations
      double precision valT,trT,trD,trD0, PI, p_mean,ec,betta, LambdaG
      double precision tmp, KG,alpha, m,f,g, valR,h,valu33
      double precision R(3,3),R0(3,3),eMat(3,3),eHeMat(3,3), Einh(3,3),zeros(3,3),T0(3,3)


      zeros = reshape( (/ 0.0,0.0,0.0 , 0.0,0.0,0.0 , 0.0,0.0,0.0 /), (/3,3/))
      PI=3.141592653589793

! set material parameters
!------------------------------
      phic = props(1)*PI/180
      eN   = props(2)
      lambda   = props(3)
      kappa   = props(4)
      Sigst  = props(5)


! calculate abbreviations
!------------------------------
      trD = D(1,1)+D(2,2)+D(3,3)
      trT = T(1,1)+T(2,2)+T(3,3)
      valT = valu33(T)
      call normalize33(T,valT,T0)

      trD0 = (D0(1,1)+D0(2,2)+D0(3,3))
      p_mean = -1.0/3.0*trT


! calculate material constants:
!------------------------------

      Kc = (1.0-sin(phic))/(1.0+sin(phic))
      K0 = 1.0-sin(phic)

      c2= -(3.0*sqrt(2.0)+3.0)/2.0
      c1 = (Kc)  /  (c2**2.0 *(1.0-Kc))
      c4 = 1.0
      c5 = 1.0/Kc
      c3 = (-1.0*3**0.5/lambda + 3**0.5/kappa)  /  (2.0**(c5*lambda) + (0.002)**(c5*lambda) -2)
      c6 = (1.0)/(2.0*(-1*3**0.5/kappa/c3-(1-2**(lambda*c5))))


      LambdaG =  (kappa-lambda)/2/3**0.5 * trD0    +    (lambda+kappa)*0.5
      betta = -1.0/c3/LambdaG + 2**(c5*lambda)/3**0.5  -  1/(3**0.5)


      m = 0.0
      tmp=6.0-2.0*trD0**2.0
      if(tmp.gt.0.0) then
         m = (-3.0*trD0)  /  ((6.0-2.0*trD0**2.0)**0.5)
      endif
      KG = 1.0 - (1.0)/(1+c1*(m-c2)**2)

      alpha = 0.0
      tmp=1.5-trD0**2/2
      if(tmp.gt.0.0) then
         alpha = (log(KG))  /  ((1.5-trD0**2/2)**0.5)
      endif

! critical void ratio
      ec=exp(eN-lambda*log(2.0* (p_mean)/Sigst))-1


! calculate functions f & g
!----------------------------
      f= c6*betta*trD0 - 0.5
      g =(1.0-c6) * betta * trD0 + ((1+voidr)/(1+ec))**c5 - 0.5


! calculate R-Funktion:
!----------------------
      eMat = zeros
      eMat = alpha*D0

! e^Matrix
      call expmbaro(eMat,eHeMat)

      R =   -1.0*eheMat

      valR = valu33(R)
      call normalize33(R,valR,R0)


! calculate h-Funktion
!-----------------------
      h = c3*sigst*(valT/sigst)**c4

! get Timereate TR
      TR = h*(f*R0+g*T0)*valD

   end ! getObjTR_BaroCL


!----------------------------------------------------------------------------------------------
   subroutine getObjTR_IS(TR,T, D, voidr,valD,D0, delta,rho, props, nprops, dt,hprop, error,&
      ISmod,AddVar, nAddVar,valdelta,delta0,delta0_D0)
!----------------------------------------------------------------------------------------------
! Check IS timestep size
! Switch between IEM and EEM approach
! Calculate IS stress rate
! M. Bode 2019
!----------------------------------------------------------------------------------------------
! Input Variables
! ----------------
      integer error, nprops,nAddVar, ISmod
      double precision T(3,3), D(3,3),TR(3,3),delta(3,3),delta0(3,3)
      double precision props(nprops),AddVar(nAddVar)
      double precision valD,dt,tolabs,hprop,valdelta,delta0_D0,voidr


      integer i,j
! Intergranular Strain
      double precision rho,eta_1,eta_2, f_A,f_B1,f_B2,f_C
      double precision R_ID, beta_ID,maxD,maxeps, m_R,m_Tm,chi_ID
      double precision TR_m(3,3),TR_e(3,3),TR_ed(3,3)

! set material parameters
!---------------------------
      m_R = props(6)
      m_T = props(7)
      R_ID = props(8)
      beta_ID = props(9)
      chi_ID = props(10)


! control stepsize for intergranular strain calculation
!------------------------------------------------------
      maxD = (max(maxval(D),maxval(-D))*dt)
      if (delta0_D0.gt.0.0) then
         if (rho.lt.1.0) then
            maxeps = (0.2*R_ID)/(1-rho**beta_ID)
            hprop = abs( maxeps/maxD)
         else
            hprop = dt
         endif  !rho<1

      else
         maxeps = 0.2*R_ID
         hprop = maxeps/maxD
      endif ! d0:D0>0


      if (hprop.lt.dt) then
! restart substep with new full substepsize hprop
         error = 2
         return
      endif ! control stepsize


      select case (ISmod)
       case(1) !EEM approach (external elastic model)
         call getObjTR_EEM(TR_m,TR_e,TR_ed,T, D, voidr,valD,D0, props, nprops,&
            AddVar, nAddVar,delta0,delta0_D0)
         eta_1 = delta0_D0
         eta_2 = eta_1

       case(2) ! IEM approach (internal elastic model)
         call getObjTR_IEM(TR_m,TR_e,TR_ed,T, D, voidr,valD,D0, props, nprops,&
            AddVar, nAddVar,delta0,delta0_D0)
         ! interpolation function eta
         eta_1 = delta0_D0**2
         eta_2 = -eta_1
      end select !ISmod



! scalar stiffness factors for intergranular strain concept
      f_A = (rho**chi_ID*m_T+(1-rho**chi_ID)*m_R)
      f_B1 = rho**chi_ID*(1-m_T)
      f_B2 = rho**chi_ID*(m_R-m_T)
      f_C =  rho**chi_ID


! Calculate Stressrate using the intergranular strain concept
      if (delta0_D0.gt.0) then
         TR = f_A*TR_e+f_B1*TR_ed*eta_1*valD+f_C*(TR_m-TR_e)*eta_1
      else
         TR = f_A*TR_e+f_B2*TR_ed*eta_2*valD
      endif

   end !getObjTR_IS

!----------------------------------------------------------------------------------------------
   subroutine getObjTR_EEM(TR_m,TR_e,TR_ed,T, D, voidr,valD,D0, props, nprops, &
      AddVar, nAddVar,delta0,delta0_D0)

!----------------------------------------------------------------------------------------------
! calcualte stress rates for the Intergranular Strain  concept
! EEM approach with external elastic model
!
! M. Bode 2020
!-----------------------------------------------------------------------------------------------

! Input Variables
! ----------------
      integer nprops,nAddVar
      double precision T(3,3), D(3,3),delta0(3,3),D0(3,3)
      double precision props(nprops),AddVar(nAddVar)
      double precision valD,delta0_D0,voidr


      integer i,j ,k ,l


! EEM approach
      double precision TR_e(3,3),TR_ed(3,3), TR_m(3,3)


! calculate stress rate of barodesy
      call getObjTR_BaroCL(TR_m,T, D,voidr,valD,D0, props, nprops,&
         AddVar, nAddVar)

! calculate elastic stressrate with respect to D
      call getObjTR_elast(TR_e,T, D, voidr, props, nprops,&
         AddVar, nAddVar,0)

! calculate elastic stressrate with respect to delta0
      call getObjTR_elast(TR_ed,T, delta0, voidr, props, nprops,&
         AddVar, nAddVar,0)

   end !getObjTR_EEM

!----------------------------------------------------------------------------------------------
   subroutine getObjTR_IEM(TR_m,TR_e,TR_ed,T, D, voidr,valD,D0, props, nprops,&
      AddVar, nAddVar,delta0,delta0_D0)

!----------------------------------------------------------------------------------------------
! calcualte stress rates for the Intergranular Strain  concept
! IEM approach using the odd and even part of the stress rate
!
! M. Bode 2020
!-----------------------------------------------------------------------------------------------

! Input Variables
! ----------------
      integer  nprops,nAddVar
      double precision T(3,3), D(3,3),delta0(3,3),D0(3,3)
      double precision props(nprops),AddVar(nAddVar)
      double precision valD,valdelta,delta0_D0,voidr
      double precision TR_e(3,3),TR_ed(3,3), TR_m(3,3)
      integer i,j

! IEM approach
      double precision TR_pos(3,3),TR_neg(3,3),TR_el(3,3),TR_ne(3,3),TR_d_pos(3,3)
      double precision TR_d_neg(3,3),TR_d_el(3,3)
      double precision valdelta0


! even and odd part of stress rate
!---------------------------------

! in direction of D
      call getObjTR_BaroCL(TR_pos,T, D,voidr,valD,D0, props, nprops,&
         AddVar, nAddVar)
      call getObjTR_BaroCL(TR_neg,T, -1.0*D,voidr,valD,-1.0*D0, props, nprops,&
         AddVar, nAddVar)

! in direction of delta
      valdelta0 = 1.0
      call getObjTR_BaroCL(TR_d_pos,T, delta0,voidr,valdelta0,delta0, props, nprops,&
         AddVar, nAddVar)
      call getObjTR_BaroCL(TR_d_neg,T, -1.0*delta0,voidr,valdelta0,-1.0*delta0, props, nprops,&
         AddVar, nAddVar)

! calculate stress rates
      TR_e = 1.0/2.0*(TR_pos-TR_neg)
      TR_m = TR_pos
      TR_ed = 1.0/2.0*(TR_d_pos-TR_d_neg)



   end !getObjTR_IEM


!------------------------------------------------------------------------------
   SUBROUTINE expmbaro(a,eHm)
!------------------------------------------------------------------------------
! calculates a matrix exponential function
! only for a symmetric 3*3 matrix
!-----------------------------------------------------------------------------
      implicit none
      integer i,j,r,n,np
      double precision a(3,3),VT(3,3),V(3,3),eHm(3,3),&
         expD(3,3),h(3,3),hh(3,3), zeros(3,3)
      double precision d(3),e(3)
      double precision S(6),xN1(3),xN2(3),xN3(3), SNew(6),S1,S2,S3,P,Q
      zeros = reshape( (/ 0,0,0 , 0,0,0 , 0,0,0 /), (/3,3/))
      h = h-h
      hh =zeros
      h = zeros
      n = 3
      np = 3

! eigenvalues d and eigenvectors V
      S(1) = a(1,1)
      S(2) = a(2,2)
      S(3) = a(3,3)
      S(4) = a(1,2)
      S(5) = a(2,3)
      S(6) = a(1,3)



      call Eig_3(1,S,xN1,xN2,xN3,S1,S2,S3,P,Q)
      d(1) = S1
      d(2) = S2
      d(3) = S3

      do i = 1,3
         V(i,1) = xN1(i)
         V(i,2) = xN2(i)
         V(i,3) = xN3(i)
      end do


!V = a
! transposed matrix:
      do i = 1,3
         do j = 1,3
            VT(j,i) = V(i,j)
         end do
      end do

! exponential function
      expD(1,1) = exp(d(1))
      expD(2,2) = exp(d(2))
      expD(3,3) = exp(d(3))


! roll back in initial direction  expA = V*expD*VT (z*cc*zz)
!-----------------------------------------------------------
      do i=1,3
         do j=1,3
            do r=1,3
               h(i,j)=h(i,j)+V(i,r)*expD(r,j)
            end do
         end do
      end do

      do i=1,3
         do j=1,3
            do r=1,3
               hh(i,j)=hh(i,j)+h(i,r)*VT(r,j)
            end do
         end do
      end do


      eHm = hh
   END !expmbaro


!------------------------------------------------------------------------------
   subroutine elaststiff(ddsdde,ameanstress_inp,props,nprops,ntens,tolabs,AddVar,nAddVar)
!------------------------------------------------------------------------------
! calculate elastic stiffness tensor with compression modulus and shear modulus
! use EEM parameters for barodesy
!------------------------------------------------------------------------------
      implicit none

      integer i,j,ntens,nprops,nAddVar
      double precision ddsdde(ntens,ntens),II(6,6),krondelta(6),youngel,nuel
      double precision one, half, zero,props(nprops), ameanstress,tolabs,isMod
      double precision G_EEM,nu_EEM, m_R
      double precision ameanstress_inp
      double precision AddVar(nAddVar),p_min,p_t

      parameter(one=1.0d0,half=0.5d0,zero=0.d0)
      p_min = AddVar(11)
      p_t = AddVar(4)
      ameanstress = ameanstress_inp+p_t
! only nonzero stresses
      if(ameanstress.le.p_min) then
         ameanstress = p_min

      endif

! Check intergranular strain
      isMod = AddVar(3)


      m_R = 1.0
      if (isMod.ge.1.0) then
         m_R = props(6)
      endif


      G_EEM= AddVar(5)
      nu_EEM  = AddVar(6)

      nuel = nu_EEM
      youngel = 2*G_EEM*(1+nu_EEM)*ameanstress*max(m_R,1.d0)

! fourth order identity tensors in Voigt notation

      do i = 1,ntens
         do j=1,ntens
            II(i,j)=zero
         end do
      end do

      II(1,1)=one
      II(2,2)=one
      II(3,3)=one
      II(4,4)=half
      II(5,5)=half
      II(6,6)=half
!
      krondelta(1)=one
      krondelta(2)=one
      krondelta(3)=one
      krondelta(4)=zero
      krondelta(5)=zero
      krondelta(6)=zero
!
! Elastic stiffness tensor

      do i = 1,ntens
         do j=1,ntens
            ddsdde(i,j)=(youngel/(1+nuel))*(II(i,j) +&
               nuel/(1-2*nuel)*krondelta(i)*krondelta(j))
         end do
      end do

      return
   end !elaststiff

!------------------------------------------------------------------------------
   subroutine strength_red(temp,dtemp,props,nprops,statev,nstatv)
!------------------------------------------------------------------------------
! calculate reduced material parameters phi_c and N and save them in the props array
!
! save Initial Values of props before allpying strength reduction
! and reset parameters at the End of UMAT calculation
! M.Bode, B. Schneider-Muntau, 2018
!------------------------------------------------------------------------------
      integer nprops, nstatev
      double precision temp, props(nprops),statev(nstatv)
      double precision phi,phi_red,lambda, PI,eta,dtemp
      PI = 3.14159265358979
      lambda = props(3)
      phi = props(1) /180*PI
      phi_red = phi*180/PI


      if (temp.le.1.0) then
         temp = 1.0
      endif
      eta = temp+dtemp

! calculate only if reduction is needed (temp>1)
      if (temp.gt.1.) then
         phi_red = tan(phi)/eta
         phi_red = atan(phi_red) *180/PI
         props(1) = phi_red
      endif
! save reduced Parameters as Output State variable
      statev(24) = phi_red
      statev(23) = eta


   end !strength_red


!------------------------------------------------------------------------------
   subroutine get_voidr(initvoid,ameanstress,props,nprops,pmin)
!------------------------------------------------------------------------------
! initialise void ratio for initial OCR:
! M.Bode, 2017
!------------------------------------------------------------------------------
! adapted from David Masin
! if initial void ratio is greater than 100:
! calculate initial void ratio with OCR = initvoid-100
      integer nprops
      double precision initvoid, initOCR, ameanstress, props(nprops),pmin,p_temp

      p_temp = ameanstress
      if(initvoid .gt. 100.0) then
         if(ameanstress.lt.pmin) ameanstress=pmin
         initOCR=initvoid-100.0d0

         if (initOCR.lt.0.5) then
            write(6,*) 'UMAT: OCR less than 1, program terminated'
            call XIT
         end if

         initvoid=dexp(props(2)-props(3)*&
            dlog(ameanstress)&
            -props(3)*dlog(initOCR))-1

      endif
      ameanstress = p_temp
   end !get_voidr


!------------------------------------------------------------------------------
   subroutine ini_IS(statev,nstatv,props,nprops,drot,AddVar,nAddVar)
! Initialise values for intergranular strain
! M. Bode 2019
!------------------------------------------------------------------------------
! Initialise intergranular strain:
!    rotate Tensor
!    if valdelta>R_ID: scale to rho=1

      integer nprops, nstatv,i,j,p,q, nAddVar,error
      double precision statev(nstatv), props(nprops), valdelta, zeros(3,3)
      double precision delta(3,3), drot(3,3), eye(3,3),delta_rot(3,3),drot_T(3,3)
      double precision AddVar(nAddVar)
      double precision m_T,beta_r,chi,R_ID, m_R,isMod
      R_ID = props(8)
      eye = reshape( (/ 1.0,0.0,0.0 , 0.0,1.0,0.0 , 0.0,0.0,1.0 /), (/3,3/))
      zeros = reshape( (/ 0.0,0.0,0.0 , 0.0,0.0,0.0 , 0.0,0.0,0.0 /), (/3,3/))


      if (AddVar(10).gt.1.0) then
         AddVar(3) = 0.0
      endif
      isMod = AddVar(3)
      error = 0
      delta = zeros
      if (isMod.gt.0.0)then
! set material paramters
         m_R = props(6)
         m_T = props(7)
         R_ID = props(8)
         beta_r = props(9)
         chi = props(10)


! Check material parameters
         if (m_T.lt.1.0) then
            write(6,*) 'UMAT: Intergranular strain factor m_T < 1'
            write(6,*) 'm_T = ',m_T
            write(6,*) 'Programm terminated.'
            error = 10
         endif
         if (R_ID.le.0.0) then
            write(6,*) 'UMAT: Intergranular strain parameter R <= 0'
            write(6,*) 'R = ',R_ID
            write(6,*) 'Programm terminated.'
            error = 10
         endif
         if (m_R.lt.1.0) then
            write(6,*) 'UMAT: Intergranular strain parameter m_R < 1 '
            write(6,*) 'm_R = ',m_R
            write(6,*) 'Programm terminated.'
            error = 10
         endif
         if (beta_r.le.0.0) then
            write(6,*) 'UMAT: Intergranular strain parameter beta_r <= 0'
            write(6,*) 'beta_r = ',beta_ID
            write(6,*) 'Programm terminated.'
            error = 10
         endif
         if (chi.le.0.0) then
            write(6,*) 'UMAT: Intergranular strain parameter chi <= 0'
            write(6,*) 'Programm terminated.'
            error = 10
         endif


         if (error.gt.0) call xit()

         valdelta = 0.0
         delta = zeros
         do i = 1,3
            delta(i,i) = statev(i+1)
         enddo
         delta(1,2) = 0.5*statev(5)
         delta(2,1) = 0.5*statev(5)
         delta(1,3) = 0.5*statev(6)
         delta(3,1) = 0.5*statev(6)
         delta(2,3) = 0.5*statev(7)
         delta(3,2) = 0.5*statev(7)

         valdelta = sqrt(delta(1,1)**2+delta(2,2)**2+delta(3,3)**2+&
            2.0d0*(delta(1,2)**2+delta(1,3)**2+delta(2,3)**2))

         if (valdelta.gt.R_ID) then
            delta = delta/valdelta*R_ID

         endif

! Rotate tensor SDV with drot
         delta_rot = zeros

         do i = 1,3
            do j = 1,3
               do p = 1,3
                  do q = 1,3
                     delta_rot(i,j) = delta_rot(i,j) + drot(i,p)*drot(j,q)*delta(p,q)
                  enddo
               enddo
            enddo
         enddo

         delta = delta_rot

      else

         do i = 1,3
            do j = 1,3
               delta(i,j) = zeros(i,j)
            enddo
         enddo

      endif !ISmod


      do i = 1,3
         statev(i+1) = delta(i,i)
      enddo
      statev(5) = 2.0*delta(1,2)
      statev(6) = 2.0*delta(1,3)
      statev(7) = 2.0*delta(2,3)

   end !ini_IS



!------------------------------------------------------------------------------
   subroutine check_ASBS(isPeak,isASBS,voidr,stress_ini,ntens,props,nprops,&
      AddVar,nAddVar,phi_p,phi_a)
!------------------------------------------------------------------------------
! check if the actual state is a peak or softening state according to Medicus2019
! material parameters reduced due to strength reduction
! M. Bode 2019
!------------------------------------------------------------------------------
! Numerical cohesion is not considered within this function and has to be considered in the
! input stress
!------------------------------------------------------------------------------
      implicit none
      INTEGER nprops, nAddVar, i,ntens
      double precision isPeak, isASBS
      double precision trD0
      double precision Einh(3,3), zeros(3,3)
      double precision phic,eN,lambda, kappa, Sigst,p_t,phic_rad
      double precision alpha, betta, Kc, K0, KG, m, LambdaG
      double precision c1,c2,c3,c4,c5,c6
      double precision voidr,PI, zero, tolabs,realzero
      double precision props(nprops),stress_ini(ntens),stress(6),emax,tmp,theta
      double precision valu33, tol_pt,ameanstress,temp,stress_out(6)

      double precision  AddVar(nAddVar)


!      - State Boundary Surface
!----------------------------------------
      double precision p_e
      double precision trD0_pa_1, trD0_pa_2, p_p,phi_m,phi_p,sin_phi_p
      double precision phi_p_phi_m, phi_c_phi_m, p_p_p
      double precision a1,a2,a3,a4,a5, AA,BB,CC, trD0_pp_1,trD0_pp_2
      double precision phi_a, sin_phi_a, phi_a_phi_m

      Einh = reshape( (/ 1.0,0.0,0.0 , 0.0,1.0,0.0 , 0.0,0.0,1.0 /), (/3,3/))
      zeros = reshape( (/ 0.0,0.0,0.0 , 0.0,0.0,0.0 , 0.0,0.0,0.0 /), (/3,3/))

      zero=1d-10
      realzero=0.0d0
      PI=3.141592653589793


      do i = 1,6
         stress(i) = 0.0
      enddo
      do i = 1,ntens
         stress(i) = stress_ini(i)
      enddo

! initialise Variables
      isPeak = 0.0
      isASBS = 0.0

      temp = AddVar(10)


! get properties
      phic = props(1)
      phic_rad = phic*PI/180
      eN   = props(2)
      lambda   = props(3)
      kappa   = props(4)
      Sigst  = props(5)


      ameanstress = -(stress(1)+stress(2)+stress(3))*1.0/3.0

      p_t = AddVar(4)



! Calculate mobilized friction angle
      call phimob(3,3,6,stress,phi_m,0.0)



      if (ameanstress.gt.zero) then
! calculate material constants:
!------------------------------
         Kc = (1.0-sin(phic_rad))/(1.0+sin(phic_rad))
         K0 = 1.0-sin(phic_rad)

         c2= (  3.0*(Kc*(1.0-Kc)*K0*(1.0-K0))**0.5  +  3.0*kc*(1.0-K0)  )  /  (2.0*(Kc-K0))
         c1 = (Kc)  /  (c2**2 *(1-Kc))
         c4 = 1.0
         c5 = 1.0/Kc
         c3 = (-1.0*3**0.5/lambda + 3.0**0.5/kappa)  /  (2.0**(c5*lambda) + (0.002)**(c5*lambda) -2.0)
         c6 = (1.0)/(2.0*(-1.0*3.0**0.5/kappa/c3-(1.0-2.0**(lambda*c5))))


! Calculation of state boundarys for actual p/p_e = OCR
!------------------------------------------------------

         p_e = exp((eN-log(1+voidr))/(lambda))
         a1 = (2.0d0*ameanstress/p_e)**(lambda*c5);
         a2 = (-1.0)/(c3*lambda)
         a3 = -(lambda-kappa)/(2.0d0*sqrt(3.d0))*c3;
         a4 = (lambda+kappa)/2.d0*c3;
         a5 = 1.d0/sqrt(3.d0)*2.d0**(c5*lambda)-1/sqrt(3.d0);


! Calculation of PSBS state (peak states)
!----------------------------

         AA = -a3*a5;
         BB = a3+1-a4*a5-a1*a3;
         CC = a4-a1*a4;
! Solve quadratic equation for trD0
!-----------------------------------
         trD0_pp_1=(-BB+(BB**2.d0-4.d0*AA*CC)**0.5)/(2.d0*AA)
         trD0_pp_2=(-BB-(BB**2.d0-4.d0*AA*CC)**0.5)/(2.d0*AA)
!
!Seperate different cases of the solution
         trD0 = -3.0**0.5 !initialisation as isotropic compression
         if (trD0_pp_1.ge.-3.0**0.5.and.trD0_pp_1.le.0)then
            trD0 = trD0_pp_1
         elseif (trD0_pp_2.ge.-3.0**0.5.and.trD0_pp_2.le.0)then
            trD0 = trD0_pp_2
         endif


! Peak friction Angle for calculated trD0
!------------------------------------
         LambdaG =  (kappa-lambda)/2.0/3**0.5 * trD0    +    (lambda+kappa)*0.5
         betta = -1.0/c3/LambdaG + 2.0**(c5*lambda)/3**0.5  -  1.0/(3.0**0.5)

         m=realzero
         tmp=6.0-2.0*trD0**2.0
         if(tmp.gt.realzero) then
            m = (-3.0*trD0)  /  ((6.0-2.0*trD0**2.0)**0.5)
         endif
         KG = 1.0 - (1.0)/(1.0+c1*(m-c2)**2.0)


! peak friction angle for calculated trD0
         sin_phi_p = (1.0-KG)/(1+KG)
         phi_p = abs(asin(sin_phi_p))*180/PI
         phi_p = max(phic,phi_p)


! Calculation of ASBS state
!------------------------------
         AA = a2*a3-a3*a5;
         BB = a3+1-a4*a5-a1*a3+a2*a4;
         CC = a4-a1*a4;


! Solve quadratic equation for trD0
!-----------------------------------
         trD0_pa_1=(-BB+(BB**2.d0-4.d0*AA*CC)**0.5)/(2.d0*AA)
         trD0_pa_2=(-BB-(BB**2.d0-4.d0*AA*CC)**0.5)/(2.d0*AA)


!Seperate different cases of the solution
         trD0 = -3.0**0.5 !initialisation as isotropic compression
         if (trD0_pa_1.ge.-3.0**0.5.and.trD0_pa_1.le.3.0**0.5)then
            trD0 = trD0_pa_1
         elseif (trD0_pa_2.ge.-3.0**0.5.and.trD0_pa_2.le.3.0**0.5)then
            trD0 = trD0_pa_2
         endif


!ASBS friction Angle for calculated trD0
!------------------------------------
         LambdaG =  (kappa-lambda)/2.0/3**0.5 * trD0    +    (lambda+kappa)*0.5
         betta = -1.0/c3/LambdaG + 2.0**(c5*lambda)/3**0.5  -  1.0/(3.0**0.5)

         m=realzero
         tmp=6.0-2.0*trD0**2.0
         if(tmp.gt.realzero) then
            m = (-3.0*trD0)  /  ((6.0-2.0*trD0**2.0)**0.5)
         endif
         KG = 1.0 - (1.0)/(1.0+c1*(m-c2)**2.0)


! asymptotic friction angle for calculated trD0
         sin_phi_a = (1.0-KG)/(1+KG)
         phi_a = abs(asin(sin_phi_a))*180/PI



! Check isPeak und isASBS conditions
!-------------------------------------
         phi_p_phi_m = phi_p-phi_m
         phi_c_phi_m = phic-phi_m
         phi_a_phi_m = phi_a-phi_m


! Check peak states with a toleranze of 0.05 in friction angle.
         if (phi_p_phi_m.le.0.05)then
            if (phi_c_phi_m.le.0.05)then
               isPeak = 1.0

            endif
         endif
! Check ASBS states with a toleranze of 0.005 in friction angle.
         if (phi_a_phi_m.le.0.005)then
            isASBS = 1.0
         endif
      endif ! check ameanstress


      do i = 1,ntens
         stress_ini(i) = stress(i)
      enddo


      return

   end ! check_ASBS


!------------------------------------------------------------------------------
   subroutine Stress_correction(stress_ini,voidr,ntens,isSBScorrect,SRmod,valD,AddVar,nAddVar,props,nprops)
! Check for out of ASBS states and perform stress correction according to Bode et al 2020 (IACMAG)
! Consider numerical cohesion
! M. Bode 2020
!------------------------------------------------------------------------------
      integer ntens,nAddVar,i,nprops
      double precision stress_ini(ntens),isPeak,isASBS,isSBScorrec, SRmod
      double precision valD, AddVar(nAddVar),props(nprops),voidr

      double precision stress(6), temp, ameanstress,tol_pt, p_t, stress_out(6),realzero
      double precision phi_p,phi_a,phi_c,phi_soll

! Initialize Variables
      realzero=0.0d0

      do i = 1,6
         stress(i) = 0.0
      enddo
      phi_c = props(1)

! SRmod:
!---------
! 0 ... Project only Dry Side of Critical
! 1 ... no Projection
! 2 ... Project to ASBS
! 3 ... Project before and after Time Integration (Dry Side)

      if (valD.eq.realzero) then
         return
      endif

! consider numerical cohesion as isotropic stress
!------------------------------------------------
      tol_pt = AddVar(11) ! minimum pressure for barodesy
      ameanstress = -(stress(1)+stress(2)+stress(3))*1.0/3.0

      p_t = AddVar(4)

      do i=1,3
         stress_ini(i) = stress_ini(i) - p_t
      enddo

! swap iput stress (1xntens) to 1x6 vector
      do i = 1,ntens
         stress(i) = stress_ini(i)
      enddo

      ! Check for out of SBS states
      call check_ASBS(isPeak,isASBS,voidr,stress_ini,ntens,props,nprops,&
         AddVar,nAddVar,phi_p,phi_a)

      temp = AddVar(10)

! For Strength Reduction only:
      if (temp.gt.1.0) then

! Switch for use of correction
         if (SRmod.eq.0.0) then! Project to ASBS
            if (isSBScorrect.ge.1.0) then
               isSBScorrect = 2.0 ! Set Flag to 'was corrected before'
            else
               isSBScorrect = 0.0
            endif

            if (isASBS.eq.1.0) then
               phi_soll = phi_a
               call project_SBS(stress,phi_soll,stress_out,isSBScorrect)
               stress = stress_out
            endif !isASBS

         endif ! SRmod
      endif !temp.gt.1


! reset numerical cohesion
!------------------------------------------------
      do i=1,3
         stress(i) = stress(i) + p_t
      enddo
      do i = 1,ntens
         stress_ini(i) = stress(i)
      enddo

      return

   end !stress correction

!------------------------------------------------------------------------------
   subroutine project_SBS(stress,phi_red,stress_out,isSBScorrect)
!------------------------------------------------------------------------------
! Pressure equivalent Projection of Stress state to ASBS
! Scale deviatoric stress for a given ASBS friction angle
! M. Bode 2019
!------------------------------------------------------------------------------
      double precision phi_red,T_test,isSBScorrect
      double precision stress(6), T(3,3), ameanstress,fvec(3)
      double precision xN1(3),xN2(3),xN3(3),P,Q, V(3,3),h(3,3), hh(3,3),VT(3,3)

      integer n,i,j, tension


      double precision T_out(3,3),Einh(3,3),stress_out(6),deviator(3,3),S(3),x,x0
      Einh = reshape( (/ 1.0,0.0,0.0 , 0.0,1.0,0.0 , 0.0,0.0,1.0 /), (/3,3/))
      call getT(T,3,3,6,stress)

! Rotate T to principla directions
      call Eig_3(1,stress,xN1,xN2,xN3,S(1),S(2),S(3),P,Q)

!Transformation Matrix for principal directions
      do i = 1,3
         V(i,1) = xN1(i)
         V(i,2) = xN2(i)
         V(i,3) = xN3(i)
      end do
! transposed transformation matrix:
      do i = 1,3
         do j = 1,3
            VT(j,i) = V(i,j)
         end do
      end do

      p = -1.0/3.0*(S(1)+S(2)+S(3))
      x0 = 1.0

      call newton_project(x,x0,S,phi_red)


      deviator = T+p*Einh

! Check for tension in output stress
      do i = 1,3
         T_test = (x-1)*p+x*S(i)
         if (T_test.gt.0) then
            x = 1
            isSBScorrect=-1
            write(6,*) 'UMAT: Tension stress after stress correction'
         else
            isSBScorrect=1
         endif
      enddo

      T_out = (x-1)*p*Einh+x*T
      stress_out(1)=T_out(1,1)
      stress_out(2)=T_out(2,2)
      stress_out(3)=T_out(3,3)
      stress_out(4)=T_out(2,1)
      stress_out(5)=T_out(1,3)
      stress_out(6)=T_out(2,3)

   end ! project_SBS

!------------------------------------------------------------------------------
   subroutine newton_project(x,x0,S,phi_red)
! One dimensional Newton solver for stress projection with p=const to a given friction angle
! Calculate deviatoric stress for reduced friction angle
! M. Bode 2019
!------------------------------------------------------------------------------
      double precision x, x0,S(3),phi_red
      integer niter,maxniter
      double precision tol, theta,fx,fxs,dfx,x_new

      tol = 0.0001
      theta = 0.0001
      niter = 0
      maxniter = 50
      x = x0

      do i = 1,maxniter
         call project_fcn(x+theta,fxs,S,phi_red)

         dfx = (fxs-fx)/theta

         x_new = x - fx/dfx

         x = x_new
         call project_fcn(x,fx,S,phi_red)

         if (fx.le.tol) then

            return

         endif
      enddo
      ! no stress correction x = 1
      write(6,*) 'UMAT: Iteration Prozedure for stress correction failed. No correction performed.'
      x= 1.0
   end
!------------------------------------------------------------------------------
   subroutine project_fcn(x,feval,S,phi_red)
! Calculate mobilised friction angle for given deviatoric stress and correction factor
! Used for Strength reduction
! M. Bode 2019
!------------------------------------------------------------------------------
      double precision x,feval, S(3),p,S_corr(3)
      double precision  PI, zero
      integer i

! Variablen für Phi_mob
      double precision k_corr, k_m,phim,sin_phi_m,phi_red
      parameter(PI=3.141592653589793d0, zero=1d-10)

!    meanstress is 1/3 traceS (pressure negative

      p = -1.0/3.0*(S(1)+S(2)+S(3))

      do i = 1,3
         S_corr(i) = (x-1)*p+x*S(i)
      enddo


! deviatoric stress factor
      k_corr = (log(-S_corr(1)*((-S_corr(1)*S_corr(2)*S_corr(3))**(-1.0/3.0))))**2.0 + (log(-S_corr(2)/((-S_corr(1)*S_corr(2)*S_corr(3))**(1.0/3.0))))**2.0 + (log(-S_corr(3)/((-S_corr(1)*S_corr(2)*S_corr(3))**(1.0/3.0))))**2.0

! Mobilised friction angle for Barodesy
      k_m = exp(sqrt(3.0/2.0*k_corr))
      if (1+k_m.le.zero) then
         phim = 0.d0
      else
         sin_phi_m = (1-k_m)/(1+k_m)
         phim = abs(asin(sin_phi_m))*180/PI
      endif
      feval = phim-phi_red

   end !fcn

end module MOD_ESM_Bardesie
