module mod_SRMC_funcs
   use kind_precision_module, only: real_type => dp
   use mod_voigt_functions, only : TwoNormTensor, TwoNormTensor_strain
   use mod_stress_invariants, only : Get_invariants, calc_inc_driver_J3_invariant, calc_dev_stess, calc_mean_stress
   use mod_yield_function, only: calc_dF_to_dtheta
   use mod_state_params, only: update_strain_depend_param
   implicit none

contains



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
      implicit none

      real(real_type),intent(in) :: TensorA(N), TensorB(N)
      real(real_type),intent(out) :: Re
      integer, intent(in) :: N
      !***********************************************************************

      ! Local variables
      integer :: X, I

      X=N/2
      Re=0.0d0
      Do I=1,X
         Re=Re+TensorA(I)*TensorB(I)
      end Do
      Do I=X+1,N
         Re=Re+2*(TensorA(I)*TensorB(I))
      end do
   end subroutine TensorInnerProduct




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
      implicit none
      real(real_type), intent(in)  :: xMat(N, N), Vec(N)
      integer, intent(in)          :: IM, N
      real(real_type), intent(out) :: VecR(N)

      !***********************************************************************
      ! Local variables
      integer :: I, J
      real(real_type) :: X

      Do I=1,N
         X=0
         Do J=1,N
            X=X+xMat(I,J)*Vec(J)
         End Do
         VecR(I)=X
      End Do
      Return
   End Subroutine MatVec

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
      implicit none
      real(real_type), intent(in)  :: VecA(N), VecB(N)
      integer, intent(in)          :: N
      real(real_type), intent(out) :: Dp

      !***********************************************************************
      ! Local variables
      integer :: I
      Dp=0.0d0
      Do I=1,N
         Dp=Dp+VecA(I)*VecB(I)
      end do

   end subroutine DotProduct_2


end module mod_SRMC_funcs
