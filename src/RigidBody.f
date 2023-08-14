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
    !	Copyright (C) 2022  Members of the Anura3D MPM Research Community 
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
	  
	  
	  module ModRigidBody
      !**********************************************************************
      !
      !    Function:  
      !                     
      !     $Revision: 9707 $
      !     $Date: 2022-04-14 14:56:02 +0200 (do, 14 apr. 2022) $
      !
      !**********************************************************************
      use ModGlobalConstants
      use ModCounters 
      use ModMPMData
      use ModParticle
      use ModMeshInfo
	  use ModMPMDynContact

      implicit none

        logical, dimension(:), allocatable :: RigdBodyDOF ! Array storing the rigid body DOF
        logical, dimension(:), allocatable :: RigdBodyInterface ! Array storing the rigid body interface nodes

      contains

        subroutine InitialiseRigidBody()
        !**********************************************************************
        !
        !    Function: initialise rigid body data    
        !
        !**********************************************************************
        implicit none

          ! Local variables
          integer(INTEGER_TYPE) :: IError

          allocate(RigdBodyDOF(Counters%N), stat = IError)
          RigdBodyDOF = .false.
          
          allocate(RigdBodyInterface(Counters%NodTot), stat = IError)
          RigdBodyInterface = .false.

          if (.not.CalParams%RigidBody%IsRigidBody) RETURN

          call GetRigidBodyTotalMassAndExternalForces()
          call IdentifyRigidInterface()
		  call InitializeRigidBodyVelocity()

		end subroutine InitialiseRigidBody
		
		subroutine InitializeRigidBodyVelocity()
        !**********************************************************************
        !
        ! Function:  Sets the initial velocity of the Rigid Body if calc is restarted
        !
        !**********************************************************************

        implicit none

          ! Local variables
          integer(INTEGER_TYPE) :: IParticle, ParticleIndex, IAEl, IEl, NElemPart
          integer(INTEGER_TYPE) :: IPatch
          logical :: VelHasBeenFound
		  
		  if (.not.IsFollowUpPhase()) RETURN   !Only for restarted calculation		  

		  VelHasBeenFound= .false.
          !do IAEl = 1, Counters%NAEl ! loop over active elements
              do IPatch = 1, NumberOfPatches ! Loop over patches
                do IAEl = 1, nael_NURBS(IPatch)!Counters%NEl ! Loop over all elements 
                    
            IEl = ActiveElement(IAEl, IPatch)
            NElemPart = NPartEle(IEl, IPatch)        

            do IParticle = 1, NElemPart ! loop over material points of the element
              ParticleIndex = GetParticleIndex(IParticle, IEl, IPatch)
              if (EntityIDArray(ParticleIndex)== CalParams%RigidBody%RigidEntity) then
                ! material point belongs to the rigid body
                CalParams%RigidBody%Velocity=VelocityArray(ParticleIndex, :) ! sets the velocity
				VelHasBeenFound= .true.
				EXIT
              end if
            end do ! loop over material points           
			if (VelHasBeenFound) EXIT
                end do ! loop over active elements
                
                end do ! loop over patches

        end subroutine InitializeRigidBodyVelocity

        subroutine GetRigidBodyTotalMassAndExternalForces()
        !**********************************************************************
        !
        ! Function:  To calculate the forces and mass of a rigid body
        !
        !**********************************************************************

        implicit none

          ! Local variables
          integer(INTEGER_TYPE) :: IParticle, ParticleIndex, IAEl, IEl, NElemPart, INode, iDofOffset
          logical :: RigidEntityElm
          
          ! Multipatch variables 
          integer(INTEGER_TYPE) :: IPatch_Temporary
          integer(INTEGER_TYPE) :: IPatch

          CalParams%RigidBody%Mass = 0.0
          CalParams%RigidBody%TractionForce = 0.0
          CalParams%RigidBody%GravityForce = 0.0		  

          do IParticle = 1, Counters%NParticles ! loop over material points
            ParticleIndex = GetParticleIndexFromList(IParticle)
            if (EntityIDArray(ParticleIndex) == CalParams%RigidBody%RigidEntity) then ! particle belongs to the rigid body
              CalParams%RigidBody%Mass = CalParams%RigidBody%Mass + MassArray(ParticleIndex)
              CalParams%RigidBody%TractionForce = CalParams%RigidBody%TractionForce + Particles(ParticleIndex)%FExt(:,1)
              CalParams%RigidBody%GravityForce = CalParams%RigidBody%GravityForce + Particles(ParticleIndex)%FBody
            end if ! rigid entity particles
          end do ! loop over material points

          !do IAEl = 1, Counters%NAEl ! loop over active elements
              
              do IPatch = 1, NumberOfPatches ! Loop over patches
                do IAEl = 1, nael_NURBS(IPatch)!Counters%NEl ! Loop over all elements 
                    
            IEl = ActiveElement(IAEl, IPatch)
            NElemPart = NPartEle(IEl, IPatch)
            RigidEntityElm = .false.

            do IParticle = 1, NElemPart ! loop over material points of the element
              ParticleIndex = GetParticleIndex(IParticle, IEl, IPatch)
              if (EntityIDArray(ParticleIndex)== CalParams%RigidBody%RigidEntity) then ! material point belongs to the rigid body
                RigidEntityElm = .true.
              end if
            end do ! loop over material points

            if (RigidEntityElm) then
              do INode = 1,ELEMENTNODES ! loop over element nodes
                iDofOffset = Multipatch_Connecting_Local_To_Global_ControlPoints(ReducedDof(ElementConnectivities(INode,IEl,IPatch)), IPatch)
                ! All DOFs
                RigdBodyDOF(iDofOffset+1: iDofOffset + NDOFL) = .true.
              end do ! loop over element nodes
            end if

                end do ! loop over active elements
                
                end do ! patches

        end subroutine GetRigidBodyTotalMassAndExternalForces

        subroutine OverWriteRateofMomentumRigidBody()
        !**********************************************************************
        !
        ! Function:  To calculate the incremental nodal accelerations from nodal mass and rate of momentum
        !
        !**********************************************************************

        implicit none
        
          ! Local variables
          integer(INTEGER_TYPE) :: IDOF, RigidEntity, IAEl, iEl, iDofOffset,i, iNode
          real(REAL_TYPE) :: Mass
          real(REAL_TYPE), dimension(NVECTOR) :: TractionForce, GravityForce, RMom, InternalForce, Acc, Acc0

          ! Multipatch variables 
          integer(INTEGER_TYPE) :: IPatch_Temporary = 1
          integer(INTEGER_TYPE) :: IPatch
          
          
          if (.not.CalParams%ApplyContactAlgorithm) RETURN
          if (.not.CalParams%RigidBody%IsRigidBody) RETURN

          Mass = CalParams%RigidBody%Mass

          TractionForce = CalParams%RigidBody%CurrentTraction
          GravityForce = CalParams%RigidBody%CurrentGravity
          InternalForce = CalParams%RigidBody%InternalForce
          RMom = TractionForce + GravityForce - InternalForce
          Acc = CalParams%RigidBody%Acceleration

          CalParams%RigidBody%InAcceleration = RMom /Mass
		do i=1, NVECTOR
		    if (CalParams%RigidBody%Constrains(i)==1) &
				RMom (i)=0
		end do

          RigidEntity = CalParams%RigidBody%RigidEntity
          do IDOF = 1, Counters%N ! loop over degrees of freedom
            RateofMomentum(IDOF,RigidEntity) = 0.0 ! reset momentum of the rigid entity
            TotalVelocitySoil(IDOF,RigidEntity) = 0.0
          end do ! loop over degrees of freedom

          !do IAEl = 1, Counters%NAEl ! loop over active elements
              do IPatch = 1, NumberOfPatches ! Loop over patches
                do IAEl = 1, nael_NURBS(IPatch)!Counters%NEl ! Loop over all elements 
            IEl = ActiveElement(IAEl, IPatch)
            do INode = 1,ELEMENTNODES ! loop over element nodes
              iDofOffset = ReducedDof(Multipatch_Connecting_Local_To_Global_ControlPoints(ElementConnectivities(INode,IEl,IPatch), IPatch))
              do i = 1, NVECTOR
                iDof = iDofOffset+i
                if (RigdBodyDOF(iDof)) then
				  RateofMomentum(IDOF,RigidEntity) = CalParams%RigidBody%InAcceleration(i)*LumpedMassDry(IDOF,RigidEntity)
                  TotalVelocitySoil(IDOF, RigidEntity) = CalParams%RigidBody%Velocity(i)				  
                endif
              enddo
            enddo ! loop over element nodes
                enddo ! active elements
                end do ! patches

        end subroutine OverWriteRateofMomentumRigidBody

        subroutine GetInternalLoadRigidBody()
        !**********************************************************************
        !
        ! Function:  To calculate the incremental nodal accelerations from nodal mass and rate of momentum
        !
        !**********************************************************************

        implicit none

          integer(INTEGER_TYPE) :: INode, iDofOffset, iDim
          integer(INTEGER_TYPE) :: IPatch
          integer(INTEGER_TYPE) :: GlobalNodeID

         if (.not.CalParams%ApplyContactAlgorithm) RETURN
         if (.not.CalParams%RigidBody%IsRigidBody) RETURN
		 CalParams%RigidBody%InternalForce=0
         
         do IPatch = 1, NumberOfPatches ! loop over patches
          do INode = 1, NControlPoints(IPatch)!Counters%NodTot
              
              GlobalNodeID = Multipatch_Connecting_Local_To_Global_ControlPoints(INode,IPatch)
              
            if (RigdBodyInterface(INode)) then ! interface node
              iDofOffset = Multipatch_Connecting_Local_To_Global_ControlPoints(ReducedDof(GlobalNodeID), IPatch)
              do iDim = 1, NVECTOR
                ! reaction force on the rigid body
                CalParams%RigidBody%InternalForce(iDim) =  CalParams%RigidBody%InternalForce(iDim) + IntLoad(iDofOffset+iDim, SOFT_ENTITY)
              end do
            end if
          end do	
          end do ! loop over patches

		do iDim=1, NVECTOR
		    if (CalParams%RigidBody%Constrains(iDim)==1) &
				CalParams%RigidBody%InternalForce(iDim)=0
		end do
		
        end subroutine GetInternalLoadRigidBody

        
        subroutine IdentifyRigidInterface()
        !**********************************************************************
        !
        ! Function:  To identify the nodes where rigid body algorithm should be applied
        !
        !**********************************************************************
        
        implicit none
        
          ! Local variables
          character :: FilNME*1023
          integer(INTEGER_TYPE) :: NumNodes, INode, NodeID
          logical :: RigdInter (Counters%NodTot)
          real(REAL_TYPE) :: IDum			
		  if ( CalParams%ApplyContactAlgorithm ) then
			  RigdBodyInterface=InterfaceNodes
		  else
			  call GiveError('Contact is not enabled')
		  end if
		  

        end subroutine IdentifyRigidInterface

        subroutine GetRigidBodyAverageAcceleration()

        !**********************************************************************
        !
        !    Function:  To update particles total velocities and accelerations
        !
        !
        !    Note : This subroutine works only if NumbOfLayers = 1
        !
        !**********************************************************************

        implicit none
        
          ! Local variables
          integer(INTEGER_TYPE) :: RigidEntity, IDOF, INode, iDofOffset, i, iAEl, iEl, N
          real(REAL_TYPE) :: Mass
          real(REAL_TYPE), dimension(NVECTOR) :: TotLMass, TotAcc, Acc, TractionForce, GravityForce

          ! Multipatch variables
          integer(INTEGER_TYPE) :: IPatch_Temporary
          integer(INTEGER_TYPE) :: IPatch
          
          if (.not.(NFORMULATION==1)) RETURN
          if (.not.CalParams%ApplyContactAlgorithm) RETURN
          if (.not.CalParams%RigidBody%IsRigidBody) RETURN

          RigidEntity = CalParams%RigidBody%RigidEntity
          TotAcc = 0.0
		  TotLMass= 0.0


          Mass = CalParams%RigidBody%Mass
          TractionForce = CalParams%RigidBody%CurrentTraction
          GravityForce = CalParams%RigidBody%CurrentGravity
          Acc=CalParams%RigidBody%InAcceleration
		
		do i=1, NVECTOR
		    if (CalParams%RigidBody%Constrains(i)==1) &
				Acc(i)=0
		end do

          AccelerationSoil(1:Counters%N, RigidEntity) = 0.0
        CalParams%RigidBody%Acceleration = Acc
		CalParams%RigidBody%Velocity = CalParams%RigidBody%Velocity + CalParams%RigidBody%Acceleration * CalParams%TimeIncrement

          !do IAEl = 1, Counters%NAEl ! loop over active elements
              do IPatch = 1, NumberOfPatches ! Loop over patches
                do IAEl = 1, nael_NURBS(IAEl)!Counters%NEl ! Loop over all elements 
            IEl = ActiveElement(IAEl, IPatch)
            do INode = 1,ELEMENTNODES ! loop over element nodes
              iDofOffset = ReducedDof(Multipatch_Connecting_Local_To_Global_ControlPoints(ElementConnectivities(INode,IEl,IPatch_Temporary), IPatch) )
              do i = 1, NVECTOR
                iDof = iDofOffset+i
                if (RigdBodyDOF(iDof)) then
                  AccelerationSoil(iDof, RigidEntity) = Acc(i)
                endif
              enddo
            enddo ! loop over element nodes
                enddo ! loop over active elements
                end do ! loop over patches
        end subroutine GetRigidBodyAverageAcceleration


        subroutine OverWriteParticleAndNodalAccAndVeloRigidBody()
        !**********************************************************************
        !
        !    Function:  To calculate the forces and mass of a rigid body
        !
        !
        !    Note : this subroutine works if NumbOfLayers = 1
        !
        !**********************************************************************

        implicit none

          ! Local variables
          integer(INTEGER_TYPE) :: IParticle, ParticleIndex, EintityID, RigidEntity, IDOF, iDofOffset, i, iEl, iNode, iAEl

          
          ! Multipatch variables 
          integer(INTEGER_TYPE) :: IPatch_Temporary
          integer(INTEGER_TYPE) :: IPatch
          
          if (.not.(NFORMULATION==1)) RETURN
          if (.not.CalParams%ApplyContactAlgorithm) RETURN
          if (.not.CalParams%RigidBody%IsRigidBody) RETURN

          RigidEntity = CalParams%RigidBody%RigidEntity

          do IParticle = 1, Counters%NParticles
            ! Loop over particles
            ParticleIndex = GetParticleIndexFromList(IParticle)
            EintityID = EntityIDArray(ParticleIndex)
            if (EintityID == RigidEntity) then ! particle belongs to the rigid body
              VelocityArray(ParticleIndex,:) = CalParams%RigidBody%Velocity
              AccelerationArray(ParticleIndex, :) = CalParams%RigidBody%Acceleration
            end if ! rigid entity particles
          end do ! loop over particles

          TotalVelocitySoil(1:Counters%N,RigidEntity) = 0.0

          !do IAEl = 1, Counters%NAEl ! loop over active elements
              do IPatch = 1, NumberOfPatches ! Loop over patches
                do IAEl = 1, nael_NURBS(IPatch)!Counters%NEl ! Loop over all elements 
                    
            IEl = ActiveElement(IAEl, IPatch)
            do INode = 1,ELEMENTNODES ! loop over element nodes
              iDofOffset = ReducedDof(Multipatch_Connecting_Local_To_Global_ControlPoints(ElementConnectivities(INode,IEl,IPatch), IPatch))
              do i = 1, NVECTOR
                iDof = iDofOffset+i
                if (RigdBodyDOF(iDof)) then
                  TotalVelocitySoil(iDof,RigidEntity) = CalParams%RigidBody%Velocity(i)
                endif
              enddo
            enddo ! loop over element nodes
                enddo
                
                end do 

        end subroutine OverWriteParticleAndNodalAccAndVeloRigidBody

      end module ModRigidBody