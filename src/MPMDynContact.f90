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

    module ModMPMDynContact
    !**********************************************************************
    !
    !    Function:  This module contains all routines which are used for the Dynamic MPM contact model
    !
    !     $Revision: 10002 $
    !     $Date: 2023-06-19 12:44:04 +0200 (ma, 19 jun 2023) $
    !
    !**********************************************************************
    use ModCounters
    use ModReadCalculationData
    use ModMPMData
    use ModParticle
    use ModMeshInfo
    use ModRotBoundCond
    use ModFileIO
    use ModGeometryMath, only: rotate_2d_vector, calc_2d_line_normal
    use ModGlobalConstants, only: INTEGER_TYPE, REAL_TYPE
    use MOD2DCONTACT, only: read_2d_contact_data, ContactSurfaceNodes_2D
    
    implicit none

    character(len=MAX_FILENAME_LENGTH) :: FileVersion                       ! FileVersion for contact data
    character(len=64), dimension(:,:), allocatable :: ContactMaterialName   ! array including the Contact Material name
    integer(INTEGER_TYPE) :: PileEntityMaterial, NContactElm
    integer(INTEGER_TYPE) :: NInterfaceNodes
    integer(INTEGER_TYPE), dimension(:, :), allocatable :: NodeElement      ! Node-element relation - for each node a list of the elements containing this global node

    integer(INTEGER_TYPE), dimension(:, :), allocatable :: ContactSurfaceElmNodes

    real(REAL_TYPE), dimension(:,:), allocatable :: FrictionVector          ! friction coeff for multiple materials
    real(REAL_TYPE), dimension(:,:), allocatable :: InterfaceNodesAdhesion  ! Adhesion for multiple materials
    real(REAL_TYPE), dimension(:),   allocatable :: LumpedMassSys             ! Lumped mass vector of the entire system
    real(REAL_TYPE), dimension(:),   allocatable :: TotalVelocitySys          ! Total velocity vector of the entire system
    real(REAL_TYPE), dimension(:),   allocatable :: AccelerationSys           ! Acceleration vector of the entire system
    real(REAL_TYPE), dimension(:,:), allocatable :: ContactNodeNormals, ContactSurfaces ! Normals at nodes used with the contact algorithm
    
    logical, dimension(:), allocatable :: ContactSurfaceSoilElements
    logical, dimension(:), allocatable :: GlobContElement, InterfaceNodes ! GlobContElement - Mask for which elements are contact
                                                                          ! InterfaceNodes - Mask for which nodes are contact



    contains ! Routines of this module

    !*************************************************************************************
    !    SUBROUTINE: InitialiseContactData
    !
    !    DESCRIPTION:
    !>   Destroys the contact data from the precious time step and initiliases the new contact data
    !
    !>   @note : Notes
    !*************************************************************************************
    subroutine InitialiseContactData()
    implicit none

    ! Remove the contact data from memory
    call DestroyContactData()

    ! Initialise the contacts arrays
    call InitialiseContactArrays()

    end subroutine InitialiseContactData

    !*************************************************************************************
    !    SUBROUTINE: DestroyContactData
    !
    !    DESCRIPTION:
    !>   Deallocates the arrays used in the contact algorithm. These variables are used in a ton of modules
    !
    !>   @note : Notes
    !
    !>   @param[in] N/A: N/A
    !>   @param[out] N/A: N/A
    !>   @param[inout] N/A: N/A
    !
    !*************************************************************************************
    subroutine DestroyContactData()
    implicit none

    ! Local variables
    integer(INTEGER_TYPE) :: IError

    if (allocated(LumpedMassSys)) then
        deallocate(LumpedMassSys, stat = IError)
    end if

    if (allocated(TotalVelocitySys)) then
        deallocate(TotalVelocitySys, stat = IError)
    end if

    if (allocated(AccelerationSys)) then
        deallocate(AccelerationSys, stat = IError)
    end if

    if (allocated(ContactNodeNormals) ) then
        deallocate(ContactNodeNormals, stat = IError)
    end if

    if (allocated(InterfaceNodes) ) then
        deallocate(InterfaceNodes, stat = IError)
    end if

    if (allocated(InterfaceNodesAdhesion) ) then
        deallocate(InterfaceNodesAdhesion, stat = IError)
    end if

    if (allocated(NodeElement)) then
        deallocate(NodeElement, stat = IError)
    end if

    if (allocated(ContactSurfaceSoilElements)) then
        deallocate(ContactSurfaceSoilElements, stat = IError)
    end if

    end subroutine DestroyContactData

    !*************************************************************************************
    !    SUBROUTINE: ReadContactData
    !
    !    DESCRIPTION:
    !>   Reads and processes the contact data to node normals, depending on GOM file version 
    !
    !>   @note : Notes
    !
    !>   @param[in] N/A : N/A
    !>   @param[out] N/A : N/A
    !>   @param[inout] N/A : N/A
    !
    !*************************************************************************************
    subroutine ReadContactData()

        implicit none

        ! local variables
        character(len=MAX_FILENAME_LENGTH) :: FileName
        integer(INTEGER_TYPE) :: FileUnit, NContVolElem, nrows, ncols
        character(len=255) :: file_line
        integer(INTEGER_TYPE) :: ios
        integer(INTEGER_TYPE) :: StructureElement ! used for error control
        real(REAL_TYPE),     allocatable, dimension(:,:) :: contact_volumes
        character(len = 64), allocatable, dimension(:,:) :: material_names
    
        ! Check that the contact algorithm is being used if not return
        if (.not.CalParams%ApplyContactAlgorithm) RETURN

        FileName = trim(CalParams%FileNames%ProjectName)//GOM_FILE_EXTENSION ! Construct the file name (GOM file name from the project name and the GOM extension)
        FileUnit = TMP_UNIT                                                  ! Use a temporary file unit
    
        if ( FExist(trim(FileName)) ) then      ! check if GOM file exists in project folder, otherwise give error and stop execution
            call GiveMessage('Reading contact data from GOM file: ' // trim(FileName) )
        else
            call GiveError('GOM file does not exist!' // NEW_LINE('A') // 'required GOM file: ' // trim(FileName) )
        end if

        call FileOpen(FileUnit, trim(FileName)) ! Open the file

        ! determine current version of GOM file
        read(FileUnit, '(A)', iostat=ios) file_line ! NB: if no version is specified in the header of the GOM file, the default case will be chosen automatically below
        call Assert( ios == 0, "GOM file: Can't read flag from GOM file." )
        FileVersion = trim(file_line)
    
        ! Read each line of the file until the desired flag is found
        do
            read(FileUnit,*) file_line
            if ((trim(file_line)=='$$START_CONTACT_VOLUME').or.(trim(file_line)=='$$START_BODY_CONTACT_2D')) then
                read(FileUnit, *) NContVolElem  ! Number of element sides forming the contact surface
                if (NContVolElem>0) then

                    read(FileUnit,*) StructureElement ! Read first contact element

                    ! Get the material id of that first element and store the material id as the StructureMaterialId
                    CalParams%MovingMesh%StructureMaterialID = ElementMaterialID(StructureElement)
                    EXIT
                end if
            else if (trim(file_line)=='$$FINISH') then
                EXIT
            end if
        end do
        close(FileUnit)


    
        if (NDIM==2) then
            ! Read the 2d contact data
            call read_2d_contact_data(FileUnit, FileName, contact_volumes, ContactMaterialName)
            
            ! Get the contact nodes and the contact normals
            call ContactSurfaceNodes_2D(FileUnit, contact_volumes, ContactMaterialName, GlobContElement, ContactSurfaces,&
                                        ContactSurfaceElmNodes, InterfaceNodes, ContactNodeNormals)
        else
            select case (FileVersion) ! read GOM data depending on file version
            case (Anura3D_v2023)
                call ContactSurfaceNodes_v2021(FileUnit)
            case (Anura3D_v2022)
                call ContactSurfaceNodes_v2021(FileUnit)
            case (Anura3D_v2021)
                call ContactSurfaceNodes_v2021(FileUnit)
            end select
        end if

        ! close GOM file
        close(FileUnit)

    end subroutine ReadContactData

    

    subroutine ContactSurfaceNodes_v2021(FileUnit)
        !**********************************************************************
        !
        !    Function:  Contains code for defining contact nodes and computing the node normals
        !
        !    Note: This function is only for 3D case
        !
        !    GOM version 2021 
        !
        !**********************************************************************
        
        implicit none

          integer(INTEGER_TYPE), intent(in) :: FileUnit

          ! local variables
          character :: BName*255
          integer(INTEGER_TYPE) :: IError, ISurf, J, NContSur, NodID, I, IEl, NPileelements, NFoundNodes, CheckNodeID, &
                    ElementID, SideID, AdjacentElement, IVol, ISide, ContElem, INode, JNode, NSideNodes, &
                    CheckNodeI, II, CountSurf, DummyElem, IMat, JMat, NeighbourElement, MaterialID1, MaterialID2, AI, &
                    NContVolElem
          real(REAL_TYPE), dimension(:,:), allocatable :: TotSrfaceAndLoad, ContactVolumes, ContactVolumesDummy
          integer(INTEGER_TYPE), dimension(:,:), allocatable :: ContSurfaceNodes, ContVolSurf
          integer(INTEGER_TYPE), dimension(:), allocatable :: ContctSurfaceElmID, BoundaryElements, Sides
          logical, dimension(:), allocatable :: Element, CheckedLoadedSurface, IsContactSurface
          integer(INTEGER_TYPE) :: NFoundSides, NContactElm

          allocate(Element(Counters%NEl), stat=IError)
          allocate(GlobContElement(Counters%NEl), stat=IError) 
        
          if (allocated(ContactSurfaces) ) then
            deallocate(ContactSurfaces, stat = IError)
          end if

          do
            read(FileUnit, '(A)') BName
            if ((trim(BName)=='$$START_CONTACT_SURFACE')) then
              read(FileUnit, *) NContSur  ! no of element sides forming the contact surface
              if (NContSur>0) then
                allocate(ContactSurfaces(NContSur, 16), stat=IError)
                ContactSurfaces = 0.0 
                 do ISurf= 1, NContSur
                   read(FileUnit,*)  (ContactSurfaces(ISurf,J), J = 1,16)  ! (ElementID, 6 NodeIDs, no. of slaves, friction and adhesion for each slave)
                 end do
              end if
            else if (trim(BName)=='$$START_CONTACT_VOLUME') then 
                read(FileUnit, *) NContVolElem  
                allocate(ContactVolumes(NContVolElem,10), stat=IError)
                allocate(ContactVolumesDummy(NContVolElem,10), stat=IError)
                allocate(ContactMaterialName(NContVolElem,4), stat=IError)
                ContactVolumes = 0
                ContactVolumesDummy = 0
                do IVol = 1, NContVolElem
                    read(FileUnit,*) ContactVolumes(IVol,1:2), ContactMaterialName(Ivol,1), ContactVolumes(IVol,3:4), &
                    ContactMaterialName(Ivol,2), ContactVolumes(IVol,5:6), ContactMaterialName(Ivol,3), &
                    ContactVolumes(IVol,7:8), ContactMaterialName(Ivol,4), ContactVolumes(IVol,9:10)
                end do 
            else if (trim(BName)=='$$FINISH') then
                  EXIT
             end if 
          end do 
          ! rearrange contact properties due to MaterialIndex
          do IMat = 1, 4 !loop over "sliding" materials
            do JMat = 1, CalParams%NumberOfMaterials  !loop over all materials  
              if (ContactMaterialName(1,IMat)==MatParams(JMat)%MaterialName) then
                ContactVolumesDummy(:, MatParams(JMat)%MaterialIndex * 2 + 1 : MatParams(JMat)%MaterialIndex * 2 + 2) &      
                = ContactVolumes(:,(IMat*2+1):(IMat*2+2))     
              end if
            end do
          end do  
          ContactVolumes(:,3:10) = ContactVolumesDummy(:,3:10)
           
          allocate(ContVolSurf(4*NContVolElem,6), stat=IError) 
          ContVolSurf = 0
          Element = .False.
          GlobContElement = .False.
          NPileelements = 0
          NContSur = 0
          DummyElem = 0
          do IEl = 1, NContVolElem
             NPileelements = NPileelements+1 
             Element(ContactVolumes(IEl,1)) = .True. ! pile elemen      
             GlobContElement(ContactVolumes(IEl,1)) = .True.
               if (abs(ElementMaterialID(ContactVolumes(IEl,1)))/=CalParams%MovingMesh%StructureMaterialID) then
                 call GiveWarning('Body contact can not be assigned to more than one material')
               end if
             PileEntityMaterial = abs(ElementMaterialID(ContactVolumes(IEl,1)))
             if  (IsActiveElement(ContactVolumes(IEl,1))) then ! element belongs to the contact volume and is active
               DummyElem = ContactVolumes(IEl,1)    
               do ISide = 1, 4
                 DummyElem = ContactVolumes(IEl,1)  
                 ContElem = BoundaryElementSurface(DummyElem, ISide, &
                                                IsActiveElement, Counters%NEl)
        !Returns 1 if the element adjacent to side ISide of IElement is deactivated.
                 if (ContElem==1) then 
                   do INode = 1, 6 ! Loop over nodes of side I
                     CheckNodeI = ElementConnectivities10Node(DetermineSideNodesTetrahedronHOE(ISide, INode), ContactVolumes(IEl,1))
                     ContVolSurf((IEl+(ISide-1)+3*(IEl-1)), INode) = CheckNodeI
                   end do
                   NContSur = NContSur + 1
        !Returns 0 if side ISide of IElement lies inside a group of activated elements.  
                 else if (ContElem==0) then
                   NeighbourElement = GetAdjacentElement(DummyElem,ISide)  
                   if ((NeighbourElement/=0).or.(NeighbourElement/=-999)) then
                     MaterialID1 = ElementMaterialID(NeighbourElement)
                     MaterialID2 = ElementMaterialID(DummyElem)! get ParticleID from first particle in DummyElem
                   end if  
                      if (MaterialID1/=MaterialID2) then
                        do JNode = 1, 6  
                          CheckNodeI = ElementConnectivities10Node(DetermineSideNodesTetrahedronHOE(ISide, JNode), ContactVolumes(IEl,1))
                          ContVolSurf((IEl+(ISide-1)+3*(IEl-1)),JNode) = CheckNodeI  
                        end do
                        NContSur = NContSur + 1
                      end if
                 end if   
               end do
             end if
          end do ! loop over elements
          
          allocate(ContactSurfaces(NContSur, 16), stat=IError)
          ContactSurfaces = 0
          CountSurf = 0
          
          do II = 1, 4*NContVolElem
              if (ContVolSurf(II,1)/=0) then
                  AI = INT(II/4)+1
                  CountSurf = CountSurf + 1
                  ContactSurfaces(CountSurf,3:8) = ContVolSurf(II,1:6)
                  ContactSurfaces(CountSurf,9:16) = ContactVolumes(AI,3:10)
              end if
          end do        
             
          allocate(IsContactSurface(NContSur), stat=IError)
          IsContactSurface = .true.
          
          allocate(ContSurfaceNodes(NContSur, 6), stat=IError)
          ContSurfaceNodes= 0.0
          
          do I = 1, NContSur
              do J = 1, 6
                  ContSurfaceNodes(I,J) = ContactSurfaces(I,J+2)
              end do
          end do
          
          
          allocate(ContctSurfaceElmID(NContSur), stat=IError)
          ContctSurfaceElmID= 0
          
          allocate(Sides(NContSur), stat=IError)
          Sides= 0         
                    
          
         if (CalParams%ApplyContactMeshBoundary) then

            allocate(CheckedLoadedSurface(NContSur), stat = IError)
            CheckedLoadedSurface = .false.
            NFoundSides = 0
           do IEl = 1, Counters%NEl
             do I = 1, 4
               AdjacentElement = GetAdjacentElement(IEl, I) 
               if (AdjacentElement==0) then  
                 do ISurf= 1, NContSur
                   if (.not.CheckedLoadedSurface(ISurf)) then
                     NFoundNodes = 0
                     do J = 1, 6 ! Loop over nodes of side I
                       CheckNodeID = ElementConnectivities10Node(DetermineSideNodesTetrahedronHOE(I, J), IEl)
                       if ( (CheckNodeID== ContSurfaceNodes (ISurf, 1)).or. &
                         (CheckNodeID== ContSurfaceNodes (ISurf, 2)).or. &
                         (CheckNodeID== ContSurfaceNodes (ISurf, 3)).or. &
                         (CheckNodeID== ContSurfaceNodes (ISurf, 4)).or. &
                         (CheckNodeID== ContSurfaceNodes (ISurf, 5)).or. &
                         (CheckNodeID== ContSurfaceNodes (ISurf, 6)))then
                         NFoundNodes = NFoundNodes + 1
                       end if
                     end do
                     if (NFoundNodes==6) then
                       ContctSurfaceElmID(ISurf) = IEl
                       CheckedLoadedSurface(ISurf) = .true.
                       Sides(ISurf) = I
                       NFoundSides = NFoundSides + 1
                       EXIT
                     end if
                   end if
                 end do
               end if
               if (NFoundSides==NContSur) then
                 EXIT
               end if
             end do
             if (NFoundSides==NContSur) then
               EXIT
             end if
           end do
           
           deallocate(CheckedLoadedSurface, stat = IError)
       else
         do ISurf= 1, NContSur
           do IEl = 1, Counters%NEl
             if (Element(IEl)) then 
               do I = 1, 4 ! Loop over sides
                 NFoundNodes = 0
                  do J = 1, 6 ! Loop over nodes of side I
                    CheckNodeID = ElementConnectivities10Node &
                                  (DetermineSideNodesTetrahedronHOE(I, J), IEl)
                if ( (CheckNodeID== ContSurfaceNodes (ISurf, 1)).or. &
                     (CheckNodeID== ContSurfaceNodes (ISurf, 2)).or. &
                     (CheckNodeID== ContSurfaceNodes (ISurf, 3)).or. &
                     (CheckNodeID== ContSurfaceNodes (ISurf, 4)).or. &
                     (CheckNodeID== ContSurfaceNodes (ISurf, 5)).or. &
                     (CheckNodeID== ContSurfaceNodes (ISurf, 6)))then
                   NFoundNodes = NFoundNodes + 1
                 end if
                 end do
                 if (NFoundNodes==6) then
                  ContctSurfaceElmID(ISurf) = IEl
                  Sides(ISurf) = I
                  EXIT
                 end if
               end do ! sides
              end if ! pile
            end do ! elements
          end do ! contact surfaces
       end if   
       
       allocate(ContactSurfaceElmNodes(NContSur, 8), stat=IError)
       ContactSurfaceElmNodes= 0
        
       NContactElm = NContSur
       
       do ISurf = 1, NContSur
          ContactSurfaceElmNodes (ISurf, 1)= &
                                    INT(ContctSurfaceElmID(ISurf)) ! element ID
            ContactSurfaceElmNodes (ISurf, 2)= Sides (ISurf)
          do J = 1, 6
            ContactSurfaceElmNodes (ISurf, J+2)= &
                                              ContSurfaceNodes(ISurf, J)
          end do
       end do
       
       
          if (allocated(TotSrfaceAndLoad) ) then
            deallocate(TotSrfaceAndLoad, stat = IError)
          end if
          
          if (allocated(Sides) ) then
            deallocate(Sides, stat = IError)
          end if
          
          if (allocated(IsContactSurface) ) then
            deallocate(IsContactSurface, stat = IError)
          end if
          
          if (allocated(ContSurfaceNodes) ) then
            deallocate(ContSurfaceNodes, stat = IError)
          end if
          
          if (allocated(Element) ) then
            deallocate(Element, stat = IError)
          end if
          
          if (allocated(ContctSurfaceElmID) ) then
            deallocate(ContctSurfaceElmID, stat = IError)
          end if

            if (allocated(ContactNodeNormals)) then
              deallocate(ContactNodeNormals, stat = IError)
            end if
            if (allocated(InterfaceNodes)) then
              deallocate(InterfaceNodes, stat = IError)
            end if

            if (NContactElm==0) RETURN

            allocate(ContactNodeNormals(3, Counters%NodTot), stat = IError)
            ContactNodeNormals = 0.0
            allocate(InterfaceNodes(Counters%NodTot), stat = IError)
            InterfaceNodes = .false.
       
               do I = 1, NContactElm
                do J = 1, 3 ! corner nodes
                  NodID = ContactSurfaceElmNodes (I, J+2)
                  InterfaceNodes(NodID) = .true.
                end do
               end do
            NSideNodes = GetNSideNodes(ELEMENTNODES)

            allocate(BoundaryElements(Counters%NEl), stat = IError)
            BoundaryElements = 0
              do I = 1, NContactElm
               ElementID = ContactSurfaceElmNodes (I, 1)
               SideID = ContactSurfaceElmNodes (I, 2)
               BoundaryElements(ElementID) = SideID
              end do

            call DetermineNodeLocalCS(Counters%NodTot, NDIM, &
                                      ELEMENTNODES, Counters%NEl, &
                                      NSideNodes, &
                                      NodalCoordinates, ElementConnectivities, &
                                      BoundaryElements, &
                                      ContactNodeNormals)

            NInterfaceNodes = 0
            do I = 1, Counters%NodTot
              if (InterfaceNodes(I)) then
                NInterfaceNodes = NInterfaceNodes + 1
                if (CalParams%ApplyContactMeshBoundary) then
                  if ((ContactNodeNormals(1, I)>0.0).and.((NodalCoordinates(I, 2))>1E-5)) then
                    ContactNodeNormals(2, I) = 0.0
                    ContactNodeNormals(1:3, I) = VectorNorm(ContactNodeNormals(1:3, I), 3) ! generalise dimension
                  end if
                
                  ContactNodeNormals(1:3, I) = -1.0 * ContactNodeNormals(1:3, I)
                end if
              end if
            end do

            deallocate(BoundaryElements, stat = IError)               
               
    end subroutine ContactSurfaceNodes_v2021

    

    
    !*************************************************************************************
    !    SUBROUTINE: InitialiseContactArrays
    !
    !    DESCRIPTION:
    !>   To initialise the arrays relate to contact
    !
    !>   @note : Notes
    !
    !>   @param[in] ParameterName : ParameterDescription
    !>   @param[out] ParameterName : ParameterDescription
    !>   @param[inout] ParameterName : ParameterDescription
    !
    !*************************************************************************************
    subroutine InitialiseContactArrays()

    implicit none

    integer(INTEGER_TYPE) :: IError

    if (CalParams%ApplyContactAlgorithm) then
        allocate(LumpedMassSys(Counters%N), stat = IError)                                           ! System Lumped mass array has dimension number of degrees of freedom
        allocate(TotalVelocitySys(Counters%N), stat = IError)                                        ! System Velocity array
        allocate(AccelerationSys(Counters%N), stat = IError)                                         ! System Acceleration array
        allocate(ContactNodeNormals(NVECTOR, Counters%NodTot), stat = IError)                        ! Array to store the node normals
        allocate(InterfaceNodes(Counters%NodTot), stat = IError)                                     ! Mask (true and false values) False for nodes not in the contact interface and true for the nodes in the contact interface
        allocate(NodeElement(Counters%NodTot, 100), stat = IError)                                   ! Given a node the array returns the elements that contain that node ! Has dimension (total number of nodes by 100). Not sure why it's 100
        allocate(InterfaceNodesAdhesion(Counters%NodTot,CalParams%NumberOfMaterials), stat = IError) ! Stores each of the material adhesion properties for each node
        allocate(ContactSurfaceSoilElements(Counters%NEl), stat = IError)                            ! Mask that is true for elements that are on the contact surface (line in 2d??)
    else
        allocate(LumpedMassSys(1), stat = IError)
        allocate(TotalVelocitySys(1), stat = IError)
        allocate(AccelerationSys(1), stat = IError)
        allocate(ContactNodeNormals(1, 1), stat = IError)
        allocate(NodeElement(1, 1), stat = IError)
        allocate(InterfaceNodes(1), stat = IError)
        allocate(InterfaceNodesAdhesion(1,1), stat = IError)
        allocate(ContactSurfaceSoilElements(1), stat = IError)
    end if

    LumpedMassSys = 0.0
    TotalVelocitySys = 0.0
    AccelerationSys = 0.0
    ContactNodeNormals = 0.0
    InterfaceNodes = .false.
    NodeElement = 0
    InterfaceNodesAdhesion = 0.0
    ContactSurfaceSoilElements = .false.

    end subroutine InitialiseContactArrays

    !**********************************************************************
    !   SUBROUTINE: ApplyMPMDynamicContact
    !   DESCRIPTION:
    !   Do velocity correction to entity boundary nodes
    !
    ! Note : This function works only if Number of Layers = 1
    !
    !**********************************************************************
    subroutine ApplyMPMDynamicContact()

    implicit none

    ! local variables
    integer(INTEGER_TYPE), dimension(Counters%nEntity,Counters%NodTot) :: EntityNodes      !for each entity a list of the global nodes belonging to elements containing particles belonging to this entity
    integer(INTEGER_TYPE), dimension(CalParams%NumberOfMaterials,Counters%NodTot) :: MaterialNodes
    integer(INTEGER_TYPE), dimension(NVECTOR) :: IDof
    integer(INTEGER_TYPE) :: I, nn, ient, n

    real(REAL_TYPE), dimension(NVECTOR) :: Unormal      !unit normal vector at node for specific entity
    real(REAL_TYPE), dimension(NVECTOR) :: VelDiff      !Difference in entity and system velocity
    real(REAL_TYPE), dimension(3) :: CrossP      !cross product
    real(REAL_TYPE), dimension(3) :: Omega      !
    real(REAL_TYPE) :: AdhesionFactor, FrictionCoef, AdhesionCoef
    real(REAL_TYPE) :: DiffDotUn, Traction, vsize, Dmu, emax, ScaleFactor

    logical :: ApplyContact

    if (.not.(NFORMULATION==1)) RETURN
    if (.not.CalParams%ApplyContactAlgorithm) RETURN

    ! set up EntityNodes matrix for current particle configuration
    call SetEntityNod(EntityNodes)
    call SetMaterialNod(MaterialNodes)

    ! loop through all nodes in model, in global coordinate system
    DO nn=1,Counters%NodTot

        if (InterfaceNodes(nn)) then  ! in case of pile driving the contact surface is defined

            !get storage position of node degree-of-freedoms
            do I = 1, NVECTOR
                IDof(I) = ReducedDof(nn) + I
            end do

            !loop through all entities
            do ient=1,Counters%nEntity
                !
                if (EntityNodes(ient,nn) == 1) then             !this node (nn) belongs to the entity (ient)

                    ApplyContact = .false.
                    do I = 1, NVECTOR
                        if (TotalVelocitySoil(IDof(I), ient) /= TotalVelocitySys(IDof(I))) then
                            ApplyContact = .true.
                            exit
                        end if
                    end do

                    if (ApplyContact) then   !there is contact at this node because the entity velocity is different from system velocity

                        Unormal(:) = ContactNodeNormals(:, nn)

                        if (ient == SOFT_ENTITY) then ! if it is a soil node (entity 1) then change direction
                            Unormal = - Unormal
                        end if

                        DiffDotUn = 0.0
                        do I = 1, NVECTOR
                            ! Calculate velocity difference between entity and system
                            VelDiff(I) = TotalVelocitySoil(IDof(I), ient) - TotalVelocitySys(IDof(I))
                            ! VelDiff.dot.Unormal
                            DiffDotUn = DiffDotUn + VelDiff(I) * Unormal(I)
                        end do

                        !get the traction at the node for this entity if the modified contact model is used
                        if (CalParams%ApplyTractionContact) then               !MODIFIED (TRACTION) CONTACT MODEL IS USED
                            call EntityNodeTraction(nn, ient, Unormal, Traction)
                        else                                                   !STANDARD CONTACT MODEL IS USED
                            Traction = 0.0                                       !set equal to zero so that (DiffDotUn > 0) is only criteria used below
                        end if
                        !
                        !If the modified contact model is used, correction is done when Traction < 0 (compression)
                        ! The exeption is when it is the first contact and all stresses are still zero:
                        !    - this will give Traction = 0
                        !    - In this case the standard contact model is used: (TotalVelocitySoil-TotalVelocitySys).dot. Unormal > 0
                        !
                        !If the standard contact model is used, Traction is set to zero (above) and the only condition used...
                        ! ...to check for contact is:  (TotalVelocitySoil-TotalVelocitySys).dot. Unormal > 0
                        !
                        do I = 1, CalParams%NumberOfMaterials-1
                            if (MaterialNodes(I,nn)==1) then
                                AdhesionCoef = InterfaceNodesAdhesion(nn,I)
                                exit
                            else
                                AdhesionCoef=1000
                            end if
                        end do
                        ! Adhesion
                        AdhesionFactor = AdhesionCoef * CalParams%TimeIncrement / (LumpedMassDry(IDof(1), IEnt) * (DiffDotUn + TINY))

                        if (DiffDotUn>0.0) then ! Standard formulation: DiffDotUn greater than 0 means approaching bodies

                            if (NVECTOR == 2) then ! 2D case
                                ! VelDiff x Unormal
                                CrossP(1) = 0.0
                                CrossP(2) = 0.0
                                CrossP(3) = VelDiff(1) * Unormal(2) - VelDiff(2) * Unormal(1)
                            elseif (NVECTOR == 3) then ! 3D case
                                !VelDiff x Unormal
                                CrossP(1) = VelDiff(2) * Unormal(3) - VelDiff(3) * Unormal(2)
                                CrossP(2) = VelDiff(3) * Unormal(1) - VelDiff(1) * Unormal(3)
                                CrossP(3) = VelDiff(1) * Unormal(2) - VelDiff(2) * Unormal(1)
                            else
                                call GiveError('Dimension is not correct. It must be 2 or 3.')
                            end if

                            vsize = sqrt( CrossP(1)*CrossP(1) + CrossP(2)*CrossP(2) + CrossP(3)*CrossP(3) )

                            !CrossP/ |CrossP|
                            Omega = 0.0
                            if (vsize /= 0) then
                                Omega = CrossP/vsize
                            end if

                            do I = 1, CalParams%NumberOfMaterials-1
                                if (MaterialNodes(I,nn)==1) then
                                    FrictionCoef = FrictionVector(nn,I)
                                    exit
                                else !node belongs to contact but no material is in neighbour element
                                    FrictionCoef = 10 ! A high value will enforce that the solution is vsys
                                end if
                            end do

                            Dmu = min(FrictionCoef + AdhesionFactor, vsize/DiffDotUn) ! Slope larger than friction coefficient, then bring back to limit (max tangential force)

                            !Unormal X Omega - tangential direction normal to n-t-plane
                            if (NVECTOR == 2) then ! 2D case
                                !VelDiff x Unormal
                                CrossP(1) = Unormal(2) * Omega(3)
                                CrossP(2) = - Unormal(1) * Omega(3)
                                CrossP(3) = Unormal(1) * Omega(2) - Unormal(2) * Omega(1)
                            elseif (NVECTOR == 3)  then ! 3D case
                                !VelDiff x Unormal
                                CrossP(1) = Unormal(2) * Omega(3) - Unormal(3) * Omega(2)
                                CrossP(2) = Unormal(3) * Omega(1) - Unormal(1) * Omega(3)
                                CrossP(3) = Unormal(1) * Omega(2) - Unormal(2) * Omega(1)
                            else
                                call GiveError('Dimension is not correct. It must be 2 or 3.')
                            end if

                            do I = 1, NVECTOR
                                ! correct the velocity for soil nodes
                                TotalVelocitySoil(IDof(I),ient) =  TotalVelocitySoil(IDof(I),ient) - DiffDotUn*( Unormal(I) + Dmu*CrossP(I) )
                                ! with the corrected, get the final acceleration
                                AccelerationSoil(IDof(I),ient) = (TotalVelocitySoil(IDof(I),ient) - TotalVelocitySoilPrevious(IDof(I),ient)) / CalParams%TimeIncrement
                            end do
                            if ( CalParams%ApplyContactVelocityScaling ) then
                                call MaxStrainRate(nn, TotalVelocitySoil(IDof(1):IDof(NVECTOR),ient), &
                                    TotalVelocitySoilPrevious(IDof(1):IDof(NVECTOR),ient), emax, &
                                    LumpedMassDry(IDof(1):IDof(NVECTOR), IEnt), LumpedMassSys(IDof(1):IDof(NVECTOR)))
                                if ( (emax*CalParams%TimeIncrement) > CalParams%ContactScalingFactor ) then
                                    ScaleFactor= CalParams%ContactScalingFactor/ (emax*CalParams%TimeIncrement)
                                    do I = 1, NVECTOR
                                        ! correct the velocity for soil nodes
                                        TotalVelocitySoil(IDof(I),ient) =  TotalVelocitySoilPrevious(IDof(I),ient) - ScaleFactor*DiffDotUn*( Unormal(I) + Dmu*CrossP(I) )
                                        ! with the corrected, get the final acceleration
                                        AccelerationSoil(IDof(I),ient) = (TotalVelocitySoil(IDof(I),ient) - TotalVelocitySoilPrevious(IDof(I),ient)) / CalParams%TimeIncrement
                                    end do
                                end if
                            end if
                        else

                            ! Adhesion factor takes into account that in case of separation, the adhesive force might
                            ! still hold back the body. In this case, slip / no-slip are still relevant.
                            if (abs(AdhesionFactor)>=1.0) then
                                ! Bodies stick to each other
                                ! Set velocity of contacting entities to the system velocity
                            end if
                        end if  !(Traction == 0 .and. DiffDotUn > 0 .or......)
                    end if     !(TotalVelocitySoil(IDof,ient) /= TotalVelocitySys(IDof)...)
                end if       !(EntityNodes(ient,nn) = 1)
            end do         !ient=1,nEntity
        end if ! interface node
    end do             !nn=1,NodTot

    if (IsFEMComputation()) then
        TotalVelocitySoil = TotalVelocitySoilPrevious ! global coordinate system (not changed), not used for MPM
    end if

    end subroutine ApplyMPMDynamicContact


    !**********************************************************************
    !   SUBROUTINE: SetEntityNod
    !   DESCRIPTION:
    !   To set the entity data for the current particle configuration
    !
    !
    ! O   EntityNodes   :  For each entity, gives the active nodes
    !
    !**********************************************************************
    subroutine SetEntityNod(EntityNodes)

    implicit none

    integer(INTEGER_TYPE), dimension(Counters%nEntity,Counters%NodTot), intent(inout) :: EntityNodes

    ! local variables
    integer(INTEGER_TYPE) :: iel, iael, iEnt, inode, nn, I

    EntityNodes = 0

    !EntityNodes
    do iael=1,Counters%Nael                       !loop through all elements
        iel = ActiveElement(iael)
        do ient=1,Counters%nEntity                    !loop through all entities
            if (EntityElements(ient,iel) == 1) then !entity contains this element
                do inode=1,ELEMENTNODES              !loop through element nodes
                    nn = ElementConnectivities(inode,iel)                  !get global node number
                    EntityNodes(ient,nn) = 1              !set node as active node for entity
                end do
            end if
        end do
    end do

    if (CalParams%ApplyContactMeshBoundary) then
        do I = 1, Counters%NodTot
            if (InterfaceNodes(I)) then
                EntityNodes(HARD_ENTITY, I) = 1
            end if
        end do
    end if

    end subroutine SetEntityNod

    !**********************************************************************
    !   SUBROUTINE: SetMaterialNod
    !   DESCRIPTION:
    !   To set the material data for the current particle configuration
    !
    !
    ! O   MaterialNodes   :  For each material, gives the active nodes
    !
    !**********************************************************************
    subroutine SetMaterialNod(MaterialNodes)

    implicit none

    integer(INTEGER_TYPE), dimension(CalParams%NumberOfMaterials,Counters%NodTot), intent(inout) :: MaterialNodes

    ! local variables
    integer(INTEGER_TYPE) :: global_element_id, i, material_id, inode, global_node_id, num_active_elements

    MaterialNodes = 0

    num_active_elements = Counters%Nael                     ! Store the number of active elements for clarity
    do i=1,num_active_elements                              ! loop through all active elements
        global_element_id = ActiveElement(i)                ! Get the global element id
        do material_id=1,CalParams%NumberOfMaterials        ! loop through all materials

            if (MaterialElements(material_id,global_element_id) == 1) then          ! material contains this element
                do inode=1,ELEMENTNODES                                             ! loop through element nodes
                    global_node_id = ElementConnectivities(inode,global_element_id) ! get global node number
                    MaterialNodes(material_id, global_node_id) = 1                  ! set node as active node for material
                end do
            end if
        end do
    end do

    end subroutine SetMaterialNod

    !**********************************************************************
    !
    !   SUBTOUTINE:  EntityNodeTraction
    !   DESCRIPTION:
    !   For a given node (inode) and entity (ient), calculate (interpolate from material points) the normal traction
    !
    ! I: inode         - global node number
    !   ient          - entity number
    !   Unormal       - the unit normal vector at node (inode) for entity (ient)
    !
    ! O: Traction       - Traction (scalar)
    !
    !**********************************************************************
    subroutine EntityNodeTraction(inode, ient, Unormal, Traction)

    implicit none

    integer(INTEGER_TYPE), intent(in) :: inode, ient
    real(REAL_TYPE), dimension(NVECTOR),intent(in) :: Unormal
    real(REAL_TYPE), intent(out) :: Traction

    !local variables
    integer(INTEGER_TYPE) i,iel,nelempart,j,nn,iPart,IntGlo,nix
    real(REAL_TYPE), dimension(NTENSOR) :: Stress

    Traction = 0.0   !reset
    ! Note: Maybe this 99 has something to do with the size of the 100 column matrix?? WaveHello
    do i=1,99                                                                 !loop through all elements containing this global node (inode)
        iel = NodeElement(inode,i)                                              !element global number
        if (iel /= 0) then                                                    !NodeElement might contain zero's if node has little elements
            if (EntityElements(ient,iel) == 1) then                             !the enity contains this element
                nelempart = NPartEle(iel)                                           !number of mass points in element
                do j=1,ELEMENTNODES                                              !loop through element nodes
                    nn=ElementConnectivities(j,iel)                                                    !get global node number
                    if (nn == inode) then                                           !this is the element node corresponding to the global node
                        do iPart=1,nelempart                                            !loop through all material points in element
                            IntGlo = GetParticleIndex(iPart,iel)                          !particle global number
                            if (EntityIDArray(IntGlo) == ient ) then               !check if particle belongs to entity
                                Stress = SigmaEffArray(IntGlo,:)

                                if (NTENSOR == 4) then ! 2D case
                                    Traction = Traction + (Stress(1) * Unormal(1) * Unormal(1) + &!...the mass is used as interpolation weight
                                        Stress(2) * Unormal(2) * Unormal(2) + &
                                        Stress(4) * Unormal(1) * Unormal(2)) * &
                                        MassArray(IntGlo) * ShapeValuesArray(IntGlo,j)
                                elseif (NTENSOR == 6) then ! 3D case
                                    Traction = Traction + (Stress(1)*Unormal(1)*Unormal(1) + & !...the mass is used as interpolation weight
                                        Stress(4)*Unormal(1)*Unormal(2) + &
                                        Stress(6)*Unormal(1)*Unormal(3) + &
                                        Stress(4)*Unormal(2)*Unormal(1) + &
                                        Stress(2)*Unormal(2)*Unormal(2) + &
                                        Stress(5)*Unormal(2)*Unormal(3) + &
                                        Stress(6)*Unormal(3)*Unormal(1) + &
                                        Stress(5)*Unormal(3)*Unormal(2) + &
                                        Stress(3)*Unormal(3)*Unormal(3) )* &
                                        MassArray(IntGlo) * ShapeValuesArray(IntGlo,j)
                                else
                                    call GiveError('Dimension is not correct. It must be 2 or 3.')
                                end if
                            end if
                        end do
                        exit        !element node corresponding to global node already found (can only be one) - go to next element
                    end if
                end do
            end if
        end if
    end do

    nix=ReducedDof(inode)+1                      !global storage coordinate of x-val dof of the global node
    Traction = Traction/LumpedMassDry(nix,ient)     !divide the Traction by the nodal mass

    end subroutine EntityNodeTraction

    !*************************************************************************************
    !    SUBROUTINE: CalculateNodeElement
    !
    !    DESCRIPTION:
    !>   Generate the NodeElement array. This array stores which elements are associated with
    !    with an input node. Eg. node=1 is connected to elements 2, 3,4 NodeElement(node) = (2, 3, 4)
    !
    !>   @note : Notes
    !
    !>   @param[in] ParameterName : ParameterDescription
    !>   @param[out] ParameterName : ParameterDescription
    !>   @param[inout] ParameterName : ParameterDescription
    !
    !*************************************************************************************
    subroutine CalculateNodeElement()

    implicit none

    ! Local variables
    integer(INTEGER_TYPE) :: IEl, INode, global_node_id, I

    ! if not applying contact don't need to generate this array
    if (.not.CalParams%ApplyContactAlgorithm) RETURN

    ! Loop over all of the elements
    do IEl = 1, Counters%NEl
        ! Loop over the local node numbering
        do INode=1, ELEMENTNODES
            global_node_id = ElementConnectivities(INode, IEl) ! Get the global node id

            I = NodeElement(global_node_id, 100) + 1 ! Index something??

            NodeElement(global_node_id, I) = IEl

            NodeElement(global_node_id, 100) = I
        end do
    end do

    end subroutine CalculateNodeElement


    !*************************************************************************************
    !    SUBROUTINE: ComputeInterfaceNodesAdhesion
    !
    !    DESCRIPTION:
    !>   Function: Calculates the adhesion and friction at interface nodes
    !
    !>   @note : Wrapper for the 2D and 3D cases
    !
    !>   @param[in] ParameterName : ParameterDescription
    !>   @param[out] ParameterName : ParameterDescription
    !>   @param[inout] ParameterName : ParameterDescription
    !
    !*************************************************************************************
    subroutine ComputeInterfaceNodesAdhesion()
    ! TODO: Generalize the ComputeInterfaceNodesAdhesion for 2d and 3d for any element type and include that here. Get rid of the split between the 2D and 3D case

    implicit none

    ! Get the number of corner nodes and the local index of them
    ! Corner Node mask is IsCornerNode
    ! FIXME: At the momenet these subroutines only work for triangle element in 2D and tetrahederal elements in 3D

    ! If 2D problem
    if (NDIM == 2) then
        call ComputeInterfaceNodesAdhesion2D()

        ! If 3D problem
    elseif (NDIM == 3) then
        call ComputeInterfaceNodesAdhesion3D()
    end if

    end subroutine ComputeInterfaceNodesAdhesion


    !*************************************************************************************
    !    SUBROUTINE: ComputeInterfaceNodesAdhesion2D
    !
    !    DESCRIPTION:
    !>   Calculates the adhesion and friction at interface nodes in 2D.
    !
    !>   @note : Uses ContactSurfaceElmNodes
    !
    !>   @param[in] ParameterName : ParameterDescription
    !>   @param[out] ParameterName : ParameterDescription
    !>   @param[inout] ParameterName : ParameterDescription
    !
    !*************************************************************************************
    subroutine ComputeInterfaceNodesAdhesion2D()

    implicit none

    ! Local variables
    integer(INTEGER_TYPE) :: NodeID
    integer(INTEGER_TYPE) :: IElement, I, J, K, IError
    integer(INTEGER_TYPE), dimension(NVECTOR) :: SurfaceNodes
    integer(INTEGER_TYPE) :: CornerNodes

    real(REAL_TYPE), dimension(NVECTOR) :: Vec1
    real(REAL_TYPE) :: ElementSurface
    real(REAL_TYPE), dimension(:,:), allocatable :: AdhesionVector



    CornerNodes = NDIM ! this holds only for surface element of tetrahedron and triangle

    if (.not.CalParams%ApplyContactAlgorithm) RETURN

    allocate(AdhesionVector(Counters%NodTot, CalParams%NumberOfMaterials), stat=IError)
    allocate(FrictionVector(Counters%NodTot, CalParams%NumberOfMaterials), stat=IError)

    FrictionVector = 0.0
    AdhesionVector = 0.0
    InterfaceNodesAdhesion = 0.0

    do IElement = 1, NContactElm
        do I = 1, CornerNodes
            SurfaceNodes(I) = ContactSurfaceElmNodes(IElement, I + 2)
        end do

        do I = 1, CornerNodes
            Vec1(I) = NodalOriginalCoord(SurfaceNodes(2), I) - NodalOriginalCoord(SurfaceNodes(1), I)
        end do
        ! TODO: Generalize this for general elements and not just triangles WaveHello
        ElementSurface = sqrt(Vec1(1)*Vec1(1)+Vec1(2)*Vec1(2)) !2.0 is number of corner nodes of a line


        if (CalParams%NumberOfMaterials>1) then
            do J = 1, CalParams%NumberOfMaterials-1
                do I = 1, CornerNodes
                    InterfaceNodesAdhesion(SurfaceNodes(I),J) = InterfaceNodesAdhesion(SurfaceNodes(I),J) + ElementSurface
                    if ((CalParams%FricCoef==0).and.(allocated(ContactSurfaces))) then
                        FrictionVector(SurfaceNodes(I),J) = ContactSurfaces(IElement,3+J*2)     !!!!!!!!!!!!!!!!!!!!!!!!!!!check numbers
                    else
                        FrictionVector(SurfaceNodes(I),J) = CalParams%FricCoef
                    end if
                    if ((CalParams%Adhesion==0).and.(allocated(ContactSurfaces))) then
                        AdhesionVector(SurfaceNodes(I),J) = ContactSurfaces(IElement,4+J*2)  !!!!!!!!!!!!!!!!!!!!!!!!!!!check numbers
                    else
                        AdhesionVector(SurfaceNodes(I),J) = CalParams%Adhesion
                    end if
                end do
            end do
        else
            do I = 1, CornerNodes
                InterfaceNodesAdhesion(SurfaceNodes(I),1) = InterfaceNodesAdhesion(SurfaceNodes(I),1) + ElementSurface
                if ((CalParams%FricCoef==0).and.(allocated(ContactSurfaces))) then
                    FrictionVector(SurfaceNodes(I),1) = ContactSurfaces(IElement,7+1*2)       !!!!!!!!!!!!!!!!!!!!!!!!!!!check numbers
                else
                    FrictionVector(SurfaceNodes(I),1) = CalParams%FricCoef
                end if
                if ((CalParams%Adhesion==0).and.(allocated(ContactSurfaces))) then
                    AdhesionVector(SurfaceNodes(I),1) = ContactSurfaces(IElement,8+1*2)     !!!!!!!!!!!!!!!!!!!!!!!!!!!check numbers
                else
                    AdhesionVector(SurfaceNodes(I),1) = CalParams%Adhesion
                end if
            end do
        end if
    end do

    do K = 1, CalParams%NumberOfMaterials
        do NodeID = 1, Counters%NodTot
            if (InterfaceNodes(NodeID)) then
                InterfaceNodesAdhesion(NodeID,K) = InterfaceNodesAdhesion(NodeID,K) * AdhesionVector(NodeID,K)
            end if
        end do
    end do

    end subroutine ComputeInterfaceNodesAdhesion2D


    !**********************************************************************
    !    SUBROUTINE: ComputeInterfaceNodesAdhesion3D
    !    DESCRIPTION:
    !    Calculates the adhesion and friction at interface nodes in 3D.
    !
    !**********************************************************************
    subroutine ComputeInterfaceNodesAdhesion3D()

    implicit none

    ! Local variables
    integer(INTEGER_TYPE) :: NodeID
    integer(INTEGER_TYPE) :: IElement, I, J, K, IError
    integer(INTEGER_TYPE), dimension(NVECTOR) :: SurfaceNodes
    real(REAL_TYPE), dimension(NVECTOR) :: Vec1, Vec2
    real(REAL_TYPE), dimension(NVECTOR) :: Normal
    real(REAL_TYPE) :: ElementSurface
    real(REAL_TYPE), dimension(:,:), allocatable :: AdhesionVector
    integer(INTEGER_TYPE) :: CornerNodes

    ! FIXME: This configuration isn't going to work because the number
    CornerNodes = NDIM ! this holds only for surface element of tetrahedron and triangle

    if (.not.CalParams%ApplyContactAlgorithm) RETURN

    allocate(AdhesionVector(Counters%NodTot, CalParams%NumberOfMaterials), stat=IError)
    allocate(FrictionVector(Counters%NodTot, CalParams%NumberOfMaterials), stat=IError)

    FrictionVector = 0.0
    AdhesionVector = 0.0
    InterfaceNodesAdhesion = 0.0

    do IElement = 1, NContactElm
        do I = 1, CornerNodes
            SurfaceNodes(I) = ContactSurfaceElmNodes(IElement, I + 2)
        end do

        if (NDIM == 3) then ! 3D case
            do I = 1, CornerNodes
                Vec1(I) = NodalOriginalCoord(SurfaceNodes(2), I) - NodalOriginalCoord(SurfaceNodes(1), I)
                Vec2(I) = NodalOriginalCoord(SurfaceNodes(3), I) - NodalOriginalCoord(SurfaceNodes(1), I)
            end do
            Normal = CrossProduct(Vec1, Vec2)
            ElementSurface = Length(Normal, 3) / 2.0 / 3.0 !3.0 is number of vertices (corner nodes) of a triangle
        else ! 2D case
            do I = 1, CornerNodes
                Vec1(I) = NodalOriginalCoord(SurfaceNodes(2), I) - NodalOriginalCoord(SurfaceNodes(1), I)
            end do
            ElementSurface = sqrt(Vec1(1)*Vec1(1)+Vec1(2)*Vec1(2)) / 2.0 !2.0 is number of corner nodes of a line
        endif

        if (CalParams%NumberOfMaterials>1) then
            do J = 1, (CalParams%NumberOfMaterials-1)
                do I = 1, CornerNodes
                    InterfaceNodesAdhesion(SurfaceNodes(I),J) = InterfaceNodesAdhesion(SurfaceNodes(I),J) + ElementSurface
                    if ((CalParams%FricCoef==0).and.(allocated(ContactSurfaces))) then
                        FrictionVector(SurfaceNodes(I),J) = ContactSurfaces(IElement,7+J*2)     !!!!!!!!!!!!!!!!!!!!!!!!!!!check numbers
                    else
                        FrictionVector(SurfaceNodes(I),J) = CalParams%FricCoef
                    end if
                    if ((CalParams%Adhesion==0).and.(allocated(ContactSurfaces))) then
                        AdhesionVector(SurfaceNodes(I),J) = ContactSurfaces(IElement,8+J*2)  !!!!!!!!!!!!!!!!!!!!!!!!!!!check numbers
                    else
                        AdhesionVector(SurfaceNodes(I),J) = CalParams%Adhesion
                    end if
                end do
            end do
        else
            do I = 1, CornerNodes
                InterfaceNodesAdhesion(SurfaceNodes(I),1) = InterfaceNodesAdhesion(SurfaceNodes(I),1) + ElementSurface
                if ((CalParams%FricCoef==0).and.(allocated(ContactSurfaces))) then
                    FrictionVector(SurfaceNodes(I),1) = ContactSurfaces(IElement,7+1*2)       !!!!!!!!!!!!!!!!!!!!!!!!!!!check numbers
                else
                    FrictionVector(SurfaceNodes(I),1) = CalParams%FricCoef
                end if
                if ((CalParams%Adhesion==0).and.(allocated(ContactSurfaces))) then
                    AdhesionVector(SurfaceNodes(I),1) = ContactSurfaces(IElement,8+1*2)     !!!!!!!!!!!!!!!!!!!!!!!!!!!check numbers
                else
                    AdhesionVector(SurfaceNodes(I),1) = CalParams%Adhesion
                end if
            end do
        end if
    end do

    do K = 1, CalParams%NumberOfMaterials
        do NodeID = 1, Counters%NodTot
            if (InterfaceNodes(NodeID)) then
                InterfaceNodesAdhesion(NodeID,K) = InterfaceNodesAdhesion(NodeID,K) * AdhesionVector(NodeID,K)
            end if
        end do
    end do

    end subroutine ComputeInterfaceNodesAdhesion3D


    !*************************************************************************************
    !    SUBROUTINE: MaxStrainRate
    !
    !    DESCRIPTION:
    !>   Determine the max strain rate for scaling the wrong velocities in contact algorithm.
    !
    !>   @note : luisez after Bardenhagen et al (2001)
    !
    !>   @param[in] NodeID : Node ID of interest
    !>   @param[in] EntID : Entity ID of interest
    !>   @param[in] uncorrected_velocity : Previous velocity of soil (without contact algorithm)
    !>   @param[out] emax : maximun strain rate used for corretion
    !>   @param[in] CorrectedVelocity : Corrected velocity with the contact algoritm
    !
    !*************************************************************************************
    subroutine MaxStrainRate(NodeID, CorrectedVelocity, uncorrected_velocity, emax, lumped_masss, SysMass)

    implicit none
    integer(INTEGER_TYPE)              , intent(in):: NodeID
    real(REAL_TYPE), dimension(NVECTOR), intent(in):: uncorrected_velocity, lumped_masss, SysMass, CorrectedVelocity

    real(REAL_TYPE)                    , intent(out):: emax

    ! Local varaibles
    integer(INTEGER_TYPE)::  NelmofNode, ielm, j
    real(REAL_TYPE):: VolElm, Dx
    real(REAL_TYPE), dimension(NVECTOR, 2) :: ej


    !1 step is to find the average volume of elements to find an average length
    NelmofNode = GetNElmOfNode(NodeID)

    VolElm=0

    do ielm= 1, NelmofNode
        VolElm=VolElm + ElementSpace(GetElmIOfNode(NodeID, ielm))
    end do

    VolElm=VolElm/real(NelmofNode)

    !TODO: Need to check if this works for quad elemnts?
    if (NVECTOR==3) then
        Dx= VolElm**(0.333333333333) !Average size of cell sourronfing the element
    else
        Dx = VolElm**(0.5)
    end if

    !2 compute the strain rate
    do j=1, NVECTOR
        ej(j, 1)=(  CorrectedVelocity(j)-uncorrected_velocity(j) )/Dx
        ej(j, 2)= (lumped_masss(j)/(SysMass(j)-lumped_masss(j)))*ej(j, 1)
        ej(j, 1)=abs(ej(j, 1))
        ej(j, 2)=abs(ej(j, 2))
    end do

    emax=0
    do j=1, NVECTOR
        if   (ej(j, 1) > ej(j, 2)) then
            if (ej(j, 1) > emax ) emax= ej(j, 1)
        else
            if (ej(j, 2) > emax ) emax= ej(j, 2)
        end if
    end do
    end subroutine MaxStrainRate
    end module ModMPMDynContact