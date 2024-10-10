! Read in the contact data


! Check which elements need to have the


!*************************************************************************************
    !    SUBROUTINE: ContactSurfaceNodes_2D
    ! 
    !    DESCRIPTION:
    !>   Determines the contact nodes for 2D probelems after 2021
    !
    !>   @note : Notes
    !
    !>   @param[in] FileUnit : The file unit of the file containing the contact information (GOM file)
    !>   @param[out] ParameterName : ParameterDescription
    !>   @param[inout] ParameterName : ParameterDescription
    !
    !*************************************************************************************
    subroutine ContactSurfaceNodes_2D(FileUnit)
    implicit none

    integer(INTEGER_TYPE), intent(in) :: FileUnit

    ! local variables
    character :: file_line*255
    character(len=64), dimension(:,:), allocatable :: contact_mat_names   ! array including the Contact Material name

    integer(INTEGER_TYPE) :: IError, ISurf, I, J, num_contact_surfaces, NodID, IEl, NFoundNodes, CheckNodeID,           &
                             ElementID, SideID, AdjacentElement, IVol, ISide, is_side_on_free_surf, vertex_id, JNode,   &
                             CheckNodeI, global_material_id, adjacent_element_material_id, current_element_material_id, &
                             local_side_id, IRef

    integer(INTEGER_TYPE) :: global_element_id, num_nodes_on_side = 2, element_id_column = 1
    integer(INTEGER_TYPE) :: max_num_contact_materials = 4, & ! Do to hard coding and the current formating the max number of materials a body can have contact with is 4
                             width_ContactVolumes = 10 ! Has  width of 10 because it stores the Global element id (column 1), number of contacting materials for that element (column 2),
    
    ! the friction coeff and the adhesion for each of the 4 possible materials (2 * 4 = 8 values - other 8 columns)
    integer(INTEGER_TYPE) :: num_element_sides, num_element_vertices, num_vertex_per_side, &
                             contact_material_index, material_name_column, JMAT, num_contact_elements, &
                             local_vertex_id, global_vertex_id, node_1_coord, node_2_coord
    
    integer(INTEGER_TYPE), allocatable :: vertex_nodes(:, :)
    
    integer(INTEGER_TYPE), dimension(:,:), allocatable :: ContSurfaceNodes, ContVolSurf
    integer(INTEGER_TYPE), dimension(:)  , allocatable :: contact_surf_element_id, Sides, NNodes, MatNodes

    real(REAL_TYPE)      , dimension(:,:), allocatable :: TotSrfaceAndLoad, ContactVolumes, ContactVolumesDummy, NNormals, NNormals2,MatNodesCood, DVectors
    real(REAL_TYPE) :: dummy ,dummy2, dummy3

    logical              , dimension(:)  , allocatable :: element_mask, CheckedLoadedSurface, IsContactSurface
    integer(INTEGER_TYPE), dimension(:)  , allocatable :: local_corner_node_ids

        do
            read(FileUnit, '(A)') file_line
            if (trim(file_line)=='$$START_BODY_CONTACT_2D') then
                ! Read the number of elements in the contact surface/Volume ( This is all of the elements that have the contact material assigned to them)
                read(FileUnit, *) num_contact_elements ! Read the number of elements that are in the contact surface


                allocate(ContactVolumes( num_contact_elements, width_ContactVolumes )        , stat=IError) ! Storage for contact
                allocate(ContactVolumesDummy(num_contact_elements, width_ContactVolumes )    , stat=IError) ! dummy of the above
                allocate(contact_mat_names(num_contact_elements, max_num_contact_materials ), stat=IError)  ! Storage for the contact material names
                print *, "Contact material name", contact_mat_names
            
                ContactVolumes = 0
                ContactVolumesDummy = 0

                do IVol = 1, num_contact_elements
                    ! Read the following information...
                    read(FileUnit,*) ContactVolumes(IVol,1:2)   , & ! Global element id (1) and the number of contacting materials for that element (2)
                        contact_mat_names(Ivol,1), & ! Name of contacting material 1
                        ContactVolumes(IVol,3:4)   , & ! Contacting material 1 friction coeff and adhesion
                        contact_mat_names(Ivol,2), & ! Name of contacting material 2
                        ContactVolumes(IVol,5:6)   , & ! Contacting material 2 friction coeff and adhesion
                        contact_mat_names(Ivol,3), & ! Name of contacting material 3
                        ContactVolumes(IVol,7:8)   , & ! Contacting material 3 friction coeff and adhesion
                        contact_mat_names(Ivol,4), & ! Name of contacting material 4
                        ContactVolumes(IVol,9:10)      ! Contacting material 4 friction coeff and adhesion
                end do

            else if (trim(file_line)=='$$FINISH') then
                EXIT
            end if
        end do
        
        ! Allocate variables
        allocate(element_mask(Counters%NEl), stat=IError)         ! Allocate an array to hold the element ids?? TODO: Check that this note is right WaveHello
        allocate(GlobContElement(Counters%NEl), stat=IError) ! Global copy of element

        if (allocated(ContactSurfaces) ) then
            deallocate(ContactSurfaces, stat = IError)
        end if

        if (allocated(ContactSurfaceElmNodes) ) then
            deallocate(ContactSurfaceElmNodes, stat = IError)
        end if

        ! TODO: Generalize this in the future so multiple element types in the same model can be used
        ! Get the number of element sides (storing in local variables so multiple element types can be used in the same model)
        num_element_sides = ELEMENTSIDES

        ! Get the number of element corner nodes (storing in local variable so multiple element types can be used in the same model)
        num_element_vertices = ELEMENTVERTICES

        ! Store the information on the element side ordering
        vertex_nodes = Get_Vertex_Nodes_2d(ELEMENTTYPE)

        ! Get the number of nodes on a side (Number of rows in the vertex_nodes matrix
        ! In general will be 2 but doing this in case there's an element that has more than two vertices per side somehow
        num_vertex_per_side = size(vertex_nodes, 1)
        print *, num_vertex_per_side
        
        ! Store the number of contact elements (just doing this until
        ! Find the material index corresponding to the contact material name
        ! FIXME: Change this so that loops over the number of contact materials and not the max number of contact materials
        do contact_material_index = 1, max_num_contact_materials  !loop over contacting materials
            material_name_column = contact_material_index * 2 +1  ! Get the index for the column that contains the material name ContactVolumes for the iterated material

            do JMat = 1, CalParams%NumberOfMaterials  !loop over all materials
                if (contact_mat_names(1,contact_material_index)==MatParams(JMat)%MaterialName) then
                    ContactVolumesDummy(:, material_name_column) = ContactVolumes(:,material_name_column)
                end if

            end do
        end do

        ContactVolumes(:,3:10) = ContactVolumesDummy(:,3:10)  !until here 2D equals 3D

        ! Come back to this

        ! Stores the nodes that need to be checked for contact
        ! Each side of the element is a row and the global nodes ids that make up that side are stored as the column values
        allocate(ContVolSurf(num_element_sides * num_contact_elements, num_vertex_per_side), stat=IError)

        ContVolSurf = 0
        element_mask = .False.
        GlobContElement = .False.
        num_contact_surfaces = 0
        global_element_id = 0
        
        ! Loop over the elements in base material
        do IEl = 1, num_contact_elements
            global_element_id = ContactVolumes(IEl, element_id_column) ! Get the contact element global element id
            element_mask(global_element_id) = .True. ! Set the  mask to true for this elemn
            GlobContElement(global_element_id ) = .True. ! Set the global contact element mask to true for this element



            ! Check that the material of global_element is the assigned base material (StructureMaterial)
            if (ElementMaterialID(global_element_id)/=CalParams%MovingMesh%StructureMaterialID) then
                call GiveWarning('Body contact can not be assigned to more than one material')
            end if
            
            ! Allocate arr to store local node ids for a side
            
            ! If the element is active (ie. there's a sufficient amount of material points in the element)...
            if  (IsActiveElement(global_element_id)) then
                ! Loop over the element sides
                do ISide = 1, num_element_sides

                    ! Check whether ISide of element global_element_id is on a free surface
                    ! ie. that the only active element it is next two is global_element_id
                    ! Return 1: if on free surface (ie. the side is on the current material boundary of the model or if the adjacent element it's not active))
                    ! Return 0: if ISide is within active elements (not on the free surface)
                    is_side_on_free_surf = BoundaryElementSurface(global_element_id, ISide, IsActiveElement, Counters%NEl)

                    ! The purpose of this control statement is to store the nodes, that need to be checked for contact?
                    ! There's two conditions:
                    ! 1) The side of iEL is exposed to the free surface
                    ! 2) The side of iEl is in contact with an element that contains another material
                    ! Get the global node ids that make up that ISide
                    if (is_side_on_free_surf==1) then

                        ! Get the global node number for those nodes
                        do i = 1, num_vertex_per_side ! loop over the vertex that makes up the current side
                            ! Get the local node id
                            local_vertex_id = vertex_nodes(i, ISide)

                            ! Get the global vertex id
                            global_vertex_id = ElementConnectivities(local_vertex_id, global_element_id)
                            
                            ! Get the row the element id is in
                            j = num_element_sides * IEL - (num_element_sides -1)

                            ! Offset from element row to get the row ISide lives on
                            ContVolSurf( j + ISide-1, i) = global_vertex_id
                        end do

                        ! Increment the number of contact surfaces (ie. the number of sides that need to be checked for contact)
                        num_contact_surfaces = num_contact_surfaces + 1
                        !Returns 0 if side ISide of IElement lies inside a group of activated elements.
                    else if (is_side_on_free_surf==0) then
                        ! Get the global element id of the adjacent element to side ISide
                        AdjacentElement = GetAdjacentElement(global_element_id,ISide)

                        ! Check that there wasn't an error or if there actually wasn't an adjacent element from GetAdjacentElement
                        if ((AdjacentElement/=0).or.(AdjacentElement/=-999)) then

                            ! This values are used to check whether the neighboring element is a different material
                            ! Which means the nodes connecting the current and neighboring element should be checked for contact
                            adjacent_element_material_id = ElementMaterialID(AdjacentElement)
                            current_element_material_id = ElementMaterialID(global_element_id)! get ParticleID from first particle in global_element_id
                        end if

                        ! If the element
                        if (adjacent_element_material_id/=current_element_material_id) then  ! check if
                            ! Loop over the node
                            do i = 1, num_vertex_per_side
                                ! Get the local vertex id
                                local_vertex_id = vertex_nodes(i, ISide)

                                ! Get the global vertex id
                                global_vertex_id = ElementConnectivities(local_vertex_id, global_element_id)

                                ! Get the row index that corresponds to the first side of element IEL in ContVolSurf
                                ! There are num_element_sides per element so this index skips to the firt row index corresponding to the current element
                                j = num_element_sides * IEL - (num_element_sides -1)

                                ! Offset from start element to get the row ISide lives on
                                ContVolSurf( j + ISide-1, i) = global_vertex_id
                            end do
                            ! Increment the number of contact surfaces (ie. the number of sides that need to be checked for contact)
                            num_contact_surfaces = num_contact_surfaces + 1
                        end if
                    end if
                end do
            end if
        end do ! loop over elements
        
        ! Matrix stores
        ! rows are a contact surface
        ! columns are
        ! 1) element id and number of contacting materials
        ! 2) 
        ! 3) 1st global vertex id of side i
        ! 4) 2nd global vertex id of side i
        
        ! 5) Friction of contact material 1
        ! 6) adhesion of contact material 1
        
        ! 7) Friction of contact material 2
        ! 8) adhesion of contact material 2
        
        ! 9) Friction of contact material 3
        ! 10) adhesion of contact material 3
        
        ! 11) Friction of contact material 4
        ! 12) adhesion of contact material 4
         allocate(ContactSurfaces(num_contact_surfaces, 12), stat=IError)
    
        ! Zero the matrix
        ContactSurfaces = 0
    
        ! init variable to store the number of contact surfaces
        num_contact_surfaces = 0

        ! Loop over the element sides
        do i = 1, num_element_sides * num_contact_elements
            ! If the surface has been stored in ContVolSurf then the element id value won't be zero
            if (ContVolSurf(i,1)/=0) then            
                ! Increment the number of contact surfaces
                num_contact_surfaces = num_contact_surfaces + 1
                
                ! Get the row in ContactVolumes that corresponds to this element
                j = i/num_element_sides
                
                ! Store the global element id and the number of contacting materials for that element
                ContactSurfaces(num_contact_surfaces,1:2) = ContactVolumes(j, 1:2)
                
                ! Store the global vertex (node) id
                ContactSurfaces(num_contact_surfaces,3:4) = ContVolSurf(i, 1:num_vertex_per_side)
                
                ! Store all of the contact material properties
                ContactSurfaces(num_contact_surfaces,5:12) = ContactVolumes(j, 3:10)
            end if
        end do
        
        ! Allocate variables
        allocate(IsContactSurface(num_contact_surfaces), stat=IError)
        allocate(ContSurfaceNodes(num_contact_surfaces, num_vertex_per_side), stat=IError)
        allocate(contact_surf_element_id(num_contact_surfaces), stat=IError)
        allocate(Sides(num_contact_surfaces), stat=IError)

        ! Assign initial values
        IsContactSurface = .true. ! Array with length Number of Contact
        ContSurfaceNodes= 0.0
        contact_surf_element_id= 0
        Sides= 0

        ! Store the contact surface nodes (3:4 because there are two nodes
        ContSurfaceNodes(1:num_contact_surfaces, 1:2)  = ContactSurfaces(1:num_contact_surfaces, 3:4)
        
        ! update the module variables
        ! Allocate and store the names in the module varaible
        allocate(ContactMaterialName(num_contact_elements, max_num_contact_materials ), stat=IError)  ! Storage for the contact material names
        ContactMaterialName = contact_mat_names   ! array including the Contact Material name
        NContVolElem = num_contact_elements
        
        do ISurf= 1, num_contact_surfaces          ! Loop local surface ids
            do IEl = 1, Counters%NEl   ! Loop over all the elements
                if (element_mask(IEl)) then ! If element is in the contact elements
                    do I = 1, num_element_sides ! Loop over sides
                        NFoundNodes = 0
                        
                        do j = 1, num_vertex_per_side
                            ! Get the local node id
                            local_vertex_id = vertex_nodes(j, I)

                            global_vertex_id = ElementConnectivities(local_vertex_id, IEl)
                            
                            if ( (CheckNodeID== ContSurfaceNodes (ISurf, 1)) .or. (CheckNodeID== ContSurfaceNodes (ISurf, 2))) then
                                NFoundNodes = NFoundNodes + 1
                            end if
                        end do
                    
                        ! if both corner vertices where found store...
                        if (NFoundNodes==2) then
                            ! The global element id 
                            contact_surf_element_id(ISurf) = IEl
                            
                            ! And the side id
                            Sides(ISurf) = I
                            EXIT
                        end if
                    end do ! sides
                end if ! pile
            end do ! elements
        end do ! contact surfaces
    
        ! Stores:
        ! 1) Element id
        ! 2) Side id
        ! 3 and 4) Node 1 and Node 2
        allocate(ContactSurfaceElmNodes(num_contact_surfaces, 4), stat=IError)
        ContactSurfaceElmNodes= 0
    
        do ISurf = 1, num_contact_surfaces
            ContactSurfaceElmNodes (ISurf, 1)= contact_surf_element_id(ISurf) ! element ID
            ContactSurfaceElmNodes (ISurf, 2)= Sides(ISurf) ! Store the side ids
            
            ! Store the global node ids that are in the side
            do J = 1, 2
                ContactSurfaceElmNodes (ISurf, J+2)= ContSurfaceNodes(ISurf, J)
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

        if (allocated(element_mask) ) then
            deallocate(element_mask, stat = IError)
        end if

        if (allocated(contact_surf_element_id) ) then
            deallocate(contact_surf_element_id, stat = IError)
        end if

        if (allocated(ContactNodeNormals)) then
            deallocate(ContactNodeNormals, stat = IError)
        end if
        if (allocated(InterfaceNodes)) then\
            deallocate(InterfaceNodes, stat = IError)
        end if

        if (num_contact_surfaces==0) RETURN
    
        allocate(ContactNodeNormals(2, Counters%NodTot), stat = IError) ! used to store the node normals - each row stores the x_i component of the normal
        allocate(NNodes(3), stat = IError) ! ??
        allocate(DVectors(3,2), stat = IError) ! ??
        allocate(NNormals(3,2), stat = IError) ! Node Normals?
        allocate(NNormals2(2,2), stat = IError) ! Node Normals 2d ?
        allocate(MatNodes(num_contact_surfaces), stat = IError) ! Material Nodes
        allocate(MatNodesCood(num_contact_surfaces,2), stat = IError) ! Material Node Coordinates

        ContactNodeNormals = 0.0 ! Zero matrix
        
        ! Allocate array to make a mask for which nodes are in the interface
        allocate(InterfaceNodes(Counters%NodTot), stat = IError)
        InterfaceNodes = .false. ! Set the intial allocation to false

        ! Loop over the contact surface
        do I = 1, num_contact_surfaces
            do J = 1, 2 ! corner nodes
                NodID = ContactSurfaceElmNodes (I, J+2)
                InterfaceNodes(NodID) = .true.
            end do
        end do
       
        ! Loop over the contact surfaces
        do I = 1, num_contact_surfaces
            ! Get the element id related to that surface
            ElementID = ContactSurfaceElmNodes (I, 1)
            
            ! Store the side id
            SideID = ContactSurfaceElmNodes (I, 2)
        end do
        
        ! Loop over the contact surfaces
        do I=1,num_contact_surfaces
            ! I think this is the side id. i need to check this
            side_id = ContactSurfaceElmNodes(I,2)
            
            ElementID = ContactSurfaceElmNodes (I, 1)
            
            ! Store the global node ids for
            NNodes(1:3) = ElementConnectivities(1:3, ElementID)
            
            ! What is this doing 
            do J=1,num_side_ver
                ! This is checking if each node for the selected element is in the contact surface element nodes
                if (NNodes(J) == ContactSurfaceElmNodes(I,3)) then
                    ! Store that the equal value was found in column 3
                    NNodes(J) = 3
                elseif (NNodes(J) == ContactSurfaceElmNodes(I,4)) then
                    ! Store that the equal values was found in column 4
                    NNodes(J) = 4
                else
                    ! Store a zero to be used later so that 
                    NNodes(J) = 0
                end if
            end do
                        
            ! This node normal calculation assumes that there's only one material in an element at a time. 
            ! bardenhagen using the density gradient to determine the node normals because there could be more than one material in an element
            ! This assumption is reasonable for penetration problems and the tutorial contact problem
            
            ! Get the tangent vector to the element side
            do J=1,3
                if (NNodes(J)==0) then
                    if (J==1) then
                        node_1_coord = NodalCoordinates(ContactSurfaceElmNodes(I,NNodes(3)),1:2)
                        node_2_coord = NodalCoordinates(ContactSurfaceElmNodes(I,NNodes(2)),1:2)
                        NNormals(side_id,1:2)= node_1_coord - node_2_coord
                        EXIT
                    elseif (J==2) then
                        node_1_coord = NodalCoordinates(ContactSurfaceElmNodes(I,NNodes(1)),1:2) 
                        node_2_coord = NodalCoordinates(ContactSurfaceElmNodes(I,NNodes(3)),1:2)
                        NNormals(side_id,1:2)= node_1_coord - node_2_coord
                        EXIT
                    else
                        node_1_coord = NodalCoordinates(ContactSurfaceElmNodes(I,NNodes(2)),1:2) 
                        node_2_coord = NodalCoordinates(ContactSurfaceElmNodes(I,NNodes(1)),1:2)
                        NNormals(side_id,1:2) = node_1_coord - node_2_coord
                        EXIT
                    end if
                end if
            end do
            
            ! Apply the rotation matrix to tangent vector to get the normal
            ! Flip the first and second row
            DVectors(side_id,1) = NNormals(side_id,2)
            ! Negate the new second row
            DVectors(side_id,2) = -1.0 * NNormals(side_id,1)
            
            ! Store the rotated vector
            NNormals(side_id,1:2)=DVectors(side_id,1:2)
            
            ! Get the unit vector
            NNormals(side_id,1:2)=VectorNorm(NNormals(side_id,1:2),2)
            
            ! Adjacent Element ID
            IRef=GetAdjacentElement(ContactSurfaceElmNodes (I, 1), side_id)
            
            ! Check that the adjacenet element isn't empty
            
            if (IRef/=0) then
                IMat=ElementMaterialID(IRef) !Material Index
                if (IMat>0) then ! adjacent element is soil, then give preference
                    do J=1,2
                        ContactNodeNormals(1:2, ContactSurfaceElmNodes(I,2+J))=&
                               NNormals(ContactSurfaceElmNodes(I,2),1:2)
                    end do
                else !If the normal have been assigned then avoid writting the normal
                    do J=1,2
                        if (ContactNodeNormals(1, ContactSurfaceElmNodes(I,(2+J)))==0.0 .and. &
                            ContactNodeNormals(2, ContactSurfaceElmNodes(I,(2+J)))==0.0) then
                            ContactNodeNormals(1:2, ContactSurfaceElmNodes(I,(2+J)))=&
                            NNormals(ContactSurfaceElmNodes(I,2),1:2)
                        end if
                    end do
                end if
            end if
            
        ! Insert manual normals if given in GOM
        if (CalParams%ApplyContactNormalCorrection) then ! if there are manually defined normal vectores
            do J=1, CalParams%NCorrectedNormals !loop number of manually defined vectors (max is 10)
                NodID = CalParams%ManuallyDefinedVectors(J,1) ! Get the node ID
            
                do I=1, NVECTOR !Loop directions
                    ContactNodeNormals(I,NodID) = CalParams%ManuallyDefinedVectors(J,1+I)
                end do
            
                ContactNodeNormals(1:2,NodID)=VectorNorm(ContactNodeNormals(1:2,NodID),2)! normalizes the vector
            end do
        end if
        end subroutine ContactSurfaceNodes_2D
        
        

    function calc_side_normal(point_1, point_2) Result(side_normal)
    ! Purpose: Calc the normal vector to a point
    
    ! Assume that the point
    end function calc_side_normal
    
        