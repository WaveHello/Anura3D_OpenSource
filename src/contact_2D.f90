! Read in the contact data

module MOD2DCONTACT
   ! Select the variables/subroutines/functions to use from other modules
   use ModGlobalConstants,     only: INTEGER_TYPE, REAL_TYPE, ELEMENTSIDES, ELEMENTVERTICES, ELEMENTTYPE, NVECTOR
   use ModFileIO,              only: FileOpen
   use ModCounters,            only: Counters
   use ModMeshAdjacencies,     only: get_vertex_nodes_2d
   use ModReadCalculationData, only: CalParams
   use ModReadMaterialData,    only: MatParams
   use ModMeshInfo,            only: ElementMaterialID, IsActiveElement, ElementConnectivities, NodalCoordinates
   use ModMeshAdjacencies,     only: BoundaryElementSurface, GetAdjacentElement
   use ModGeometryMath,        only: calc_2d_line_normal, vectorNorm
   implicit none
    
   private
   public :: read_2d_contact_data, ContactSurfaceNodes_2D
   
contains

   subroutine read_2d_contact_data(FileUnit, FileName, contact_volumes, material_names)
      implicit none

      integer(INTEGER_TYPE),              intent(in)  :: FileUnit
      character(len =*),                  intent(in)  :: FileName 
      character(len = 1024), allocatable, intent(out) :: material_names(:, :)  ! Store the material names associated with each element
      real(REAL_TYPE),       allocatable, intent(out) :: contact_volumes(:, :) ! Store the material properties associated with element

      ! Local variables
      character(len = 1024) :: file_line
      integer(INTEGER_TYPE) :: num_contact_cols = 10, max_contact_materials = 4
      integer(INTEGER_TYPE) :: i, num_contact_elements, IError
    
      call FileOpen(FileUnit, trim(FileName)) ! Open the file
      do
         read(FileUnit, '(A)') file_line
         if (trim(file_line)=='$$START_BODY_CONTACT_2D') then

            read(FileUnit, *) num_contact_elements ! Read the number of elements that are in the contact surface

            allocate(contact_volumes( num_contact_elements, num_contact_cols )   , stat=IError) ! Storage for contact
            allocate(material_names(num_contact_elements, max_contact_materials ), stat=IError)  ! Storage for the contact material names

            ! Zero the contact volumes matrix
            contact_volumes = 0

            do i = 1, num_contact_elements
               ! Read the following information...
               read(FileUnit,*) contact_volumes(i, 1:2), & ! Global element id (1) and the number of contacting materials for that element (2)
                  material_names(i, 1),    & ! Name of contacting material 1
                  contact_volumes(i, 3:4), & ! Contacting material 1 friction coeff and adhesion
                  material_names(i, 2),    & ! Name of contacting material 2
                  contact_volumes(i, 5:6), & ! Contacting material 2 friction coeff and adhesion
                  material_names(i, 3),    & ! Name of contacting material 3
                  contact_volumes(i, 7:8), & ! Contacting material 3 friction coeff and adhesion
                  material_names(i, 4),    & ! Name of contacting material 4
                  contact_volumes(i, 9:10)   ! Contacting material 4 friction coeff and adhesion
            end do
            Exit ! Break out of the outer loop
         else if (trim(file_line)=='$$FINISH') then ! Break out of the file if the flags are never found
            EXIT
         end if
      end do
      
      ! Close the file
      close(FileUnit)
   end subroutine read_2d_contact_data

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
   subroutine ContactSurfaceNodes_2D(FileUnit, contact_volumes, material_names, element_mask, contact_surfaces,&
                                     contact_surf_elm_node, interface_nodes, node_normals)
      implicit none
   
      character(len=64),      dimension(:,:),              intent(in)   :: material_names ! matrix containing the Contact Material name
      integer(INTEGER_TYPE),                              intent(in)    :: FileUnit
      integer(INTEGER_TYPE), dimension(:,:), allocatable, intent(out)   :: contact_surf_elm_node
      real(REAL_TYPE),       dimension(:,:), allocatable, intent(out)   :: contact_surfaces, node_normals
      real(REAL_TYPE),       dimension(:,:),              intent(inout) :: contact_volumes
      logical,               dimension(:),  allocatable,  intent(out)   :: element_mask, interface_nodes
      
      ! local variables
      integer(INTEGER_TYPE) :: IError, ISurf, I, J, num_contact_surfaces, NodID, IEl, NFoundNodes,            &
         ElementID, side_id, AdjacentElement, IVol, ISide, is_side_on_free_surf, vertex_id, JNode,   &
         adjacent_element_material_id, current_element_material_id, &
         local_side_id, global_element_id , node_1, node_2
   
      integer(INTEGER_TYPE) :: global_material_id, num_nodes_on_side = 2, element_id_column = 1
      integer(INTEGER_TYPE) :: max_num_contact_materials = 4, & ! Do to hard coding and the current formating the max number of materials a body can have contact with is 4
         width_contact_volumes = 10 ! Has  width of 10 because it stores the Global element id (column 1), number of contacting materials for that element (column 2),
   
      ! the friction coeff and the adhesion for each of the 4 possible materials (2 * 4 = 8 values - other 8 columns)
      integer(INTEGER_TYPE) :: num_element_sides, num_element_vertices, num_vertex_per_side, &
         contact_material_index, material_name_column, JMAT, num_contact_elements, &
         local_vertex_id, global_vertex_id
   
      integer(INTEGER_TYPE), allocatable :: vertex_nodes(:, :)
   
      integer(INTEGER_TYPE), dimension(:,:), allocatable :: ContSurfaceNodes, ContVolSurf
      integer(INTEGER_TYPE), dimension(:)  , allocatable :: contact_surf_element_id, Sides, NNodes
   
      real(REAL_TYPE)      , dimension(:,:), allocatable ::  contact_volumes_dummy
      real(REAL_TYPE)                                    :: vector_2d(2), node_1_coord(2), node_2_coord(2)
   
      logical              , dimension(:)  , allocatable :: IsContactSurface
      integer(INTEGER_TYPE), dimension(:)  , allocatable :: local_corner_node_ids
     
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
      i = SIZE(contact_volumes,1)
      j = SIZE(contact_volumes,2)
      allocate(contact_volumes_dummy(i,j )) ! Make a dummy matrix so reorginization
      num_contact_elements = i ! Store the number of contact elements
      ! rearrange contact properties due to Material Index
      do contact_material_index = 1, max_num_contact_materials  !loop over contacting materials
         material_name_column = contact_material_index * 2 +1  ! Get the index for the column that contains the material name contact_volumes for the iterated material
         do JMat = 1, CalParams%NumberOfMaterials  !loop over all materials
            if (material_names(1,contact_material_index)==MatParams(JMat)%MaterialName) then
               contact_volumes_dummy(:, material_name_column) = contact_volumes(:,material_name_column)
            end if
      
         end do
      end do
      
      !contact_volumes(:,3:10) = contact_volumes_dummy(:,3:10)  !until here 2D equals 3D
      
      ! Come back to this
      
      ! Stores the nodes that need to be checked for contact
      ! Each side of the element is a row and the global nodes ids that make up that side are stored as the column values
      allocate(ContVolSurf(num_element_sides * num_contact_elements, num_vertex_per_side), stat=IError)
      allocate(element_mask(Counters%NEl))
      ContVolSurf = 0
      element_mask = .False.
      num_contact_surfaces = 0
      global_element_id = 0
      
      ! Loop over the elements in base material
      do IEl = 1, num_contact_elements
         global_element_id = contact_volumes(IEl, element_id_column) ! Get the contact element global element id
         element_mask(global_element_id) = .True. ! Set the  mask to true for this elemn  
      
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
               ! Return -1: If on mesh boundary
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
      allocate(contact_surfaces(num_contact_surfaces, 12), stat=IError)
      
      ! Zero the matrix
      contact_surfaces(:, :) = 0
      
      ! init variable to store the number of contact surfaces
      num_contact_surfaces = 0  
      ISide =0
      do i = 1, num_contact_elements
          do j = 1, num_element_sides
              ISide = ISide +1
              if (ContVolSurf(ISide,1)/=0) then
                num_contact_surfaces = num_contact_surfaces + 1 ! increment number of contact surface
                
                ! Store the global element id and the number of contacting materials for that element
                contact_surfaces(num_contact_surfaces, 1:2) = contact_volumes(i, 1:2)
                
                ! Store the global vertex (node) id
                contact_surfaces(num_contact_surfaces, 3:4) = ContVolSurf(ISide, 1:num_vertex_per_side)
                
                ! Store all of the contact material properties
                contact_surfaces(num_contact_surfaces, 5:12) = contact_volumes(i, 3:10)
              end if         
          end do
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
      ContSurfaceNodes(1:num_contact_surfaces, 1:2)  = contact_surfaces(1:num_contact_surfaces, 3:4)
            
      do ISurf= 1, num_contact_surfaces          ! Loop local surface ids
         do IEl = 1, Counters%NEl   ! Loop over all the elements
            if (element_mask(IEl)) then ! If element is in the contact elements
               do I = 1, num_element_sides ! Loop over sides
                  NFoundNodes = 0
      
                  do j = 1, num_vertex_per_side
                     ! Get the local node id
                     local_vertex_id = vertex_nodes(j, I)
      
                     global_vertex_id = ElementConnectivities(local_vertex_id, IEl)
      
                     if ( (global_vertex_id== ContSurfaceNodes (ISurf, 1)) .or. (global_vertex_id== ContSurfaceNodes (ISurf, 2))) then
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
      allocate(contact_surf_elm_node(num_contact_surfaces, 4), stat=IError)
      contact_surf_elm_node= 0
      
      do ISurf = 1, num_contact_surfaces
         contact_surf_elm_node (ISurf, 1)= contact_surf_element_id(ISurf) ! element ID
         contact_surf_elm_node (ISurf, 2)= Sides(ISurf) ! Store the side ids
      
         ! Store the global node ids that are in the side
         do J = 1, 2
            contact_surf_elm_node (ISurf, J+2)= ContSurfaceNodes(ISurf, J)
         end do
      end do
      
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
      
      
      if (num_contact_surfaces==0) RETURN
      
      allocate(node_normals(2, Counters%NodTot), stat = IError) ! used to store the node normals - each row stores the x_i component of the normal
      
      node_normals = 0.0 ! Zero matrix
      
      ! Allocate array to make a mask for which nodes are in the interface
      allocate(interface_nodes(Counters%NodTot), stat = IError)
      interface_nodes = .false. ! Set the intial allocation to false
      
      ! Loop over the contact surface
      do I = 1, num_contact_surfaces
         do J = 1, 2 ! corner nodes
            NodID = contact_surf_elm_node (I, J+2)
            interface_nodes(NodID) = .true.
         end do
      end do
      
      ! This node normal calculation assumes that there's only one material in an element at a time.
      ! bardenhagen using the density gradient to determine the node normals because there could be more than one material in an element
      ! This assumption is reasonable for penetration problems and the tutorial contact problem
      
      ! Loop over the contact surfaces
      do I=1,num_contact_surfaces
         ! I think this is the side id. i need to check this
         side_id = contact_surf_elm_node(I,2)
      
         ElementID = contact_surf_elm_node (I, 1)
      
         ! Use the side id to get the local node ids
         node_1 = vertex_nodes(1, side_id)
         node_2 = vertex_nodes(2, side_id)
      
         ! Get the global node ids
         node_1 = ElementConnectivities(node_1, ElementID)
         node_2 = ElementConnectivities(node_2, ElementID)
      
         node_1_coord = NodalCoordinates(node_1, 1:2)
         node_2_coord = NodalCoordinates(node_2, 1:2)
      
         ! Calc the normal to the line segment formed by node 1 and 2
         vector_2d = calc_2d_line_normal(node_1_coord, node_2_coord)
      
         ! Check if there's an adjacent element, returns 0 if no adjacent element exists
         global_element_id =GetAdjacentElement(contact_surf_elm_node (I, 1), side_id) ! Element id
      
         if (global_element_id /=0) then
            global_material_id=ElementMaterialID(global_element_id ) !Material Index
            if (global_material_id>0) then ! adjacent element is soil, then give preference
               do J=1,2 ! Write the normal for both nodes
                  node_normals(1:2, contact_surf_elm_node(I,2+J))= vector_2d(:)
               end do
            else !If the normal have been assigned then avoid writting the normal. Assumes the first node normal is the one that should be used for the calculation
               ! TODO: Might be worth averaging the node normals, if the node appears more than once
               do J=1,2
                  if (node_normals(1, contact_surf_elm_node(I,(2+J)))==0.0 .and. &
                     node_normals(2, contact_surf_elm_node(I,(2+J)))==0.0) then
                     node_normals(1:2, contact_surf_elm_node(I,(2+J)))= vector_2d
                  end if
               end do
            end if
         end if
      end do
      
      ! Insert manual normals if given in GOM
      if (CalParams%ApplyContactNormalCorrection) then ! if there are manually defined normal vectores
         do J=1, CalParams%NCorrectedNormals !loop number of manually defined vectors (max is 10)
            NodID = CalParams%ManuallyDefinedVectors(J,1) ! Get the node ID
      
            do I=1, NVECTOR !Loop directions
               node_normals(I,NodID) = CalParams%ManuallyDefinedVectors(J,1+I)
            end do
      
            node_normals(1:2,NodID)=VectorNorm(node_normals(1:2,NodID),2)! normalizes the vector
         end do
      end if
   end subroutine ContactSurfaceNodes_2D

! Read in the material parameters
end module MOD2DCONTACT

