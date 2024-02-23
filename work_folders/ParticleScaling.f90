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

module MODParticleScaling
   !**********************************************************************
   !
   ! Module: Contains the procedures required to scale particle properties
   ! Currently that only includes velocity, stress,
   !
   ! Note: This module
   !
   ! TODO:
   !
   !     $Revision: ????? $
   !     $Date:  2023-12-29 1:32 +0500 (WaveHello, 29 Dec 2023) $
   !
   !**********************************************************************
   use ModGlobalConstants, only : INTEGER_TYPE, REAL_TYPE
                              ! , CalParams%ApplyConvContactStressScaling, &
                              !    CalParams%ApplyConvContactVelocityScaling
   use ModMeshInfo, only : NodalCoordinatesUpd
   use ModGeometryMath, only : check_points_in_box, CheckMinMax

   implicit none
   private ! Makes all functions and subroutines private in the module
   public apply_contact_velocity_scaling, determine_elements_in_box ! Make single subroutine public
contains

  ! Make a global flag that keeps track if the elements in the box have already been calculated
   if (CalParams%ApplyConvContactVelocityScaling .or. CalParams%ApplyConvContactStressScaling) then
     
     !! Check if the elements that need scaling are already calculated
     !if (CalParams%CalcElementsForScaling) then !TODO: Need to create this flag
     !   ! Calculate the elements that need scaling and store in global variable
     !   ElementsForScaling_global = determine_elements_in_box(corner_nodes, node_ids, node_coords, element_ids, element_connectivity)
     !end if
     !
     !if (CalParams%ApplyConvContactVelocityScaling) then
     !   ! Apply contact velocity scaling
     !   ! TODO: Update the name of the variables here 
     !   call apply_contact_velocity_scaling(element_ids = , min_element_dim_arr, particle_connectivity, particle_velocities, time_step, &
     !                                       velocity_scale_factor)
     !end if
   
     ! if (CalParams%ApplyConvContactStressScaling) then
     !    ! Apply stress scaling
        
     ! end if 
   end if

   subroutine apply_particle_scaling(corner_nodes, node_ids, node_coords, element_ids, element_connectivity, &
      vel_flag, stress_flag)
      !! Note: Apparently fortran passes matrices by reference so it doesn't actually take that
      !           much memory to do this

      integer(INTEGER_TYPE), dimension(:), intent(in) :: node_ids, element_ids
      integer(INTEGER_TYPE), dimension(:,:), intent(in) :: element_particle_connectivity

      ! Local variables
      real(REAL_TYPE), dimension(:), allocatable :: min_box_coord, max_box_coord
      integer(INTEGER_TYPE) :: i, num_nodes

      ! Init local variables
      num_nodes = size(node_ids)
      allocate(min_box_coord(num_nodes))
      allocate(max_box_coord(num_nodes))

      ! First determine the box dimensions (set dummy values)
      min_box_coord(:) = 1e30
      max_box_coord(:) = -1e30

      ! Loop through the nodes and determine the bounding min, max values
      do i =1, num_nodes
         CheckMinMax(node_coords(node_ids(i), :), min_box_coord, max_box_coord)
      end do

      ! Determine the nodes in the box
      nodes_in_box = determine_nodes_in_box(node_ids, nodal_coordinates, min_box_coord, max_box_coord, offset_in)

      ! Determine the element_ids in those nodes
      elements_in_box = determine_elements_in_nodes(nodes_in_box, element_ids, element_connectivity)

      ! Pass those elements and the element MP list to the following functions
      apply_particle_scaling(elements, element_particle_connectivity)
      if (CalParams%ApplyConvContactStressScaling) call apply_stress_scaling_to_region()
      if (CalParams%ApplyConvContactVelocityScaling) call apply_contact_velocity_scaling()

   end subroutine apply_particle_scaling


   subroutine apply_contact_velocity_scaling(element_ids, min_element_dim_arr, particle_connectivity, particle_velocities, time_step, &
      velocity_scale_factor)
      integer(INTEGER_TYPE) ,dimension(:)    ,intent(in) :: element_ids
      real(REAL_TYPE)       ,dimension(:)    ,intent(in) :: min_element_dim_arr
      real(REAL_TYPE)       ,dimension(:, :) ,intent(in) :: particle_connectivity
      real(REAL_TYPE)       ,dimension(:, :) ,intent(in) :: particle_velocities
      real(REAL_TYPE) :: time_step, velocity_scale_factor

      ! Local variables
      integer(INTEGER_TYPE) :: i, j, iElem, particle_id
      real(REAL_TYPE)       :: speed_criteria, min_elem_length,
      real(REAL_TYPE)       ,dimension(:), allocatable :: element_particles

      do i = 1, size(element_ids)
         ! Get element id of current element
         iElem = element_ids(i)

         ! Get the min length of the element
         min_elem_length = min_element_dim_arr(iElem)

         ! Calc velocity criteria
         speed_criteria = min_elem_length/time_step

         ! Select the MPs in iElem
         num_particles_in_elem = size(particle_connectivity(iElem, :))
         allocate(element_particle(num_particles_in_elem))

         !TODO - MAKE SURE IM SELECTING THE CORRECT INDEX HERE
         element_particles = particle_connectivity(iElem,:) !- select the right indexes of this matrix

         do j = 1, size(element_MPs)
            ! Get the current MP global id
            particle_id = element_MPs(j)

            ! Get the MP velocity
            particle_velocity = particle_velocities(particle_id)

            ! Get MP speed
            particle_speed = norm2(particle_velocities)

            if (particle_speed > speed_criteria) then
               ! Calc unit norm velocity - need to keep direction correct
               unit_norm_velocity = particle_velocity/particle_speed

               ! set particle velocity
               particle_velocities(particle_id) = velocity_scale_factor * speed_criteria * unit_norm_velocity
            end if

            deallocate(element_particles)
         end do
      end do
   end subroutine apply_contact_velocity_scaling

   function determine_elements_in_box(corner_nodes, node_ids, node_coords, element_ids, element_connectivity) result(elements_in_box)
      implicit none

      integer(INTEGER_TYPE) ,dimension(:)    ,intent(in) :: corner_nodes, node_ids, element_ids
      real(REAL_TYPE)       ,dimension(:, :) ,intent(in) :: node_coords, element_connectivity

      ! Local variables
      real(REAL_TYPE), dimension(:) :: max_box_coord, min_box_coord
      real(REAL_TYPE), dimension(:) :: nodes_in_box, elements_in_box
      integer(INTEGER_TYPE) :: i


      ! Determine the size of the box
      call find_min_max_box_corners(corner_nodes, node_coords, min_box_coord, max_box_coord)

      ! Determine the nodes in the box
      nodes_in_box = determine_nodes_in_box(node_ids, node_coords, min_box_coord, max_box_coord)

      ! Determine the elements in the box
      elements_in_box = determine_elements_in_nodes(nodes_in_box, element_ids, element_connectivity)
   end function determine_elements_in_box


   ! Apply stress scaling

   ! TODO: Move the element and nodes functions to another module (I'm just not sure which one they should go into at this time)
   pure function check_values_in_integer_arr(arr, values) result(inside)
      !!author: WaveHello
      !!date: 12/25/2023
      !! Determines if an array of values are in an array

      !! @todo
      !! Check to see if this can be a pure function
      !! might want to to overload this function to handle over integer types
      !! Adding floats doesn't make too much sense because the floating point procedure
      !! Will throw stuff off
      !! @endtodo

      integer(INTEGER_TYPE), dimension(:), intent(in) :: arr, values

      !Local Variables
      integer(INTEGER_TYPE) :: i
      logical :: inside

      !!True if all elements of values are in arr, false otherwise
      inside= .true.

      do i = 1, size(values)
         if (inside .eqv. .true.) then
            inside = inside .and. any(arr == values(i))
         else
            ! If one value not in arr, return don't need to check others
            return
         end if
      end do

   end function check_values_in_integer_arr

   pure function determine_nodes_in_box(node_ids, nodal_coordinates, min_box_coord, max_box_coord, offset_in) result(nodes_in_box)
      !!author: WaveHello
      !!date: 12/26/2023
      !! From an matrix of nodes and the corresponding coordiantes function determines
      !! the node ids that are in the input box
      !!
      !! nodal_coordinates JUST the nodal coordinates - (Global location)
      real(REAL_TYPE), dimension(:,:), intent(in) :: nodal_coordinates
      integer(INTEGER_TYPE), dimension(:), intent(in):: node_ids
      real(REAL_TYPE), dimension(:), intent(in) :: min_box_coord, max_box_coord
      real(REAL_TYPE), optional, intent(in) :: offset_in

      !Local variables
      real(REAL_TYPE) :: offset
      integer(INTEGER_TYPE), dimension(:), allocatable :: nodes_in_box
      logical, dimension(:), allocatable :: nodes_in_box_mask

      !Set the dimension of nodes_in_box_mask

      ! Check a offset value was passed
      if (present(offset_in)) then
         offset = offset_in
      else
         offset = 1.e-10
      endif

      ! Get the nodes in the box
      nodes_in_box_mask = check_points_in_box(nodal_coordinates, min_box_coord, max_box_coord, offset)

      ! Apply the mask and get nodes in the box
      ! Converting nodal ids to integer
      nodes_in_box = pack(node_ids, nodes_in_box_mask)
   end function determine_nodes_in_box

   pure function determine_elements_in_nodes(node_ids, element_ids, element_connectivity) result(contained_elements)
      !!author: WaveHello
      !!date: 12/26/2023
      !! Searches through a list of nodes to determine the elements that have all there nodes in the list

      integer(INTEGER_TYPE), dimension(:),    intent(in) :: node_ids, element_ids

      !> Rows are nodes and columns are global element ids
      integer(INTEGER_TYPE), dimension(:, :), intent(in) :: element_connectivity

      !> Elements that have their nodes contained in node_ids (Output)
      integer(INTEGER_TYPE), dimension(:), allocatable :: contained_elements

      ! Local variables
      integer(INTEGER_TYPE) :: i, num_elements, current_element_id
      integer(INTEGER_TYPE), dimension(:), allocatable :: current_element_nodes

      logical, dimension(:), allocatable :: elements_mask

      ! Get array size and allocate variables
      num_elements = size(element_ids, 1)
      allocate(elements_mask(num_elements))

      do i =1, num_elements
         !Get the nodes for selected element
         current_element_id = element_ids(i)
         current_element_nodes = element_connectivity(:, current_element_id)

         ! Determine if all element nodes are in nodes_ids
         elements_mask(i)  = check_values_in_integer_arr(node_ids, current_element_nodes)
      end do

      ! Store contained elements
      contained_elements = pack(element_ids, elements_mask)
   end function determine_elements_in_nodes

   ! Determines the min and max corners of a box encompassing an array of nodes
   !TODO: Add subroutine description, calls CheckMinMax()
   !TODO: Move away from CheckMinMax use a pure function to get the min and max box values
   subroutine find_min_max_box_corners(point_ids, point_coords, min_box_coord, max_box_coord)
      implicit none
      integer(INTEGER_TYPE) ,dimension(:)    ,intent(in)  :: point_ids
      real(REAL_TYPE)       ,dimension(:, :) ,intent(in)  :: point_coords
      real(REAL_TYPE)       ,dimension(:)    ,intent(out) :: min_box_coord, max_box_coord

      ! Local variables
      integer(INTEGER_TYPE) :: i, num_nodes

      ! Init variables
      num_nodes = size(point_ids)
      ! Store largest possible number to make sure the min and max evaluates correctly in CheckMinMax
      !! TODO: i don't think this is necessary but need to check
      min_box_coord = huge(min_box_coord)
      max_box_coord = - huge(max_box_coord)

      do i= 1, num_nodes
         call CheckMinMax(point_coords(point_ids(i), :), min_box_coord, max_box_coord)
      end do
   end subroutine find_min_max_box_corners


! Subroutine for modifying the stresses




end module MODParticleScaling
