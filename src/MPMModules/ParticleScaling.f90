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
    use ModGlobalConstants, only :: INTEGER_TYPE, REAL_TYPE
    ! Add files to the geometry module
    use ModMatrixMath, only 
    use ModGeometryMath, only :: check_points_in_box
    implicit none
    
contains
    ! Subroutine for modifying the stresses
    ! Subroutine for modifying the velocities

end module MODParticleScaling