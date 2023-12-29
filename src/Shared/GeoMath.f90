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


!*********************************************************************************************
!
!    MODULE:    ModGeometryMath
!
!    DESCRIPTION:
!>   This module provides basic mathematical routines for dealing with
!>   geometric objects.
!
!     REVISION HISTORY
!>    $Revision: 9707 $
!>    $Date: 2022-04-14 14:56:02 +0200 (do, 14 apr 2022) $
!
!*********************************************************************************************

module ModGeometryMath

   use ModGlobalConstants, only : INTEGER_TYPE, REAL_TYPE, NDIM, TINY

   implicit none

contains
   !! TODO: Update Geomath to use intrinsic functions and use functions from matrix math

   real(REAL_TYPE) function PlanePointDistance(PlaneNormal, PlanePoint, DistantPoint) ! 3D function
      !*************************************************************************************
      !    FUNCTION:     PlanePointDistance
      !
      !    DESCRIPTION:
      !>   Returns the distance between DistantPoint and the plane
      !>   defined by PlaneNormal and PlanePoint. \n
      !>   Distance = <PlaneNormal, (DistantPoint - PlanePoint)>
      !
      !>   @note: only required for 3D
      !
      !>   @param[in] PlaneNormal : Vector normal to the plane of unit length
      !>   @param[in] PlanePoint : Point on the plane
      !>   @param[in] DistantPoint : Point whose distance to the plane is to be returned
      !
      !>   @return PlanePointDistance : Distance between the plane and the point
      !
      !*************************************************************************************
      implicit none

      real(REAL_TYPE), dimension(3), intent(in) :: PlaneNormal, PlanePoint, DistantPoint ! dimension(3) as only required for 3D
      ! local variables
      real(REAL_TYPE), dimension(3) :: Vector

      Vector = DistantPoint - PlanePoint
      PlanePointDistance = DotProduct(PlaneNormal, Vector, 3 ) ! < PN, (DP - PP) >

   end function PlanePointDistance

   !-----------------------------------------------------------------
   function VectorPlaneProjection(PlaneV1, PlaneV2, Vector)
      !*************************************************************************************
      !    FUNCTION:     VectorPlaneProjection
      !
      !    DESCRIPTION:
      !>   Projects Vector onto the plane spanned by PlaneV1 and PlaneV2.
      !
      !>   @note: only required for 3D
      !
      !>   @param[in] PlaneV1 : 3D vector spanning a plane
      !>   @param[in] PlaneV2 : 3D vectors spanning a plane
      !>   @param[in] Vector : Vector whose projection onto the plane is to be determined
      !
      !>   @return VectorPlaneProjection : Projection of Vector on to a plane
      !
      !*************************************************************************************
      implicit none

      real(REAL_TYPE), dimension(3), intent(in) :: PlaneV1, PlaneV2, Vector ! dimension(3) as only required for 3D
      real(REAL_TYPE), dimension(3) :: VectorPlaneProjection ! dimension(3) as only required for 3D
      ! Local variables
      real(REAL_TYPE), dimension(3) :: Normal, UnitVector ! dimension(3) as only required for 3D
      real(REAL_TYPE) :: VectorLength


      Normal = VectorNorm(CrossProduct(PlaneV1, PlaneV2), 3)
      UnitVector = VectorNorm(Vector, 3)
      VectorLength = Length(Vector, 3)

      VectorPlaneProjection = CrossProduct(Normal, CrossProduct(UnitVector, Normal) )

      VectorPlaneProjection = ScalarMultiplication(VectorLength, VectorPlaneProjection, 3)

   end function VectorPlaneProjection

   !-----------------------------------------------------------------
   real(REAL_TYPE) function LinePointDistance(LineNormal, LinePoint, DistantPoint)
      !**********************************************************************
      !
      !    Function:  Returns the distance between DistantPoint and the line
      !               defined by LineNormal and LinePoint (Note: 2D).
      !               Distance = < LineNormal, (DistantPoint - LinePoint) >
      !
      ! I    LineNormal : Vector normal to the line of unit length
      ! I    LinePoint : Point on the line
      ! I    DistantPoint : Point whose distance to the line is to be returned
      !
      ! O   LinePointDistance : Distance between the line and the point
      !
      !**********************************************************************

      implicit none

      real(REAL_TYPE), dimension(:), intent(in) :: LineNormal, LinePoint, DistantPoint

      LinePointDistance = DotProduct(LineNormal, DistantPoint - LinePoint, NDIM) ! < PN, (DP - PP) >

   end function LinePointDistance

   !-----------------------------------------------------------------
   function ScalarMultiplication(Value, Vector, IDim)
      !*************************************************************************************
      !    FUNCTION:     ScalarMultiplication
      !
      !    DESCRIPTION:
      !>   Multiplies Vector with the scalar Value.
      !
      !>   @param[in] Value : Scalar
      !>   @param[in] Vector : Vector
      !>   @param[in] IDim : Dimension of Vector
      !
      !>   @return ScalarMultiplication : Value * Vector
      ! !! TODO: Depricate this. Vectors can be scaled in place no need for a function
      !*************************************************************************************
      implicit none

      integer(INTEGER_TYPE), intent(in) :: IDim
      real(REAL_TYPE), intent(in) :: Value
      real(REAL_TYPE), dimension(IDim), intent(in) :: Vector
      real(REAL_TYPE), dimension(IDim) :: ScalarMultiplication
      ! Local variables
      integer(INTEGER_TYPE) :: I

      do I = 1, IDim
         ScalarMultiplication(I) = Value * Vector(I)
      end do

   end function ScalarMultiplication

   !-----------------------------------------------------------------
   function PolyhedronCentrePoint(Vertices)
      !*************************************************************************************
      !    FUNCTION : PolyhedronCentrePoint
      !
      !    DESCRIPTION :
      !>   Returns the centre point of the polyhedron spanned by the Vertices.
      !
      !>   @note : for 2D and 3D polyhedrons
      !
      !>   @param[in] Vertices : Vertices of the polyhedron
      !
      !>   @return PolyhedronCentrePoint : Centre point of the polyhedron
      !
      !*************************************************************************************
      implicit none

      real(REAL_TYPE), dimension(:, :), intent(in) :: Vertices
      real(REAL_TYPE), dimension(:), allocatable :: PolyhedronCentrePoint
      ! Local variables
      integer(INTEGER_TYPE) :: IVertex, NVertices, IDim

      NVertices = size(Vertices, 1)
      IDim = size(Vertices, 2)
      allocate( PolyhedronCentrePoint(IDim) )

      PolyhedronCentrePoint = 0.0
      do IVertex = 1, NVertices  ! loop vertice nodes of polyhedron
         PolyhedronCentrePoint(:) = PolyhedronCentrePoint(:) + Vertices(IVertex, :)
      end do

      if (NVertices > 0) then
         PolyhedronCentrePoint = PolyhedronCentrePoint / NVertices
      else
         call GiveError("Number of vertices <= 0 in function PolyhedronCentrePoint().")
      end if

   end function PolyhedronCentrePoint

   !-----------------------------------------------------------------
   logical function CheckInsideTetrahedron(A, B, C, D, P)
      !*************************************************************************************
      !    FUNCTION:     CheckInsideTetrahedron
      !
      !    DESCRIPTION:
      !>   Returns .true. if P lies inside the tetrahedron spanned by the vertices A, B, C and D.\n
      !>   (From: Presentation of Prof. Jarek Rossignac, College of Computing, Georgia Institute of Technology, 2008)
      !
      !>   @note: only required for 3D
      !
      !>   @param[in] A : 1st Vertex of a tetrahedron
      !>   @param[in] B : 2nd Vertex of a tetrahedron
      !>   @param[in] C : 3rd Vertex of a tetrahedron
      !>   @param[in] D : 4th Vertex of a tetrahedron
      !>   @param[in] P : Global coordinates of a checked point P
      !
      !>   @note: A, B, C, O, P : dimension(3) as only required for 3D
      !>   @note: A, B, C, D : assumed not to be all co-planar
      !
      !>   @return CheckInsideTetrahedron : True, if P lies inside the tetrahedron
      !
      !*************************************************************************************
      implicit none

      real(REAL_TYPE), dimension(3), intent(in) :: A, B, C, D, P
      ! Local variables
      real(REAL_TYPE), dimension(3) :: DA, DB, DC, DP, PA, PB, PC ! dimension(3) as only required for 3D
      real(REAL_TYPE), dimension(5) :: STerms
      logical, dimension(5) :: Mask
      integer(INTEGER_TYPE) :: TrueCounter, FalseCounter

      DA = A - D
      DB = B - D
      DC = C - D
      DP = P - D
      PA = A - P
      PB = B - P
      PC = C - P

      ! SABCD = < DA, (DB x DC) > greater than zero
      STerms(1) = DotProduct(DA, CrossProduct(DB, DC), 3)
      ! SPBCD = < DP, (DB x DC) > greater than zero
      STerms(2) = DotProduct(DP, CrossProduct(DB, DC), 3)
      ! SAPCD = < DA, (DP x DC) > greater than zero
      STerms(3) = DotProduct(DA, CrossProduct(DP, DC), 3)
      ! SABPD = < DA, (DB x DP) > greater than zero
      STerms(4) = DotProduct(DA, CrossProduct(DB, DP), 3)
      ! SABCP = < PA, (PB x PC) > greater than zero
      STerms(5) = DotProduct(PA, CrossProduct(PB, PC), 3)

      ! Check whether conditions are all true or all false
      Mask = STerms > 0.0
      TrueCounter = Count(Mask)
      Mask = STerms < 0.0
      FalseCounter = Count(Mask)

      CheckInsideTetrahedron = (TrueCounter == 5) .or. (FalseCounter == 5)

   end function CheckInsideTetrahedron

   !-----------------------------------------------------------------
   integer function CheckInsideSubTetrahedron(A, B, C, O, P) ! 3D function
      !*************************************************************************************
      !    FUNCTION:     CheckInsideSubTetrahedron
      !
      !    DESCRIPTION:
      !>   Returns .true. if P lies inside the tetrahedron spanned by the vertices A, B, C and
      !>   center point O of the tetrahedron.
      !
      !>   @note: only required for 3D
      !
      !>   @param[in] A : 1st Vertex of a tetrahedron
      !>   @param[in] B : 2nd Vertex of a tetrahedron
      !>   @param[in] C : 3rd Vertex of a tetrahedron
      !>   @param[in] O : Center point of a tetrahedron
      !>   @param[in] P : Global coordinates of a checked point P
      !
      !>   @note: A, B, C, O, P  : dimension(3) as only required for 3D
      !>   @note: A, B, C : assumed not to be all co-planar
      !
      !>   @return CheckInsideSubTetrahedron : 0, if the track of P does not lie inside the tetrahedron\n
      !>                                           1, if the track of P lies inside the tetrahedron\n
      !>                                           2, if particle P lies inside the tetrahedron
      !*************************************************************************************
      implicit none

      real(REAL_TYPE), dimension(:), intent(in) :: A, B, C, O, P ! dimension(3) as only required for 3D
      ! Local variables
      real(REAL_TYPE), dimension(3) :: OA, OB, OC, OP, PA, PB, PC, temp ! dimension(3) as only required for 3D
      real(REAL_TYPE) :: SOABC, SOPAB, SOPBC, SOPCA, SPABC
      logical, dimension(4) :: Mask

      OA = A - O
      OB = B - O
      OC = C - O
      OP = P - O
      PA = A - P
      PB = B - P
      PC = C - P

      ! SOABC = < OA, (OB x OC) >
      temp = CrossProduct(OB, OC)
      SOABC = DotProduct(OA, temp, 3)
      ! SOPAB = < OP, (OA x OB) >
      temp = CrossProduct(OA, OB)
      SOPAB = DotProduct(OP, temp, 3)
      Mask(1) = ( (SOPAB / SOABC) >= 0)
      ! SOPBC = < OP, (OB x OC) >
      temp = CrossProduct(OB, OC)
      SOPBC = DotProduct(OP, temp, 3)
      Mask(2) = ( (SOPBC / SOABC) >= 0)
      ! SOPCA = < OP, (OC x OA) >
      temp = CrossProduct(OC, OA)
      SOPCA = DotProduct(OP, temp, 3)
      Mask(3) = ( (SOPCA / SOABC) >= 0)
      ! SPABC = < PA, (PB x PC) >
      temp = CrossProduct(PB, PC)
      SPABC = DotProduct(PA, temp, 3)
      Mask(4) = ( (SPABC / SOABC) >= 0)

      CheckInsideSubTetrahedron = 0
      if (Mask(1).and.Mask(2).and.Mask(3) ) then
         ! Track of P lies inside the sub-tetrahedron
         CheckInsideSubTetrahedron = 1
         if (Mask(4) ) then
            ! P itself lies inside the sub-tetrahedron
            CheckInsideSubTetrahedron = 2
         end if
      end if

   end function CheckInsideSubTetrahedron

   !-----------------------------------------------------------------
   logical function CheckInsideTriangle(A, B, C, P) result(res)
      !**********************************************************************
      !
      !    Function:  Returns .true. if P lies inside the triangle
      !               spanned by the vertices A, B and C (Note: 2D).
      !               (Adapted from: Presentation of Prof. Jarek Rossignac, College of Computing, Georgia Institute of Technology, 2008)
      !
      ! I    A, B, C : Vertices of a triangle (they are assumed not to be all co-linear!)
      ! I    P : Global coordinates of a checked point P
      !
      ! O   CheckInsideTriangle : True, if P lies inside the triangle
      !
      !**********************************************************************

      implicit none

      real(REAL_TYPE), dimension(:), intent(in) :: A, B, C, P
      ! local variable
      integer(INTEGER_TYPE) :: IDim

      IDim = size(A)
      res = isInAngle(P, A, B, C, IDim) .and. isInAngle(P, B, C, A, IDim) .and. isInAngle(P, C, A, B, IDim)

   end function CheckInsideTriangle

   !-----------------------------------------------------------------
   logical function CheckInsideQuadrilateral(A, B, C, D, P) result(res)
      !**********************************************************************
      !
      !    Function:  Returns .true. if P lies inside the quadrilateral
      !               spanned by the vertices A, B, C and D (Note: 2D).
      !               (Adapted from: Presentation of Prof. Jarek Rossignac, College of Computing, Georgia Institute of Technology, 2008)
      !
      ! I   A, B, C, D : Vertices of a quadrilateral (they are assumed not to be all co-linear!)
      ! I   P : Global coordinates of a checked point P
      !
      ! O   CheckInsideQuadrilateral : True, if P lies inside the quadrilateral
      !
      !**********************************************************************

      implicit none

      real(REAL_TYPE), dimension(:), intent(in) :: A, B, C, D, P
      ! local variable
      integer(INTEGER_TYPE) :: IDim

      res = isInAngle(P, A, B, C, IDim) .and. isInAngle(P, B, C, D, IDim) .and. isInAngle(P, C, D, A, IDim) .and. isInAngle(P, D, A, B, IDim)

   end function CheckInsideQuadrilateral

   !-----------------------------------------------------------------
   real(REAL_TYPE) function DotProduct(A, B, IDim)
      !*************************************************************************************
      !    FUNCTION:     DotProduct
      !
      !    DESCRIPTION:
      !>   Multiplies Vector with the scalar Value.
      !
      !>   @param[in] A : Vector A
      !>   @param[in] B : Vector B
      !>   @param[in] IDim : Dimensions of Vectors A and B
      !
      !>   @return DotProduct : Dot Product of A, B (< A, B >)
      !   !! TODO: Deprecate this, use intrinsic dot_product function
      !*************************************************************************************
      implicit none

      integer(INTEGER_TYPE), intent(in) :: IDim
      real(REAL_TYPE), dimension(IDim), intent(in) :: A, B
      ! Local variables
      integer(INTEGER_TYPE) :: I

      DotProduct = 0.0
      do I = 1, IDim
         DotProduct = DotProduct + A(I) * B(I)
      end do

   end function DotProduct

   !-----------------------------------------------------------------
   function CrossProduct(A, B) ! 3D function
      !*************************************************************************************
      !    FUNCTION:     Cross Product
      !
      !    DESCRIPTION:
      !>   Multiplies Vector with the scalar Value.
      !
      !>   @note: only required for 3D
      !
      !>   @param[in] A : Vector A
      !>   @param[in] B : Vector B
      !
      !>   @note: A, B : dimension(3) as only required for 3D
      !
      !>   @return CrossProduct : Cross Product of A, B ( A X B )
      !
      !>   @note: CrossProduct : dimension(3) as only required for 3D
      !*************************************************************************************
      implicit none

      real(REAL_TYPE), dimension(3), intent(in) :: A, B
      real(REAL_TYPE), dimension(3) :: CrossProduct

      CrossProduct(1) = A(2) * B(3) - A(3) * B(2)
      CrossProduct(2) = A(3) * B(1) - A(1) * B(3)
      CrossProduct(3) = A(1) * B(2) - A(2) * B(1)

   end function CrossProduct

   !-----------------------------------------------------------------
   real(REAL_TYPE) function Length(Vector, IDim)
      !*************************************************************************************
      !    FUNCTION:     Length
      !
      !    DESCRIPTION:
      !>   Returns the length of Vector.
      !
      !>   @param[in] Vector : Vector
      !>   @param[in] IDim : Dimension of Vector
      !
      !>   @return Length : Length of Vector i.e. SQRT(V1 * V1 + ... + Vn * Vn)
      ! !! TODO: Deprecate this use intrinsic norm2 function
      !*************************************************************************************
      implicit none

      integer(INTEGER_TYPE), intent(in) :: IDim
      real(REAL_TYPE), dimension(IDim), intent(in) :: Vector

      Length = sqrt( DotProduct(Vector, Vector, IDim) )

   end function Length

   !-----------------------------------------------------------------
   function VectorNorm(Vector, IDim)
      !*************************************************************************************
      !    FUNCTION:     VectorNorm
      !
      !    DESCRIPTION:
      !>   Returns Vector with a length of 1.
      !
      !>   @param[in] Vector : Vector
      !>   @param[in] IDim : Dimension of Vector
      !
      !>   @return VectorNorm : Normalizes a vector i.e. Vector / |Vector|
      ! !! TODO: Rename this function. Too easily cofused with getting the norm of a vector
      !           Change to sometinge like GetUnitVector
      !*************************************************************************************
      implicit none

      integer(INTEGER_TYPE), intent(in) :: IDim
      real(REAL_TYPE), dimension(IDim), intent(in) :: Vector
      real(REAL_TYPE), dimension(IDim) :: VectorNorm
      ! Local variables
      integer(INTEGER_TYPE) :: I
      real(REAL_TYPE) :: VectorLength

      VectorLength = Length(Vector, IDim)
      do I = 1, IDim
         VectorNorm(I) = Vector(I) / (VectorLength + TINY)
      end do

   end function VectorNorm

   !-----------------------------------------------------------------
   real(REAL_TYPE) function Distance(PointA, PointB, IDim)
      !*************************************************************************************
      !    FUNCTION:     Distance
      !
      !    DESCRIPTION:
      !>   Returns the distance between PointA and PointB.
      !
      !>   @param[in] PointA : location of first point
      !>   @param[in] PointB : location of second point
      !>   @param[in] IDim : Dimensions of Points A and B
      !
      !>   @return Distance : Distance between two points
      !
      !*************************************************************************************
      implicit none

      integer(INTEGER_TYPE), intent(in) :: IDim
      real(REAL_TYPE), dimension(IDim), intent(in) :: PointA, PointB
      ! Local variables
      real(REAL_TYPE), dimension(IDim) :: Vector
      integer(INTEGER_TYPE) :: I

      do I = 1, IDim
         Vector(I) = PointA(I) - PointB(I)
      end do

      Distance = Length(Vector, IDim)

   end function Distance

   !-----------------------------------------------------------------
   real(REAL_TYPE) function VectorAngle(CrossPoint, A1, A2, IDim)
      !*************************************************************************************
      !    FUNCTION:     Distance
      !
      !    DESCRIPTION:
      !>   Returns the angle between the vectors CrossPoint-A1 and CrossPoint-A2
      !>   in radians.
      !
      !>   @param[in] CrossPoint : location of angle
      !>   @param[in] A1 : Point 1
      !>   @param[in] A2 : Point 2
      !>   @param[in] IDim : Dimensions of Points CrossPoint, A1, and A2
      !
      !>   @return Distance : Angle between two vectors
      !
      !*************************************************************************************
      implicit none

      integer(INTEGER_TYPE), intent(in) :: IDim
      real(REAL_TYPE), dimension(IDim), intent(in) :: CrossPoint, A1, A2
      ! Local variables
      real(REAL_TYPE), dimension(IDim) :: Vector1, Vector2
      real(REAL_TYPE) :: Denominator

      Vector1 = A1 - CrossPoint
      Vector2 = A2 - CrossPoint

      Denominator = Length(Vector1, IDim) * Length(Vector2, IDim)

      if (abs(Denominator) > 1d-10) then
         VectorAngle = acos(DotProduct(Vector1, Vector2, IDim) / Denominator)

         if (IsNaN(VectorAngle)) then
            VectorAngle = 0.0
         end if

      else
         VectorAngle = -999
      end if

   end function VectorAngle

   !-----------------------------------------------------------------
   logical function isInAngle(pointP, pointA, pointB, pointC, iDim) result(res)
      ! is in an angle:
      !  C |
      !    |
      !    |
      !    |____________
      !   A             B
      implicit none
      integer(INTEGER_TYPE), parameter :: Dim3D = 3

      integer(INTEGER_TYPE), intent(in):: iDim
      real(REAL_TYPE) :: dotProduct12, dotProduct13, dotProduct1, dotProduct2
      real(REAL_TYPE) :: crossProdLength
      real(REAL_TYPE), dimension(:) :: pointP, pointA, pointB, pointC
      real(REAL_TYPE), dimension(Dim3D) :: crossProd1, crossProd2, crossProd3
      real(REAL_TYPE), dimension(Dim3D) :: AB, AP, AC

      AB = 0.0
      AP = 0.0
      AC = 0.0

      AB(1:iDim) = pointB(1:iDim) - pointA(1:iDim)
      AP(1:iDim) = pointP(1:iDim) - pointA(1:iDim)
      AC(1:iDim) = pointC(1:iDim) - pointA(1:iDim)
      crossProd1 = crossProduct(AB, AP)
      crossProd2 = crossProduct(AP, AC)
      crossProd3 = crossProduct(AB, AC)

      dotProduct12 = dotProduct(crossProd1, crossProd2, Dim3D)
      dotProduct13 = dotProduct(crossProd1, crossProd3, Dim3D)

      res = dotProduct12 > 0.0 .and. dotProduct13 > 0.0

      if (res) then
         RETURN
      else
         if (abs(dotProduct12) < TINY .or. abs(dotProduct13) < TINY) then
            ! point lies on on one the edges
            crossProdLength = Length(crossProd1, Dim3D)
            if (abs(crossProdLength) < TINY) then
               ! the point lies on line AB
               dotProduct1 = dotProduct(AB, AP, Dim3D)
               if (dotProduct1 > 0.0) then
                  res = .true.
               else
                  res = .false.
               end if
               RETURN
            end if

            crossProdLength = Length(crossProd2, Dim3D)
            if (abs(crossProdLength) < TINY) then
               ! the point lies on line AC
               dotProduct2 = dotProduct(AC, AP, Dim3D)
               if (dotProduct2 > 0.0) then
                  res = .true.
               else
                  res = .false.
               end if
               RETURN
            end if
         end if
      end if

   end function isInAngle

   !-----------------------------------------------------------------
   integer function CheckInsideSubTriangle(A, B, C, P) result(res)
      !**********************************************************************
      !
      !    Function:  Returns .true. if P lies inside the triangle
      !               spanned by the vertices A, B and C (Note: 2D).
      !               (Adapted from: Presentation of Prof. Jarek Rossignac, College of Computing, Georgia Institute of Technology, 2008)
      !
      !     A, B, C : Vertices of a triangle (they are assumed not to be all co-linear!)
      !     P : Global coordinates of a checked point P
      !
      ! O   CheckInsideTriangle : True, if P lies inside the triangle
      !
      !**********************************************************************

      implicit none

      integer(INTEGER_TYPE), parameter :: IS_OUTSIDE  = 0
      integer(INTEGER_TYPE), parameter :: IS_IN_ANGLE = 1
      integer(INTEGER_TYPE), parameter :: IS_INSIDE   = 2

      real(REAL_TYPE), dimension(:), intent(in) :: A, B, C, P
      logical :: isInside
      logical :: isInAngleCAB
      logical :: isInAngleABC
      logical :: isInAngleBCA
      integer(INTEGER_TYPE) :: IDim

      IDim = size(A)
      res = IS_OUTSIDE

      isInAngleCAB = isInAngle(P, A, B, C, IDim)
      isInAngleABC = isInAngle(P, B, C, A, IDim)
      isInAngleBCA = isInAngle(P, C, A, B, IDim)

      isInside = isInAngleCAB .and. isInAngleABC .and. isInAngleBCA

      if (isInside) then
         res = IS_INSIDE
      else
         if (isInAngleBCA) then
            res = IS_IN_ANGLE
         end if
      end if

   end function CheckInsideSubTriangle

   !*************************************************************************************
   !    FUNCTION: check_point_in_box
   ! 
   !    DESCRIPTION:
   !>   Checks whether point lies inside or on the box spanned by <Minimum, Maximum>
   !
   !>   @param[in] point(:) : coordinates of the point that is checked
   !>   @param[in] min_box_coord(:) : minimum coordinates of a bounding box
   !>   @param[in] max_box_coord(:) : maximum coordinates of a bounding box
   !>   @param[in] (optional) offset_in: allowed offset from the box to account for float precision
   !
   !>   @return check_point_in_box : .true. if inside or on boundary, .false. if outside
   !
   !*************************************************************************************
   logical function check_point_in_box(point, min_box_coord, max_box_coord, offset_in)
        
   implicit none

     real(REAL_TYPE), intent(in), dimension(:) :: point
     real(REAL_TYPE), intent(in), dimension(:) :: min_box_coord, max_box_coord
     real, optional :: offset_in ! Could not define type at run time for optional argument

     ! local variables
     integer(INTEGER_TYPE) :: I
     real(REAL_TYPE) :: offset
     
     if(present(offset_in)) then
       offset = offset_in 
     else
       offset = 1.e-10
     endif
     
     if ((size(point) /= size(min_box_coord)) .or. size(point) /= size(max_box_coord)) then
       ! Throw an error
       !TODO: Create an error condition here
     endif

     check_point_in_box = .true.
     do I = 1, size(point)
       check_point_in_box = check_point_in_box .and. ( ( point(I) >= (min_box_coord(I) - offset) ) .and. ( point(I) <= (max_box_coord(I) + offset) ) )                
     end do
     
   end function check_point_in_box

   function check_points_in_box(points, min_box_coord, max_box_coord, offset_in) result(inside)
      !!author: WaveHello
      !!date: 12/25/2023
      !! Determines if a matrix of points (where each point is stored row-wise) is inside a box.
      !! Use the error module to throw size error in this module

      real(REAL_TYPE), dimension(:,:) ,intent(in)  :: points
      !! input min and max values act as an initial guess
      real(REAL_TYPE), dimension(:)   ,intent(in)  :: min_box_coord(:), max_box_coord(:)
      real(REAL_TYPE), optional       ,intent(in)  :: offset_in
      logical,         dimension(:)   ,allocatable :: inside

      ! Local variables
      real(REAL_TYPE), dimension(:), allocatable :: one_point
      real(REAL_TYPE) :: offset
      integer(INTEGER_TYPE) :: num_points, i, point_dimension

      ! Check that the points and box coords are the same dimension
      if ((size(points, 2) /= size(min_box_coord)) .or. size(points, 2) /= size(max_box_coord)) then
         error stop "Incorrect box or point size"
      endif

      ! Check a offset value was passed
      if (present(offset_in)) then
         offset = offset_in
      else
         offset = 1.e-10
      endif

      ! store dimensions for array allocation and looping
      num_points = size(points, 1)
      point_dimension = size(points, 2)

      allocate(inside(num_points))
      allocate(one_point(point_dimension))

      ! Determine the points inside the box
      do i = 1, num_points
         ! Store single point
         one_point = points(i,:)

         ! Condition on truth
         inside(i) = .true.

         ! Check if point is inside the box
         inside(i) = inside(i) .and. all(one_point >= (min_box_coord - offset)) .and. all(one_point <= (max_box_coord+offset))
      end do
   end function check_points_in_box
end module ModGeometryMath
