! Module contains contains functions and subroutines to create and modify arrays
! This is an attempt to recreate some of the functionality available in python
module ModArray
    use ModGlobalConstants, only: INTEGER_TYPE, REAL_TYPE
   
    implicit none
contains
    !*************************************************************************************
    !    FUNCTION: create_int_array
    !
    !    DESCRIPTION:
    !>   Recreates the the functionallity of python range function
    !
    !>   @note : Notes
    !
    !>   @param[in] n1  : First value of the integer array
    !>   @param[in] n2  : Last value of the integer array
    !>   @param[in] dn_ : Step size
    !
    !>   @return array  : integer array
    !
    !*************************************************************************************
    pure function create_int_array(n1,n2,dn_) result(array)
       integer,           intent(in) :: n1,n2
       integer, optional, intent(in) :: dn_
 
       ! Local variables
       integer, allocatable :: array(:)
       integer ::dn, i ! Local step size
 
       dn=1; if(present(dn_))dn=dn_
 
       if(dn<=0)then
          allocate(array(0))
       else
          allocate(array(1+(n2-n1)/dn))
          array=[(i,i=n1,n2,dn)]
       endif
    end function create_int_array
 
    !*************************************************************************************
    !    FUNCTION: linspace
    !
    !    DESCRIPTION:
    !>   Recreates the the functionallity of numpy linspace function
    !
    !>   @note : Notes
    !
    !>   @param[in] start : The starting value of the sequence
    !>   @param[in] end   : The end value of the sequence (End value is inclusive)
    !>   @param[in] num   : number of elements in the array
    !
    !>   @return array  : array of linearly spaced values
    !
    !*************************************************************************************
    pure function linspace(start, end, num) result(array)
       real(REAL_TYPE), intent(in) :: start, end
       integer(INTEGER_TYPE), intent(in) :: num
 
       ! Local Variables
       real(REAL_TYPE), dimension(num) :: array
       real(REAL_TYPE)       :: array_range
       integer(INTEGER_TYPE) :: n, i
 
       array_range = end - start
 
       if (num == 0) return
 
       if (num == 1) then
          array(1) = start
          return
       end if
 
       do i=1, num
          array(i) = start + array_range * (i - 1) / (num - 1)
       end do
    end function linspace
 
 end module ModArray
 