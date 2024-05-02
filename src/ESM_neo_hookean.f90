
      !*************************************************************************************
      !    SUBROUTINE: ESM_NeoHookean
      ! 
      !    DESCRIPTION:
      !>   Description of subroutine
      !
      !>   @note : Notes
      !
      !>   @param[in] ParameterName : ParameterDescription
      !>   @param[out] ParameterName : ParameterDescription
      !>   @param[inout] ParameterName : ParameterDescription
      !
      !*************************************************************************************
subroutine ESM_NeoHookean(F,  arg2)
    real(REAL_TYPE), dimension(:, :) intent(in) :: F
    type2, intent(out) ::  arg2

    call umat_NeoHookean()
end subroutine ESM_NeoHookean

      !*************************************************************************************
      !    SUBROUTINE: umat_NeoHookean
      ! 
      !    DESCRIPTION:
      !>   Calculates the stress change for a NeoHookean material. Follows the format from
      !>   https://osupdocs.forestry.oregonstate.edu/index.php/Neo-Hookean_Material
      !
      !>   @note : The subroutine refers to 
      !
      !>   @param[in] ParameterName : ParameterDescription
      !>   @param[out] ParameterName : ParameterDescription
      !>   @param[inout] ParameterName : ParameterDescription
      !
      !*************************************************************************************
subroutine umat_NeoHookean(lame_modulus, F,  UJOption, stress)

    real(REAL_TYPE)                 , intent(in) :: lame_modulus
    real(REAL_TYPE), dimension(:, :), intent(in) :: F
    real(REAL_TYPE), dimension(:), intent(inout) :: stress, stress_increment
    real, intent(out) ::  arg2

    integer(INTEGER_TYPE), intent(in) :: UJOption


    ! Local Variables
    integer(INTEGER_TYPE) :: 
    real(REAL_TYPE) ::

    ! $$F = \frac{\partial}{\partial X}(X + u) = Identity + \frac{\partial u}{\partial X}$$
    ! Calc the derminant of F (This relative volume change of the material)

    ! Get the UJ_term depending on the input option

    ! Calc the stress increment

    ! Calc the determinant 
    
end subroutine umat_NeoHookean