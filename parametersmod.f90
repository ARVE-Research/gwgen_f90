module parametersmod
    ! Simple module defining some types and parameters

    implicit none

    integer, parameter :: dp = selected_real_kind(9)      ! 8 byte real
    integer, parameter :: sp = selected_real_kind(4)      ! 4 byte real
    integer, parameter :: i4 = selected_int_kind(9)       ! 10^9 fits into 4 bytes

    integer, parameter, dimension(12) :: ndaymonth = (/ 31,28,31,30,31,30,31,31,30,31,30,31 /) ! number of days in each month

    real(sp), parameter :: Tfreeze = 273.15      ! freezing temperature of freshwater (K)

end module parametersmod
