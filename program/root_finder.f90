! Interface to compute roots using (F)GSL

module root_finder
  use FGSL
  implicit none
  private
  public :: CubicRootClose

contains
  double precision function CubicRootClose(a3, a2, a1, a0, refpoint)
    ! Finds the roots of the cubic equation: a3*x^3 + a2*x^2 + a1*x + a0 = 0
    ! Returns the root closest to refpoint
    implicit none

    double precision, intent(in) :: a1, a2, a3, a0, refpoint
    double precision, dimension(3) :: roots
    double precision :: distance, previous_distance
    real(fgsl_double) :: A, B, C
    integer :: number_of_solutions, i

    A = a2/a3
    B = a1/a3
    C = a0/a3
    
    ! This GSL routine will solve:  x^3 + A*x^2 + B*x + C = 0
    number_of_solutions = fgsl_poly_solve_cubic(A, B, C,   &
                                                roots(1),roots(2),roots(3))
    previous_distance = huge(previous_distance)
    do i=1,number_of_solutions
      distance = (roots(i)-refpoint)**2
      if (distance < previous_distance) then
        previous_distance = distance
        CubicRootClose = roots(i)
      endif
    end do
    return
  end function CubicRootClose
end module root_finder