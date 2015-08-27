! Contains the subroutines which compute initial profiles used in the dynamo calculations
module profiles
  use global_input_parameters
  use calc_params
  use grid
!
  implicit none
!
  double precision, dimension(nx) :: h, h_kpc
  double precision, dimension(nx) :: om, G, om_kmskpc, G_kmskpc
  double precision, dimension(nx) :: Uz, Uz_kms
  double precision, dimension(nx) :: Ur, dUrdr, d2Urdr2, Ur_kms
  double precision, dimension(nx) :: n, n_cm3
  double precision, dimension(nx) :: l, l_kpc, dldr
  double precision, dimension(nx) :: v, v_kms, dvdr
  double precision, dimension(nx) :: etat, etat_cm2s, etat_kmskpc
  double precision, dimension(nx) :: tau, tau_Gyr, tau_s
  double precision :: tau_sol, tau_sol_Gyr, tau_sol_s
  double precision, dimension(nx) :: Beq, Beq_mkG
  integer :: ialp_k
  double precision, dimension(nx) :: alp_k, alp_k_kms

  private :: disk_rotation_curve
  private :: bulge_rotation_curve
  private :: halo_rotation_curve

!
contains
  subroutine construct_profiles
!     SCALE HEIGHT PROFILE
    if (Flaring) then
      h= h_sol*dexp((r-r_sol)/r_h)
    else
      h= h_sol
    endif
    h_kpc=h*h0_kpc/h0

    ! ROTATION CURVE
    if (Om_Brandt) then
      om=om0/((1.d0+(r/r_om)**2))**(1.d0/2)  !Brandt angular velocity profile
      if (Shear) then
        G=-om0*(r/r_om)**2/(1.d0+(r/r_om)**2)**(3.d0/2)
      else
        G=0.
      endif
!       else if (three_components_rotation_curve)
!         call compute_rotation_curve(M_bulge,r_bulge,M_disk,r_disk,Mhalo,cs, om, G)
    endif
    om_kmskpc=om*h0_km/h0_kpc/t0_s*t0
    G_kmskpc=  G*h0_km/h0_kpc/t0_s*t0
!
!     VERTICAL VELOCITY PROFILE
    if (.not.Var_Uz) then
      Uz= Uz_sol  !No variation of Uz
    else
      Uz= Uz_sol*dexp(-(r-r_sol)/r_Uz)  !Decreasing with radius according to exponential
    endif
    Uz_kms=Uz*h0_km/h0/t0_s*t0
!
!     RADIAL VELOCITY PROFILE
    Ur=Ur_sol
    dUrdr=0.d0
    d2Urdr2=0.d0
    Ur_kms= Ur*h0_km/h0/t0_s*t0
!
!     NUMBER DENSITY PROFILE
    if (.not.Var_n) then
      n=n_sol
    else
      n=n_sol*exp(-(r-r_sol)/r_n)
    endif
    n_cm3=n*n0_cm3/n0
!
!     TURBULENT SCALE PROFILE
    if (.not.Var_l) then
      l=l_sol
      dldr=0.d0
    else
      l=l_sol*exp((r-r_sol)/r_l)
      dldr=l/r_l
    endif
    l_kpc=l*h0_kpc/h0
!
!     RMS TURBULENT VELOCITY PROFILE
    if (.not.Var_v) then
      v=v_sol
      dvdr=0.d0
    else
      v=v_sol*exp(-(r-r_sol)/r_v)
      dvdr=v/r_v
    endif
    v_kms=v*h0_km/h0/t0_s*t0
!
!     TURBULENT DIFFUSIVITY PROFILE
    etat=1.d0/3*l*v  !Formula for etat from mixing length theory
    etat_cm2s=etat*h0_cm**2/h0**2/t0_s*t0
    etat_kmskpc=etat*h0_km*h0_kpc/h0**2/t0_s*t0
!
!     TURBULENT CORRELATION TIME PROFILE
    tau=        ctau*l/v  !Formula for tau from mixing length theory
    tau_Gyr=    ctau*tau*t0_Gyr/t0
    tau_s=      ctau*tau*t0_s/t0
    tau_sol=    ctau*l_sol/v_sol
    tau_sol_Gyr=ctau*tau_sol*t0_Gyr/t0
    tau_sol_s=  ctau*tau_sol*t0_s/t0
!
!     EQUIPARTITION MAGNETIC FIELD STRENGTH PROFILE
    Beq=dsqrt(4*pi*n)*v  !Formula for equiparition field strength
    Beq_mkG=Beq*B0_mkG/B0
!
!     KINETIC ALPHA PROFILE
    if (.not.Krause) then
      alp_k= C_alp  !No variation of alpha
    else
      alp_k= C_alp*l**2/h*om  !Decreasing with radius
    endif
    if (Alp_ceiling) then
      do ialp_k=1,nx
        if (alp_k(ialp_k)>alpceil*v(ialp_k)) then
          alp_k(ialp_k)=alpceil*v(ialp_k)
        endif
      enddo
    endif
    alp_k_kms=alp_k*h0_km/h0/t0_s*t0
  end subroutine construct_profiles

  function disk_rotation_curve(r,r_disk,v_disk) result(v)
    ! Computes the rotation curve associated with an exponential disk
    ! Input: r -> the radii where the rotation curve will be computed
    !        r_disk ->  the half mass radius of the disk
    !        v_disk -> the circular velocity at r_disk
    !
    ! Info:  V^2 \propto y^2 [I0(y)K0(y)-I1(y)K1(y)] where y=r/r_s
    !        r_s --> the scale radius
    ! Ref:   Binney & Tremaine or  Mo, Bosch & White

    use Bessel_Functions
    implicit none
    double precision, intent(in) :: r_disk, v_disk
    double precision, dimension(:), intent(in)  :: r
    double precision, dimension(size(r)) :: v
    double precision, parameter :: rs_to_r50 = 1.678346990d0
    double precision, dimension(size(r)) :: y
    integer :: i

    y = r / (r_disk/rs_to_r50)
    do i=1,nx
      v(i) = ( v_disk * y(i)/rs_to_r50 )**2                                   &
             * ( Bessel_Function_I0(y(i)) * Bessel_Function_K0(y(i))          &
               - Bessel_Function_I1(y(i)) * Bessel_Function_K1(y(i)) )        &
             / ( Bessel_Function_I0(rs_to_r50) * Bessel_Function_K0(rs_to_r50)&
               - Bessel_Function_I1(rs_to_r50) * Bessel_Function_K1(rs_to_r50))
      v(i) = sqrt(v(i))
    end do
    return
  end function disk_rotation_curve

  function bulge_rotation_curve(r,r_bulge,v_bulge) result(v)
    ! Computes the rotation curve associated with an Hernquist profile
    ! Input: r -> the radii where the rotation curve will be computed
    !        r_bulge ->  the half mass radius of the spheroid
    !        v_bulge -> the circular velocity at r_bulge
    !
    ! Info:  V^2 \propto y/(y+1)^2  where y=r/a
    ! Ref:   http://adsabs.harvard.edu/abs/1990ApJ...356..359H

    implicit none
    double precision, intent(in) :: r_bulge, v_bulge
    double precision, dimension(:), intent(in)  :: r
    double precision, dimension(size(r)) :: v
    double precision, parameter :: a_to_r50 = 1.0d0/(sqrt(2.0d0)-1.0d0)
    double precision, dimension(size(r)) :: y
    integer :: i

    y = r / (r_bulge/a_to_r50)

    v =  (y/a_to_r50) * (a_to_r50 +1d0)**2 * (y +1d0)**(-1)
    v = v_bulge * sqrt(v)
    return
  end function bulge_rotation_curve

  function halo_rotation_curve(r,r_halo, v_halo, cs1) result(v)
    ! Computes the rotation curve associated with an NFW halo
    ! Input: r -> the radii where the rotation curve will be computed
    !        r_halo -> the virial radius of the halo
    !        v_halo -> the circular velocity at r_halo
    !        cs1 -> 1/c_s -> inverse of the NFW concentration parameter
    !                        i.e. (NFW scale radius)/(virial radius)
    ! Warning: This ignores effects of adiabatic contraction!
    !
    ! Info:  V^2 \propto {ln[(cs1+y)/cs1] - y/(cs1+y)} /
    !                    {ln[(cs1+1)/cs1] - 1/(cs1+1)} / y
    ! Ref: NFW profile

    implicit none
    double precision, intent(in) :: r_halo, v_halo, cs1
    double precision, dimension(:), intent(in)  :: r
    double precision, dimension(size(r)) :: v
    double precision, dimension(size(r)) :: y
    integer :: i

    y = r / r_halo

    v = (log((cs1+y)/cs1) - y/(cs1+y)) / &
                       (log((cs1+1d0)/cs1) - 1d0/(cs1+1d0)) / y
    v = v_halo * sqrt(v)
    return
  end function halo_rotation_curve

end module profiles
