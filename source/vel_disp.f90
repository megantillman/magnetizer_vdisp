module vel_disp
  use root_finder
  use grid
  
  implicit none
  private

  ! defining constants to be used in this calculation
  double precision, parameter :: epsff=0.005d0
  double precision, parameter :: alps=-1.5d0
  double precision, parameter :: phia=2.0d0
  
  double precision, parameter :: fsf=0.5d0
  double precision, parameter :: fgP=0.5d0
  double precision, parameter :: fgQ=0.5d0
  double precision, parameter :: Q=1.0d0
  
  double precision, parameter :: beta=0.0d0
  double precision, parameter :: eta=1.5d0
  double precision, parameter :: pmstar=3.0d3
  double precision, parameter :: sigma_sf_fix=0.0d0
  double precision, parameter :: phimp=1.4d0
  double precision, parameter :: phint=0.0d0
  double precision, parameter :: phiQ=2.0d0
  
  double precision, parameter :: pc=3.0856776d18
  double precision, parameter :: yr=365.25d0*24.0d0*3600.0d0
  double precision, parameter :: Myr=1.0d6*yr
  double precision, parameter :: kmps = 1.0d5
  double precision, parameter :: G=6.6743d-8
  double precision, parameter :: Msun=1.9891d33
  double precision, parameter :: pi = 3.14159265358979d0
  
  ! make sure this can be called in the actual code
  public :: calc_v_kms  
  
contains 
  subroutine calc_v_kms(Mgas_disk, r_disk, v_disk, SFR, p_ISM_sound_speed_km_s, v_kms)

    ! subroutine inputs from the main code
    double precision, intent(in) :: p_ISM_sound_speed_km_s
    double precision, intent(in) :: Mgas_disk, r_disk, v_disk, SFR 

    ! local variables
    double precision :: t_orb, t_disk, tsfmax
    double precision, dimension (nx) :: Sig_g, Sig_sf
    double precision, dimension(nx) :: sigma_g, sigma_sf
    double precision :: maxfac_g
    
    ! subroutine output, v_kms is an input to be consistent with the main code but is not affected by what is actually given to the subroutine
    double precision, dimension(nx), intent(out) :: v_kms
    
    ! calc orbital time in Myrs
    t_orb = (2.0d0*pi*(r_disk*1000.0d0*pc)/(v_disk*kmps))/Myr
    ! calc Sigma gas in Msun / pc^2
    Sig_g = Mgas_disk/(pi*(r_disk*1000.0d0)**2)
    ! calc t disk in seconds
    t_disk = (r_disk*1000.0d0*pc)/(v_disk*kmps)
    ! calc maximum star formation time scale in Myrs
    tsfmax = (t_disk*(v_disk/200.0d0)**alps/epsff)/Myr
    
    if ((sqrt(2.0d0*(1.0d0+beta)/(3.0d0*fgP*phimp))*8.0d0*epsff*fgQ/Q) .gt. (t_orb/tsfmax)) then
      maxfac_g = sqrt(2.0d0*(1.0d0+beta)/(3.0d0*fgP*phimp))*8.0d0*epsff*fgQ/Q
    else
      maxfac_g = t_orb/tsfmax
    endif
    
    ! calc the turbulent velocity for the non SF supported case
    sigma_g = ((SFR*Msun/yr)*(pi*G*Q)/sqrt(2.0d0/(1.0d0+beta))/maxfac_g/ (phia*fsf*fgQ*(v_disk*kmps)**2))/kmps
    
    ! calc the turbulent velocity for the SF supported case
    call calc_sig_sf(p_ISM_sound_speed_km_s, sigma_sf, t_orb, tsfmax)
    
    ! calc approximately the gas surface density to see which turb vel to use
    Sig_sf = sqrt(2.0d0*(1.0d0+beta)) * (2.0d0*pi/(t_orb*Myr)) * sigma_sf * kmps / (pi*G*(Q/fgQ)) / (Msun/pc**2)
    
    ! check if we can support via SF or through SF and mass transport
    if (Sig_g(nx) .gt. Sig_sf(nx) .and. sigma_g(1) .gt. sigma_sf(1)) then
      v_kms = sigma_g
    else
      v_kms = sigma_sf
    endif
    
    ! this prevents zeros which make a nan later and causes the code to fail
    ! this probably isn't needed anymore since I implemented the sig_g > sig_sf condition
    if (v_kms(1) .lt. 10.0d0) then
      v_kms = 10.0d0
    endif
    
    ! cap the velocity to somthing physical?
    if (v_kms(1) .gt. 200.0d0) then
      v_kms = 200.0d0
    endif
    
  end subroutine calc_v_kms
  
  
  
  
  subroutine calc_sig_sf(p_ISM_sound_speed_km_s, sigma_sf, t_orb, tsfmax)

    ! subroutine inputs
    double precision, intent(in) :: p_ISM_sound_speed_km_s
    double precision, intent(in) :: t_orb, tsfmax

    ! local variables
    double precision :: fac1, fac
    double precision :: a3, a2, a1, a0
    double precision :: refpoint
    double precision :: sigma_sf2
    double precision, dimension(3) :: roots

    !subroutine output
    double precision, dimension(nx), intent(out):: sigma_sf
    
    fac1=sqrt(3.0d0*fgP/(8.0d0*(1.0d0+beta)))*Q*phimp/(4.0d0*fgQ*epsff)*t_orb/tsfmax
    if (fac1 .lt. 1.0d0) then
      fac1=1.0d0
    endif
    fac=4.0d0*fac1*fsf*epsff*pmstar/(sqrt(3.0d0*fgP)*pi*eta*phimp*phiQ)
    
    a3=1.0d0
    a2=-(3.0d0*p_ISM_sound_speed_km_s**2+fac**2)
    a1=3.0d0*p_ISM_sound_speed_km_s**4
    a0=-p_ISM_sound_speed_km_s**6
    refpoint = 230.0d0
  
    sigma_sf2 = CubicRootClose(a3, a2, a1, a0, refpoint, roots)
    sigma_sf = sqrt(sigma_sf2)
    
  end subroutine calc_sig_sf




end module vel_disp
