module vel_disp
  use root_finder
  use grid
  
  implicit none
  private

  ! defining constants to be used in this calculation
  double precision, parameter :: epsff=0.005
  double precision, parameter :: alps=-1.5
  double precision, parameter :: phia=2.0
  
  double precision, parameter :: fsf=0.5
  double precision, parameter :: fgP=0.5
  double precision, parameter :: fgQ=0.5
  double precision, parameter :: Q=1.0
  
  double precision, parameter :: beta=0.0
  double precision, parameter :: eta=1.5
  double precision, parameter :: pmstar=3d3
  double precision, parameter :: sigma_sf_fix=0.0
  double precision, parameter :: phimp=1.4
  double precision :: phint=0.0
  double precision, parameter :: phiQ=2.0
  
  double precision, parameter :: pc=3.0856776d18
  double precision, parameter :: yr=365.25*24.*3600.
  double precision, parameter :: Myr=1e6*yr
  double precision, parameter :: kmps = 1d5
  double precision, parameter :: G=6.6743d-8
  double precision, parameter :: Msun=1.9891d33
  
  ! make sure this can be called in the actual code
  public :: calc_v_kms  
  
contains 
  subroutine calc_v_kms(Mgas_disk, r_disk, v_disk, SFR, p_ISM_sound_speed_km_s, v_kms)
  
    double precision, intent(in) :: p_ISM_sound_speed_km_s
    double precision, intent(in) :: Mgas_disk, r_disk, v_disk, SFR 
    double precision :: t_orb, t_disk, tsfmax
    double precision, dimension (nx) :: Sig_g, Sig_sf
    double precision, dimension(nx), intent(inout) :: v_kms
    
    
    double precision, dimension(nx) :: sigma_g, sigma_sf
    double precision ::  maxfac_g
    
    
    ! calc orbital time in Myrs
    t_orb = (2*3.14*(r_disk*1000.*pc)/(v_disk*kmps))/Myr
    ! calc Sigma gas in Msun / pc^2
    Sig_g = Mgas_disk/(3.14*(r_disk*1000.)**2)
    ! calc t disk in seconds
    t_disk = (r_disk*1000*pc)/(v_disk*kmps)
    ! calc maximum star formation time scale in Myrs
    tsfmax = (t_disk*(v_disk/200.)**alps/epsff)/Myr
    
    if ((sqrt(2.0*(1.0+beta)/(3.0*fgP*phimp))*8.0*epsff*fgQ/Q) .gt. (t_orb/tsfmax)) then
      maxfac_g = sqrt(2.0*(1.0+beta)/(3.0*fgP*phimp))*8.0*epsff*fgQ/Q
    else
      maxfac_g = t_orb/tsfmax
    endif
    
    ! calc the turbulent velocity for the non SF supported case
    sigma_g = ((SFR*Msun/yr)*(3.14*G*Q)/sqrt(2.0/(1.0+beta))/maxfac_g/ (phia*fsf*fgQ*(v_disk*kmps)**2))/kmps
    
    ! calc the turbulent velocity for the SF supported case
    call calc_sig_sf(p_ISM_sound_speed_km_s, sigma_sf, t_orb, tsfmax)
    
    ! calc approximately the gas surface density to see which turb vel to use
    Sig_sf = sqrt(2.0*(1.0+beta)) * (2.0*3.14/(t_orb*Myr)) * sigma_sf * kmps / (3.14*G*(Q/fgQ)) / (Msun/pc**2)
    
    ! check if we can support via SF or through SF and mass transport
    if (Sig_g(nx) .gt. Sig_sf(nx) .and. sigma_g(1) .gt. sigma_sf(1)) then
      v_kms = sigma_g
    else
      v_kms = sigma_sf
    endif
    
    ! this prevents zeros which make a nan later and causes the code to fail
    ! this probably isn't needed anymore since I implemented the sig_g > sig_sf condition
    if (v_kms(1) .lt. 10.) then
      v_kms = v_kms*0.+10.
    endif
    
    ! cap the velocity to somthing physical?
    if (v_kms(1) .gt. 200.) then
      v_kms = v_kms*0.+200.
    endif
    
  end subroutine calc_v_kms
  
  
  
  
  subroutine calc_sig_sf(p_ISM_sound_speed_km_s, sigma_sf, t_orb, tsfmax)
    
    double precision :: fac1, fac
    double precision :: a3, a2, a1, a0
    double precision :: refpoint
    double precision :: sigma_sf2
    double precision, dimension(3) :: roots
    double precision, intent(in) :: p_ISM_sound_speed_km_s
    double precision, intent(in) :: t_orb, tsfmax
    
    double precision, dimension(nx), intent(out):: sigma_sf
    
    fac1=sqrt(3.0*fgP/(8.0*(1.0+beta)))*Q*phimp/(4*fgQ*epsff)*t_orb/tsfmax
    if (fac1 .lt. 1.0) then
      fac1=1.0
    endif
    fac=4.0*fac1*fsf*epsff*pmstar/(sqrt(3.0*fgP)*3.14*eta*phimp*phiQ)
    
    a3=1.0
    a2=-(3.0*p_ISM_sound_speed_km_s**2+fac**2)
    a1=3.0*p_ISM_sound_speed_km_s**4
    a0=-p_ISM_sound_speed_km_s**6
    refpoint = 230.
  
    sigma_sf2 = CubicRootClose(a3, a2, a1, a0, refpoint, roots)
    sigma_sf = real(sqrt(sigma_sf2))
    
  end subroutine calc_sig_sf




end module vel_disp
