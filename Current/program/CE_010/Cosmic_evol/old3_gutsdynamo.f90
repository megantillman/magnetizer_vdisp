!*****************************************************
!This code contains modules, settings and parameters used by the galactic dynamo code dynamo.f90.
!*****************************************************
module direc_names  !Specifies paths
  implicit none
!
  character (len=6), parameter :: s0= 'CE_010'
  character (len=42), parameter :: s1= '/Users/luke/fortran_pde/1D/telegraph/r_noz'
end module direc_names
!*****************************************************
module math_constants  !Contains math constants, unit conversions
  implicit none
! CONSTANTS
  double precision, parameter :: pi=     3.14156295358d0
  double precision, parameter :: mp_g=   1.67262158d-24  !proton mass in grams
  double precision, parameter :: cm_kpc= 3.08567758d21  !number of cm in 1 kpc
  double precision, parameter :: km_kpc= 3.08567758d16  !number of km in 1 kpc
  double precision, parameter :: cm_km=  1.d5  !number of cm in 1 km
  double precision, parameter :: s_Gyr=  1.d9*365.25d0*24*3600  !number of s in 1 Gyr
  double precision, parameter :: mkG_G=  1.d6  !number of mkG in 1 G
end module math_constants
!*****************************************************
module input_constants  !Contains math constants, unit conversions
  use math_constants
  implicit none
!   TURBULENCE
    double precision, parameter :: ctau= 1.d0 !Ratio of tau (eddy turnover time) to correlation time of the turbulence
!   ADVECTION NO-z
    double precision, parameter :: C_U=0.25d0  !No-z correction term for terms involving Uz*Br, Uz*Bp
    double precision, parameter :: C_a=1.d0!0.3d0  !No-z correction term for terms involving Uz*alp_m
!   SEED FIELD
    double precision, parameter :: r1_kpc= 15.0d0 !Only relevant if Rand_seed=F
    double precision, parameter :: frac_seed= 0.0001d0 !Seed field amplitude as a fraction of equipartition magnetic field strength
    integer, parameter :: nn=2 !Only relevant if Rand_seed=F
!   DIFFUSIVE MAGNETIC HELICITY FLUX
    double precision, parameter :: C_d= -pi**2  !No-z correction term for terms involving diffusive flux
!   RESISTIVITY
    double precision, parameter :: Rm_inv= 0.d0  !1.e-5 !Inverse magnetic Reynolds number
!   ALPHA EFFECT
    double precision, parameter :: C_alp=   1.0d0 !Factor determining strength of alpha effect
    double precision, parameter :: alpceil= 0.5d0  !Relevant only if module Alp_ceiling=1; Limits maximum value of alpha to ceiling*v_kms
!   RANDOM MAGNETIC FIELD
    double precision, parameter :: fmag= 0.5d0 !Strength of the rms random magnetic field brms in units of Beq
end module input_constants
!*****************************************************
module modules  !Contains switches
  implicit none
!
! PARAMETER INPUTS
  logical, parameter :: Read_param=   .false.  !Set to T to read in parameters from file, set to F to use defaults
! PARAMETER INPUTS
  logical, parameter :: Time_evol=    .true.  !Set to T to read in parameters at several timesteps; set to F to read in parameters at the start only
! ALPHA QUENCHING
  logical, parameter :: Dyn_quench=   .false.  !Works with alg_quench=F; Set to T for dynamical quenching (d_alpha_m/dt eqn incl in sim)
  logical, parameter :: Alg_quench=   .false. !Works with dyn_quench=F; Set to T for algebraic alpha quenching; F for no quenching
! RANDOM MAGNETIC FIELD
  logical, parameter :: Floor=        .true.  !Strength of the rms random magnetic field brms in units of Beq
! CLOSURE APPROXIMATION
  logical, parameter :: Damp=         .false. !Set to T for FOSA, F for minimal tau approximation
! OUTFLOW VELOCITY
  logical, parameter :: Var_Uz=       .false.  !Set to T to make Uz decrease with radius; set to F to keep Uz constant with radius
! TURBULENT SCALE
  logical, parameter :: Var_l=        .false. !Set to T to make l decrease with radius; set to F to keep l constant with radius
! TURBULENT VELOCITY
  logical, parameter :: Var_v=        .false. !Set to T to make v decrease with radius; set to F to keep v constant with radius
! DENSITY
  logical, parameter :: Var_n=        .true.  !Set to T to make n decrease with radius; set to F to keep n constant with radius (Beq^2=4pi*n*mp*v^2)
! SCALE HEIGHT
  logical, parameter :: Flaring=      .true.  !Set to T to include variable scale height. 
! SEED MAGNETIC FIELD
  logical, parameter :: Rand_seed=    .false.  !Set to F for random seed magnetic field, T for some other seed
! CEILING ON ALPHA EFFECT
  logical, parameter :: Alp_ceiling=  .true. !Set to T to put a ceiling for alpha at alpceil*v
! BRANDT ROTATION CURVE
  logical, parameter :: Om_Brandt =   .true.  !Set to T to use a Brandt rotation curve
! ALPHA^2 EFFECT
  logical, parameter :: Alp_squared=  .true.  !Set to T to include alpha^2 effect; set to F to use alpha-omega approximation equations
! KRAUSE'S LAW
  logical, parameter :: Krause=       .true.  !Set to T for alpha effect to decrease in proportion to omega (Krause's formula)
! OMEGA EFFECT
  logical, parameter :: Shear=        .true.  !Set to T to include Omega effect in dynamo, 0 for alpha^2 dynamo
! ADVECTION
  logical, parameter :: Advect=       .true.  !Set to F to turn off advection of the magnetic field
! TURBULENT DIFFUSION
  logical, parameter :: Turb_dif=     .true.  !Set to F to turn off turbulent diffusion
end module modules
!*****************************************************
module units  !Specifies code units
  use math_constants
!
  implicit none
!
! DIMENSIONLESS UNITS
  double precision, parameter :: etat0=  1.0d0
  double precision, parameter :: h0=     1.0d0
  double precision, parameter :: r_disk= 1.0d0  !Note that radii have different length units than other quantities
  double precision, parameter :: B0=     1.0d0
  double precision, parameter :: t0=     h0**2/etat0
  double precision, parameter :: n0=     B0**2/h0**2*t0**2
! DIMENSIONAL UNITS
  double precision, parameter :: r_sol_kpc=    8.0d0  !Chosen radius at which to specify param values (e.g. radius of solar neighbourhood)
  double precision, parameter :: etat0_kmskpc= 1.d0/3*0.1d0*10.d0  !Unit of etat, corresp to a typical value, in km*kpc/s
  double precision, parameter :: etat0_cm2s=   etat0_kmskpc*cm_km*cm_kpc  !Unit of etat, corresp to a typical value, in cm^2/s
  double precision, parameter :: h0_kpc=       0.5d0  !Unit of length, corresp to a typical half-disk thickness, in kpc
  double precision, parameter :: h0_km=h0_kpc*km_kpc  !Unit of length, corresp to a typical half-disk thickness, in km
  double precision, parameter :: h0_cm=h0_km*cm_km  !Unit of length, corresp to a typical half-disk thickness, in cm
  double precision, parameter :: B0_mkG=       1.0d0  !Unit of magnetic field, in mkG
  double precision, parameter :: t0_kpcskm= h0_kpc**2/etat0_kmskpc  !Unit of time, corresp to a typical vert diff time, in s*kpc/km
  double precision, parameter :: t0_Gyr=t0_kpcskm/s_Gyr*km_kpc  !Unit of time in Gyr
  double precision, parameter :: t0_s=t0_kpcskm*km_kpc  !Unit of time in s
  double precision, parameter :: n0_cm3= (B0_mkG/mkG_G)**2/h0_cm**2*t0_s**2/mp_g  !Unit of number density, in cm^-3
end module units
!*****************************************************
module grid  !Contains grid parameters and subroutine for constructing grid
  use math_constants
  use units
!
  implicit none
!
  integer, parameter :: nxphys= 201!61!91  !Resolution in r (excluding ghost zones) (for convenience should be N*r_disk_kpc+1)
  integer, parameter :: nxghost= 3  !Number of ghost cells at each end in r
  integer, parameter :: nx= nxphys +2*nxghost  !Resolution in r
  double precision, parameter :: len_phi= 2*pi, r_in=0.0d0  !phi domain, min radius of disk, max radius of disk
  double precision :: dx
  double precision, dimension(nx) :: x
  double precision, dimension(nx) :: r
!
  contains
    subroutine construct_grid
  !
      integer :: i
  !
      dx= (r_disk-r_in)/(nxphys-1)  !x corresponds to r coordinate
      do i=1,nx
        x(i)=r_in -nxghost*dx +(dble(i)-1.)*dx
      enddo
      r=x !use explicit name
    endsubroutine construct_grid
end module grid
!*****************************************************
module input_params  !Contains default input parameters and a subroutine for reading input parameters from files
  use modules
  use math_constants
  use direc_names
!
  implicit none
!
  character (len=200) :: header_time_indep, header_time_dep
  integer :: iread=0
  double precision :: r_disk_kpc,R_kappa
  double precision :: l_sol_kpc, r_l_kpc, v_sol_kms, r_v_kpc, n_sol_cm3, r_n_kpc, Uz_sol_kms, r_Uz_kpc, &
                      h_sol_kpc, r_h_kpc, Uphi_sol_kms, r_om_kpc
  contains
    subroutine set_input_params(gal_id_string,info)
!
      character (len=8), intent(in) :: gal_id_string
      integer, intent(in) :: info
!
!
!     DIMENSIONAL INPUT PARAMETERS THAT MUST BE SPECIFIED AND MAY VARY FROM GALAXY TO GALAXY AND FROM TIMESTEP TO TIMESTEP
      if (Read_param) then
        open(29,file= s1 // '/input/' // s0 // '/time_indep_params_' &
                      // gal_id_string // '.in',status="old")
        read(29,*)header_time_indep  !Header gives time-independent parameters in order
        read(29,*)r_disk_kpc  !Read time-independent parameter values
        close(29)  !Close file containing time-independent parameter values
        open(30,file= s1 // '/input/' // s0 // '/time_dep_params_'   &
                      // gal_id_string // '.in',status="old")
        if (iread == 0) then
          read(30,*)header_time_dep  !Read header if first iteration only
        endif
        read(30,*)l_sol_kpc,r_l_kpc,v_sol_kms,r_v_kpc,n_sol_cm3,r_n_kpc,Uz_sol_kms,r_Uz_kpc, &
                  h_sol_kpc,r_h_kpc,Uphi_sol_kms,r_om_kpc  !Read time-dependent parameter values
        iread=iread+1
        if (info> 0) then
          print*,''
          print*,'Reading new parameters, iread=',iread, ', gal_id_string=',gal_id_string
        endif
      else
!       NUMERICAL (DEFAULTS)
        r_disk_kpc=     15.0d0  !The maximum radius of the domain in kpc (unit of length along r)
!
!       TURBULENCE (DEFAULTS)
        l_sol_kpc=       0.1d0  !Typical size of largest turbulent eddies, in kpc
        r_l_kpc=         7.5d0  !Exponential scale radius of l in kpc; Relevant only if Var_l=T
        v_sol_kms=      10.0d0  !Typical turbulent velocity in km/s; together with h0_kpc it determ the unit of time
        r_v_kpc=         7.5d0  !Exponential scale radius of l in kpc; Relevant only if Var_v=T
!
!       EQUIPARTITION FIELD (DEFAULTS)
        n_sol_cm3=       1.0d0  !Number density at r=r_sol in cm^-3
        r_n_kpc=         7.5d0  !Exponential scale radius of density in kpc n=n_sol*exp((r-r_sol)/r_n), note that r_Beq=2*r_n if r_v=infty
!
!       VERTICAL WIND (DEFAULTS)
        Uz_sol_kms=     5.0d0  !Vertical mean velocity at r=r_sol in km/s (if Var_Uz=F then Uz_kms=Uz_sol_kms)
        r_Uz_kpc=       15.0d0 !Uz_kms=Uz_sol_kms*exp((r_kpc-r_sol_kpc)/r_Uz_kpc) if Var_Uz=T
!
!       SCALE HEIGHT (DEFAULTS)
        h_sol_kpc=       0.5d0 !amplitude of h_kpc=h_sol_kpc*exp((r_kpc-r_sol_kpc)/r_h_kpc)
        r_h_kpc=         9.8d0 !Relevant only if Flaring=1; charac radius of the hyperb/exp varying scale height in kpc
!
!       ROTATION CURVE (DEFAULTS)
        Uphi_sol_kms=  220.0d0  !Circular rotation speed at r=r_sol in km/s
        r_om_kpc=        2.0d0  !Relevant only if Om_Brandt =1; r_om is the charac radius of Brandt profile in kpc
      endif
!
!     DIMENSIONLESS PARAMETERS THAT MUST BE SPECIFIED BUT WILL NOT NORMALLY VARY FROM GALAXY TO GALAXY
!     DIFFUSIVE MAGNETIC HELICITY FLUX (DEFAULTS)
      R_kappa=         1.0d0 !Ratio kappa_t/eta_t of turbulent diffusivities of alpha_m and B
!
!      print*,'r_disk_kpc= ',r_disk_kpc ,'R_kappa=   ',R_kappa
!      print*,'l_sol_kpc=  ',l_sol_kpc  ,'r_l_kpc=   ',r_l_kpc   ,'v_sol_kms=    ',v_sol_kms    ,'r_v_kpc=    ',r_v_kpc
!      print*,'n_sol_cm3=  ',n_sol_cm3  ,'r_n_kpc=   ',r_n_kpc   ,'Uz_sol_kms=   ',Uz_sol_kms   ,'r_Uz_kpc=   ',r_Uz_kpc
!      print*,'h_sol_kpc=  ',h_sol_kpc  ,'r_h_kpc=   ',r_h_kpc   ,'Uphi_sol_kms= ',Uphi_sol_kms ,'r_om_kpc=   ',r_om_kpc
!      print*,''
    endsubroutine set_input_params
end module input_params
!*****************************************************
module calc_params  !Contains parameters that are calculated from the input parameters
  use units
  use math_constants
  use input_params
  use grid
  use input_constants
!
  implicit none
!
  double precision :: etat_sol,etat_sol_kmskpc,etat_sol_cm2s,td_sol,td_sol_kpcskm,td_sol_Gyr,td_sol_s,om0_kmskpc, &
                      Ur_sol_kms,r_sol,lambda,n_sol,r_n,h_sol,r_h,l_sol,r_l,v_sol,r_v,Uz_sol,r_Uz, &
                      Ur_sol,Uphi_sol,om0,r_om,r1
  double precision, dimension(nx) :: r_kpc
!
  contains
    subroutine set_calc_params
!     DIMENSIONAL PARAMETERS THAT CAN BE CALCULATED FROM SPECIFIED DIMENSIONAL PARAMETERS OR THAT ARE NOT NORMALLY VARIED:
!     TURBULENCE
      etat_sol_kmskpc= 1.d0/3*l_sol_kpc*v_sol_kms  !Typical value of etat in units of km kpc/s
      etat_sol_cm2s= 1.d0/3*l_sol_kpc*cm_kpc*v_sol_kms*cm_km  !Typical turbulent diffusivity in units of cm^2/s
      td_sol_kpcskm= h0_kpc**2/etat_sol_kmskpc  !Typical vertical turbulent diffusion timescale in units of kpc s/km
      td_sol_Gyr= h0_kpc**2/etat_sol_cm2s/s_Gyr*cm_kpc*cm_kpc  !Typical vertical turbulent diffusion timescale in units of Gyr
      td_sol_s= td_sol_Gyr*s_Gyr  !Typical vertical turbulent diffusion timescale in units of seconds
!
!     ROTATION CURVE
      om0_kmskpc= dsqrt(1.d0+(r_sol_kpc/r_om_kpc)**2)*Uphi_sol_kms/r_sol_kpc  !om0 of Brandt profile in physical units
!
!     RADIAL FLOW
      Ur_sol_kms=   0.0d0  !Radial mean velocity at r=r_sol in km/s
!
!     DIMENSIONLESS PARAMETERS THAT CAN BE CALCULATED OR THAT ARE NOT NORMALLY VARIED:
!     NUMERICAL
      r_sol=r_sol_kpc/r_disk_kpc !Chosen radius at which to specify param values (e.g. radius of solar neighbourhood)
      lambda=h0_kpc/r_disk_kpc  !Typical aspect ratio of disk
!
!     TURBULENCE
      h_sol= h_sol_kpc/h0_kpc*h0  !Disk thickness at r=r_sol in units of h0
      r_h= r_h_kpc/r_disk_kpc*r_disk  !Exponential scale radius of disk scale height
      n_sol= n_sol_cm3/n0_cm3*n0  !Equipartition magnetic field at r=r_sol
      r_n= r_n_kpc/r_disk_kpc*r_disk  !Exponential scale radius of equipartition magnetic field
      l_sol= l_sol_kpc/h0_kpc*h0  !Size of largest turbulent eddies
      r_l= r_l_kpc/r_disk_kpc*r_disk  !Exponential scale radius of turbulent scale
      v_sol=v_sol_kms/h0_km*t0_s*h0/t0  !Turbulent velocity
      r_v= r_v_kpc/r_disk_kpc*r_disk  !Exponential scale radius of rms turbulent velocity
      etat_sol=1.d0/3*l_sol*v_sol  !Typical turbulent diffusivity
      td_sol= h_sol**2/etat_sol  !Typical vertical turbulent diffusion timescale
!
!     VERTICAL WIND
      Uz_sol=Uz_sol_kms/h0_km*t0_s*h0/t0  !Vertical mean velocity
      r_Uz= r_Uz_kpc/r_disk_kpc*r_disk  !Exponential scale radius of vertical mean velocity
!
!     RADIAL FLOW
      Ur_sol=Ur_sol_kms/h0_km*t0_s*h0/t0  !Radial mean velocity at r=r_sol
!
!     ROTATION CURVE
      Uphi_sol=Uphi_sol_kms/h0_km*t0_s*h0/t0  !Dimensionless circular velocity at r=r_sol
      om0=om0_kmskpc/km_kpc*t0_s/t0  !om0 of Brandt profile
      r_om=r_om_kpc/r_disk_kpc*r_disk  !Characteristic radius of Brandt profile
!
!     SEED FIELD
      r1= r1_kpc/r_disk_kpc*r_disk !Only relevant if Rand_seed=F
!
!     PHYSICAL GRID
      r_kpc=r*r_disk_kpc
  endsubroutine set_calc_params
end module calc_params
!*****************************************************
module ts_params  !Contains time-stepping parameters
  use input_params
  use calc_params
!
  implicit none
!
  integer, parameter :: n1= 50!420   !Number of snapshots
  integer, parameter :: nread= 99999  !Read in new input parameters every nread timesteps
  integer, parameter :: nscreen= 100  !Print output to screen every nscreen*n2 timesteps
  double precision, parameter :: tsnap= 0.025d0/t0_Gyr !Time between successive snapshots
  integer :: n2  !Number of timesteps in between snapshots
  double precision :: dt,t=0.d0,first=0.d0, eps_t=0.005d0
!
  contains
    subroutine set_ts_params
      n2= nint(tsnap/t0_Gyr/eps_t)  !Change n2 by changing eps_t
      dt= tsnap/n2  !Timestep in units of t0=h0^2/etat0
    endsubroutine set_ts_params
end module ts_params
!*****************************************************
module var  !Contains subroutine for setting the number of variables for which to solve
  use modules
!
  implicit none
  integer :: nvar
!
  contains
    subroutine init_var
!     DEFINE ARRAYS FOR PHYSICAL VARIABLES
      if (.not.Dyn_quench) then
        if (.not.Damp) then
          nvar=2  !vars are Br, Bp
        else
          nvar=4  !vars are Br, Bp, Fr, Fp
        endif
      else
        if (.not.Damp) then
          nvar=3  !vars are Br, Bp, alp_m
        else
          nvar=7  !vars are Br, Bp, Fr, Fp, Er, Ep, alp_m
        endif
      endif
    end subroutine init_var
end module var
!*****************************************************
module ts_arrays  !Contains subroutine that stores time series data (n1 snapshots, separated by n2 timesteps)
!
  use modules
  use var
  use grid
  use ts_params
!
  implicit none
!
  double precision, dimension(n1+1) :: ts_t, ts_rmax, ts_delta_r  !time array
  double precision, dimension(n1+1,nx) :: ts_Br, ts_Bp, ts_alp_m, ts_Bzmod, ts_h, ts_om, ts_G, ts_l, ts_v, &
                                          ts_etat, ts_tau, ts_alp_k, ts_Uz, ts_Ur, ts_n, ts_Beq
!
  contains
    subroutine make_ts_arrays(it,t,f,Bzmod,h,om,G,l,v,etat,tau,alp_k,Uz,Ur,n,Beq,rmax,delta_r)
!
      integer, intent(in) :: it
      double precision, intent(in) :: t
      double precision, intent(in) :: rmax
      double precision, intent(in) :: delta_r
      double precision, dimension(nx,nvar), intent(in) :: f
      double precision, dimension(nx), intent(in) :: Bzmod
      double precision, dimension(nx), intent(in) :: h
      double precision, dimension(nx), intent(in) :: om
      double precision, dimension(nx), intent(in) :: G
      double precision, dimension(nx), intent(in) :: l
      double precision, dimension(nx), intent(in) :: v
      double precision, dimension(nx), intent(in) :: etat
      double precision, dimension(nx), intent(in) :: tau
      double precision, dimension(nx), intent(in) :: alp_k
      double precision, dimension(nx), intent(in) :: Uz
      double precision, dimension(nx), intent(in) :: Ur
      double precision, dimension(nx), intent(in) :: n
      double precision, dimension(nx), intent(in) :: Beq
!
      ts_t(it+1)=       t
      ts_rmax(it+1)=    rmax
      ts_delta_r(it+1)= delta_r
      ts_Br(it+1,:)=    f(:,1)
      ts_Bp(it+1,:)=    f(:,2)
      if (Dyn_quench) then   
        if (.not.Damp) then
          ts_alp_m(it+1,:)=   f(:,3)
        else
          ts_alp_m(it+1,:)=   f(:,7)
        endif
      endif
      ts_h(it+1,:)=     h(:)
      ts_om(it+1,:)=    om(:)
      ts_G(it+1,:)=     G(:)
      ts_l(it+1,:)=     l(:)
      ts_v(it+1,:)=     v(:)
      ts_etat(it+1,:)=  etat(:)
      ts_tau(it+1,:)=   tau(:)
      ts_alp_k(it+1,:)= alp_k(:)
      ts_Uz(it+1,:)=    Uz(:)
      ts_Ur(it+1,:)=    Ur(:)
      ts_n(it+1,:)=     n(:)
      ts_Beq(it+1,:)=   Beq(:)
      ts_Bzmod(it+1,:)= Bzmod(:)
    end subroutine make_ts_arrays
end module ts_arrays
!*****************************************************
module profiles
  use modules
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
!
!     ROTATION CURVE
      if (Om_Brandt) then
        om=om0/((1.d0+(r/r_om)**2))**(1.d0/2)  !Brandt angular velocity profile
        if (Shear) then
          G=-om0*(r/r_om)**2/(1.d0+(r/r_om)**2)**(3.d0/2)
        else
          G=0.
        endif
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
end module profiles
!*****************************************************
module boundary_conditions  !Specify and implement boundary conditions at r=0, r=R
  use modules
  use grid
  use var
!  
  implicit none
!
  contains
    subroutine impose_bc(f)
!
      double precision, dimension(nx,nvar), intent(inout) :: f
      integer :: ix
!
      if (nxghost/=0) then
        f(nxghost+1 ,1:2)=  0.d0  !Set Br=Bp=0 at r=0 	Dirichlet BC on Br, Bp
        f(nx-nxghost,1:2)=  0.d0  !Set Br=Bp=0 at r=R 	Dirichlet BC on Br, Bp
        do ix=1,nxghost
          f(ix     ,1:2)= -f(2*(nxghost+1)-ix     ,1:2)  !Antisymmetric about r=0    Dirichlet BC on Br, Bp: Br=Bp=0 at r=0
          f(nx+1-ix,1:2)= -f(nx+1-2*(nxghost+1)+ix,1:2)  !Antisymmetric about r=R    Dirichlet BC on Br, Bp: Br=Bp=0 at r=R
          if (Dyn_quench) then
            f(ix     ,nvar  )=  f(2*(nxghost+1)-ix     ,nvar  )  !Symmetric     about r=0    Neumann   BC on alp_m    : dalp_mdr=0 at r=0
            f(nx+1-ix,nvar  )=  f(nx+1-2*(nxghost+1)+ix,nvar  )  !Symmetric     about r=R    Neumann   BC on alp_m    : dalp_mdr=0 at r=R
          endif
        enddo
      endif
    end subroutine impose_bc
end module boundary_conditions  
!*****************************************************
module initial_conditions  !Set initial conditions
  use calc_params
  use grid
  use var
  use boundary_conditions
  use random
  use profiles
!  
  implicit none
!
  contains
    subroutine init_seed(f)
!
      integer :: iseed,var
      double precision, dimension(nx) :: Bseed
      double precision, dimension(nx,nvar), intent(inout) :: f
!
      Bseed=frac_seed*Beq
!
      f(:,:)=0.d0 !Initialize
!
      if (Rand_seed) then  
        do var=1,2 !Seed r and phi components of B
        do iseed=1,nx
          f(iseed,var)=Bseed(iseed)*random_normal()
        enddo
        enddo
      else
        f(:,1)=-Bseed*r/(r_disk-dx)*(1.d0-r/(r_disk-dx))**nn*dexp(-r/r1)
        f(:,2)= Bseed*r/(r_disk-dx)*(1.d0-r/(r_disk-dx))**nn*dexp(-r/r1)
      endif
      call impose_bc(f)
!      print*,'Seed field: Br        =',f(:,1)
!      print*,'            Bphi      =',f(:,2)
!      print*,'Seed field: Br   (mkG)=',f(:,1)*B0_mkG/B0
!      print*,'            Bphi (mkG)=',f(:,2)*B0_mkG/B0
    end subroutine init_seed
end module initial_conditions
!*****************************************************
module deriv  !Finite differencing routine for spatial differentiation (same as PENCIL Code)
  use grid
!
  contains
    function xder(f)
!
      implicit none
!
      integer :: ix
      double precision, dimension(nx), intent(in) :: f
      double precision :: fac
      double precision, dimension(nx) :: xder
!
!     assume uniform mesh
!
      fac=1.d0/(60.*(x(2)-x(1)))
!
      do ix=4,nx-3
        xder(ix)=fac*(+45.*(f(ix+1)-f(ix-1)) & 
                      - 9.*(f(ix+2)-f(ix-2)) &
                      +    (f(ix+3)-f(ix-3)) &
                     )
      enddo
!
      xder(1)=fac*(-147.d0*f(1)+360.d0*f(2)-450.d0*f(3)+400.d0*f(4)-225.d0*f(5)+72.d0*f(6)-10.d0*f(7))
      xder(2)=fac*( -10.d0*f(1)- 77.d0*f(2)+150.d0*f(3)-100.d0*f(4)+ 50.d0*f(5)-15.d0*f(6)+ 2.d0*f(7))
      xder(3)=fac*(   2.d0*f(1)- 24.d0*f(2)- 35.d0*f(3)+ 80.d0*f(4)- 30.d0*f(5)+ 8.d0*f(6)-      f(7))
!
!     outer points
!
      xder(nx  )=fac*(147.d0*f(nx)-360.d0*f(nx-1)+450.d0*f(nx-2)-400.d0*f(nx-3)+225.d0*f(nx-4)-72.d0*f(nx-5) &
                      +10.d0*f(nx-6))
      xder(nx-1)=fac*( 10.d0*f(nx)+ 77.d0*f(nx-1)-150.d0*f(nx-2)+100.d0*f(nx-3)- 50.d0*f(nx-4)+15.d0*f(nx-5) &
                      - 2.d0*f(nx-6))
      xder(nx-2)=fac*( -2.d0*f(nx)+ 24.d0*f(nx-1)+ 35.d0*f(nx-2)- 80.d0*f(nx-3)+ 30.d0*f(nx-4)- 8.d0*f(nx-5) &
                      +      f(nx-6))
    end function xder
!
    function xder2(f)
!
      implicit none
!
      integer :: ix
      double precision, dimension(nx), intent(in) :: f
      double precision :: fac
      double precision, dimension(nx) :: xder2
!
!     assume uniform mesh
!
      fac=1.d0/(180.d0*(x(2)-x(1))**2)
!
      do ix=4,nx-3
        xder2(ix)=fac*(-490.d0*f(ix)             &
                       +270.d0*(f(ix+1)+f(ix-1)) &
                       - 27.d0*(f(ix+2)+f(ix-2)) &
                       +  2.d0*(f(ix+3)+f(ix-3)) &
                      )
      enddo
!
      xder2(1)=fac*(812.d0*f(1)-3132.d0*f(2)+5265.d0*f(3)-5080.d0*f(4)+2970.d0*f(5)-972.d0*f(6)+137.d0*f(7))
      xder2(2)=fac*(137.d0*f(1)- 147.d0*f(2)- 255.d0*f(3)+ 470.d0*f(4)- 285.d0*f(5)+ 93.d0*f(6)- 13.d0*f(7))
      xder2(3)=fac*(-13.d0*f(1)+ 228.d0*f(2)- 420.d0*f(3)+ 200.d0*f(4)+  15.d0*f(5)- 12.d0*f(6)+  2.d0*f(7))
!
!     outer points
!
      xder2(nx  )=fac*( 812.d0*f(nx)-3132.d0*f(nx-1)+5265.d0*f(nx-2)-5080.d0*f(nx-3)+2970.d0*f(nx-4) &
                       -972.d0*f(nx-5) +137.d0*f(nx-6))
      xder2(nx-1)=fac*( 137.d0*f(nx)- 147.d0*f(nx-1)- 255.d0*f(nx-2)+ 470.d0*f(nx-3)- 285.d0*f(nx-4) &
                       + 93.d0*f(nx-5) - 13.d0*f(nx-6))
      xder2(nx-2)=fac*( -13.d0*f(nx)+ 228.d0*f(nx-1)- 420.d0*f(nx-2)+ 200.d0*f(nx-3)+  15.d0*f(nx-4) &
                        - 12.d0*f(nx-5) + 2.d0*f(nx-6))
    end function xder2
end module deriv
!*****************************************************
module galaxy_model
  use input_params
  use calc_params
  use profiles
!
  implicit none
!
  contains
    subroutine construct_galaxy_model(gal_id_string,info)
!
      character (len=8), intent(in) :: gal_id_string
      integer, intent(in) :: info
!
!     SET UP NUMERICS AND PLOTTING
!
      if (Time_evol) then
        call set_input_params(gal_id_string,info)  !read in model parameters
        call set_calc_params  !set other parameters
      endif
!
      call construct_profiles  !construct profiles of h(r), om(r), Uz(r), Ur(r), l(r), v(r), etat(r), tau(r), n(r), Beq(r), alp_k(r)
    end subroutine construct_galaxy_model
end module galaxy_model
!*****************************************************
module bzcalc  !Calculates |Bz| using Div B=0 in the no-z approximation
  use math_constants
  use calc_params
  use var
  use grid
  use ts_params
  use deriv
  use profiles
!
  implicit none
!
  double precision, dimension(nx) :: Bzmod
!
  contains
    subroutine estimate_Bzmod(f)
!
      double precision, dimension(nx,nvar), intent(in) :: f
!
      if (r(nxghost+1)==0.d0) then
        Bzmod= lambda*h*abs(xder(f(:,1)))
      else
        Bzmod= lambda*h*abs(f(:,1)/r +xder(f(:,1)))
      endif
    end subroutine estimate_Bzmod
end module bzcalc
!*****************************************************
module equ  !Contains the partial differential equations to be solved
  use galaxy_model
  use deriv
  use boundary_conditions
  use bzcalc
!  
  implicit none
!
  double precision, dimension(nx) :: Br, Bp, Fr, Fp, Er, Ep, dBrdr, d2Brdr2, dBpdr, d2Bpdr2
  double precision, dimension(nx) :: alp_m, alp, dalp_mdr, d2alp_mdr2, dalpdr, detatdr
  double precision, dimension(nx) :: Bsqtot, DivVishniac, Dyn_gen, sqrt_Dyn_gen, d2Bsqtotdr2, &
                                     d2Dyn_gendr2, d2sqrt_Dyn_gendr2
  double precision, dimension(nx) :: brms, Bp_floor, gamma_anal, Dcrit_anal, d2gamma_analdr2
  integer :: i, imax, ifloor1, ifloor2
  double precision :: rmax, delta_r, hmax, lmax, Ncells, rmax_B, delta_r_B, rmax_D, delta_r_D
!
  contains
    subroutine pde(f,dfdt)
!
      double precision, dimension(nx,nvar) :: f, dfdt
      intent(inout) :: f
      intent(out) :: dfdt
!
      call impose_bc(f)
!
!     USE EXPLICIT NAMES
      Br=f(:,1)
      Bp=f(:,2)
!
      if (Damp) then
        Fr=f(:,3)
        Fp=f(:,4)
        if (Dyn_quench) then
          Er=f(:,5)
          Ep=f(:,6)
          alp_m=f(:,7)
        endif
      else
        Fr=0*Br
        Fp=0*Bp
        if (Dyn_quench) then
          alp_m=f(:,3)
        endif
      endif
!
      dBrdr=xder(Br)
      d2Brdr2=xder2(Br)
      dBpdr=xder(Bp)
      d2Bpdr2=xder2(Bp)
!
      if (Dyn_quench) then
        dalp_mdr=xder(alp_m)
        d2alp_mdr2=xder2(alp_m)
      endif
!
!     CALCULATE MAGNETIC ENERGY (WITHOUT THE FACTOR 1/(8PI))
      Bsqtot=Br**2 +Bp**2 +Bzmod**2
      d2Bsqtotdr2= xder2(Bsqtot)
!
      if (Alg_quench) then
        alp= alp_k/(1.d0 +Bsqtot/Beq**2)  !Formula for simple alpha quenching
      elseif (Dyn_quench) then
        alp= alp_k +alp_m  !Total alpha is equal to the sum of kinetic and magnetic parts
      else
        alp= alp_k
      endif
!
      detatdr= 1.d0/3*l*dvdr+1.d0/3*v*dldr
!
!     IMPOSE MINIMUM (FLOOR) ON B_PHI DUE TO SMALL-SCALE TURBULENT FLUCTUATING MAGNETIC FIELD
!
!     CALCULATE DYNAMO NUMBER
      Dyn_gen=G*alp*h**3/etat**2
      sqrt_Dyn_gen=dsqrt(abs(Dyn_gen))
      d2Dyn_gendr2=xder2(Dyn_gen)
      d2sqrt_Dyn_gendr2=xder2(sqrt_Dyn_gen)
!
      if (Floor) then
        Dcrit_anal= -pi/32*(pi**2 +12*C_U*h*Uz/l/v)**2
        gamma_anal= dsqrt(2.d0/pi)*(dsqrt(abs(Dyn_gen))-dsqrt(abs(Dcrit_anal)))*etat/h**2
        d2gamma_analdr2= xder2(gamma_anal)
        do i=nxghost+1,nx-nxghost-1
          if (abs(Bp(i))==maxval(abs(Bp))) then
!          if (abs(Bsqtot(i))==maxval(abs(Bsqtot))) then
            rmax=r(i)  !radius at max of Bsqtot(r)
!            rmax_B=r(i)  !radius at max of Bsqtot(r)
            delta_r=2*dsqrt(abs(Bp(i)/d2Bpdr2(i)))  !width of Gaussian approx to Bp(r)
!            delta_r_B=2*dsqrt(abs(Bp(i)/d2Bpdr2(i)))  !width of Gaussian approx to Bp(r)
!            delta_r_B=2*dsqrt(abs(Bsqtot(i)/d2Bsqtotdr2(i)))  !width of Gaussian approx to Bsqtot(r)
            hmax=h(i)  !scale height at max of Bp(r)
            lmax=l(i)  !turbulent scale at max of Bp(r)
            imax=i
          endif
!          if (abs(Dyn_gen(i))==maxval(abs(Dyn_gen))) then
!            rmax_D=r(i)  !radius at max of Dyn_gen(r)
!            delta_r_D=2*dsqrt(abs(sqrt_Dyn_gen(i)/d2sqrt_Dyn_gendr2(i)))  !width of Gaussian approx to Dyn_gen(r)
!          endif
!          if (gamma_anal(i)==maxval(gamma_anal(nxghost+1:nx-nxghost-1))) then
!            rmax=r(i)  !radius at max of gamma_anal
!            delta_r=2*dsqrt(abs(gamma_anal(i)/d2gamma_analdr2(i)))  !width of Gaussian approx to Bp(r)
!            hmax=h(i)  !scale height at max of Bp(r)
!            lmax=l(i)  !turbulent scale at max of Bp(r)
!            imax=i
!          endif
        enddo
        print*, 'gamma_anal=', gamma_anal
        print*, '-Dcrit_anal=', -Dcrit_anal
        print*, '-Dyn_gen=', -Dyn_gen
        print*, 'sqrt(-Dyn_gen)-sqrt(-Dcrit_anal)=', sqrt(-Dyn_gen)-sqrt(-Dcrit_anal)
        print*, 'td=', h**2/etat
        print*, 'rmax=',rmax*r_disk_kpc
        print*, 'rmax_D=',rmax_D*r_disk_kpc
        print*, 'rmax_B=',rmax_B*r_disk_kpc
        print*, 'delta_r=',delta_r*r_disk_kpc
        print*, 'delta_r_D=',delta_r_D*r_disk_kpc
        print*, 'delta_r_B=',delta_r_B*r_disk_kpc
        Ncells= 3.d0*rmax*delta_r*hmax/lmax**3/lambda**2
        brms= fmag*Beq
        Bp_floor= brms(imax)/dsqrt(Ncells)*lmax/delta_r*lambda !brms/dsqrt(Ncells)*lmax/delta_r*lambda
print*, Bp_floor
        do i=nxghost+1,nx-nxghost
          if (abs(r(i)-(rmax-delta_r/2)) == minval(abs(r-(rmax-delta_r/2)))) then
            ifloor1=i
          endif
          if (abs(r(i)-(rmax+delta_r/2)) == minval(abs(r-(rmax+delta_r/2)))) then
            ifloor2=i
          endif
        enddo
        Bp(ifloor1:ifloor2)=Bp(ifloor1:ifloor2)*(1.d0 +abs(Bp_floor(ifloor1:ifloor2)/Bp(ifloor1:ifloor2)))  !to leave the sign of Bp unchanged
        do i=ifloor1,ifloor2
          if (abs(Bp(i))<Bp_floor(i)) then
            Bp(i)=sign(Bp_floor(i),Bp(i))
          endif
        enddo
      endif
!      print*, 'r(ifloor1:ifloor2)='
!      print*, r(ifloor1:ifloor2)
!      print*, 'Bp_floor(ifloor1:ifloor2)='
!      print*, Bp_floor(ifloor1:ifloor2)
!      print*, 'Bp(ifloor1:ifloor2)='
!      print*, Bp(ifloor1:ifloor2)
!      print*, 'Bp(ifloor1:ifloor2)*(1.d0 +abs(Bp_floor(ifloor1:ifloor2)/Bp(ifloor1:ifloor2)))='
!      print*, Bp(ifloor1:ifloor2)*(1.d0 +abs(Bp_floor(ifloor1:ifloor2)/Bp(ifloor1:ifloor2)))
!      print*,'rmax (kpc)=',rmax*r_disk_kpc,'delta_r (kpc)=',delta_r*r_disk_kpc,'hmax (kpc)=',hmax*h0_kpc,'lmax (kpc)=',lmax*h0_kpc
!      print*,'Ncells=',Ncells
!      print*,'Bp_floor(imax)=',Bp_floor(imax)
!      print*, 'Bp_floor_with_sign_of_Bp=',sign(Bp_floor,Bp)
!      print*, 'ifloor1, ifloor2,delta_r,rmax,r(ifloor1),r(ifloor2)=', &
!              ifloor1,ifloor2,delta_r*r_disk_kpc,rmax*r_disk_kpc, &
!              r(ifloor1)*r_disk_kpc,r(ifloor2)*r_disk_kpc
!      print*, 'r(ifloor1)+delta_r/2,r(ifloor2)-delta_r/2,dx', &
!              (r(ifloor1)+delta_r/2)*r_disk_kpc,(r(ifloor2)-delta_r/2)*r_disk_kpc,dx
!      print*, 'r(ifloor1-1)+delta_r/2,r(ifloor2+1)-delta_r/2,dx', &
!              (r(ifloor1+1)+delta_r/2)*r_disk_kpc,(r(ifloor2-1)-delta_r/2)*r_disk_kpc
!
!     LIST OF VARIABLE NAMES FOR f ARRAY
!   
!     UNDER FOSA          UNDER TAU APPROXIMATION
!     f(:,1)=Br           f(:,1)=Br
!     f(:,2)=Bp           f(:,2)=Bp
!     f(:,3)=alp_m        f(:,3)=Fr
!                         f(:,4)=Fp
!                         f(:,5)=Er
!                         f(:,6)=Ep
!                         f(:,7)=alp_m
!   
!     METHOD: ALL ARRAYS ARE 2-DIMENSIONAL
!
      r(nxghost+1)=0.000001d0
      dalp_mdr(nxghost+1)=0.d0
!
      if (.not.Damp) then
!       CASE 1: FOSA (tau-->0 LIMIT)--NOTE: dfdt BLOWS UP AT ORIGIN BUT SET IT TO 0 ANYWAY
        dfdt(:,1)=      -C_U*Uz*Br/h -lambda*Ur*Br/r -lambda*Ur*dBrdr         &
                   -2.d0/pi/h*ctau*alp*Bp -pi**2/4/h**2*(ctau+Rm_inv)*etat*Br &
                   +(ctau+Rm_inv)*etat*lambda**2*(-Br/r**2 +dBrdr/r +d2Brdr2)                     
!                   +ctau*( detatdz*dBrdz -detatdz*lambda*dBzdr)  !Contains detatdz terms
        dfdt(:,2)= G*Br -C_U*Uz*Bp/h -lambda*dUrdr*Bp -lambda*Ur*dBpdr        &
                   -2.d0/pi/h*ctau*alp*Br -pi**2/4/h**2*(ctau+Rm_inv)*etat*Bp &
                   +(ctau+Rm_inv)*etat*lambda**2*(-Bp/r**2 +dBpdr/r +d2Bpdr2) &
                   +ctau*lambda**2*( detatdr*Bp/r +detatdr*dBpdr)  !Contains detatdr terms 
!                   +ctau*(detatdz*dBpdz)  !Contains detatdz terms
        if (.not.Alp_squared) then
          dfdt(:,2)=dfdt(:,2) +2.d0/pi/h*ctau*alp*Br
        endif
        if (Dyn_quench) then
          dfdt(:,3)= -2*(h0_kpc/l_kpc)**2*etat*(ctau*alp*(Br**2+Bp**2+Bzmod**2)/Beq**2                      &
                     +ctau*3*etat/pi**(3.d0/2)/h*abs(Dyn_gen)**(1.d0/2)*Br*Bp/Beq**2+Rm_inv*alp_m) &
                     -C_a*alp_m*Uz/h -lambda*alp_m*Ur/r -lambda*alp_m*dUrdr -lambda*Ur*dalp_mdr    &
                     +R_kappa*etat*(lambda**2*d2alp_mdr2 +lambda**2/r*dalp_mdr +C_d/h**2*alp_m)
        endif
      else
!       CASE 2: MTA (FINITE tau)--NOTE: dfdt BLOWS UP AT ORIGIN BUT SET IT TO 0 ANYWAY
        dfdt(:,1)=       -C_U*Uz*Br/h -lambda*Ur*Br/r -lambda*Ur*dBrdr +Fr                       &
                   +Rm_inv*(-pi**2/4/h**2*etat*Br +etat*lambda**2*(-Br/r**2 +dBrdr/r +d2Brdr2))
        dfdt(:,2)=  G*Br -C_U*Uz*Bp/h                 -lambda*dUrdr*Bp   -lambda*Ur*dBpdr +Fp    &
                   +Rm_inv*(-pi**2/4/h**2*etat*Bp +etat*lambda**2*(-Bp/r**2 +dBpdr/r +d2Bpdr2))
        dfdt(:,3)=  tau**(-1)*(-2.d0/pi/h*ctau*alp*Bp -pi**2/4/h**2*ctau*etat*Br                 &
                   +ctau*etat*lambda**2*(-Br/r**2 +dBrdr/r +d2Brdr2) -Fr)                          
!                   +ctau*( detatdz*dBrdz -detatdz*lambda*dBzdr)  !Contains detatdz terms
        dfdt(:,4)=  tau**(-1)*(-2.d0/pi/h*ctau*alp*Br -pi**2/4/h**2*ctau*etat*Bp                 &
                   +ctau*etat*lambda**2*(-Bp/r**2 +dBpdr/r +d2Bpdr2) -Fp                         &
                   +ctau*lambda**2*( detatdr*Bp/r +detatdr*dBpdr))  !Contains detatdr terms
!                   +ctau*detatdz*dBpdz  !Contains detatdz terms
        if (.not.Alp_squared) then
          dfdt(:,4)=dfdt(:,4) +tau**(-1)*2./pi/h*ctau*alp*Br
        endif
        if (Dyn_quench) then
          dfdt(:,5)=  tau**(-1)*(ctau*alp*Br -ctau*etat*pi/2/h                                                 *Bp -Er)
          dfdt(:,6)=  tau**(-1)*(ctau*alp*Bp +ctau*etat*pi/2/h*(1. +1.d0/2/pi**(3.d0/2)*abs(Dyn_gen)**(1.d0/2))*Br -Ep)
          dfdt(:,7)= -2*(h0_kpc/l_kpc)**2*etat*((Er*Br +Ep*Bp)/Beq**2 +Rm_inv*alp_m)             &
                     -C_a*alp_m*Uz/h -lambda*alp_m*Ur/r -lambda*alp_m*dUrdr -lambda*Ur*dalp_mdr  &
                     +R_kappa*etat*(lambda**2*d2alp_mdr2 +lambda**2/r*dalp_mdr +C_d/h**2*alp_m)
        endif
      endif
    end subroutine pde
end module equ
!*****************************************************
module timestep  !Contains time-stepping routine
  use var
  use grid
  use ts_params
  use equ
!
contains
  subroutine rk(f)
!
  implicit none
!
  double precision :: gam1, gam2, gam3, zet1, zet2
  double precision, dimension(nx,nvar) :: f, dfdt, pdef, ftmp
  double precision :: ttmp
!
  intent(inout) :: f
!
!  Runge Kutta 3rd order time advance
!  f = f(exact) - 0.0046*dt**3*d3f/dt3
!
    gam1=8.d0/15 !gam1, gam2 AND gam3 ARE THE COEFFICIENTS OF THE TIMESTEPS AT WHICH dfdt IS CALCULATED
    gam2=5.d0/12
    gam3=3.d0/4
    zet1=-17.d0/60
    zet2=-5.d0/12
!
    call pde(f,dfdt)
    pdef=dfdt !pde FUNCTION CALCULATES VALUE OF TIME DERIVATIVE ACCORDING TO P.D.E.s
    f=f+dt*gam1*pdef !FIRST GET INTERMEDIATE VALUE OF f (AT t=t+dt*gam1) USING dfdt AT t_0
    t=t+dt*gam1 !THEN GO TO THAT TIMESTEP
    ftmp=f+dt*zet1*pdef !NOW CALCULATE A TEMPORARY f (AT t=t_0+dt*gam1+dt*zet1) USING dfdt AT t_0
    ttmp=t+dt*zet1 !DEFINE TEMPORARY t AS THAT TIMESTEP (t=t_0+dt*gam1+dt*zet1)
!
    call pde(f,dfdt)
    pdef=dfdt !NOW CALCULATE dfdt AT THE NEW TIMESTEP t_0+dt*gam1
    f=ftmp+dt*gam2*pdef !USE THIS TO GET ANOTHER INTERMEDIATE VALUE OF f (AT t=t_0+dt*gam1+dt*zet1+dt*gam2)
    t=ttmp+dt*gam2 !THEN GO TO THAT TIMESTEP
    ftmp=f+dt*zet2*pdef !NOW CALCULATE A TEMPORARY f (AT t=t_0+dt*gam1+dt*zet1+dt*gam2+dt*zet2) USING dfdt AT t=t_0+dt*gam1
    ttmp=t+dt*zet2 !DEFINE TEMPORARY t AS THAT TIMESTEP (t=t_0+dt*gam1+dt*zet1+dt*gam2+dt*zet2)
!
    call pde(f,dfdt)
    pdef=dfdt !CALCULATE dfdt AT THE NEW TIMESTEP t_0+dt*gam1+dt*zet1+dt*gam2
    f=ftmp+dt*gam3*pdef !USE THIS TO GET THE FINAL VALUE OF f (AT t=t_0+dt*gam1+dt*zet1+dt*gam2+dt*zet2+dt*gam3)
    t=ttmp+dt*gam3 !THEN GO TO THAT TIMESTEP
!
    first=first+1. !COUNTS THE NUMBER OF TIMES RUNGA-KUTTA ROUTINE IS EXCUTED
  end subroutine rk
end module timestep
!*****************************************************
module diagnostic  !Writes diagnostic information to screen, file
  use direc_names
  use math_constants
  use units
  use input_params
  use calc_params
  use ts_params
  use grid
  use modules
  use var
  use galaxy_model
  use initial_conditions
  use equ
  use bzcalc
!
  implicit none
!
  contains
    subroutine print_info(info)
      
      integer, intent(in) :: info
       
      if (info>1) then
        print*,''
        print*,'UNITS'
        print*,'t0 (Gyr)=',t0_Gyr,'t0 (kpc s/km)=',t0_kpcskm,'h0 (kpc)=',h0_kpc,'etat0 (cm2s)= ',etat0_cm2s
        print*,''
        print*,'NUMERICS:'
        print*,'nvar=     ',nvar     
        print*,'dt=       ',dt       ,'   n1=       ',n1       ,'   n2=     ',n2
        print*,'dx=       ',dx       ,'   nxphys=   ',nxphys   ,'   nxghost=',nxghost,'   nx=',nx
        print*,'minval(x)=',minval(x),'   maxval(x)=',maxval(x)
        print*,''
        print*,'PHYSICAL GRID:'
        print*,'r_kpc=',r_kpc
        print*,''
        print*,'MODULES:' 
        print*,'Read_param= ',Read_param ,'Dyn_quench=',Dyn_quench,'Alg_quench=   ',Alg_quench   ,'Damp=       ' ,Damp
        print*,'Var_Uz=     ',Var_Uz     ,'Var_l=     ',Var_l     ,'Var_v=        ',Var_v        ,'Var_n=      ',Var_n
        print*,'Flaring=    ',Flaring    ,'Rand_seed= ',Rand_seed ,'Alp_ceiling=  ',Alp_ceiling  ,'Om_Brandt=  ',Om_Brandt
        print*,'Alp_squared=',Alp_squared,'Krause=    ',Krause    ,'Shear=        ',Shear        ,'Advect=     ',Advect
        print*,'Turb_dif=   ',Turb_dif
        print*,''
        print*,'INPUT PARAMETERS:'
        print*,'r_disk_kpc= ',r_disk_kpc ,'R_kappa=   ',R_kappa
        print*,'l_sol_kpc=  ',l_sol_kpc  ,'r_l_kpc=   ',r_l_kpc   ,'v_sol_kms=    ',v_sol_kms    ,'r_v_kpc=    ',r_v_kpc
        print*,'n_sol_cm3=  ',n_sol_cm3  ,'r_n_kpc=   ',r_n_kpc   ,'Uz_sol_kms=   ',Uz_sol_kms   ,'r_Uz_kpc=   ',r_Uz_kpc
        print*,'h_sol_kpc=  ',h_sol_kpc  ,'r_h_kpc=   ',r_h_kpc   ,'Uphi_sol_kms= ',Uphi_sol_kms ,'r_om_kpc=   ',r_om_kpc
        print*,''
        print*,'OTHER FREE PARAMETERS:'
        print*,'r_in= ',r_in ,'r_disk_kpc=',r_disk_kpc,'r_sol_kpc=',r_sol_kpc,'r1_kpc=',r1_kpc
        print*,'ctau= ',ctau ,'nn=        ',nn        ,'lambda=   ',lambda
        print*,'C_alp=',C_alp,'alpceil=   ',alpceil   ,'Rm_inv=   ',Rm_inv
        print*,''
        print*,'IMPORTANT CALCULATED PARAMETERS'
        print*,'etat_sol_cm2s=',etat_sol_cm2s,'td_sol_Gyr=',td_sol_Gyr,'tau_sol_Gyr=',tau_sol_Gyr
        print*,'etat_sol=     ',etat_sol     ,'td_sol=    ',td_sol    ,'tau_sol     ',tau_sol
        print*,''
        print*,'INPUT/OUTPUT INFORMATION'
        print*,'Directory ending   =',s0
        print*,'home directory     =',s1
      endif
!
!     WRITE INFO TO FILE "diagnostic.out"
      open(20,file= 'diagnostic.out',status="replace")
!
      write(20,*),'UNITS'
      write(20,*),'t0 (Gyr)=',t0_Gyr,'t0 (kpc s/km)=',t0_kpcskm,'h0 (kpc)=',h0_kpc,'etat0 (cm2s)= ',etat0_cm2s
      write(20,*),''
      write(20,*),'NUMERICS:'
      write(20,*),'nvar=     ',nvar     
      write(20,*),'dt=       ',dt       ,'   n1=       ',n1       ,'   n2=     ',n2
      write(20,*),'dx=       ',dx       ,'   nxphys=   ',nxphys   ,'   nxghost=',nxghost,'   nx=',nx
      write(20,*),'minval(x)=',minval(x),'   maxval(x)=',maxval(x)
      write(20,*),''
      write(20,*),'PHYSICAL GRID:'
      write(20,*),'r_kpc=',r_kpc
      write(20,*),''
      write(20,*),'MODULES:' 
      write(20,*),'Read_param= ',Read_param ,'Dyn_quench=',Dyn_quench,'Alg_quench=   ',Alg_quench   ,'Damp=       ' ,Damp
      write(20,*),'Var_Uz=     ',Var_Uz     ,'Var_l=     ',Var_l     ,'Var_v=        ',Var_v        ,'Var_n=      ',Var_n
      write(20,*),'Flaring=    ',Flaring    ,'Rand_seed= ',Rand_seed ,'Alp_ceiling=  ',Alp_ceiling  ,'Om_Brandt=  ',Om_Brandt
      write(20,*),'Alp_squared=',Alp_squared,'Krause=    ',Krause    ,'Shear=        ',Shear        ,'Advect=     ',Advect
      write(20,*),'Turb_dif=   ',Turb_dif
      write(20,*),''
      write(20,*),'INPUT PARAMETERS:'
      write(20,*),'r_disk_kpc= ',r_disk_kpc ,'R_kappa=   ',R_kappa
      write(20,*),'l_sol_kpc=  ',l_sol_kpc  ,'r_l_kpc=   ',r_l_kpc   ,'v_sol_kms=    ',v_sol_kms    ,'r_v_kpc=    ',r_v_kpc
      write(20,*),'n_sol_cm3=  ',n_sol_cm3  ,'r_n_kpc=   ',r_n_kpc   ,'Uz_sol_kms=   ',Uz_sol_kms   ,'r_Uz_kpc=   ',r_Uz_kpc
      write(20,*),'h_sol_kpc=  ',h_sol_kpc  ,'r_h_kpc=   ',r_h_kpc   ,'Uphi_sol_kms= ',Uphi_sol_kms ,'r_om_kpc=   ',r_om_kpc
      write(20,*),''
      write(20,*),'OTHER FREE PARAMETERS:'
      write(20,*),'r_in= ',r_in ,'r_disk_kpc=',r_disk_kpc,'r_sol_kpc=',r_sol_kpc,'r1_kpc=',r1_kpc
      write(20,*),'ctau= ',ctau ,'nn=        ',nn        ,'lambda=   ',lambda
      write(20,*),'C_alp=',C_alp,'alpceil=   ',alpceil   ,'Rm_inv=   ',Rm_inv
      write(20,*),''
      write(20,*),'IMPORTANT CALCULATED PARAMETERS'
      write(20,*),'etat_sol_cm2s=',etat_sol_cm2s,'td_sol_Gyr=',td_sol_Gyr,'tau_sol_Gyr=',tau_sol_Gyr
      write(20,*),'etat_sol=     ',etat_sol     ,'td_sol=    ',td_sol    ,'tau_sol     ',tau_sol
      write(20,*),''
      write(20,*),'INPUT/OUTPUT INFORMATION'
      write(20,*),'Directory ending   =',s0
      write(20,*),'home directory     =',s1
      write(20,*),''
      close(20)
    end subroutine print_info
end module diagnostic
!*****************************************************
module initial_data_dump
  use direc_names
  use math_constants
  use units
  use input_params
  use calc_params
  use ts_params
  use grid
  use modules
  use var
  use galaxy_model
  use initial_conditions
  use equ
  use bzcalc
!
  contains
    subroutine write_initial_data(f,gal_id_string,info)
!
      double precision, dimension(nx,nvar), intent(in) :: f
      character (len=8), intent(in) :: gal_id_string
      integer, intent(in) :: info
!
!     WRITE DATA TO FILE "param.out"
      if (info> 0) then
        print*,''
        print*,'Writing parameters to file param',gal_id_string,'.out'
      endif
      open(20,file= 'diagnostic.out',status="old",position="append")
      write(20,*)''
      write(20,*)'Writing parameters to file param',gal_id_string,'.out'
      close(20)
      open(10,file= s1 // '/output/' // s0 &
                    // '/param_' // gal_id_string // '.out',status="replace")
      write(10,*)t0_Gyr,t0_kpcskm,h0_kpc,etat0_cm2s,n0_cm3,B0_mkG
      write(10,*)nvar,dt,n1,n2,dx,nxphys,nxghost,nx
      write(10,*)l_sol_kpc,r_l_kpc,v_sol_kms,r_v_kpc,n_sol_cm3,r_n_kpc,Uz_sol_kms,r_Uz_kpc, &
                 h_sol_kpc,r_h_kpc,Uphi_sol_kms,r_om_kpc,R_kappa
      write(10,*)r_in,r_disk_kpc,r_sol_kpc,r1_kpc,ctau,nn,lambda,C_alp,alpceil,Rm_inv
      write(10,*)etat_sol_cm2s,td_sol_Gyr,tau_sol_Gyr,etat_sol,td_sol,tau_sol
      close(10)
      if (info> 0) then
        print*,'Finished writing parameters to file param',gal_id_string,'.out'
        print*,'Writing initial output to file init',gal_id_string,'.out'
      endif
      open(20,file= 'diagnostic.out',status="old",position="append")
      write(20,*)'Finished writing parameters to file param',gal_id_string,'.out'
      write(20,*)'Writing initial output to file init',gal_id_string,'.out'
      write(20,*)''
      close(20)
!
!     WRITE DATA TO FILE "init.out"
      open(11,file= s1 // '/output/' // s0 &
                    // '/init_' // gal_id_string // '.out' ,status="replace")
      write(11,*)r
      write(11,*)h
      write(11,*)om
      write(11,*)G
      write(11,*)Uz
      write(11,*)Ur
      write(11,*)l
      write(11,*)v
      write(11,*)etat
      write(11,*)tau
      write(11,*)alp_k
      write(11,*)n
      write(11,*)Beq
      write(11,*)f(:,1)
      write(11,*)f(:,2)
      write(11,*)Bzmod
      close(11)
      if (info> 0) then
        print*,'Finished writing initial output to file init',gal_id_string,'.out'
      endif
      open(20,file= 'diagnostic.out',status="old",position="append")
      write(20,*)'Finished writing initial output to file init',gal_id_string,'.out'
      close(20)
    end subroutine write_initial_data
end module initial_data_dump
!*****************************************************
module start  !Contains initialization routine for simulation
  use direc_names
  use math_constants
  use units
  use input_params
  use calc_params
  use ts_params
  use grid
  use modules
  use galaxy_model
  use initial_conditions
  use var
  use equ
  use bzcalc
  use diagnostic
  use initial_data_dump
!
  implicit none
  double precision, allocatable, dimension(:,:) :: f, dfdt
!
  contains
    subroutine init_start(gal_id_string,info)
!
      character (len=8), intent(in) :: gal_id_string
      integer, intent(in) :: info
!
!     INITIALIZE VARIABLE ARRAY
      call init_var
!
      if (.not. allocated(f)) allocate(f(nx,nvar)) !Copied from Pitch/run.f90
      if (.not. allocated(dfdt)) allocate(dfdt(nx,nvar)) !Copied from Pitch/run.f90
!
!     SET INPUT PARAMETERS
!      call set_input_params(gal_id_string,info)
!
!     SET UP NUMERICS AND PLOTTING
      call construct_grid
!
!     SET UP GALAXY MODEL
      call construct_galaxy_model(gal_id_string,info)
!
!     SEED FIELD
      call init_seed(f)
!  
!     CALCULATE |Bz|
      call estimate_Bzmod(f)
!
      if (.not.Turb_dif) then
        etat=0.d0
      endif
!
      if (.not.Advect) then
        om=0.d0
      endif
!
!     PRINT INFO TO SCREEN AND DIAGNOSTIC FILE "diagnostic.out"
      call print_info(info)
!
!     WRITE INITIAL DATA TO FILES "param.out" and "init.out" and continue appending "diagnostic.out"
      call write_initial_data(f,gal_id_string,info)
    end subroutine init_start
end module start
!*****************************************************
module output_dump
!
  contains
    subroutine write_final_output(f,gal_id_string,info)
!
      use direc_names
      use galaxy_model
      use bzcalc
      use ts_arrays
!  
      double precision, dimension(nx,nvar), intent(in) :: f
      character (len=8), intent(in) :: gal_id_string
      integer, intent(in) :: info
  
!     WRITE DATA FOR FINAL TIMESTEP TO FILE "run.out"
      if (info> 0) then
        print*,'Writing output for final timestep to file run',gal_id_string,'.out'
      endif
      open(20,file= 'diagnostic.out',status="old",position="append")
      write(20,*)'Writing output for final timestep to file run',gal_id_string,'.out'
      close(20)
      open(12,file= s1 // '/output/' // s0 &
                    // '/run_' // gal_id_string // '.out',status="replace")
      write(12,*)t
      write(12,*)f(:,1)
      write(12,*)f(:,2)
      if (Dyn_quench) then
        if (.not.Damp) then
          write(12,*)f(:,3)
        else
          write(12,*)f(:,7)
        endif
      endif
      write(12,*)Bzmod
      write(12,*)h
      write(12,*)om
      write(12,*)G
      write(12,*)l
      write(12,*)v
      write(12,*)etat
      write(12,*)tau
      write(12,*)alp_k
      write(12,*)Uz
      write(12,*)Ur
      write(12,*)n
      write(12,*)Beq
      close(12)
!  
!     WRITE DATA FOR ALL TIMESTEPS TO FILE "ts.out"
      if (info> 0) then
        print*,'Writing time series output to file ts',gal_id_string,'.out'
      endif
      open(13,file= s1 // '/output/' // &
                    s0 // '/ts_' // gal_id_string // '.out',status="replace")
      write(13,*)ts_t
      write(13,*)ts_Br
      write(13,*)ts_Bp
      if (Dyn_quench) then
        if (.not.Damp) then
          write(13,*)ts_alp_m
        else
          write(13,*)ts_alp_m
        endif
      endif
      write(13,*)ts_Bzmod
      write(13,*)ts_h
      write(13,*)ts_om
      write(13,*)ts_G
      write(13,*)ts_l
      write(13,*)ts_v
      write(13,*)ts_etat
      write(13,*)ts_tau
      write(13,*)ts_alp_k
      write(13,*)ts_Uz
      write(13,*)ts_Ur
      write(13,*)ts_n
      write(13,*)ts_Beq
      write(13,*)ts_rmax
      write(13,*)ts_delta_r
      close(13)
print*,'ts_Br(1,1)'   ,ts_Br(1,1)
print*,'ts_Br(1,37)'  ,ts_Br(1,37)
print*,'ts_Br(n1+1,1)' ,ts_Br(n1+1,1)
print*,'ts_Br(n1+1,37)',ts_Br(n1+1,37)
print*,'ts_Br(10,10)',ts_Br(10,10)
    end subroutine write_final_output
end module output_dump
!*****************************************************
