subroutine XYZJ(J2,r,a,n,e,Oe,Oq,Oh,OO,X,Y,Z)
!
! This subroutine calculates the tidal factors associated with 
! the permanent J2 quadrupole moment, as derived by Ragozzine &
! Brown (2009). 
!
! Sep-28-2019/wlyra
!
  real :: fac
  real, intent(in):: J2,r,a,n,e,Oe,Oq,Oh,OO
  real, intent(out):: X,Y,Z
  
  fac = 1.5 * J2 * (r/a)**2 * n/(1-e**2)**2/OO**2
  X = fac * Oe*Oh
  Y = fac * Oq*Oh
  Z = fac * (2*Oh**2 - Oe**2 - Oq**2)*.5
!
endsubroutine XYZJ
!***********************************************************************  
subroutine XYZVW(Qtidal,klove,m1,m2,GNewton,sqrtGM1,r,a,mu,n,e,Oe,Oq,Oh,X,Y,Z,V,W)
!
! This subroutine calculates the tidal factors associated with 
! the induced quadrupole moment, as in the model of Eggleton et al. (2001),
! adapted by Fabrycky & Tremaine (2007) for the planetary case of tidal dissipation. 
!
! Sep-28-2019/wlyra
!  
  real :: tf
  real, intent(in):: Qtidal,klove,m1,m2,GNewton,sqrtGM1,r,a,mu,n,Oe,Oq,Oh,e
  real, intent(out):: X,Y,Z,V,W
!
  tf = 1./6 * Qtidal/klove * m1/m2 * sqrtGM1 * 1./r**5 * a**(13./2)
  V = 9./tf * ((1+15./4*e**2 + 15./8*e**4 + 5./64*e**6)/(1-e**2)**(13./2) - 11./18*Oh/n*(1+3./2*e**2+1./8*e**4)/(1-e**2)**5)
  W = 1./tf * ((1+15./2*e**2 + 45./8*e**4 + 5./16*e**6)/(1-e**2)**(13./2) - Oh/n*(1+3.*e**2+3./8*e**4)/(1-e**2)**5)
!       
  X = -m2*klove*(r/a)**5/(2*mu*n) * Oh*Oe/(1-e**2)**2 - Oq*(1+9./2*e**2+5./8*e**4)/(2*n*tf*(1-e**2)**5)
  Y = -m2*klove*(r/a)**5/(2*mu*n) * Oh*Oq/(1-e**2)**2 + Oe*(1+3./2*e**2+1./8*e**4)/(2*n*tf*(1-e**2)**5)
  Z =  m2*klove*(r/a)**5/(2*mu*n) * ((2*Oh**2-Oe**2-Oq**2)/(2*(1-e**2)**2) &
       + 15*GNewton*m2/a**3 * (1+3./2*e**2+1./8*e**4)/(1-e**2)**5)
!
endsubroutine XYZVW
!***********************************************************************  
subroutine get_orbit(e,h,sqrtGM,GM,a,n)
!
! This subroutine takes eccentricity and angular momentum, outputting 
! semimajor axis, period, and pericenter distance. 
!
! Sep-28-2019/wlyra
!
  real, intent(in):: e,h,sqrtGM,GM
  real, intent(out):: a,n
!  
  a = h**2/(GM*(1-e**2))
  n = sqrtGM*a**(-1.5) 
!
endsubroutine get_orbit
!***********************************************************************    
subroutine cross(a,b,c)
!  
! Cross product.
!
! Sep-28-2019/wlyra
!
  real, dimension(3), intent(in):: a,b
  real, dimension(3), intent(out):: c
!  
  c(1) = a(2)*b(3) - a(3)*b(2)
  c(2) = a(3)*b(1) - a(1)*b(3)
  c(3) = a(1)*b(2) - a(2)*b(1)
!  
endsubroutine cross
!***********************************************************************   
subroutine dot(a,b,c)
!  
! Dot product.
!
! Sep-28-2019/wlyra
!
  real, dimension(3), intent(in):: a,b
  real, intent(out) :: c
!  
  c = a(1)*b(1) + a(2)*b(2) + a(3)*b(3)
!  
endsubroutine dot
!***********************************************************************      
subroutine get_eqh_components(A,x,y,z,Ax,Ay,Az)
!
! This subroutine calculates the components of a given vector A in the 
! coordinate system defined by orthogonal vectors x, y, and z. 
!
! Sep-28-2019/wlyra
!  
  real :: Asize, Asize1, A_dot_x, A_dot_y, A_dot_z
  real :: alpha, sina, cosa
  real :: beta,  sinb, cosb
!  
  real, dimension(3), intent(in) :: A
  real, dimension(3), intent(in) :: x,y,z
  real, intent(out):: Ax,Ay,Az
! 
  Asize  = sqrt(sum(A**2))
  Asize1 = 1./Asize
!        
  call dot(A*Asize1,x,A_dot_x)
  call dot(A*Asize1,y,A_dot_y)
  call dot(A*Asize1,z,A_dot_z)
!
  if (A_dot_z > 1.) A_dot_z=1.0  !ugly patch
  alpha = acos(A_dot_z)
  cosa  = cos(alpha)
  sina  = sin(alpha)
!
  beta = atan2(A_dot_y,A_dot_x)
  cosb = cos(beta)
  sinb = sin(beta)
!        
  Ax = Asize*sina*cosb
  Ay = Asize*sina*sinb
  Az = Asize*cosa
!
endsubroutine get_eqh_components
!***********************************************************************
subroutine tidalJ(a,n,e,&
           r1,O1e,O1q,O1h,Omega1,J21,&
           r2,O2e,O2q,O2h,Omega2,J22,&
           XX1,XX2,YY1,YY2,ZZ1,ZZ2)
!
! This subroutine adds the contribution of the permanent quadrupole to 
! the X, Y, and Z tidal parameters. 
!
! Sep-28-2019/wlyra
!
  real, intent(in) :: a,n,e
  real, intent(in) :: r1,O1e,O1q,O1h,Omega1,J21
  real, intent(in) :: r2,O2e,O2q,O2h,Omega2,J22
  real, intent(inout) :: XX1,XX2,YY1,YY2,ZZ1,ZZ2
  real :: XJ1,YJ1,ZJ1,XJ2,YJ2,ZJ2
!
  call XYZJ(J21,r1,a,n,e,O1e,O1q,O1h,Omega1,XJ1,YJ1,ZJ1)
  call XYZJ(J22,r2,a,n,e,O2e,O2q,O2h,Omega2,XJ2,YJ2,ZJ2)
!  
  XX1 = XX1 + XJ1; XX2 = XX2 + XJ2
  YY1 = YY1 + YJ1; YY2 = YY2 + YJ2
  ZZ1 = ZZ1 + ZJ1; ZZ2 = ZZ2 + ZJ2
!
endsubroutine tidalJ
!***********************************************************************
subroutine dissipation(a,mu,n,e,GNewton,sqrtGM1,&
           Qtidal1,klove1,m1,r1,O1e,O1q,O1h,&
           Qtidal2,klove2,m2,r2,O2e,O2q,O2h,&
           XX1,XX2,YY1,YY2,ZZ1,ZZ2,VV1,VV2,WW1,WW2)
!
! This subroutine adds the contribution of the induced quadrupole to 
! the V and W dissipation functions, as well as to the 
! X, Y, and Z tidal parameters. 
!
! Sep-28-2019/wlyra
!
  real, intent(in) :: GNewton,sqrtGM1,a,mu,n,e
  real, intent(in) :: Qtidal1,klove1,m1,r1,O1e,O1q,O1h
  real, intent(in) :: Qtidal2,klove2,m2,r2,O2e,O2q,O2h
  real, intent(inout) :: XX1,XX2,YY1,YY2,ZZ1,ZZ2,VV1,VV2,WW1,WW2
  real :: X1,Y1,Z1,V1,W1,X2,Y2,Z2,V2,W2
!
  call XYZVW(Qtidal1,klove1,m1,m2,GNewton,sqrtGM1,r1,a,mu,n,e,O1e,O1q,O1h,X1,Y1,Z1,V1,W1)
  call XYZVW(Qtidal2,klove2,m2,m1,GNewton,sqrtGM1,r2,a,mu,n,e,O2e,O2q,O2h,X2,Y2,Z2,V2,W2)
!
  XX1 = XX1 + X1; XX2 = XX2 + X2
  YY1 = YY1 + Y1; YY2 = YY2 + Y2
  ZZ1 = ZZ1 + Z1; ZZ2 = ZZ2 + Z2
  VV1 = VV1 + V1; VV2 = VV2 + V2            
  WW1 = WW1 + W1; WW2 = WW2 + W2
!  
endsubroutine dissipation
!***********************************************************************
subroutine thirdbody(Mbody,Msun,Pout,pi,esun,HHH,&
           n,e,ev,qv,hv,See,Sqq,Shh,Sqh,Seh,Seq)
!            
! This subroutine calculates the coupling of the solar orbit to the
! binary orbit, using the model of Eggleton et al. (2001), accurate to
! the quadrupole term. 
!
! Sep-28-2019/wlyra
!
  real, intent(in) :: Mbody,Msun,Pout,pi,esun,n,e
  real, dimension(3), intent(in) :: HHH,ev,qv,hv
  real, intent(out) :: See,Sqq,Shh,Sqh,Seh,Seq
!  
  real :: Pin,tau,C
  real :: He,Hq,Hh
!
  Pin=2*pi/n
  tau = 2.*Pout**2/(3*pi*Pin) * (Mbody+Msun)/Msun * (1-esun**2)**1.5
  C = 1./(3*tau) * (1-e**2)**(-0.5)
!
  call get_eqh_components(HHH,ev,qv,hv,He,Hq,Hh)
!        
! Calculate the components of the Sij orbit coupling tensor. 
!
  See=C*(1-3*He**2)
  Sqq=C*(1-3*Hq**2)
  Shh=C*(1-3*Hh**2)
  Sqh=-3*C*Hq*Hh
  Seh=-3*C*He*Hh
  Seq=-3*C*He*Hq
!
endsubroutine thirdbody
!***********************************************************************
subroutine evolve_quantities(e,h,ev,hv,O1,O2,t,&
     dedt,dhdt,devdt,dhvdt,dO1dt,dO2dt,ds,dt_,levolve_spin)
!
! This subroutine adds the derivatives to the evolution quantities, 
! and advances the timestep. 
!
! Sep-28-2019/wlyra
!
  real, intent(in) :: dt_,ds,dedt,dhdt
  real, dimension(3), intent(in) :: devdt,dhvdt,dO1dt,dO2dt
  real, intent(inout) :: e,h,t
  real, dimension(3), intent(inout) :: ev,hv,O1,O2
  logical, intent(in) :: levolve_spin
  integer :: i 
!
!
!  Break if any quantity is a NaN
!  
  call check_NaN(t,'t')
  call check_NaN(dt_,'dt') 
  call check_NaN(h,'h')
  call check_NaN(e,'e')
  do i=1,3
    call check_NaN(ev(i),'ev')
    call check_NaN(hv(i),'hv')
    call check_NaN(O1(i),'O1')
    call check_NaN(O2(i),'O2')
  enddo
!
  e = e + dt_*dedt
  h = h + dt_*dhdt
  ev = ev + dt_*devdt
  hv = hv + dt_*dhvdt
!        
  if (levolve_spin) then 
    O1 = O1 + dt_*dO1dt    
    O2 = O2 + dt_*dO2dt
  endif
!                
  t = t + dt_*ds
!  
endsubroutine evolve_quantities
!***********************************************************************
subroutine check_NaN(q,str)

  real :: q
  character(len=*) :: str
  
  if (q /= q) then
    print*,'variable ',str,' has a NaN'
    stop
  endif
  
endsubroutine check_NaN
!***********************************************************************
!***********************************************************************          
program kozai_tidal
!
! Main program. This program solves the fourteen equations of the
! Kozai cycles plus tidal friction model, with the addition of a planetary 
! permanent J2, and our modest addition of gas drag.
!
! Sep-28-2019/wlyra
!
!***********************************************************************          
!
! Constants (in cgs) 
!
  real :: pi      = 3.141592653589793238
  real :: yr      = 3.15576008d7
  real :: GNewton = 6.67430d-8
  real :: Msun    = 1.98847d33
  real :: AU      = 1.49597871d13
  real :: km      = 1d5
  real :: hour    = 3600.
!
! Default parameters: MU69  
! 
  real :: x1km         = 20.6      ! km
  real :: y1km         = 19.9      ! km
  real :: z1km         = 9.4       ! km
  real :: x2km         = 15.4      ! km
  real :: y2km         = 13.8      ! km
  real :: z2km         = 9.8       ! km
  real :: rho_bullet   = 0.5       ! g/cm3
  real :: mu_body      = 4d10      ! rigidity 4e9 N/m2; 4e10 g/cm/s2
  real :: Qtidal1      = 100.
  real :: Qtidal2      = 100.
  real :: asun_AU      = 45.       ! AU
  real :: frac_hill    = 0.1
  real :: e            = 0.4
  real :: t            = 0.0
  real :: O1_hour      = 15.       ! hr
  real :: O2_hour      = 15.       ! hr
  real :: esun         = 0.0417249
  real :: I0_degree    = 95.       ! degree
  real :: omega        = 0.0
  real :: tau_drag_Myr = 10.       ! Myr
!  
! orbital vectors
!  
  real, dimension(3) :: ev=(/1.,0.,0./)
  real, dimension(3) :: qv=(/0.,1.,0./)  
  real, dimension(3) :: hv=(/0.,0.,1./)
  real, dimension(3) :: O1=(/0.,0.,0./)
  real, dimension(3) :: O2=(/0.,0.,0./)
  real, dimension(3) :: HHH =(/0.,0.,0./)
!
  real :: twopi
  real :: Myr,tmaxMyr=10,tmax
  real :: Rbody,Mbody,GM,GM1,sqrtGM,sqrtGM1
  real :: asun,rhill
  real :: J21,J22
  real :: a,h,n,mu,nsun,Hsun,Pout
  real :: II,sinI,cosI
  real :: beta,cosb,sinb
  real :: Hcte,tau_drag,tau1_drag
  real :: klove1,klove2,Qtidal,klove
  real :: x1_body,y1_body,z1_body,vol1,r1,I1,m1
  real :: x2_body,y2_body,z2_body,vol2,r2,I2,m2
!
  integer :: itmax     = 100000000
  integer :: itsave    = 1000000
  integer :: itdiagnos = 1000
  logical :: lJ2_tidal=.true.
  logical :: ldissipation=.true.
  logical :: lthirdbody=.true.
  logical :: lorbit_drag=.true.
  logical :: levolve_spin
  logical :: file_exists
  real    :: dt_period=1d-2
!  
! Runtime parameters
!
  integer :: it,itsub,j
  real :: ds=0.0
  real :: period,dt
  real :: dhdt,dedt
  real, dimension(3) :: devdt,dhvdt,dO1dt,dO2dt
  real :: O1e,O1q,O1h,O2e,O2q,O2h
  real :: Omega1,Omega2
  real :: XX1,XX2,YY1,YY2,ZZ1,ZZ2   
  real :: VV1,VV2,WW1,WW2
  real :: See,Sqq,Shh,Sqh,Seh,Seq
  real :: Wd
  real :: muhI1,muhI2
  real :: min_distance_km=0.0,min_distance  
!
! Parameters for Runge-Kutta
!
  real, dimension(3) :: alpha_ts
  real, dimension(3) :: beta_ts
  real, dimension(3) :: dt_beta_ts
!
! Namelist allows for input:
!
!  IO_degree        : Initial inclination, in degrees.
!  frac_hill        : Initial semimajor axis, in units of Hill radius (i.e, frac_hill=0.1 equals 10% of Hill radius).
!  e                : Initial eccentricity.
!  itmax            : Maximum number of iterations.
!  itsave           : Frequency that snapshots of dynamical variables (t,e,h,ev,hv,O1 and O2) are written to disk.
!  itdiagnos        : Frequency that the timeseries is written (itdiagnos=1000 means will be written every 1000 steps).
!  lthirdbody       : Include the kozai-inducing 3rd body.
!  ldissipation     : Include the effect of tidally-induced quadrupole.
!  lJ2_tidal        : Include the effect of J2 permanent quadrupole.
!  lorbit_drag      : Include orbital drag (1/tau exponential decay in angular momentum).
!  dt_period        : Timestep, given in fraction of the orbital period of the internal binary.
!  x1km             : Length in km of the principal axis (not semiaxis!) of the primary body.
!  y1km             : Length in km of the y-axis of the primary body.
!  z1km             : Length in km of the z-axis of the primary body.
!  x2km             : Length in km of the principal axis of the secondary body.
!  y2km             : Length in km of the y-axis of the secondary body.
!  z2km             : Length in km of the z-axis of the secondary body.
!  rho_bullet       : Density of the bodies (assumed equal).
!  mu_body          : Rigidity of the bodies (assumed equal).
!  Qtidal1          : Tidal dissipation quality of the primary.
!  Qtidal2          : Tidal dissipation quality of the secondary.
!  asun_AU          : Semimajor axis of the orbit of the center of mass of the binary around the 3rd body.
!  esun             : Eccentricity of the orbit of the center of mass of the binary around the 3rd body.
!  O1_hour          : Spin period of the primary, in hours.
!  O2_hour          : Spin period of the secondary, in hours.
!  tmaxMyr          : Maximum integration time, in millions of years.
!  min_distance_km  : The minimum distance (in km) that the bodies can come close; interrupts simulation.
!  
  namelist /input/ I0_degree, frac_hill, e, itmax, itdiagnos, itsave, lthirdbody, ldissipation, lJ2_tidal,&
       lorbit_drag, dt_period,x1km,y1km,z1km,x2km,y2km,z2km,rho_bullet,mu_body,Qtidal1,Qtidal2,asun_AU,esun,&
       O1_hour,O2_hour,tmaxMyr,min_distance_km
!
! Parameters for the time-integrator (3rd order Runge-Kutta scheme by Williamson, 1980).
! So far the only integrator used. Implement different ones in the future.   
!
  alpha_ts = (/0.   , -5./9.  ,-153./128./)
  beta_ts  = (/ 1./3., 15./16. ,   8./15. /)
!  
!  Read the input namelist with the user-defined parameters. 
!  
  open(20,file='input.in')
  read(20,nml=input)
  close(20)
!
!  Define some shorthands.
!
  twopi     = 2*pi
  Myr       = yr*1d6
  tmax      = tmaxMyr*Myr
  tau_drag  = tau_drag_Myr * Myr
  tau1_drag = 1.0/tau_drag
!
! Principal semiaxes of the ellipsoids, in cgs; volume, and their sphere-equivalent radii.
! 
  x1_body = x1km/2 * km
  y1_body = y1km/2 * km
  z1_body = z1km/2 * km
  x2_body = x2km/2 * km
  y2_body = y2km/2 * km
  z2_body = z2km/2 * km
!
  vol1 = 4./3 *pi*x1_body*y1_body*z1_body
  vol2 = 4./3 *pi*x2_body*y2_body*z2_body
!
  r1 = (x1_body*y1_body*z1_body)**(1./3)
  r2 = (x2_body*y2_body*z2_body)**(1./3)
!
! Set the minimum distance to interrupt the simulation.
! If not user-specified, use the sum of the radii of the bodies,
! which is the closest that the bodies can come to contact. 
!
  if (min_distance_km /= 0.0) then 
    min_distance = min_distance_km*km
  else
    min_distance = r1+r2
  endif
!
! Calculate the effective quadrupole J2 moment for a
! rotating homogeneous tri-axial ellipsoid (Scheeres 1994). 
!  
  J21 = 0.1 * (x1_body**2 + y1_body**2 - 2*z1_body**2) / r1**2
  J22 = 0.1 * (x2_body**2 + y2_body**2 - 2*z2_body**2) / r2**2
!
! Masses of the bodies, and reduced mass. 
!
  m1 = vol1 * rho_bullet
  m2 = vol2 * rho_bullet
  mu = m1*m2/(m1+m2)
!
! Inertia moments of the bodies. 
!
  I1 = 0.2*m1*(x1_body**2+y1_body**2)
  I2 = 0.2*m2*(x2_body**2+y2_body**2)
!
!  One-body equivalent
!
  Rbody   = (r1**3+r2**3)**(1./3)
  Mbody   = m1+m2
  GM      = GNewton*Mbody
  GM1     = 1./GM
  sqrtGM  = sqrt(GM)
  sqrtGM1 = sqrt(GM1)
!
!  Love number and tidal factor.
!  
  klove1 = 3./2 / (1. + 19.*mu_body*r1/(2.*GNewton*m1*rho_bullet))
  klove2 = 3./2 / (1. + 19.*mu_body*r2/(2.*GNewton*m2*rho_bullet))
  klove  = .5*(klove1+klove2)
  Qtidal = .5*(Qtidal1+Qtidal2)
!
! Parameters of the orbit around the 3rd body:
! semimajor axis, Hill radius, mean motion, angular momentum, 
! orbital period. These do not change during the integration.
!
  asun  = asun_AU*AU
  rhill = asun*((m1+m2)/(3*Msun))**(1./3)
  nsun  = sqrt(GNewton*Msun)/asun**1.5
  Hsun  = sqrt(GNewton*Msun*asun*(1-esun**2))
  Pout  = twopi/nsun
!  
! Parameters of the internal binary. 
! Initial semimajor axis, angular momentum,
! mean motion, and inclination
!
  a    = frac_hill*rhill
  n    = sqrtGM/a**1.5
  h    = sqrt(a*(1-e**2)*GM)
  II   = I0_degree * pi/180.
  sinI = sin(II)
  cosI = cos(II) 
!
! Spin period
!
  O1(3)  = twopi/(O1_hour*hour) 
  O2(3)  = twopi/(O2_hour*hour) 
  Omega1 = sqrt(sum(O1**2))
  Omega2 = sqrt(sum(O2**2))
  levolve_spin = lJ2_tidal.or.ldissipation
!
! Sanity check.
! 
  print*,''
  print*,'Equivalent spherical radii of the bodies (km)=',r1/km,r2/km
  print*,'Masses of the bodies (g) =',m1,m2
  print*,'Inertia moments about the spin axis =',I1,I2
  print*,'Love numbers =',klove1,klove2
  print*,'J2=',J21,J22
  print*,'Semimajor axis (km)=',a/km
  print*,'Mean motion (yr)=',twopi/n/yr
  print*,''
!
! Project the angular momentum of the external orbit into the axes defined by the orbit of the binary. 
!
  beta = 3*pi/2 - omega
  cosb = cos(beta)
  sinb = sin(beta)  
  HHH(1) = sinI*cosb
  HHH(2) = sinI*sinb
  HHH(3) = cosI
!
  inquire(file="snapshot.dat", exist=file_exists)
  if (file_exists) then
    print*,''
    print*,'re-starting from snapshot'
    open(40,file="snapshot.dat",status='old')
    read(40, FMT=*) t,e,h,&
           ev(1),ev(2),ev(3),&
           hv(1),hv(2),hv(3),&
           O1(1),O1(2),O1(3),&
           O2(1),O2(2),O2(3)
    close(40)
    call get_orbit(e,h,sqrtGM,GM,a,n)
    print*,'t (Myr)=',t/Myr
    print*,'e,h,a (km),n (yr)=',e,h,a/km,twopi/n/yr
    print*,'ev=',ev(1),ev(2),ev(3)
    print*,'hv=',hv(1),hv(2),hv(3)
    print*,'O1=',O1(1),O1(2),O1(3)
    print*,'O2=',O2(1),O2(2),O2(3)
    print*,''
  endif
!
  open(10,file="timeseries.dat",status="replace",action='write')
  print*,         'it -- t (Myr) -- dt(Myr) -- a(km) -- e -- q -- Q -- h -- i -- hconst -- Pbin (yr) -- Spin1(hr) -- Spin2(hr) '
  write(10,FMT=*) 'it -- t (Myr) -- dt(Myr) -- a(km) -- e -- q -- Q -- h -- i -- hconst -- Pbin (yr) -- Spin1(hr) -- Spin2(hr) '
!
  timestepping: do it=1,itmax
!
    period     = twopi*sqrtGM1*a**1.5
!
! Scale the timestep by the velocity at pericenter.
!
    dt         = dt_period*period*sqrt((1-e)/(1+e))
    dt_beta_ts = dt*beta_ts
!
    stageing: do itsub=1,3
      if (itsub == 1) then
        dhdt   = 0.0 ; dedt   = 0.0
        dhvdt  = 0.0 ; devdt  = 0.0
        dO1dt  = 0.0 ; dO2dt  = 0.0
        ds     = 0
      else
        dhdt   = alpha_ts(itsub)*dhdt ;  dedt  = alpha_ts(itsub)*dedt
        dhvdt  = alpha_ts(itsub)*dhvdt; devdt  = alpha_ts(itsub)*devdt
        dO1dt  = alpha_ts(itsub)*dO1dt; dO2dt  = alpha_ts(itsub)*dO2dt
        ds     = alpha_ts(itsub)*ds
      endif
!      
      ds=ds+1
!
!  Calculate the auxiliary variables semimajor axis and mean motion from the
!  dynamically evolved variables eccentricity and angular momentum.       
!
      call get_orbit(e,h,sqrtGM,GM,a,n)
!      
!  Get the auxiliary unit vector, \hat{q}, from the dynamically evolved unit vectors \hat{h} and \hat{e}.
!      
      call cross(hv,ev,qv)
!
!  If the spin is evolved, the spin-related variables (spin angular momenta and inertia moments)
!  are computed. Otherwise, don't waste time. 
!
      if (levolve_spin) then 
        call get_eqh_components(O1,ev,qv,hv,O1e,O1q,O1h)
        call get_eqh_components(O2,ev,qv,hv,O2e,O2q,O2h)      
        Omega1=sqrt(sum(O1**2))
        Omega2=sqrt(sum(O2**2))
        muhI1=mu*h/I1
        muhI2=mu*h/I2
      endif
!
!  Initialize the variables needed for the evolution equations. 
!
      XX1 = 0.; XX2 = 0.; YY1 = 0.; YY2 = 0.; ZZ1 = 0.; ZZ2 = 0.
      VV1 = 0.; VV2 = 0.; WW1 = 0.; WW2 = 0.; Wd  = 0.
      !See = 0.; Sqq = 0.; Shh = 0.; Sqh = 0.; Seh = 0.; Seq = 0.
!
!  Enter the KTJD model. K is for Kozai ...
!
      if (lthirdbody) call thirdbody(Mbody,Msun,Pout,pi,esun,HHH,&
           n,e,ev,qv,hv,See,Sqq,Shh,Sqh,Seh,Seq)
!
! T is for tides ...
!
      if (ldissipation) call dissipation(a,mu,n,e,&
           GNewton,sqrtGM1,&
           Qtidal1,klove1,m1,r1,O1e,O1q,O1h,&
           Qtidal2,klove2,m2,r2,O2e,O2q,O2h,&
           XX1,XX2,YY1,YY2,ZZ1,ZZ2,VV1,VV2,WW1,WW2)
!      
! J is for J2 ...
!
      if (lJ2_tidal) call tidalJ(a,n,e,&
           r1,O1e,O1q,O1h,Omega1,J21,&
           r2,O2e,O2q,O2h,Omega2,J22,&
           XX1,XX2,YY1,YY2,ZZ1,ZZ2)
!
!  and D is for drag (the last call).             
!
      if (lorbit_drag) Wd=tau1_drag
!
!  Evolution equations from Eggleton et al. (2001), modified by Fabrycky & Tremaine (2007) for the
!  planetary case, with the addition of permament J2 from Ragozzine & Brown (2009), and our modest
!  addition of gas drag.       
!
      dedt = dedt - e*(VV1 + VV2      + 5*(1-e**2)*Seq)
      dhdt = dhdt - h*(WW1 + WW2 + Wd - 5*   e**2 *Seq)
!        
      do j=1,3
        devdt(j) = devdt(j) + (ZZ1 + ZZ2 + (1-e**2)*(4*See-Sqq))*qv(j) - (YY1 + YY2 + (1-e**2)*Sqh)*hv(j)
        dhvdt(j) = dhvdt(j) + (YY1 + YY2 + (1-e**2)*Sqh)*ev(j) - (XX1 + XX2 + (4*e**2+1)*Seh)*qv(j)
!
        if (levolve_spin) then 
          dO1dt(j) = dO1dt(j) + muhI1 * (-YY1*ev(j) + XX1*qv(j) + WW1*hv(j)) 
          dO2dt(j) = dO2dt(j) + muhI2 * (-YY2*ev(j) + XX2*qv(j) + WW2*hv(j))
        endif
!        
      enddo
!
!  Add the derivatives to the variables and advance time. 
!      
      call evolve_quantities(e,h,ev,hv,O1,O2,t,dedt,dhdt,devdt,dhvdt,dO1dt,dO2dt,ds,dt_beta_ts(itsub),levolve_spin)
!
    enddo stageing
!    
!  On-the-fly diagnostic output: iteration, time in Myr, timestep in Myr, semimajor axis in km, eccentricity,
!  pericenter (km), apocenter (rhill), angular momentum in cgs, inclination (degrees), Kozai constant,
!  orbital period in years, and spin periods in hours.    
!
    if (mod(it,itdiagnos) == 0) then 
      call dot(HHH,hv,cosI)
      Hcte = sqrt(1-e**2)*cosI
      print*,             it,t/Myr, dt/Myr, a/km, e, a*(1-e)/km, a*(1+e)/rhill, h, acos(cosI)*180/pi, &
           Hcte, twopi/n/yr, twopi/Omega1/hour, twopi/Omega2/hour
      write(unit=10,FMT=*) it,t/Myr, dt/Myr, a/km, e, a*(1-e)/km, a*(1+e)/rhill, h, acos(cosI)*180/pi, &
           Hcte, twopi/n/yr, twopi/Omega1/hour, twopi/Omega2/hour
    endif

    if (mod(it,itsave) == 0) then
      open(30,file="snapshot.dat",status="replace",action='write')
      write(30,FMT=*) t,e,h,&
           ev(1),ev(2),ev(3),&
           hv(1),hv(2),hv(3),&
           O1(1),O1(2),O1(3),&
           O2(1),O2(2),O2(3)
      close(30)
    endif
!
!  Break the code if
!    (1) maximum number of iterations if reached,
!    (2) maximum integration time is reached,
!    (3) pericenter is smaller than the sum of the radii of the bodies. 
!
    if ((it == itmax) .or. (t > tmax) .or. (a*(1-e) < min_distance) .or. (a*(1+e) > rhill)) then
      print*,''
      print*,'Simulation finished.'
      print*,'Number of iterations it, and maximum allowed itmax=',it,itmax
      print*,'Time t (Myr), and maximum allowed tmax (Myr)=',t/Myr,tmax/Myr
      print*,'Pericenter q (km) =',a*(1-e)/km
      print*,'Apocenter Q (rhill) =',a*(1+e)/rhill
      print*,''
!
      if (it == itmax)            print*,"Maximum number of iterations reached."
      if (t > tmax)               print*,"Maximum integration time reached."
      if (a*(1-e) < min_distance) print*,"Pericenter smaller than minimum distance: bodies came into contact."
      if (a*(1+e) > rhill)        print*,"Apocenter larger than Hill radius: binary ionized."
!      
      close(10)
      stop
    endif
!
  enddo timestepping
!   
endprogram kozai_tidal
