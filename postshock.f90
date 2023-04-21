!$Id: postshock.f90,v 1.11 2007-08-14 12:05:58 bgiacoma Exp $

!!$ Copyright (C) 2005  B. Giacomazzo, L. Rezzolla

subroutine postshock(Vs,Byb,Bzb,vxb,ahead,behind)
  use type
  use global
  use wavecheck !contains wave_error
  use eos_param
  implicit none
  real(DP),intent(IN)::Vs,Byb,Bzb,vxb
  real(DP),intent(IN),dimension(7)::ahead !value of primitives ahead the shock
  real(DP),intent(OUT),dimension(7)::behind

  !Given the values of the post-shock By, Bz and v^x and the value of the shock velocity Vs,
  !returns the values of the remaining postshock quantities using the Bt-method.

  real(DP)::vyb,vzb,Db,pb,Wb,v2b,Vb,Ea,normBa
  real(DP)::rhoa,pa,vxa,vya,vza,Bya,Bza,etaa,b2a,ha,v2a,Wa,Ws,j,Da,taua,Va
  real(DP)::rhob,etab,b2b,F,B2


  wave_error=.false.

  !values ahead the shock
  rhoa = ahead(1) !mass density
  pa   = ahead(2) !Totale pressure (gas pressure+ magnetic pressure)
  vxa  = ahead(3) !velocity
  vya  = ahead(4)
  vza  = ahead(5)
  Bya  = ahead(6) !magnetic field
  Bza  = ahead(7)

  normBa=sqrt(Bya**2+Bza**2) !Bt

  etaa=Bx*vxa+Bya*vya+Bza*vza !B^j v_j

  v2a=vxa**2+vya**2+vza**2 !v^2
  if (v2a>=1.0e0_dp) stop 'error in postshock.f90: v^2a>c'
  Wa=1.0e0_dp/sqrt(1.0e0_dp-v2a) !Lorentz factor

  b2a=(Bx**2+Bya**2+Bza**2)/Wa**2+etaa**2 !(b_\mu b^\mu)

  !EoS
  call eos_enthalpy(pa-0.5e0_dp*b2a,rhoa,gamma,ha)

  !ha=1.0e0_dp+gamma*(pa-0.5e0_dp*b2a)/(rhoa*(gamma-1.0e0_dp)) !specifc gas enthalpy

  taua=(rhoa*ha+b2a)*Wa**2-pa-Da

  Da=rhoa*Wa

  Ea=taua-Wa**2*etaa**2

  Va=1/Da
  
  Ws=1.0e0_dp/sqrt(1-Vs**2) !Lorents factor of the shock wave

  j=Ws*rhoa*Wa*(Vs-vxa) !J


  F=j/Ws

  !-----------------------------------------------------!
  !--------------COMPUTE THE POST-SHOCK VALUES----------!
  !-----------------------------------------------------!
  B2=Bx**2+Byb**2+Bzb**2 !B^2

  Vb = Va + (vxa - vxb)/F

  
  if (Bx==0.0e0_dp) then
     stop 'riemann.f90: postshock: no solution for Bx=0 within the Bt-method'
  else
     vyb=(Bya*F*Va - Byb*(F*Va + vxa - vxb) + Bx*vya)/Bx
     vzb=(Bza*F*Va - Bzb*(F*Va + vxa - vxb) + Bx*vza)/Bx
  end if

  if (((vxb**2+vyb**2+vzb**2)>1.0e0_dp).OR.(abs(Vs)>1.0e0_dp)) then
     wave_error=.true.
     return
  end if

  Wb=1.0e0_dp/sqrt(1.0e0_dp-vxb**2-vyb**2-vzb**2)

  Db=1.0e0_dp/Vb

  rhob=Db*sqrt(1.0e0_dp-vxb**2-vyb**2-vzb**2)


  etab=vxb*Bx+vyb*Byb+vzb*Bzb

  b2b=B2/Wb**2+etab**2


  !This equation is obtained using the invariance of (iii) (see Anile, page 260)
  if (eos_ideal) then
     
     pb = (b2b*gamma*rhoa*Wa*(etab*j*Wb - Bx*rhob*Ws) + &
          2*(-1 + gamma)*rhob*(-(etab*j*rhoa*Wa*Wb) + &
          etaa*ha*j*rhob*Wa*Wb + Bx*rhoa*rhob*(Wa - ha*Wb)*Ws))/ &
          (2.*gamma*rhoa*Wa*(etab*j*Wb - Bx*rhob*Ws))

  elseif (eos_meliani) then
     
     pb = (b2b*rhoa*Wa*(etab*j*Wb - Bx*rhob*Ws)* &
          (j*(2*etaa*ha + etab*(-1 + gamma)*rhoa)*Wa*Wb - &
          Bx*rhoa*((-1 + gamma)*rhob*Wa + 2*ha*Wb)*Ws) - &
          2*(-1 + gamma)*(j*(etab*rhoa + etaa*ha*rhob)*Wa*Wb - &
          Bx*rhoa*rhob*(Wa + ha*Wb)*Ws)* &
          (etab*j*rhoa*Wa*Wb - &
          rhob*(etaa*ha*j*Wa*Wb + Bx*rhoa*(Wa - ha*Wb)*Ws)))/ &
          (2.*rhoa*Wa*(etab*j*Wb - Bx*rhob*Ws)* &
          (j*(2*etaa*ha + etab*(-1 + gamma)*rhoa)*Wa*Wb - &
          Bx*rhoa*((-1 + gamma)*rhob*Wa + 2*ha*Wb)*Ws))

  else
     write(*,*) 'postshock.f90:postshock: Error in the EOS'
     STOP
  end if

  behind(1)=rhob
  behind(2)=pb
  behind(3)=vxb
  behind(4)=vyb
  behind(5)=vzb
  behind(6)=Byb
  behind(7)=Bzb

end subroutine postshock



!---------------------------------------------------------------!
!-------------------------VELOCITY------------------------------!
!---------------------------------------------------------------!
subroutine velocity(unk1,unk2,ahead,switchLR,vx,Vs,switchPB,behind)
  use type
  use global
  use accmod !contains accuracy
  use wavecheck
  use interfaces,only:velocity_eqn,velocity_df,postshock
  implicit none
  real(DP),intent(IN)::unk1,unk2
  real(DP),dimension(7),intent(IN)::ahead
  character(len=2),intent(IN)::switchLR
  real(DP),intent(OUT)::vx,Vs
  character(len=1),intent(IN)::switchPB
  real(DP),dimension(7),intent(OUT),OPTIONAL::behind

  !given the post-shock value of By and Bz it computes 
  ! the postshock v^x and the shock velocity Vs
  !(It's used only in the Bt-method)

  real(DP)::vx1,vx2,f1,f2,vxl,vxh,dvx,dvxold,Vs1,Vs2,d,df,f,temp,ftol,step,vxinit
  real(DP)::xmin,xmax,step_min,cf,vx_max,vx_min
  integer(I4B)::j
  logical::upperbound,lowerbound
  real(DP),dimension(7)::plotstate

  !This routine is needed only if you use the norm of Bt
  if (switchPB=='P') stop 'postshock.f90: velocity: it has no sense to use this routine with the p-method'

  cf=1.0e0_dp !relaxation factor

  ftol=1.0e-10_dp !accuracy used in this subroutine

  if (accuracy<ftol) ftol=accuracy !to be consistent with the solution of the equations at the CD


  !It will search for the value of vx in the interval [xmin,xmax]

  vx=ahead(3)
  xmin=-0.9e0_dp
  xmax=0.9e0_dp


  select case(initial_data)
  case (13) 
     !For the generic Alfven test (from HLLE)
     if (switchLR=='LS') then
        vx=0.03885e0_dp
        xmin=0.0388e0_dp
        xmax=0.0395e0_dp
     else if (switchLR=='RS') then
        vx=0.03885e0_dp
        xmin=0.0388e0_dp
        xmax=0.0395e0_dp
     end if

  case (7)
     !For KO: Collision (from HLLE)
     if (switchLR=='LS') then
        vx=0.0e0_dp
        xmin=-1.0e-3_dp
        xmax=1.0e-3_dp
     else if (switchLR=='RS') then
        vx=0.0e0_dp
        xmin=-1.0e-3_dp
        xmax=1.0e-3_dp
     end if

  case (8) 
     !For Balsara 1 (from HLLE)
     if (switchLR=='LS') then
        vx=0.2555e0_dp
        xmin=0.25e0_dp
        xmax=0.26e0_dp
     else if (switchLR=='RS') then
        vx=0.25545e0_dp
        xmin=0.25540e0_dp
        xmax=0.25546e0_dp
     end if

  case (9)
     !For Balsara 2
     if (switchLR=='RS') then
        vx=0.676980015e0_dp
        xmin=0.675
        xmax=0.677
     end if

  case (10)
     !For Balsara 3 (from HLLE)
     if (switchLR=='RS') then
        vx=0.953e0_dp
        xmin=0.95e0_dp
        xmax=0.954e0_dp
     end if

  case(11)
     !For Balsara 4 (from HLLE)
     if (switchLR=='LS') then
        vx=0.0e0_dp
        xmin=-1.0e-4_dp
        xmax=5.0e-4_dp
     else if (switchLR=='RS') then
        vx=0.0e0_dp
        xmin=-2.0e-4_dp
        xmax=5.0e-4_dp
     end if

  case (12)
     !For Balsara 5 (from HLLE)
     if (switchLR=='RS') then
        vx=-0.454303238E-01
        xmin=-0.0455
        xmax=-0.0454
     end if

  case default
     !change these values if you are using your own initial condition
     ! (and if Bx is different from zero)
     if (switchLR=='LS') then
     vx=0.0e0_dp
     xmin=-0.99e0_dp
     xmax=0.99e0_dp
     else if (switchLR=='RS') then
     vx=0.0e0_dp
     xmin=-0.99e0_dp
     xmax=0.99e0_dp        
     end if

  end select


  vxinit=vx

  !step used to search the right interval in which the solution is contained
  ! (i.e. it must be f(xmin)*f(xmax)<0, and the function f must be defined in that interval)
  step=1.0e-1_dp*(xmax-xmin)
  step_min=1.0e-5*step 



  !-------------------------------------------------------------------------------------!
  !-----------------------ONLY FOR DEBUG------------------------------------------------!
  !---------------------------PLOT------------------------------------------------------!
  !-------------------------------------------------------------------------------------!

  !This plot the function whose roots give the solution for the velocity
  !Try this to get an initial guess about the interval in which the root is present

  !switchLR can be set equal to RS or LS
  if ((switchLR=="noplot").AND.(niter>-1)) then
     open(UNIT=1234,FILE='vel_func.dat',STATUS='REPLACE',ACTION='WRITE')

     print *
     print *,'vmax     = ',xmax
     print *,'vmin     = ',xmin
     print *,'switchLR = ',switchLR

     wave_error=.false.


     vx=xmin
     write(1234,*)

     do while (vx<=xmax)
        call velocity_eqn(unk1,unk2,ahead,vx,switchLR,f1,Vs1,plotstate)
        print *,vx,wave_error
        !if wave_error==true then f is not defined for that value
        if (.NOT.(wave_error)) write(1234,FMT='(13E30.18)') vx,f1,df,Vs1,unk1,unk2,plotstate
        wave_error=.false.
        vx=vx+1.0D-4*(xmax-xmin)
     end do
     close (1234)
     stop 'function plotted in vel_func.dat'
  end if

  !-------------------------------------------------------------------------------------!
  !-------------------------------------------------------------------------------------!
  !-------------------------------------------------------------------------------------!


    !-----------------------!
    !Bracketing the solution!
    !-----------------------!
2   if(upperbound .and. lowerbound) then
       step=step/2.0e0_dp
       if (step<step_min) then !STEP TOO SMALL
          wave_error=.true.
          stop 'error in velocity'
          !return
       end if
       if (verbose) then
          print *,'STEP=',step,'switchLR=',switchLR,'vx=',vx
          print *
       end if
       lowerbound=.false.
       upperbound=.false.
       vx=vxinit
    end if

    call velocity_eqn(unk1,unk2,ahead,vx,switchLR,f2,Vs)

    if ((abs(f2)<ftol).AND.(.NOT.(wave_error))) goto 1000

    if (wave_error) then
       vx=vx-step !Go to the left
       if (vx<xmin) then
          vx=xmax-(xmin-vx)
          if (lowerbound) then
             !reduce the step
             !It cannot find a valid interval
             upperbound=.true.
             goto 2
          end if
          lowerbound=.true.
       end if
       goto 2
    end if
    vx2=vx !In this point the function is defined 
    lowerbound=.false.
    upperbound=.false.
    if (verbose) then
       print *
       print *,'vx2=',vx2,'f2=',f2
    end if

    !Now start from vx2 and go to the left, if the function changes its sign
    !then we have found the right interval
    !else we have to reduce the step

    do
21     vx=vx-step !Go to the left
       if (vx<xmin) then
          wave_error=.true.
          goto 22
       end if

       call velocity_eqn(unk1,unk2,ahead,vx,switchLR,f1,Vs)

       if ((abs(f1)<ftol).AND.(.NOT.(wave_error))) goto 1000


22     if (wave_error) then !f is not defined go to the right
          vx1=vx
          vx=vx2
4         vx=vx+step
          if (vx>xmax) then
             wave_error=.true.
             goto 5
          end if
          call velocity_eqn(unk1,unk2,ahead,vx,switchLR,f1,Vs)

          if ((abs(f1)<ftol).AND.(.NOT.(wave_error))) goto 1000


5         if (wave_error) then 
             step=step/2.0e0_dp
             if (step<step_min) then !STEP TOO SMALL
                wave_error=.true.
                stop 'error in velocity'                
                return
             end if

             if (verbose) then
                print *,'reducing the step and tring again...'
                print *,'STEP=',step
                print *
             end if
             vx=vx2
             goto 21
          end if
          if ((f1*f2)<0.0e0_dp) then
             vx1=vx
             goto 10
          else
             goto 4
          end if
       end if
       if ((f1*f2)<0.0e0_dp) then
          vx1=vx
          goto 10
       end if
    end do



10  if (verbose) then
     print *
     print *,'Byb=',unk1,'Bzb=',unk2
     print *,'I will try to find the solution in the following interval:'
     print *,'vx_min=',vx1,' f1=',f1
     print *,'vx_max=',vx2,' f2=',f2
     print *
  end if


  !Find the post-shock velocity
  !Newton-Raphson and bisection method

  if (abs(f1)<=ftol) then
     if (veryverbose) print *,'postshock.f90: pressure: f1 is zero'
     vx=vx1
     Vs=Vs1
     if (present(behind)) call postshock(Vs,unk1,unk2,vx,ahead,behind)
     return
  end if

  if (abs(f2)<=ftol) then
     if (veryverbose) print *,'postshock.f90: pressure: f2 is zero'
     vx=vx2
     Vs=Vs2
     if (present(behind)) call postshock(Vs,unk1,unk2,vx,ahead,behind)
     return
  end if

  if (f1*f2>0.0e0_dp) stop 'f1 and f2 have the same sign; there is no solution!'

  if (f1 < 0.0) then !Orient the search so that f(Vsl)<0. 
     vxl=vx1 
     vxh=vx2 
  else 
     vxh=vx1 
     vxl=vx2 
  end if

  vx=0.5e0_dp*(vx1+vx2) !Initialize the guess for root
  dvxold=abs(vx2-vx1) !the "stepsize before last",  
  dvx=dvxold !and the last step. 

  call velocity_df(unk1,unk2,ahead,vx,switchLR,f,df,Vs)

  do j=1,200 !200 is the maximum number of iterations

     if (((vx-vxh)*df-f)*((vx-vxl)*df-f) > 0.0 .or. & 
          abs(2.0_DP*f) > abs(dvxold*df) ) then 
        !Bisect if Newton out of range, or not decreasing fast enough. 
        dvxold=dvx 
        dvx=0.5_DP*(vxh-vxl) 
        vx=vxl+dvx 
        if (vxl == vx) goto 1000 !Change in root is negligible. 
     else !Newton step acceptable. Take it. 
        dvxold=dvx 
        dvx=f/df 
        temp=vx 
        vx=vx-dvx 
        if (temp == vx) goto 1000
     end if
     if (abs(dvx) < ftol) goto 1000 !Convergence criterion. 
     call velocity_df(unk1,unk2,ahead,vx,switchLR,f,df,Vs)
     !One new function evaluation per iteration. 

     if (f < 0.0) then !Maintain the bracket on the root. 
        vxl=vx 
     else 
        vxh=vx
     end if
  end do
  stop 'I cannot find post-shock velocity'

999 call velocity_eqn(unk1,unk2,ahead,vx,switchLR,f1,Vs)
  if (wave_error) stop 'postshock.f90: velocity: no solution for the value vx used'

1000 if (present(behind)) call postshock(Vs,unk1,unk2,vx,ahead,behind)


end subroutine velocity



!---------------------------------------------------------------!
!-------------------------VELOCITY_DF---------------------------!
!---------------------------------------------------------------!
subroutine velocity_df(unk1,unk2,ahead,vx,switchLR,f,df,Vs)
  use type
  use global
  use wavecheck
  use interfaces,only:velocity_eqn
  implicit none
  real(DP),intent(IN)::unk1,unk2,vx
  real(DP),dimension(7),intent(IN)::ahead
  character(len=2),intent(IN)::switchLR
  real(DP),intent(OUT)::f,df
  real(DP),intent(OUT),OPTIONAL::Vs

  !Compute the first derivative of f(vx)

  real(DP)::delta,f2,dump

  if (veryverbose) then
     print *
     print *,'postshock.f90: velocity_df: I''m computing the derivative of f(vx)'
     print *
  end if

  call velocity_eqn(unk1,unk2,ahead,vx,switchLR,f,Vs)
  if (wave_error) return 
  
  !First derivative
  delta=1.0D-8*vx
  if (delta<=1.0D-10) delta=1.0D-10

  call velocity_eqn(unk1,unk2,ahead,vx-delta,switchLR,f2,dump)  
  if (wave_error) then
     print *,'vx=',vx,'delta=',delta,'switchLR=',switchLR
     stop 'postshock.f90: velocity_df: error computing the derivative'
  end if
  
  df=(f-f2)/delta
  
end subroutine velocity_df


!---------------------------------------------------------------!
!-------------------------VELOCITY_EQN--------------------------!
!---------------------------------------------------------------!

subroutine velocity_eqn(unk1,unk2,ahead,vxb,switchLR,f,Vs,stateb)
  use type
  use global
  use wavecheck!This module contains wave_error
  use interfaces,only:ContactVelocity,shockvelocity,postshock
  implicit none
  real(DP),intent(IN)::unk1,unk2
  real(DP),dimension(7),intent(IN)::ahead
  real(DP),intent(IN)::vxb
  character(len=2),intent(IN)::switchLR
  real(DP),intent(OUT)::f
  real(DP),intent(OUT)::Vs
  real(DP),dimension(7),OPTIONAL,intent(OUT)::stateb

  !The zero of the equation (4.16) gives the right value for vxb
  !Used only within the Bt-method

  real(DP)::rhoa,pa,vxa,vya,vza,Bya,Bza,etaa,v2a,Wa,b2a,ha,taua,Da,Ws,j
  real(DP)::rhob,pb,vyb,vzb,Byb,Bzb,etab,v2b,Wb,b2b,hb,taub,Db
  real(DP)::Ea,Eb,dump,Va,Vb,normBa,B2,bigF
  real(DP),dimension(7)::behind

  Byb=unk1
  Bzb=unk2

  rhoa = ahead(1) !mass density
  pa   = ahead(2) !Totale pressure (gas pressure+ magnetic pressure)
  vxa  = ahead(3) !velocity
  vya  = ahead(4)
  vza  = ahead(5)
  Bya  = ahead(6) !magnetic field
  Bza  = ahead(7)

  if (abs(vxa-vxb)<1.0e-15_dp) then
     !It's not a shock
     wave_error=.true.
     return
  end if


  normBa=sqrt(Bya**2+Bza**2) !Bt ahead the shock

  etaa=Bx*vxa+Bya*vya+Bza*vza !B^j v_j

  v2a=vxa**2+vya**2+vza**2 !v^2
  if (v2a>=1.0e0_dp) stop 'error in postshock.f90: pressure_eqn: v^2a>c'
  Wa=1.0e0_dp/sqrt(1.0e0_dp-v2a) !Lorentz factor

  b2a=(Bx**2+Bya**2+Bza**2)/Wa**2+etaa**2 !(b_\mu b^\mu)

  !EoS
  call eos_enthalpy(pa-0.5e0_dp*b2a,rhoa,gamma,ha)
  !ha=1.0e0_dp+gamma*(pa-0.5e0_dp*b2a)/(rhoa*(gamma-1.0e0_dp)) !specific gas enthalpy

  Da=rhoa*Wa

  Va=1.0e0_dp/Da

  taua=(rhoa*ha+b2a)*Wa**2-pa-Da

  Ea=taua-Wa**2*etaa**2


  !Compute the shock velocity from normB, vxb and the postshock values
  wave_error=.false.
  call ContactVelocity(Byb,Bzb,vxb,Vs,ahead,dump,switchLR,'B',behind)
  if (wave_error) then
     if (veryverbose) print *,'postshock.f90: velocity_eqn: wave_error'
     return
  end if
  Ws=1.0e0_dp/sqrt(1-Vs**2) !Lorentz factor for the shock
  j=Ws*rhoa*Wa*(Vs-vxa) !J

  bigF=j/Ws

  rhob = behind(1) !mass density
  pb   = behind(2) !total pressure
  vyb  = behind(4) !tangential components of velocity
  vzb  = behind(5)
  Byb  = behind(6) !magnetic field
  Bzb  = behind(7)

  if (present(stateb)) stateb=behind


  B2=Bx**2+Byb**2+Bzb**2

  etab=Bx*vxb+Byb*vyb+Bzb*vzb

  v2b=vxb**2+vyb**2+vzb**2

  if (v2b>=1.0e0_dp) then
     if (verbose) print *,'velocity_eqn: v2b=',v2b
     wave_error=.true.
     return
  end if

  Wb=1.0e0_dp/sqrt(1.0e0_dp-v2b)

  !Compare j before and after the shock (it's only a check)
  if (veryverbose) print *,'ja=',j,'jb=',Ws*rhob*Wb*(Vs-vxb)
  

  b2b=(Bx**2+Byb**2+Bzb**2)/Wb**2+etab**2 !this depends on normB

  !EoS
  call eos_enthalpy(pb-0.5e0_dp*b2b,rhob,gamma,hb)
  !hb=1.0e0_dp+gamma*(pb-0.5e0_dp*b2b)/(rhob*(gamma-1.0e0_dp))

  Db=rhob*Wb

  Vb=1.0e0_dp/Db

  taub=(rhob*hb+b2b)*Wb**2-pb-Db

  Eb=taub-Wb**2*etab**2

  !equation (4.16)
  f = bigF*(((Ea+pa+Da)*vya-etaa*Bya)/Da-((Eb+pb+Db)*vyb-etab*Byb)/Db) + Bx*(etaa*vya-etab*vyb)+ &
       Bx*(Bya/Wa**2-Byb/Wb**2)

  if (.NOT.(abs(f)>=0.0e0_dp)) then
     if (verbose) print *,'velocity_eqn: f=',f
     wave_error=.true.
     return
  end if


end subroutine velocity_eqn



!---------------------------------------------------------------!
!---------------------CONTACT VELOCITY--------------------------!
!---------------------------------------------------------------!

subroutine ContactVelocity(unk1,unk2,vxb,Vs,ahead,vx,switchLR,switchPB,solution)
  use type
  use accmod !This module contains accuracy
  use global !This module contains Bx and gamma
  use wavecheck!This module contains wave_error
  use odeswitch!This module contains switch
  use interfaces,only: zeroeqn,postshock,eos_cs2,eos_enthalpy,postshock_pressure
  implicit none
  real(DP),intent(IN)::unk1,unk2,vxb
  real(DP),intent(OUT)::Vs
  real(DP),dimension(7),intent(IN)::ahead
  real(DP),intent(OUT)::vx
  character(len=2),intent(IN)::switchLR
  character(len=1),intent(IN)::switchPB
  real(DP),dimension(7),intent(OUT),OPTIONAL::solution

  !Given vxb this subroutine returns Vs and
  !the value of the postshock vx

  !It can be used with both the Bt-method (switchPB=B) and the p-method (switchPB=P)

  integer(I4B)::j
  real(DP)::dVs,f,df,Vsinit,temp,Vs1,Vs2,f1,f2
  real(DP),dimension(7)::behind
  logical::errcheck,secondtry,lowerbound,upperbound
  real(DP)::Vsl,Vsh,dVsold,step
  real(DP)::b2,pgas,wtot,leftAlfven,rightAlfven,cdvel,xmax,xmin,B2big,W2inv,eta,enthalpy

  real(DP)::step_min,ftol,cf !Minimum step allowed, accuracy and relaxation factor

  ftol=1.0e-10_dp !accuracy used in this subroutine

  if (accuracy<ftol) ftol=accuracy !to be consistent with the solution of the equations at the CD

  cf=0.5e0_dp

  wave_error=.false.

  switch=switchLR


  !Compute Alfven velocities (left and right-going)
  B2big=Bx**2+ahead(6)**2+ahead(7)**2 !B^2
  W2inv=1.0e0_dp-ahead(3)**2-ahead(4)**2-ahead(5)**2
  eta=Bx*ahead(3)+ahead(6)*ahead(4)+ahead(7)*ahead(5) !B_j v^j
  b2=B2big*W2inv+eta**2 !b_mu b^\mu
  pgas=ahead(2)-0.5e0_dp*b2 !gas pressure
  call eos_enthalpy(pgas,ahead(1),gamma,enthalpy)
  wtot=ahead(1)*enthalpy + b2 !total relativistic enthalpy
  leftAlfven =ahead(3)+Bx*W2inv/(eta-sqrt(wtot)) !left-going Alfven velocity
  rightAlfven=ahead(3)+Bx*W2inv/(eta+sqrt(wtot)) !right-going Alfven velocity
  !vx
  if ((switchLR=='LF').OR.(switchLR=='LS'))then
     cdvel=min(ahead(3),vxb) !because J should be always < 0
     if (switchPB=='P') cdvel=ahead(3)
  else
     cdvel=max(ahead(3),vxb) !because J should be always > 0
     if (switchPB=='P') cdvel=ahead(3)
  end if

  if (switchPB=='B') then
     if(abs(vxb-ahead(3))<epsilon(vxb)) then !No change in v^x (there is no discontinuity)
        if (present(solution)) solution=ahead
        vx=ahead(3)
        if ((switchLR=='LS').OR.(switchLR=='RS')) then
           Vs=vx
        else if (switchLR=='LF') then
           Vs=leftAlfven
        else if (switchLR=='RF') then
           Vs=rightAlfven
        else
           stop 'riemann.f90: ContactVelocity: error in switchLR'
        end if
        return
     end if
  end if


  if (switchLR=='LF') then
     xmax=leftAlfven
     xmin=-0.9999999e0_dp
     !if (initial_data==0) xmin=0.99 !for user defined Riemann problem (if needed)
  else if (switchLR=='LS') then
     xmax=cdvel
     xmin=leftAlfven

     !This is due to the presence of more the one root in the interval (xmin,xmax)
     if (initial_data==13) xmin=-0.2e0_dp   !for test alfven 1
     if (initial_data==8)  xmin=0.02e0_dp   !for Balsara 1
     if (initial_data==11) xmin=-0.150e0_dp !for Balsara 4 

     !if (initial_data==0) xmin=???   !for user defined Riemann problem (if needed)

  else if (switchLR=='RS') then
     xmax=rightAlfven
     xmin=cdvel

     !This is due to the presence of more the one root in the interval (xmin,xmax)     
     if (initial_data==9)  xmin=0.8e0_dp  !for Balsara 2
     if (initial_data==8)  xmin=0.38e0_dp !for Balsara 1
     if (initial_data==11) xmax=0.16e0_dp !for Balsara 4
     if (initial_data==12) xmin=0.36e0_dp !for Balsara 5
     !if (initial_data==10) xmax=0.995e0_dp !for Balsara 4 !Only for Meliani EOS
     !if (initial_data==0) xmin=???   !for user defined Riemann problem (if needed)

  else if (switchLR=='RF') then
     xmax=0.9999999e0_dp
     xmin=rightAlfven
     !if (initial_data==0) xmin=0.999999   !for user defined Riemann problem (if needed)
     !This is due to the presence of more than one root in the interval (xmin,xmax)
     if (initial_data==10)  xmin=0.994e0_dp !for Balsara 3
  else
     stop 'riemann.f90: ContactVelocity: error in switchLR'
  end if


  call zeroeqn(vxb,unk1,unk2,xmax,ahead,f,df,switchLR,switchPB,errcheck)
  if (veryverbose) print *,'ContactVelocity: f(xmax)=',f,'xmax=',xmax,'errcheck=',errcheck
  if ((abs(f)<=ftol).AND.(.NOT.(errcheck))) then
     Vs=xmax
     write(*,*) vxb,unk1,unk2,xmax,ahead,f,df,switchLR,switchPB,errcheck
     stop 'postshock.f90: contact_velocity: Vs=vmax; this is strange!'
     if (veryverbose) print *,'Vs=',Vs
     return
  end if
  call zeroeqn(vxb,unk1,unk2,xmin,ahead,f,df,switchLR,switchPB,errcheck)
  if (veryverbose) print *,'ContactVelocity: f(xmin)=',f ,'xmin=',xmin,'errcheck=',errcheck
  if ((abs(f)<=ftol).AND.(.NOT.(errcheck))) then
     Vs=xmin
     write(*,*) vxb,unk1,unk2,xmin,ahead,f,df,switchLR,switchPB,errcheck
     stop 'postshock.f90: contact_velocity: Vs=vmin; this is strange!'
     if (veryverbose) print *,'Vs=',Vs
     return
  end if

  secondtry=.false.


  Vs=0.5e0_dp*(xmax+xmin)!First guess

  Vsinit=Vs


  step=1.0D-3*abs(xmax-xmin) !Initial step

  step_min=1.0e-4_dp*step !Minimum step allowed

  upperbound=.false.
  lowerbound=.false.


  !-----------------------ONLY FOR DEBUG-----------------------------!
  !--------------------------PLOT------------------------------------!

  !This plots the function whose roots give the value of the shock velocity
  !This can be used to have an idea of what is going wrong if the routine doesn't work and
  ! also to check for the presence of more than one root in the interval

  !substitute noplot with one of the following: RF, RS, LF, LS 
  ! (respectively right-fast, right-slow, left-fast, left-slow shock)
  if (switchLR=='noplot') then
     open(UNIT=123,FILE='funcd.dat',STATUS='REPLACE',ACTION='WRITE')

     print *
     print *,'cdvel =',cdvel
     print *,'xmax  =',xmax
     print *,'xmin  =',xmin


!!$     xmin=-0.99999e0_dp
!!$     xmax=0.99999e0_dp

     errcheck=.false.
     Vs=xmin

     write(123,*)

     do while (Vs<=xmax)
        call zeroeqn(vxb,unk1,unk2,Vs,ahead,f,df,switchLR,switchPB,errcheck)
        !if errcheck==true then f is not defined for that value
        if (.NOT.(errcheck)) write(123,FMT='(3E30.18)') Vs,f,vxb
!!$        write(*,*) Vs,f
        errcheck=.false.
        Vs=Vs+1.0D-4*(xmax-xmin)
     end do
     close(123)
     stop 'function plotted in funcd.dat'
  end if
  !------------------------------------------------------------------!
  !------------------------------------------------------------------!


  Vs=Vsinit

  !-----------------------!
  !Bracketing the solution!
  !-----------------------!
2 if(upperbound .and. lowerbound) then
     step=step/2.0e0_dp
     if (step<step_min) then !STEP TOO SMALL
        wave_error=.true.
        !rarefaction_pmin=ahead(2)
        shock=.true.
        return
     end if
     if (veryverbose) then
        print *,'STEP=',step,'switchLR=',switchLR,'vxb=',vxb
        print *
     end if
     lowerbound=.false.
     upperbound=.false.
     Vs=Vsinit
  end if

  wave_error=.false.


  call zeroeqn(vxb,unk1,unk2,Vs,ahead,f2,df,switchLR,switchPB,errcheck)

  if ((abs(f2)<ftol).AND.(.NOT.(errcheck))) goto 1000

  if (errcheck) then
     Vs=Vs-step !Go to the left
     if (Vs<xmin) then
        Vs=xmax-(xmin-Vs)
        if (lowerbound) then
           !reduce the step
           !It cannot find a valid interval
           upperbound=.true.
           goto 2
        end if
        lowerbound=.true.
     end if
     goto 2
  end if
  Vs2=Vs !In this point the function is defined 
  lowerbound=.false.
  upperbound=.false.
  if (veryverbose) then
     print *
     print *,'Vs2=',Vs2,'f2=',f2,'vxb=',vxb
  end if

  !Now start from Vs2 and go to the left, if the function changes its sign
  !then we have found the right interval
  !else we have to change interval

  do
     Vs=Vs-step !Go to the left
     if (Vs<xmin) then
        errcheck=.true.
        goto 22
     end if
     call zeroeqn(vxb,unk1,unk2,Vs,ahead,f1,df,switchLR,switchPB,errcheck)

     if ((abs(f1)<ftol).AND.(.NOT.(errcheck))) goto 1000


22   if (errcheck) then !f is not defined go to the right
        Vs1=Vs
        Vs=Vs2
4       Vs=Vs+step
        if (Vs>xmax) then
           errcheck=.true.
           goto 5
        end if
        call zeroeqn(vxb,unk1,unk2,Vs,ahead,f1,df,switchLR,switchPB,errcheck)

        if ((abs(f1)<ftol).AND.(.NOT.(errcheck))) goto 1000


5       if (errcheck) then !wrong interval. chose a different starting point
           step=step/2.0e0_dp
           if (step<step_min) then !STEP TOO SMALL
              wave_error=.true.
              !rarefaction_pmin=ahead(2)
              shock=.true.
              return
           end if

           if (veryverbose) then
              print *,'wrong interval chose a different starting point'
              print *,'STEP=',step
              print *
           end if
           Vs=Vs1-step
           if (Vs<xmin) Vs=xmax-(xmin-Vs)
           goto 2
        end if
        if ((f1*f2)<0.0e0_dp) then
           Vs1=Vs
           goto 10
        else
           goto 4
        end if
     end if
     if ((f1*f2)<0.0e0_dp) then
        Vs1=Vs
        goto 10
     end if
  end do


  !Find the shock velocity
  !Newton-Raphson and bisection method
10 if (veryverbose) then
     print *,'Finding Root for'
     print *,'vxb=',vxb
     print *,'in the following interval'
     print *,'Vs1=',Vs1,'f1=',f1
     print *,'Vs2=',Vs2,'f2=',f2
  end if

  if (f1==0.0e0_dp) then
     Vs=Vs1
     stop 'Vs=Vs1; seems strange!'
     goto 1000 !return
  end if

  if (f2==0.0e0_dp) then
     Vs=Vs2
     stop 'Vs=Vs2; seems strange!'
     goto 1000 !return
  end if

  if (f1 < 0.0) then !Orient the search so that f(Vsl)<0. 
     Vsl=Vs1 
     Vsh=Vs2 
  else 
     Vsh=Vs1 
     Vsl=Vs2 
  end if

  Vs=0.5e0_dp*(Vs1+Vs2) !Initialize the guess for root
  dVsold=abs(Vs2-Vs1) !the "stepsize before last",  
  dVs=dVsold !and the last step. 

  !This returns f(Vs) and df/dVs
  call zeroeqn(vxb,unk1,unk2,Vs,ahead,f,df,switchLR,switchPB,errcheck)


  if (errcheck) then
     wave_error=.true.
     print *,'something wrong: a point in the interval returns an error!'
     stop 'postshock.f90: ContactVelocity: something wrong 1'
  end if

  do j=1,200 !200 is the maximum number of iterations

     if (((Vs-Vsh)*df-f)*((Vs-Vsl)*df-f) > 0.0 .or. & 
          abs(2.0_DP*f) > abs(dVsold*df) ) then 
        !Bisect if Newton out of range, or not decreasing fast enough. 
        dVsold=dVs 
        dVs=0.5_DP*(Vsh-Vsl) 
        Vs=Vsl+dVs 
        if (Vsl == Vs) goto 1000 !Change in root is negligible. 
     else !Newton step acceptable. Take it. 
        dVsold=dVs 
        dVs=f/df 
        temp=Vs 
        Vs=Vs-cf*dVs !use a relaxation factor 
        if (temp == Vs) goto 1000
     end if
     if (abs(dVs) < ftol) goto 1000 !Convergence criterion. 

113  call zeroeqn(vxb,unk1,unk2,Vs,ahead,f,df,switchLR,switchPB,errcheck)



     !One new function evaluation per iteration. 
     if (errcheck) then
        wave_error=.true.
        print *,'something wrong: a point in the interval returns an error!'
        stop 'postshock.f90: ContactVelocity: something wrong 2'
     end if
     if (f < 0.0) then !Maintain the bracket on the root. 
        Vsl=Vs 
     else 
        Vsh=Vs
     end if
  end do
  stop 'I cannot find shock velocity'


1000 if (switchPB=='B') then
     !we are using the Bt-method
     call postshock(Vs,unk1,unk2,vxb,ahead,behind)
  else
     !we are using the p-method
     call postshock_pressure(Vs,unk1,ahead,behind,.false.)
  end if

  if (wave_error) stop 'postshock.f90: error in postshock called by ContactVelocity'
  if (veryverbose) then
     print *
     print *, 'This is the value I found:'
     print *,'Vs=',Vs
     print *,'post-shock values=',behind
     print *,'exiting from ContactVelocity'
     print *
  end if
  vx=behind(3)

  if (present(solution)) solution=behind


end subroutine ContactVelocity


!---------------------------------------------------------------!
!--------------------------SHOCKFUNC----------------------------!
!---------------------------------------------------------------!


subroutine shockfunc(Vs,unk1,unk2,vxbin,ahead,f,switchLR,switchPB,errcheck) 
  USE type
  use global !this module contains Bx and gamma
  use wavecheck !contains wave_error
  use interfaces,only: postshock, postshock_pressure
  IMPLICIT NONE 
  REAL(DP),INTENT(IN) :: Vs,unk1,unk2,vxbin
  real(DP),dimension(7),intent(IN) ::ahead
  REAL(DP),intent(OUT):: f
  character(len=2),intent(IN)::switchLR
  character(len=1),intent(IN)::switchPB
  logical,intent(OUT)::errcheck

  !This subroutine computes equation (4.26)

  real(DP)::rhob,pb,vxb,vyb,vzb,Byb,Bzb,v2b,Wb,etab,b2b,Ws,j,hb
  real(DP)::rhoa,vxa,vya,vza,Bya,Bza,v2a,Wa,etaa,b2a,ha,pa
  real(DP),dimension(7)::behind

  errcheck=.false.

  !STATE AHEAD
  rhoa=ahead(1)
  pa=ahead(2)
  vxa=ahead(3)
  vya=ahead(4)
  vza=ahead(5)
  Bya=ahead(6)
  Bza=ahead(7)

  v2a=vxa**2+vya**2+vza**2 !v^2
  Wa=1.0e0_dp/sqrt(1.0e0_dp-v2a) !Lorentz factor
  etaa=vxa*Bx+vya*Bya+vza*Bza !B_j v^j
  b2a=(Bx**2+Bya**2+Bza**2)/Wa**2+etaa**2 !b_mu b^mu

  !EoS
  call eos_enthalpy(pa-0.5e0_dp*b2a,rhoa,gamma,ha)
  !ha=1.0e0_dp+gamma*(pa-0.5e0_dp*b2a)/(rhoa*(gamma-1.0e0_dp)) !specific gas enthalpy

  !STATE BEHIND

  if (((Vs-vxb)==0.0e0_dp).OR.((Vs-vxa)==0.0e0_dp)) then
     errcheck=.true.
     f=-1.0d0
     f=sqrt(f)!NaN
     return
  else
     if (switchPB=='B') then
        call postshock(Vs,unk1,unk2,vxbin,ahead,behind)
        if (wave_error) then
           errcheck=.true.
           f=-1.0d0
           f=sqrt(f)!NaN
           return
        end if
     else
        call postshock_pressure(Vs,unk1,ahead,behind,.false.)
     end if
  end if

  rhob = behind(1)
  pb   = behind(2)
  vxb  = behind(3)
  vyb  = behind(4)
  vzb  = behind(5)
  Byb  = behind(6)
  Bzb  = behind(7)

  v2b=vxb**2+vyb**2+vzb**2 !v^2

  if ((v2b>1.0e0_dp).OR.(pb<=0.0e0_dp)) then 
     errcheck=.true.
     f=-1.0d0
     f=sqrt(f)!NaN
     return
  end if

  Wb=1.0e0_dp/sqrt(1.0e0_dp-v2b) !Lorentz factor
  etab=vxb*Bx+vyb*Byb+vzb*Bzb !B_j v^j
  b2b=(Bx**2+Byb**2+Bzb**2)/Wb**2+etab**2 !b_mu b^mu

  !EoS
  call eos_enthalpy(pb-0.5e0_dp*b2b,rhob,gamma,hb)
  !hb=1.0e0_dp+gamma*(pb-0.5e0_dp*b2b)/(rhob*(gamma-1.0e0_dp)) !specific gas enthalpy


  Ws=1.0e0_dp/sqrt(1.0e0_dp-Vs**2) !Lorentz factor of the shock
  j=rhob*Wb*Ws*(Vs-vxb) !J

  !equation (4.26)
  f=(pa-pb)+j**2*(ha/rhoa-hb/rhob)


  if ((hb<1.0e0_dp).OR.(rhob<0.0e0_dp).OR.((pb-0.5e0_dp*b2b)<0.0e0_dp)) then 
     !the results are not physical
     !the function is considered as not defined for this value of Vs
     errcheck=.true.
     return
  end if

  if(.NOT.(abs(f)>=0)) then
     !NaN
     errcheck=.true.
     return
  end if

END subroutine shockfunc



!---------------------------------------------------------------!
!-------------------------ZEROEQN-------------------------------!
!---------------------------------------------------------------!


SUBROUTINE zeroeqn(vxb,unk1,unk2,Vs,ahead,f,df,switchLR,switchPB,errcheck) 
  USE type
  use interfaces, only: shockfunc
  IMPLICIT NONE 
  REAL(DP),INTENT(IN) :: vxb,unk1,unk2,Vs
  real(DP),dimension(7),intent(IN)::ahead
  REAL(DP),INTENT(OUT) :: f 
  REAL(DP),INTENT(OUT) :: df 
  character(len=2),intent(IN)::switchLR
  character(len=1),intent(IN)::switchPB
  logical,intent(OUT)::errcheck

  !This subroutine is called by contact velocity and 
  ! compute the value of f(Vs) and its first derivative
  ![ f(Vs)=0 --> exact value of Vs ]

  real(DP)::delta,f2

  errcheck=.false.

  call shockfunc(Vs,unk1,unk2,vxb,ahead,f,switchLR,switchPB,errcheck)
5 if (errcheck) return

  !First derivative
  delta=1.0D-5*(1.0e0_dp-abs(Vs)) !avoid abs(Vs+delta)>1
  if (delta==0.0e0_dp) delta=1.0D-5

10 call shockfunc(Vs+delta,unk1,unk2,vxb,ahead,f2,switchLR,switchPB,errcheck)

  if (errcheck) then
     !delta is too big, reduce it
     delta=delta/10
     if (delta<1.0e-15) goto 5
     errcheck=.false.
     goto 10
  end if

  df=(f2-f)/delta

END SUBROUTINE zeroeqn


!---------------------------------------------------------------!
!-------------------POSTSHOCK_PRESSURE--------------------------!
!---------------------------------------------------------------!


subroutine postshock_pressure(Vs,pb,ahead,behind,alfcheck)
  use type
  use global    !values of Bx and gamma
  implicit none
  real(DP),intent(IN)::Vs,pb
  real(DP),intent(IN),dimension(7)::ahead !value of primitives ahead the shock
  real(DP),intent(OUT),dimension(7)::behind
  logical,intent(IN),OPTIONAL::alfcheck

  !Given the values of the post-shock pressure and of the shock velocity,
  !returns the values of the postshock quantities using the p-method

  real(DP)::vxb,vyb,vzb,Db,Wb,v2b,alpha
  real(DP)::rhoa,pa,vxa,vya,vza,Bya,Bza,etaa,b2a,ha,v2a,Wa,Ws,j,Da,taua,wtota


  rhoa = ahead(1) !mass density
  pa   = ahead(2) !Total pressure (gas pressure + magnetic pressure)
  vxa  = ahead(3) !velocity
  vya  = ahead(4)
  vza  = ahead(5)
  Bya  = ahead(6) !magnetic field
  Bza  = ahead(7)

  etaa=Bx*vxa+Bya*vya+Bza*vza !B_j v^j

  v2a=vxa**2+vya**2+vza**2 !v^2
  if (v2a>=1.0e0_dp) stop 'error in postshock.f90: v^2a>c'
  Wa=1.0e0_dp/sqrt(1.0e0_dp-v2a) !Lorentz factor

  b2a=(Bx**2+Bya**2+Bza**2)/Wa**2+etaa**2 !b_mu b^mu

  !EoS
  call eos_enthalpy(pa-0.5e0_dp*b2a,rhoa,gamma,ha)
  !ha=1.0e0_dp+gamma*(pa-0.5e0_dp*b2a)/(rhoa*(gamma-1.0e0_dp)) !specific gas enthalpy

  wtota=rhoa*ha+b2a !total relativistic enthalpy

  Ws=1.0e0_dp/sqrt(1-Vs**2) !Lorentz factor of the shock

  j=Ws*rhoa*Wa*(Vs-vxa) !J

  Da=rhoa*Wa

  taua=(rhoa*ha+b2a)*Wa**2-pa-Da


  !-----------------------------------------------------!
  !--------------------VELOCITY-------------------------!
  !-----------------------------------------------------!


  vzb =  -((j ** 3*Wa **4 *(vza*(Da + pa + taua - etaa ** 2*Wa ** 2)*(-Bx ** 2 - &
       Bya ** 2 + Da + pb + taua - etaa ** 2*Wa ** 2) + &
       Bza*(-(etaa*(pb + taua)) + (pa + taua)*(Bx*vxa + Bya*vya) + &
       Da*(-etaa + Bx*vxa + Bya*vya) + etaa ** 3*Wa ** 2 - &
       etaa ** 2*(Bx*vxa + Bya*vya)*Wa ** 2)) - &
       Da*j**2*Wa**2*(-(Bx*Bza*(Da+pb+taua))+(-(Bx*Bza*(-3*etaa**2+2*Bya*etaa*vya + &
       Da*(-1 + vxa ** 2 + vya ** 2) + (pa + taua)*(-1 + vxa ** 2 + vya ** 2))) + &
       2 *Bx ** 3*etaa*vza + Bx*(2*Bya ** 2*etaa -etaa*(3*Da + 2*pa + pb + 3*taua) + &
       Bya*(Da + pa + taua)*vya)*vza +Bx ** 2*vxa*(-2*Bza*etaa +(Da+pa+taua)*vza) + &
       2 *(pa - pb)*vxa*(-(Bza*etaa) + (Da + pa + taua)*vza))* Wa ** 2 + &
       etaa ** 2 *(Bx*Bza*(-1 + vxa ** 2 + vya ** 2) - (-3*Bx*etaa + &
       Bx ** 2*vxa + 2*pa*vxa - 2*pb*vxa + Bx*Bya*vya)*vza)*Wa ** 4)*Ws - &
       Da ** 2*j* Wa ** 2*(-(Bx ** 4*vza) +2 *Bx*(pa - pb)*vxa*(Bza+etaa*vza*Wa**2) + &
       Bx ** 3*vxa*(Bza + 2*etaa*vza*Wa ** 2) +Bx ** 2*(vza*(-Bya**2+Da+pa+taua - &
       etaa*(3*etaa - 2*Bya*vya)*Wa ** 2) +Bza*(-3*etaa + Bya*vya - &
       2 *etaa*(-1 + vxa ** 2 + vya ** 2)*Wa ** 2)) + (pa - pb)*(-1 + vxa ** 2)* &
       Wa ** 2*((Da + pa + taua)*vza -etaa*(Bza + etaa*vza*Wa ** 2)))*Ws ** 2 + &
       Bx*Da **3 *(Bx ** 3*vxa*vza*Wa ** 2 - (pa - pb)*(-1 + vxa ** 2)* &
       Wa ** 2*(Bza + etaa*vza*Wa ** 2) + &
       Bx ** 2*(-Bza - (Bza*(-1 + vxa ** 2 + vya ** 2) + (etaa - &
       Bya*vya)*vza)*Wa ** 2))*Ws ** 3)/(Wa **2 *(-(Bx**4*Da**2*Ws**2*(j+Da*vxa*Ws)) + &
       Bx*Da*Wa ** 2*  Ws*(j **2 *(2*Bya ** 2*etaa + 2*Bza ** 2*etaa - &
       3 *etaa*(Da + pb + taua) + 3*etaa ** 3*Wa ** 2 + &
       Bya*vya*(Da - pa + 2*pb + taua - etaa ** 2*Wa ** 2) + &
       Bza*vza*(Da - pa + 2*pb + taua - etaa ** 2*Wa ** 2)) + &
       2 *Da*j*(pa - pb)*vxa*(2*etaa - Bya*vya - Bza*vza)*Ws + &
       Da ** 2*etaa*(pa - pb)*(-1 + vxa ** 2)*Ws ** 2) - &
       Wa ** 2*(j*(Da + pb + taua - etaa ** 2*Wa ** 2) + &
       Da*(-pa + pb)*vxa*Ws)*(-(j **2 *(Bya ** 2 + Bza ** 2 - Da - pb - taua + &
       etaa ** 2*Wa ** 2)) + 2 *Da*j*(-pa + pb)*vxa*Ws + &
       Da ** 2*(-pa + pb)*(-1 + vxa ** 2)*Ws ** 2) + &
       Bx ** 3*Da*Ws*(2*etaa*j ** 2*Wa ** 2 + 2*Da*etaa*j*vxa*Wa ** 2*Ws + &
       Da ** 2*(etaa - Bya*vya - Bza*vza)*Ws ** 2) + &
       Bx ** 2*(j ** 3*Wa ** 2*(Da + pb + taua - etaa ** 2*Wa ** 2) + &
       Da*j ** 2*vxa*  Wa ** 2*(Da - 2*pa + 3*pb + taua - etaa ** 2*Wa ** 2)*Ws + &
       Da ** 2*j*(-Bya ** 2 - Bza ** 2 + Da + pb + &
       taua + (-3*etaa ** 2 + 2 *etaa*(Bya*vya + Bza*vza) - (pa - pb)*(-1 + &
       3 *vxa ** 2 + vya ** 2 + vza ** 2))*Wa ** 2)*Ws ** 2))))




  vyb=  -((j ** 3*Wa ** 4 *(vya*(Da + pa + taua - etaa ** 2*Wa ** 2)*(-Bx ** 2 + Da + &
       pb + taua - etaa ** 2*Wa ** 2) + &
       Bya*(-(etaa*(pb + taua)) + Bx*(pa + taua)*vxa + & 
       Bza*(pb + taua)*vzb + Da*(-etaa + Bx*vxa + Bza*vzb) + & 
       etaa ** 3*Wa ** 2 - etaa ** 2*(Bx*vxa + Bza*vzb)*Wa ** 2)) - & 
       Da*j ** 2* Wa ** 2*(-(Bx*Bya*(Da + pb + taua)) + (2*Bx ** 3*etaa*vya + & 
       Bx ** 2*vxa*(-2*Bya*etaa + (Da + pa + taua)*vya) + &
       (pa - pb)*vxa*(2*(Da + pa + taua)*vya + Bya*(-2*etaa + Bza*vzb)) - &
       Bx*(etaa*(3*Da + 2*pa + pb + 3*taua)*vya - & 
       Bza*(pa - pb)*vya*vzb +Bya*(-3*etaa ** 2 + (pa + taua)*(-1 + vxa ** 2) + & 
       2 *Bza*etaa*vzb + (pb + taua)*vza*vzb + Da*(-1 + vxa ** 2 + vza*vzb))))*Wa ** 2 + &
       etaa ** 2 *(-(Bx ** 2*vxa*vya) + 2*(-pa + pb)*vxa*vya + & 
       Bx*(3*etaa*vya + Bya*(-1 + vxa ** 2 + vza*vzb)))*Wa ** 4)*Ws - & 
       Da ** 2*j*Wa ** 2*(-(Bx ** 4*vya) +Bx ** 3*vxa*(Bya + 2*etaa*vya*Wa ** 2) + &
       Bx*(pa - pb)*vxa*(2*Bya + (2*etaa*vya + Bza*vya*vzb + Bya*vza*vzb)*Wa ** 2) + &
       (pa - pb)*(-1 + vxa ** 2)*Wa ** 2*((Da + pa + taua)*vya - &
       etaa*(Bya + etaa*vya*Wa ** 2)) + Bx ** 2*(Bya*(-3*etaa + Bza*vzb - &
       2 *etaa*(-1 + vxa ** 2 + vza*vzb)*Wa ** 2) +&
       vya*(Da + pa + taua - (3*etaa ** 2 + (-pa + pb)*vza*vzb)* Wa ** 2)))*Ws ** 2 + &
       Bx*Da **3 *(Bx ** 3*vxa*vya*Wa ** 2 - & 
       Bx*(pa - pb)*vxa*vya*vza*vzb* Wa ** 4 - (pa - pb)*(-1 + vxa ** 2)* &
       Wa ** 2*(Bya + etaa*vya*Wa ** 2) + &
       Bx ** 2*(-Bya - (etaa*vya + Bya*(-1 + vxa ** 2 + vza*vzb))* &
       Wa ** 2))*Ws ** 3))!!!!!


  vyb=vyb/(Wa **2 *(-(Bx ** 4*Da ** 2*Ws ** 2*(j + Da*vxa*Ws)) + &
       Bx*Da*Wa ** 2* Ws*(j **2 *(2*Bya ** 2*etaa - 3*etaa*(Da + pb + taua) + &
       3 *etaa ** 3*Wa ** 2 + Bya*vya*(Da - pa + 2*pb + taua - etaa ** 2*Wa ** 2)) + &
       2 *Da*j*(pa - pb)*vxa*(2*etaa - Bya*vya)*Ws +&
       Da ** 2*etaa*(pa - pb)*(-1 + vxa ** 2)*Ws ** 2) - &
       Wa ** 2*(j*(Da + pb + taua - etaa ** 2*Wa ** 2) + &
       Da*(-pa + pb)*vxa*Ws)*(-(Bya ** 2*j ** 2) + &
       j ** 2*(Da + pb + taua - etaa ** 2*Wa ** 2) +2 *Da*j*(-pa + pb)*vxa*Ws + &
       Da ** 2*(-pa + pb)*(-1 + vxa ** 2)*Ws ** 2) + &
       Bx ** 3*Da* Ws*(2*etaa*j ** 2*Wa ** 2 + 2*Da*etaa*j*vxa*Wa ** 2*Ws + &
       Da ** 2*(etaa - Bya*vya)*Ws ** 2) + &
       Bx ** 2*(j ** 3*Wa ** 2*(Da + pb + taua - etaa ** 2*Wa ** 2) + &
       Da*j ** 2*vxa* Wa ** 2*(Da - 2*pa + 3*pb + taua - etaa ** 2*Wa ** 2)* Ws + &
       Da ** 2* j*(-Bya ** 2 + Da + pb + taua + (-3*etaa ** 2 +2 *Bya*etaa* &
       vya - (pa - pb)*(-1 + 3*vxa ** 2 + vya ** 2))*Wa ** 2)*Ws ** 2 + &
       Da**3*(-pa+pb)*vxa*(1+(-1+vxa**2+vya**2)*Wa ** 2)*Ws ** 3)))  


  vxb=(j*Wa**2*(-((Da+pa+taua)*vxa)+Bx*(etaa-Bya*vyb-Bza*vzb)+ & 
       etaa**2*vxa*Wa**2)- &
       Da*(Bx**2+(-pa+pb+Bx*etaa*vxa+Bx**2*(-1+vya*vyb+vza*vzb))* &
       Wa**2)*Ws)/(Wa**2*(-(j*(Da+pb+taua-etaa**2*Wa**2))- &
       Bx*Da*etaa*Ws+Da*(pa-pb)*vxa*Ws+Bx**2*(j+Da*vxa*Ws)))


  behind(3)=vxb !vxb
  behind(4)=vyb !vyb
  behind(5)=vzb !vzb


  !-----------------------------------------------------!
  !----------MAGNETIC FIELD COMPONENTS------------------!
  !-----------------------------------------------------!

  Db=1.0e0_dp/((vxa-vxb)*Ws/j+1.0e0_dp/Da)

  !Byb
  behind(6)=Db*(Bya/Da+Ws*Bx*vya/j-Ws*Bx*vyb/j)
  !Bzb
  behind(7)=Db*(Bza/Da+Ws*Bx*vza/j-Ws*Bx*vzb/j)


  !-----------------------------------------------------!
  !----------------MATTER DENSITY-----------------------!
  !-----------------------------------------------------!

  v2b=vxb**2+vyb**2+vzb**2 !v^2

  Wb=1.0e0_dp/sqrt(1.0e0_dp-v2b) !lorentz factor behind the shock

  !rhob
  behind(1)=Db/Wb

  !total pressure
  behind(2)=pb


  !Alfven discontinuity
  if(present(alfcheck).OR.(abs(pa-pb)<=epsilon(pb))) then
     if (present(alfcheck)) then
        if (.NOT.(alfcheck)) return
     end if

     stop 'no analytic solution for the Alfven dicontinuity at the moment!'

  end if

end subroutine postshock_pressure


!---------------------------------------------------------------!
!----------------------RAREFACTION_PRESSURE---------------------!
!---------------------------------------------------------------!

subroutine Rarefaction_pressure(pb,Vs,init,solution,switchLR)
  use type
  use interfaces,only:RHS_pressure,xi
  use wavecheck
  use global !This module contains verbose
  implicit none
  real(DP),intent(IN)::pb
  real(DP),intent(OUT)::Vs
  real(DP),dimension(7),intent(IN)::init
  real(DP),dimension(7),intent(OUT)::solution
  character(len=2),intent(IN)::switchLR

  !Solve the system of ODEs for rarefaction waves using the total pressure p as the self-similar variable

  real(DP),dimension(7)::y,k1,k2,k3,k4,f
  real(DP)::x,h
  integer(I4B)::j
  integer(I4B),PARAMETER::nsteps=1000 !number of steps used by the 4th order Runge-Kutta



  !Solve the ODE from init(2) to pb
  x=init(2)
  y=init

  !Head velocity
  call xi(y,switchLR,Vs)

  if (veryverbose) print *,'RAREFACTION: dp=',abs(pb-init(2))

  if (abs(pb-init(2))<=epsilon(pb)) then
     if (verbose) print *,'Pressure doesn''t change!'
     solution=init
     return
  end if

  !NOTE: pb<init(2)
  !(pb is the pressure at the Contact Discontinuity)

  h=(pb-init(2))/nsteps !use nsteps


  !Fourth order Runge-Kutta with a fixed step
  do j=1,nsteps
     call RHS_pressure(x,y,f,switchLR)
     if (wave_error) return 
     k1=h*f
     call RHS_pressure(x+0.5e0_dp*h,y+0.5e0_dp*k1,f,switchLR)
     if (wave_error) return 
     k2=h*f
     call RHS_pressure(x+0.5e0_dp*h,y+0.5e0_dp*k2,f,switchLR)
     if (wave_error) return 
     k3=h*f
     call RHS_pressure(x+h,y+k3,f,switchLR)
     k4=h*f

     y=y+(k1/6.0e0_dp)+(k2/3.0e0_dp)+(k3/3.0e0_dp)+(k4/6.0e0_dp)
     x=x+h

     !write(*,*) 'x=',x,'pb=',pb,'init(2)=',init(2), 'y(2)=',y(2) 

     if (wave_error) STOP 'wave error' !return 

  end do

  !write(*,*)
  !write(*,*) '----------------------------------------------'
  !write(*,*)

  solution=y


end subroutine Rarefaction_pressure


!---------------------------------------------------------------!
!----------------------RHS_PRESSURE-----------------------------!
!---------------------------------------------------------------!
subroutine RHS_pressure(x,y,f,switchLR)
  use type
  use interfaces,only:xi
  use global !This module contains Bx and gamma
  use wavecheck !This module contains wave_check
  implicit none
  real(DP),intent(IN)::x
  real(DP),dimension(7),intent(IN)::y
  real(DP),dimension(7),intent(OUT)::f
  character(len=2),intent(IN)::switchLR

  !subroutine called by Rarefaction
  !It returns the Right Hand Side of the ODE equations for rarefaction waves
  ! when using the p-method

  real(DP)::Vs,rho,P,vx,vy,vz,By,Bz,cs2,v2
  real(DP)::wtot,b2,Pgas,h,W,B2big,eta,lam,lap

  !NOTE
  !y(1) = rho
  !y(2) = total pressure P
  !y(3) = vx
  !y(4) = vy
  !y(5) = vz
  !y(6) = By
  !y(7) = Bz

  !x    = pb

  rho=y(1)
  P=x
  vx=y(3)
  vy=y(4)
  vz=y(5)
  By=y(6)
  Bz=y(7)


  v2=vx**2+vy**2+vz**2 !v^2

  if (v2>=1.0e0_dp) then
     print *
     print *,'RHS:v2>=1...'
     print *,'switchLR= ',switchLR
     wave_error=.true.
     return
     !stop 'riemann.f90: RHS: v2>1!'
  end if

  B2big=Bx**2+By**2+Bz**2 !B^2
  eta=Bx*vx+By*vy+Bz*vz !B^j v_j
  W=1.0e0_dp/sqrt(1.0e0_dp-v2) !Lorentz factor

  !compute Vs=\xi
  call xi(y,switchLR,Vs)


  !-----------------------EoS---------------------------!
  !EoS (sound velocity)
  b2=B2big/W**2+eta**2 !b_mu b^mu
  Pgas=P-0.5e0_dp*b2 !gas pressure
  call eos_enthalpy(Pgas,rho,gamma,h)!specific gas enthalpy
  wtot=rho*h + b2 !total relativistic enthalpy
  call eos_cs2(Pgas,rho,gamma,cs2)
  !cs2=gamma*Pgas/(wtot-b2) !c_s^2=\Gamma*P/w, squared sound speed
  !-----------------------------------------------------!


  !-----------------Alfven velocities-------------------!
  lap=vx+Bx*(1.0e0_dp-v2)/(eta+sqrt(wtot))!Right going
  lam=vx+Bx*(1.0e0_dp-v2)/(eta-sqrt(wtot))!Left going
  !-----------------------------------------------------!

  !dP/dP
  f(2)=1.0e0_dp

  !dv^x/dP
  f(3)=-(rho*h*W**2+Bx**2)*(vx-Vs)*(vx*Vs-1.0e0_dp)/(rho*h*(lap-Vs)*(lam-Vs))+&
       Bx**2*(Vs**2-1.0e0_dp)/(rho*h*W**2*(vx-Vs)*(lap-Vs)*(lam-Vs))+&
       Bx**2*Vs*(1-vx**2)/(rho*h*(lap-Vs)*(lam-Vs))+&
       Bx*(eta*(Vs**2-1.0e0_dp)-Bx*vx*(1.0e0_dp-2.0e0_dp*vx*Vs+Vs**2))/(rho*h*(lap-Vs)*(lam-Vs))
  f(3)=f(3)/(W**4*(eta**2-wtot))

  !dv^y/dP
  f(4)=-Bx**2*(vy*Vs*(Vs+vx))/(rho*h*(lap-Vs)*(lam-Vs))+&
       Bx*2.0e0_dp*vy*(eta-Bz*vz)*Vs/(rho*h*(lap-Vs)*(lam-Vs))+&
       (vy*(Bz**2+W**2*(eta**2-wtot))*(vx-Vs)*Vs+By**2*vy*(-1.0e0_dp+vx*Vs)+By*Bz*vz*(-1.0e0_dp+Vs**2))/&
       (rho*h*(lap-Vs)*(lam-Vs))+ &
       Bx*By*((-1.0e0_dp+vy**2+vz**2)+(vx-2.0e0_dp*vx*vy**2)*Vs+(1.0e0_dp+vy**2-vz**2)*Vs**2-vx*Vs**3)/&
       (rho*h*(vx-Vs)*(lap-Vs)*(lam-Vs))
  f(4)=f(4)/(W**4*(eta**2-wtot))

  !dv^z/dP
  f(5)=-Bx**2*(vz*Vs*(Vs+vx))/(rho*h*(lap-Vs)*(lam-Vs))+&
       Bx*2.0e0_dp*vz*(eta-By*vy)*Vs/(rho*h*(lap-Vs)*(lam-Vs))+&
       (vz*(By**2+W**2*(eta**2-wtot))*(vx-Vs)*Vs+Bz**2*vz*(-1.0e0_dp+vx*Vs)+By*Bz*vy*(-1.0e0_dp+Vs**2))/&
       (rho*h*(lap-Vs)*(lam-Vs))+ &
       Bx*Bz*((-1.0e0_dp+vy**2+vz**2)+(vx-2.0e0_dp*vx*vz**2)*Vs+(1.0e0_dp+vz**2-vy**2)*Vs**2-vx*Vs**3)/&
       (rho*h*(vx-Vs)*(lap-Vs)*(lam-Vs))
  f(5)=f(5)/(W**4*(eta**2-wtot))


  !drho/dP
  f(1)=-rho*(W**2*vx+1.0e0_dp/(vx-Vs))*f(3)-rho*W**2*vy*f(4)-rho*W**2*vz*f(5)

  !dBy/dP
  f(6)=-W**2*(By-By*vx*Vs+Bx*vy*Vs)/(Bx**2+2.0e0_dp*Bx*eta*W**2*(vx-Vs)+W**4*(eta**2-wtot)*(vx-Vs)**2)
  !dBz/dP
  f(7)=-W**2*(Bz-Bz*vx*Vs+Bx*vz*Vs)/(Bx**2+2.0e0_dp*Bx*eta*W**2*(vx-Vs)+W**4*(eta**2-wtot)*(vx-Vs)**2)


end subroutine RHS_pressure


subroutine alfven(Vs,ahead,behind)
  use type
  use global !contains Bx
  use interfaces,only:alfven_func,alfven_df,ludcmp,lubksb
  use wavecheck !This module contains wave_check
  implicit none
  real(DP),intent(IN)::Vs
  real(DP),dimension(7),intent(IN)::ahead
  real(DP),dimension(7),intent(INOUT)::behind

  !find the state behind an Alfven discontinuity using a Newton-Raphson method

  real(DP)::temp,accuracy,d,tolx,tolf,relax
  real(DP),dimension(5,5)::df
  real(DP),dimension(5)::f,p
  integer(I4B)::i,ntrial,nit
  integer(I4B),dimension(5)::indx

  accuracy=1.0e-13_dp !accuracy used in this subroutine

  relax=1.0e0_dp !relaxation factor

  ntrial=200 !Maximum number of steps
  tolx=epsilon(behind)
  tolf=accuracy
  nit=0


  call alfven_df(ahead,Vs,behind,f,df)

  do i=1,ntrial 
     !The subroutine alfven_df has supplied function values in f and Jacobian
     !matrix in df.

     nit=nit+1 !iterations counter

     if (veryverbose) print *,'Alfevn: I''m using the Newton-Raphson Method'
     if (veryverbose) then
        print *,'behind=',behind
        print *,'f=',f
     end if

     !Newton-Raphson routine
     if (sum(abs(f)) <= tolf) goto 1000 !Check function convergence. 
     p=-f !Right-hand side of linear equations. 


     !NR routine
     !Solve linear equations using LU decomposition.
     call ludcmp(df,indx,d) 
     call lubksb(df,indx,p) 

     p=p*relax !use a relaxation factor
     behind(3:7)=behind(3:7)+p !Update solution. 

     if (sum(abs(p)) <= tolx) goto 1000 !Check root convergence. 

     if (veryverbose) then
        print *
        print *,'++++ These are the values suggested by the algorithm ++++'
        print *,'behind=',behind
     end if

     call alfven_df(ahead,Vs,behind,f,df)

     if (veryverbose) then
        print *,'---------------------------'
        print *,'ITERATION END'
        print *,'desired accuracy=',tolf
        print *,'fvec =',f
        print *,'unknowns=',behind
        print *,'********************************************************'
        print *,'# of iterations=',niter
     end if

  end do

  wave_error=.true.
  return


1000 if (verbose) then
     print *,'*************ALFVEN WAVE***********************************' 
     print *,'accuracy',tolf,'reached!'
     print *,'number of iterations=',nit
     print *,'fvec =',f
     print *,'unknowns=',behind
  end if

end subroutine alfven



subroutine alfven_df(ahead,Vs,behind,f,df)
  use type
  use interfaces,only:alfven_func
  implicit none
  real(DP),intent(IN)::Vs
  real(DP),intent(IN),dimension(7)::behind,ahead
  real(DP),intent(OUT),dimension(5)::f
  real(DP),intent(OUT),dimension(5,5)::df

  !compute the Jacobian

  real(DP),dimension(5)::x,xsav,h,xph,f2
  real(DP),dimension(7)::behind_tmp
  real(DP)::EPS
  integer(I4B)::j
  
  EPS=1.0e-5_dp

  behind_tmp=behind

  call alfven_func(ahead,Vs,behind,f)
  
  x=behind(3:7)
  xsav=x 
  h=EPS*abs(xsav) 
  where (h == 0.0) h=EPS 
  xph=xsav+h !Trick to reduce finite precision error. 
  h=xph-xsav 
  do j=1,5
     x(j)=xph(j) 
     behind_tmp(3:7)=x
     call alfven_func(ahead,Vs,behind_tmp,f2)
     df(:,j)=(f2(:)-f(:))/h(j) !Forward difference formula 
     x(j)=xsav(j) 
  end do

end subroutine alfven_df


subroutine alfven_func(ahead,Vs,behind,f)
  use type
  use global
  use interfaces,only:eos_enthalpy
  implicit none
  real(DP),intent(IN)::Vs
  real(DP),dimension(7),intent(IN)::ahead,behind
  real(DP),intent(OUT),dimension(5)::f

  !compute equations (4.15)-(4.17) and (4.19)-(4.20)

  real(DP)::rhoa,pa,vxa,vya,vza,Bya,Bza,etaa,Wa,Da,b2a,wtota,taua
  real(DP)::rhob,pb,vxb,vyb,vzb,Byb,Bzb,etab,Wb,Db,b2b,wtotb,taub,tauOVD
  real(DP)::j,Ws,enthalpy


  !values ahead the Alfven discontinuity
  rhoa=ahead(1)
  pa=ahead(2)
  vxa=ahead(3)
  vya=ahead(4)
  vza=ahead(5)
  Bya=ahead(6)
  Bza=ahead(7)

  etaa=Bx*vxa+Bya*vya+Bza*vza !B^j v_j
  Wa=1.0e0_dp/sqrt(1.0e0_dp-vxa**2-vya**2-vza**2) !Lorentz factor
  Da=rhoa*Wa
  b2a=(Bx**2+Bya**2+Bza**2)/Wa**2+etaa**2 !b^mu b_mu
  call eos_enthalpy(pa-0.5e0_dp*b2a,rhoa,gamma,enthalpy)!specific gas enthalpy
  wtota=rhoa*enthalpy + b2a !total relativistic enthalpy
  taua=wtota*Wa**2-pa-Da

  Ws=1.0e0_dp/sqrt(1.0e0_dp-Vs**2) !Lorentz factor of the discontinuity
  j=rhoa*Wa*Ws*(Vs-vxa) !J
  
  !values behind the Alfven discontinuity
  rhob=rhoa
  pb=pa
  vxb = behind(3)
  vyb = behind(4)
  vzb = behind(5)
  Byb = behind(6)
  Bzb = behind(7)
  etab=Bx*vxb+Byb*vyb+Bzb*vzb !B^j v_j
  Wb=1.0e0_dp/sqrt(1.0e0_dp-vxb**2-vyb**2-vzb**2) !Lorentz factor
  Db=rhob*Wb
  b2b=b2a !b_mu b^mu (it's constant across Alfven discontinuities)
  wtotb=wtota !total relativistic enthalpy (it's constant across Alfven discontinuities)
  taub=wtotb*Wb**2-pb-Db
  tauOVD=taub/Db

  !eqnvx  !equation (4.15)
  f(1)= pa - pb - Bx*(etaa*vxa - etab*vxb) - Bx**2*(Wa**(-2) - Wb**(-2)) + (Bx*(etaa/Da - etab/Db)*j)/Ws - &
       (j*(vxa + (pa*vxa)/Da + (taua*vxa)/Da - vxb - (pb*vxb)/Db - tauOVD*vxb))/Ws + &
       (j*((etaa**2*vxa*Wa**2)/Da - (etab**2*vxb*Wb**2)/Db))/Ws

  !eqnvy  !equation (4.16)
  f(2)= -(Bx*(etaa*vya - etab*vyb)) - Bx*(Bya/Wa**2 - Byb/Wb**2) + (((Bya*etaa)/Da - (Byb*etab)/Db)*j)/Ws - &
       (j*(vya + (pa*vya)/Da + (taua*vya)/Da - vyb - (pb*vyb)/Db - tauOVD*vyb))/Ws + &
       (j*((etaa**2*vya*Wa**2)/Da - (etab**2*vyb*Wb**2)/Db))/Ws

  !eqnvz  !equation (4.17)
  f(3)= -(Bx*(etaa*vza - etab*vzb)) - Bx*(Bza/Wa**2 - Bzb/Wb**2) + (((Bza*etaa)/Da - (Bzb*etab)/Db)*j)/Ws - &
       (j*(vza + (pa*vza)/Da + (taua*vza)/Da - vzb - (pb*vzb)/Db - tauOVD*vzb))/Ws + &
       (j*((etaa**2*vza*Wa**2)/Da - (etab**2*vzb*Wb**2)/Db))/Ws

  !eqnBy  !equation (4.19)
  f(4)= Bx*(vya - vyb) + ((Bya/Da - Byb/Db)*j)/Ws

  !eqnBz  !equation (4.20)
  f(5)= Bx*(vya - vyb) + ((Bya/Da - Byb/Db)*j)/Ws

end subroutine alfven_func



!---------------------------------------------------------------!
!----------------------RAREFACTION_PSI--------------------------!
!---------------------------------------------------------------!

subroutine Rarefaction_psi(psib,Vs,init,solution,switchLR)
  use type
  use interfaces,only:RHS,xi
  use odeswitch
  use wavecheck
  use global !This module contains verbose
  implicit none
  real(DP),intent(IN)::psib !the value of psi behind the rarefaction wave
  real(DP),intent(OUT)::Vs
  real(DP),dimension(7),intent(IN)::init
  real(DP),dimension(7),intent(OUT)::solution
  character(len=2),intent(IN)::switchLR

  !Solve the system of ODEs for slow rarefaction waves 
  ! when the ratio By/Bz (i.e. the angle psi) varies across the wave
  !(used only within the Bt-method)

  real(DP),dimension(7)::y,k1,k2,k3,k4,f
  real(DP)::x,h
  integer(I4B)::j
  integer(I4B),PARAMETER::nsteps=1000 !number of steps used by the 4th order Runge-Kutta

  switch=switchLR !Contained in module odeswitch and used in odeint subroutines

  !Solve the ODE from psi(ahead) to psi(behind)
  x=atan(init(7)/init(6))

  if ((init(7)*init(6)<0.0e0_dp) .AND. (init(7)>0.0e0_dp)) then
     x=acos(init(6)/sqrt(init(6)**2+init(7)**2))
  end if

  y=init

  !Head velocity (Vs=xi)
  call xi(y,switchLR,Vs)

  if (veryverbose) print *,'dpsi=',abs(psib-x)

  if (abs(psib-x)<=epsilon(psib)) then
     if (veryverbose) print *,'psi doesn''t change!'
     solution=init
     return
  end if


  h=(psib-x)/nsteps !use nsteps


  !Fourth order Runge-Kutta with a fixed step
  do j=1,nsteps
     call RHS_psi(x,y,f,switchLR)
     if (wave_error) return 
     k1=h*f
     call RHS_psi(x+0.5e0_dp*h,y+0.5e0_dp*k1,f,switchLR)
     if (wave_error) return 
     k2=h*f
     call RHS_psi(x+0.5e0_dp*h,y+0.5e0_dp*k2,f,switchLR)
     if (wave_error) return 
     k3=h*f
     call RHS_psi(x+h,y+k3,f,switchLR)
     k4=h*f

     y=y+(k1/6.0e0_dp)+(k2/3.0e0_dp)+(k3/3.0e0_dp)+(k4/6.0e0_dp)
     x=x+h

     if (wave_error) return 

  end do

  solution=y


end subroutine Rarefaction_psi


!---------------------------------------------------------------!
!------------------------------RHS_PSI--------------------------!
!---------------------------------------------------------------!
subroutine RHS_psi(x,y,f,switchLR)
  use type
  use interfaces,only:xi
  use global !This module contains Bx and gamma
  use wavecheck !This module contains wave_check
  implicit none
  real(DP),intent(IN)::x
  real(DP),dimension(7),intent(IN)::y
  real(DP),dimension(7),intent(OUT)::f
  character(len=2),intent(IN)::switchLR
  !subroutine called by Rarefaction_psi
  !It returns the Right Hand Side of the ODE equations for rarefaction waves
  ! when the angle psi is used as independent variable

  real(DP)::psi,dV,Vs,rho,p,vx,vy,vz,By,Bz,normB,dBy,dBz,W,eta,B2,b2small,wtot,h,cs2

  real(DP)::denBz,psi1,psi2,drho,dvx,dvy,dvz,dptot

  !NOTE
  !y(1) = rho
  !y(2) = gas pressure P
  !y(3) = vx
  !y(4) = vy
  !y(5) = vz
  !y(6) = By
  !y(7) = Bz

  !x    = psi
  psi=x

  normB=sqrt(y(6)**2+y(7)**2) !Bt

  rho = y(1)
  P   = y(2)
  vx  = y(3)
  vy  = y(4)
  vz  = y(5)
  By  = y(6)
  Bz  = y(7)



  if ((vx**2+vy**2+vz**2)>1.0e0_dp) then
     wave_error=.true.
     return
  end if

  W = 1.0e0_dp/sqrt(1.0e0_dp-vx**2-vy**2-vz**2) !Lorentz factor

  eta = Bx*vx + By*vy + Bz*vz !B_j v^j

  B2 = Bx**2 + By**2 + Bz**2 !B^2

  b2small = B2/W**2 + eta**2 !b_mu b^mu

  !relativistic specific gas enthalpy
  call eos_enthalpy(p-0.5e0_dp*b2small,rho,gamma,h)
  wtot = rho*h + b2small !total relativistic enthalpy

  !sound velocity squared
  call eos_cs2(p-0.5e0_dp*b2small,rho,gamma,cs2)


  !compute xi (Vs=xi)
  call xi(y,switchLR,Vs)

  !dvz/dpsi
  dvz = (Bz*sin(psi)*(Bx**2*vz*W**2*(vx - Vs)*Vs*(vx + Vs) + &
       Bx*(-2*(eta - By*vy)*vz*W**2*(vx - Vs)*Vs + &
       Bz*(1 + vx**2*W**2 - vx*(1 + (vx**2 + vy**2 - vz**2)*W**2)*Vs + &
       (-1 + vy**2 - vz**2)*W**2*Vs**2 + vx*W**2*Vs**3)) - &
       W**2*(vx - Vs)*(By**2*vz*(vx - Vs)*Vs + By*Bz*vy*(-1 + Vs**2) + &
       vz*(W**2*(eta**2 - wtot)*(vx - Vs)*Vs + Bz**2*(-1 + vx*Vs)))))/ &
       ((sin(psi))**2*W**2*(2*Bx*(-eta + By*vy + Bz*vz) - &
       (By**2 + Bz**2 + W**2*(eta**2 - wtot))*(vx - Vs) + Bx**2*(vx + Vs))* &
       (Bx*(-(sin(psi)*vy) + cos(psi)*vz)*Vs))

  !dvx/dpsi
  dvx = (Bz*sin(psi)*(-(W**2*(By**2 + Bz**2 + W**2*(eta**2 - wtot))*(vx - Vs)**2* &
       (-1 + vx*Vs)) + Bx**2* &
       (1 - vx*(1 + (vy**2 + vz**2)*W**2)*Vs + (vy**2 + vz**2)*W**2*Vs**2) - &
       Bx*W**2*(vx - Vs)*(2*eta*(-1 + vx*Vs) - &
       (By*vy + Bz*vz)*(-1 + 2*vx*Vs - Vs**2))))/ &
       ((sin(psi))**2*W**2*(2*Bx*(-eta + By*vy + Bz*vz) - &
       (By**2 + Bz**2 + W**2*(eta**2 - wtot))*(vx - Vs) + Bx**2*(vx + Vs))* &
       (Bx*(-(sin(psi)*vy) + cos(psi)*vz)*Vs))

  !dp/dpsi
  dptot = (Bz*sin(psi)*(Bx**2 + 2*Bx*eta*W**2*(vx - Vs) + &
       W**4*(eta**2 - wtot)*(vx - Vs)**2))/ &
       ((sin(psi))**2*W**2*(Bx*(sin(psi)*vy - cos(psi)*vz)*Vs))

  !dvy/dpsi
  dvy = (By*dvx + (-(Bz/(sin(psi))**2) + (cos(psi)*(-(Bz*dvx) + Bx*dvz))/(sin(psi)*(vx - Vs)))* &
       (vx - Vs))/Bx

  !dBy/dpsi
  dBy = -(Bz/(sin(psi))**2) + (cos(psi)*(-(Bz*dvx) + Bx*dvz))/(sin(psi)*(vx - Vs))

  !dBz/dpsi
  dBz = (-(Bz*dvx) + Bx*dvz)/(vx - Vs)

  !drho/dpsi
  drho = -rho*(W**2*vx+1.0e0_dp/(vx-Vs))*dvx-rho*W**2*vy*dvy-rho*W**2*vz*dvz



  f(1)=drho
  f(2)=dptot
  f(3)=dvx
  f(4)=dvy
  f(5)=dvz
  f(6)=dBy
  f(7)=dBz


end subroutine RHS_psi
