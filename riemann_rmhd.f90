!$Id: riemann_rmhd.f90,v 1.8 2007-08-14 12:05:58 bgiacoma Exp $

program riemann_rmhd
  use type
  use global !This module contains Bx, gamma, degen, niter and initial_data
  use accmod !This module contains accuracy
  use interfaces,only:initialdata,contact,output,fullcontact,xi,eos_enthalpy
  use output_grid !This module contains x1, x2, t, nx
  use alfven_wave !contains an initial guess for the alfven waves
  use eos_param !Set the EOS
  implicit none
  integer(I4B)::verbo,init
  real(DP),dimension(7)::left,right
  real(DP)::pb,VsLeft,VsRight,dump
  real(DP),dimension(7,2)::solution
  real(DP),dimension(7,6)::fullsolution
  real(DP),dimension(3)::VsLv3,VsRv3
  !This program find the exact solution of the Riemann Problem
  !in Relativistic Magneto-HydroDynamics with an ideal Equation of State.
  !For a detailed description of the method see:
  !Giacomazzo & Rezzolla 2006, J. Fluid Mech. 562, 223-259
  !(arxiv.org: gr-qc/0507102)

  !Copyright (C) 2005  B. Giacomazzo, L. Rezzolla

  integer(I4B)::degen_case,eos_param_value
  real(DP)::leftalfven_left,rightalfven_left,leftalfven_right,rightalfven_right,eta,b2,wtot,enthalpy
  real(DP),dimension(4)::all_left,all_right
  real(DP),dimension(4)::unk
  real(DP),PARAMETER::small_vel=1.0d-15

  degen_case=0

  print *,'Select initial condition:'
  print *,' 0) User defined'
  print *,' 1) Marti-Muller figure 7 (hydro-test)'
  print *,' 2) Marti-Muller figure 6 (hydro-test)'
  print *,' 3) Marti-Muller figure 5 (hydro-test)'
  print *,' 4) Generic Shock-Tube Test (Bx=0)'
  print *,' 5) Komissarov: Shock-Tube Test 2 (Bx=0)'
  print *,' 6) Komissarov: Shock-Tube Test 1'
  print *,' 7) Komissarov: Collision'
  print *,' 8) Balsara Test 1'
  print *,' 9) Balsara Test 2'
  print *,'10) Balsara Test 3'
  print *,'11) Balsara Test 4'
  print *,'12) Balsara Test 5'
  print *,'13) Generic Alfven Test'



  read(*,*) initial_data

  print *,'EOS:'
  print *,'1) Ideal Fluid'
  print *,'2) Meliani et al.'
  read(*,*) eos_param_value

  if (eos_param_value==1) then
     eos_ideal = .true.    !we use an ideal EOS
     eos_meliani = .false. !we do not use the Meliani et al. EOS
  elseif (eos_param_value==2) then
     eos_ideal = .false.   !we do not use an ideal EOS
     eos_meliani = .true.  !we use the Meliani et al. EOS
  else
     write(*,*)'riemann_rmhd: wrong choice for the EOS'
     STOP
  end if


  print *,'left boundary (used only for output)='
  read(*,*) x1
  print *,'right boundary (used only for output)='
  read(*,*) x2
  print *,'time (used only for output)='
  read(*,*) t
  print *,'number of grid points to be used to plot the solution ='
  read(*,*) nx
  print *,'accuracy (usually 1.0e-10)='
  read(*,*) accuracy
  print *,'verbose (2=a lot of output, 1=normal, 0=only essential things)='
  read(*,*) verbo


  if (verbo==1) then
     veryverbose=.false.
     verbose=.true.
  else if (verbo==2) then
     veryverbose=.true.
     verbose=.true.
  else
     veryverbose=.false.
     verbose=.false.
  end if

  !initial conditions
  init=initial_data
  call initialdata(init,left,right)


  !Initial guess
  !NOTE:
  !unk(1)= p  in regions R2-R3 (total pressure)
  !unk(2)= By in regions R4-R5
  !unk(3)= Bz in regions R4-R5
  !unk(4)= p  in regions R6-R7 (total pressure)
  !unk is used only if Bx is different from zero and the hybrid method is used,
  ! otherwise there is no need for an initial guess from an approximate Riemann solver.
  select case(initial_data)
  case (0) !User defined
     unk=0.0e0_dp
     if (abs(Bx)>epsilon(Bx)) then
        print *
        print *,'!!!!!!!!!!!!!!!!!!WARNING!!!!!!!!!!!!!!!!!!!!!!!'
        print *,'If Bx is different from zero'
        print *,'you should provide a good initial guess'
        print *,'taken from an approximate Riemann solver.'
        print *,'If slow shocks are present then you probably need also'
        print *,'to have a look at subroutines velocity (lines 230-239)'
        print *,'and contact_velocity (lines 754 and 766)'
        print *,'contained in postshock.f90'
        print *
        print *,'press enter to continue...'
        read(*,*)
        print *,'approximate value of the total pressure in regions R2-R3 ='
        read(*,*) unk(1)
        print *,'approximate value of By in regions R4-R5 ='
        read(*,*) unk(2)
        print *,'approximate value of Bz in regions R4-R5 ='
        read(*,*) unk(3)
        print *,'approximate value of the total pressure in regions R6-R7 ='
        read(*,*) unk(4)
        print *
        print *,'Are alfven discontinuities present in the approximate solution? (y/n)'
        read(*,*) alfvenwave
        if (alfvenwave=='y') then
           print *,'give me an approximate value for the following'
           print *,'quantities behind the LEFT GOING Alfven discontinuity:'
           print *,'vx='
           read(*,*) leftalfven(1)
           print *,'vy='
           read(*,*) leftalfven(2)
           print *,'vz='
           read(*,*) leftalfven(3)
           print *,'By='
           read(*,*) leftalfven(4)
           print *,'Bz='
           read(*,*) leftalfven(5)
           print *
           print *,'give me an approximate value for the following'
           print *,'quantities behind the RIGHT GOING Alfven discontinuity:'
           print *,'vx='
           read(*,*) rightalfven(1)
           print *,'vy='
           read(*,*) rightalfven(2)
           print *,'vz='
           read(*,*) rightalfven(3)
           print *,'By='
           read(*,*) rightalfven(4)
           print *,'Bz='
           read(*,*) rightalfven(5)
        end if
     end if
  case (6)
     unk(1)=0.292650423E+02_dp   !For Komissarov: ShockTube1 
     unk(2)=0.0e0_dp             !For Komissarov: ShockTube1
     unk(3)=0.0e0_dp             !For Komissarov: ShockTube1
     unk(4)=0.292650423E+02_dp   !For Komissarov: ShockTube1
  case (7)
     unk(1)=257.0e0_dp   !For Komissarov: Collision 
     unk(2)=0.0e0_dp     !For Komissarov: Collision
     unk(3)=0.0e0_dp     !For Komissarov: Collision
     unk(4)=257.0e0_dp   !For Komissarov: Collision
  case (8)
     unk(1)=0.698933452e0_dp  !For BALSARA 1
     unk(2)=-0.428492941e0_dp !For BALSARA 1
     unk(3)=0.0e0_dp          !For BALSARA 1
     unk(4)=0.697622885e0_dp  !For BALSARA 1
  case (9)
     unk(1)=0.232143586E+02_dp  !For BALSARA 2
     unk(2)=0.320517414E+01_dp  !For BALSARA 2
     unk(3)=0.320517414E+01_dp  !For BALSARA 2
     unk(4)=0.207206998E+02_dp  !For BALSARA 2
  case (10)
     unk(1)=0.860438975E+02_dp  !For BALSARA 3
     unk(2)=0.466965496E+01_dp  !For BALSARA 3
     unk(3)=0.466965496E+01_dp  !For BALSARA 3
     unk(4)=0.636328552E+02_dp  !For BALSARA 3
  case (11)
     unk(1)=1183.5e0_dp  !For BALSARA 4
     unk(2)=0.0e0_dp     !For BALSARA 4
     unk(3)=0.0e0_dp     !For BALSARA 4
     unk(4)=1183.5e0_dp  !For BALSARA 4
  case (12)
     unk(1)=5.908e0_dp  !For BALSARA 5
     unk(2)=-1.175e0_dp !For BALSARA 5
     unk(3)=0.585e0_dp  !For BALSARA 5
     unk(4)=5.488e0_dp  !For BALSARA 5
  case (13)
     unk(1)=20.832e0_dp !For Generic Alfven Test
     unk(2)=5.1302e0_dp !For Generic Alfven Test
     unk(3)=0.7675e0_dp !For Generic Alfven Test
     unk(4)=20.853e0_dp !For Generic Alfven Test
  case default
     !cases that don't need an intial guess from HLLE (i.e. cases in which Bx=0)
     unk(1)=0.0e0_dp
     unk(2)=0.0e0_dp
     unk(3)=0.0e0_dp
     unk(4)=0.0e0_dp
  end select


  open(UNIT=1223,FILE='solution.sol',STATUS='REPLACE',ACTION='WRITE')
  open(UNIT=100,FILE='solution.dat',STATUS='REPLACE',ACTION='WRITE')
  open(UNIT=200,FILE='funcv.dat',STATUS='REPLACE',ACTION='WRITE')

  !looking for the presence of degeneracies
  degen=.false.
  if(.NOT.(Bx==0.0e0_dp)) then
     !Check that the alfven velocities are different from contact discontinuity's velocity
     !If they are equal the case is similar to Bx=0 case:
     !only 3 waves: two fast waves and the tangential discontinuity

     !Compute Alfven Velocities from the left state
     eta=Bx*left(3)+left(6)*left(4)+left(7)*left(5)
     b2=(Bx**2+left(6)**2+left(7)**2)*(1.0e0_dp-left(3)**2-left(4)**2-left(5)**2)+eta**2
     call eos_enthalpy(left(2)-0.5e0_dp*b2,left(1),gamma,enthalpy)
     wtot=left(1)*enthalpy + b2
     !right-going
     rightalfven_left=left(3)+Bx*(1.0e0_dp-left(3)**2-left(4)**2-left(5)**2)/(eta+sqrt(wtot))
     !left-going
     leftalfven_left=left(3)+Bx*(1.0e0_dp-left(3)**2-left(4)**2-left(5)**2)/(eta-sqrt(wtot))
     !check if the alfven velocities are equal to vx
     if ((abs(rightalfven_left-left(3))<epsilon(rightalfven_left)).OR.(abs(leftalfven_left-left(3))<epsilon(leftalfven_left))) &
          degen=.true.

     !Compute Alfven Velocities from the right state
     eta=Bx*right(3)+right(6)*right(4)+right(7)*right(5)
     b2=(Bx**2+right(6)**2+right(7)**2)*(1.0e0_dp-right(3)**2-right(4)**2-right(5)**2)+eta**2
     call eos_enthalpy(right(2)-0.5e0_dp*b2,right(1),gamma,enthalpy)
     wtot=right(1)*enthalpy + b2
     !right-going
     rightalfven_right=right(3)+Bx*(1.0e0_dp-right(3)**2-right(4)**2-right(5)**2)/(eta+sqrt(wtot))
     !left-going
     leftalfven_right=right(3)+Bx*(1.0e0_dp-right(3)**2-right(4)**2-right(5)**2)/(eta-sqrt(wtot))
     !check if the alfven velocities are equal to vx
     if ((abs(rightalfven_right-right(3))<epsilon(rightalfven_right)).OR.(abs(leftalfven_right-right(3))<epsilon(leftalfven_right))) &
          degen=.true.


     call xi(left ,'LF',dump,all_left) !compute slow and fast eigenvalues from the left state
     call xi(right,'LF',dump,all_right)!compute slow and fast eigenvalues from the right state

     !EXPERIMENTAL: this is true in cases such as ShockTube1 of Komissarov,
     !  but still need to be tested in other cases
     if ((abs(leftalfven_left-all_left(1))<=small_vel).AND.&
          (abs(rightalfven_left-all_left(4))<=small_vel)) then
        print *,'degeneracy ok kind 3, only SLOW WAVES'
        degen_case=3
        if (verbose) then
           print *,'left and right-going alfven velocity and vx at the left state'
           print *,leftalfven_left,rightalfven_left,left(3)
           print *,'left and right-going alfven velocity and vx at the right state'
           print *,leftalfven_right,rightalfven_right,right(3)
           print *,'fast and slow characteristic velocities from the left state'
           print *,all_left
           print *,'fast and slow characteristic velocities from the right state'
           print *,all_right
        end if
        goto 3
     else if ((abs(leftalfven_right-all_right(1))<=small_vel).AND.&
          (abs(rightalfven_right-all_right(4))<=small_vel)) then
        print *,'degeneracy ok kind 3, only SLOW WAVES'
        degen_case=3
        if (verbose) then
           print *,'left and right-going alfven velocity and vx at the left state'
           print *,leftalfven_left,rightalfven_left,left(3)
           print *,'left and right-going alfven velocity and vx at the right state'
           print *,leftalfven_right,rightalfven_right,right(3)
           print *,'fast and slow characteristic velocities from the left state'
           print *,all_left
           print *,'fast and slow characteristic velocities from the right state'
           print *,all_right
        end if
        goto 3
     else if ((abs(leftalfven_left-all_left(2))<=small_vel).AND.&
          (abs(rightalfven_left-all_left(3))<=small_vel)) then
        print *,'degeneracy ok kind 2, only FAST WAVES'
        degen_case=2
        degen=.true.
        goto 3
     else if ((abs(leftalfven_right-all_right(2))<=small_vel).AND.&
          (abs(rightalfven_right-all_right(3))<=small_vel)) then
        print *,'degeneracy ok kind 2, only FAST WAVES'
        degen_case=2
        degen=.true.
        goto 3
     end if

  end if

  !Find the solution at the contact discontinuity
3 if ((Bx==0.0e0_dp).OR.(degen)) then
     if (degen) then
        print *,'-----------------------'
        print *,'-----------------------'
        print *,'NO DIFFERENCE WITH Bx=0'
        if (verbose) then
           print *,'left and right-going alfven velocity and vx at the left state'
           print *,leftalfven_left,rightalfven_left,left(3)
           print *,'left and right-going alfven velocity and vx at the right state'
           print *,leftalfven_right,rightalfven_right,right(3)
           print *,'fast and slow characteristic velocities from the left state'
           print *,all_left
           print *,'fast and slow characteristic velocities from the right state'
           print *,all_right
        end if
        print *,'-----------------------'
        print *,'-----------------------'
     end if
     
     !Only two fast waves and a tangential discontinuity
     write(200,*) '#iteration, [[vx]] at the TD'
     write(200,*) !blank line
     call contact(left,right,pb,VsLeft,VsRight,solution,0) !p-method
     fullsolution=0.0e0_dp
     fullsolution(:,3)=solution(:,1)
     fullsolution(:,4)=solution(:,2)
     VsLv3=VsLeft
     VsRv3=VsRight
     unk=pb
  else if (degen_case==0) then
     !All the waves are present:
     !two fast and two slow waves,
     !two Alfven discontinuities and a contact discontinuity
     write(200,*) '#iteration, [[vx]],[[vy]],[[vz]],[[By]],[[Bz]],[[p]] at the CD'
     write(200,*) !blank line
     call fullcontact(left,right,unk,VsLv3,VsRv3,fullsolution) !hybrid method
  else if (degen_case==2) then !Only fast waves
     write(200,*) '#iteration, [[vx]] at the TD'
     write(200,*) !blank line
     call contact(left,right,pb,VsLeft,VsRight,solution,degen_case) !p-method
     fullsolution=0.0e0_dp
     fullsolution(:,3)=solution(:,1)
     fullsolution(:,4)=solution(:,2)
     VsLv3=VsLeft
     VsRv3=VsRight
  else if (degen_case==3) then !Only slow waves
     write(200,*) '#iteration, [[vx]] at the TD'
     write(200,*) !blank line
     call contact(left,right,pb,VsLeft,VsRight,solution,degen_case) !p-method
     fullsolution=0.0e0_dp
     fullsolution(:,3)=solution(:,1)
     fullsolution(:,4)=solution(:,2)
     VsLv3=VsLeft
     VsRv3=VsRight
  else
     stop 'riemann_rmhd.f90: riemann: error in degen_case'
  end if

  write(1223,*)
  write(1223,*)
  write(1223,*)
  write(1223,*) 'Exact solution found with accuracy',accuracy
  write(1223,*) 'Initial Condition     = ',init
  write(1223,*) 'left boundary (x1)    = ',x1
  write(1223,*) 'right boundary (x2)   = ',x2
  write(1223,*) 'time                  = ',t
  write(1223,*) 'number of points (nx) = ',nx
  write(1223,*)
  write(1223,*) 'Bx        = ',Bx
  write(1223,*) 'EOS Gamma = ',gamma
  write(1223,*)
  write(1223,*)
  write(1223,*) 'SOLUTION (at the LEFT side of the CD)'
  write(1223,*) 'rho = ',fullsolution(1,3)
  write(1223,*) 'p   = ',fullsolution(2,3)
  write(1223,*) 'vx  = ',fullsolution(3,3)
  write(1223,*) 'vy  = ',fullsolution(4,3)
  write(1223,*) 'vz  = ',fullsolution(5,3)
  write(1223,*) 'By  = ',fullsolution(6,3)
  write(1223,*) 'Bz  = ',fullsolution(7,3)
  write(1223,*)
  write(1223,*) 'SOLUTION (at the RIGHT side of the CD)'
  write(1223,*) 'rho = ',fullsolution(1,4)
  write(1223,*) 'p   = ',fullsolution(2,4)
  write(1223,*) 'vx  = ',fullsolution(3,4)
  write(1223,*) 'vy  = ',fullsolution(4,4)
  write(1223,*) 'vz  = ',fullsolution(5,4)
  write(1223,*) 'By  = ',fullsolution(6,4)
  write(1223,*) 'Bz  = ',fullsolution(7,4)
  write(1223,*)

  !write the exact values obtained in the different regions
  open(UNIT=4000,FILE='exact.sol',STATUS='REPLACE',ACTION='WRITE')
  write(4000,*) '#rho, total pressure p, vx, vy, vz, By, Bz'
  write(4000,*)
  if ((Bx==0.0e0_dp).OR.(degen)) then
     write(4000,'(7E20.9)') left
     write(4000,'(7E20.9)') fullsolution(:,3)
     write(4000,'(7E20.9)') fullsolution(:,4)
     write(4000,'(7E20.9)') right
  else
     write(4000,'(7E20.9)') left
     write(4000,'(7E20.9)') fullsolution(:,1)
     write(4000,'(7E20.9)') fullsolution(:,2)
     write(4000,'(7E20.9)') fullsolution(:,3)
     write(4000,'(7E20.9)') fullsolution(:,4)
     write(4000,'(7E20.9)') fullsolution(:,5)
     write(4000,'(7E20.9)') fullsolution(:,6)
     write(4000,'(7E20.9)') right
  end if
  close(4000)

  !Produce the output on a non-uniform grid.
  ! This is due to the fact that we want to draw the exact positions of shocks.
  !The number of points used can be slightly larger than the number of points nx chosen by the user
  call output(x1,x2,t,nx,left,right,VsLv3,VsRv3,unk,fullsolution)

  close(100)
  close(200)
  close(1223)

  print *
  print *,'- The eight states of the exact solution are written in file ''exact.sol'''
  print *
  print *,'- The solution is written in file ''solution.dat'' AT EACH ITERATION. This means that the exact solution computed with the desired accuracy is written at the end of the file.'
  print *
  print *,'- The jumps at the CD of vx,vy,vz,By,Bz and p at each iteration are saved in ''funcv.dat'''
  print *
  print *,'- The values at the contact discontinuity and the waves'' velocities are saved in ''solution.sol'' at each iteration'
  print *
  
end program riemann_rmhd



!---------------------------------------------------------------!
!-------------------------CONTACT-------------------------------!
!---------------------------------------------------------------!


subroutine contact(left,right,pb,VsL,VsR,solution,degen_case)
  use type
  use global
  use accmod !This module contains accuracy
  use interfaces,only:funcd,funcv
  implicit none
  real(DP),dimension(7),intent(IN)::left,right
  real(DP),intent(OUT)::pb
  real(DP),intent(OUT)::VsL,VsR
  real(DP),dimension(7,2),intent(OUT)::solution
  integer(I4B),intent(IN)::degen_case
  !This subroutine solves the riemann problem at the contact discontinuity using the p-method
  integer(I4B)::j,ILOOP
  real(DP)::temp,dpb,tolx,d,relax,pmin,pmax,fmin,fmax,check
  real(DP)::p,f,df,xl,xh,dpbold
  character(len=2)::left_wave,right_wave

  tolx=epsilon(pb)

  relax=1.0e0_dp !relaxation factor

  solution=0.0e0_dp


  !Are we considering a case with Bx different from zero but only one kind of waves, 
  ! i.e. only fast or slow waves?
  if ((degen_case==0).OR.(degen_case==1)) then
     left_wave='LF'
     right_wave='RF'
  else if (degen_case==2) then
     left_wave='LF'
     right_wave='RF'
  else if (degen_case==3) then
     left_wave='LS'
     right_wave='RS'
  end if



  !Braketing the function
  pmin=0.5e0_dp*(left(2) + right(2))
  pmax=pmin

  print *
  print *,'Bracketing the solution...'
  print *

  !----------ONLY FOR DEBUG--------------------------!
  if(.false.) then
     open(UNIT=12001,FILE='funcv_debug.dat',STATUS='REPLACE',ACTION='WRITE')
     pmin=1.0e-15
     pmax=1.0e-10
     dpbold=(pmax-pmin)/1001
     p=pmin
     do j=1,1001
        call funcv(left,right,p,VsL,VsR,fmin,left_wave,right_wave)
        write(12001,*)p,fmin
        p=p+dpbold
     end do
     CLOSE(12001)
     STOP 'data plotted in funcv_debug.dat'
  end if
  !----------ONLY FOR DEBUG--------------------------!

  do
     ILOOP = ILOOP + 1

     call funcv(left,right,pmin,VsL,VsR,fmin,left_wave,right_wave)
     call funcv(left,right,pmax,VsL,VsR,fmax,left_wave,right_wave)
     
     if ((abs(fmin)<accuracy).and.(abs(fmin)>0)) then !check abs(fmin)>0 to avoid NaN
        !We have the solution
        print *,'pmin=',pmin,'fmin=',fmin
        pb=pmin
        goto 1000
     end if

     if ((abs(fmax)<accuracy).and.(abs(fmax)>0)) then !check abs(fmax)>0 to avoid NaN
        !We have the solution
        print *,'pmax=',pmax,'fmax=',fmax
        pb=pmax
        goto 1000
     end if

     CHECK = fmin*fmax 
     print *,'p1=',pmin,'p2=',pmax
     print *,'f(p1)=',fmin,'f(p2)=',fmax
     print *,'f(p1)*f(p2)=',CHECK
     IF (CHECK<0.0e0_dp) GOTO 5 !the root is in the interval (pmin,pmax)
     
     !The function does not change sign in this interval, move the boundaries
     if (pmin>1e-15_dp) pmin = 0.9e0_dp*pmin
     pmax = 1.1e0_dp*pmax
  end do

5  if (fmin < 0.0) then !Orient the search so that f(xl)<0. 
     xl=pmin 
     xh=pmax
  else 
     xh=pmin 
     xl=pmax
  end if

  print * !blank line

  !Find the pressure using the continuity of v^x at the contact discontinuity
  !Use Newton-Raphson and bisection method
  pb=0.5e0_dp*(pmin+pmax)
  dpbold=abs(pmax-pmin)
  dpb=dpbold
  call funcd(left,right,pb,VsL,VsR,f,df,left_wave,right_wave)
  if (f==0.0e0_dp) goto 1000
  do j=1,200 !the maximum number of steps is 200
     if (((pb-xh)*df-f)*((pb-xl)*df-f) > 0.0 .or. & 
          abs(2.0_DP*f) > abs(dpbold*df) ) then 
        !Bisect if Newton out of range, or not decreasing fast enough. 
        dpbold=dpb 
        dpb=0.5e0_dp*(xh-xl) 
        pb=xl+dpb 
        if (xl == pb) goto 1000 !Change in root is negligible. 
     else !Newton step acceptable. Take it. 
        dpbold=dpb
        dpb=f/df
        temp=pb
        pb=pb-dpb
        if (temp==pb) goto 1000
     end if
     
     if (abs(dpb)<tolx) goto 1000 !Convergence
     call funcd(left,right,pb,VsL,VsR,f,df,left_wave,right_wave)!One new function evaluation per iteration. 
     print *
     print *,'pb=', pb,'f=',f,'desired accuracy=',accuracy
     print *
     niter=j !number of iterations
     write(200,*) j,f
     
     if (abs(f)<accuracy) goto 1000 !Convergence
     if (f < 0.0) then !Maintain the bracket on the root. 
        xl=pb 
     else 
        xh=pb
     end if

     

  end do
  print *, 'I cannot find the post-shock pressure'
  print *, 'Try to reduce the accuracy'
  stop

1000 call funcv(left,right,pb,VsL,VsR,f,left_wave,right_wave,solution)

  print *
  print *,'left  side of TD',solution(:,1)
  print *,'right side of TD',solution(:,2)
  print *

end subroutine contact



!---------------------------------------------------------------!
!-----------------------FULLCONTACT-----------------------------!
!---------------------------------------------------------------!

subroutine fullcontact(left,right,unk,VsLv3,VsRv3,fullsolution)
  use type
  use global !This module contains Bx, the EOS gamma, degen, niter and solmethod
  use accmod!This module contains accuracy
  use interfaces, only:fullfuncd,fullfuncv,ludcmp,lubksb
  use wavecheck
  use odeswitch
  use output_grid !This module contains t,x1,x2,nx
  implicit none

  real(DP),dimension(7),intent(IN)::left,right
  real(DP),dimension(4),intent(INOUT)::unk
  real(DP),dimension(3),intent(OUT)::VsLv3,VsRv3
  real(DP),dimension(7,6),intent(OUT)::fullsolution
  !Given the left and right states, this subroutine finds the solution at
  !the contact discontinuity in the case Bx different from 0.
  !It uses the pressure as an unknown between the fast and slow waves, and
  ! the tangential components of the magnetic field between the slow waves (i.e.
  ! at the contact discontinuity)

  !unk(1)= p  in regions R2-R3 (total pressure)
  !unk(2)= By in regions R4-R5
  !unk(3)= Bz in regions R4-R5
  !unk(4)= p  in regions R6-R7 (total pressure)

  integer(I4B) :: i,ntrial
  real(DP) :: d,tolx,tolf
  real(DP), dimension(size(unk)) :: fvec,p
  integer(I4B),dimension(size(unk)) ::indx
  real(DP), dimension(size(unk),size(unk)) :: fjac
  real(DP), parameter :: EPS=epsilon(x1)
  real(DP)::relax

  ntrial=200 !maximum number of steps
  tolx=epsilon(unk)
  tolf=accuracy
  niter=0

  !-------------------------------------------------------------------!
  !Check that the function is defined in the first point given by unk !
  !-------------------------------------------------------------------!

  shock=.false.
  wave_error=.false.

  call fullfuncv(left,right,unk,VsLv3,VsRv3,fvec,fullsolution)
  !Test if the starting point is a valid one

  if (wave_error) stop 'riemann.f90: fullcontact: chose a different initial guess. This isn''t valid!'

  !Plot the results
  call output(x1,x2,t,nx,left,right,VsLv3,VsRv3,unk,fullsolution)

  if (maxval(abs(fvec)) <= tolf) goto 1000 !Check function convergence. 

  relax=1.0 !RELAXATION FACTOR


  !compute the Jacobian
  call fullfuncd(left,right,unk,VsLv3,VsRv3,fvec,fjac)


  do i=1,ntrial 
     !The subroutine fullfuncd has supplied function values at unk in fvec and Jacobian
     ! matrix in fjac.

     niter=niter+1 !iterations counter

     if (verbose) then
        print *,'unk=',unk
        print *,'fvec=',fvec
        print *,'fjac(1,:)=',fjac(1,:)
        print *,'fjac(2,:)=',fjac(2,:)
        print *,'fjac(3,:)=',fjac(3,:)
        print *,'fjac(4,:)=',fjac(4,:)
        print *
     end if

     !Newton-Raphson routine
     if (maxval(abs(fvec)) <= tolf) goto 1000 !Check function convergence. 
     p=-fvec !Right-hand side of linear equations. 

     !NR routine
     !Solve linear equations using LU decomposition.
     call ludcmp(fjac,indx,d) 
     call lubksb(fjac,indx,p) 

     p=p*relax !use relaxation factor
     unk=unk+p !Update solution. 

     if (sum(abs(p)) <= tolx) goto 1000 !Check root convergence. 

     if (verbose) then
        print *
        print *,'++++ These are the values suggested by the algorithm ++++'
        print *,'[ptot(2), By(4), Bz(4), ptot(7)]=',unk
     end if

     !---------------------------------------------------------------------!
     !Tring to avoid negative values for the pressure
     if(unk(1)<=1.0e-4_dp) unk(1)=1.0e-4_dp
     if(unk(4)<=1.0e-4_dp) unk(4)=1.0e-4_dp
     !---------------------------------------------------------------------!

     if (verbose) then
        print *, 'and these are the values I will use:'
        print *,'[ptot(2), By(4), Bz(4), ptot(7)]=',unk
        print *
     end if


     shock=.false.
     wave_error=.false.

     print *,'---------------------------'
     print *,'ITERATION END'
     print *,'desired accuracy=',tolf
     print *,'fvec =',fvec
     print *,'new guess=',unk
     print *,'********************************************************'
     print *,'# of iterations=',niter

     !Compute fvec and its Jacobian with the new guess
     call fullfuncd(left,right,unk,VsLv3,VsRv3,fvec,fjac,fullsolution)

     !----------------------------------------------!
     !At the end of each iteration, plot the results!
     !----------------------------------------------!

     call output(x1,x2,t,nx,left,right,VsLv3,VsRv3,unk,fullsolution)

  end do

  print *, 'riemann.f90: fullcontact: unable to find the solution'
  print *, 'try to reduce the accuracy'
  stop

1000 print *,'***********************************************************' 
  print *,'accuracy',tolf,'reached!'
  print *,'number of iterations=',niter
  print *,'fvec =',fvec
  print *,'unknowns=',unk
  print *
  print *,'Computing the fullsolution...'
  print *

  call fullfuncv(left,right,unk,VsLv3,VsRv3,fvec,fullsolution)

end subroutine fullcontact


!---------------------------------------------------------------!
!---------------------------FUNCD-------------------------------!
!---------------------------------------------------------------!

subroutine funcd(left,right,pb,VsL,VsR,f,df,left_wave,right_wave)
  use type
  use wavecheck
  use interfaces,only: funcv
  implicit none
  real(DP),dimension(7),intent(IN)::left,right
  real(DP),intent(IN)::pb
  real(DP),intent(OUT)::VsL,VsR
  real(DP),intent(OUT)::f,df
  character(len=2),intent(IN)::left_wave,right_wave

  !This subroutine computes the first derivative of f

  real(DP)::f2,dpb,dumpVsL,dumpVsR

  call funcv(left,right,pb,VsL,VsR,f,left_wave,right_wave)

  dpb=1.0D-5*abs(pb)
  if (dpb==0.0e0_dp) dpb=1.0D-5

  call funcv(left,right,pb+dpb,dumpVsL,dumpVsR,f2,left_wave,right_wave)
  
  !Compute the first derivative
  df=(f2-f)/dpb

end subroutine funcd


!---------------------------------------------------------------!
!---------------------------FUNCV-------------------------------!
!---------------------------------------------------------------!

subroutine funcv(left,right,pb,VsL,VsR,f,left_wave,right_wave,solution)
  use type
  use global
  use wavecheck
  use interfaces,only:Rarefaction_pressure,ContactVelocity
  implicit none
  real(DP),dimension(7),intent(IN)::left,right
  real(DP),intent(IN)::pb
  real(DP),intent(OUT)::VsL,VsR
  real(DP),intent(OUT)::f
  character(len=2),intent(IN)::left_wave,right_wave
  real(DP),dimension(7,2),intent(OUT),OPTIONAL::solution

  !This subroutine computes the jump in vx at the tangential discontinuity

  real(DP),dimension(7)::solutionL,solutionR
  real(DP)::vx1,vx2,dump
  

  if (pb>left(2)) then
     !Shock
     call ContactVelocity(pb,0.0e0_dp,0.0e0_dp,VsL,left,dump,left_wave,'P',solutionL)
     if (wave_error) stop 'funcv: error solving the LF shock equations'
  elseif (pb==left(2)) then
     solutionL=left
  else
     !Rarefaction
     call Rarefaction_pressure(pb,VsL,left,solutionL,left_wave)
     if (wave_error) stop 'funcv: error solving the LF rarefaction equations'
  end if
  vx1=solutionL(3)!vx at the left of the TD


  if (pb>right(2)) then
     !Shock
     call ContactVelocity(pb,0.0e0_dp,0.0e0_dp,VsR,right,dump,right_wave,'P',solutionR)
     if (wave_error) stop 'funcv: error solving the RF shock equation'
  elseif (pb==right(2)) then
     solutionR=right
  else
     !Rarefaction
     call Rarefaction_pressure(pb,VsR,right,solutionR,right_wave)
     if (wave_error) stop 'funcv: error solving the LF rarefaction equations'
  end if
  vx2=solutionR(3)!vx at the right of the TD

  f=vx1-vx2 !This is zero if the solution is exact

  if (present(solution)) then
     solution(:,1)=solutionL
     solution(:,2)=solutionR
     !output the accuracy f reached at iteration niter to funcv.dat
     write(200,'(I3,2E20.8)') niter,f 
  end if

end subroutine funcv

!---------------------------------------------------------------!
!-------------------------FULLFUNCD-----------------------------!
!---------------------------------------------------------------!

subroutine fullfuncd(left,right,unk,VsLv3,VsRv3,fvec,fjac,fullsolution)
  use type
  use interfaces,only:fullfuncv
  use wavecheck
  use global !It contains verbo
  implicit none
  real(DP),dimension(7),intent(IN)::left,right
  real(DP),dimension(4),intent(IN)::unk
  real(DP),dimension(3),intent(OUT)::VsLv3,VsRv3
  real(DP),dimension(4),intent(OUT)::fvec
  real(DP),dimension(4,4),intent(OUT)::fjac
  real(DP),dimension(7,6),intent(OUT),OPTIONAL::fullsolution

  !Compute the Jacobian of fvec

  real(DP):: EPS
  real(DP),dimension(4)::x,xsav,xph,h,fvec2,dumpv3,dumpv3b
  integer(I4B)::j

  EPS=1.0D-5

  if(present(fullsolution)) then
     call fullfuncv(left,right,unk,VsLv3,VsRv3,fvec,fullsolution)
  else
     call fullfuncv(left,right,unk,VsLv3,VsRv3,fvec)
  end if
  
  
  if (wave_error) return !The function is not defined in these points


  !Now computes the Jacobian fjac
  x=unk
  xsav=x 
  h=EPS*abs(xsav) 
  where (h == 0.0) h=EPS 
  xph=xsav+h !Trick to reduce finite precision error. 
  h=xph-xsav 
  do j=1,4
     x(j)=xph(j) 
     call fullfuncv(left,right,x,dumpv3,dumpv3b,fvec2)
     if (wave_error) stop 'fullfuncd: error computing the Jacobian matrix'
     fjac(:,j)=(fvec2(:)-fvec(:))/h(j) !Forward difference formula 
     x(j)=xsav(j) 
  end do


end subroutine fullfuncd


!---------------------------------------------------------------!
!-------------------------FULLFUNCV-----------------------------!
!---------------------------------------------------------------!

subroutine fullfuncv(left,right,unk,VsLv3,VsRv3,fvec,fullsolution)
  use type
  use global !This module contains Bx and gamma (the EOS' gamma)
  use interfaces,only: velocity,Rarefaction,alfven,eos_enthalpy
  use interfaces,only: Rarefaction_pressure,ContactVelocity
  use wavecheck !This module contains wave_check
  use accmod !It contains accuracy
  use odeswitch !It contains switch
  use alfven_wave !guess for the alfven discontinuities (if present)
  implicit none
  real(DP),dimension(7),intent(IN)::left,right
  real(DP),dimension(4),intent(IN)::unk
  real(DP),dimension(3),intent(OUT)::VsLv3,VsRv3
  real(DP),dimension(4),intent(OUT)::fvec
  real(DP),dimension(7,6),intent(OUT),OPTIONAL::fullsolution

  !compute vx, vy, vz and p at the left and right boundaries of the
  !Contact Discontinuity; then takes the differences which should be zero

  real(DP),dimension(7)::zone1a,zone1b,zone2a,zone2b,zone3a,zone3b,zonedump
  real(DP)::vx1,vy1,vz1,ptot1,vx2,vy2,vz2,ptot2,vxb,B2big,eta,W2inv,b2,wtot,normB,psi,dump,enthalpy

  real(DP),parameter::small_psi=1.0e-4_dp

  wave_error=.false.

  !LEFT GOING FAST WAVE
  if(verbose) then
     print *
     print *,'--------------------------------------'
     print *,'I''m solving equations for the'
     print *,'Left Going Fast Wave'
     print *,'--------------------------------------'
     print *
  end if

  !Is it a shock?
  !If it's a shock pressure increases
  if (abs(unk(1)-left(2))<epsilon(unk(1))) then
     if (verbose) print *,'no shock or rarefaction!'
     zone1a=left
     goto 10
  end if

  if (unk(1)>=left(2)) then
     !SHOCK
     if (verbose) print *,'It''s a shock!'
     !compute the velocity of the left going fast shock
     call ContactVelocity(unk(1),0.0e0_dp,0.0e0_dp,VsLv3(1),left,dump,'LF','P',zone1a)
     if (wave_error) stop 'fullfuncv: error in the left-going fast shock'
  else
     !RAREFACTION
     if (verbose) print *,'It''s a rarefaction!'
     call Rarefaction_pressure(unk(1),VsLv3(1),left,zone1a,'LF')
     if (wave_error) stop 'fullfuncv: error in the left-going fast rarefaction'
  end if




  !LEFT GOING ALFVEN DISCONTINUITY
10 if (verbose) then
     print *
     print *,'--------------------------------------'
     print *,'I''m solving equations for the'
     print *,'Left Going Alfven Wave'
     print *,'--------------------------------------'
     print *
  end if

  !Alfven velocity
  B2big=Bx**2+zone1a(6)**2+zone1a(7)**2 !Bx**2+By**2+Bz**2
  eta=Bx*zone1a(3)+zone1a(6)*zone1a(4)+zone1a(7)*zone1a(5) !Bx*vx+By*vy+Bz*vz
  W2inv=1.0e0_dp-zone1a(3)**2-zone1a(4)**2-zone1a(5)**2 !Inverse of squared Lorentz factor
  b2=B2big*W2inv+eta**2 !b_mu b^mu
  
  call eos_enthalpy(zone1a(2)-0.5e0_dp*b2,zone1a(1),gamma,enthalpy)
  wtot=zone1a(1)*enthalpy + b2 !Total relativistic enthalpy
  VsLv3(2)= zone1a(3)+Bx*W2inv/(eta-sqrt(wtot))!Alfven velocity = vx-Bx/(W**2*(sqrt(wtot)-eta))


  !First guess for the Alfven subroutine
  zone1b=zone1a


  if (initial_data==13) then
     zone1b(3)= 0.0710516 !for the generic Alfven test (LAW)
     zone1b(4)= 0.36694   !for the generic Alfven test (LAW)
     zone1b(5)= 0.242851  !for the generic Alfven test (LAW)
     zone1b(6)= 5.691     !for the generic Alfven test (LAW)
     zone1b(7)= 0.850178  !for the generic Alfven test (LAW)
  elseif (initial_data==12) then
     zone1b(3)=-0.121516834e0_dp !for Balsara5 (LAW)
     zone1b(4)=0.126358623e0_dp  !for Balsara5 (LAW)
     zone1b(5)=0.115811357e0_dp  !for Balsara5 (LAW)
     zone1b(6)=-0.118226506e0_dp !for Balsara5 (LAW)
     zone1b(7)=0.230194489e0_dp  !for Balsara5 (LAW)
  elseif ((initial_data==0).AND.(alfvenwave=='y')) then
     zone1b(3)= leftalfven(1) !for user defined Riemann problem (LAW)
     zone1b(4)= leftalfven(2) !for user defined Riemann problem (LAW)
     zone1b(5)= leftalfven(3) !for user defined Riemann problem (LAW)
     zone1b(6)= leftalfven(4) !for user defined Riemann problem (LAW)
     zone1b(7)= leftalfven(5) !for user defined Riemann problem (LAW)
  end if

  !solve the system of equations (4.15)-(4.17),(4.19)-(4.20)
  call alfven(VsLv3(2),zone1a,zone1b)

  if (wave_error) stop 'fullfuncv: error in the left-going alfven discontinuity'

  !LEFT GOING SLOW WAVE
  if (verbose) then
     print *
     print *,'--------------------------------------'
     print *,'I''m solving equations for the'
     print *,'Left Going Slow Wave'
     print *,'--------------------------------------'
     print *
  end if


  ! Is it a shock?
  ! If it's a slow shock normB decreases
  normB=sqrt(unk(2)**2+unk(3)**2)

  if (abs(normB-sqrt(zone1b(6)**2+zone1b(7)**2))<epsilon(unk(2))) then
     if (verbose) print *,'no shock or rarefaction!'
     zone2a=zone1b
     goto 20
  end if

  if (normB<=sqrt(zone1b(6)**2+zone1b(7)**2)) then
     !SHOCK
     if (verbose) print *,'It''s a shock!'
     !compute the velocity of the left going slow shock
     call velocity(unk(2),unk(3),zone1b,'LS',vxb,VsLv3(3),'B',zone2a)
     if (wave_error) stop 'error in the left-going slow shock' 
  else
     !RAREFACTION
     if (verbose) print *,'It''s a rarefaction!'

     if (abs((zone1b(7)/zone1b(6))-(unk(3)/unk(2)))<=small_psi) then
        !psi is constant
        if (verbose) print *,'psi is constant'
        call Rarefaction(normB,VsLv3(3),zone1b,zone2a,'LS')
     else
        !psi is not constant
        if (verbose) print *,'psi is NOT constant'
        psi=atan(unk(3)/unk(2))
        if (((unk(3)*unk(2))<0.0e0_dp) .AND. (unk(3)>0.0e0_dp)) then
           psi=acos(unk(2)/sqrt(unk(2)**2+unk(3)**2))
        end if
        call Rarefaction_psi(psi,VsLv3(3),zone1b,zone2a,'LS')
     end if
     if (wave_error) stop 'error in the left-going slow rarefaction' 
  end if



  !Values of velocity and pressure at the left side of the contact discontinuity
20 vx1=zone2a(3)
  vy1=zone2a(4)
  vz1=zone2a(5)
  ptot1=zone2a(2)

  !RIGHT GOING FAST WAVE
  if (verbose) then
     print *
     print *,'--------------------------------------'
     print *,'I''m solving equations for the'
     print *,'Right Going Fast Wave'
     print *,'--------------------------------------'
     print *
  end if


  ! Is it a shock?
  ! If it's a fast shock pressure increases
  if (abs(unk(4)-right(2))<epsilon(unk(4))) then
     if (verbose) print *,'no shock or rarefaction!'
     zone3b=right
     goto 11
  end if
  
  if (unk(4)>=right(2)) then
     !SHOCK
     if (verbose) print *,'It''s a shock!'
     !compute the velocity of the right going fast shock
     call ContactVelocity(unk(4),0.0e0_dp,0.0e0_dp,VsRv3(1),right,dump,'RF','P',zone3b)
     if (wave_error) stop 'fullfuncv: error at the right-going fast shock'
  else
     !RAREFACTION
     if (verbose) print *,'It''s a rarefaction!'
     call Rarefaction_pressure(unk(4),VsRv3(1),right,zone3b,'RF')
     if (wave_error) stop 'fullfuncv: error at the right-going fast rarefaction'
  end if



  !RIGHT GOING ALFVEN DISCONTINUITY
11 if (verbose) then
     print *
     print *,'--------------------------------------'
     print *,'I''m solving equations for the'
     print *,'Right Going Alfven Wave'
     print *,'--------------------------------------'
     print *
  end if


  !Alfven velocity (in lagrangian mass coordinates)
  B2big=Bx**2+zone3b(6)**2+zone3b(7)**2 !Bx**2+By**2+Bz**2
  eta=Bx*zone3b(3)+zone3b(6)*zone3b(4)+zone3b(7)*zone3b(5) !Bx*vx+By*vy+Bz*vz
  W2inv=1.0e0_dp-zone3b(3)**2-zone3b(4)**2-zone3b(5)**2 !Inverse of squared Lorentz factor
  b2=B2big*W2inv+eta**2 !b_{\mu} b^{\mu}

  call eos_enthalpy(zone3b(2)-0.5e0_dp*b2,zone3b(1),gamma,enthalpy)
  wtot=zone3b(1)*enthalpy + b2 !Total relativistic enthalpy
  VsRv3(2)= zone3b(3)+Bx*W2inv/(eta+sqrt(wtot))!Alfven velocity = vx+Bx/(W**2*(sqrt(wtot)+eta))


  !First guess for the Alfven subroutine
  zone3a=zone3b


  if (initial_data==13) then  
     zone3a(3)= 0.016e0_dp  !for the generic Alfven test (RAW)
     zone3a(4)= -0.048e0_dp !for the generic Alfven test (RAW)
     zone3a(5)= 0.181e0_dp  !for the generic Alfven test (RAW)
     zone3a(6)= 5.50e0_dp   !for the generic Alfven test (RAW)
     zone3a(7)= 0.8195e0_dp !for the generic Alfven test (RAW)
  elseif (initial_data==12) then
     zone3a(3)= -0.113e0_dp   !for Balsara5 (RAW)
     zone3a(4)= -0.0462e0_dp  !for Balsara5 (RAW)
     zone3a(5)= 0.16e0_dp     !for Balsara5 (RAW)
     zone3a(6)= -1.429e0_dp   !for Balsara5 (RAW)
     zone3a(7)= 0.732e0_dp    !for Balsara5 (RAW)
  elseif ((initial_data==0).AND.(alfvenwave=='y')) then
     zone3a(3)= rightalfven(1) !for user defined Riemann problem (RAW)
     zone3a(4)= rightalfven(2) !for user defined Riemann problem (RAW)
     zone3a(5)= rightalfven(3) !for user defined Riemann problem (RAW)
     zone3a(6)= rightalfven(4) !for user defined Riemann problem (RAW)
     zone3a(7)= rightalfven(5) !for user defined Riemann problem (RAW)
  end if

  !solve the system of equations (4.15)-(4.17),(4.19)-(4.20)
  call alfven(VsRv3(2),zone3b,zone3a)
  
  if (wave_error) stop 'error at the right-going alfven discontinuity'

  !RIGHT GOING SLOW WAVE
  if (verbose) then
     print *
     print *,'--------------------------------------'
     print *,'I''m solving equations for the'
     print *,'Right Going Slow Wave'
     print *,'--------------------------------------'
     print *
  end if

  normB=sqrt(unk(2)**2+unk(3)**2)
  ! Is it a shock?
  ! If it's a slow shock normB decreases
  if (abs(normB-sqrt(zone3a(6)**2+zone3a(7)**2))<epsilon(unk(3))) then
     if (verbose) print *,'no shock or rarefaction!'
     zone2b=zone3a
     goto 21
  end if

  if (normB<=sqrt(zone3a(6)**2+zone3a(7)**2)) then
     if (verbose) print *,'It''s a shock!'
     !SHOCK
     !compute the velocity of the right going slow shock
     call velocity(unk(2),unk(3),zone3a,'RS',vxb,VsRv3(3),'B',zone2b)
     if (wave_error) stop 'fullfuncv: error at the right-going shock wave'
  else
     !RAREFACTION
     if (verbose) print *,'It''s a rarefaction!'

     if(abs(zone3a(7)/zone3a(6)-unk(3)/unk(2))<small_psi) then
        !psi is ocnstant
        if (verbose) print *,'psi is constant'
        call Rarefaction(normB,VsRv3(3),zone3a,zone2b,'RS')
     else
        !psi is not constant
        if (verbose) print *,'psi is NOT constant'
        psi=atan(unk(3)/unk(2))
        if (((unk(3)*unk(2))<0.0e0_dp) .AND. (unk(3)>0.0e0_dp)) then
           psi=acos(unk(2)/sqrt(unk(2)**2+unk(3)**2))
        end if
        call Rarefaction_psi(psi,VsRv3(3),zone3a,zone2b,'RS')
     end if
     if (wave_error) stop 'fullfuncv: error at the right-going rarefaction wave'
  end if



  !Values of the velocity and pressure at the right side of the contact discontinuity
21 vx2=zone2b(3)
  vy2=zone2b(4)
  vz2=zone2b(5)
  ptot2=zone2b(2)
  
  fvec(1)=vx1-vx2
  fvec(2)=vy1-vy2
  fvec(3)=vz1-vz2
  fvec(4)=ptot1-ptot2
  
  if (present(fullsolution)) then
     fullsolution(:,1)=zone1a !region R2
     fullsolution(:,2)=zone1b !region R3
     fullsolution(:,3)=zone2a !region R4
     fullsolution(:,4)=zone2b !region R5
     fullsolution(:,5)=zone3a !region R6
     fullsolution(:,6)=zone3b !region R7

     write(200,'(I3,6E20.8)') niter,vx1-vx2,vy1-vy2,vz1-vz2,zone2a(6)-zone2b(6),zone2a(7)-zone2b(7),ptot1-ptot2
  end if


end subroutine fullfuncv


!---------------------------------------------------------------!
!-------------------------OUTPUT--------------------------------!
!---------------------------------------------------------------!

subroutine output(x1,x2,t,nx,left,right,VsLv3,VsRv3,unk,fullsolution)
  use type
  use global !This module contains Bx, the EOS' gamma, degen and niter
  use interfaces,only:xi,Rarefaction
  implicit none
  real(DP),intent(IN)::x1,x2,t
  real(DP),dimension(4),intent(IN)::unk
  real(DP),dimension(3),intent(IN)::VsLv3,VsRv3
  integer(I4B),intent(IN)::nx
  real(DP),dimension(7),intent(IN)::left,right
  real(DP),dimension(7,6),intent(IN)::fullsolution

  !This subroutine writes the solution to file 'solution.dat'
  !NOTE: the initial discontinuity is located at (x2-x1)/2

  real(DP),dimension(nx)::x
  real(DP)::xc,vxc,pbc,h,pb,xir,xp,VsLeft,VsRight,normb,psi
  real(DP),dimension(7,2)::solution
  real(DP),dimension(7)::solutionr
  integer(I4B)::i,j,nxr
  real(DP),parameter::small_psi=1.0e-4_dp

  nxr=nx/5 !Number of points used for rarefaction waves

  xc=0.5e0_dp*(x2+x1)!position of the initial interface

  if ((Bx==0.0e0_dp).OR.(degen)) then
     pbc=unk(1)
     solution(:,1)=fullsolution(:,3)
     solution(:,2)=fullsolution(:,4)
     vxc=0.5e0_dp*(solution(3,1)+solution(3,2)) !Velocity of the contact discontinuity
  else
     vxc=0.5e0_dp*(fullsolution(3,3)+fullsolution(3,4)) !Velocity of the contact discontinuity
     pbc=0.5e0_dp*(fullsolution(2,3)+fullsolution(2,4)) !pressure at the contact discontinuity
  end if

  write(1223,*) '#iteration number',niter
  write(1223,*) 'These are the values I have found:'
  write(1223,*) 'vx at the CD     =',vxc
  write(1223,*) 'p at the CD      =',pbc
  write(1223,*)
  write(1223,*) 'left-going  fast wave velocity =',VsLv3(1)
  if (.NOT.((Bx==0.0e0_dp).OR.(degen))) write(1223,*) 'left-going  alfven velocity    =',VsLv3(2)
  if (.NOT.((Bx==0.0e0_dp).OR.(degen))) write(1223,*) 'left-going  slow wave velocity =',VsLv3(3)
  if (.NOT.((Bx==0.0e0_dp).OR.(degen))) write(1223,*) 'right-going slow wave velocity =',VsRv3(3)
  if (.NOT.((Bx==0.0e0_dp).OR.(degen))) write(1223,*) 'right-going alfven velocity    =',VsRv3(2)
  write(1223,*) 'right-going fast wave velocity =',VsRv3(1)
  write(1223,*)
  if (.NOT.((Bx==0.0e0_dp).OR.(degen))) then
     write(1223,*) 'total pressure in regions:'
     write(1223,*) 'R1    =',left(2)
     write(1223,*) 'R2-R3 =',fullsolution(2,1)
     write(1223,*) 'R4-R5 =',fullsolution(2,3)
     write(1223,*) 'R6-R7 =',fullsolution(2,5)
     write(1223,*) 'R8    =',right(2)
     write(1223,*)
     write(1223,*) 'Values used to get the solution:'
     write(1223,*) 'total pressure in regions R2-R3  =',unk(1)
     write(1223,*) 'By             in regions R4-R5  =',unk(2)
     write(1223,*) 'Bz             in regions R4-R5  =',unk(3)
     write(1223,*) 'total pressure in regions R6-R7  =',unk(4)
  else
     write(1223,*) 'total pressure in regions:'
     write(1223,*) 'R1    =',left(2)
     write(1223,*) 'R2-R3 =',pbc
     write(1223,*) 'R4    =',right(2)
     write(1223,*)
     write(1223,*) 'Value used to get the solution:'
     write(1223,*) 'total pressure in regions R2-R3  =',unk(1)
  end if

  write(1223,*)
  write(1223,*) 'Solution plotted in file ''solution.dat'''
  write(1223,*)
  write(1223,*)


  !grid
  write(100,*)
  write(100,*) '#iteration number',niter
  write(100,*) '# x, rho, total pressure p, vx, vy, vz, By, Bz'

  do i=1,nx
     x(i)=x1+(i-1)*(x2-x1)/(nx-1)
  end do



  if ((Bx==0.0e0_dp).OR.(degen)) then
     !Only two fast waves and a tangential discontinuity
     VsLeft=VsLv3(1)
     VsRight=VsRv3(1)


     i=1
     do while(x(i)-xc<=VsLeft*t)
        write(100,FMT='(8E20.9)') x(i),left
        i=i+1
        xp=x(i)
        if (i>nx) then
           return
        end if
     end do
     write(100,FMT='(8E20.9)') xc+VsLeft*t,left !EXACT POSITION OF THE SHOCK

     if (pbc>left(2)) then !Shock
        write(100,FMT='(8E20.9)') xc+VsLeft*t,solution(:,1) !EXACT POSITION OF THE SHOCK
        do while(x(i)-xc<=vxc*t)
           write(100,FMT='(8E20.9)') x(i),solution(:,1)
           i=i+1
           xp=x(i)
           if (i>nx) then
              return
           end if
        end do
        write(100,FMT='(8E20.9)') xc+vxc*t,solution(:,1) !EXACT POSITION OF THE SHOCK
     else !Rarefaction
        pb=left(2)
        h=(pbc-left(2))/nxr
        solutionr=left
        pb=pb+h
        do j=1,nxr
           call Rarefaction_pressure(pb,xir,left,solutionr,'LF')
           !Rarefaction returns always the value of xi at left(2)
           call xi(solutionr,'LF',xir)!Computes the value of xi at pb
           write(100,FMT='(8E20.9)') xc+xir*t,solutionr
           xp=xc+xir*t
           pb=pb+h
           if (xp>x2) return
23         if (x(i+1)<=xp) then
              i=i+1
              goto 23
           end if
        end do
        i=i+1
        do while(x(i)-xc<=vxc*t)
           write(100,FMT='(8E20.9)') x(i),solutionr
           i=i+1
           xp=x(i)
           if (i>nx) then
              return
           end if
        end do
        write(100,FMT='(8E20.9)') xc+vxc*t,solution(:,1) !EXACT POSITION OF THE SHOCK
     end if

     if (pbc>right(2)) then !Shock
        write(100,FMT='(8E20.9)') xc+vxc*t,solution(:,2) !EXACT POSITION OF THE SHOCK
        do while(x(i)-xc<=VsRight*t)
           write(100,FMT='(8E20.9)') x(i),solution(:,2)
           xp=x(i)
           i=i+1
           if (i>nx) then
              return
           end if
        end do
        write(100,FMT='(8E20.9)') xc+VsRight*t,solution(:,2) !EXACT POSITION OF THE SHOCK
        write(100,FMT='(8E20.9)') xc+VsRight*t,right !EXACT POSITION OF THE SHOCK
     else !Rarefaction
        write(100,FMT='(8E20.9)') xc+vxc*t,solution(:,2) !EXACT POSITION OF THE SHOCK
        call xi(solution(:,2),'RF',xir)!compute the value of xi at pb
        do while(x(i)<xc+xir*t)
           write(100,FMT='(8E20.9)') x(i),solution(:,2)
           xp=x(i)
           i=i+1
           if (i>nx) then
              return
           end if
        end do

        pb=pbc
        h=(right(2)-pb)/nxr
        solutionr=solution(:,2)      
        pb=pb+h
        do j=1,nxr
           call Rarefaction_pressure(pb,xir,solution(:,2),solutionr,'RF')
           !Rarefaction returns always the value of xi at right(2)
           call xi(solutionr,'RF',xir)!compute the value of xi at pb
           write(100,FMT='(8E20.9)') xc+xir*t,solutionr
           pb=pb+h
           xp=xc+xir*t
           if (xp>x2) return
        end do

     end if

     do 
        if (x(i)>=xc+VsRight*t) write(100,FMT='(8E20.9)') x(i),right
        i=i+1
        if (i>nx) then
           return
        end if
     end do

  else

     !Bx different from zero
     !ALL THE WAVES ARE PRESENT

     i=1
     do while(x(i)-xc<VsLv3(1)*t)
        write(100,FMT='(8E20.9)') x(i),left
        i=i+1
        xp=x(i)
        if (i>nx) then
           return
        end if
     end do
     if ((xc+VsLv3(1)*t)>x1) then
        write(100,FMT='(8E20.9)') xc+VsLv3(1)*t,left !EXACT POSITION OF THE SHOCK
     end if


     !---------------------------------!
     !------LEFT GOING FAST WAVE-------!
     !---------------------------------!

     print *,'PLOTTING... Left Going Fast Wave'

     if (fullsolution(2,1)>=left(2)) then !Shock
        if ((xc+VsLv3(1)*t)>x1) then
           write(100,FMT='(8E20.9)') xc+VsLv3(1)*t,fullsolution(:,1) !EXACT POSITION OF THE SHOCK
        end if

        do while(x(i)-xc<=VsLv3(2)*t)
           write(100,FMT='(8E20.9)') x(i),fullsolution(:,1)
           i=i+1
           xp=x(i)
           if (i>nx) then
              return
           end if
        end do
        write(100,FMT='(8E20.9)') xc+VsLv3(2)*t,fullsolution(:,1) !EXACT POSITION OF THE SHOCK
     else !Rarefaction
        pb=left(2)
        h=(unk(1)-pb)/nxr
        solutionr=left
        pb=pb+h
        do j=1,nxr
           call Rarefaction_pressure(pb,xir,left,solutionr,'LF')
           !Rarefaction returns always the value of xi at left(2)
           call xi(solutionr,'LF',xir)!Computes the value of xi at pb
           write(100,FMT='(8E20.9)') xc+xir*t,solutionr
           xp=xc+xir*t
           pb=pb+h
           if (xp>x2) return
233        if (x(i+1)<=xp) then
              i=i+1
              goto 233
           end if
        end do
        i=i+1
        do while(x(i)-xc<=VsLv3(2)*t)
           write(100,FMT='(8E20.9)') x(i),solutionr
           i=i+1
           xp=x(i)
           if (i>nx) return
        end do
        write(100,FMT='(8E20.9)') xc+VsLv3(2)*t,solutionr !EXACT POSITION OF THE SHOCK
     end if




     !---------------------------------!
     !------LEFT GOING ALFVEN WAVE-----!
     !---------------------------------!

     print *,'PLOTTING... Left Going Alfven Wave'

     write(100,FMT='(8E20.9)') xc+VsLv3(2)*t,fullsolution(:,2) !EXACT POSITION OF THE SHOCK
     do while(x(i)-xc<=VsLv3(3)*t)
        write(100,FMT='(8E20.9)') x(i),fullsolution(:,2)
        i=i+1
        xp=x(i)
        if (i>nx) return
     end do
     write(100,FMT='(8E20.9)') xc+VsLv3(3)*t,fullsolution(:,2) !EXACT POSITION OF THE SHOCK


     !---------------------------------!
     !------LEFT GOING SLOW WAVE-------!
     !---------------------------------!
     
     print *,'PLOTTING... Left Going Slow Wave'
     
     if (fullsolution(2,3)>=fullsolution(2,2)) then !Shock
        write(100,FMT='(8E20.9)') xc+VsLv3(3)*t,fullsolution(:,3) !EXACT POSITION OF THE SHOCK
        do while(x(i)-xc<=vxc*t)
           write(100,FMT='(8E20.9)') x(i),fullsolution(:,3)
           i=i+1
           xp=x(i)
           if (i>nx) return
        end do
        write(100,FMT='(8E20.9)') xc+vxc*t,fullsolution(:,3) !EXACT POSITION OF THE SHOCK
     else !Rarefaction
        solutionr=fullsolution(:,2)
        if (abs((fullsolution(7,2)/fullsolution(6,2))-(unk(3)/unk(2)))<=small_psi) then
           !psi is constant
           normb=sqrt(fullsolution(6,2)**2+fullsolution(7,2)**2)
           h=(sqrt(unk(2)**2+unk(3)**2)-normb)/nxr
           normb=normb+h
           if (veryverbose) print *,'psi is constant'
        else
           !psi is not constant
           psi=atan(fullsolution(7,2)/fullsolution(6,2))
           h=(atan(unk(3)/unk(2))-psi)/nxr
           if (((fullsolution(7,2)*fullsolution(6,2))<0.0e0_dp) .AND. (fullsolution(7,2)>0.0e0_dp)) then
              psi=acos(fullsolution(6,2)/sqrt(fullsolution(6,2)**2+fullsolution(7,2)**2))
              h=(acos(unk(2)/sqrt(unk(2)**2+unk(3)**2))-psi)/nxr
           end if
           psi=psi+h
        end if
        if (verbose) print*,'Rarefaction being called'
        do j=1,nxr
           !  Modified by M.A. Aloy (20-01-2009):
           if (abs((fullsolution(7,2)/fullsolution(6,2))-(unk(3)/unk(2)))<=small_psi) then
              !psi is constant
              if (veryverbose) print *,'psi is constant'
              call Rarefaction(normb,xir,fullsolution(:,2),solutionr,'LS')
              normb=normb+h
           else
              !psi is not constant
              if (veryverbose) print *,'psi is NOT constant'
              call Rarefaction_psi(psi,xir,fullsolution(:,2),solutionr,'LS')
              psi=psi+h
           end if
           !
           !  End modification
           !
           !Rarefaction returns always the value of xi at left(2)           
           call xi(solutionr,'LS',xir)!Computes the value of xi at pb
           write(100,FMT='(8E20.9)') xc+xir*t,solutionr
           xp=xc+xir*t
           
           if (xp>x2) return
2334       if (x(i+1)<=xp) then
              i=i+1
              goto 2334
           end if
        end do
        i=i+1
        do while(x(i)-xc<=vxc*t)
           write(100,FMT='(8E20.9)') x(i),fullsolution(:,3)
           i=i+1
           xp=x(i)
           if (i>nx) return
        end do
        write(100,FMT='(8E20.9)') xc+vxc*t,fullsolution(:,3) !EXACT POSITION OF THE SHOCK
     end if


     !---------------------------------!
     !------RIGHT GOING SLOW WAVE------!
     !---------------------------------!

     print *,'PLOTTING... Right Going Slow Wave'

     if (fullsolution(2,4)>=fullsolution(2,5)) then !Shock
        write(100,FMT='(8E20.9)') xc+vxc*t,fullsolution(:,4) !EXACT POSITION OF THE SHOCK
        do while(x(i)-xc<=VsRv3(3)*t)
           write(100,FMT='(8E20.9)') x(i),fullsolution(:,4)
           i=i+1
           xp=x(i)
           if (i>nx) return
        end do
        write(100,FMT='(8E20.9)') xc+VsRv3(3)*t,fullsolution(:,4) !EXACT POSITION OF THE SHOCK
     else !Rarefaction
        call xi(fullsolution(:,4),'RS',xir)!Computes the value of xi at pb
        do while(x(i)-xc<=xir*t)
           write(100,FMT='(8E20.9)') x(i),fullsolution(:,4)
           i=i+1
           xp=x(i)
           if (i>nx) then
              return
           end if
        end do

        normb=sqrt(unk(2)**2+unk(3)**2)
        h=(sqrt(fullsolution(6,5)**2+fullsolution(7,5)**2)-normb)/nxr
        solutionr=fullsolution(:,4)
        normb=normb+h
        do j=1,nxr
           call Rarefaction(normb,xir,fullsolution(:,4),solutionr,'RS')
           !Rarefaction returns always the value of xi at left(2)
           call xi(solutionr,'RS',xir)!Computes the value of xi at pb
           write(100,FMT='(8E20.9)') xc+xir*t,solutionr
           xp=xc+xir*t
           normb=normb+h
           if (xp>x2) return
2335       if (x(i+1)<=xp) then
              i=i+1
              goto 2335
           end if
        end do
     end if


     !---------------------------------!
     !------RIGHT GOING ALFVEN WAVE----!
     !---------------------------------!

     print *,'PLOTTING... Right Going Alfven Wave'
     write(100,FMT='(8E20.9)') xc+VsRv3(3)*t,fullsolution(:,5) !EXACT POSITION OF THE SHOCK
     do while(x(i)-xc<=VsRv3(2)*t)
        write(100,FMT='(8E20.9)') x(i),fullsolution(:,5)
        i=i+1
        xp=x(i)
        if (i>nx) return
     end do
     write(100,FMT='(8E20.9)') xc+VsRv3(2)*t,fullsolution(:,5) !EXACT POSITION OF THE SHOCK

     !---------------------------------!
     !------RIGHT GOING FAST WAVE------!
     !---------------------------------!

     print *,'PLOTTING... Right Going Fast Wave'


     if (fullsolution(2,5)>=right(2)) then !Shock
        write(100,FMT='(8E20.9)') xc+VsRv3(2)*t,fullsolution(:,6) !EXACT POSITION OF THE SHOCK
        do while(x(i)-xc<=VsRv3(1)*t)
           write(100,FMT='(8E20.9)') x(i),fullsolution(:,6)
           xp=x(i)
           i=i+1
           if (i>nx) then
              return
           end if
        end do
        write(100,FMT='(8E20.9)') xc+VsRv3(1)*t,fullsolution(:,6) !EXACT POSITION OF THE SHOCK
        write(100,FMT='(8E20.9)') xc+VsRv3(1)*t,right !EXACT POSITION OF THE SHOCK
     else !Rarefaction
        call xi(fullsolution(:,6),'RF',xir)!compute the value of xi at pb
        do while(x(i)<xc+xir*t)
           write(100,FMT='(8E20.9)') x(i),fullsolution(:,6)
           xp=x(i)
           i=i+1
           if (i>nx) then
              return
           end if
        end do

        pb=unk(4)
        h=(right(2)-pb)/nxr
        solutionr=fullsolution(:,6)      
        pb=pb+h
        do j=1,nxr
           call Rarefaction_pressure(pb,xir,fullsolution(:,6),solutionr,'RF')
           !Rarefaction returns always the value of xi at right(2)
           call xi(solutionr,'RF',xir)!compute the value of xi at pb
           write(100,FMT='(8E20.9)') xc+xir*t,solutionr
           pb=pb+h
           xp=xc+xir*t
           if (xp>x2) return
        end do
     end if

     do 
        if(x(i)>xp) write(100,FMT='(8E20.9)') x(i),right
        i=i+1
        if (i>nx) return
     end do

  end if


end subroutine output


!---------------------------------------------------------------!
!----------------------RAREFACTION------------------------------!
!---------------------------------------------------------------!

subroutine Rarefaction(normB,Vs,init,solution,switchLR)
  use type
  use interfaces,only:RHS,xi
  use odeswitch
  use wavecheck
  use global !This module contains verbose
  implicit none
  real(DP),intent(IN)::normB !This is sqrt(By+Bz) at the tail of the rarefaction (i.e. at the CD)
  real(DP),intent(OUT)::Vs
  real(DP),dimension(7),intent(IN)::init
  real(DP),dimension(7),intent(OUT)::solution
  character(len=2),intent(IN)::switchLR

  !Solve the system of ODEs for rarefaction waves in the Bt-method if By/Bz is constant across the wave
  !Use a fourth order Runge-Kutta with fixed stepsize

  real(DP),dimension(7)::y,k1,k2,k3,k4,f
  real(DP)::x,h
  integer(I4B)::j
  integer(I4B),PARAMETER::nsteps=10000 !number of steps

  switch=switchLR !Contained in module odeswitch and used in odeint subroutines

  !Solve the ODE from normB(ahead) to normB
  x=sqrt(init(6)**2+init(7)**2) !sqrt(Bya**2+Bza**2)
  y=init

  !Head velocity
  call xi(y,switchLR,Vs)

  if (veryverbose) print *,'dnormB=',abs(normB-x)

  if (abs(normB-x)<=epsilon(normB)) then
     if (veryverbose) print *,'normB doesn''t change!'
     solution=init
     return
  end if

  !step
  h=(normB-x)/nsteps !use nsteps


  !Fourth order Runge-Kutta with a fixed step
  do j=1,nsteps
     call RHS(x,y,f,switchLR)
     if (wave_error) return 
     k1=h*f
     call RHS(x+0.5e0_dp*h,y+0.5e0_dp*k1,f,switchLR)
     if (wave_error) return 
     k2=h*f
     call RHS(x+0.5e0_dp*h,y+0.5e0_dp*k2,f,switchLR)
     if (wave_error) return 
     k3=h*f
     call RHS(x+h,y+k3,f,switchLR)
     k4=h*f

     y=y+(k1/6.0e0_dp)+(k2/3.0e0_dp)+(k3/3.0e0_dp)+(k4/6.0e0_dp)
     x=x+h

     if (wave_error) return 

  end do

  solution=y


end subroutine Rarefaction


!---------------------------------------------------------------!
!------------------------------RHS------------------------------!
!---------------------------------------------------------------!
subroutine RHS(x,y,f,switchLR)
  use type
  use interfaces,only:xi,eos_enthalpy,eos_cs2
  use global !This module contains Bx and gamma
  use wavecheck !This module contains wave_check
  implicit none
  real(DP),intent(IN)::x
  real(DP),dimension(7),intent(IN)::y
  real(DP),dimension(7),intent(OUT)::f
  character(len=2),intent(IN)::switchLR

  !subroutine called by Rarefaction
  !It returns the Right Hand Side of the ODE equations for rarefaction waves
  ! (used for slow rarefactions when the angle psi=atan(Bz/By) is constant)

  real(DP)::psi,dV,Vs,rho,p,vx,vy,vz,By,Bz,normB,dBy,dBz,W,eta,B2,b2small,wtot,h,cs2


  !NOTE
  !y(1) = rho
  !y(2) = total pressure P
  !y(3) = vx
  !y(4) = vy
  !y(5) = vz
  !y(6) = By
  !y(7) = Bz

  !x    = normB
  normB=x

  rho = y(1)
  P   = y(2)
  vx  = y(3)
  vy  = y(4)
  vz  = y(5)
  By  = y(6)
  Bz  = y(7)

  !Lorentz factor
  W = 1.0e0_dp/sqrt(1.0e0_dp-vx**2-vy**2-vz**2)

  eta = Bx*vx + By*vy + Bz*vz !B^j v_j

  B2 = Bx**2 + By**2 + Bz**2 !B^2

  b2small = B2/W**2 + eta**2 !b^mu b_mu

  !relativistic specific enthalpy
  call eos_enthalpy(p-0.5e0_dp*b2small,rho,gamma,h)
  !total relativistic enthalpy
  wtot = rho*h + b2small


  !sound velocity squared
  call eos_cs2(p-0.5e0_dp*b2small,rho,gamma,cs2)


  !angle psi
  psi=atan(Bz/By)

  if ((normB>0).AND.(Bz==0.0e0_dp)) then
     psi=acos(By/normB)
  end if

  if ((By*Bz<0.0e0_dp).AND.(Bz>0)) then
     psi=acos(By/normB)
  end if

  !compute \xi
  call xi(y,switchLR,Vs)

  !I'm supposing psi constant
  dBy=cos(psi)
  dBz=sin(psi)


  !drho/dnormB
  f(1) = -((rho*(-((Bz*dBy - By*dBz)*eta*(Bz*vy - By*vz)*W**4*(vx - Vs)**2*(-1 + vx*Vs)) + & 
       Bx**3*(dBy*vy + dBz*vz)*W**2*(vx - Vs)*(-1 + W**2*(vx - Vs)*Vs) + &
       Bx**2*(eta*(dBy*vy + dBz*vz)*W**2* &
       (-1 - W**2*Vs**2 + vx**2*W**2*(-2 + Vs**2) + & 
       vx*(Vs + 3*W**2*Vs - W**2*Vs**3)) + &
       By*(dBz*vy*vz*W**4*(vx - Vs)**2 + &
       dBy*(-1 + vx**2*(vx**2 + vy**2)*W**4 - &
       vx*(-1 + W**2*(1 + 2*vy**2*W**2 + vx**2*(-1 + 3*W**2)))* &
       Vs + W**2*(1 + vy**2*W**2 + vx**2*(-1 + 3*W**2))*Vs**2 - &
       vx*W**4*Vs**3)) + &
       Bz*(dBy*vy*vz*W**4*(vx - Vs)**2 + &
       dBz*(-1 + vx**2*(vx**2 + vz**2)*W**4 - &
       vx*(-1 + W**2*(1 + 2*vz**2*W**2 + vx**2*(-1 + 3*W**2)))* &
       Vs + W**2*(1 + vz**2*W**2 + vx**2*(-1 + 3*W**2))*Vs**2 - &
       vx*W**4*Vs**3))) + &
       Bx*W**2*(vx - Vs)*(-((dBy*vy + dBz*vz)* &
       (W**4*(eta**2 - wtot)*(vx - Vs)**2 + B2*(-1 + vx*Vs))) + &
       By**2*(-(dBz*vz*W**2*(vx - Vs)**2) + dBy*vy*(-1 + vx*Vs)) + &
       Bz**2*(-(dBy*vy*W**2*(vx - Vs)**2) + dBz*vz*(-1 + vx*Vs)) + &
       Bz*eta*(-1 + vx*Vs)*(dBy*vy*vz*W**2 + dBz*(1 + W**2*(vx**2 + vz**2 - vx*Vs))) + & 
       By*(dBz*vy*(eta*vz*W**2*(-1 + vx*Vs) + &
       Bz*(-1 + vx**2*W**2 + W**2*Vs**2 + vx*(Vs - 2*W**2*Vs))) &
       + dBy*(eta*(-1 + vx*Vs)* &
       (1 + W**2*(vx**2 + vy**2 - vx*Vs)) + &
       Bz*vz*(-1 + vx**2*W**2 + W**2*Vs**2 + vx*(Vs - 2*W**2*Vs)))))))/ &
       (Bx*W**2*(Bx**3*(vx - Vs)*(vx + Vs) - &
       (By**2*eta - By*vy*(B2 + cs2*h*rho*W**2) + &
       Bz*(Bz*eta - vz*(B2 + cs2*h*rho*W**2)))*(vx - Vs)*(-1 + vx*Vs) &
       - Bx**2*(vx - Vs)*(eta - 2*(By*vy + Bz*vz) + eta*vx*Vs) + &
       Bx*(cs2*h*rho*(1 + vx*W**2*(vx - Vs))*(-1 + vx*Vs) + &
       (vx - Vs)*(-((By**2 + Bz**2 + W**2*(eta**2 - wtot))* &
       (vx - Vs)) + B2*vx*(-1 + vx*Vs))))))

  !dp/dnormB
  f(2) = ((vx - Vs)*(W**2*(B2*(Bz**2*dBy*vy + By**2*dBz*vz - &
       By*Bz*(dBz*vy + dBy*vz) + &
       (dBy*vy + dBz*vz)*W**2*(eta**2 - wtot)) + &
       W**2*(Bz**2*cs2*dBy*h*rho*vy + By**2*cs2*dBz*h*rho*vz + &
       cs2*h*rho*(dBy*vy + dBz*vz)*W**2*(eta**2 - wtot) + &
       Bz*dBz*eta*(-eta**2 + wtot) - &
       By*(Bz*cs2*dBz*h*rho*vy + &
       dBy*(eta**3 + Bz*cs2*h*rho*vz - eta*wtot))))*(vx - Vs)**2 + &
       Bx**3*(2*eta*(dBy*vy + dBz*vz)*W**2*Vs + By*dBy*(vx + Vs) + &
       Bz*dBz*(vx + Vs)) + &
       Bx**2*(2*By**2*dBy*vy + 2*Bz**2*dBz*vz + &
       (dBy*vy + dBz*vz)*W**2* &
       (-2*eta**2 + B2*Vs*(-vx + Vs) + &
       cs2*h*rho*(1 + W**2*Vs*(-vx + Vs))) + &
       2*Bz*eta*(dBy*vy*vz*W**2 + &
       dBz*(-1 + W**2*(vx**2 + vz**2 - vx*Vs))) + &
       2*By*(dBz*vy*(Bz + eta*vz*W**2) + &
       dBy*(Bz*vz + eta*(-1 + W**2*(vx**2 + vy**2 - vx*Vs))))) - &
       Bx*(vx - Vs)*(By**3*dBy + Bz**3*dBz + 2*Bz**2*dBy*eta*vy*W**2 + &
       By**2*dBz*(Bz + 2*eta*vz*W**2) + &
       eta*(dBy*vy + dBz*vz)*W**2* &
       (-2*B2 + W**2*(eta**2 - 2*cs2*h*rho - wtot)) + & 
       Bz*W**2*(dBy*vy*vz*(B2 + cs2*h*rho*W**2) + &
       dBz*(3*eta**2 - wtot + B2*(vx**2 + vz**2 - vx*Vs) + &
       cs2*h*rho*(1 + W**2*(vx**2 + vz**2 - vx*Vs)))) + &
       By*(Bz**2*dBy - 2*Bz*eta*(dBz*vy + dBy*vz)*W**2 + &
       W**2*(dBz*vy*vz*(B2 + cs2*h*rho*W**2) + &
       dBy*(3*eta**2 - wtot + B2*(vx**2 + vy**2 - vx*Vs) + &
       cs2*h*rho*(1 + W**2*(vx**2 + vy**2 - vx*Vs))))))))/ &
       (W**2*(Bx**3*(vx - Vs)*(vx + Vs) - &
       (By**2*eta - By*vy*(B2 + cs2*h*rho*W**2) + &
       Bz*(Bz*eta - vz*(B2 + cs2*h*rho*W**2)))*(vx - Vs)*(-1 + vx*Vs) - &
       Bx**2*(vx - Vs)*(eta - 2*(By*vy + Bz*vz) + eta*vx*Vs) + &
       Bx*(cs2*h*rho*(1 + vx*W**2*(vx - Vs))*(-1 + vx*Vs) + &
       (vx - Vs)*(-((By**2 + Bz**2 + W**2*(eta**2 - wtot))*(vx - Vs)) + &
       B2*vx*(-1 + vx*Vs)))))

  !dvx/dnormB
  f(3) = ((vx - Vs)*(Bx**2*(dBy*vy + dBz*vz)*W**2*(-vx + Vs) - &
       W**2*(-(By*dBy*eta) - Bz*dBz*eta + &
       (dBy*vy + dBz*vz)*(B2 + cs2*h*rho*W**2))*(vx - Vs)*(-1 + vx*Vs) + &
       Bx*(eta*(dBy*vy + dBz*vz)*W**2*(-1 + vx*Vs) + &
       Bz*dBz*(-1 + vx**2*W**2 + vx*Vs - 2*vx*W**2*Vs + W**2*Vs**2) + &
       By*dBy*(-1 + vx**2*W**2 + W**2*Vs**2 + vx*(Vs - 2*W**2*Vs)))))/ &
       (W**2*(Bx**3*(vx - Vs)*(vx + Vs) - &
       (By**2*eta - By*vy*(B2 + cs2*h*rho*W**2) + &
       Bz*(Bz*eta - vz*(B2 + cs2*h*rho*W**2)))*(vx - Vs)*(-1 + vx*Vs) - &
       Bx**2*(vx - Vs)*(eta - 2*(By*vy + Bz*vz) + eta*vx*Vs) + &
       Bx*(cs2*h*rho*(1 + vx*W**2*(vx - Vs))*(-1 + vx*Vs) + &
       (vx - Vs)*(-((By**2 + Bz**2 + W**2*(eta**2 - wtot))*(vx - Vs)) + &
       B2*vx*(-1 + vx*Vs)))))
  !dvy/dnormB
  f(4) = ((vx - Vs)*(dBy - (By*(Bx**2*(dBy*vy + dBz*vz)*W**2*(-vx + Vs) - &
       W**2*(-(By*dBy*eta) - Bz*dBz*eta + &
       (dBy*vy + dBz*vz)*(B2 + cs2*h*rho*W**2))*(vx - Vs)* &
       (-1 + vx*Vs) + &
       Bx*(eta*(dBy*vy + dBz*vz)*W**2*(-1 + vx*Vs) + &
       Bz*dBz*(-1 + vx**2*W**2 + vx*Vs - 2*vx*W**2*Vs + &
       W**2*Vs**2) + &
       By*dBy*(-1 + vx**2*W**2 + W**2*Vs**2 + vx*(Vs - 2*W**2*Vs)))))/ &
       (W**2*((By**2*eta - By*vy*(B2 + cs2*h*rho*W**2) + &
       Bz*(Bz*eta - vz*(B2 + cs2*h*rho*W**2)))*(vx - Vs)* &
       (-1 + vx*Vs) + &
       Bx**2*(vx - Vs)*(eta - 2*(By*vy + Bz*vz) + eta*vx*Vs) + &
       Bx**3*(-vx**2 + Vs**2) + &
       Bx*(-(cs2*h*rho*(1 + vx*W**2*(vx - Vs))*(-1 + vx*Vs)) - &
       (vx - Vs)*(-((By**2 + Bz**2 + W**2*(eta**2 - wtot))* &
       (vx - Vs)) + B2*vx*(-1 + vx*Vs)))))))/Bx

  !dvz/dnormB
  f(5) = ((vx - Vs)*(-(dBz*(By*(By*eta - vy*(B2 + cs2*h*rho*W**2))*(vx - Vs)* &
       (-1 + vx*Vs) + Bx**2*(vx - Vs)*(eta - 2*By*vy + eta*vx*Vs) + &
       Bx**3*(-vx**2 + Vs**2) + &
       Bx*(-(cs2*h*rho*(1 + vx*W**2*(vx - Vs))*(-1 + vx*Vs)) - &
       (vx - Vs)*(-((By**2 + W**2*(eta**2 - wtot))*(vx - Vs)) + &
       B2*vx*(-1 + vx*Vs))))) + &
       Bz*((Bx*(By*dBy*(-1 + vx*Vs) + Bz*dBz*(-1 + vx*Vs) + &
       (dBy*vy + dBz*vz)*W**2*(-eta + Bx*vx - Bx*Vs + eta*vx*Vs)))/ &
       W**2 - dBy*(vx - Vs)* &
       (By*(eta - eta*vx*Vs + Bx*(-vx + Vs)) + &
       vy*(2*Bx**2 + B2*(-1 + vx*Vs) + cs2*h*rho*W**2*(-1 + vx*Vs))))) &
       )/(Bx*(Bx**3*(vx - Vs)*(vx + Vs) - &
       (By**2*eta - By*vy*(B2 + cs2*h*rho*W**2) + &
       Bz*(Bz*eta - vz*(B2 + cs2*h*rho*W**2)))*(vx - Vs)*(-1 + vx*Vs) - &
       Bx**2*(vx - Vs)*(eta - 2*(By*vy + Bz*vz) + eta*vx*Vs) + &
       Bx*(cs2*h*rho*(1 + vx*W**2*(vx - Vs))*(-1 + vx*Vs) + &
       (vx - Vs)*(-((By**2 + Bz**2 + W**2*(eta**2 - wtot))*(vx - Vs)) + &
       B2*vx*(-1 + vx*Vs)))))


  !dBy/dnormB (I'm supposing that psi doesn't change)
  f(6)=dBy


  !dBz/dnormB (I'm supposing that psi doesn't change)
  f(7)=dBz

end subroutine RHS


!---------------------------------------------------------------!
!------------------------------XI-------------------------------!
!---------------------------------------------------------------!

subroutine xi(state,switchLR,Vs,allwaves)
  use type
  use interfaces,only:quartic
  use global !This module contains Bx and gamma
  implicit none
  real(DP),dimension(7),intent(IN)::state
  character(len=2),intent(IN)::switchLR
  real(DP),intent(OUT)::Vs !this is xi
  real(DP),dimension(4),intent(OUT),OPTIONAL::allwaves !output all the eigenvalues

  !Compute the eigenvalues of the system
  
  real(DP)::rho,P,vx,vy,vz,By,Bz,v2,W,b2,Pgas,wtot,h,cs2,B2big,eta,B
  real(DP)::leftAlfven,rightAlfven
  real(DP),dimension(0:4)::a
  real(DP),dimension(4)::rtr
  real(DP),dimension(0:1)::btilde
  real(DP)::eps2,btilde2,temp
  integer(I4B)::i,j


  rho=state(1)
  P=state(2)
  vx=state(3)
  vy=state(4)
  vz=state(5)
  By=state(6)
  Bz=state(7)


  v2=vx**2+vy**2+vz**2

  if (v2>1.0e0_dp) then
     write(*,*) 'vx=',vx
     write(*,*) 'vy=',vy
     write(*,*) 'vz=',vz
     write(*,*) 'v2=',v2
     stop 'riemann.f90: xi: v2>1.0!'
  end if

  B2big=Bx**2+By**2+Bz**2 !B^2
  B=sqrt(B2big)
  eta=Bx*vx+By*vy+Bz*vz !B^j v_j
  W=1.0e0_dp/sqrt(1.0e0_dp-v2) !Lorentz factor

  !-----------------------EoS---------------------------!
  !EoS (sound velocity)
  b2=B2big/W**2+eta**2 !b^mu b_mu
  Pgas=P-0.5e0_dp*b2 !gas pressure
  call eos_enthalpy(Pgas,rho,gamma,h)!specific relativistic enthalpy
  wtot=rho*h + b2 !total relativistic enthalpy
  call eos_cs2(Pgas,rho,gamma,cs2)
  !-----------------------------------------------------!

  !left-going alfven velocity
  leftAlfven  = vx+Bx*(1.0e0_dp-vx**2-vy**2-vz**2)/(eta-sqrt(wtot))
  !right-going alfven velocity
  rightAlfven = vx+Bx*(1.0e0_dp-vx**2-vy**2-vz**2)/(eta+sqrt(wtot))


  !The values for the rarefaction velocities are given by the solution
  !of the following quartic equation:
  !a(0)+a(1)*Vs+a(2)*Vs**2+a(3)*Vs**3+a(4)*Vs**4

  btilde(0)=W*eta
  btilde(1)=Bx/W+W*eta*vx
  btilde=btilde/sqrt(wtot)
  btilde2=(Bx**2+By**2+Bz**2)*(1.0e0_dp-vx**2-vy**2-vz**2)+eta**2
  btilde2=btilde2/wtot

  eps2=cs2+btilde2-cs2*btilde2

  a(0)=btilde(1)**2*cs2-eps2*W**2*vx**2-(eps2-1.0e0_dp)*W**4*vx**4
  a(1)=-2.0e0_dp*btilde(0)*btilde(1)*cs2+2.0e0_dp*eps2*W**2*vx+4.0e0_dp*(eps2-1.0e0_dp)*W**4*vx**3
  a(2)=(btilde(0)**2-btilde(1)**2)*cs2+W**2*(6.0e0_dp*W**2*vx**2+ &
       eps2*(-1.0e0_dp+(1.0e0_dp-6.0e0_dp*W**2)*vx**2))
  a(3)=2.0e0_dp*btilde(0)*btilde(1)*cs2+2.0e0_dp*W**2*(-eps2+2.0e0_dp*(-1.0e0_dp+eps2)*W**2)*vx
  a(4)=-btilde(0)**2*cs2+W**4+eps2*(W**2-W**4)

  if (Bx<=epsilon(Bx)) then
     !if Bx=0 the quartic reduces to a quadratic equation with the following coefficients
     !a(2)+a(3)*Vs+a(4)*Vs**2

     a(2)=-By ** 4 - Bz ** 4 - Bz*eta*vz*W ** 4*(eta ** 2 - wtot) - &
          W ** 4*(cs2*h*rho*(1 + vx ** 2*W ** 2) + &
          vx ** 2*(B ** 2 + W ** 2*(eta ** 2 - wtot)))*(eta ** 2 - wtot) - &
          Bz ** 2*W **2*(B ** 2*(vx ** 2 + vz ** 2) + eta ** 2*(1 + vx ** 2*W ** 2) + &
          cs2*h*rho*(1 + (vx ** 2 + vz ** 2)*W ** 2) - (1 + vx ** 2*W ** 2)* &
          wtot) - By*vy* W ** 2*(2*B ** 2*Bz*vz + &
          W ** 2*(eta ** 3 + 2*Bz*cs2*h*rho*vz - eta*wtot)) - &
          By ** 2*(2*Bz ** 2 + W ** 2*(B ** 2*(vx ** 2 + vy ** 2) + eta ** 2*(1 + vx ** 2*W ** 2) + &
          cs2*h*rho*(1 + (vx ** 2 + vy ** 2)*W ** 2) - (1 + vx ** 2*W ** 2)*wtot))

     a(3)=vx*(By ** 2 + Bz ** 2 + W ** 2*(B ** 2*(1 + vx ** 2 + vy ** 2 + vz ** 2) + &
          cs2*h*rho*(1 + (1 + vx ** 2 + vy ** 2 + vz ** 2)*W ** 2) + &
          2*W ** 2*(eta ** 2 - wtot)))*(By ** 2 + Bz ** 2 + &
          W ** 2*(eta ** 2 - wtot))

     a(4)=W ** 2*(B ** 2*(-(Bz ** 2*(vx ** 2 + vy ** 2)) + 2*By*Bz*vy*vz - &
          By ** 2*(vx ** 2 + vz ** 2) + (vx ** 2 + vy ** 2 + vz ** 2)* &
          W ** 2*(-eta ** 2 + wtot)) + W ** 2*(Bz*eta*vz*(eta ** 2 - wtot) - &
          Bz ** 2*(eta ** 2 + cs2*h*rho*(vx ** 2 + vy ** 2) - wtot) - &
          By ** 2*(eta ** 2 + cs2*h*rho*(vx ** 2 + vz ** 2) - wtot) - &
          W ** 2*(eta ** 2 - wtot)*(eta ** 2 + &
          cs2*h*rho*(vx ** 2 + vy ** 2 + vz ** 2) - wtot) + &
          By*vy*(eta ** 3 + 2*Bz*cs2*h*rho*vz - eta*wtot)))

     a(0)=0.0e0_dp

     a(1)=0.0e0_dp
  end if

  !write(*,*) 'a=',a !!!!!!!!!!!!!!!!

  call quartic(a,rtr) !Analytic solution of quartic equation

  !write(*,*) 'rtr=',rtr !!!!!!!!!!!!
  !write(*,*)!!!!!!!!!!!!!!!!!!!!!!!!


  !Sort the roots in ascending order
  do i=1,3
     do j=i+1,4
        if(rtr(i) > rtr(j)) then
           temp = rtr(j) 
           rtr(j) = rtr(i) ! Swap the contents 
           rtr(i) = temp
        end if
     end do
  end do

  !Note: rtr(1) is the left going fast magnetosonic wave
  !      rtr(2) is the left going slow magnetosonic wave
  !      rtr(3) is the right going slow magnetosonic wave
  !      rtr(4) is the right going fast magnetosonic wave

  if (present(allwaves)) then
     allwaves=rtr
  end if

  if (switchLR=='LF') then
     Vs=rtr(1) !Left going fast rarefaction
     if (Vs>leftAlfven) then
        print *,'xi=',Vs
        print *,'Left Alfven=',leftAlfven
        stop 'XI: xi>leftAlfven'
     end if
  else if (switchLR=='LS') then
     Vs=rtr(2) !Left going slow rarefaction
     if ((Vs<leftAlfven).and.(abs(Vs-leftAlfven)>1.0e-10)) then
        print *,'xi=',Vs
        print *,'Left Alfven=',leftAlfven
        stop'XI: xi<leftAlfven'
     end if

     if (Vs>vx) then
        print *,'xi=',Vs
        print *,'vx=',vx
        stop 'XI: xi>vx'
     end if
  else if (switchLR=='RS') then
     Vs=rtr(3) !Right going slow rarefaction
     if(Vs>rightAlfven) then
        print *,'xi=',Vs
        print *,'Right Alfven=',rightAlfven
        stop 'XI: xi>rightAlfven'
     end if
     if (Vs<vx) then
        print *,'xi=',Vs
        print *,'vx=',vx
        stop 'XI: xi<vx'
     end if
  else if(switchLR=='RF') then
     Vs=rtr(4) !Right going fast rarefaction
     if ((Vs<rightAlfven).and.(abs(Vs-rightAlfven)>1.0e-10)) then
        print *,'xi=',Vs
        print *,'Right Alfven=',rightAlfven
        stop 'XI: xi<rightAlfven'
     end if
  else
     stop 'riemann.f90: xi: wrong switchLR'
  endif


end subroutine xi


subroutine eos_enthalpy(pgas,rho,gamma,enthalpy)
  use type
  use eos_param
  implicit none
  real(DP), intent(IN)::pgas,rho,gamma
  real(DP), intent(OUT)::enthalpy

  if (eos_ideal) then
     enthalpy = 1.0e0_dp+ gamma*pgas/((gamma-1.0e0_dp)*rho)
  elseif (eos_meliani) then
     enthalpy = pgas/(gamma-1.0e0_dp) + rho*sqrt((pgas/((gamma-1.0e0_dp)*rho))**2 + 1.0e0_dp) + pgas
     enthalpy = enthalpy/rho
  else
     write(*,*) 'eos_enthalpy: Error in the EOS'
     STOP
  end if


end subroutine eos_enthalpy


subroutine eos_cs2(pgas,rho,gamma,cs2)
  use type
  use eos_param
  implicit none
  real(DP), intent(IN)::pgas,rho,gamma
  real(DP), intent(OUT)::cs2
  
  real(DP)::w,emel,h

  if (eos_ideal) then

     w = rho + gamma*pgas/(gamma-1.0e0_dp)
     cs2 = gamma*pgas/w !c_s^2=\Gamma*Pgas/w

  elseif (eos_meliani) then

     h = pgas/(gamma-1.0e0_dp) + &
          rho*sqrt((pgas/((gamma-1.0e0_dp)*rho))**2 + 1.0e0_dp) + pgas
     h = h/rho
     emel = pgas/((gamma-1.0)*rho) + sqrt((pgas/((gamma-1)*rho))**2+1.0)
     cs2 = pgas/(rho*h)*(2.0*gamma + (gamma+1.0)*(emel**2-1.0))/(2.0*emel**2)
  
  else
     write(*,*) 'eos_cs2: Error in the EOS'
     STOP
  end if


end subroutine eos_cs2
