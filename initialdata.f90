!$Id: initialdata.f90,v 1.2 2007-07-16 10:42:34 bgiacoma Exp $

!!$ Copyright (C) 2005  B. Giacomazzo, L. Rezzolla

subroutine initialdata(init,left,right)
  use type
  use global !This module contains gamma and Bx
  implicit  none
  integer(I4B),intent(IN)::init
  real(DP),dimension(7),intent(OUT)::left,right

  !Initial conditions for the Riemann problem

  real(DP),dimension(16)::value
  character(len=22)::dump
  integer(I4B)::i

  select case(init)
  case (0) !User defined

     !Read from file RInput.txt
     open(UNIT=200,FILE='RInput.txt',STATUS='OLD',ACTION='READ')
     do i = 1,16
        read (200,280) dump,value(i)
        write(*,*)dump,value(i)
     enddo
     write(*,*)

280  format(a22,f40.20)
     close (200)

     !B^x
     Bx=value(1)
     !EoS, ideal gas: P=\rho\eps*(\Gamma-1.0)
     gamma=value(2)

     !Left state
     left(1)=value(3) !rho
     left(2)=value(4) !gas pressure
     left(3)=value(5) !vx
     left(4)=value(6) !vy
     left(5)=value(7) !vz
     left(6)=value(8) !By
     left(7)=value(9) !Bz

     !Right state
     right(1)=value(10)  !rho
     right(2)=value(11)  !gas pressure
     right(3)=value(12)  !vx
     right(4)=value(13)  !vy
     right(5)=value(14)  !vz
     right(6)=value(15)  !By
     right(7)=value(16)  !Bz


  case (1) !Marti-Muller figure 7
     !B^x
     Bx=0.0D0
     !EoS, ideal gas: P=\rho\eps*(\Gamma-1.0)
     gamma=4.0D0/3.0D0

     !Left state
     left(1)=1.0D0 !rho
     left(2)=1.0D0 !gas pressure
     left(3)=0.9D0 !vx
     left(4)=0.0D0 !vy
     left(5)=0.0D0 !vz
     left(6)=0.0D0 !By
     left(7)=0.0D0 !Bz

     !Right state
     right(1)=1.0D0  !rho
     right(2)=10.0D0 !gas pressure
     right(3)=0.0D0  !vx
     right(4)=0.0D0  !vy
     right(5)=0.0D0  !vz
     right(6)=0.0D0  !By
     right(7)=0.0D0  !Bz


  case (2) !Marti-Muller figure 6
     !B^x
     Bx=0.0D0
     !EoS, ideal gas: P=\rho\eps*(\Gamma-1.0)
     gamma=5.0D0/3.0D0

     !Left state
     left(1)=1.0D0  !rho
     left(2)=10.0D0 !gas pressure
     left(3)=-0.6D0 !vx
     left(4)=0.0D0  !vy
     left(5)=0.0D0  !vz
     left(6)=0.0D0  !By
     left(7)=0.0D0  !Bz

     !Right state
     right(1)=10.0D0  !rho
     right(2)=20.0D0  !gas pressure
     right(3)=0.5D0   !vx
     right(4)=0.0D0   !vy
     right(5)=0.0D0   !vz
     right(6)=0.0D0   !By
     right(7)=0.0D0   !Bz


  case (3) !Marti-Muller figure 5
     !B^x
     Bx=0.0D0
     !EoS, ideal gas: P=\rho\eps*(\Gamma-1.0)
     gamma=5.0D0/3.0D0

     !Left state
     left(1)=1.0D0 !rho
     left(2)=1.0D3 !gas pressure
     left(3)=0.0D0 !vx
     left(4)=0.0D0 !vy
     left(5)=0.0D0 !vz
     left(6)=0.0D0 !By
     left(7)=0.0D0 !Bz

     !Right state
     right(1)=1.0D0   !rho
     right(2)=1.0D-2  !gas pressure
     right(3)=0.0D0   !vx
     right(4)=0.0D0   !vy
     right(5)=0.0D0   !vz
     right(6)=0.0D0   !By
     right(7)=0.0D0   !Bz


  case (4) !Generic Shock-Tube Test
     !B^x
     Bx=0.0D0
     !EoS, ideal gas: P=\rho\eps*(\Gamma-1.0)
     gamma=5.0D0/3.0D0

     !Left state
     left(1)=1.0D0  !rho
     left(2)=1.0D-2 !gas pressure
     left(3)=0.1D0  !vx
     left(4)=0.3D0  !vy
     left(5)=0.4D0  !vz
     left(6)=6.0D0  !By
     left(7)=2.0D0  !Bz

     !Right state
     right(1)=1.0D-2  !rho
     right(2)=5.0D3   !gas pressure
     right(3)=0.5D0   !vx
     right(4)=0.4D0   !vy
     right(5)=0.3D0   !vz
     right(6)=5.0D0   !By
     right(7)=20.0D0  !Bz


  case (5) !Komissarov: Shock Tube 2
     !B^x
     Bx=0.0D0
     !EoS, ideal gas: P=\rho\eps*(\Gamma-1.0)
     gamma=4.0D0/3.0D0

     !Left state
     left(1)=1.0D0  !rho
     left(2)=30.0D0 !gas pressure
     left(3)=0.0D0  !vx
     left(4)=0.0D0  !vy
     left(5)=0.0D0  !vz
     left(6)=20.0D0 !By
     left(7)=0.0D0  !Bz

     !Right state
     right(1)=0.1D0   !rho
     right(2)=1.0D0   !gas pressure
     right(3)=0.0D0   !vx
     right(4)=0.0D0   !vy
     right(5)=0.0D0   !vz
     right(6)=0.0D0   !By
     right(7)=0.0D0   !Bz

  case (6) !Komissarov: Shock Tube 1
     !B^x
     Bx=1.0D0
     !EoS, ideal gas: P=\rho\eps*(\Gamma-1.0)
     gamma=4.0D0/3.0D0

     !Left state
     left(1)=1.0D0  !rho
     left(2)=1.0D3  !gas pressure
     left(3)=0.0D0  !vx
     left(4)=0.0D0  !vy
     left(5)=0.0D0  !vz
     left(6)=0.0D0  !By
     left(7)=0.0D0  !Bz

     !Right state
     right(1)=0.1D0   !rho
     right(2)=1.0D0   !gas pressure
     right(3)=0.0D0   !vx
     right(4)=0.0D0   !vy
     right(5)=0.0D0   !vz
     right(6)=0.0D0   !By
     right(7)=0.0D0   !Bz

  case (7) !Komissarov: Collision
     !B^x
     Bx=10.0D0
     !EoS, ideal gas: P=\rho\eps*(\Gamma-1.0)
     gamma=4.0D0/3.0D0

     !Left state
     left(1)=1.0D0               !rho
     left(2)=1.0D0               !gas pressure
     left(3)=5.0D0/sqrt(26.0D0)  !vx
     left(4)=0.0D0               !vy
     left(5)=0.0D0               !vz
     left(6)=10.0D0              !By
     left(7)=0.0D0               !Bz

     !Right state
     right(1)=1.0D0              !rho
     right(2)=1.0D0              !gas pressure
     right(3)=-5.0D0/sqrt(26.0D0)!vx
     right(4)=0.0D0              !vy
     right(5)=0.0D0              !vz
     right(6)=-10.0D0            !By
     right(7)=0.0D0              !Bz

  case (8) !Balsara: figure1 (also known as relativistic Brio-Wu shock tube test)
     !B^x
     Bx=0.5D0
     !EoS, ideal gas: P=\rho\eps*(\Gamma-1.0)
     gamma=2.0D0 
     

     !Left state
     left(1)=1.0D0  !rho
     left(2)=1.0D0  !gas pressure
     left(3)=0.0D0  !vx
     left(4)=0.0D0  !vy
     left(5)=0.0D0  !vz
     left(6)=1.0D0  !By
     left(7)=0.0D0  !Bz

     !Right state
     right(1)=0.125D0 !rho
     right(2)=0.1D0   !gas pressure
     right(3)=0.0D0   !vx
     right(4)=0.0D0   !vy
     right(5)=0.0D0   !vz
     right(6)=-1.0D0  !By
     right(7)=0.0D0   !Bz

  case (9) !Balsara: figure2
     !B^x
     Bx=5.0D0
     !EoS, ideal gas: P=\rho\eps*(\Gamma-1.0)
     gamma=1.6666D0

     !Left state
     left(1)=1.0D0  !rho
     left(2)=30.0D0 !gas pressure
     left(3)=0.0D0  !vx
     left(4)=0.0D0  !vy
     left(5)=0.0D0  !vz
     left(6)=6.0D0  !By
     left(7)=6.0D0  !Bz

     !Right state
     right(1)=1.0D0 !rho
     right(2)=1.0D0 !gas pressure
     right(3)=0.0D0 !vx
     right(4)=0.0D0 !vy
     right(5)=0.0D0 !vz
     right(6)=0.7D0 !By
     right(7)=0.7D0 !Bz

  case (10) !Balsara: figure3
     !B^x
     Bx=10.0D0
     !EoS, ideal gas: P=\rho\eps*(\Gamma-1.0)
     gamma=1.6666D0

     !Left state
     left(1)=1.0D0  !rho
     left(2)=1.0D3  !gas pressure
     left(3)=0.0D0  !vx
     left(4)=0.0D0  !vy
     left(5)=0.0D0  !vz
     left(6)=7.0D0  !By
     left(7)=7.0D0  !Bz

     !Right state
     right(1)=1.0D0 !rho
     right(2)=0.1D0 !gas pressure
     right(3)=0.0D0 !vx
     right(4)=0.0D0 !vy
     right(5)=0.0D0 !vz
     right(6)=0.7D0 !By
     right(7)=0.7D0 !Bz

  case (11) !Balsara: figure 4
     !B^x
     Bx=10.0D0
     !EoS, ideal gas: P=\rho\eps*(\Gamma-1.0)
     gamma=1.6666D0

     !Left state
     left(1)=1.0D0  !rho
     left(2)=0.1D0  !gas pressure
     left(3)=0.999D0!vx
     left(4)=0.0D0  !vy
     left(5)=0.0D0  !vz
     left(6)=7.0D0  !By
     left(7)=7.0D0  !Bz

     !Right state
     right(1)=1.0D0   !rho
     right(2)=0.1D0   !gas pressure
     right(3)=-0.999D0!vx
     right(4)=0.0D0   !vy
     right(5)=0.0D0   !vz
     right(6)=-7.0D0  !By
     right(7)=-7.0D0  !Bz

  case (12) !Balsara: figure 5
     !B^x
     Bx=2.0D0
     !EoS, ideal gas: P=\rho\eps*(\Gamma-1.0)
     gamma=1.6666D0

     !Left state
     left(1)=1.08D0 !rho
     left(2)=0.95D0 !gas pressure
     left(3)=0.4D0  !vx
     left(4)=0.3D0  !vy
     left(5)=0.2D0  !vz
     left(6)=0.3D0  !By
     left(7)=0.3D0  !Bz

     !Right state
     right(1)=1.0D0  !rho
     right(2)=1.0D0  !gas pressure
     right(3)=-0.45D0!vx
     right(4)=-0.2D0 !vy
     right(5)=0.2D0  !vz
     right(6)=-0.7D0 !By
     right(7)=0.5D0  !Bz

  case (13) !generic Alfven test
     !B^x
     Bx=1.0e0_dp
     !EoS, ideal gas: P=\rho\eps*(\Gamma-1.0)
     gamma=5.0e0_dp/3.0e0_dp

     !Left state
     left(1)=1.0e0_dp  !rho
     left(2)=5.0e0_dp  !gas pressure
     left(3)=0.0e0_dp  !vx
     left(4)=0.3e0_dp  !vy
     left(5)=0.4e0_dp  !vz
     left(6)=6.0e0_dp  !By
     left(7)=2.0e0_dp  !Bz

     !Right state
     right(1)=0.9e0_dp   !rho
     right(2)=5.3e0_dp   !gas pressure
     right(3)=0.0e0_dp   !vx
     right(4)=0.0e0_dp   !vy
     right(5)=0.0e0_dp   !vz
     right(6)=5.0e0_dp   !By
     right(7)=2.0e0_dp   !Bz

  case default
     stop 'Wrong choice'

  end select

  !Computing the total pressure
  !Total pessure is used instead of gas pressure 
  !in shocks and rarefactions equations
  !P_{tot} = P_{gas} + 0.5 * (b^\mu b_\mu)

  left(2)=left(2)+0.5D0*((Bx**2+left(6)**2+left(7)**2)* &
       (1.0D0-left(3)**2-left(4)**2-left(5)**2)+ &
       (left(3)*Bx+left(4)*left(6)+left(5)*left(7))**2)

  right(2)=right(2)+0.5D0*((Bx**2+right(6)**2+right(7)**2)* &
       (1.0D0-right(3)**2-right(4)**2-right(5)**2)+ &
       (right(3)*Bx+right(4)*right(6)+right(5)*right(7))**2)

end subroutine initialdata
