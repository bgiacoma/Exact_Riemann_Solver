!$Id: Interfaces.f90,v 1.3 2007-08-13 09:56:44 bgiacoma Exp $

!!$ Copyright (C) 2005  B. Giacomazzo, L. Rezzolla

MODULE type
  INTEGER, PARAMETER :: I4B = SELECTED_INT_KIND(9)
  INTEGER, PARAMETER :: DP = KIND(1.0D0)
  INTEGER, PARAMETER :: DPC = KIND(1.0D0)!,1.0D0)
  INTEGER, PARAMETER :: LGT = KIND(.true.)
END MODULE type

module eos_param
  implicit none
  logical::eos_ideal, eos_meliani
end module eos_param


module global
  use type
  implicit none
  real(DP)::Bx,gamma
  logical::degen,verbose,veryverbose
  integer(I4B)::niter,initial_data
end module global

module accmod
  use type
  implicit none
  real(DP):: accuracy
end module accmod

module odeswitch
  use type
  implicit none
  character(len=2)::switch  
end module odeswitch

module wavecheck
  use type
  implicit none
  logical::wave_error,shock
  !real(DP)::rarefaction_pmin
end module wavecheck

module output_grid
  use type
  implicit none
  real(DP)::x1,x2,t
  integer(I4B)::nx
end module output_grid

module alfven_wave
  use type
  implicit none
  character(len=1)::alfvenwave
  real(DP),dimension(5)::leftalfven,rightalfven
end module alfven_wave


module interfaces
  
  interface
     subroutine initialdata(init,left,right)
       use type
       implicit  none
       integer(I4B),intent(IN)::init
       real(DP),dimension(7),intent(OUT)::left,right
     end subroutine initialdata
  end interface

  interface
     subroutine postshock(Vs,Byb,Bzb,vxb,ahead,behind)
       use type
       implicit none
       real(DP),intent(IN)::Vs,Byb,Bzb,vxb
       real(DP),intent(IN),dimension(7)::ahead
       real(DP),intent(OUT),dimension(7)::behind
     end subroutine postshock
  end interface

  interface
     subroutine velocity(unk1,unk2,ahead,switchLR,vxb,Vs,switchPB,behind)
       use type
       implicit none
       real(DP),intent(IN)::unk1,unk2
       real(DP),dimension(7),intent(IN)::ahead
       character(len=2),intent(IN)::switchLR
       real(DP),intent(OUT)::vxb,Vs
       character(len=1),intent(IN)::switchPB
       real(DP),dimension(7),intent(OUT),OPTIONAL::behind
     end subroutine velocity
  end interface

  interface
     subroutine velocity_df(unk1,unk2,ahead,vx,switchLR,f,df,Vs)
       use type
       implicit none
       real(DP),intent(IN)::unk1,unk2,vx
       real(DP),dimension(7),intent(IN)::ahead
       character(len=2),intent(IN)::switchLR
       real(DP),intent(OUT)::f,df
       real(DP),intent(OUT),OPTIONAL::Vs
     end subroutine velocity_df
  end interface

  interface
     subroutine velocity_eqn(unk1,unk2,ahead,vx,switchLR,f,Vs,stateb)
       use type
       implicit none
       real(DP),intent(IN)::unk1,unk2
       real(DP),dimension(7),intent(IN)::ahead
       real(DP),intent(IN)::vx
       character(len=2),intent(IN)::switchLR
       real(DP),intent(OUT)::f
       real(DP),intent(OUT)::Vs
       real(DP),dimension(7),OPTIONAL,intent(OUT)::stateb
     end subroutine velocity_eqn
  end interface

  interface
     subroutine ContactVelocity(unk1,unk2,vxb,Vs,ahead,vx,switchLR,switchPB,solution)
       use type
       implicit none
       real(DP),intent(IN)::unk1,unk2,vxb
       real(DP),intent(OUT)::Vs
       real(DP),dimension(7),intent(IN)::ahead
       real(DP),intent(OUT)::vx
       character(len=2),intent(IN)::switchLR
       character(len=1),intent(IN)::switchPB
       real(DP),dimension(7),intent(OUT),OPTIONAL::solution
     end subroutine ContactVelocity
  end interface

  interface
     subroutine shockfunc(Vs,unk1,unk2,pb,ahead,f,switchLR,switchPB,errcheck) 
       USE type
       IMPLICIT NONE 
       REAL(DP),INTENT(IN) :: Vs,unk1,unk2,pb
       real(DP),dimension(7),intent(IN) ::ahead
       REAL(DP),intent(OUT):: f
       character(len=2),intent(IN)::switchLR
       character(len=1),intent(IN)::switchPB
       logical,intent(OUT)::errcheck
     end subroutine shockfunc
  end interface

  interface
     SUBROUTINE zeroeqn(vxb,unk1,unk2,Vs,ahead,f,df,switchLR,switchPB,errcheck) 
       USE type
       IMPLICIT NONE 
       REAL(DP),INTENT(IN) :: vxb,unk1,unk2,Vs
       real(DP),dimension(7),intent(IN)::ahead
       REAL(DP),INTENT(OUT) :: f 
       REAL(DP),INTENT(OUT) :: df 
       character(len=2),intent(IN)::switchLR
       character(len=1),intent(IN)::switchPB
       logical,intent(OUT)::errcheck
     end SUBROUTINE zeroeqn
  end interface

  interface
     subroutine contact(left,right,pb,VsL,VsR,solution,degen_case)
       use type
       implicit none
       real(DP),dimension(7),intent(IN)::left,right
       real(DP),intent(OUT)::pb
       real(DP),intent(OUT)::VsL,VsR
       real(DP),dimension(7,2),intent(OUT)::solution
       integer(I4B),intent(IN)::degen_case
     end subroutine contact
  end interface

  interface
     subroutine fullcontact(left,right,unk,VsLv3,VsRv3,fullsolution)
       use type
       implicit none
       real(DP),dimension(7),intent(IN)::left,right
       real(DP),dimension(4),intent(INOUT)::unk
       real(DP),dimension(3),intent(OUT)::VsLv3,VsRv3
       real(DP),dimension(7,6),intent(OUT)::fullsolution
     end subroutine fullcontact
  end interface

  interface
     subroutine funcv(left,right,pb,VsL,VsR,f,left_wave,right_wave,solution)
       use type
       implicit none
       real(DP),dimension(7),intent(IN)::left,right
       real(DP),intent(IN)::pb
       real(DP),intent(OUT)::VsL,VsR
       real(DP),intent(OUT)::f
       character(len=2),intent(IN)::left_wave,right_wave
       real(DP),dimension(7,2),intent(OUT),OPTIONAL::solution
     end subroutine funcv
  end interface

  interface
     subroutine funcd(left,right,pb,VsL,VsR,f,df,left_wave,right_wave)
       use type
       implicit none
       real(DP),dimension(7),intent(IN)::left,right
       real(DP),intent(IN)::pb
       real(DP),intent(OUT)::VsL,VsR
       real(DP),intent(OUT)::f,df
       character(len=2),intent(IN)::left_wave,right_wave
     end subroutine funcd
  end interface

  interface
     subroutine fullfuncd(left,right,unk,VsLv3,VsRv3,fvec,fjac,fullsolution)
       use type
       implicit none
       real(DP),dimension(7),intent(IN)::left,right
       real(DP),dimension(4),intent(IN)::unk
       real(DP),dimension(3),intent(OUT)::VsLv3,VsRv3
       real(DP),dimension(4),intent(OUT)::fvec
       real(DP),dimension(4,4),intent(OUT)::fjac
       real(DP),dimension(7,6),intent(OUT),OPTIONAL::fullsolution
     end subroutine fullfuncd
  end interface

  interface
     subroutine fullfuncv(left,right,unk,VsLv3,VsRv3,fvec,fullsolution)
       use type
       implicit none
       real(DP),dimension(7),intent(IN)::left,right
       real(DP),dimension(4),intent(IN)::unk
       real(DP),dimension(3),intent(OUT)::VsLv3,VsRv3
       real(DP),dimension(4),intent(OUT)::fvec
       real(DP),dimension(7,6),intent(OUT),OPTIONAL::fullsolution
     end subroutine fullfuncv
  end interface

  interface
     subroutine output(x1,x2,t,nx,left,right,VsLv3,VsRv3,unk,fullsolution)
       use type
       implicit none
       real(DP),intent(IN)::x1,x2,t
       real(DP),dimension(4),intent(IN)::unk
       real(DP),dimension(3),intent(IN)::VsLv3,VsRv3
       integer(I4B),intent(IN)::nx
       real(DP),dimension(7),intent(IN)::left,right
       real(DP),dimension(7,6),intent(IN)::fullsolution
     end subroutine output
  end interface

  interface
     subroutine Rarefaction(pb,Vs,init,solution,switchLR)
       use type
       implicit none
       real(DP),intent(IN)::pb
       real(DP),intent(OUT)::Vs
       real(DP),dimension(7),intent(IN)::init
       real(DP),dimension(7),intent(OUT)::solution
       character(len=2),intent(IN)::switchLR
     end subroutine Rarefaction
  end interface

  interface
     subroutine RHS(x,y,f,switchLR)
       use type
       implicit none
       real(DP),intent(IN)::x
       real(DP),dimension(7),intent(IN)::y
       real(DP),dimension(7),intent(OUT)::f
       character(len=2),intent(IN)::switchLR
     end subroutine RHS
  end interface

  interface
     subroutine xi(state,switchLR,Vs,allwaves)
       use type
       implicit none
       real(DP),dimension(7),intent(IN)::state
       character(len=2),intent(IN)::switchLR
       real(DP),intent(OUT)::Vs
       real(DP),dimension(4),intent(OUT),OPTIONAL::allwaves
     end subroutine xi
  end interface

  interface
     subroutine quartic(a,rtr)
       use type
       implicit none
       real(DP),dimension(0:4),intent(IN)::a
       real(DP),dimension(4),intent(OUT)::rtr
     end subroutine quartic
  end interface

  interface
     subroutine cubic(c,W,check)
       use type
       implicit none
       real(DP),dimension(0:2),intent(IN)::c
       real(DP),dimension(1:3),intent(OUT)::W
       logical(LGT),intent(OUT),optional::check
     end subroutine cubic
  end interface

  interface
     subroutine shockvelocity(ahead,normB,wavetype,W)
       use type
       implicit none
       real(DP),dimension(7),intent(IN)::ahead
       real(DP),intent(IN)::normB
       character(len=2),intent(IN)::wavetype
       real(DP),intent(OUT)::W
     end subroutine shockvelocity
  end interface

  interface
     subroutine postshock_pressure(Vs,pb,ahead,behind,alfcheck)
       use type
       implicit none
       real(DP),intent(IN)::Vs,pb
       real(DP),intent(IN),dimension(7)::ahead !value of primitives ahead the shock
       real(DP),intent(OUT),dimension(7)::behind
       logical,intent(IN),OPTIONAL::alfcheck
     end subroutine postshock_pressure
  end interface

  interface
     subroutine Rarefaction_pressure(pb,Vs,init,solution,switchLR)
       use type
       implicit none
       real(DP),intent(IN)::pb
       real(DP),intent(OUT)::Vs
       real(DP),dimension(7),intent(IN)::init
       real(DP),dimension(7),intent(OUT)::solution
       character(len=2),intent(IN)::switchLR
     end subroutine Rarefaction_pressure
  end interface

  interface
     subroutine RHS_pressure(x,y,f,switchLR)
       use type
       implicit none
       real(DP),intent(IN)::x
       real(DP),dimension(7),intent(IN)::y
       real(DP),dimension(7),intent(OUT)::f
       character(len=2),intent(IN)::switchLR
     end subroutine RHS_pressure
  end interface

  interface
     subroutine alfven(Vs,ahead,behind)
       use type
       implicit none
       real(DP),intent(IN)::Vs
       real(DP),dimension(7),intent(IN)::ahead
       real(DP),dimension(7),intent(INOUT)::behind
     end subroutine alfven
  end interface

  interface
     subroutine alfven_df(ahead,Vs,behind,f,df)
       use type
       implicit none
       real(DP),intent(IN)::Vs
       real(DP),intent(IN),dimension(7)::ahead,behind
       real(DP),intent(OUT),dimension(5)::f
       real(DP),intent(OUT),dimension(5,5)::df
     end subroutine alfven_df
  end interface

  interface
     subroutine alfven_func(ahead,Vs,behind,f)
       use type
       implicit none
       real(DP),intent(IN)::Vs
       real(DP),dimension(7),intent(IN)::ahead,behind
       real(DP),intent(OUT),dimension(5)::f
     end subroutine alfven_func
  end interface

  interface
     subroutine Rarefaction_psi(normB,Vs,init,solution,switchLR)
       use type
       implicit none
       real(DP),intent(IN)::normB
       real(DP),intent(OUT)::Vs
       real(DP),dimension(7),intent(IN)::init
       real(DP),dimension(7),intent(OUT)::solution
       character(len=2),intent(IN)::switchLR
     end subroutine Rarefaction_psi
  end interface

  interface
     subroutine RHS_psi(x,y,f,switchLR)
       use type
       implicit none
       real(DP),intent(IN)::x
       real(DP),dimension(7),intent(IN)::y
       real(DP),dimension(7),intent(OUT)::f
       character(len=2),intent(IN)::switchLR
     end subroutine RHS_psi
  end interface

  interface
     subroutine eos_enthalpy(pgas,rho,gamma,enthalpy)
       use type
       implicit none
       real(DP), intent(IN)::pgas,rho,gamma
       real(DP), intent(OUT)::enthalpy
     end subroutine eos_enthalpy
  end interface

  interface
     subroutine eos_cs2(pgas,rho,gamma,cs2)
       use type
       implicit none
       real(DP), intent(IN)::pgas,rho,gamma
       real(DP), intent(OUT)::cs2
     end subroutine eos_cs2
  end interface

  !interfaces to NR routines   
  INTERFACE
     SUBROUTINE lubksb(a,indx,b)
       USE type
       REAL(DP), DIMENSION(:,:), INTENT(IN) :: a
       INTEGER(I4B), DIMENSION(:), INTENT(IN) :: indx
       REAL(DP), DIMENSION(:), INTENT(INOUT) :: b
     END SUBROUTINE lubksb
  END INTERFACE
  INTERFACE
     SUBROUTINE ludcmp(a,indx,d)
       USE type
       REAL(DP), DIMENSION(:,:), INTENT(INOUT) :: a
       INTEGER(I4B), DIMENSION(:), INTENT(OUT) :: indx
       REAL(DP), INTENT(OUT) :: d
     END SUBROUTINE ludcmp
  END INTERFACE
  

end module interfaces
