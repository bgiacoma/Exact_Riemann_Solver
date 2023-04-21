!$Id: quartic.f90,v 1.3 2007-08-13 09:56:44 bgiacoma Exp $

!!$ Copyright (C) 2005  B. Giacomazzo, L. Rezzolla

subroutine quartic(a,rtr)
  use type
  use interfaces,only:cubic
  implicit none
  
  real(DP),dimension(0:4),intent(IN)::a
  real(DP),dimension(4),intent(OUT)::rtr

  !This subroutine finds the analytic solutions of a quartic equation:
  !a4*x^4+a3*x^3+a2*x^2+a1*x+a0 == 0

  real(DP)::a0,a1,a2,a3,c0,c1,c2

  real(DP),dimension(3)::W
  real(DP),dimension(0:2)::c
  real(DP)::p,q
  real(DP)::alpha,beta,f,g,h
  
  logical(LGT)::check
  
  if (a(4)==0.0D0) stop 'a4=0: this is not a quartic equation'
  
  a0=a(0)/a(4)
  a1=a(1)/a(4)
  a2=a(2)/a(4)
  a3=a(3)/a(4)

  !x^4+a3*x^3+a2*x^2+a1*x+a0==0

  if(a0==0.0D0) then
     !It can be reduced to a cubic equation
     if (a1==0.0D0) then
        !It can be reduced to a square equation
        if (a2==0.0D0) then
           !It can be reduced to a linear equation
           rtr(1)=-a3 !0.0D0
           rtr(2)=-a3 !0.0D0
           rtr(3)=-a3 !0.0D0
           rtr(4)=-a3
           return
        end if
        rtr(1)=0.5D0*(-a3+sqrt(a3**2-4.0D0*a2))
        rtr(2)=0.5D0*(-a3-sqrt(a3**2-4.0D0*a2))
        rtr(3)=rtr(2)!0.0D0
        rtr(4)=rtr(1)!0.0D0
        return
     end if
     !solve x^4+a3*x^3+a2*x^2+a1*x==0
     c(0)=a1
     c(1)=a2
     c(2)=a3
     call cubic(c,W)
     rtr(1)=W(1)
     rtr(2)=W(2)
     rtr(3)=W(3)
     rtr(4)=rtr(1) !0.0D0
     return
     
  end if


  !First step: find the real solution of the following cubic equation
  !W^3+c2*W^2+c1*W+c0==0

  !these are the coefficients of the cubic equation
  f=a2-3.0D0*a3**2/8.0D0
  g=a1+a3**3/8.0D0-a3*a2/2.0D0
  h=a0-3.0D0*a3**4/256+a3**2*a2/16.0D0-a3*a1/4.0D0
  
  c(0)=-g**2/64.0D0                 !-a1**2-a0*a3**2+4.0D0*a0*a2
  c(1)=(f**2-4.0D0*h)/16.0D0        !a1*a3-4.0D0*a0
  c(2)=f/2.0D0                      !-a2

  check=.false.
  
  call cubic(c,W,check) 

  !write(*,*) 'cubic roots=',W,check !!!!!!!!!!!!!!!!!!!!!!!!!

  if(check) goto 1000 !only one real solution

  !p and q are the sqrt of non zero roots
  if (W(1)>0.0D0) then 
     p=sqrt(W(1))
     if (W(2)>0.0D0) then
        q=sqrt(W(2))
     else if (W(3)>0.0D0) then
        q=sqrt(W(3))
     else
        stop 'error in quartic: I cannot find two positive roots in the cubic resolvent'
        goto 1000
     end if
  else if (W(2)>0.0D0) then
     p=sqrt(W(2))
     if (W(3)>0.0D0) then
        q=sqrt(W(3))
     else
        stop 'error in quartic: I cannot find two positive roots in the cubic resolvent'
        goto 1000
     end if
  else
     stop 'error in quartic: I cannot find two positive roots in the cubic resolvent'
     goto 1000
  end if
  
  alpha=-g/(8.0D0*p*q)
  beta=a3/(4.0D0)
  rtr(1)=p+q+alpha-beta
  rtr(2)=p-q-alpha-beta
  rtr(3)=-p+q-alpha-beta
  rtr(4)=-p-q+alpha-beta
  return

  !alternative method to be used if only one root is real

1000 p=0.5D0*a3-sqrt(0.25D0*a3**2+W(1)-a2)
  q=0.5D0*W(1)-sqrt((0.5D0*W(1))**2-a0)
  rtr(1)=0.5D0*(-p+sqrt(p**2-4.0D0*q))
  rtr(2)=0.5D0*(-p-sqrt(p**2-4.0D0*q))

  !write(*,*) 'radice=',0.25D0*a3**2+W(1)-a2 !!!!!!!!!!
  !write(*,*) 'radice2=',(0.5D0*W(1))**2-a0 !!!!!!!!!!!!
  !write(*,*) p,q !!!!!!!!!!!!!!!!!!!!!!!

  p=0.5D0*a3+sqrt(0.25D0*a3**2+W(1)-a2)
  q=0.5D0*W(1)+sqrt((0.5D0*W(1))**2-a0)
  rtr(3)=0.5D0*(-p+sqrt(p**2-4.0D0*q))
  rtr(4)=0.5D0*(-p-sqrt(p**2-4.0D0*q))

  !write(*,*) p,q !!!!!!!!!!!!!!!!!!!!!!!

end subroutine quartic


subroutine cubic(c,W,check)
  use type
  implicit none
  real(DP),dimension(0:2),intent(IN)::c
  real(DP),dimension(1:3),intent(OUT)::W
  logical(LGT),intent(OUT),optional::check
  !This subroutine find the solution to a cubic edquation
  !check is true if only one real solution exists
  complex(DPC)::q,r,s1,s2
  complex(DPC),parameter::I=(0.0D0,1.0D0)
  

  if (c(0)==0.0D0) then
     if (c(1)==0.0D0) then
        !It's a linear equation
        W(1)=-c(2)
        W(2)=0.0D0
        W(3)=0.0D0
        return
     end if
     !It's a quadratic equation
     if(present(check)) then
        check=.false.
        if ((c(2)**2-4.0D0*c(1))<0.0D0) check=.true. !only one real solution
     end if
     W(1)=0.0D0
     W(2)=0.5D0*(-c(2)+sqrt(c(2)**2-4.0D0*c(1)))
     W(3)=0.5D0*(-c(2)-sqrt(c(2)**2-4.0D0*c(1)))
     return
  else

     !Solve the cubic equation (see Abramowitz and Stegun
     !"Handbook of Mathematical Functions")
     q=1.0D0/3.0D0*c(1)-1.0D0/9.0D0*c(2)**2
     r=1.0D0/6.0D0*(c(1)*c(2)-3.0D0*c(0))-1.0D0/27.0D0*c(2)**3
     
     !Trick to avoid a compiler problem
     !Compiler has some problem to do s**1/3 if s is a negative real number
     s1=sqrt(q**3+r**2)
     s1=(r+s1)
     if ((real(s1)<0.0D0).and.(aimag(s1)==0.0D0)) then
        s1=-((-s1)**(1.0D0/3.0D0))
     else
        s1=(s1)**(1.0D0/3.0D0)
     end if
     s2=sqrt(q**3+r**2)
     s2=(r-s2)
     if ((real(s2)<0.0D0).and.(aimag(s2)==0.0D0)) then
        s2=-((-s2)**(1.0D0/3.0D0))
     else
        s2=(s2)**(1.0D0/3.0D0)
     end if

     !write(*,*) '(q**3+r**2)', (q**3+r**2)!!!!!!!!!!!!!!!!!!!

     if ((real((q**3+r**2))<=0.0D0).or.(abs(real((q**3+r**2)))<=1.0e-25)) then
        !ALL the solutions are real
        if(present(check)) check=.false.
        W(1)=real((s1+s2)-c(2)/3.0D0)
        W(2)=real(-0.5D0*(s1+s2)-c(2)/3.0D0+0.5D0*I*sqrt(3.0D0)*(s1-s2))
        W(3)=real(-0.5D0*(s1+s2)-c(2)/3.0D0-0.5D0*I*sqrt(3.0D0)*(s1-s2))
     else
        !ONLY one solution is real
        if(present(check)) check=.true.
        W(1)=real((s1+s2)-c(2)/3.0D0)
        W(2)=W(1)!0.0/0.0!NaN
        W(3)=W(1)!0.0/0.0!NaN
     end if
  end if
     
end subroutine cubic
