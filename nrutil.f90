!(C) 1986-1996 Numerical Recipes Software

MODULE nrutil
  USE type
  IMPLICIT NONE
  INTERFACE swap
     MODULE PROCEDURE swap_i,swap_r,swap_rv, &
          masked_swap_rs,masked_swap_rv,masked_swap_rm
  END INTERFACE

  INTERFACE imaxloc
     MODULE PROCEDURE imaxloc_r,imaxloc_i
  END INTERFACE

  INTERFACE assert_eq
     MODULE PROCEDURE assert_eq2,assert_eq3,assert_eq4,assert_eqn
  END INTERFACE

  INTERFACE outerprod
     MODULE PROCEDURE outerprod_r
  END INTERFACE
  
contains
  
  SUBROUTINE swap_i(a,b)
    INTEGER(I4B), INTENT(INOUT) :: a,b
    INTEGER(I4B) :: dum
    dum=a
    a=b
    b=dum
  END SUBROUTINE swap_i

  SUBROUTINE swap_r(a,b)
    REAL(DP), INTENT(INOUT) :: a,b
    REAL(DP) :: dum
    dum=a
    a=b
    b=dum
  END SUBROUTINE swap_r

  SUBROUTINE swap_rv(a,b)
    REAL(DP), DIMENSION(:), INTENT(INOUT) :: a,b
    REAL(DP), DIMENSION(SIZE(a)) :: dum
    dum=a
    a=b
    b=dum
  END SUBROUTINE swap_rv

  SUBROUTINE masked_swap_rs(a,b,mask)
    REAL(DP), INTENT(INOUT) :: a,b
    LOGICAL(LGT), INTENT(IN) :: mask
    REAL(DP) :: swp
    if (mask) then
       swp=a
       a=b
       b=swp
    end if
  END SUBROUTINE masked_swap_rs

  SUBROUTINE masked_swap_rv(a,b,mask)
    REAL(DP), DIMENSION(:), INTENT(INOUT) :: a,b
    LOGICAL(LGT), DIMENSION(:), INTENT(IN) :: mask
    REAL(DP), DIMENSION(size(a)) :: swp
    where (mask)
       swp=a
       a=b
       b=swp
    end where
  END SUBROUTINE masked_swap_rv

  SUBROUTINE masked_swap_rm(a,b,mask)
    REAL(DP), DIMENSION(:,:), INTENT(INOUT) :: a,b
    LOGICAL(LGT), DIMENSION(:,:), INTENT(IN) :: mask
    REAL(DP), DIMENSION(size(a,1),size(a,2)) :: swp
    where (mask)
       swp=a
       a=b
       b=swp
    end where
  END SUBROUTINE masked_swap_rm

  FUNCTION imaxloc_r(arr)
    REAL(DP), DIMENSION(:), INTENT(IN) :: arr
    INTEGER(I4B) :: imaxloc_r
    INTEGER(I4B), DIMENSION(1) :: imax
    imax=maxloc(arr(:))
    imaxloc_r=imax(1)
  END FUNCTION imaxloc_r

  FUNCTION imaxloc_i(iarr)
    INTEGER(I4B), DIMENSION(:), INTENT(IN) :: iarr
    INTEGER(I4B), DIMENSION(1) :: imax
    INTEGER(I4B) :: imaxloc_i
    imax=maxloc(iarr(:))
    imaxloc_i=imax(1)
  END FUNCTION imaxloc_i

  FUNCTION assert_eq2(n1,n2,string)
    CHARACTER(LEN=*), INTENT(IN) :: string
    INTEGER, INTENT(IN) :: n1,n2
    INTEGER :: assert_eq2
    if (n1 == n2) then
       assert_eq2=n1
    else
       write (*,*) 'nrerror: an assert_eq failed with this tag:', &
            string
       STOP 'program terminated by assert_eq2'
    end if
  END FUNCTION assert_eq2

  FUNCTION assert_eq3(n1,n2,n3,string)
    CHARACTER(LEN=*), INTENT(IN) :: string
    INTEGER, INTENT(IN) :: n1,n2,n3
    INTEGER :: assert_eq3
    if (n1 == n2 .and. n2 == n3) then
       assert_eq3=n1
    else
       write (*,*) 'nrerror: an assert_eq failed with this tag:', &
            string
       STOP 'program terminated by assert_eq3'
    end if
  END FUNCTION assert_eq3

  FUNCTION assert_eq4(n1,n2,n3,n4,string)
    CHARACTER(LEN=*), INTENT(IN) :: string
    INTEGER, INTENT(IN) :: n1,n2,n3,n4
    INTEGER :: assert_eq4
    if (n1 == n2 .and. n2 == n3 .and. n3 == n4) then
       assert_eq4=n1
    else
       write (*,*) 'nrerror: an assert_eq failed with this tag:', &
            string
       STOP 'program terminated by assert_eq4'
    end if
  END FUNCTION assert_eq4

  FUNCTION assert_eqn(nn,string)
    CHARACTER(LEN=*), INTENT(IN) :: string
    INTEGER, DIMENSION(:), INTENT(IN) :: nn
    INTEGER :: assert_eqn
    if (all(nn(2:) == nn(1))) then
       assert_eqn=nn(1)
    else
       write (*,*) 'nrerror: an assert_eq failed with this tag:', &
            string
       STOP 'program terminated by assert_eqn'
    end if
  END FUNCTION assert_eqn

  SUBROUTINE nrerror(string)
    CHARACTER(LEN=*), INTENT(IN) :: string
    write (*,*) 'nrerror: ',string
    STOP 'program terminated by nrerror'
  END SUBROUTINE nrerror

  FUNCTION outerprod_r(a,b)
    REAL(DP), DIMENSION(:), INTENT(IN) :: a,b
    REAL(DP), DIMENSION(size(a),size(b)) :: outerprod_r
    outerprod_r = spread(a,dim=2,ncopies=size(b)) * &
         spread(b,dim=1,ncopies=size(a))
  END FUNCTION outerprod_r

END MODULE nrutil
