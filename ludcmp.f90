SUBROUTINE ludcmp(a,indx,d)
  !Given an N×N input matrix a, 
  !this routine replaces it by the LU decomposition
  !of a rowwise permutation of itself.
  !On output, a is arranged as in equation (2.3.14);
  !indx is an output vector of length N that records the row permutation
  !effected by the partial pivoting;
  !d is output as ±1 depending on whether the number of row 
  !interchanges was even or odd, respectively.
  !This routine is used in combination with lubksb to solve linear equations 
  !or invert a matrix.

  !(C) 1986-1996 Numerical Recipes Software

  USE type; USE nrutil, ONLY : assert_eq,imaxloc,nrerror,outerprod,swap 
  IMPLICIT NONE 
  REAL(DP), DIMENSION(:,:), INTENT(INOUT) :: a 
  INTEGER(I4B), DIMENSION(:), INTENT(OUT) :: indx 
  REAL(DP), INTENT(OUT) :: d
  REAL(DP), DIMENSION(size(a,1)) :: vv 
  !vv stores the implicit scaling of each row. 
  REAL(DP), PARAMETER :: TINY=1.0e-20_dp !A small number. 
  INTEGER(I4B) :: j,n,imax 
  n=assert_eq(size(a,1),size(a,2),size(indx),'ludcmp')
  d=1.0 !No row interchanges yet. 
  vv=maxval(abs(a),dim=2) 
  !Loop over rows to get the implicit scaling information
  if (any(vv == 0.0)) call nrerror('singular matrix in ludcmp') !There is a row of zeros. 
  vv=1.0_dp/vv !Save the scaling.
  do j=1,n 
     imax=(j-1)+imaxloc(vv(j:n)*abs(a(j:n,j))) !Find the pivot row. 
     if (j /= imax) then !Do we need to interchange rows? 
        call swap(a(imax,:),a(j,:)) !Yes,do so... 
        d=-d !...and change the parity of d 
        vv(imax)=vv(j) !Also interchange the scale factor. 
     end if
     indx(j)=imax 
     if (a(j,j) == 0.0) a(j,j)=TINY 
     !If the pivot element is zero the matrix is singular 
     !(at least to the precision of the algorithm).
     !For some applications on singular matrices, it is desirable 
     !to substitute TINY for zero. 
     a(j+1:n,j)=a(j+1:n,j)/a(j,j) !Divide y the pivot element. 
     a(j+1:n,j+1:n)=a(j+1:n,j+1:n)-outerprod(a(j+1:n,j),a(j,j+1:n)) 
     !Reduce remaining submatrix. 
  end do
END SUBROUTINE ludcmp
