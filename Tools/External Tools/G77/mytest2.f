      program mytest
c=================================================================
c  Solves the linear system defined as:
c
c    | 1    2    3    4....    N    | |X1|   |1.0|
c    | N    1    2    3....    N-1  | |X2|   |1.0|
c    | N-1  N    1    2....    N-2  |*|X3| = |1.0|
c    |    ................          | |. |   |1.0|
c    | 2    3    4    5....    1    | |XN|   |1.0|
c
c  which has the known exact solution: Xi = 2 /(n*(n+1), i=1,2,..N
c
c  NOTE: This version uses the Linpack routines DGEFA & DGESL
c        to perform the LU decomposition and solution of the
c        linear system.  These depend on some BLAS routines.
c	 All the necessary routines are contained in the SLATEC 
c        library, thus the user has to retrieve the library
c	 (see http://www.geocities.com/Athens/Olympus/5564/)
c	 and install it beforehand.
c
c        To compile:
c
c         g77 -O2 -mpentium mytest2.f -lslatec -o mytest2.exe
c                                      
c        or, alternatively (with explicit supply of the source code 
c 	 of the relevant LINPACK and BLAS routines:
c
c     	  g77 -O2 -mpentium mytest2.f dgefa.f dgesl.f daxpy.f dscal.f
c                         ddot.f idamax.f -o mytest2.exe
c
c        It's suggested that the user compares the solution speed
c        of mytest2.exe to the speed of mytest1.exe for, say, 1000 
c	 equations.
c
c=================================================================
      implicit none
      integer nmax
      parameter(nmax = 2000)		! max allowable number of equations
      real*8 aa(nmax,nmax),bb(nmax)
      character*1 affirm
      integer nn
      real*4 t, secnds	! secnds is a G77 intrinsic timer function


  5   write(*,10) 'Enter number of equations (up to ', nmax,'): '
 10   format(1x,a,i5,a,$)
      read(*,*) nn
      if ((nn.gt.nmax) .or. (nn.lt.1)) go to 5

      call prepare(aa, bb, nn, nmax)
C      call test_print(aa, bb, nn, nmax)  ! only for small test arrays
      t = secnds(0.0)
      call solve(aa, bb, nn, nmax)
      t = secnds(0.0) -t
      call results(bb, nn)
      write(*,15) 'Solution of the linear system of',nn,
     &           ' equations took ',t,' seconds'
 15   format(1x,a,i4,a,f7.2,a)
      if (affirm('Want to run again').eq.'Y') goto 5
      stop
      end

c--------------------------------------------
      character*1 function affirm(s)
c--------------------------------------------
      implicit none
      character *(*) s
      character*1    c

  5   write(*,10) s,' (y/n) ? : '
 10   format(1x,a,a,$)
      read(*,15) c
 15   format(a1)
      if (c.ne.'Y'.and.c.ne.'y'.and.c.ne.'N'.and.c.ne.'n') goto 5
      if (c.eq.'y') c='Y'
      affirm=c
      return
      end

c---------------------------------------
	subroutine prepare(a, b, n, nmax)
c---------------------------------------
      implicit none
      integer n, nmax
      real*8 a(nmax,*),b(*)
      integer i, j

      do 10 j=1,n
          a(1,j)=dble(j)*1.0D0
 10       b(j)=1.0D0
      do 20 i=2,n
         a(i,1)=a(i-1,n)
         do 20 j=2,n
 20         a(i,j)=a(i-1,j-1)
      return
      end


c-----------------------------------------
	subroutine test_print(a, b, n, nmax)
c-----------------------------------------
      implicit none
      integer n, nmax
      real*8 a(nmax,*),b(*)	
      integer i, j

      do 10 i=1,n
 10      print *,(a(i,j), j=1,n), b(i)

      return
      end

c-------------------------------------------------------
      subroutine solve(a, b, n, nmax)
c-------------------------------------------------------
c This version uses the LINPACK routines DGEFA, DGESL
c to effect the solution.  These routines are supplied
c by the SLATEC library, see the SLATEC documentation
c for details.
c-------------------------------------------------------
      implicit none
      integer  n, nmax
      real*8   a(nmax,*),b(*)
      integer  info
      integer  ipvt(nmax)	! dynamic creation of an auxiliary
				! array, needed for the solution!
				! This is not standard Fortran 77,
				! but it's supported by G77
	
      call dgefa( a, nmax, n, ipvt, info)
c	 on return, info = 0 implies everything's OK
c                   info = x implies a tiny pivot on equation x

      print *, ' on return from DGEFA, info= ', info

      call dgesl(a, nmax, n, ipvt, b, 0)
	
      return
      end


	subroutine results(b, n)
c---------------------------------------
      implicit none
      integer  n
      real*8   b(*)	
      real*8   diff, exact, temp
      integer  i, idiff

      diff = 0.0D0	
      idiff = 0
      exact=2.0D0/DBLE(n*(n+1))     
      do 10 i=1,n
c         print *,'x[',i,']=',b(i) ! activate either one only
c         write(*,20) i, b(i)	  ! for very small arrays.	
      temp = abs(b(i) - exact)
      if (temp .gt. diff) then
          diff = temp
          idiff = i
      endif	
 10   continue
c      print *,'exact=',exact
      write(*,30) exact
      write(*,40) 100.0D0*diff/exact, idiff
 20   format(1x,'x(',i5,')=',f16.10)
 30   format(1x,'exact =',f16.10)
 40   format(1x,'max diff was ',f20.12,'% of the exact ', 
     &        ' in eqn. No.', i5)
      return
      end





