      program mytest1
c=================================================================
c  Solves the linear system of equations defined as:
c
c    | 1    2    3    4  ....    N    | |X1|   |1.0|
c    | N    1    2    3  ....    N-1  | |X2|   |1.0|
c    | N-1  N    1    2  ....    N-2  |*|X3| = |1.0|
c    |    ..................          | |..|   |1.0|
c    | 2    3    4    5  ....    1    | |XN|   |1.0|
c
c  which has the known exact solution: Xi = 2 /(n*(n+1)), i=1,2,..N
c
c  by in-situ decomposition of the coefficient array into upper
c  and lower triagular parts ( A = L * U, see subroutine "solve".)
c
c  To compile:  g77 -O2 -mpentium mytest1.f -o mytest1.exe
c=================================================================
      implicit none
      integer nmax, nn
      parameter (nmax=1000)       ! Max allowable number of equations
      real*8 aa(nmax,nmax),bb(nmax)
      character*1 affirm
      real*4 secnds               ! an intrinsic timer function in G77
      real*4 t

 10    write(*,20) 'Enter number of equations (up to ', nmax,'): '
 20   format(1x,a,i4,a,$)
      read(*,*) nn
      if ((nn.gt.nmax) .or. (nn.lt.1)) go to 10

      call prepare(aa,bb,nn)
c      call testprint(aa,bb,nn)     ! activate for very small arrays
      t = secnds(0.0)
      call solve(aa,bb,nn)
      t = secnds(0.0) -t
      call results(bb,nn)
      write(*,30) 'Solution of the linear system of ',nn,
     1           ' equations took ',t,' seconds'
 30   format(1x,a,i4,a,f7.2,a)
      if (affirm('Want to run again').eq.'Y') goto 10  ! repeat
	stop
	end

      character*1 function affirm(s)
c--------------------------------------------
      implicit none
      character *(*) s
      character*1    c

 10   write(*,20) s,' (y/n) ? : '
 20   format(1x,a,a,$)
      read(*,30) c
 30   format(a1)
      if (c.ne.'Y'.and.c.ne.'y'.and.c.ne.'N'.and.c.ne.'n') goto 10
      if (c.eq.'y') c='Y'
      affirm=c
      return
      end


      subroutine prepare(a,b,n)
c---------------------------------------
      implicit none
      integer n, j, i
      real*8 a(n,n), b(n)

      do 10 j=1,n
         a(1,j)=dble(j)*1.0D0
 10      b(j)=1.0D0
      do 20 i=2,n
         a(i,1)=a(i-1,n)
         do 20 j=2,n
 20         a(i,j)=a(i-1,j-1)
	return
	end


      subroutine testprint(a,b,n)
c-----------------------------------------
      implicit none
      integer n, i, j
      real*8 a(n,n), b(n)

      do 10 i=1,n
 10      print *,(a(i,j), j=1,n), b(i)

      return
      end

	
	subroutine solve(a,b,n)
c-----------------------------------------
      implicit none
      integer n, k, i, m, j
      real*8 a(n,n), b(n), sum

c  First decompose in-situ into LU
	
	do 50 k=2,n	
         do 30 i=1,k-1
            do 10 m=1,i-1
               a(k,i)=a(k,i)-a(k,m)*a(m,i)
               a(i,k)=a(i,k)-a(m,k)*a(i,m)
 10         continue
c           Here fix the "lower" part term (divide it by
c           the diagonal coefficient Uii)
            if (abs(a(i,i)).ge.1.e-10) go to 20 ! check for near zero
            write(*,*) 'Fatal ERROR - Tiny pivot in eq.',i
            stop
 20         a(k,i)=a(k,i)/a(i,i)
 30      continue
c        Now calculate the diagonal term of the "upper" part
         do 40 m=1,k-1
            a(k,k)=a(k,k)-a(m,k)*a(k,m)
 40      continue  ! next m
 50     continue  ! next k

c  Now solve for the lower triangular part

	do 70 i=2,n
         sum=0.0
         do 60 j=1,i-1
            sum=sum+a(i,j)*b(j)
 60      continue
         b(i)=b(i)-sum
 70     continue

c  Now solve for the upper triangular part

      do 100 i=n,1,-1
         sum=0.0
         do 80 j=i+1,n
            sum=sum+a(i,j)*b(j)
 80      continue

c             Check the diagonal coefficient prior to division

         if (abs(a(i,i)).ge.1.e-10) go to 90
         write(*,*) 'Fatal ERROR - Tiny pivot in eq.',i
         stop

 90      b(i)=(b(i)-sum)/a(i,i)
100   continue                             

      return
      end


      subroutine results(b,n)
c---------------------------------------
      implicit none
      integer n, i, idiff
      real*8 b(n), exact, diff, diffmax

      exact=2.0D0/dble(n*(n+1))
      diffmax = 0.0D0
      idiff = 1
      do 10 i=1,n
c        write(*,20) i, b(i) ! only for very small arrays
         diff = abs(b(i) - exact)
         if (diff .gt. diffmax) then
            diffmax = diff
            idiff = i
            endif
 10      continue
      write(*,30) 'exact =', exact
      write(*,40) 'Max. difference was ', diffmax*100.0D0/exact,
     1            '% of the exact in eq. No.',idiff
 20   format(1x,'x(',i5,')=',f16.10)
 30   format(1x,a,f16.10)
 40   format(1x,a,f20.12,a,i5)
      return
      end

