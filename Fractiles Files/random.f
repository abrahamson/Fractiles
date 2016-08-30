
c ----------------------------
      subroutine GetRandom1 ( iseed, n, wt, i1, iSave, n1, iflag )

      implicit none

      integer iseed, i1, isave, n1, n, i, iflag
      real ran1, wt(n1, n1), x
      
c     Get random number
      x = ran1( iseed )

      do i=1,n
        if ( x .le. wt(i1,i) ) then
          iSave = i
          return
        endif
      enddo
      
      write (*,*) ' Get Random Number 1'
      write (*,'( 2x,''Error - bad ran number or weights'')')
      write (*,*) ' Random Number        = ', x
      write (*,*) ' Fixed Parameter      = ', i1
      write (*,*) ' Number of Parameters = ', n
      do i=1,n
         write (*,*) wt(i1,i)
      enddo
      write (*,'( 2x,''iflag ='',i5)') iflag
      stop 99
      end
      
c ----------------------------
      subroutine GetRandom1b ( iseed, n, wt, i1, iSave, n1, n2 )

      implicit none

      integer iseed, i1, isave, n1, n, n2, i
      real ran1, wt(n2, n1), x

c     Get random number
      x = ran1( iseed )

      do i=1,n
        if ( x .le. wt(i1,i) ) then
          iSave = i
          return
        endif
      enddo
      
      write (*,*) ' Get Random Number 1b'
      write (*,'( 2x,''Error - bad ran number or weights'')')
      write (*,*) x
      write (*,*) (wt(i1,i),i=1,n)
      stop 99
      end
      
c ----------------------------
      subroutine GetRandom2 ( iseed, n, wt, i1, i2, iSave, n1, n2 )

      implicit none
      include 'fract.h'

      integer n1, n2, n, iseed, i, i1, i2, iSave
      real ran1, wt(n1, n2, MAXPARAM), x

c     Get random number
      x = ran1( iseed )
      do i=1,n
        if ( x .le. wt(i1,i2,i) ) then
          iSave = i
          return
        endif
      enddo
      
      write (*,*) ' Get Random Number 2'
      write (*,'( 2x,''Error - bad ran number 2 or weights'')')

      write (*,*) 'Random Number = ', x
      write (*,'(2x,''wts:'',10f10.4)') (wt(i1,i2,i),i=1,n)
      write (*,'( 3i5)') n,i1, i2

      stop 99
      end

c ----------------------------

      function Ran1 ( idum )

      implicit none

c     Random number generator, From numerical recipes
      integer idum, ia, im, iq, ir, ntab, ndiv
      real ran1, am, eps, rnmx
      parameter (ia=16807, im=2147483647, am=1./im,iq=127773,ir=2836,
     1      ntab=32,ndiv=1+(im-1)/ntab,eps=1.2e-7,rnmx=1.-eps)
      integer j, k, iv(ntab), iy
      save iv, iy
      data iv /ntab*0/, iy /0/
      
      if (idum .le. 0 .or. iy .eq. 0 ) then
        idum=max(-idum,1)
        do j=ntab+8,1,-1
          k=idum/iq
          idum=ia*(idum-k*iq)-ir*k
          if( idum .lt. 0) idum=idum+im
          if (j .le. ntab) iv(j)=idum
        enddo
        iy = iv(1)
      endif
      k = idum/iq
      idum=ia*(idum-k*iq) - ir*k
      if ( idum .lt. 0 ) idum = idum + im
      j = 1 + iy/ndiv
      iy = iv(j)
      iv (j) = idum
      ran1 = min(am*iy,rnmx)
      return
      end

c ----------------------------------------------------------------------

      subroutine CheckDim ( n, nMax, name )
      
      implicit none
     
      integer n, nMax
      character*80 name
      
      if ( n .gt. nMax ) then
        write (*,'( 2x,''Array Dimension Too Small'')')
        write (*,'( 2x,''Increase '',a20,'' to '',i5)') name, n
        stop 99
      endif
      return
      end

c --------------------------

      subroutine CheckWt ( x, n, fName, name )
      
      implicit none
      
      integer i, n
      real x(1), sum
      character*80 name, fName
      
      sum = 0.
      do i=1,n
        sum = sum + x(i)
      enddo
      if ( sum .ne. 1. ) then
        write (*,*) ' CheckWt Subroutine.'
        write (*,'( 2x,''Error -- Weights do not sum to unity'')')
        write (*,'( 2x,a80)') name
        write (*,'( 2x,a80)') fName
        stop 99
      endif
      return
      end

c --------------------------
C      ________________________________________________________
C     |                                                        |
C     |            SORT AN ARRAY IN INCREASING ORDER           |
C     |                                                        |
C     |    INPUT:                                              |
C     |                                                        |
C     |         X     --ARRAY OF NUMBERS                       |
C     |                                                        |
C     |         Y     --WORKING ARRAY (LENGTH  AT LEAST N)     |
C     |                                                        |
C     |         N     --NUMBER OF ARRAY ELEMENTS TO SORT       |
C     |                                                        |
C     |    OUTPUT:                                             |
C     |                                                        |
C     |         X     --SORTED ARRAY                           |
C     |________________________________________________________|
C
      SUBROUTINE SORT(X,Y,N)
      
      implicit none

      INTEGER I,J,K,L,M,N      
      REAL X(1),Y(1),S,T

      I = 1
10    K = I
20    J = I
      I = I + 1
      IF ( J .EQ. N ) GOTO 30
      IF ( X(I) .GE. X(J) ) GOTO 20
      Y(K) = I
      GOTO 10
30    IF ( K .EQ. 1 ) RETURN
      Y(K) = N + 1
40    M = 1
      L = 1
50    I = L
      IF ( I .GT. N ) GOTO 120
      S = X(I)
      J = Y(I)
      K = J
      IF ( J .GT. N ) GOTO 100
      T = X(J)
      L = Y(J)
      X(I) = L
60    IF ( S .GT. T ) GOTO 70
      Y(M) = S
      M = M + 1
      I = I + 1
      IF ( I .EQ. K ) GOTO 80
      S = X(I)
      GOTO 60
70    Y(M)= T
      M = M + 1
      J = J + 1
      IF ( J .EQ. L ) GOTO 110
      T = X(J)
      GOTO 60
80    Y(M) = T
      K = M + L - J
      I = J - M
90    M = M + 1
      IF ( M .EQ. K ) GOTO 50
      Y(M) = X(M+I)
      GOTO 90
100   X(I) = J
      L = J
110   Y(M) = S
      K = M + K - I
      I = I - M
      GOTO 90
120   I = 1
130   K = I
      J = X(I)
140   X(I) = Y(I)
      I = I + 1
      IF ( I .LT. J ) GOTO 140
      Y(K) = I
      IF ( I .LE. N ) GOTO 130
      IF ( K .EQ. 1 ) RETURN
      GOTO 40
      END


c -------------

 

