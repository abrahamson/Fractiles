
      subroutine GetRandom0 ( iseed, n, wt, iSave )

      integer iseed,isave,n
      real wt(1)
      real x

c     Get random number
      x = ran1( iseed )

c      write (*,*) 'N=', n

      do i=1,n
        if ( x .le. wt(i) ) then
          iSave = i
          return
        endif
      enddo
      
      write (*,*) ' Get Random Number 0'
      write (*,*) 'Weights = ',wt
      write (*,*) 'Random Number = ', x
      write (*,'( 2x,''Error - bad ran number or weights'')')
      stop 99
      end

c ----------------------------
      subroutine GetRandom1 ( iseed, n, wt, i1, iSave, n1 )

      integer iseed,i1,isave,n1,n
      real wt(n1, n1)
      real x

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
      stop 99
      end
      
c ----------------------------
      subroutine GetRandom1b ( iseed, n, wt, i1, iSave, n1, n2 )

      integer iseed,i1,isave,n1,n, n2
      real wt(n2, n1)
      real x

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

      include 'FRACT.H'

      real wt(n1, n2, MAXPARAM), x
      integer n1, n2
      
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
      subroutine GetRandom3 ( iseed, n, wt, i1, i2, i3, iSave, n1, n2, n3 )

      include 'FRACT.H'

      real wt(n1, n2, n3, MAXPARAM), x
      integer n1, n2, n3
      
c     Get random number
      x = ran1( iseed )

      do i=1,n
        if ( x .le. wt(i1,i2,i3,i) ) then
          iSave = i
          return
        endif
      enddo
      
      write (*,*) ' Get Random Number 3'
      write (*,'( 2x,''Error - bad ran number 3 or weights'')')

      write (*,*) 'Random Number = ', x
      write (*,'(2x,''wts:'',10f10.4)') (wt(i1,i2,i3,i),i=1,n)
      write (*,'( 5i5)') n, i1, i2, i3

      stop 99
      end

c ----------------------------

      function Ran1 ( idum )

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

