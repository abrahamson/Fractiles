
      subroutine CheckDim ( n, nMax, name )
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
      real x(1)
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

      subroutine CheckWt1 ( x, n, j, n1, fName, name  )
      real x(n1,1), delta
      character*80 fName, name
      
      sum = 0.
      do i=1,n
        sum = sum + x(j,i)
      enddo
      delta = abs(sum - 1.0)
      if ( delta .gt. 0.01 ) then
        write (*,*) ' CheckWt1 Subroutine.'
        write (*,'( 2x,''Error -- Weights do not sum to unity'')')
        write (*,'( 2x,a80)') name
        write (*,'( 2x,a80)') fName
        write (*,*) ' Sum = ', sum
        do k=1,n
           write (*,*) k,x(j,k)
        enddo
        stop 99
      endif
      return
      end
      
      
