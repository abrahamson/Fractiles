      subroutine read_logicRisk ( isite, nSite, nFlt, nfiles, ix, risk, nParamVar, nAtten, nFtype)

      include 'fract.h'

      real risk(MAX_INTEN, MAX_ATTEN, MAX_FLT, MAX_WIDTH, MAXPARAM,MAX_FTYPE)
      real temp(MAX_INTEN)
      integer isite, nInten, nFlt,  nSite, nAtten(MAX_FLT),ix(MAX_FILES), iProb,
     1        nParamVar(MAX_FLT,MAX_WIDTH), nWidth(MAX_FLT), nfiles, iA1, ntotal
      integer nFtype(MAX_FLT)
      character*80 fName, file1
      nwr = 12

      nfiles = 1

C     Loop over the number of files.
      do 100 i1=1,nfiles
                     
c     Open output file
      read (5,'( a80)') file1

      write (*,*) 'Opening output file from the hazard runs.'
      write (*,*) file1
      open (nwr,file=file1,status='old')

C      Loop over nSite taken out since each site is now given its own output file
c           from the hazard code.
c      do jSite=1,nSite
            
      do iFlt=1,nFlt      
c         read (nwr,*) jFlt, nProb, nAtten, nWidth(iFlt), nFtype(iFlt), 
c     1        (nParamVar(iFlt,i,1),i=1,nWidth(iFlt)),
c     2        nInten

         read (nwr,*) jFlt, nProb, nAtten(iFlt), nWidth(iFlt), nFtype(iFlt), 
     1         (nParamVar(iFlt,i),i=1,nWidth(iFlt)),
     2         nInten

C     Check for Array Dimension of Risk Array.
        if (nFlt .gt. MAX_FLT) then
           write (*,*) 'MAX_FLT needs to be increased to ', nFlt
           write (*,*) 'Change Array Parameter in FRACT.H and recompile.'
           stop 99
        endif
        if (nWidth(iFlt) .gt. MAX_WIDTH) then
           write (*,*) 'MAX_WIDTH needs to be increased to ', nWidth(iFlt)
           write (*,*) 'Change Array Parameter in FRACT.H and recompile.'
           stop 99
        endif
        if (nProb .gt. MAX_PROB) then
           write (*,*) 'MAX_PROB needs to be increased to ', nProb
           write (*,*) 'Change Array Parameter in FRACT.H and recompile.'
           stop 99
        endif
        if (nFtype(iFlt) .gt. MAX_FTYPE) then
           write (*,*) 'MAX_FTYPE needs to be increased to ', nFtype(iFlt)
           write (*,*) 'Change Array Parameter in FRACT.H and recompile.'
           stop 99
        endif

        do iProb=1,nProb
          do iAtten=1,nAtten(iFlt)
            do iWidth=1,nWidth(iFlt)
               do iFtype=1,nFtype(iFlt)
                  do i=1,nParamVar(iFlt,iWidth)
 	             read (nwr,*) (temp(j),j=1,nInten)

	             do j=1,nInten
	                Risk(j,iAtten,iFlt,iWidth,i,iFtype) = temp(j)
	             enddo

                  enddo
               enddo
	    enddo
          enddo
        enddo

      enddo

      write (*,*) 'In logic Read Loop = ', nInten

c        if (jSite .eq. iSite ) then
c          close ( nwr )
c          return
c        endif

      close (nwr)
c      enddo

  100 continue

      return
  
      write (*,'( 2x,''bad site number'')')
      stop 98

      end

