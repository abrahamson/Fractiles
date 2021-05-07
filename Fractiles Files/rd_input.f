c -------------------------------------------------------------------------

      subroutine RdInput (nInten, testInten, nGM_model, cumWt_GM, nAttenType,
     1                    attenType, nProb, iPer, SpecT1, version, PCflag )

      implicit none
      include 'fract.h'

      integer nInten, ntotal, attentype(MAX_FLT), nfiles, PCflag(MAX_PROB),
     1        jCalc(MAX_ATTENTYPE,MAX_ATTEN), nProb, nattentype,
     2        nGM_Model(MAX_ATTENTYPE), j, jj, iprob, iPer, nwr
      real testInten(MAX_INTEN), dummy, specT, dirflag,
     2     checkwt, c1, c2, wtgm(4,MAX_ATTEN), sigtrunc, Varadd,
     4     cumWt_GM(MAX_ATTEN,MAX_ATTEN), SpecT1, version
      character*80 filein, title

c     Set Data file units
      nwr = 11

      ntotal = 0

c     Program no longer allowed to read from multiple files.
      nFiles = 1

c     Open PSHA Run Input File
      read (31,'( a80)') filein
      write (*,'( a80)') filein
      open (20,file=filein,status='old')

c     Open Input PSHA Source/Fault file
      read (20,'( a80)') filein
      open (10,file=filein,status='old')

C     Check for version compatibility with hazard code
        read (20,*) version
         if (version .ne. 45.3 .and. version .ne. 45.2 .and. version .ne. 45.1) then
         write (*,*) 'Hazard faultfile format incompatible, use Haz45.3, Haz45.2, or Haz45.1'
         stop 99
        endif

c     Read in parameters for background grid.
      read (20,*) dummy
      read (20,*) dummy

c     Input Title (not used)
      read(20,'( a80)') title

c     Number of Spectral Periods and Number of attenuation relations types
      read(20,*) nProb, nattentype
      call checkDim ( nProb, MAX_PROB, 'MAX_PROB' )
      call checkDim ( nattentype, MAX_ATTENTYPE, 'MAX_ATTENTYPE' )

      do iprob=1,nProb

C       Read period, maxeps, dir flag, PC flag, and gm intensities
        read (20,*) specT, sigtrunc, dirflag, PCflag(iProb)
        read (20,*) nInten
        call CheckDim ( nInten, MAX_INTEN, 'MAX_INTEN' )
        backspace (20)
        read (20,*) nInten, (testInten(j), j=1,nInten)


C       Read in the suite of attenution models and wts for each attentype
        do j=1,nattentype
          checkwt = 0.0
          read (20,*) nGM_model(j)
          call checkDim ( nGM_model(j), MAX_ATTEN, 'MAX_ATTEN' )

          do jj=1,nGM_model(j)
            read (20,*) jcalc(j,jj), c1, c2, wtgm(j,jj), Varadd
          enddo
        enddo

c       keep the weights for the selected problem only
        if ( iper .eq. iProb) then
         specT1 = specT
C         Set up cum wts for attenuation models.

          do j=1,nattentype
            do jj=1,nGM_model(j)
              if (jj .eq. 1) then
                cumWt_GM(j,jj) = Wtgm(j,jj)
              else
                cumWt_GM(j,jj) = Wtgm(j,jj) + cumWt_GM(j,jj-1)
              endif
            enddo
          enddo
        endif

      enddo

      close (20)

       return
       end
