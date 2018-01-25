c -------------------------------------------------------------------------

      subroutine RdInput ( nFlt, nFlt0, f_start, f_num, faultFlag, nInten,  
     1           testInten, probAct, al_segWt, cumWt_SegModel, cumWt_Width, 
     2           cumWt_param, cumWt_ftype, cumWt_GM, nSegModel, nWidth, 
     3           nParamVar, nFtype, nGM_model, nAttenType, attenType, 
     4           nProb, iPer, SpecT1, version )

      implicit none
      include 'fract.h'

      integer nInten, ntotal, attentype(MAX_FLT), nfiles, ix(MAX_FILES),
     1        jCalc(MAX_ATTENTYPE,MAX_ATTEN), nParamVar(MAX_FLT,MAX_WIDTH),
     2        nWidth(MAX_FLT), nProb, nattentype, nGM_Model(MAX_ATTENTYPE), 
     3        nFtype(MAX_FLT), faultflag(MAX_FLT,MAX_SEG,MAX_FLT), j,
     4        f_start(MAX_FLT), f_num(MAX_FLT), nSegModel(MAX_FLT), jj
      integer iprob, iPer, nFlt, nFlt0, nwr
      real testInten(MAX_INTEN), probAct(MAX_FLT), minlat, maxlat,
     1     minlong, maxlong, maxdist, specT, dirflag, al_segWt(MAX_FLT),
     2     checkwt, c1, c2, wtgm(4,MAX_ATTEN), cumWt_SegModel(MAX_FLT,MAX_SEG),
     3     cumWt_param(MAX_FLT,MAX_WIDTH,MAXPARAM), sigtrunc, Varadd,
     4     cumWt_Width(MAX_FLT,MAX_WIDTH), cumWt_GM(MAX_ATTEN,MAX_ATTEN)
      real cumwt_FTYPE(MAX_FLT,MAX_FTYPE), SpecT1, version
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
      read (20,*) minlat,maxlat,minlong,maxlong

      read (20,*) maxdist

c     Input Title (not used) 
      read(20,'( a80)') title

c     Number of Spectral Periods and Number of attenuation relations types
      read(20,*) nProb, nattentype

      do iprob=1,nProb

C       Read period, maxeps dir flag and gm intensities
        read (20,*) specT, sigtrunc, dirflag 
        read (20,*) nInten, (testInten(j), j=1,nInten)
        call CheckDim ( nInten, MAX_INTEN, 'MAX_INTEN' )

C       Read in the suite of attenution models and wts for each attentype
        do j=1,nattentype
          checkwt = 0.0
          read (20,*) nGM_model(j)
c         Check for Max number of attenuation model
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

      ix(1) = 0


       return
       end
