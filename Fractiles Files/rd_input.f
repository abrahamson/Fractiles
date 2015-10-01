c -------------------------------------------------------------------------

      subroutine RdInput ( nFlt, nFlt0, f_start, f_num, faultFlag, 
     1     nInten,  testInten, lgTestInten, probAct, al_segWt, 
     3     cumWt_SegModel, cumWt_Width, cumWt_param, cumWt_ftype, cumWt_GM, 
     4     nSegModel,      nWidth,      nParamVar,   nFtype,      nGM_model, nAttenType, attenType, nProb, iPer )

      include 'fract.h'

      real testInten(MAX_INTEN)
      real segModelWt(1), probAct(1), minlat,maxlat,minlong,maxlong,maxdist
      integer nInten, ntotal, attentype(MAX_FLT)
      integer jCalc(MAX_ATTENTYPE,MAX_ATTEN),  nParamVar(MAX_FLT,MAX_WIDTH),iii
      character*80 filein, title, fname(1), dummy
      integer nFlt1(1), nFlt2, nfiles, ix(MAX_FILES), ndip(MAX_FLT), n_dip(MAX_FLT)
      integer nWidth(MAX_FLT)
      real cumWt_SegModel(MAX_FLT,1),cumWt_param(MAX_FLT,MAX_WIDTH,MAX_DIP,MAXPARAM),
     1     cumWt_Width(MAX_FLT,1), cumWt_GM(MAX_ATTEN,MAX_ATTEN),
     2     cumwt_FTYPE(MAX_FLT,MAX_FTYPE)
      real x(100), specT, dirflag, vs30, depth10, depth15, SegWt1(MAX_FLT), al_segWt(1)

c      integer TreeIndex(MAX_FLT,MAX_DIP,MAX_WIDTH,MAX_N2,MAX_N2,MAX_N2,MAX_N2)
      integer nMaxMag(MAX_FLT,MAXPARAM), nSlipRate(MAX_FLT,MAXPARAM)
      integer iOverRide, nbvalue(MAX_FLT), nRecur(MAX_FLT)
      integer nMagLength, nMagArea

      integer nMagBins, nDistBins, nEpsBins, nXcostBins, soilampflag
      real magBins(MAX_MAG), distBins(MAX_DIST), epsBins(MAX_EPS)
      real XcostBins(MAX_XCOST)
      integer nProb, nattentype, nGM_Model(MAX_ATTENTYPE), nFtype(MAX_FLT)
      real testwt, checkwt, c1, c2, wtgm(4,MAX_ATTEN)
      integer faultflag(MAX_FLT,MAX_SEG,MAX_FLT)
      integer f_start(1), f_num(1), nSegModel(1)
 


c     Set Data file units
      nwr = 11


      ntotal = 0

c     Read in the number of data files.
c      read (5,*) nfiles
c     Program no longer allowed to read from multiple files.
      nFiles = 1

c     Loop over the number of files.
c      do 111 iii=1,nfiles

c     Open PSHA Run Input File
      read (31,'( a80)') filein
      write (*,'( a80)') filein
      open (20,file=filein,status='old')

c     Open Input PSHA Source/Fault file
      read (20,'( a80)') filein
      open (10,file=filein,status='old')

c     Read in parameters for background grid.
      read (20,*) minlat,maxlat,minlong,maxlong

Cnjg  Added back in read of single maxdist
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


c     Read Fault Data
           call Rd_Fault_Data ( nFlt, nFlt0,
     7     cumWt_segModel, cumWt_param, cumWt_width, probAct,
     2     nParamVar, attentype, cumwt_ftype, 
     3     nFtype, nWidth, nSegModel,  f_start, f_num, faultFlag, al_segwt)
       

       return
       end
