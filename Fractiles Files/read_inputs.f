
      subroutine RdInput ( nFlt, nParamVar, jCalc, nAtten, nInten,  
     1     testInten, nSite, lgTestInten, findex, probAct, nWidth, 
     2     nFlt0, nFlt1, nFlt2, cumWt_flt1, cumWt_param, cumWt_Thick, 
     3     nfiles, ix, cumWt_Dip, ndip, cumwt_atten, attentype, cumWt_Ftype, 
     4     nFtype, n_Dip, nThick1, nmodel, segwt1 )       

      include 'fract.h'

      real testInten(MAX_INTEN)
      real segModelWt(1), probAct(1), minlat,maxlat,minlong,maxlong,maxdist
      integer nInten, nWidth(1), ntotal, attentype(MAX_FLT)
      integer jCalc(MAX_ATTENTYPE,MAX_ATTEN),  nParamVar(MAX_FLT,MAX_WIDTH),iii,
     1        findex(MAX_FLT,MAX_FLT,MAX_FLT)
      character*80 filein, title, fname(1)
      integer nFlt1(1), nFlt2(MAX_FLT,1), nfiles, ix(MAX_FILES), ndip(MAX_FLT), n_dip(MAX_FLT)
      integer nThick1(MAX_FLT)
cnjg      real cumWt_flt1(MAX_FLT,1),cumWt_param(MAX_FLT,MAX_WIDTH,MAX_DIP,MAXPARAM),
      real cumWt_flt1(MAX_FLT),cumWt_param(MAX_FLT,MAX_WIDTH,MAX_DIP,MAXPARAM),
     1     cumWt_Thick(MAX_FLT,1), cumWt_Dip(MAX_FLT,MAX_DIP), cumWt_atten(MAX_ATTEN,MAX_ATTEN),
     2     cumwt_FTYPE(MAX_FLT,MAX_FTYPE)
      real x(100), specT, dirflag, vs30, depth10, depth15, SegWt1(MAX_FLT)

c      integer TreeIndex(MAX_FLT,MAX_DIP,MAX_WIDTH,MAX_N2,MAX_N2,MAX_N2,MAX_N2)
      integer nMaxMag(MAX_FLT,MAXPARAM), nSlipRate(MAX_FLT,MAXPARAM)
      integer iOverRide, nbvalue(MAX_FLT), nRecur(MAX_FLT)
      integer nMagLength, nMagArea

      integer nMagBins, nDistBins, nEpsBins, nXcostBins, soilampflag
      real magBins(MAX_MAG), distBins(MAX_DIST), epsBins(MAX_EPS)
      real XcostBins(MAX_XCOST)
      integer nProb, nattentype, nModel(MAX_ATTENTYPE), nFtype(MAX_FLT)
      real testwt, checkwt, c1, c2, wtgm(4,MAX_ATTEN)

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
      read (5,'( a80)') filein
      open (20,file=filein,status='old')

c     Open Input PSHA Source/Fault file
      read (20,'( a80)') filein
      open (10,file=filein,status='old')

c     Read in parameters for background grid.
      read (20,*) minlat,maxlat,minlong,maxlong
      read (20,*) maxdist

c     Input Title (not used) 
      read(20,'( a80)') title

c     Number of Spectral Periods and Number of attenuation relations types
      read(20,*) nProb, nattentype

C     Read period, maxeps dir flag and gm intensities
      read (20,*) specT, sigtrunc, dirflag 
      read (20,*) nInten, (testInten(j), j=1,nInten)
      call CheckDim ( nInten, MAX_INTEN, 'MAX_INTEN' )

C     Read in the suite of attenution models and wts for each attentype
      do j=1,nattentype
         checkwt = 0.0
         read (20,*) nmodel(j)
c     Check for Max number of attenuation model
         call checkDim ( nmodel(j), MAX_ATTEN, 'MAX_ATTEN' )

         do jj=1,nmodel(j)
cnjg            read (20,*) jcalc(jj,j), c1, c2, wtgm(j,jj)
            read (20,*) jcalc(j,jj), c1, c2, wtgm(j,jj)
            checkwt = checkwt + wtgm(j,jj)
         enddo

C     Check weight for each attentype group
         testwt = abs(checkwt-1.0)
         if (testwt .gt. 0.02) then
            write (*,*) 'Ground motion model weights do not sum to unity!!!'
            write (*,*) 'Need to check PSHA input file and PSHA may need'
            write (*,*) 'to be re-run with corrected weights.'
            write (*,*) 'Attenuation Type = ', j
            write (*,*) 'Number of Attenuation models = ', nmodel(j)
            write (*,*) 'Ground Motion Model Weight sum = ', checkwt
            stop 99
         endif
      enddo 

C     Now set up cum wts for attenuation models. 
      do j=1,nattentype
         do jj=1,nmodel(j)
            if (jj .eq. 1) then
               cumWt_atten(j,jj) = Wtgm(j,jj)
            else
               cumWt_atten(j,jj) = Wtgm(j,jj) + cumWt_atten(j,jj-1)
            endif
         enddo
      enddo

C*** Stuff not used and hence not needed ***
c      read (20,*) psCorFlag
cC     Read over the magnitude, distance, epsilon and Xcostheta bins.
c      read (20,*) nMagBins
c      read (20,*) (magBins(k),k=1,nMagBins)
c      read (20,*) nDistBins
c      read (20,*) (distBins(k),k=1,nDistBins)
c      read (20,*) nEpsBins
c      read (20,*) (epsBins(k),k=1,nEpsBins)
c      read (20,*) nXcostBins
c      read (20,*) (XcostBins(k),k=1,nXcostBins)
cC     Read over the soilampflag
c      read (20,*) soilAmpFlag
cC     Read in the number of sites to process.
c      read (20,*) nSite
c      call CheckDim ( nSite, MAX_SITE, 'MAX_SITE' )
      close (20)

      ix(1) = 0

c     Read Fault Data
      
      call Rd_Fault_Data ( nFlt, nFlt0,nFlt1,nFlt2,
     7     cumWt_flt1,cumWt_param,cumWt_Thick, nWidth, probAct,
     2     nParamVar, fIndex, cumWt_Dip, ndip, attentype, cumwt_ftype, 
     3     nFtype, n_DIP, nThick1, segwt1)

       return
       end
       
c ----------------------------------------------------------------------

      subroutine Rd_Fault_Data (nFlt,nFlt0,nFlt1,nFlt2,
     1     cumWt_flt1,cumWt_param,cumWt_Thick, nWidth, probAct,
     2     nParamVar, fIndex, cumWt_dip, nDip, AttenType, cumwt_Ftype, 
     3     nFtype, n_DIP, nThick1, segwt1 )

      include 'fract.h'
      integer fIndex(MAX_FLT,MAX_FLT,MAX_FLT), nWidth(1), ndip(MAX_FLT)
      integer  nParamVar(MAX_FLT,MAX_WIDTH), nRT(MAXPARAM)
      integer nMaxMag(MAX_FLT,MAXPARAM), nSlipRate(MAXPARAM), Attentype(MAX_FLT)
      real segwt(MAX_FLT,MAX_FLT), wt, al_segWt(MAX_FLT)
      real minmag, magstep, minDepth, coef_area(2), coef_width(2)
      real sliprateWt(MAXPARAM,MAXPARAM),  bValueWt(MAXPARAM),
     1     magRecurWt(MAXPARAM), faultWidthWt(MAXPARAM),
     2     maxMagWt(MAXPARAM,MAXPARAM), rt_wt(10,10),
     3     dipWt(MAXPARAM)
      integer nFlt, nFlt0, nFlt1(1), nFlt2
cnjg      real cumWt_flt1(MAX_FLT,1),cumWt_param(MAX_FLT,MAX_WIDTH,MAX_DIP,MAXPARAM),
      real cumWt_flt1(MAX_FLT),cumWt_param(MAX_FLT,MAX_WIDTH,MAX_DIP,MAXPARAM),
     1     cumWt_Thick(MAX_FLT,1), cumWt_Dip(MAX_FLT,MAX_DIP), cumwt_Ftype(MAX_FLT,MAX_FTYPE)
      real magRecur1(MAXPARAM)
      real x(MAXPARAM), x2(MAXPARAM,MAXPARAM)
      real ProbAct(1), topdepth
      character*80 fName1, fName
      character*1 tempname
      integer nFm, nFtypeModels, nFtype1(MAX_FLT)
      real ftmodelwt(MAX_FLT), ftype1(MAX_FLT,MAX_FLT), Ftype_wt1(MAX_FLT,MAX_FLT)

c      integer TreeIndex(MAX_FLT,MAX_DIP,MAX_WIDTH,MAX_N2,MAX_N2,MAX_N2,MAX_N2)
c      real tree_RecurWt(MAX_FLT,MAXPARAM), tree_bvalueWt(MAX_FLT,MAXPARAM), 
c     1     tree_SlipWt(MAX_FLT,MAXPARAM,MAXPARAM) 
c      real treeMagWt(MAX_FLT,MAX_WIDTH,MAX_FIXED_MAG,MAXRUP,2,MAX_MAGDIM)
c      integer treeMag(MAX_FLT,MAX_WIDTH,MAX_FIXED_MAG,MAXRUP,2,MAX_MAGDIM)
cnjg      real rupLength_wt(MAX_N1), magApproach_wt(MAX_FLT,2)
      real rupLength_wt(MAX_N1), magApproach_wt(MAX_FLT)
      real magLength_wt(MAX_N1), magArea_wt(MAX_N1)
      integer iOverRide, nbvalue(MAX_FLT), nRecur(MAX_FLT), nThick1(MAX_FLT)
      integer nMagLength, nMagArea
      integer iOverRideMag, nruplength

      real magsyn, rupsyn, jbsyn, seismosyn, hyposyn, wtsyn
      integer directflag, synflag

      integer nSR, nActRAte, nRecInt, n_Dip(MAX_FLT)
      real wt_srBranch, wt_ActBranch, wt_recIntBranch
      real wt_SR(MAXPARAM), wt_ActRate(MAXPARAM), wt_RecInt(MAXPARAM)
      real bValue2(MAXPARAM), actRate(MAXPARAM), actRateWt(MAXPARAM)
      real wt_MoRate(MAXPARAM), RateWt1(MAXPARAM), faultThickWt(MAX_WIDTH)
      real refMagWt(MAX_Width,MAX_Width), MagDisp_Wt(MAX_FLT), Disp_Wt(MAX_FLT)
      integer nRefMag(MAX_WIDTH), nFtype(MAX_FLT)
      real ftype(MAX_FLT,MAX_FLT), ftype_wt(MAX_FLT,MAX_FLT), SegWt1(MAX_FLT)
      integer faultflag(MAX_FLT,MAX_SEG,MAX_FLT)

C     Input Fault Parameters

      read (10,*) iCoor
      read (10,*) NFLT

      call CheckDim ( NFLT, MAX_FLT, 'MAX_FLT' )
  
      iflt = 0
      ifsytem = 1

C.....Loop over each fault in the source file....
      DO iFlt0=1,NFLT
        read (10,'( a80)') fName1
        read (10,*) probAct0

c       Read number of segmentation models for this fault system       
        read (10,*) nSegModel
        read (10,*) (segWt(iFlt0,k),k=1,nSegModel)

c       Read total number of fault segments 
        read (10,*) nFlt2
        do i=1,nSegModel
         read (10,*) (faultFlag(iFlt0,i,k),k=1,nFlt2)
        enddo

C.......Loop over number of individual fault segments....                        
        do iflt2=1,nflt2
          iFlt = iFlt + 1
          call CheckDim ( iflt, MAX_FLT, 'MAX_FLT   ' )

c       Set fault indexes
c          fIndex(1,iflt) = iflt0
c          fIndex(2,iflt) = iflt2

c       Find total weight for this segment
          sumseg = 0.
          do i=1,nSegModel
            sumseg = sumseg + segWt(iFlt0,i) * faultFlag(iFlt0,i,iFlt2)
          enddo
          segWt1(iFlt) = sumseg

C     Set up cum wt for fault segments...
          if (iFlt .eq. 1) then
             cumWt_flt1(iFlt) = segWt1(iFlt)
          else
             cumWt_flt1(iFlt) = segWt1(iFlt) + cumWt_flt1(iFlt-1)
          endif

c         write (*,*) 'Loading Flt Wts: ', iFlt, segwt1(iFlt), cumwt_flt1(iFlt)
c         
           
c           Set fault indexes to keep track of epistemic variability vs multple faults
            fIndex(iflt0,iflt1,iflt2) = iflt

c           Read past name of this segment
            read(10,'( a80)') fname
            write(*,'(2i6, a80)') iflt, iflt2, fName
            read (10,*) isourceType, attenType(iFlt), sampleStep, directflag, synflag

c           Read synchronous Rupture parameters
            if (synflag .gt. 0) then
               read (10,*) nsyn, attensyn
                    do insyn=1,nsyn
                       read (10,*) magsyn, rupsyn, jbsyn, seismosyn, hwsyn,
     1                             ftypesyn, hyposyn, wtsyn
                    enddo
            endif

c           Read aleatory segmentation wts
            read (10,*) al_segWt(iFlt)

c       Check for standard fault source or areal source
        if ( isourceType .eq. 1 .or. isourceType .eq. 2) then
          read (10,*) dip1, top
          read(10,*) nfp     
          do ipt=1,nfp
            read (10,*) fLong, fLat
          enddo
        endif

c        Check for grid source (w/o depth)
        if ( isourceType .eq. 3 ) then
           read (10,*)  dip1, top
c     Read over grid filename...
           read (10,*)
        endif

c        Check for grid source (w/ depth)
        if ( isourceType .eq. 4 ) then
              read (10,*)  dip1
c     Read over grid filename...
              read (10,*)
        endif

c        Check for custom fault source
         if ( isourceType .eq. 5) then
           read(10,*) nDownDip, nfp
c........Only read in the first downdip case - rest is not needed...
           do ipt=1,nfp
              read (10,*) fLong, fLat, fZ
           enddo
         endif

c        Read dip Variation
         if ( isourceType .ne. 5 ) then
           read (10,*) n_Dip(iflt)
           read (10,*) (x(i),i=1,n_Dip(iflt))
           read (10,*) (dipWt(i),i=1,n_Dip(iflt))
         else
           n_Dip(iflt) = 1
           x(1) = 0.
           dipWt(1) = 1.
         endif

            do iDip=1,n_Dip(iflt)
              if ( iDip .eq. 1) then
                cumWt_dip(iFlt,iDip) = dipWt(iDip)
              else
                cumWt_dip(iFlt,iDip) = cumWt_dip(iFlt,iDip-1)
     1                   + DipWt(iDip)
              endif
            enddo

c        Read b-values (not for activity rate cases)
         read (10,*) n_bValue
         call CheckDim ( n_bValue, MAX_N1, 'MAX_N1    ' )
         if ( n_bValue .gt. 0 ) then
           read (10,*) (x(i),i=1,n_bValue)
           read (10,*) (bValueWt(i),i=1,n_bValue)
         endif
                
c        Read activity rate - b-value pairs
         read (10,*) nActRate 
         if ( nActRate .ne. 0 ) then
           do ii=1,nActRate
             read (10,*) bValue2(ii), actRate(ii), actRateWt(ii)
           enddo
         endif

c        Read weights for rate methods
         read (10,*) wt_srBranch, wt_ActRateBranch, wt_recIntBranch, wt_MoRateBranch
         sum = wt_srBranch + wt_ActRateBranch + wt_recIntBranch + wt_MoRateBranch
         if ( sum .lt. 0.999 .or. sum .gt. 1.001 ) then
              write (*,'( 2x,''rate method weights do not sum to unity for fault, '',a30)') fName
              stop 99
         endif

c        Read slip-rates
         read (10,*) nSR
         if ( nSR .gt. 0 ) then
           read (10,*) (x(k),k=1,nSR)
           read (10,*) (wt_sr(k),k=1,nSR)
         endif

c        Read recurrence intervals
         read (10,*) nRecInt
         if ( nRecInt .gt. 0 ) then
           read (10,*) (x(k),k=1,nRecInt)
           read (10,*) (wt_recInt(k),k=1,nRecInt)
         endif

c        Read moment-rates
         read (10,*) nMoRate
         if ( nMoRate .gt. 0 ) then
           read (10,*) (x(k),k=1,nMoRate)
           read (10,*) (x(k),k=1,nMoRate)
           read (10,*) (wt_MoRate(k),k=1,nMoRate)
         endif
              
c        Load into single array called "rate_param"
         nRate = nSR + nActRate + nRecInt + nMoRate
         call CheckDim ( nRate, MAX_N1, 'MAX_N1    ' )
         do k=1,nSR
c            rateParam1(k) = sr(k)
            rateWt1(k) = wt_sr(k)*wt_srbranch              
c            rateType1(k) = 1
c            temp_BR(k) = k
c            temp_BR_wt(k) = wt_sr(k)
         enddo
         do k=1,nActRate
c           rateParam1(k+nSR) = actRate(k)
           rateWt1(k+nSR) = actRateWt(k) * wt_actRateBranch             
c           rateType1(k+nSR) = 2
c           temp_BR(k+nSR) = k
c           temp_BR_wt(k+nSR) = actRateWt1(k)
         enddo
         do k=1,nRecInt
c           rateParam1(k+nSR+nActRate) = rec_Int(k)
           rateWt1(k+nSR+nActRate) = wt_recInt(k) *wt_recIntBranch              
c           rateType1(k+nSR+nActRate) = 3
c           temp_BR(k+nSR+nActRate) = k
c           temp_BR_wt(k+nSR+nActRate) = wt_recInt(k)
         enddo
         do k=1,nMoRate
c           rateParam1(k+nSR+nActRate+nRecInt) = MoRate(k)
           rateWt1(k+nSR+nActRate+nRecInt) = wt_MoRate(k) *wt_MoRateBranch              
c           rateType1(k+nSR+nActRate+nRecInt) = 4
c           temp_BR(k+nSR+nActRate+nRecInt) = k
c           temp_BR_wt(k+nSR+nActRate+nRecInt) =  wt_MoRate(k)
         enddo
                     
c        Read Mag recurrence weights (char and exp)
         read (10,*) nMagRecur
         call CheckDim ( nMagRecur, MAX_N1, 'MAX_N1    ' )
         read (10,*) (x(i),i=1,nMagRecur)
         read (10,*) (magRecurWt(i),i=1,nMagRecur)

c        Read in corresponding magnitude parameters. 
         do iRecur=1,nMagRecur
           if (x(iRecur) .eq. 4) then
             read (10,*) x(iRecur), x(iRecur), x(iRecur), x(iRecur), x(iRecur)
           else              
             read (10,*) x(iRecur), x(iRecur), x(iRecur)
             x(iRecur) = 0.0
             x(iRecur) = 0.0
           endif
         enddo

c        Read seismogenic thickness
         if ( isourceType .ne. 5) then
           read (10,*) nThick1(iFlt)
           call CheckDim ( nThick1(iFlt), MAX_WIDTH, 'MAX_WIDTH ' )
           read (10,*) (x(i),i=1,nThick1(iFlt))
           read (10,*) (faultThickWt(i),i=1,nThick1(iFlt))
         else
           nThick1(iFlt) = 1
           x(1) = -99.
           faultThickWt(1) = 1.
         endif
         
            do iThick=1,nThick1(iFlt)
              if ( iThick .eq. 1) then
                cumWt_Thick(iFlt,iThick) = faultThickWt(iThick)
              else
                cumWt_Thick(iFlt,iThick) = cumWt_Thick(iFlt,iThick-1)
     1                   + faultThickWt(iThick)
              endif
            enddo

c        Read depth pdf
         read (10,*) iDepthModel, (x2(iflt,k),k=1,3)       

c        Read Mag method (scaling relations or set values)
         read (10,*) iOverRideMag

         if ( iOverRideMag .eq. 1 ) then
c         Read reference mags for each fault thickness
          iThick = 1
          do iThick1=1,nThick1(iFlt)
            read (10,*) nRefMag(iThick)
            read (10,*) (x2(iThick,i),i=1,nRefMag(iThick))
            read (10,*) (refMagWt(iThick,i),i=1,nRefMag(iThick))
            if ( nRefMag(iThick) .ne. 0 ) then
            endif
              iThick = iThick + 1
              
c           Copy these ref magnitudes for each addtional dip (no correlation allowed)
            do iDip=2,n_Dip(iflt)
              nRefMag(iThick) = nRefMag(iThick-1)
              do i=1,nRefMag(iThick)
                x2(iThick,i) = x2(iThick-1,i)
                refMagWt(iThick,i) = refMagWt(iThick-1,i)
              enddo
              iThick = iThick + 1
            enddo
          enddo

         else
c   *** New read for Mean Char Eqk computed in Code 	
c       (i.e. OverridgeMag .ne. 1) 
              read (10,*) nRupLength
              read (10,*) (x(k),k=1,nRupLength)
              read (10,*) (rupLength_wt(k),k=1,nRupLength)
              read (10,*) (magApproach_wt(k),k=1,3)

c             read mag-length models
              read (10,*) nMaglength
              if ( nMagLength .ne. 0 ) then
                do k=1,nMagLength
                  read (10,*) (x2(kk,k), kk=1,4)
                enddo
                read (10,*) (magLength_wt(k),k=1,nMagLength)
              endif

c             read mag-area models
              read (10,*) nMagArea
              if ( nMagArea .ne. 0 ) then
                do k=1,nMagArea
                  read (10,*) (x2(kk,k), kk=1,4), aspectMin
                enddo
                read (10,*) (magArea_wt(k),k=1,nMagArea)
              endif
               
c             read mag-displacement models
              read (10,*) nMagDisp
              if ( nMagDisp .ne. 0 ) then
                do k=1,nMagDisp
                  read (10,*) (x2(kk,k), kk=1,4)
                enddo
                read (10,*) (magDisp_wt(k),k=1,nMagDisp)
                read (10,*) nDisp
                read (10,*) (x(k), k=1,nDisp)
                read (10,*) (Disp_wt(k),k=1,nDisp)
              endif
         endif

c           Read Past remaining input for this fault
         read (10,*) minMag, magStep, hxStep, 
     1              hyStep, nRupArea, nRupWidth, minDepth
         read (10,*) (x2(k,iFlt),k=1,2), sigArea
         read (10,*) (x2(k,iFlt),k=1,2), sigWidth

c        Read ftype Models
         read (10,*) nFtypeModels
         do iFM=1,nFtypeModels
            read (10,*) ftmodelwt(iFM)
            read (10,*) nFtype1(iFM)
            read (10,*) (ftype1(iFM,k),k=1,nFtype1(iFM))
            read (10,*) (ftype_wt1(iFM,k), k=1,nFtype1(iFM))
         enddo

c        Load up Ftype Models and weights.
         nFM = 1
         do iFM=1,nFtypeModels
            do k=1,nFtype1(iFM)
               ftype(iFlt,nFM) = ftype1(iFM,k)
               ftype_wt(iFlt,nFM) = ftype_wt1(iFM,k)*ftmodelwt(iFM)
               nFm = nFm + 1
            enddo
         enddo
         nFm = nFm - 1
         nFtype(iFlt) = nFm

C     Set up Cum Wts for Ftype
         do iFtype=1,nFtype(iFlt)
           if (iFtype .eq. 1) then
             cumWt_Ftype(iFlt,iFtype) = Ftype_Wt(iFlt,iFtype)
           else
             cumWt_Ftype(iFlt,iFtype) = cumWt_Ftype(iFlt,iFtype-1)
     1                + Ftype_Wt(iFlt,iFtype)
           endif
c           write (*,'(a12,2(i4,2x),f8.4)') 'Fault Type: ', iFtype,nftype(iFlt),cumwt_Ftype(iFlt,iFtype)
         enddo

c     Load up parameter variations into large single dimension arrays        
          i1 = 0
	  do iDip=1,n_Dip(iFlt)
            do iThick=1,nThick1(iFlt)

c             do iftype=1,nFtype(iFlt)
              i1 = i1 + 1
              i = 0

              do iRecur=1,nMagRecur
                do i_bValue=1,n_bValue
                  do iSlipRate=1,nSr+nActRate+nMoRate+nRecInt
                    do iRefMag=1,nRefMag(iThick)
                         i = i + 1
                         call CheckDim ( i, MAXPARAM, 'MAXPARAM' )

c                     Set the weight for this set of parameters
                      wt = RateWt1(iSlipRate)
     1                    * bValueWt(i_bValue)
     2                    * magRecurWt(iRecur) 
     3                    * refMagWt(iThick,iRefMag) 

                      if ( i .eq. 1 ) then
                         cumWt_param(iFlt,iThick,iDip,i) = wt
                      else
                        cumWt_param(iFlt,iThick,iDip,i) = wt 
     1               +  cumWt_param(iFlt,iThick,iDip,i-1)
                      endif
    
c     End of Loop over iRefMag  
                  enddo
c     End of Loop over ib_value
                enddo
c     End of Loop over iRate
	      enddo

c     End of Loop over iRecur
	    enddo
            nParamVar(iFlt,iThick) = i
            if ( cumwt_param(iFlt,iThick,iDip,nParamVar(iFlt,iThick)) .gt. 0.99 ) then
	      cumwt_param(iFlt,iThick,iDip,nParamVar(iFlt,iThick)) = 1.0
            endif
c	    write (*,'( 2x,''nParam='',5i5,f10.4)') iFlt, iThick,iDip,iFtype,
c     1            nParamVar(iFlt,iThick), CUMWT_PARAM(IFLT,iThick,iDip,nParamVar(iFlt,iThick))

c     End of Loop over iFtype    
c             enddo

c     End of Loop over iThick1
            enddo
            probAct(iFlt) = probAct0
c     End of Loop over iDip
          enddo
c     End of Loop over iFlt2 - number of segments    
        enddo
c     End of Loop over iFlt
      enddo
      nFlt = iFlt

      return
      end

