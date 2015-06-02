      program Fract_Haz46

c     This program will compute the fractiles from a seismic hazard
c     run. This version of the program works with the seismic hazard
c     code version 42 or later. The fractiles are computed based on 
c     a Monte Carlo simulation. This version will be able to read fault
c     files in which synchronous rupture is modeled. 
 
c     compatible with Haz42
c     Last modified: 7/10

      include 'fract.h'
c      implicit none
      
      real risk(MAX_INTEN,MAX_ATTEN, MAX_FLT, MAX_WIDTH, MAXPARAM,MAX_FTYPE), 
     1     testInten(MAX_INTEN)
c      real probAct(MAX_FLT)
c      real cum_wt(5)
      integer nInten, jcalc(MAX_ATTENTYPE,MAX_ATTEN), nAtten(MAX_ATTEN), nsite, nFlt,
     1        isite, iflt, 
     1        iAtten, iInten, nWidth(MAX_FLT)
      integer iParam, nParamVar(MAX_FLT,MAX_WIDTH)
      integer fIndex(MAX_FLT,MAX_FLT,MAX_FLT), iFltWidth, iX(MAX_FILES)
      character*80 filein, file1, fileinPSHA
      integer nFlt1(MAX_FLT), nFlt2(MAX_FLT,MAX_FLT), attentype(MAX_FLT), nmodel(MAX_ATTENTYPE)
cnjg      real cumWt_flt1(MAX_FLT,MAX_FLT),
      real cumWt_flt1(MAX_FLT),
     1     cumWt_param(MAX_FLT,MAX_WIDTH,MAX_DIP,MAX_FTYPE,MAXPARAM),
     1     cumWt_Thick(MAX_FLT,MAX_WIDTH), cumWt_Ftype(MAX_FLT,MAX_FTYPE)
      real ran1, risk1(MAX_SAMPLE,MAX_INTEN)
      real sortarray(MAX_SAMPLE), attenWt(MAX_PROB)
      real cumWt_atten(MAX_ATTEN,MAX_ATTEN), cumwt_dip(MAX_FLT,MAX_DIP), SegWt1(MAX_FLT)
      integer Fflag 
      integer iperc, nfiles, n_Dip(MAX_FLT)
      real perc(nPerc,MAX_INTEN),mean(MAX_INTEN) 
      integer nProb1(MAX_ATTENTYPE), jProb(MAX_PROB), kProb, nDip(MAX_FLT)
      integer iseed, nSample, nFtype(MAX_FLT), nThick1(MAX_FLT)

c      integer TreeIndex(MAX_FLT,MAX_DIP,MAX_WIDTH,MAX_N2,MAX_N2,MAX_N2,MAX_N2)      

c      open (33,file='debug.out')
    
      write (*,*) '*************************'
      write (*,*) '* Fractile Code for use *'
      write (*,*) '*  Hazard_v42.exe code  *'
      write (*,*) '*      July, 2010       *'
      write (*,*) '*************************'

      write (*,*) 'Enter the input filename.'
      
      read (*,'(a80)') filein
      open (5,file=filein,status='old')
      pause 'test -1'

      read (5,*) iseed
      read (5,*) nSample
      pause 'test 0'

c     Enter the input filename from the PSHA run.
c      read (5,*) fileinPSHA

      read (5,*) nAttentype
      read (5,*) (nProb1(j), j=1,nattentype)
      pause 'test 1'

c      read (5,*) (attenWt(k),k=1,nProb1)

c      cumWt_atten(1) = attenWt(1)
c      do i=2,nProb1
c        cumWt_atten(i) = cumWt_atten(i-1) + attenWt(i)
c      enddo
      
      call CheckDim ( nSample, MAX_SAMPLE, 'MAX_SAMPLE ' )

c     Read Input File
      call RdInput ( nFlt, nParamVar, jcalc, nAtten, nInten,  
     1     testInten, nSite, lgTestInten, findex, probAct, nWidth, 
     2     nFlt0, nFlt1, nFlt2, cumWt_flt1, cumWt_param, cumWt_Thick,
     3     nfiles, ix, cumWt_Dip, nDip, CumWt_atten, attentype, cumwt_Ftype, 
     4     nFtype, n_dip, nThick1, nmodel, segwt1 )

C     Write out the intensity values if requested.
c       if (outflag .eq. 1) then
c          write (17,*) nInten, (testInten(j2), j2=1,nInten) 
c       endif

c     Loop Over Number of Sites
      nsite = 1
      do 1000 iSite = 1, nSite

c       Initialize risk1 and mean
        do i=1,nSample
          do j=1,MAX_INTEN
            risk1(i,j) = 0.
          enddo
        enddo

c       Initialize mean
        do i=1,MAX_INTEN
          mean(i) = 0.
        enddo

c       Initialize random number generator
        do i=1,500
          tempx = Ran1 (iseed)
        enddo

c       Each "set" is the file (mixure and mno mixture)
        do 900 iSet=1,nSet

        write (*,'(2x,a18,i4,a5,2x,i4)') 'Looping over site ', isite,'  of ', nsite
c       Read site location information from hazard input file (data not used)
c        read (20,*) SX, SY, vs30, depth10, depth15
c        do j=1,4
c           read (20,*)
c        enddo

        write (*,'( 2x,''reading input file'')')
        call read_logicRisk ( isite, nSite, nFlt, nfiles, ix, risk, nParamVar, natten, nFtype )
        call read_logicRisk ( isite, nSite, nFlt, nfiles, ix, risk_nomix, nParamVar, natten, nFtype )

c       write (*,*) 'Done reading logic Risk Data.'

c       Monte Carlo Sampling of Hazard
        do iSample=1,nSample
          write (*,*) 'Looping over samples:        ', iSample, ' of ', nSample

c         Sample the hazard from each geometrically independent source region (or fault)
c         (a source region may include multiple subsources)
          do iFlt0 = 1, nFlt
c             write (*,*) 'Looping over Faults:        ', iFlt0, ' of ', nFlt

c           Select the source geometry model for this source region
c            write (*,*) 'Calling Fault wts'
c            call GetRandom0 ( iseed, nFlt, cumWt_Flt1, mcFlt)
             mcFlt = 1
c            write (*,*) 'Selected Fault Number ', mcFlt, ' out of ', nFlt

c           Sample the attenuation relation
c            write (*,*) 'Calling Atten wts'
            call GetRandom1 ( iseed, nmodel(attentype(iFlt0)), cumWt_atten, attentype(iFlt0), 
     1                        mcAtten, MAX_ATTEN ) 
c            write (*,*) 'Selected Attenuation Number', mcAtten, ' out of ', nAtten(iFlt0)

c           Select the fault dip for this subsource 
c            write (*,*) 'Calling Dip wts'
            call GetRandom1b ( iseed, n_Dip(iFlt0), cumWt_dip,
     1                        iflt0, mcDip, MAX_DIP, MAX_FLT )
c            write (*,*) 'Selected ', mcDip, ' out of ', n_Dip(iFlt0)

c           Select the fault mechanism type 
c            write (*,*) 'Calling Ftype wts'
            call GetRandom1 ( iseed, nFtype(iFlt0), cumWt_Ftype,
     1                        iflt0, mcFTYPE, MAX_FLT )
c            write (*,*) 'Selected mcFtype', mcFtype, ' out of ', nFtype(iFlt0)

c           Select the fault width for this subsource 
c            write (*,*) 'Calling Width wts'
            call GetRandom1 ( iseed, nThick1(iFlt0), cumWt_Thick,
     1                      iflt0, mcWidth, MAX_FLT )
c            write (*,*) 'Selected ', mcWidth, ' out of ', nThick1(iFlt0)

c           Select the parameter variation for this subsource and fault width
c            write (*,*) 'Calling Param wts'
c            call GetRandom3 (iseed,nParamVar(mcFlt,mcWidth), 
c     1          cumWt_param,mcFlt,mcWidth,mcDip,mcParam,MAX_FLT,MAX_WIDTH,MAX_DIP)
            call GetRandom3 (iseed,nParamVar(iFlt0,mcWidth), 
     1          cumWt_param,iFlt0,mcWidth,mcDip,mcParam,MAX_FLT,MAX_WIDTH,MAX_DIP)
c            write (*,*) 'Selected ', mcParam, ' out of ', nParamVar(mcFlt,mcWidth)

c      write (*,*) 'Sampled Cases:'
c      write (*,'(7i8,f8.3)') iSample, iFlt0, mcAtten, mcFlt, mcWidth, mcParam, mcFtype, segwt1(iFlt0)

c             Add the sampled hazard curve to the total hazard
              do iInten=1,nInten
                risk1(iSample,iInten) = risk1(iSample,iInten) 
     1                    + risk(iInten,mcAtten,iFlt0,mcWidth,mcParam,mcFtype)*segwt1(iFlt0)*0.8
     1                    + risk_nomix(iInten,mcAtten,iFlt0,mcWidth,mcParam,mcFtype)*segwt1(iFlt0)*0.2

     1        
                mean(iInten)=mean(iInten)
     1	 	          + risk(iInten,mcAtten,iFlt0,mcWidth,mcParam,mcFtype)*segwt1(iFlt0)*0.8
     1	 	          + risk_nomix(iInten,mcAtten,iFlt0,mcWidth,mcParam,mcFtype)*segwt1(iFlt0)*0.2
              enddo
            enddo

c        write (*,*) 'Done with sample: ', iSample
        enddo
  900   continue

c       Compute mean hazard as check between fractile and haz code
         do iInten = 1,nInten
           mean(iInten) = mean(iInten)/nSample
        enddo

c       Sort risk data so that cummulative empirical distribution function
c       can be computed  (sorts in place)
c       sortarray is a working array
        do iInten = 1,nInten
           call sort(risk1(1,iInten),sortarray,nSample)
        enddo

c       Compute fractile values
        step = 1./(nperc+1)
        do iInten=1,nInten
          do iperc = 1,nPerc
            i1 = int(iperc*step*nSample)
	  perc(iperc,iInten) = risk1(i1,iInten)
	enddo
        enddo

c       Write output
        read (5,'( a80)') file1
        write(*,'( a32,a80)') 'Writing results to output file: ',file1

cnjg        open (30,file=file1,status='new')
        open (30,file=file1,status='unknown')
        write(30,'(3x,25f12.4)') 
     1      (testInten(J2),J2=1,nInten)

        do i=1,nPerc
          write (30,'(2x,f4.2,30e12.4)') i*step,
     1	  (perc(i,k),k=1,nInten)
        enddo

        write(30,'(/,2x,a4,30e12.4)') 'mean',
     1         (mean(j),j=1,nInten)
        close (30)

 1000 continue

      write (*,*) 
      write (*,*) '*** Fractile Code (42) Completed with Normal Termination ***'

      stop
      end

c -------------------------------------------------------------------------
