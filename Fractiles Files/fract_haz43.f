      program Fract_Haz43

c     This program will compute the fractiles from a seismic hazard
c     run. This version of the program works with the seismic hazard
c     code version 43 or later. The fractiles are computed based on 
c     a Monte Carlo simulation. This version will be able to read fault
c     files in which synchronous rupture is modeled. 
 
c     compatible with Haz43
c     Last modified: 2/11

      include 'fract.h'
c      implicit none
      
      real risk(MAX_INTEN,MAX_ATTEN, MAX_FLT, MAX_WIDTH, MAXPARAM,MAX_FTYPE), 
     1     testInten(MAX_INTEN)
      real al_segwt(MAX_FLT), cumwt_probact(2)
      real risk_noMix(MAX_INTEN,MAX_ATTEN, MAX_FLT, MAX_WIDTH, MAXPARAM,MAX_FTYPE)
       real probAct(MAX_FLT)
c      real cum_wt(5)
      integer nInten, jcalc(MAX_ATTENTYPE,MAX_ATTEN), nAtten(MAX_FLT), nsite, nFlt,
     1        isite, iflt, 
     1        iAtten, iInten, nWidth(MAX_FLT)
      integer iParam, nParamVar(MAX_FLT,MAX_WIDTH)
      integer  iFltWidth, iX(MAX_FILES)
      character*80 filein, file1, fileinPSHA
      integer nFlt1(MAX_FLT), nFlt2, attentype(MAX_FLT), nGM_model(MAX_ATTENTYPE)
      real cumWt_segModel(MAX_FLT,MAX_FLT),
     1     cumWt_param(MAX_FLT,MAX_WIDTH,MAXPARAM),
     1     cumWt_Width(MAX_FLT,MAX_WIDTH), cumWt_Ftype(MAX_FLT,MAX_FTYPE)
      real ran1, risk1(MAX_SAMPLE,MAX_INTEN)
      real sortarray(MAX_SAMPLE), attenWt(MAX_PROB)
      real cumWt_GM(MAX_ATTEN,MAX_ATTEN)
      integer Fflag 
      integer iperc, nfiles, n_Dip(MAX_FLT)
      real perc(100,MAX_INTEN),mean(MAX_INTEN) 
      integer nProb1(MAX_ATTENTYPE), jProb(MAX_PROB), kProb, nDip(MAX_FLT)
      integer iseed, nSample, nFtype(MAX_FLT)
      real timeDepfactor(10), timeDepwt(10), cumwt_TDF(10), cumwt_mix(2)
      integer mcatten(MAX_SAMPLE), iMix(MAX_SAMPLE), c
      integer f_start(MAX_FLT), f_num(MAX_FLT), nSegModel(MAX_FLT)
      integer faultflag(MAX_FLT,MAX_SEG,MAX_FLT)
      real haz_other(100), haz_SA(100)
      real timeDepfactor1(10), timeDepwt1(10), cumwt_TDF1(10)

*** need to fix: treating all ftype as epistemic

c      integer TreeIndex(MAX_FLT,MAX_DIP,MAX_WIDTH,MAX_N2,MAX_N2,MAX_N2,MAX_N2)      

c      open (33,file='debug.out')
    
      write (*,*) '*************************'
      write (*,*) '* Fractile Code for use *'
      write (*,*) '*  Hazard_v43.exe code  *'
      write (*,*) '*    March, 2015, NAA     *'
      write (*,*) '*************************'

      write (*,*) 'Enter the input filename.'
      
      read (*,'(a80)') filein
      open (31,file=filein,status='old')

      read (31,*) iseed
      read (31,*) nSample

      read (31,*) nAttentype
      read (31,*) (nProb1(j), j=1,nattentype)

      read (31,*) nSet

      read (31,*) nn
      read (31,*) (Haz_other(k),k=1,nn)
      read (31,*) (Haz_SA(k),k=1,nn)
          read (31,*) nTimeDep1
          read (31,*) (timedepFactor1(k),k=1,nTimeDep1)
          read (31,*) ( timeDepWt1(k), k=1,nTimeDep1)
          do k=1,nTimeDep1
            if (k .eq. 1 ) then
              cumwt_TDF1(k) = timeDepwt1(k)
            else
              cumwt_TDF1(k) = cumwt_TDF1(k-1) + timeDepwt1(k)
            endif
          enddo

      call CheckDim ( nSample, MAX_SAMPLE, 'MAX_SAMPLE ' )

C     Write out the intensity values if requested.
c       if (outflag .eq. 1) then
c          write (17,*) nInten, (testInten(j2), j2=1,nInten) 
c       endif

c     Loop Over Number of Sites
      nsite = 1
      do 1000 iSite = 1, nSite

c       Initialize risk1 and mean to the other faults 
        do i=1,nSample
          do j=1, nn
            risk1(i,j) = haz_other(j)
          enddo
        enddo

c       Initialize mean
        do i=1, nn
          mean(i) = haz_other(i)*nSample
        enddo

c       Initialize random number generator
        do i=1,500
          tempx = Ran1 (iseed)
        enddo

        write (*,'(2x,a18,i4,a5,2x,i4)') 'Looping over site ', isite,'  of ', nsite

        do iSet=1,nSet

c         Read Input File
          call RdInput ( nFlt, nFlt0, f_start, f_num, faultFlag, 
     1     nInten,  testInten, lgTestInten, probAct, al_segWt, 
     3     cumWt_SegModel, cumWt_Width, cumWt_param, cumWt_ftype, cumWt_GM, 
     4     nSegModel,      nWidth,      nParamVar,   nFtype,      nGM_model )

          write (*,'( i5)') nGM_model(1)
          read (31,*) ifix

          write (*,'( 10f10.3)') (probAct(k),k=1,nFlt)

       do iflt0=1,nflt0
        write (*,'( 2i5)') iflt0, nsegModel(iflt0)
        do i=1, nSegModel(iFlt0)
          write (*,'( 20i5)') (faultFlag(iFlt0,i,k),k=1,f_num(iflt0))
        enddo
       enddo
       


          read (31,*) nTimeDep
          read (31,*) (timedepFactor(k),k=1,nTimeDep)
          read (31,*) ( timeDepWt(k), k=1,nTimeDep)
          do k=1,nTimeDep
            if (k .eq. 1 ) then
              cumwt_TDF(k) = timeDepwt(k)
            else
              cumwt_TDF(k) = cumwt_TDF(k-1) + timeDepwt(k)
            endif
          enddo

        write (*,'( 2x,''reading logictree file'')')
        call read_logicRisk ( isite, nSite, nFlt, nfiles, ix, risk, nParamVar, natten, nFtype )
        call read_logicRisk ( isite, nSite, nFlt, nfiles, ix, risk_nomix, nParamVar, natten, nFtype )

c       write (*,*) 'Done reading logic Risk Data.'

c       set fixed mixture weights
        cumwt_mix(1) = 0.8
        cumwt_mix(2) = 1.0

c       Monte Carlo Sampling of Hazard
        do iSample=1,nSample

           nn100 = (iSample/100)*100
           if ( nn100 .eq. iSample) then
             write (*,*) 'Looping over samples:        ', iSample, ' of ', nSample
           endif

c          Sample the attenuation relation (correlated for all sources)
           if ( iSet .eq. 1 ) then
             call GetRandom1 ( iseed, nGM_model(1), cumWt_GM, 1, mcAtten(iSample), MAX_ATTEN ) 
             call GetRandom0 ( iseed, 2, cumWt_Mix, iMix(iSample))

             call GetRandom0 ( iseed, nTimeDep1, cumWt_TDF1, iTDF1 ) 
             do j=1, nn
               risk1(iSample,j) =  risk1(iSample,j) + Haz_SA(j)*timedepFactor1(iTDF1)
             enddo
            
           endif
             write (44,'( ''GM model'', 2x, 10i5)') isample, mcAtten(iSample), iMix(iSample), iTDF1,
     1           nGM_model(1)


c         Sample the time dep factr
          call GetRandom0 ( iseed, nTimeDep, cumWt_TDF, iTDF)
c          write (*,'( 2i5, 2x,''iTDF, nflt0'')') iTDF, nFlt0

c         Sample the hazard from each geometrically independent source region (or fault)
c         (a source region may include multiple subsources)
 
          do iFlt0 = 1, nFlt0



c           Select the seg model for this source region
            call GetRandom1 ( iseed, nSegModel(iFlt0), cumWt_segModel, iflt0, mcSegModel, MAX_FLT )
c           write (*,'( 3i5)') nSegModel(iFlt0), iflt0, mcSegModel

            i1 = f_start(iFlt0)
            i2 = f_start(iFlt0) + f_num(iFLt0) - 1
c            write (*,'( 3i5)') mcSegModel , i1, i2

            do iflt1=i1,i2

c             reset the index to the faults in segment
              jFlt = iFlt1-i1+1
c              write (*,'( 2x''fltflag'', 10i5)') iFlt0, mcSegMOdel, jFlt, faultflag(iFlt0,mcSegModel,jFlt)

c             Skip this fault if not included in this seg model
              if (faultflag(iFlt0,mcSegModel,jFlt) .eq. 0. ) goto 50
c                write (*,'( 3i5, 2x,''isample, mcSegmodel, iflt1'' )') isample, mcSegmodel, iflt1

c              write (*,'( i5)') nFtype(iFlt1)
c              write (*,'( 10f10.4)') (cumWt_Ftype(iflt1, ii),ii=1,nFtype(iFLt1) )
c             Select the fault mechanism type 
c              call GetRandom1b ( iseed, nFtype(iFlt1), cumWt_Ftype, iflt1, mcFTYPE, MAX_FTYPE, MAX_FLT )
c   * not working, all set to 1 for now.
               mcFType = 1

c             Select the fault width for this subsource 
              call GetRandom1 ( iseed, nWidth(iFlt1), cumWt_width, iflt1, mcWidth, MAX_FLT )

c             Select the parameter variation for this subsource and fault width
              call GetRandom2 (iseed, nParamVar(iFlt1,mcWidth), 
     1          cumWt_param, iFlt1, mcWidth, mcParam, MAX_FLT, MAX_WIDTH)

c              write (*,'( 6i5)') mcAtten(iSample), mcFtype, mcWidth, mcParam, nParamVar(iFlt1,mcWidth)

c      Temp fix for IHEB
c       if ( ifix .eq. 1 ) then
c         if (iflt1 .le. 18 ) then
c           wt_fix = 0.7
c         else
c           wt_fix = 0.3
c         endif
c       else
c          wt_fix = 1.
c       endif

              If ( iSample .eq. 1 ) write (*,'( i5,f10.4, 2x,''alsegwt'')') iflt1, al_segwt(iFlt1)
c             Add the sampled hazard curve to the total hazard
              do iInten=1,nInten
                if ( imix(iSample) .eq. 1 ) then
                  risk1(iSample,iInten) = risk1(iSample,iInten) 
     1                    + risk(iInten,mcAtten(iSample),iFlt1,mcWidth,mcParam,mcFtype)
     1                    * timeDepFactor(iTDF)  * al_segwt(iflt1)
                  mean(iInten)=mean(iInten)
     1	 	          +  risk(iInten,mcAtten(iSample),iFlt1,mcWidth,mcParam,mcFtype)
     1                    * timeDepFactor(iTDF) * al_segwt(iflt1)
                else
                  risk1(iSample,iInten) = risk1(iSample,iInten) 
     1                    + risk_nomix(iInten,mcAtten(iSample),iFlt1,mcWidth,mcParam,mcFtype)
     1                    * timeDepFactor(iTDF) * al_segwt(iflt1)
                  mean(iInten)=mean(iInten)
     1	 	          +  risk_nomix(iInten,mcAtten(iSample),iFlt1,mcWidth,mcParam,mcFtype)
     1                    * timeDepFactor(iTDF) * al_segwt(iflt1)
                endif
c         if ( iSample .lt. 5 ) write (*,'( i5,12e12.4)') isample, (risk1(iSample,jj),jj=1,nInten)

              enddo
  50          continue
            enddo
  60          continue
         enddo

c        write (*,*) 'Done with sample: ', iSample
        enddo
        enddo

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

        nperc = 99
c       Compute fractile values
        step = 1./(nperc+1)
        do iInten=1,nInten
          do iperc = 1,nPerc
            i1 = int(iperc*step*nSample)
	  perc(iperc,iInten) = risk1(i1,iInten)
	enddo
        enddo

c       Write output
        read (31,'( a80)') file1
        write(*,'( a32,a80)') 'Writing results to output file: ',file1

cnjg        open (30,file=file1,status='new')
        open (30,file=file1,status='new')
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

c -------------------------------------------------------------------

      subroutine read_logicRisk ( isite, nSite, nFlt, nfiles, ix, risk, nParamVar, nAtten, nFtype)

      include 'fract.h'

      real risk(MAX_INTEN, MAX_ATTEN, MAX_FLT, MAX_WIDTH, MAXPARAM,MAX_FTYPE)
      real temp(MAX_INTEN)
      integer isite, nInten, nFlt,  nSite, nAtten(1),ix(MAX_FILES), iProb,
     1        nParamVar(MAX_FLT,MAX_WIDTH), nWidth(MAX_FLT), nfiles, iA1, ntotal
      integer nFtype(MAX_FLT)
      character*80 fName, file1
      nwr = 12

      nfiles = 1

C     Loop over the number of files.
      do 100 i1=1,nfiles
                     
c     Open output file
      read (31,'( a80)') file1

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

        write (*,'( 2i5)') iFlt, nFlt 

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

c -------------------------------------------------------------------------
