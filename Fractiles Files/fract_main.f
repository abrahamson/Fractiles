      program Fract_Haz45

c     This  program will compute the fractiles from a seismic hazard
c     run. This version of the program works with the seismic hazard
c     code version 45.3. The fractiles are computed based on
c     a Monte Carlo simulation.

c     compatible with Haz45.3

      implicit none
      include 'fract.h'

      integer iPer, nFlt0, nattentype, nProb, iFlt0, iFlt1, jAttenType, jFlt,
     1        i, j, k, jj, i1, i2, j2, nn100, iSample, mcSegModel, mcFType,
     2        mcWidth, mcParam, nperc, nInten, nAtten(MAX_FLT), nsite, nFlt,
     3        isite, iInten, nWidth(MAX_FLT), nParamVar(MAX_FLT,MAX_WIDTH),
     4        attentype(MAX_FLT), nGM_model(MAX_ATTENTYPE), nHazLevel
      integer iseed, nSample, nFtype(MAX_FLT), mcatten(MAX_SAMPLE), iHazLevel,
     1        f_start(MAX_FLT), f_num(MAX_FLT), nSegModel(MAX_FLT), iAmp, iPC,
     2        faultflag(MAX_FLT,MAX_SEG,MAX_FLT), iFract, nPer, jPer, iperc,
     3        PCflag(MAX_PROB)
      real tempx, step, testInten(MAX_INTEN), probAct(MAX_FLT),
     1     Haz(7,MAX_ATTEN, MAX_FLT, MAX_WIDTH, MAXPARAM, MAX_FTYPE),
     2     al_segwt(MAX_FLT), cumWt_segModel(MAX_FLT,MAX_SEG), ran1,
     3     cumWt_param(MAX_FLT,MAX_WIDTH,MAXPARAM), x, perc(100,MAX_INTEN),
     4     cumWt_Width(MAX_FLT,MAX_WIDTH), cumWt_Ftype(MAX_FLT,MAX_FTYPE)
      real Haz1(MAX_SAMPLE,MAX_INTEN), sortarray(MAX_SAMPLE), mean(MAX_INTEN),
     1     cumWt_GM(MAX_ATTEN,MAX_ATTEN), hazTotal(100), UHS(MAX_PROB,5,10),
     2     specT1(MAX_PROB), hazLevel(10), testHaz, version, xi, hermite(10),
     3     gasdev
      real*8 sum
      character*80 filein, file1


*** need to fix: treating all ftype as epistemic

      write (*,*) '**********************************'
      write (*,*) '*  Fractiles Code: Version 45.3  *'
      write (*,*) '*       Under Development        *'
      write (*,*) '*           May, 2021            *'
      write (*,*) '**********************************'

c     Open and read the run file
      write (*,*) 'Enter the input filename.'
      read (*,'(a80)') filein
      open (31,file=filein,status='old')
      read (31,*) iseed
      read (31,*) nSample
      read (31,*) nPer
      read (31,*) nHazLevel, (HazLevel(k),k=1,nHazLevel)

      read (31,'( a80)') file1
      write (*,'( a80)') file1
      open (30,file=file1,status='new')

      call CheckDim ( nSample, MAX_SAMPLE, 'MAX_SAMPLE ' )

c     Loop over each period
      do 1100 jPer=1,nPer
       read (31,*) iPer

c     Read Input File
      call RdInput (nInten, testInten, nGM_model, cumWt_GM, nAttenType,
     1              attenType, nProb, iPer, specT1(jPer), version, PCflag)

c     Read Fault Data (only if this is the first period)
      if ( jPer .eq. 1 ) then
           call Rd_Fault_Data ( version, nFlt, nFlt0,
     7     cumWt_segModel, cumWt_param, cumWt_width, probAct,
     2     nParamVar, attentype, cumwt_ftype,
     3     nFtype, nWidth, nSegModel, f_start, f_num, faultFlag, al_segwt)
      endif

c     Loop Over Number of Sites
      nsite = 1
      do 1000 iSite = 1, nSite

c       Initialize Haz1 and mean
        do i=1,nSample
          do j=1,MAX_INTEN
            Haz1(i,j) = 0.
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

c       Loop over Inten to save memory (read haz at only one intensity value at a time)
        do iInten=1,nInten

         write (*,'( 2x,''reading out1 file for iInten='',i5)') iInten

        call read_Out1 ( nFlt, Haz, nWidth, nGM_Model, attenType,
     1       nParamVar, nAtten, nFtype, iPer, nProb, iInten, nInten, PCflag)
         write (*,'( 2x,''out of read_Out1'')')

c        Monte Carlo Sampling of Hazard
         do iSample=1,nSample

          nn100 = (iSample/100)*100
          if ( nn100 .eq. iSample) then
             write (*,*) 'Looping over samples:        ', iSample, ' of ', nSample
          endif

c         Sample standard normal for the PC method
          if (PCflag(iPer) .eq. 1 ) then
            xi = gasdev(iseed)
            call calc_Hermite ( xi, hermite )
          endif

c         Sample the attenuation relation (correlated for all sources)
          do jj=1,nAttenType
           call GetRandom1 ( iseed, nGM_model(jj), cumWt_GM, jj, mcAtten(jj), MAX_ATTEN, 1 )
          enddo

c         Sample the hazard from each geometrically independent source region (or fault)
c         (a source region may include multiple subsources)
          do iFlt0 = 1, nFlt0

c           Select the seg model for this source region
            call GetRandom1 ( iseed, nSegModel(iFlt0), cumWt_segModel, iflt0, mcSegModel, MAX_FLT, 2 )

            i1 = f_start(iFlt0)
            i2 = f_start(iFlt0) + f_num(iFLt0) - 1

            do iflt1=i1,i2

c             set the attenuation type for this fault
              jAttenType = attenType(iflt1)

c             Reset the index to the faults in segment
              jFlt = iFlt1-i1+1

c             Skip this fault if not included in this seg model
              if (faultflag(iFlt0,mcSegModel,jFlt) .eq. 0 ) goto 50

c *** fix this - Now treats aleatory ftype variability as epistemic ***
c             Select the fault mechanism type
              call GetRandom1b ( iseed, nFtype(iFlt1), cumWt_Ftype, iflt1,
     1            mcFTYPE, MAX_FTYPE, MAX_FLT )

c             Select the fault width for this subsource
              call GetRandom1 ( iseed, nWidth(iFlt1), cumWt_width, iflt1, mcWidth, MAX_FLT, 3 )

c             Select the parameter variation for this subsource and fault width
              call GetRandom2 (iseed, nParamVar(iFlt1,mcWidth),
     1          cumWt_param, iFlt1, mcWidth, mcParam, MAX_FLT, MAX_WIDTH)

c             Add the sampled hazard curve to the total hazard
c             Standard Method
              if (PCflag(iPer) .eq. 0) then
                Haz1(iSample,iInten) = Haz1(iSample,iInten)
     1                    + Haz(1,mcAtten(jAttenType),iFlt1,mcWidth,mcParam,mcFtype)
     1                    * al_segwt(iflt1)
                mean(iInten)=mean(iInten)
     1	 	          + Haz(1,mcAtten(jAttenType),iFlt1,mcWidth,mcParam,mcFtype)
     1                    * al_segwt(iflt1)

c             PC method
              else if (PCflag(iPer) .eq. 1) then
                sum = 0.
                do iPC=1,7
                  sum = sum + hermite(iPC)* Haz(iPC,mcAtten(jAttenType),iFlt1,mcWidth,mcParam,mcFtype)
                enddo
                Haz1(iSample,iInten) = Haz1(iSample,iInten)
     1                    + sum * al_segwt(iflt1)
                mean(iInten)=mean(iInten)
     1	 	                 + sum * al_segwt(iflt1)
              endif

  50          continue

c           end loop over faults in this flt system
            enddo

c         end loop over flt systems
          enddo

c        end loop over number of monte carlo samples
         enddo

c       end loop over iInten
        enddo

c       Compute mean hazard as check between fractile and haz code
         do iInten = 1,nInten
           mean(iInten) = mean(iInten)/nSample
        enddo

c       Sort Haz data so that cummulative empirical distribution function
c       can be computed  (sorts in place)
c       sortarray is a working array
        do iInten = 1,nInten
           call sort(Haz1(1,iInten),sortarray,nSample)
        enddo

        nperc = 99
c       Compute fractile values
        step = 1./(nperc+1)
        do iInten=1,nInten
          do iperc = 1,nPerc
            i1 = int(iperc*step*nSample)
  	    perc(iperc,iInten) = Haz1(i1,iInten)
          enddo
	enddo

c       Write output
        write (30,'( f10.3,2x,''Period (sec)'')') specT1(jPer)
        write(30,'(3x,25f12.4)') (testInten(J2),J2=1,nInten)
        do i=1,nPerc
          write (30,'(2x,f4.2,30e12.4)') i*step, (perc(i,k),k=1,nInten)
        enddo
        write(30,'(/,2x,a4,30e12.4)') 'mean', (mean(j),j=1,nInten)
        write(30,'(/,''----------------------------------------------'',/)')

c       Compute the UHS value for this period for the 5th, 10th, 50th, 90th, 95th
        do iHazLevel=1,nHazLevel
         do iFract=1,5
          if ( iFract .eq. 1 ) i1 = 5
          if ( iFract .eq. 2 ) i1 = 10
          if ( iFract .eq. 3 ) i1 = 50
          if ( iFract .eq. 4 ) i1 = 90
          if ( iFract .eq. 5 ) i1 = 95

c         Copy sorted hazard to 1-d array
          do k=1,nInten
            hazTotal(k) = perc(i1,k)
          enddo

          testHaz = hazLevel(iHazLevel)
	  do iAmp=2,nInten

c          Check for zero values in hazard curve.
           if (hazTotal(iAmp) .eq. 0. ) then
            write (*,*) 'warning: Zero Values for hazard curve at desired haz level.'
            write (*,*) 'Setting UHS to last nonzero value (neg to indicate a problem)'
            UHS(jPer,iFract,iHazLevel) = -exp(x)
           endif

c          Interpolate the hazard curve.
           if ( hazTotal(iAmp) .lt. testHaz ) then
            x = ( alog(testHaz) - alog(hazTotal(iAmp-1)) )/
     1            ( alog(hazTotal(iAmp))-alog(hazTotal(iAmp-1)))
     2          * (alog(testInten(iAmp))-alog(testInten(iAmp-1)))
     3          + alog(testInten(iAmp-1))
            UHS(jPer,iFract,iHazLevel) = exp(x)
            goto 65
           endif
          enddo
 65       continue
         enddo
        enddo

 1000 continue
 1100 continue

c     Write UHS
      do iHazLevel=1,nHazLevel
        write (30,'( //,2x,'' UHS Fractiles for hazard level: '',e10.2)') hazLevel(iHazLevel)
        write (30,'(/,''    Period      5th         10th       50th        90th        95th'')')
        do jPer=1,nPer
          write (30,'( f10.3,5e12.4)') specT1(jPer), (UHS(jPer,k,iHazLevel),k=1,5)
        enddo
      enddo


      write (*,*)
      write (*,*) '*** Fractile Code (45) Completed with Normal Termination ***'

      stop
      end

c -------------------------------------------------------------------

      subroutine read_Out1 ( nFlt, haz, nWidth, nGM_Model, attenType,
     1           nParamVar, nAtten, nFtype, iPer, nProb, jInten, nInten, PCflag)

      implicit none
      include 'fract.h'

      integer nInten, nFlt, nAtten(1), iProb, nAtten1, nWidth1, jFlt,
     1        nParamVar(MAX_FLT,MAX_WIDTH), nWidth(MAX_FLT), nfiles,
     2        nFtype(MAX_FLT), iPer, nProb, nwr, i, j, iFlt, jProb,
     3        i1, iAtten, iWidth, iFtype, jInten, jFltWidth, nInten1,
     4        nFtype1, nParamVar1(100), nGM_Model(MAX_ATTENTYPE), iPC
      integer attenType(MAX_FLT), PCflag(MAX_PROB)
      real haz(7,MAX_ATTEN,MAX_FLT,MAX_WIDTH,MAXPARAM,MAX_FTYPE),
     1     temp(MAX_INTEN), version
      character*80 file1, dummy

      nwr = 12

c     Open output file
      if (jInten .eq. 1 ) then
        read (31,'( a80)') file1
        write (*,*) 'Opening out1 file from the hazard runs.'
        write (*,*) file1
        open (nwr,file=file1,status='old')
      else
        rewind (nwr)
      endif

C     Check for version compatibility with hazard code
      read (nwr,*) version
      if (version .ne. 45.3) then
        write (*,*) 'Incompatible version of Haz45, use Haz45.3'
        stop 99
      endif

c     Read out1 file
      do iFlt=1,nFlt
       do iWidth=1,nWidth(iFlt)
        do jProb=1,nProb

          read (nwr,*,err=200) jFlt, jFltWidth, iProb, nAtten1, nWidth1, nFtype1,
     1         nParamVar1(iWidth),nInten1

C        Check for Array Dimension of haz Array.
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
         if (iProb .gt. MAX_PROB) then
           write (*,*) 'MAX_PROB needs to be increased to ', iProb
           write (*,*) 'Change Array Parameter in FRACT.H and recompile.'
           stop 99
         endif
         if (nFtype(iFlt) .gt. MAX_FTYPE) then
           write (*,*) 'MAX_FTYPE needs to be increased to ', nFtype(iFlt)
           write (*,*) 'Change Array Parameter in FRACT.H and recompile.'
           stop 99
         endif

         do iAtten=1,nGM_model(attenType(iFlt))
          do iFtype=1,nFtype(iFlt)
           do i=1,nParamVar(iFlt,iWidth)
             if (PCflag(iProb) .eq. 0) then
               read (nwr,*,err=201) (temp(j),j=1,nInten)
c              Only keep if it is the desired spectral period
               if ( iProb .eq. iPer) then
                 haz(1,iAtten,iFlt,iWidth,i,iFtype) = temp(jInten)
               endif
             elseif (PCflag(iProb) .eq. 1) then
               do iPC=1,7
                 read (nwr,*,err=201) (temp(j),j=1,nInten)
c                only keep if it is the desired spectral period
                 if ( iProb .eq. iPer) then
                   haz(iPC,iAtten,iFlt,iWidth,i,iFtype) = temp(jInten)
                 endif
               enddo
             endif

           enddo
          enddo
	       enddo
        enddo
       enddo
      enddo

      if (jInten .eq. nInten ) close (nwr)

      return
 200  write (*,'( 2x,''Error reading iflt line in out1'',3i5)') iflt, iWidth, iProb
      backspace (nwr)
      read (nwr,'( a80)') dummy
      write (*,'( a80)') dummy
      stop 99
 201  write (*,'( 2x,''Error reading haz line in out1'',5i5 )') iflt, iWidth, iProb, iAtten, iFtype, i
      backspace (nwr)
      read (nwr,'( a80)') dummy
      write (*,'( a80)') dummy
      stop 99
      end

c -------------------------------------------------------------------------

      subroutine calc_Hermite ( xi, hermite)

      implicit none

      real xi, hermite(10), xi2, xi3, xi4, xi5, xi6

      xi2 = xi**2
      xi3 = xi2*xi
      xi4 = xi3*xi
      xi5 = xi4*xi
      xi6 = xi5*xi

c     Evaluate the hermit polynomials
      hermite(1) = 1.
      hermite(2) = xi
      hermite(3) = xi2 - 1
      hermite(4) = xi3 - 3.*xi
      hermite(5) = xi4 - 6.*xi2 + 3
      hermite(6) = xi5 - 10.*xi3 + 15.*xi
      hermite(7) = xi6 - 15.*xi4 + 45*xi2 - 15.

      return
      end
