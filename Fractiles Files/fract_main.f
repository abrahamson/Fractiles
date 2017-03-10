      program Fract_Haz45

c     This  program will compute the fractiles from a seismic hazard
c     run. This version of the program works with the seismic hazard
c     code version 45.2. The fractiles are computed based on 
c     a Monte Carlo simulation.  
 
c     compatible with Haz45.2

      implicit none
      include 'fract.h'
      
      integer iPer, nFlt0, nattentype, nProb, iFlt0, iFlt1, jAttenType, jFlt,
     1        i, j, k, jj, i1, i2, j2, nn100, iSample, mcSegModel, 
     2        nperc, nInten, nAtten(MAX_FLT), nsite, nFlt,
     3        isite, iInten, nWidth(MAX_FLT), nParamVar(MAX_FLT,MAX_WIDTH),
     4        attentype(MAX_FLT), nGM_model(MAX_ATTENTYPE), nHazLevel   
      integer iseed, nSample, mcatten(MAX_SAMPLE), iHazLevel,
     1        f_start(MAX_FLT), f_num(MAX_FLT), nSegModel(MAX_FLT), iAmp,
     2        faultflag(MAX_FLT,MAX_SEG,MAX_FLT), iFract, nPer, jPer, iperc   
      real tempx, step, testInten(MAX_INTEN), probAct(MAX_FLT),
     1     Haz(MAX_ATTEN, MAX_FLT, MAX_WIDTH, MAXPARAM, MAX_FTYPE), 
     2     al_segwt(MAX_FLT), cumWt_segModel(MAX_FLT,MAX_SEG), ran1,
     3     x, perc(100,MAX_INTEN)
      real Haz1(MAX_SAMPLE,MAX_INTEN), sortarray(MAX_SAMPLE), mean(MAX_INTEN),
     1     cumWt_GM(MAX_ATTEN,MAX_ATTEN), hazTotal(100), UHS(MAX_PROB,5,10), 
     2     specT1(MAX_PROB), hazLevel(10), testHaz, version      
      character*80 filein, file1, fname(MAX_FLT)
      
      integer iCorrFlag, nCorr1, iFlt10, jCOrr,
     1        corrNode(100), corrNFlt(100), corrFlt(100,MAX_FLT), iCOrr

      integer nFtype(MAX_FLT,MAXPARAM)

      real segwt(MAX_FLT,MAX_SEG)
      integer nNode_SSC, iFLt2, iNode, iThick, iFtype
      integer jBR(MAX_FLT,MAXPARAM), mcBR, iRate, nRate, iParam
      integer iWidth, iFM, isum, iFtypeIndex
      integer iBR, nBR_all(MAX_FLT,MAXPARAM,10)
      real wt_cum(MAXPARAM), wt_cum_all(MAX_FLT,12,10,10),
     1     wt_ftype(MAX_FLT,5,5)
      integer iSave, iFlag, firstCorr(100), kflt2(12,100)

      
*** need to fix: treating all ftype as epistemic

c      integer TreeIndex(MAX_FLT,MAX_DIP,MAX_WIDTH,MAX_N2,MAX_N2,MAX_N2,MAX_N2)      

c      open (33,file='debug.out')
    
*** need to fix: treating all ftype as epistemic
        
      write (*,*) '*************************'
      write (*,*) '*   Fractile Code for   *'
      write (*,*) '*       HAZ 45.2        *'
      write (*,*) '*       Mar 2017        *'
      write (*,*) '*************************'

c     Open and read the run file
      write (*,*) 'Enter the input filename.'      
      read (*,'(a80)') filein
      open (31,file=filein,status='old')
      read (31,*) iseed
      read (31,*) nSample
      read (31,*) nPer
      read (31,*) nHazLevel, (HazLevel(k),k=1,nHazLevel)

c     Read correlation
      read (31,*) nCorr1
      if ( nCorr1 .gt. 100 ) then
        write (*,'( 2x,''Maximum number of correlations is 100'')')
        stop 99
      endif
      if ( nCorr1 .ne. 0 ) then
        do i=1,nCorr1
          read (31,*) corrNode(i), corrNFlt(i)
          read (31,*) (corrFlt(i,j),j=1,corrNFlt(i))
        enddo
      endif

      read (31,'( a80)') file1
      write (*,'( a80)') file1
      open (30,file=file1,status='unknown')
      rewind (30)

      call CheckDim ( nSample, MAX_SAMPLE, 'MAX_SAMPLE ' )

c     Loop over each period
      do 1100 jPer=1,nPer
       read (31,*) iPer

c     Read Input File
      call RdInput ( nInten,  testInten, nGM_model, cumWt_GM, 
     5     nattentype, attenType, nProb, iPer, specT1(jPer), version)
     
c     Read Fault Data (only the if this is the first period)
      if ( jPer .eq. 1 ) then
c       Read fault file
           call Rd_Fault_Data ( version, nFlt, nFlt0, probAct, attentype,
     1          nSegModel, cumWt_segModel, nFtype, wt_ftype, f_start, f_num,
     2          faultFlag, al_segwt, nBR_all, wt_cum_all, fname )
        do i=1,nFlt
          isum = 0
          do ithick=1,nBR_all(i,1,1)
            isum = isum + nBR_all(i,2,ithick)
          enddo
          nWidth(i) = isum
        enddo
      write (*,'( 2x,'' out of read flt'',)')

      endif
      
c     SSC branches: 1=crustal thick, 2 = Dip, , 3=ftype, 4=magpdf, 5=refmag,
c                   6=RateType, 7=SR, 8=a-value, 9=paleo, 
c                   10=Moment, 11=b_value flt, 12=segModel
      nNode_SSC = 12

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
     1       natten, iPer, nProb, iInten, nInten )
         write (*,'( 2x,''out of read_Out1'')')  

c        Monte Carlo Sampling of Hazard
         do iSample=1,nSample

c         initialize correlation flags for this realization
          do icorr=1,nCorr1
            firstCorr(iCorr) = 0
          enddo         

          nn100 = (iSample/100)*100
          if ( nn100 .eq. iSample) then
             write (*,*) 'Looping over samples:        ', iSample, ' of ', nSample
          endif

c         Sample the attenuation relation (correlated for all sources)
          do jj=1,nAttenType
           call GetRandom1 ( iseed, nGM_model(jj), cumWt_GM, jj, mcAtten(jj), MAX_ATTEN, 1 ) 
          enddo
            
c         Sample the hazard from each geometrically independent source region (or fault)
c         (a source region may include multiple subsources) 
          iFlt2 = 0
          do iFlt0 = 1, nFlt0

c           Select the seg model for this source region
            call GetRandom1 ( iseed, nSegModel(iFlt0), cumWt_segModel, iFlt0, 
     1              mcSegModel, MAX_FLT, 2  )

            i1 = f_start(iFlt0)
            i2 = f_start(iFlt0) + f_num(iFLt0) - 1
            
            do iflt1=i1,i2
              iFlt2 = iFlt2 + 1

c             set the attenuation type for this fault
              jAttenType = attenType(iflt1)

c             Reset the index to the faults in segment
              jFlt = iFlt1-i1+1

c             Skip this fault if not included in this seg model
              if (faultflag(iFlt0,mcSegModel,jFlt) .eq. 0 ) goto 50

c             Sample each of the SSC nodes (skip seg model node) 
c             thickness is the first node           
              do iNode=1,11 

c              Check if this node has branches
               if ( nBR_all(iFlt2,iNode,1) .eq. 0 ) goto 100

c              Check if this is a correlated node
c              if so, check if this is the first sample of this correlated node
               iCorrFlag = 0
               do i=1,nCorr1
                 if ( iNode .eq. corrNode(i) ) then
                   do iflt10=1,corrNFlt(i)
                    if ( iFlt2 .eq. corrFlt(i,iflt10)) then
                      jCorr = i
                      if ( firstCorr(i) .eq. 0 ) then
                        iCorrFlag = 1
                        firstCorr(i) = 1                       
                      else
                        iCorrFlag = 2
                      endif
c                      write (*,'( 2x,''correlated node, flt, corrIndex:'',4i5)') iNode, iFlt2, i
                    endif
                   enddo
                 endif
               enddo
c               write (*,'( 5i5)') iNode, iFlt2,iCorrFlag
               
               if ( iCorrFlag .eq. 0 .or. iCorrFlag .eq. 1) then
c                copy branch weights for this node to single dimension array 
                 if (iNode .eq. 1 ) then
                   ithick=1
                 else
                   ithick = jBR(iFlt2,1)
                 endif
                 do iBR=1,nBR_all(iFlt2,iNode,ithick)
                   wt_cum(iBR) = wt_cum_all(iFlt2,iNode,ithick,iBR)
                 enddo              
                 call GetRandom0 ( iseed, nBR_all(iFlt2,iNode,ithick), wt_cum, mcBR, 
     1                MAXPARAM, iNode, iFlt2 )
                 jBR(iFlt2,iNode) = mcBR
 
c                if this is the first sample of a correlated bramch, then save for later use
                 if ( iCorrFlag .eq. 1 ) then
                   kflt2(iNode,jCorr) = iFlt2
                 endif
               else
c                Correlated, first check if number of branches match
                 if ( nBR_all(iFlt2,iNode,ithick) .ne. nBR_all(kFlt2(iNode,jCorr),iNode,ithick) ) then
                   write (*,'( 2x,''Correlated nodes do not have the same number of branches'')')
                   write (*,'( 2x,a60,i5)') fname(iFlt2), nBR_all(iFlt2,iNode,ithick)              
                   write (*,'( 2x,a60,i5)') fname(kFlt2(iNode,jCorr)), nBR_all(kFlt2(iNode,jCorr),iNode,ithick)  
                   stop 99
                 endif            
c                Set to previous selected branch 
                 jBR(iFlt2,iNode) = jBR(kFlt2(iNode,jCorr),iNode)
c                 write (*,'( 6i5)') iflt2, icorrFlag, kflt2(iNode,jCorr), jCorr, iNode, 
c     1               jBR(iFlt2,iNode)
c                 pause 'test correlated branch'
                 
               endif
 100           continue

              enddo

c             Here, all branches are set for this fault        
c             Set the parameter index as used in the hazard calc
              ithick = jBR(iFlt2,1)
              if (jBR(iFlt2,6) .eq. 1 ) then
               iRate = jBR(iFlt2,7)
              elseif (jBR(iFlt2,6) .eq. 2 ) then
               iRate = nBR_all(iFlt2,7,1) + jBR(iFlt2,8)
              elseif (jBR(iFlt2,6) .eq. 3 ) then
               iRate = nBR_all(iFlt2,7,1)+ nBR_all(iFlt2,8,1) + jBR(iFlt2,9)
              elseif (jBR(iFlt2,6) .eq. 4 ) then
               iRate = nBR_all(iFlt2,7,1)+ nBR_all(iFlt2,8,1) + nBR_all(iFlt2,9,1) + jBR(iFlt2,10)
              endif
              nRate = nBR_all(iFlt2,7,1)+ nBR_all(iFlt2,8,1) + nBR_all(iFlt2,9,1) + nBR_all(iFlt2,9,10)           
c             set iParam = (iMagRecur-1)*Nrate*n_b*nRefMag
c                          + (iRate-1)*N_b*NrefMag + (ib-1)*nRefMag + iRefMag
              iParam = (jBR(iFlt2,4)-1) * nRate * nBR_all(iflt2,11,iThick)
     1                           * nBR_all(iflt2,5,iThick)
     1            + (iRate-1) *  nBR_all(iflt2,11,ithick) * nBR_all(iflt2,5,ithick)
     2            + (jBR(iFlt2,11)-1) * nBR_all(iflt2,5,ithick)
     3            + jBR(iFlt2,5)      
     
c              iWidth = (ithick-1) * nDip + iDip
              iWidth = (ithick-1) *nBR_all(iflt2,2,ithick) + jBR(iFlt2,2)
              
c             rename ftype model selected
              iFM = jBR(iFlt2,3)            
              isum = 0
              do i2=1,iFM-1
                isum = isum + nFtype(iflt2,i2)
              enddo
              
c              write (*,'( 12i5)') mcSegModel, (jBR(iFlt2,k),k=1,11)
c              write (*,'( 12i5)') nSegModel(iflt0), (nBR_all(iflt2,k,iThick),k=1,11)
c              write (*,'( 4i5)') iRate, iParam, iWidth, isum

c             sum over aleatory ftype 
              do iFtype=1, nFtype(iflt2,iFM)
               iFtypeIndex = isum + iFtype
           
c              Add the sampled hazard curve to the total hazard
               Haz1(iSample,iInten) = Haz1(iSample,iInten) 
     1                    + Haz(mcAtten(jAttenType),iFlt2,iWidth,iParam,iFtypeIndex)
     1                    * al_segwt(iflt2) * wt_ftype(iflt2,iFM,iFtype)
     
               mean(iInten)=mean(iInten)
     1	 	          + Haz(mcAtten(jAttenType),iFlt2,iWidth,iParam,iFtypeIndex)
     1                    * al_segwt(iflt2) * wt_ftype(iflt2,iFM,iFtype)
              enddo

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
     1            nAtten, iPer, nProb, jInten, nInten)

      implicit none
      include 'fract.h'

      integer nInten, nFlt, nAtten(1), iProb, nAtten1, nWidth1, jFlt,
     1        nWidth(MAX_FLT), nfiles,
     2        iPer, nProb, nwr, i, j, iFlt, jProb,  
     3        i1, iAtten, iWidth, iFtype, jInten, jFltWidth, nInten1,
     4        nFtype1, nParamVar1(100), nGM_Model(MAX_ATTENTYPE)
      integer attenType(MAX_FLT)
      real haz(MAX_ATTEN,MAX_FLT,MAX_WIDTH,MAXPARAM,MAX_FTYPE),
     1     temp(MAX_INTEN), version
      character*80 file1, dummy

      nwr = 12
      nfiles = 1

C     Loop over the number of files.
      do 100 i1=1,nfiles
                     
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
         if (version .ne. 45.2) then
         write (*,*) 'Incompatible version of Haz45 out1 file, use Haz45.2'
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
         if (nWidth1 .gt. MAX_WIDTH) then
           write (*,*) 'MAX_WIDTH needs to be increased to ', nWidth(iFlt)
           write (*,*) 'Change Array Parameter in FRACT.H and recompile.'
           stop 99
         endif
         if (iProb .gt. MAX_PROB) then
           write (*,*) 'MAX_PROB needs to be increased to ', iProb
           write (*,*) 'Change Array Parameter in FRACT.H and recompile.'
           stop 99
         endif
         if (nFtype1 .gt. MAX_FTYPE) then
           write (*,*) 'MAX_FTYPE needs to be increased to ', nFtype1
           write (*,*) 'Change Array Parameter in FRACT.H and recompile.'
           stop 99
         endif

         do iAtten=1,nGM_model(attenType(iFlt))
          do iFtype=1,nFtype1
           do i=1,nParamVar1(iWidth)
 	    read (nwr,*,err=201) (temp(j),j=1,nInten)

c           Only keep if it is the desired spectral period
	    if ( iProb .eq. iPer) then
	      haz(iAtten,iFlt,iWidth,i,iFtype) = temp(jInten)
            endif
           enddo
          enddo
	 enddo
        enddo
       enddo
      enddo

  100 continue

      if (jInten .eq. nInten ) close (nwr)

      return
 200  write (*,'( 2x,''Error reading iflt line in out1'',3i5)') iflt, iWidth, iProb
      backspace (nwr)
      read (nwr,'( a80)') dummy
      write (*,'( a80)') dummy
      stop 99   
 201  write (*,'( 2x,''Error reading haz line in out1'',5i5 )') iflt, iWidth,
     1          iProb, iAtten, iFtype, i
      backspace (nwr)
      read (nwr,'( a80)') dummy
      write (*,'( a80)') dummy
      stop 99
      end

c -------------------------------------------------------------------------
