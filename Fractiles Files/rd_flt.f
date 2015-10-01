      subroutine Rd_Fault_Data (nFltTotal,nFlt0,
     1     cumWt_flt1, cumWt_param, cumWt_width, probAct,
     2     nParamVar, AttenType, cumwt_Ftype, 
     3     nFtype, nWidth, nSegModel, f_start, f_num, faultFlag, al_Segwt )

      include 'fract.h'
      integer ndip(MAX_FLT), nWidth(MAX_FLT)
      integer  nParamVar(MAX_FLT,MAX_WIDTH), nRT(MAXPARAM)
      integer nMaxMag(MAX_FLT,MAXPARAM), nSlipRate(MAXPARAM), Attentype(MAX_FLT)
      real segwt(MAX_FLT,MAX_FLT), wt, al_segWt(MAX_FLT)
      real minmag, magstep, minDepth, coef_area(2), coef_width(2)
      real sliprateWt(MAXPARAM,MAXPARAM),  bValueWt(MAXPARAM),
     1     magRecurWt(MAXPARAM), faultWidthWt(MAXPARAM),
     2     maxMagWt(MAXPARAM,MAXPARAM), rt_wt(10,10),
     3     dipWt(MAXPARAM)
      integer nFlt, nFlt0, nFlt1(1), nFlt2
      real cumWt_flt1(MAX_FLT,MAX_FLT),cumWt_param(MAX_FLT,MAX_WIDTH,MAXPARAM),
     1     cumWt_width(MAX_FLT,1),cumwt_Ftype(MAX_FLT,MAX_FTYPE)
      real magRecur1(MAXPARAM)
      real x(MAXPARAM), x2(MAXPARAM,MAXPARAM)
      real ProbAct(1), topdepth
      character*80 fName1, fName
      character*1 tempname
      integer nFm, nFtypeModels(MAX_FLT), nFtype1(MAX_FLT)
      real ftmodelwt(MAX_FLT), ftype1(MAX_FLT,MAX_FLT), Ftype_wt1(MAX_FLT,MAX_FLT)
      real rateType1(MAXPARAM), width_wt(MAX_FLT)

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
      integer f_start(1), f_num(1) , nSegModel(1)

C     Input Fault Parameters
      read (10,*) iCoor
      read (10,*) NFLT0

      call CheckDim ( NFLT0, MAX_FLT, 'MAX_FLT' )
  
      iflt = 0

C.....Loop over each fault in the source file....
      DO iFlt0=1,NFLT0
        read (10,'( a80)') fName1
        read (10,*) probAct0

c       Read number of segmentation models for this fault system       
        read (10,*) nSegModel(iFlt0)
        read (10,*) (segWt(iFlt0,k),k=1,nSegModel(iFlt0))

C       Set up cum wt for segment models for this fault system
        do k=1,nSegModel(iFlt0)
          if ( k .eq. 1) then
            cumWt_flt1(iFlt0,1) = segWt(iFlt0,1)
          else
            cumWt_flt1(iFlt0,k) = segWt(iFlt0,k) + cumWt_flt1(iFlt0,k-1)
          endif
        enddo

c       Read total number of fault segments for this fault system
        read (10,*) nFlt2
        do i=1, nSegModel(iFlt0)
         read (10,*) (faultFlag(iFlt0,i,k),k=1,nFlt2)
        enddo

c       Set the index for the first fault in this fault system and the number
c       This allows us to find the right fault from the large list
        f_start(iFlt0) = iFlt + 1
        f_num(iFlt0) = nFlt2 

C.......Loop over number of individual fault segments....                        
        do iflt2=1,nflt2

          iFlt = iFlt + 1
          call CheckDim ( iflt, MAX_FLT, 'MAX_FLT   ' )

c         Read past name of this segment
          read(10,'( a80)') fname
          read (10,*) isourceType, attenType(iFlt), sampleStep, directflag, synflag

c         Read past the synchronous Rupture parameters
          if (synflag .gt. 0) then
            read (10,*) nsyn, attensyn
            do insyn=1,nsyn
              read (10,*) magsyn, rupsyn, jbsyn, seismosyn, hwsyn,ftypesyn, hyposyn, wtsyn
            enddo
          endif

c         Read aleatory segmentation wts
          read (10,*) al_segWt(iFlt)

c         Check for standard fault source or areal source
          if ( isourceType .eq. 1 .or. isourceType .eq. 2) then
            read (10,*) dip1, top
            read(10,*) nfp     
            do ipt=1,nfp
              read (10,*) fLong, fLat
            enddo
          endif


c         Check for grid source (w/o depth)
          if ( isourceType .eq. 3 .or. isourceType .eq. 7 ) then
            read (10,*)  dip1, top
c           Read past grid filename...
            read (10,*)
          endif

c         Check for grid source (w/ depth)
          if ( isourceType .eq. 4 ) then
            read (10,*)  dip1
c           Read over grid filename...
            read (10,*)
          endif

c         Check for custom fault source
          if ( isourceType .eq. 5) then
            read(10,*) nDownDip, nfp
c........   Only read in the first downdip case - rest is not needed...
            do ipt=1,nfp
              read (10,*) fLong, fLat, fZ
            enddo
          endif

c         Read dip Variation
          if ( isourceType .ne. 5 ) then
            read (10,*) n_Dip(iflt)
            read (10,*) (x(i),i=1,n_Dip(iflt))
            read (10,*) (dipWt(i),i=1,n_Dip(iflt))
          else
            n_Dip(iflt) = 1
            x(1) = 0.
            dipWt(1) = 1.
          endif

c         Read b-values (not for activity rate cases)
          read (10,*) n_bValue
          call CheckDim ( n_bValue, MAX_N1, 'MAX_N1    ' )
          if ( n_bValue .gt. 0 ) then
            read (10,*) (x(i),i=1,n_bValue)
            read (10,*) (bValueWt(i),i=1,n_bValue)
          endif
                
c         Read activity rate - b-value pairs
          read (10,*) nActRate 
          if ( nActRate .ne. 0 ) then
            do ii=1,nActRate
              read (10,*) bValue2(ii), actRate(ii), actRateWt(ii)
            enddo
          endif

c         Read weights for rate methods
          read (10,*) wt_srBranch, wt_ActRateBranch, wt_recIntBranch, wt_MoRateBranch
          sum = wt_srBranch + wt_ActRateBranch + wt_recIntBranch + wt_MoRateBranch
          if ( sum .lt. 0.999 .or. sum .gt. 1.001 ) then
              write (*,'( 2x,''rate method weights do not sum to unity for fault, '',a30)') fName
              stop 99
          endif

c         Read slip-rates
          read (10,*) nSR
          if ( nSR .gt. 0 ) then
            read (10,*) (x(k),k=1,nSR)
            read (10,*) (wt_sr(k),k=1,nSR)
          endif

c         Read recurrence intervals
          read (10,*) nRecInt
          if ( nRecInt .gt. 0 ) then
            read (10,*) (x(k),k=1,nRecInt)
            read (10,*) (wt_recInt(k),k=1,nRecInt)
          endif

c         Read moment-rates
          read (10,*) nMoRate
          if ( nMoRate .gt. 0 ) then
            read (10,*) (x(k),k=1,nMoRate)
            read (10,*) (x(k),k=1,nMoRate)
            read (10,*) (wt_MoRate(k),k=1,nMoRate)
          endif
              
c         Load into single array called "rate_param"
          nRate = nSR + nActRate + nRecInt + nMoRate
          call CheckDim ( nRate, MAX_N1, 'MAX_N1    ' )
          do k=1,nSR
            rateWt1(k) = wt_sr(k)*wt_srbranch              
            rateType1(k) = 1
          enddo
          do k=1,nActRate
            rateWt1(k+nSR) = actRateWt(k) * wt_actRateBranch             
            rateType1(k+nSR) = 2
          enddo
          do k=1,nRecInt
            rateWt1(k+nSR+nActRate) = wt_recInt(k) *wt_recIntBranch              
            rateType1(k+nSR+nActRate) = 3
          enddo
          do k=1,nMoRate
            rateWt1(k+nSR+nActRate+nRecInt) = wt_MoRate(k) *wt_MoRateBranch              
            rateType1(k+nSR+nActRate+nRecInt) = 4
          enddo
                     
c         Read Mag recurrence weights (char and exp)
          read (10,*) nMagRecur
          call CheckDim ( nMagRecur, MAX_N1, 'MAX_N1    ' )
          read (10,*) (x(i),i=1,nMagRecur)
          read (10,*) (magRecurWt(i),i=1,nMagRecur)

c         Read past corresponding magnitude parameters. 
          do iRecur=1,nMagRecur
            read (10,*) x(iRecur), x(iRecur), x(iRecur)
          enddo

c         Read seismogenic thickness
          if ( isourceType .ne. 5) then
            read (10,*) nThick1(iFlt)
            call CheckDim ( nThick1(iFlt), MAX_WIDTH, 'MAX_WIDTH ' )
            read (10,*) (x(i),i=1,nThick1(iFlt))
            read (10,*) (faultThickWt(i),i=1,nThick1(iFlt))
          else
            nThick1(iFlt) = 1
            faultThickWt(1) = 1.
          endif
         
c         Read depth pdf
          read (10,*) iDepthModel       

c         Read Mag method (scaling relations or set values)
          read (10,*) iOverRideMag

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
c  ** not needed, we just use the iThick index for this weight
c            do iDip=2,n_Dip(iflt)
c              nRefMag(iThick) = nRefMag(iThick-1)
c              do i=1,nRefMag(iThick)
c                refMagWt(iThick,i) = refMagWt(iThick-1,i)
c              enddo
c              iThick = iThick + 1
c            enddo
          enddo

c         Read Past remaining input for this fault
          read (10,*) minMag, magStep, hxStep, hyStep, nRupArea, nRupWidth, minDepth
          read (10,*) (x2(k,iFlt),k=1,2), sigArea
          read (10,*) (x2(k,iFlt),k=1,2), sigWidth

c         Read ftype Models
          read (10,*) nFtypeModels(iFlt)
          do iFM=1,nFtypeModels(iFlt)
            read (10,*) ftmodelwt(iFM)
            read (10,*) nFtype1(iFM)
            read (10,*) (ftype1(iFM,k),k=1,nFtype1(iFM))
            read (10,*) (ftype_wt1(iFM,k), k=1,nFtype1(iFM))
          enddo

c        Load up Ftype Models and weights.
         nFM = 1
         do iFM=1,nFtypeModels(iFlt)
            do k=1,nFtype1(iFM)
               ftype(iFlt,nFM) = ftype1(iFM,k)
               ftype_wt(iFlt,nFM) = ftype_wt1(iFM,k)*ftmodelwt(iFM)
               nFm = nFm + 1
            enddo
         enddo
         nFm = nFm - 1
         nFtype(iFlt) = nFm

c  **** this is not correct, it treats all as epistemic
C        Set up Cum Wts for Ftype
         do iFM=1,nFtypeModels(iFlt)
c          write (*, '( i5,f10.4)') (iflt,  Ftype_Wt(iFlt,iFtype) )

           if (iFtype .eq. 1) then
             cumWt_Ftype(iFlt,iFtype) = Ftype_Wt(iFlt,iFtype)
           else
             cumWt_Ftype(iFlt,iFtype) = cumWt_Ftype(iFlt,iFtype-1) + Ftype_Wt(iFlt,iFtype)
           endif

         enddo

c     Load up parameter variations into large single dimension arrays        
         i1 = 0
	 do iDip=1,n_Dip(iFlt)
           do iThick=1,nThick1(iFlt)
             i1 = i1 + 1
             width_wt(i1) = faultThickWt(iThick) * dipWt(iDip)

             i = 0

             do iRecur=1,nMagRecur
               do iRate=1,nRate
                 if ( rateType1(iRate) .eq. 2 ) then
                   nb1 = 1
                 else
                   nb1 = n_bvalue
                 endif
                 do i_bValue=1,nb1
                   do iRefMag=1,nRefMag(iThick)
                     i = i + 1
                     call CheckDim ( i, MAXPARAM, 'MAXPARAM' )

c                    Set the weight for this set of parameters
                     if ( rateType1(iRate) .eq. 2 ) then
                       wt = RateWt1(iRate) * magRecurWt(iRecur) * refMagWt(iThick,iRefMag) 
                     else
                       wt = RateWt1(iRate) * bValueWt(i_bValue)* magRecurWt(iRecur) * refMagWt(iThick,iRefMag)                            
                     endif

c             Set cumulative weights for parameter variations     
                     if ( i .eq. 1 ) then
                        cumWt_param(iFlt,i1,i) = wt
                     else
                        cumWt_param(iFlt,i1,i) = wt +  cumWt_param(iFlt,i1,i-1)
                     endif
    
c                   End of Loop over iRefMag  
                  enddo
c                 End of Loop over ib_value
                enddo
c               End of Loop over iRate
	      enddo

c             End of Loop over iRecur
	    enddo
            nParamVar(iFlt,i1) = i
c            write (*,'( 2x,''iflt, dip-thick, nparam var ='',3i5)') iFlt, i1,  nParamVar(iFLt,i1)
c            write (*,'( 10f10.4)') cumwt_param(iFlt,i1,nParamVar(iFlt,i1))
            if ( cumwt_param(iFlt,i1,nParamVar(iFlt,i1)) .gt. 0.999 ) then
	      cumwt_param(iFlt,i1,nParamVar(iFlt,i1)) = 1.0
            endif

c         End of Loop over iThick1
          enddo

          probAct(iFlt) = probAct0

c         End of Loop over iDip
         enddo
         nWidth(iflt) = i1

         do iWidth=1,nWidth(iFlt)
            if ( iWidth .eq. 1) then
              cumWt_Width(iFlt,iWidth) = faultThickWt(iWidth)
            else
              cumWt_Width(iFlt,iWidth) = cumWt_Width(iFlt,iWidth-1) + faultThickWt(iWidth)
            endif
          enddo


c       End of Loop over iFlt2 - number of segments    
        enddo

c     End of Loop over iFlt
      enddo
      nFltTotal = iFlt

      return
      end

