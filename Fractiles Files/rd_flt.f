      
      subroutine Rd_Fault_Data (version, nFltTotal,nFlt0, probAct, AttenType, 
     2           nSegModel, cumWt_segModel, nFtype, wt_ftype, f_start, f_num,
     1           faultFlag, al_Segwt, nBR_all, wt_cum_all, fname )

      implicit none
      include 'fract.h'
      
      integer Attentype(MAX_FLT),
     1        nFlt0, nFtype(MAX_FLT), faultflag(MAX_FLT,MAX_SEG,MAX_FLT),
     2        f_start(1), f_num(1), nSegModel(1), nFltTotal
      real al_segWt(MAX_FLT), cumWt_segModel(MAX_FLT,MAX_SEG), version, 
     1     ProbAct(MAX_FLT)
      character*80 fname(MAX_FLT) 
      real wt_cum_all(MAX_FLT,12,10,10)
      integer nBR_all(MAX_FLT,MAXPARAM,10)
      real wt_ftype(MAX_FLT,5,5)

      if (version .eq. 45.1) then
          call Rd_Fault_Data_45_1 (nFltTotal, nFlt0, Wt_cum_all, nBR_all, 
     1           probAct, cumWt_segModel, nSegModel, AttenType, wt_Ftype, 
     2           nFtype, f_start, f_num, faultFlag, al_Segwt, fname )
    
      elseif (version .eq. 45.2) then
          call Rd_Fault_Data_45_2 (nFltTotal, nFlt0, Wt_cum_all, nBR_all, 
     1           probAct, cumWt_segModel, nSegModel, AttenType, wt_Ftype, 
     2           nFtype, f_start, f_num, faultFlag, al_Segwt, fname )
      else
          write (*,*) 'Incompatible fault file, use Haz45.2 or Haz45.1'
          stop 99
      endif
        
      return 
      end
      
c ----------------------------------------------------------------------

      subroutine Rd_Fault_Data_45_1 (nFltTotal, nFlt0, Wt_cum_all, nBR_all, 
     1           probAct, 
     1           cumWt_segModel, nSegModel, AttenType, wt_Ftype, 
     2           nFtype, f_start, f_num, faultFlag, al_Segwt, fname )

      implicit none
      include 'fract.h'
      
      integer Attentype(MAX_FLT),
     1        nFlt0, nFlt2, nFm, nFtype(MAX_FLT,MAXPARAM),
     2        nThick1, directflag, synflag, 
     3        nSR, nActRAte, nRecInt, n_Dip, nRefMag(MAX_WIDTH), 
     4        faultflag(MAX_FLT,MAX_SEG,MAX_FLT), nsyn, nFltTotal
      integer f_start(MAX_FLT), f_num(MAX_FLT), nSegModel(MAX_FLT),
     1        iflt, iflt0, iflt2, k, i, iCoor, isourceType,
     2        iFM, nRupArea, nRupWidth, iThick, iThick1, nMoRate,
     3        nMagRecur, n_bValue, ii, ipt, nDownDip, insyn, iRecur,
     4        iDepthModel, nfp 
      real segwt(MAX_FLT,MAX_FLT), al_segWt(MAX_FLT), minmag, magstep, 
     1     magRecurWt(MAXPARAM), dipWt(MAXPARAM), bValueWt(MAXPARAM), minDepth, 
     4     x(MAXPARAM), ProbAct(MAX_FLT) 
      real ftmodelwt(MAXPARAM), ftype1(MAXPARAM,MAXPARAM), 
     1      magsyn, rupsyn, jbsyn, 
     2     seismosyn, hyposyn, wtsyn, wt_srBranch, wt_recIntBranch, wt_SR(MAXPARAM), 
     3     wt_RecInt(MAXPARAM), bValue2(MAXPARAM), actRate(MAXPARAM), 
     4     actRateWt(MAXPARAM), wt_MoRate(MAXPARAM)
      real faultThickWt(MAX_WIDTH), refMagWt(MAX_WIDTH,MAXPARAM), fZ,
     1     hxStep, hyStep,
     2     probAct0, sampleStep, dip1, top, fLong, fLat, wt_ActRateBranch,
     3     wt_MoRateBranch, sigArea, sigWidth, sum, attensyn, hwsyn, ftypesyn
      character*80 fName1, fName(MAX_FLT)
      
      real wt_cum_all(MAX_FLT,12,10,10)
      integer nBR_all(MAX_FLT,MAXPARAM,10), iNode, iBR, iOverRideMag
      real cumWt_segModel(MAX_FLT,MAX_FLT)
      real wt_ftype(MAX_FLT,5,5)

C     Input Fault Parameters
      read (10,*) iCoor
      read (10,*) NFLT0

      call CheckDim ( NFLT0, MAX_FLT, 'MAX_FLT' )
  
      iflt = 0

C.....Loop over each fault in the source file....
      DO 900 iFlt0=1,NFLT0
        read (10,'( a80)') fName1
        read (10,*) probAct0

c       Read number of segmentation models for this fault system       
        read (10,*) nSegModel(iFlt0)
        read (10,*) (segWt(iFlt0,k),k=1,nSegModel(iFlt0))

C       Set up cum wt for segment models for this fault system
        do k=1,nSegModel(iFlt0)
          if ( k .eq. 1) then
            cumWt_segModel(iFlt0,1) = segWt(iFlt0,1)
          else
            cumWt_segModel(iFlt0,k) = segWt(iFlt0,k) + cumWt_segModel(iFlt0,k-1)
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
          read(10,'( a80)') fname(iFlt)
          read (10,*) isourceType, attenType(iFlt), sampleStep, directflag, synflag

          write (*,'( 2x,''iFlt'',i5, 2x,a80)') iFlt, fname(iFLt)

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
            read (10,*) n_Dip
            read (10,*) (x(i),i=1,n_Dip)
            read (10,*) (dipWt(i),i=1,n_Dip)
          else
            n_Dip = 1
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
              write (*,'( 2x,''rate method weights do not sum to unity for fault, '',a30)') fName(iFLt)
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
            read (10,*) nThick1
            call CheckDim ( nThick1, MAX_WIDTH, 'MAX_WIDTH ' )
            read (10,*) (x(i),i=1,nThick1)
            read (10,*) (faultThickWt(i),i=1,nThick1)
          else
            nThick1 = 1
            faultThickWt(1) = 1.
          endif

c         Read depth pdf
          read (10,*) iDepthModel       

c        Read Mag method (scaling relations or set values)
         read (10,*) iOverRideMag
         if ( iOverRideMag .ne. 1 ) then
           write (*,'( 2x,''iOverRideMag flag option not working'')') 
           stop 99
         endif

c         Read reference mags for each fault thickness
          do iThick=1,nThick1
            read (10,*) nRefMag(iThick)
            read (10,*) (x(i),i=1,nRefMag(iThick))
            read (10,*) (refMagWt(iThick,i),i=1,nRefMag(iThick))
          enddo

c         Read Past remaining input for this fault
          read (10,*) minMag, magStep, hxStep, hyStep, nRupArea, nRupWidth, minDepth
          read (10,*) (x(k),k=1,2), sigArea
          read (10,*) (x(k),k=1,2), sigWidth

c         Read ftype Models
          read (10,*) nFM
          do iFM=1,nFM
            read (10,*) ftmodelwt(iFM)
            read (10,*) nFtype(iFlt,iFM)
            read (10,*) (ftype1(iFM,k),k=1,nFtype(iFlt,iFM))
            read (10,*) (wt_ftype(iFlt,iFM,k), k=1,nFtype(iFlt,iFM))
          enddo

c         Add index on thickness because the refMag is correlated with thickness
c         For other param, just enter them for each of the thickness to make this easier to track
          do iThick=1,nThick1

c          Load seismo thickness wts into global cumulative wt array
           iNode = 1
           nBR_all(iflt,iNode,iThick) = nThick1
           wt_cum_all(iflt,iNode,iThick,1) = faultThickWt(1)
           do iBR=2,nBR_all(iflt,iNode,iThick)
            wt_cum_all(iflt,iNode,iThick,iBR) = wt_cum_all(iflt,iNode,iThick,iBR-1) + faultThickWt(iBR)
           enddo

c          Load dip wts into global cumulative wt array (node 1)
           iNode = 2
           nBR_all(iflt,iNode,iThick) = n_Dip
           wt_cum_all(iflt,iNode,iThick,1) = dipWt(1)
           do iBR=2,nBR_all(iflt,iNode,iThick)
            wt_cum_all(iflt,iNode,iThick,iBR) = wt_cum_all(iflt,iNode,iThick,iBR-1) + dipWt(iBR)
           enddo

c          Load Ftype Model wts into global cumulative wt array
           iNode = 3
           nBR_all(iflt,iNode,iThick) = nFM
           wt_cum_all(iflt,iNode,iThick,1) = ftmodelwt(1)
           do iBR=2,nBR_all(iflt,iNode,iThick)
             wt_cum_all(iflt,iNode,iThick,iBR) = wt_cum_all(iflt,iNode,iThick,iBR-1) + ftmodelwt(iBR)
           enddo

c          Load Mag Recur  wts into global cumulative wt array
           iNode = 4
           nBR_all(iflt,iNode,iThick) = nMagRecur
           wt_cum_all(iflt,iNode,iThick,1) = magRecurWt(1)
           do iBR=2,nBR_all(iflt,iNode,iThick)
             wt_cum_all(iflt,iNode,iThick,iBR) = wt_cum_all(iflt,iNode,iThick,iBR-1) + magRecurWt(iBR)
           enddo

c          Load ref mag  wts into global cumulative wt array
           iNode = 5
           nBR_all(iflt,iNode,iThick) = nRefMag(iThick)
           wt_cum_all(iflt,iNode,iThick,1) = refMagWt(ithick,1)
           do iBR=2,nBR_all(iflt,iNode,iThick)
             wt_cum_all(iflt,iNode,iThick,iBR) = wt_cum_all(iflt,iNode,iThick,iBR-1) 
     1               + refMagWt(ithick,iBR)
           enddo

c          Load rate method wts into global cumulative wt array
           iNode = 6
           nBR_all(iflt,iNode,iThick) = 4
           wt_cum_all(iflt,iNode,iThick,1) = wt_srBranch
           wt_cum_all(iflt,iNode,iThick,2) = wt_ActRateBranch +  wt_cum_all(iflt,iNode,iThick,1)
           wt_cum_all(iflt,iNode,iThick,3) = wt_recIntBranch +  wt_cum_all(iflt,iNode,iThick,2)
           wt_cum_all(iflt,iNode,iThick,4) = wt_MoRateBranch +  wt_cum_all(iflt,iNode,iThick,3)

c          Load slip-rate wts into global cumulative wt array
           iNode = 7
           nBR_all(iflt,iNode,iThick) = nSR
           if ( nSR .ne. 0 ) then
             wt_cum_all(iflt,iNode,iThick,1) = wt_sr(1)
             do iBR=2,nBR_all(iflt,iNode,iThick)
               wt_cum_all(iflt,iNode,iThick,iBR) = wt_cum_all(iflt,iNode,iThick,iBR-1) + wt_sr(iBR)
             enddo
           endif
                
c          Load activity rate - b-value pairs wts into global cumulative wt array
           iNode = 8
           nBR_all(iflt,iNode,iThick) = nActRate
           if ( nActRate .ne. 0 ) then
            wt_cum_all(iflt,iNode,iThick,1) = actRateWt(1)
            do iBR=2,nBR_all(iflt,iNode,iThick)
              wt_cum_all(iflt,iNode,iThick,iBR) = wt_cum_all(iflt,iNode,iThick,iBR-1) + actRateWt(iBR)
            enddo
           endif

c          Load recurrence intervals wts into global cumulative wt array
           iNode = 9
           nBR_all(iflt,iNode,iThick) = nRecInt
           if ( nRecInt .ne. 0 ) then
            wt_cum_all(iflt,iNode,iThick,1) = wt_recInt(1)
            do iBR=2,nBR_all(iflt,iNode,iThick)
              wt_cum_all(iflt,iNode,iThick,iBR) = wt_cum_all(iflt,iNode,iThick,iBR-1) + wt_recInt(iBR)
            enddo
           endif

c          Load recurrence intervals wts into global cumulative wt array
           iNode = 10
           nBR_all(iflt,iNode,iThick) = nMoRate
           if ( nMoRate .ne. 0 ) then
            wt_cum_all(iflt,iNode,iThick,1) = wt_MoRate(1)
            do iBR=2,nBR_all(iflt,iNode,iThick)
              wt_cum_all(iflt,iNode,iThick,iBR) = wt_cum_all(iflt,iNode,iThick,iBR-1) + wt_MoRate(iBR)
            enddo
           endif
                       
c          Load b-values wts into global cumulative wt array (node 11)
           iNode = 11
           nBR_all(iflt,iNode,iThick) = n_bValue
           wt_cum_all(iflt,iNode,iThick,1) = bValueWt(1)
           do iBR=2,nBR_all(iflt,iNode,iThick)
            wt_cum_all(iflt,iNode,iThick,iBR) = wt_cum_all(iflt,iNode,iThick,iBR-1) + bValueWt(iBR)
           enddo

c          End of Loop over iThick1
          enddo
                   
         probAct(iFlt) = probAct0
         
c       End of Loop over iFlt2 - number of segments    
        enddo

c     End of Loop over iFlt
  900 continue
      nFltTotal = iFlt
      
      return
      end

c -------------------

      subroutine Rd_Fault_Data_45_2 (nFltTotal, nFlt0, Wt_cum_all, nBR_all, 
     1           probAct, 
     1           cumWt_segModel, nSegModel, AttenType, wt_Ftype, 
     2           nFtype, f_start, f_num, faultFlag, al_Segwt, fname )

      implicit none
      include 'fract.h'
      
      integer Attentype(MAX_FLT),
     1        nFlt0, nFlt2, nFm, nFtype(MAX_FLT,MAXPARAM),
     2        nThick1, directflag, synflag, 
     3        nSR, nActRAte, nRecInt, n_Dip, nRefMag(MAX_WIDTH), 
     4        faultflag(MAX_FLT,MAX_SEG,MAX_FLT), nsyn, nFltTotal
      integer f_start(MAX_FLT), f_num(MAX_FLT), nSegModel(MAX_FLT),
     1        iflt, iflt0, iflt2, k, i, iCoor, isourceType,
     2        iFM, nRupArea, nRupWidth, iThick, iThick1, nMoRate,
     3        nMagRecur, n_bValue, ii, ipt, nDownDip, insyn, iRecur,
     4        iDepthModel, nfp 
      real segwt(MAX_FLT,MAX_FLT), al_segWt(MAX_FLT), minmag, magstep, 
     1     magRecurWt(MAXPARAM), dipWt(MAXPARAM), bValueWt(MAXPARAM), minDepth, 
     4     x(MAXPARAM), ProbAct(MAX_FLT) 
      real ftmodelwt(MAXPARAM), ftype1(MAXPARAM,MAXPARAM), 
     1      magsyn, rupsyn, jbsyn, 
     2     seismosyn, hyposyn, wtsyn, wt_srBranch, wt_recIntBranch, wt_SR(MAXPARAM), 
     3     wt_RecInt(MAXPARAM), bValue2(MAXPARAM), actRate(MAXPARAM), 
     4     actRateWt(MAXPARAM), wt_MoRate(MAXPARAM)
      real faultThickWt(MAX_WIDTH), refMagWt(MAX_WIDTH,MAXPARAM), fZ,
     1     hxStep, hyStep,
     2     probAct0, sampleStep, dip1, top, fLong, fLat, wt_ActRateBranch,
     3     wt_MoRateBranch, sigArea, sigWidth, sum, attensyn, hwsyn, ftypesyn
      character*80 fName1, fName(MAX_FLT)
      
      real wt_cum_all(MAX_FLT,12,10,10)
      integer nBR_all(MAX_FLT,MAXPARAM,10), iNode, iBR
      real cumWt_segModel(MAX_FLT,MAX_FLT)
      real wt_ftype(MAX_FLT,5,5)

C     Input Fault Parameters
      read (10,*) iCoor
      read (10,*) NFLT0

      call CheckDim ( NFLT0, MAX_FLT, 'MAX_FLT' )
  
      iflt = 0

C.....Loop over each fault in the source file....
      DO 900 iFlt0=1,NFLT0
        read (10,'( a80)') fName1
        read (10,*) probAct0

c       Read number of segmentation models for this fault system       
        read (10,*) nSegModel(iFlt0)
        read (10,*) (segWt(iFlt0,k),k=1,nSegModel(iFlt0))

C       Set up cum wt for segment models for this fault system
        do k=1,nSegModel(iFlt0)
          if ( k .eq. 1) then
            cumWt_segModel(iFlt0,1) = segWt(iFlt0,1)
          else
            cumWt_segModel(iFlt0,k) = segWt(iFlt0,k) + cumWt_segModel(iFlt0,k-1)
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
          read(10,'( a80)') fname(iFlt)
          read (10,*) isourceType, attenType(iFlt), sampleStep, directflag, synflag

          write (*,'( 2x,''iFlt'',i5, 2x,a80)') iFlt, fname(iFLt)

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
            read (10,*) n_Dip
            read (10,*) (x(i),i=1,n_Dip)
            read (10,*) (dipWt(i),i=1,n_Dip)
          else
            n_Dip = 1
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
              write (*,'( 2x,''rate method weights do not sum to unity for fault, '',a30)') fName(iFLt)
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
            read (10,*) nThick1
            call CheckDim ( nThick1, MAX_WIDTH, 'MAX_WIDTH ' )
            read (10,*) (x(i),i=1,nThick1)
            read (10,*) (faultThickWt(i),i=1,nThick1)
          else
            nThick1 = 1
            faultThickWt(1) = 1.
          endif

c         Read depth pdf
          read (10,*) iDepthModel       

c         Read reference mags for each fault thickness
          do iThick=1,nThick1
            read (10,*) nRefMag(iThick)
            read (10,*) (x(i),i=1,nRefMag(iThick))
            read (10,*) (refMagWt(iThick,i),i=1,nRefMag(iThick))
          enddo

c         Read Past remaining input for this fault
          read (10,*) minMag, magStep, hxStep, hyStep, nRupArea, nRupWidth, minDepth
          read (10,*) (x(k),k=1,2), sigArea
          read (10,*) (x(k),k=1,2), sigWidth

c         Read ftype Models
          read (10,*) nFM
          do iFM=1,nFM
            read (10,*) ftmodelwt(iFM)
            read (10,*) nFtype(iFlt,iFM)
            read (10,*) (ftype1(iFM,k),k=1,nFtype(iFlt,iFM))
            read (10,*) (wt_ftype(iFlt,iFM,k), k=1,nFtype(iFlt,iFM))
          enddo

c         Add index on thickness because the refMag is correlated with thickness
c         For other param, just enter them for each of the thickness to make this easier to track
          do iThick=1,nThick1

c          Load seismo thickness wts into global cumulative wt array
           iNode = 1
           nBR_all(iflt,iNode,iThick) = nThick1
           wt_cum_all(iflt,iNode,iThick,1) = faultThickWt(1)
           do iBR=2,nBR_all(iflt,iNode,iThick)
            wt_cum_all(iflt,iNode,iThick,iBR) = wt_cum_all(iflt,iNode,iThick,iBR-1) + faultThickWt(iBR)
           enddo

c          Load dip wts into global cumulative wt array (node 1)
           iNode = 2
           nBR_all(iflt,iNode,iThick) = n_Dip
           wt_cum_all(iflt,iNode,iThick,1) = dipWt(1)
           do iBR=2,nBR_all(iflt,iNode,iThick)
            wt_cum_all(iflt,iNode,iThick,iBR) = wt_cum_all(iflt,iNode,iThick,iBR-1) + dipWt(iBR)
           enddo

c          Load Ftype Model wts into global cumulative wt array
           iNode = 3
           nBR_all(iflt,iNode,iThick) = nFM
           wt_cum_all(iflt,iNode,iThick,1) = ftmodelwt(1)
           do iBR=2,nBR_all(iflt,iNode,iThick)
             wt_cum_all(iflt,iNode,iThick,iBR) = wt_cum_all(iflt,iNode,iThick,iBR-1) + ftmodelwt(iBR)
           enddo

c          Load Mag Recur  wts into global cumulative wt array
           iNode = 4
           nBR_all(iflt,iNode,iThick) = nMagRecur
           wt_cum_all(iflt,iNode,iThick,1) = magRecurWt(1)
           do iBR=2,nBR_all(iflt,iNode,iThick)
             wt_cum_all(iflt,iNode,iThick,iBR) = wt_cum_all(iflt,iNode,iThick,iBR-1) + magRecurWt(iBR)
           enddo

c          Load ref mag  wts into global cumulative wt array
           iNode = 5
           nBR_all(iflt,iNode,iThick) = nRefMag(iThick)
           wt_cum_all(iflt,iNode,iThick,1) = refMagWt(ithick,1)
           do iBR=2,nBR_all(iflt,iNode,iThick)
             wt_cum_all(iflt,iNode,iThick,iBR) = wt_cum_all(iflt,iNode,iThick,iBR-1) 
     1               + refMagWt(ithick,iBR)
           enddo

c          Load rate method wts into global cumulative wt array
           iNode = 6
           nBR_all(iflt,iNode,iThick) = 4
           wt_cum_all(iflt,iNode,iThick,1) = wt_srBranch
           wt_cum_all(iflt,iNode,iThick,2) = wt_ActRateBranch +  wt_cum_all(iflt,iNode,iThick,1)
           wt_cum_all(iflt,iNode,iThick,3) = wt_recIntBranch +  wt_cum_all(iflt,iNode,iThick,2)
           wt_cum_all(iflt,iNode,iThick,4) = wt_MoRateBranch +  wt_cum_all(iflt,iNode,iThick,3)

c          Load slip-rate wts into global cumulative wt array
           iNode = 7
           nBR_all(iflt,iNode,iThick) = nSR
           if ( nSR .ne. 0 ) then
             wt_cum_all(iflt,iNode,iThick,1) = wt_sr(1)
             do iBR=2,nBR_all(iflt,iNode,iThick)
               wt_cum_all(iflt,iNode,iThick,iBR) = wt_cum_all(iflt,iNode,iThick,iBR-1) + wt_sr(iBR)
             enddo
           endif
                
c          Load activity rate - b-value pairs wts into global cumulative wt array
           iNode = 8
           nBR_all(iflt,iNode,iThick) = nActRate
           if ( nActRate .ne. 0 ) then
            wt_cum_all(iflt,iNode,iThick,1) = actRateWt(1)
            do iBR=2,nBR_all(iflt,iNode,iThick)
              wt_cum_all(iflt,iNode,iThick,iBR) = wt_cum_all(iflt,iNode,iThick,iBR-1) + actRateWt(iBR)
            enddo
           endif

c          Load recurrence intervals wts into global cumulative wt array
           iNode = 9
           nBR_all(iflt,iNode,iThick) = nRecInt
           if ( nRecInt .ne. 0 ) then
            wt_cum_all(iflt,iNode,iThick,1) = wt_recInt(1)
            do iBR=2,nBR_all(iflt,iNode,iThick)
              wt_cum_all(iflt,iNode,iThick,iBR) = wt_cum_all(iflt,iNode,iThick,iBR-1) + wt_recInt(iBR)
            enddo
           endif

c          Load recurrence intervals wts into global cumulative wt array
           iNode = 10
           nBR_all(iflt,iNode,iThick) = nMoRate
           if ( nMoRate .ne. 0 ) then
            wt_cum_all(iflt,iNode,iThick,1) = wt_MoRate(1)
            do iBR=2,nBR_all(iflt,iNode,iThick)
              wt_cum_all(iflt,iNode,iThick,iBR) = wt_cum_all(iflt,iNode,iThick,iBR-1) + wt_MoRate(iBR)
            enddo
           endif
                       
c          Load b-values wts into global cumulative wt array (node 11)
           iNode = 11
           nBR_all(iflt,iNode,iThick) = n_bValue
           wt_cum_all(iflt,iNode,iThick,1) = bValueWt(1)
           do iBR=2,nBR_all(iflt,iNode,iThick)
            wt_cum_all(iflt,iNode,iThick,iBR) = wt_cum_all(iflt,iNode,iThick,iBR-1) + bValueWt(iBR)
           enddo

c          End of Loop over iThick1
          enddo
                   
         probAct(iFlt) = probAct0
         
c       End of Loop over iFlt2 - number of segments    
        enddo

c     End of Loop over iFlt
  900 continue
      nFltTotal = iFlt
      
      return
      end
