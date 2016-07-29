c     Set array dimensions
      

      integer MAX_SITE, MAX_FLT, MAX_SEG,
     1        MAX_INTEN, MAX_PROB, MAX_DIP,  
     2        MAXPARAM, MAX_MAG, MAX_DIST, 
     3        MAX_EPS, MAX_N1, MAX_Files, 
     4        MAX_N2, MAX_Xcost, MAX_WIDTH, 
     5        MAX_SAMPLE, MAX_RISK
      integer MAXRUP, MAX_FIXED_MAG, MAX_MAGDIM, 
     1        MAX_FTYPE, MAX_ATTEN, MAX_ATTENTYPE


      PARAMETER ( MAX_SITE=1, MAX_FLT=80, MAX_SEG=10,
     1            MAX_INTEN=18, MAX_PROB=20, MAX_DIP=5,  
     2            MAXPARAM=40, MAX_MAG=3, MAX_DIST=15, 
     3            MAX_EPS=15, MAX_N1=150, MAX_Files=3, 
     4            MAX_N2=6, MAX_Xcost=10, MAX_WIDTH=15, 
     5            MAX_SAMPLE=60000, MAX_RISK=4000 )
      PARAMETER ( MAXRUP=4, MAX_FIXED_MAG=3, MAX_MAGDIM=3, 
     1            MAX_FTYPE=3, MAX_ATTEN=27, MAX_ATTENTYPE=2 )

