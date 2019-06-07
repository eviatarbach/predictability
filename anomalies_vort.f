*234567**************** get_anomalies_fourieranalysis.f *************************
C     
C The CLIMATOLOGY of a time series Y of NTM time steps, with NTSPYR time steps per year
C (i.e., with NTYR(=NTM/NTSPYR) years) can be obtained from the summation of the
C mean value of the time series YM plus the first NTSPYR/2 (or (NTSPYR-1)/2 for odd NTSPYR)
C harmonics. This implies that when creating anomalies wrt the climatology, variability
C at frequencies/periods explained by these first harmonics are eliminated from
C the original time series. So it is better to use harmonic analysis to cretae anomalies
C if one is interested in frequencies/periods that may be deleted by using long-term
C mean climatologies.
C
C If memory issues arise, compile as follows     
C gfortran -O3 -mcmodel=medium -g -o a.out get_anomalies_fourieranalysis.f
C where a.out is the executable file
C     
      PARAMETER(N3264=8)
      PARAMETER(UNDEF=9999.)
      PARAMETER(NGP=88838)
      PARAMETER(NI=1, NF=14235, NTM=NF-NI+1)  !01/01/1979-12/31/2017 Daily-365
C      PARAMETER(NI=2483, NF=4818, NTM=NF-NI+1)  !01/01/1982-12/27/2013 Pentad
C      PARAMETER(NI=1, NF=732, NTM=NF-NI+1)    
      PARAMETER(NTSPYR=365) !73 for pentad
      PARAMETER(KPER=NTSPYR)
      PARAMETER(NHARM=2)!Number of harmonics to be extracted, NHARM has to be <= NTSPYR/2,
C                       !or (NTSPYR-1)/2 if NTSPYr is odd; if equal then that would be
      DIMENSION Y(NTM)  !equivalent to extract the long-term climatological mean
      DIMENSION X(NGP), XX(NGP,NTM)
      DIMENSION YY(NGP,NTM)
      DIMENSION HI(NTM,NHARM)
      DIMENSION HIS(NGP,NTM,NHARM), YMS(NGP), TVARS(NGP)
      DIMENSION H(NGP,NTM,NHARM), H2(NGP,NHARM,2)
      DIMENSION VARIH(NHARM),R2(NHARM), VARIHS(NGP,NHARM),R2S(NGP,NHARM)
C
      CHARACTER*24 PATHI
      CHARACTER*24 PATHO
      CHARACTER*60 FILEI, FILEO, FILEOV, FILEOG, FILEOVH, FILEOH
C     
      PATHI='/lustre/ebach/causality/'
      PATHO='/lustre/ebach/causality/'
C     
      FILEI  ='vort850_daily365_1979-2017.dat'
      FILEO  ='vortah850_daily365_1979-2017.dat'
      FILEOV ='vort850_daily365_1979-2017_meanvar.dat'
      FILEOH ='vort850_daily365_1979-2017_harms.dat'
      FILEOVH='vort850_daily365_1979-2017_varharms.dat'
      FILEOG ='grid_vort850_daily365.dat'
C      FILEI  ='vort850ocn_ncep_2.5x2.5_pentad_1948-2014.dat'
C      FILEO  ='vortah850ocn_ncep_2.5x2.5_pen1982-2013.dat'
C      FILEOV ='vort850ocn_ncep_2.5x2.5_pen1982-2013_meanvar.dat'
C      FILEOH ='vort850ocn_ncep_2.5x2.5_pen1982-2013_harms.dat'
C      FILEOVH='vort850ocn_ncep_2.5x2.5_pen1982-2013_varharms.dat'
C      FILEOG ='grid_vort850ocn_ncep_2.5x2.5.dat'
C     
      OPEN(11,FILE=PATHO//FILEI,FORM='UNFORMATTED',
     -       STATUS='OLD',ACCESS='DIRECT',RECL=NGP*N3264)
      OPEN(13,FILE=PATHO//FILEO,FORM='UNFORMATTED',
     -       STATUS='UNKNOWN',ACCESS='DIRECT',RECL=NGP*N3264)
      OPEN(14,FILE=PATHO//FILEOV,FORM='UNFORMATTED',
     -       STATUS='UNKNOWN',ACCESS='DIRECT',RECL=NGP*N3264)
      OPEN(15,FILE=PATHO//FILEOH,FORM='UNFORMATTED',
     -       STATUS='UNKNOWN',ACCESS='DIRECT',RECL=NGP*N3264)
      OPEN(16,FILE=PATHO//FILEOVH,FORM='UNFORMATTED',
     -       STATUS='UNKNOWN',ACCESS='DIRECT',RECL=NGP*N3264)
C    
      OPEN(22,FILE=FILEOG,STATUS='unknown')
C     
      WRITE(*,*)'GETTING TIME SERIES FROM MAPS'
      DO I=NI, NF
         J=I-NI+1
         READ(11,REC=I) X
         DO LO = 1, NGP
           XX(LO,J)=X(LO)
         ENDDO
      ENDDO
      DO J=1, NTM
         N=0
         DO LO=1, NGP
               X(LO)=XX(LO,J)
	       IF(X(LO).NE.UNDEF) THEN 
	          N = N + 1
                  YY(N,J)=X(LO)
		  IF(J.EQ.1)
     -               WRITE(22,*) J, LO, N
	          ENDIF
	 ENDDO
         WRITE(*,*)N,' DEFINED GRID POINTS'
      ENDDO
      CLOSE(22)
C
      NUMDEF=N
C     
      WRITE(*,*)'GETTING HARMONICS'
C     
      NMODD=MOD(NTM,2)
      NH=NTM/2
      IF(NMODD.NE.0) NH = (NTM-1)/2
      WRITE(*,*)'Getting',NHARM,'Harmonics of a maximum of',NH,
     -          'for ',NUMDEF,'grid points'
C  
      DO I=1, NUMDEF
	 DO J =1, NTM
             Y(J) = YY(I,J)
         ENDDO
         CALL FOURIERANALYSIS(Y,NTM,NHARM,HI,KPER,YM,TVAR,VARIH,R2)
         IF(I.EQ.4721) THEN
           DO J =1, NTM
              WRITE(*,*)'J=',J,'Y=',Y(J),'Var 1/2',VARIH(1),VARIH(2),
     -                  'R^2 1/2',R2(1),R2(2)
           ENDDO
         ENDIF
         YMS(I)=YM
         TVARS(I)=TVAR
         DO K=1, NHARM
            VARIHS(I,K)=VARIH(K)
            R2S(I,K)=R2(K)
            DO J=1, NTM
               HIS(I,J,K)=HI(J,K)
            ENDDO
         ENDDO
         DO J=1, NTM
            SH=0.0
            DO K=1, NHARM
               SH=SH+HI(J,K)
            ENDDO
            YY(I,J)=YY(I,J)-(YM+SH)
         ENDDO
      ENDDO
C
C Saving Anomalies
      WRITE(*,*)'Saving Anomalies'
      DO J=1, NTM
         DO LO=1, NGP
           H(LO,J,1)=UNDEF
         ENDDO
      ENDDO
      DO J=1, NTM
         OPEN(22,FILE=FILEOG,STATUS='UNKNOWN')
         DO I=1, NUMDEF
            READ(22,*)MES,LON,NP
            H(LON,J,1)=YY(NP,J)     
         ENDDO
         CLOSE(22)
      ENDDO
      NRECW=0
      DO J=1, NTM          !NTM time steps
         NRECW=NRECW+1
         DO LO=1, NGP ! NLON*NLAT grid points
           X(LO)=H(LO,J,1)
         ENDDO
         WRITE(13,REC=NRECW) X
      ENDDO
C
C Saving Mean and Variance of the original Raw Sample Time Series
      WRITE(*,*)'Saving Mean and Variance of the original Raw Maps'
      DO LO=1, NGP
        H2(LO,1,1)=UNDEF
        H2(LO,1,2)=UNDEF
      ENDDO
      OPEN(22,FILE=FILEOG,STATUS='UNKNOWN')
      DO I=1, NUMDEF
         READ(22,*)MES,LON,NP
         H2(LON,1,1)=YMS(NP)     
         H2(LON,1,2)=TVARS(NP)     
      ENDDO
      CLOSE(22)
      NRECW=0
      DO NR=1, 2
         NRECW=NRECW+1
         DO LO=1, NGP ! NLON*NLAT grid points
           X(LO)=H2(LO,1,NR)
         ENDDO
         WRITE(14,REC=NRECW) X
      ENDDO
C     
C Saving Individual Harmonics
      WRITE(*,*)'Saving Individual Harmonics'
      DO J=1, NTM
         DO K=1, NHARM
            DO LO=1, NGP
              H(LO,J,K)=UNDEF
            ENDDO
         ENDDO
      ENDDO
      DO J=1, NTM
         DO K=1, NHARM
            OPEN(22,FILE=FILEOG,STATUS='UNKNOWN')
            DO I=1, NUMDEF
               READ(22,*)MES,LON,NP
               H(LON,J,K)=HIS(NP,J,K)     
            ENDDO
            CLOSE(22)
        ENDDO
      ENDDO
      NRECW=0
      DO J=1, NTM          !NTM time steps
         DO K=1, NHARM  !NHARM levels
            NRECW=NRECW+1
            DO LO=1, NGP ! NLON*NLAT grid points
              X(LO)=H(LO,J,K)
            ENDDO
            WRITE(15,REC=NRECW) X
         ENDDO
      ENDDO
C
C Saving Variances and Coeff of Determination of the Harmonics
      WRITE(*,*)'Saving Variances and Coeff of Det of the Harmonics'
      DO K=1, NHARM
         DO LO=1, NGP
           H2(LO,K,1)=UNDEF
           H2(LO,K,2)=UNDEF
         ENDDO
      ENDDO
      DO K=1, NHARM
         OPEN(22,FILE=FILEOG,STATUS='UNKNOWN')
         DO I=1, NUMDEF
            READ(22,*)MES,LON,NP
            H2(LON,K,1)=VARIHS(NP,K)     
            H2(LON,K,2)=R2S(NP,K)     
         ENDDO
         CLOSE(22)
      ENDDO
      NRECW=0
      DO K=1, NHARM  !NHARM levels
         DO NR=1, 2
            NRECW=NRECW+1
            DO LO=1, NGP ! NLON*NLAT grid points
              X(LO)=H2(LO,K,NR)
            ENDDO
            WRITE(16,REC=NRECW) X
         ENDDO
      ENDDO
C     
C     
C gfortran get_anomalies_fourieranalysis.f
C ./a.out 
C     
      END     
C
      SUBROUTINE FOURIERANALYSIS(Y,NTM,NHARM,HI,KPER,YM,TVAR,VARIH,R2)
C     
      DIMENSION Y(NTM), HI(NTM,NHARM)
      DIMENSION A(NHARM), B(NHARM), C(NHARM), PHIR(NHARM)
      DIMENSION VARIH(NHARM), R2(NHARM)
C     DIMENSION YHI(NTM,NHARM)
C      DIMENSION PHI(NHARM)
C
      PI=3.1415927
      S = 0.0
      DO J =1, NTM
	 S = S + Y(J)
      ENDDO
      YM=S/FLOAT(NTM)
C
C Getting Harmonics
      DO K = 1, NHARM
	 SC=0.0
	 SS=0.0
	 DO J=1, NTM
	    ARG=2.0*PI*FLOAT(K)*FLOAT(J)/FLOAT(KPER)
	    SC = SC+Y(J)*COS(ARG)
	    SS = SS+Y(J)*SIN(ARG)
	 ENDDO
	 A(K)=2.0*SC/FLOAT(NTM)
	 B(K)=2.0*SS/FLOAT(NTM)
         C(K)=SQRT(A(K)**2+B(K)**2)
	 IF(A(K).GT.0.) THEN
	    IF (B(K).LT.0.0) PHIR(K)=ATAN(B(K)/A(K)) +2.0* PI
            IF (B(K).EQ.0.0) PHIR(K)=0.0
            IF (B(K).GT.0.0) PHIR(K)=ATAN(B(K)/A(K))
         ELSEIF(A(K).LT.0.0) THEN
	    IF (B(K).EQ.0.0) PHIR(K)=PI
	       PHIR(K) = ATAN(B(K)/A(K))+PI
	    ELSE
               PHIR(K)=PI/2.0
               IF(B(K).LT.0.0) PHIR(K)=3.0*PI/2.0
	 ENDIF
C	 PHI(K)=PHIR(K)*360.0/2.0/PI !Changing phase from radians to degrees
      ENDDO
C
C Reconstructing Original data by adding harmonics
      DO J=1, NTM
         DO K=1, NHARM
            ARGG=2.0*PI*FLOAT(K)*FLOAT(J)/FLOAT(KPER)-PHIR(K)
            HI(J,K)=C(K)*COS(ARGG)!Individual harmonics
C            YHI(J,K)=YM+HI(J,K) !Reconstructed data from individual harmonics
         ENDDO
      ENDDO
C
C Getting variance from time series
      SV = 0.0
      DO J =1, NTM
	 SV = SV + (Y(J)-YM)**2
      ENDDO
      TVAR=SV/FLOAT(NTM-1)
C
C Getting variances from individeal harmonics
      DO K=1, NHARM
         VARIH(K)=0.5*C(K)**2
         R2(K)=VARIH(K)*FLOAT(NTM)/(FLOAT(NTM-1)*TVAR)
      ENDDO
C
      RETURN     
      END     
      
