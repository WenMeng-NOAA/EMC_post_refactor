  MODULE CALRH_MODULE
!-------------------------------------------------------------------------------------
! Computes RH using various algorithms.
! The NAM v4.1.18 ALGORITHM is selected as default for the UPP 2020 unification. 
!
! program log:
!   MAY, 2020    Jesse Meng   Initial code
!-------------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------------
!
  contains

      SUBROUTINE CALRH(P1,T1,Q1,RH)

      use ctlblk_mod, only: jsta, jend, im

      REAL,dimension(IM,jsta:jend),intent(in)    :: P1,T1
      REAL,dimension(IM,jsta:jend),intent(inout) :: Q1
      REAL,dimension(IM,jsta:jend),intent(out)   :: RH

      CALL CALRH_NAM(P1,T1,Q1,RH)

      END SUBROUTINE CALRH
!
!-------------------------------------------------------------------------------------
!
      SUBROUTINE CALRH_NAM(P1,T1,Q1,RH)
!      SUBROUTINE CALRH(P1,T1,Q1,RH)
!$$$  SUBPROGRAM DOCUMENTATION BLOCK
!                .      .    .     
! SUBPROGRAM:    CALRH       COMPUTES RELATIVE HUMIDITY
!   PRGRMMR: TREADON         ORG: W/NP2      DATE: 92-12-22       
!     
! ABSTRACT:  
!     THIS ROUTINE COMPUTES RELATIVE HUMIDITY GIVEN PRESSURE, 
!     TEMPERATURE, SPECIFIC HUMIDITY. AN UPPER AND LOWER BOUND
!     OF 100 AND 1 PERCENT RELATIVE HUMIDITY IS ENFORCED.  WHEN
!     THESE BOUNDS ARE APPLIED THE PASSED SPECIFIC HUMIDITY 
!     ARRAY IS ADJUSTED AS NECESSARY TO PRODUCE THE SET RELATIVE
!     HUMIDITY.
!   .     
!     
! PROGRAM HISTORY LOG:
!   ??-??-??  DENNIS DEAVEN
!   92-12-22  RUSS TREADON - MODIFIED AS DESCRIBED ABOVE.
!   98-06-08  T BLACK      - CONVERSION FROM 1-D TO 2-D
!   98-08-18  MIKE BALDWIN - MODIFY TO COMPUTE RH OVER ICE AS IN MODEL
!   98-12-16  GEOFF MANIKIN - UNDO RH COMPUTATION OVER ICE
!   00-01-04  JIM TUCCILLO - MPI VERSION
!   02-06-11  MIKE BALDWIN - WRF VERSION
!   06-03-19  Wen Meng     - MODIFY TOP PRESSURE to 1 PA
!     
! USAGE:    CALL CALRH(P1,T1,Q1,RH)
!   INPUT ARGUMENT LIST:
!     P1     - PRESSURE (PA)
!     T1     - TEMPERATURE (K)
!     Q1     - SPECIFIC HUMIDITY (KG/KG)
!
!   OUTPUT ARGUMENT LIST: 
!     RH     - RELATIVE HUMIDITY  (DECIMAL FORM)
!     Q1     - ADJUSTED SPECIFIC HUMIDITY (KG/KG)
!     
!   OUTPUT FILES:
!     NONE
!     
!   SUBPROGRAMS CALLED:
!     UTILITIES:
!     LIBRARY:
!       NONE
!     
!   ATTRIBUTES:
!     LANGUAGE: FORTRAN
!     MACHINE : CRAY C-90
!$$$  
!
     use params_mod, only: PQ0, a2, a3, a4, rhmin
     use ctlblk_mod, only: jsta, jend, spval, im
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      implicit none
!     
!     SET PARAMETER.
!
!     DECLARE VARIABLES.
!     
      REAL,dimension(IM,jsta:jend),intent(in)    :: P1,T1
      REAL,dimension(IM,jsta:jend),intent(inout) :: Q1
      REAL,dimension(IM,jsta:jend),intent(out)   :: RH
      REAL QC
      integer I,J
!***************************************************************
!
!     START CALRH.
!
      DO J=JSTA,JEND
        DO I=1,IM
          IF (T1(I,J) < SPVAL) THEN
            IF (ABS(P1(I,J)) >= 1) THEN
              QC = PQ0/P1(I,J)*EXP(A2*(T1(I,J)-A3)/(T1(I,J)-A4))
!
              RH(I,J) = Q1(I,J)/QC
!
!   BOUNDS CHECK
!
              IF (RH(I,J) > 1.0) THEN
                RH(I,J) = 1.0
                Q1(I,J) = RH(I,J)*QC
              ENDIF
              IF (RH(I,J) < RHmin) THEN  !use smaller RH limit for stratosphere
                RH(I,J) = RHmin
                Q1(I,J) = RH(I,J)*QC
              ENDIF
!
            ENDIF
          ELSE
            RH(I,J) = SPVAL
          ENDIF
        ENDDO
      ENDDO
!
!
!      END SUBROUTINE CALRH
      END SUBROUTINE CALRH_NAM
!
!-------------------------------------------------------------------------------------
!
      SUBROUTINE CALRH_GFS(P1,T1,Q1,RH)
!$$$  SUBPROGRAM DOCUMENTATION BLOCK
!                .      .    .     
! SUBPROGRAM:    CALRH       COMPUTES RELATIVE HUMIDITY
!   PRGRMMR: TREADON         ORG: W/NP2      DATE: 92-12-22       
!     
! ABSTRACT:  
!     THIS ROUTINE COMPUTES RELATIVE HUMIDITY GIVEN PRESSURE, 
!     TEMPERATURE, SPECIFIC HUMIDITY. AN UPPER AND LOWER BOUND
!     OF 100 AND 1 PERCENT RELATIVE HUMIDITY IS ENFORCED.  WHEN
!     THESE BOUNDS ARE APPLIED THE PASSED SPECIFIC HUMIDITY 
!     ARRAY IS ADJUSTED AS NECESSARY TO PRODUCE THE SET RELATIVE
!     HUMIDITY.
!   .     
!     
! PROGRAM HISTORY LOG:
!   ??-??-??  DENNIS DEAVEN
!   92-12-22  RUSS TREADON - MODIFIED AS DESCRIBED ABOVE.
!   98-06-08  T BLACK      - CONVERSION FROM 1-D TO 2-D
!   98-08-18  MIKE BALDWIN - MODIFY TO COMPUTE RH OVER ICE AS IN MODEL
!   98-12-16  GEOFF MANIKIN - UNDO RH COMPUTATION OVER ICE
!   00-01-04  JIM TUCCILLO - MPI VERSION
!   02-06-11  MIKE BALDWIN - WRF VERSION
!   13-08-13  S. Moorthi   - Threading
!   06-03-19  Wen Meng     - MODIFY TOP PRESSURE to 1 PA
!     
! USAGE:    CALL CALRH(P1,T1,Q1,RH)
!   INPUT ARGUMENT LIST:
!     P1     - PRESSURE (PA)
!     T1     - TEMPERATURE (K)
!     Q1     - SPECIFIC HUMIDITY (KG/KG)
!
!   OUTPUT ARGUMENT LIST: 
!     RH     - RELATIVE HUMIDITY  (DECIMAL FORM)
!     Q1     - ADJUSTED SPECIFIC HUMIDITY (KG/KG)
!     
!   OUTPUT FILES:
!     NONE
!     
!   SUBPROGRAMS CALLED:
!     UTILITIES:
!     LIBRARY:
!       NONE
!     
!   ATTRIBUTES:
!     LANGUAGE: FORTRAN
!     MACHINE : CRAY C-90
!$$$  
!
      use params_mod, only: rhmin
      use ctlblk_mod, only: jsta, jend, spval, im
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      implicit none
!
      real,parameter:: con_rd      =2.8705e+2 ! gas constant air    (J/kg/K)
      real,parameter:: con_rv      =4.6150e+2 ! gas constant H2O 
      real,parameter:: con_eps     =con_rd/con_rv
      real,parameter:: con_epsm1   =con_rd/con_rv-1
!      real,external::FPVSNEW

      INTERFACE
        ELEMENTAL FUNCTION FPVSNEW (t)
          REAL  FPVSNEW
          REAL, INTENT(IN) :: t
        END FUNCTION FPVSNEW
      END INTERFACE
!
      REAL,dimension(IM,jsta:jend),intent(in)   :: P1,T1
      REAL,dimension(IM,jsta:jend),intent(inout):: Q1,RH
      REAL ES,QC
      integer :: I,J
!***************************************************************
!
!     START CALRH.
!
!$omp parallel do private(i,j,es,qc)
      DO J=JSTA,JEND
        DO I=1,IM
          IF (T1(I,J) < SPVAL .AND. P1(I,J) < SPVAL.AND.Q1(I,J)/=SPVAL) THEN
!           IF (ABS(P1(I,J)) > 1.0) THEN
!            IF (P1(I,J) > 1.0) THEN
            IF (P1(I,J) >= 1.0) THEN
              ES = MIN(FPVSNEW(T1(I,J)),P1(I,J))
              QC = CON_EPS*ES/(P1(I,J)+CON_EPSM1*ES)

!             QC=PQ0/P1(I,J)*EXP(A2*(T1(I,J)-A3)/(T1(I,J)-A4))

              RH(I,J) = min(1.0,max(Q1(I,J)/QC,rhmin))
              q1(i,j) = rh(i,j)*qc

!   BOUNDS CHECK
!
!             IF (RH(I,J) > 1.0) THEN
!               RH(I,J) = 1.0
!               Q1(I,J) = RH(I,J)*QC
!             ELSEIF (RH(I,J) < RHmin) THEN  !use smaller RH limit for stratosphere
!               RH(I,J) = RHmin
!               Q1(I,J) = RH(I,J)*QC
!             ENDIF

            ENDIF
          ELSE
            RH(I,J) = SPVAL
          ENDIF
        ENDDO
      ENDDO

      END SUBROUTINE CALRH_GFS
!
!-------------------------------------------------------------------------------------
!      
      SUBROUTINE CALRH_GSD(P1,T1,Q1,RHB)
!
! Algorithm use at GSD for RUC and Rapid Refresh                           
!------------------------------------------------------------------
!

      use ctlblk_mod, only: jsta, jend, im

      implicit none

      integer :: j, i
      real :: tx, pol, esx, es, e
      real, dimension(im,jsta:jend) :: P1, T1, Q1, RHB


      DO J=JSTA,JEND
        DO I=1,IM

! - compute relative humidity
          Tx=T1(I,J)-273.15
          POL = 0.99999683       + TX*(-0.90826951E-02 +    &
             TX*(0.78736169E-04   + TX*(-0.61117958E-06 +   &
             TX*(0.43884187E-08   + TX*(-0.29883885E-10 +   &
             TX*(0.21874425E-12   + TX*(-0.17892321E-14 +   &
             TX*(0.11112018E-16   + TX*(-0.30994571E-19)))))))))
          esx = 6.1078/POL**8

          ES = esx
          E = P1(I,J)/100.*Q1(I,J)/(0.62197+Q1(I,J)*0.37803)
          RHB(I,J) = MIN(1.,E/ES)

        ENDDO
      ENDDO

      END SUBROUTINE CALRH_GSD
!
!-------------------------------------------------------------------------------------
!
      SUBROUTINE CALRH_PW(RHPW)
!
! Algorithm use at GSD for RUC and Rapid Refresh                           
!------------------------------------------------------------------
!

      use vrbls3d, only: q, pmid, t
      use params_mod, only: g
      use ctlblk_mod, only: lm, jsta, jend, lm, im
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
       implicit none

      real,PARAMETER :: svp1=6.1153,svp2=17.67,svp3=29.65

      REAL, dimension(im,jsta:jend):: PW, PW_SAT, RHPW
      REAL deltp,sh,qv,temp,es,qs,qv_sat
      integer i,j,l,k,ka,kb

      pw     = 0.
      pw_sat = 0.
      rhpw   = 0.

      DO L=1,LM
        k=lm-l+1
       DO J=JSTA,JEND
        DO I=1,IM
! -- use specific humidity for PW calculation
           sh = q(i,j,k)
           qv = sh/(1.-sh)
           KA = MAX(1,K-1)
           KB = MIN(LM,K+1)

!   assumes that P is in mb at this point - be careful!
           DELTP = 0.5*(PMID(I,J,KB)-PMID(I,J,KA))
           PW(I,J) = PW(I,J) + sh *DELTP/G

!Csgb -- Add more for RH w.r.t. PW-sat

          temp = T(I,J,K)
! --- use saturation mixing ratio w.r.t. water here
!       for this check.
          es = svp1*exp(SVP2*(Temp-273.15)/(Temp-SVP3))
! -- get saturation specific humidity (w.r.t. total air)
          qs = 0.62198*es/(pmid(i,j,k)*1.e-2-0.37802*es)
! -- get saturation mixing ratio (w.r.t. dry air)
          qv_sat = qs/(1.-qs)

          pw_sat(i,j) = pw_sat(i,j) + max(sh,Qs)*DELTP/G

        if (i.eq.120 .and. j.eq.120 )                        &
          write (6,*)'pw-sat', temp, sh, qs, pmid(i,j,kb)    &
          ,pmid(i,j,ka),pw(i,j),pw_sat(i,j)

!sgb - This IS RH w.r.t. PW-sat.
           RHPW (i,j) = min(1.,PW(i,j) / pw_sat(i,j)) * 100.

        ENDDO
       ENDDO
      ENDDO

      END SUBROUTINE CALRH_PW
!
!-------------------------------------------------------------------------------------
!
      END MODULE CALRH_MODULE

