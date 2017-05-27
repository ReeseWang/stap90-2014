! . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
! .                                                                       .
! .                            S T A P 9 0                                .
! .                                                                       .
! .     AN IN-CORE SOLUTION STATIC ANALYSIS PROGRAM IN FORTRAN 90         .
! .     Adapted from STAP (KJ Bath, FORTRAN IV) for teaching purpose      .
! .                                                                       .
! .     Xiong Zhang, (2013)                                               .
! .     Computational Dynamics Group, School of Aerospace                 .
! .     Tsinghua Univerity                                                .
! .                                                                       .
! . . . . . . . . . . . . . .  . . .  . . . . . . . . . . . . . . . . . . .

SUBROUTINE SOLID
! . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
! .                                                                   .
! .   To set up storage and call the SOLID element subroutine         .
! .                                                                   .
! . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

  USE GLOBALS
  USE MEMALLOCATE

  IMPLICIT NONE
  INTEGER :: NUME, NUMMAT, MM, N101, N102, N103, N104, N105, N106, N107

  NUME = NPAR(2)
  NUMMAT = NPAR(3)

! Allocate storage for element group data
  IF (IND == 1) THEN
      MM = 3*NUMMAT*ITWO + 25*NUME + 24*NUME*ITWO
      CALL MEMALLOC(11,"ELEGP",MM,1)
  END IF

  NFIRST=NP(11)   ! Pointer to the first entry in the element group data array
                  ! in the unit of single precision (corresponding to A)

! Calculate the pointer to the arrays in the element group data
! N101: E(NUMMAT)
! N102: RHO(NUMMAT)
! N103: NIU(NUMMAT)
! N104: LM(24,NUME)
! N105: XYZ(24,NUME)
! N106: MTAP(NUME)
  N101=NFIRST
  N102=N101+NUMMAT*ITWO
  N103=N102+NUMMAT*ITWO
  N104=N103+NUMMAT*ITWO
  N105=N104+24*NUME
  N106=N105+24*NUME*ITWO
  N107=N106+NUME
  NLAST=N107

  MIDEST=NLAST - NFIRST

  CALL OLID (IA(NP(1)),DA(NP(2)),DA(NP(3)),DA(NP(4)),DA(NP(4)),IA(NP(5)),   &
       A(N101),A(N102),A(N103),A(N104),A(N105),A(N106))

  RETURN

END SUBROUTINE SOLID


SUBROUTINE OLID (ID,X,Y,Z,U,MHT,E,RHO,NIU,LM,XYZ,MATP)
! . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
! .                                                                   .
! .   SOLID element subroutine                                        .
! .                                                                   .
! . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

  USE GLOBALS
  USE MEMALLOCATE

  IMPLICIT NONE
  INTEGER :: ID(3,NUMNP),LM(24,NPAR(2)),MATP(NPAR(2)),MHT(NEQ)
  REAL(8) :: X(NUMNP),Y(NUMNP),Z(NUMNP),E(NPAR(3)),RHO(NPAR(3)),NIU(NPAR(3)),  &
             XYZ(24,NPAR(2)),U(NEQ)
  REAL(8) :: S(24,24),B(6,24),D(6,6),ST(6),STR(6)

  INTEGER :: NPAR1, NUME, NUMMAT, ND, I1, I2, I3, I4, I5, I6, I7, I8, L, N, I, J
  INTEGER :: MTYPE, IPRINT
  REAL(8) :: X0, Y0, Z0, XL, YL, ZL, V, XI(8), ITA(8), ZETA(8), CO, DNX(8), DNY(8), DNZ(8)

  NPAR1  = NPAR(1)
  NUME   = NPAR(2)
  NUMMAT = NPAR(3) 

  ND=24

! Read and generate element information
  IF (IND .EQ. 1) THEN

     WRITE (IOUT,"(' E L E M E N T   D E F I N I T I O N',//,  &
                   ' ELEMENT TYPE ',13(' .'),'( NPAR(1) ) . . =',I5,/,   &
                   '     EQ.1, TRUSS ELEMENTS',/,      &
                   '     EQ.8, SOLID ELEMENTS',/,  &
                   '     EQ.OTHERS, NOT AVAILABLE',//,      &
                   ' NUMBER OF ELEMENTS.',10(' .'),'( NPAR(2) ) . . =',I5,/)") NPAR1,NUME

     IF (NUMMAT.EQ.0) NUMMAT=1

     WRITE (IOUT,"(' M A T E R I A L   D E F I N I T I O N',//,  &
                   ' NUMBER OF DIFFERENT SETS OF MATERIAL',  &
                    (' .'),'( NPAR(3) ) . . =',I5,/)") NUMMAT

     WRITE (IOUT,"('  SET       YOUNG''S         DENSITY         POSSION''S',/,  &
                   ' NUMBER     MODULUS',27X,'RATIO',/,  &
                   15 X,'E',14X,'RHO',14X,'NIU')")

     DO I=1,NUMMAT
        READ (IIN,'(I5,3F10.0)') N,E(N),RHO(N),NIU(N)  ! Read material information
        WRITE (IOUT,"(I5,4X,E12.5,4X,E12.5,4X,E12.5)") N,E(N),RHO(N),NIU(N)
     END DO

     WRITE (IOUT,"(//,' E L E M E N T   I N F O R M A T I O N',//,  &
                      ' ELEMENT     NODE     NODE     NODE     NODE     NODE     NODE     NODE     NODE       MATERIAL',/,   &
                      ' NUMBER-N     I1       I2       I3       I4       I5       I6       I7       I8       SET NUMBER')")

     N=0
     DO WHILE (N .NE. NUME)
        READ (IIN,'(10I5)') N,I1,I2,I3,I4,I5,I6,I7,I8,MTYPE  ! Read in element information

!       Save element information
        XYZ(1,N)=X(I1)  ! Coordinates of the element's 1st node
        XYZ(2,N)=Y(I1)
        XYZ(3,N)=Z(I1)
        XYZ(4,N)=X(I2)  ! Coordinates of the element's 2nd node
        XYZ(5,N)=Y(I2)
        XYZ(6,N)=Z(I2)
        XYZ(7,N)=X(I3)  ! Coordinates of the element's 3rd node
        XYZ(8,N)=Y(I3)
        XYZ(9,N)=Z(I3)
        XYZ(10,N)=X(I4)  ! Coordinates of the element's 4th node
        XYZ(11,N)=Y(I4)
        XYZ(12,N)=Z(I4)
        XYZ(13,N)=X(I5)  ! Coordinates of the element's 5th node
        XYZ(14,N)=Y(I5)
        XYZ(15,N)=Z(I5)
        XYZ(16,N)=X(I6)  ! Coordinates of the element's 6th node
        XYZ(17,N)=Y(I6)
        XYZ(18,N)=Z(I6)
        XYZ(19,N)=X(I7)  ! Coordinates of the element's 7th node
        XYZ(20,N)=Y(I7)
        XYZ(21,N)=Z(I7)
        XYZ(22,N)=X(I8)  ! Coordinates of the element's 8th node
        XYZ(23,N)=Y(I8)
        XYZ(24,N)=Z(I8)

        MATP(N)=MTYPE  ! Material type

        DO L=1,24
           LM(L,N)=0
        END DO

        DO L=1,3
           LM(L,N)=ID(L,I1)     ! Connectivity matrix
           LM(L+3,N)=ID(L,I2)
           LM(L+6,N)=ID(L,I3)
           LM(L+9,N)=ID(L,I4)
           LM(L+12,N)=ID(L,I5)
           LM(L+15,N)=ID(L,I6)
           LM(L+18,N)=ID(L,I7)
           LM(L+21,N)=ID(L,I8)
        END DO

!       Update column heights and bandwidth
        CALL COLHT (MHT,ND,LM(1,N))   

        WRITE (IOUT,"(I5,6X,I5,4X,I5,4X,I5,4X,I5,4X,I5,4X,I5,4X,I5,4X,I5,7X,I5)") N,I1,I2,I3,I4,I5,I6,I7,I8,MTYPE

     END DO

     RETURN

! Assemble stucture stiffness matrix
  ELSE IF (IND .EQ. 2) THEN

     DO N=1,NUME
        MTYPE=MATP(N)
        
        X0 = (XYZ(1,N)+XYZ(4,N))/2.
        Y0 = (XYZ(2,N)+XYZ(11,N))/2.
        Z0 = (XYZ(3,N)+XYZ(15,N))/2.
        XL = ABS(XYZ(4,N)-XYZ(1,N))
        YL = ABS(XYZ(11,N)-XYZ(2,N))
        ZL = ABS(XYZ(15,N)-XYZ(3,N))
        V = XL*YL*ZL
        
        DO L = 1,8
            XI(L) = (XYZ(1+(L-1)*3,N)-X0)/XL*2.
            ITA(L) = (XYZ(2+(L-1)*3,N)-Y0)/YL*2.
            ZETA(L) = (XYZ(3+(L-1)*3,N)-Z0)/ZL*2.
        END DO
        
        CO = E(MTYPE)*V/(1+NIU(MTYPE))/(1-2*NIU(MTYPE))/16.
        DO I = 1,8
            DO J = 1,8
                S(1+(I-1)*3, 1+(J-1)*3) = (1-NIU(MTYPE))*XI(I)*XI(J)*(1+ITA(I)*ITA(J)/3.)*(1+ZETA(I)*ZETA(J)/3.)/XL/XL + &
                    & (1-2*NIU(MTYPE))*(ITA(I)*ITA(J)*(1+ZETA(I)*ZETA(J)/3.)*(1+XI(I)*XI(J)/3.)/YL/YL + &
                    & ZETA(I)*ZETA(J)*(1+XI(I)*XI(J)/3.)*(1+ITA(I)*ITA(J)/3.)/ZL/ZL)/2.
                S(2+(I-1)*3, 2+(J-1)*3) = (1-NIU(MTYPE))*ITA(I)*ITA(J)*(1+ZETA(I)*ZETA(J)/3.)*(1+XI(I)*XI(J)/3.)/YL/YL + &
                    & (1-2*NIU(MTYPE))*(XI(I)*XI(J)*(1+ITA(I)*ITA(J)/3.)*(1+ZETA(I)*ZETA(J)/3.)/XL/XL + &
                    & ZETA(I)*ZETA(J)*(1+XI(I)*XI(J)/3.)*(1+ITA(I)*ITA(J)/3.)/ZL/ZL)/2.
                S(3+(I-1)*3, 3+(J-1)*3) = (1-NIU(MTYPE))*ZETA(I)*ZETA(J)*(1+XI(I)*XI(J)/3.)*(1+ITA(I)*ITA(J)/3.)/ZL/ZL + &
                    & (1-2*NIU(MTYPE))*(XI(I)*XI(J)*(1+ITA(I)*ITA(J)/3.)*(1+ZETA(I)*ZETA(J)/3.)/XL/XL + &
                    & ITA(I)*ITA(J)*(1+XI(I)*XI(J)/3.)*(1+ZETA(I)*ZETA(J)/3.)/YL/YL)/2.
                S(1+(I-1)*3, 2+(J-1)*3) = (1+ZETA(I)*ZETA(J)/3.)*(NIU(MTYPE)*XI(I)*ITA(J)+(1-2*NIU(MTYPE))*ITA(I)*XI(J)/2.)/XL/YL
                S(1+(I-1)*3, 3+(J-1)*3) = (1+ITA(I)*ITA(J)/3.)*(NIU(MTYPE)*XI(I)*ZETA(J)+(1-2*NIU(MTYPE))*ZETA(I)*XI(J)/2.)/XL/ZL
                S(2+(I-1)*3, 1+(J-1)*3) = (1+ZETA(I)*ZETA(J)/3.)*(NIU(MTYPE)*ITA(I)*XI(J)+(1-2*NIU(MTYPE))*XI(I)*ITA(J)/2.)/XL/YL
                S(2+(I-1)*3, 3+(J-1)*3) = (1+XI(I)*XI(J)/3.)*(NIU(MTYPE)*ITA(I)*ZETA(J)+(1-2*NIU(MTYPE))*ZETA(I)*ITA(J)/2.)/YL/ZL
                S(3+(I-1)*3, 1+(J-1)*3) = (1+ITA(I)*ITA(J)/3.)*(NIU(MTYPE)*ZETA(I)*XI(J)+(1-2*NIU(MTYPE))*XI(I)*ZETA(J)/2.)/XL/ZL
                S(3+(I-1)*3, 2+(J-1)*3) = (1+XI(I)*XI(J)/3.)*(NIU(MTYPE)*ZETA(I)*ITA(J)+(1-2*NIU(MTYPE))*ITA(I)*ZETA(J)/2.)/YL/ZL
            END DO
        END DO
        
        DO I = 1,24
            DO J = 1,24
                S(I,J) = CO*S(I,J)
            END DO
        END DO
        
!        ! TEST
!        DO I = 1,24
!            WRITE(IOUT, "(24E20.8)") S(I,:)
!        END DO
        


!        XL2=0.
!        DO L=1,3
!           D(L)=XYZ(L,N) - XYZ(L+3,N)
!           XL2=XL2 + D(L)*D(L)
!        END DO
!        XL=SQRT(XL2)   ! Length of element N

!        XX=E(MTYPE)*AREA(MTYPE)*XL   !  E*A*l
!        XX = 0
!        DO L=1,3
!           ST(L)=D(L)/XL2
!           ST(L+3)=-ST(L)
!        END DO
!
!        DO J=1,ND
!           YY=ST(J)*XX
!           DO I=1,J
!              S(I,J)=ST(I)*YY
!           END DO
!        END DO

        CALL ADDBAN (DA(NP(3)),IA(NP(2)),S,LM(1,N),ND)

     END DO

     RETURN

! Stress calculations
  ELSE IF (IND .EQ. 3) THEN

     IPRINT=0
     DO N=1,NUME
        IPRINT=IPRINT + 1
        IF (IPRINT.GT.50) IPRINT=1
        IF (IPRINT.EQ.1) WRITE (IOUT,"(//,' S T R E S S  C A L C U L A T I O N S  F O R  ',  &
                                           'E L E M E N T  G R O U P',I4,//,   &
                                           '  ELEMENT',47X,'STRESS',/,'  NUMBER',8X,'Sxx',13X,'Syy',14X,'Szz',14X,'Sxy',14X,'Syz',14X,'Szx')") NG
        MTYPE=MATP(N)
        
        X0 = (XYZ(1,N)+XYZ(4,N))/2.
        Y0 = (XYZ(2,N)+XYZ(11,N))/2.
        Z0 = (XYZ(3,N)+XYZ(15,N))/2.
        XL = ABS(XYZ(4,N)-XYZ(1,N))
        YL = ABS(XYZ(11,N)-XYZ(2,N))
        ZL = ABS(XYZ(15,N)-XYZ(3,N))
        
        DO L = 1,8
            XI(L) = (XYZ(1+(L-1)*3,N)-X0)/XL*2.
            ITA(L) = (XYZ(2+(L-1)*3,N)-Y0)/YL*2.
            ZETA(L) = (XYZ(3+(L-1)*3,N)-Z0)/ZL*2.
        END DO

        
        DO I = 1,8
            DNX(I) = XI(I)/XL/4.
            DNY(I) = ITA(I)/YL/4.
            DNZ(I) = ZETA(I)/ZL/4.
            B(1,1+(I-1)*3) = DNX(I)
            B(2,2+(I-1)*3) = DNY(I)
            B(3,3+(I-1)*3) = DNZ(I)
            B(1,2+(I-1)*3) = 0
            B(1,3+(I-1)*3) = 0
            B(2,3+(I-1)*3) = 0
            B(2,1+(I-1)*3) = 0
            B(3,1+(I-1)*3) = 0
            B(3,2+(I-1)*3) = 0
            B(4,1+(I-1)*3) = DNY(I)
            B(4,2+(I-1)*3) = DNX(I)
            B(4,3+(I-1)*3) = 0
            B(5,1+(I-1)*3) = 0
            B(5,2+(I-1)*3) = DNZ(I)
            B(5,3+(I-1)*3) = DNY(I)
            B(6,1+(I-1)*3) = DNZ(I)
            B(6,2+(I-1)*3) = 0
            B(6,3+(I-1)*3) = DNX(I)
        END DO
        
        DO I = 1,6
            ST(I) = 0.
            DO J = 1,24
                L = LM(J,N)
                IF (L.GT.0) ST(I) = ST(I) + B(I,J)*U(L)
            END DO
        END DO
        
        DO I = 1,6
            DO J = 1,6
                D(I,J) = 0
            END DO
        END DO
        D(1,1) = 1
        D(2,2) = 1
        D(3,3) = 1
        D(1,2) = NIU(MTYPE)/(1-NIU(MTYPE))
        D(1,3) = NIU(MTYPE)/(1-NIU(MTYPE))
        D(2,3) = NIU(MTYPE)/(1-NIU(MTYPE))
        D(2,1) = NIU(MTYPE)/(1-NIU(MTYPE))
        D(3,1) = NIU(MTYPE)/(1-NIU(MTYPE))
        D(3,2) = NIU(MTYPE)/(1-NIU(MTYPE))
        D(4,4) = (1-2*NIU(MTYPE))/(1-NIU(MTYPE))/2.
        D(5,5) = (1-2*NIU(MTYPE))/(1-NIU(MTYPE))/2.
        D(6,6) = (1-2*NIU(MTYPE))/(1-NIU(MTYPE))/2.
        DO I = 1,6
            DO J = 1,6
                D(I,J) = D(I,J)*E(MTYPE)*(1-NIU(MTYPE))/(1+NIU(MTYPE))/(1-2*NIU(MTYPE))
            END DO
        END DO
        
        DO I = 1,6
            STR(I) = 0.
            DO J = 1,6
                STR(I) = STR(I) + D(I,J)*ST(J)
            END DO
        END DO
        

!        XL2=0.
!        DO L=1,3
!           D(L) = XYZ(L,N) - XYZ(L+3,N)
!           XL2=XL2 + D(L)*D(L)
!        END DO
!
!        DO L=1,3
!           ST(L)=(D(L)/XL2)*E(MTYPE)
!           ST(L+3)=-ST(L)
!        END DO

!        STR=0.0
!        DO L=1,3
!           I=LM(L,N)
!           IF (I.GT.0) STR=STR + ST(L)*U(I)
!
!           J=LM(L+3,N)
!           IF (J.GT.0) STR=STR + ST(L+3)*U(J)
!        END DO

!        P=STR*AREA(MTYPE)

        WRITE (IOUT,"(1X,I5,4X,E13.6,4X,E13.6,4X,E13.6,4X,E13.6,4X,E13.6,4X,E13.6)") N,STR(:)
     END DO

  ELSE 
     STOP "*** ERROR *** Invalid IND value."
  END IF

END SUBROUTINE OLID
