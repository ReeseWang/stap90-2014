      SUBROUTINE QUADS (NEL,ITYPE,NINT,THIC,YM,PR,XX,S,IOUT)            QUA00001
C                                                                       QUA00002
C . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . QUA00003
C .                                                                   . QUA00004
C .   P R O G R A M                                                   . QUA00005
C .        TO CALCULATE ISOPARAMETRIC QUADRILATERAL ELEMENT STIFFNESS . QUA00006
C .        MATRIX FOR AXISYMMETRIC, PLANE STRESS, AND PLANE STRAIN    . QUA00007
C .        CONDITIONS                                                 . QUA00008
C .                                                                   . QUA00009
C .  - - INPUT VARIABLES - -                                          . QUA00010
C .        NEL       = NUMBER OF ELEMENT                              . QUA00011
C .        ITYPE     = ELEMENT TYPE                                   . QUA00012
C .                        EQ.0 = AXISYMMETRIC                        . QUA00013
C .                        EQ.1 = PLANE STRAIN                        . QUA00014
C .                        EQ.2 = PLANE STRESS                        . QUA00015
C .        NINT      = GAUSS NUMERICAL INTEGRATION ORDER              . QUA00016
C .        THIC      = THICKNESS OF ELEMENT                           . QUA00017
C .        YM        = YOUNG'S MODULUS                                . QUA00018
C .        PR        = POISSON'S RATIO                                . QUA00019
C .        XX(2,4)   = ELEMENT NODE COORDINATES                       . QUA00020
C .        S(8,8)    = STORAGE FOR STIFFNESS MATRIX                   . QUA00021
C .        IOUT      = UNIT NUMBER USED FOR OUTPUT                    . QUA00022
C .                                                                   . QUA00023
C .  - - OUTPUT - -                                                   . QUA00024
C .        S(8,8)    = CALCULATED STIFFNESS MATRIX                    . QUA00025
C .                                                                   . QUA00026
C . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . QUA00027
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)                               QUA00028
C . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . QUA00029
C .   THIS PROGRAM IS USED IN SINGLE PRECISION ARITHMETIC ON CRAY     . QUA00030
C .   EQUIPMENT AND DOUBLE PRECISION ARITHMETIC ON IBM MACHINES,      . QUA00031
C .   ENGINEERING WORKSTATIONS AND PCS. DEACTIVATE ABOVE LINE FOR     . QUA00032
C .   SINGLE PRECISION ARITHMETIC.                                    . QUA00033
C . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . QUA00034
      DIMENSION D(4,4),B(4,8),XX(2,4),S(8,8),XG(4,4),WGT(4,4),DB(4)     QUA00035
C                                                                       QUA00036
C     MATRIX XG STORES GAUSS - LEGENDRE SAMPLING POINTS                 QUA00037
C                                                                       QUA00038
      DATA XG/   0.D0,   0.D0,   0.D0,   0.D0,   -.5773502691896D0,     QUA00039
     1 .5773502691896D0,   0.D0,   0.D0,   -.7745966692415D0,   0.D0,   QUA00040
     2 .7745966692415D0,   0.D0,   -.8611363115941D0,                   QUA00041
     3 -.3399810435849D0,   .3399810435849D0,   .8611363115941D0 /      QUA00042
C                                                                       QUA00043
C     MATRIX WGT STORES GAUSS - LEGENDRE WEIGHTING FACTORS              QUA00044
C                                                                       QUA00045
      DATA WGT /  2.D0,   0.D0,   0.D0,   0.D0,   1.D0,   1.D0,         QUA00046
     1 0.D0,   0.D0,   .5555555555556D0,   .8888888888889D0,            QUA00047
     2 .5555555555556D0,   0.D0,   .3478548451375D0,   .6521451548625D0,QUA00048
     3 .6521451548625D0,   .3478548451375D0 /                           QUA00049
C                                                                       QUA00050
C     O B T A I N  S T R E S S - S T R A I N  L A W                     QUA00051
C                                                                       QUA00052
      F=YM/(1.+PR)                                                      QUA00053
      G=F*PR/(1.-2.*PR)                                                 QUA00054
      H=F + G                                                           QUA00055
C                                                                       QUA00056
C     PLANE STRAIN ANALYSIS                                             QUA00057
C                                                                       QUA00058
      D(1,1)=H                                                          QUA00059
      D(1,2)=G                                                          QUA00060
      D(1,3)=0.                                                         QUA00061
      D(2,1)=G                                                          QUA00062
      D(2,2)=H                                                          QUA00063
      D(2,3)=0.                                                         QUA00064
      D(3,1)=0.                                                         QUA00065
      D(3,2)=0.                                                         QUA00066
      D(3,3)=F/2.                                                       QUA00067
      IF (ITYPE.EQ.1) THEN                                              QUA00068
      THIC=1.                                                           QUA00069
      GO TO 20                                                          QUA00070
      ENDIF                                                             QUA00071
C                                                                       QUA00072
C     AXISYMMETRIC ANALYSIS                                             QUA00073
C                                                                       QUA00074
      D(1,4)=G                                                          QUA00075
      D(2,4)=G                                                          QUA00076
      D(3,4)=0.                                                         QUA00077
      D(4,1)=G                                                          QUA00078
      D(4,2)=G                                                          QUA00079
      D(4,3)=0.                                                         QUA00080
      D(4,4)=H                                                          QUA00081
      IF (ITYPE.EQ.0) GO TO 20                                          QUA00082
C                                                                       QUA00083
C     FOR PLANE STRESS ANALYSIS CONDENSE STRESS-STRAIN MATRIX           QUA00084
C                                                                       QUA00085
      DO 10 I=1,3                                                       QUA00086
      A=D(I,4)/D(4,4)                                                   QUA00087
      DO 10 J=I,3                                                       QUA00088
      D(I,J)=D(I,J) - D(4,J)*A                                          QUA00089
   10 D(J,I)=D(I,J)                                                     QUA00090
C                                                                       QUA00091
C     C A L C U L A T E  E L E M E N T  S T I F F N E S S               QUA00092
C                                                                       QUA00093
   20 DO 30 I=1,8                                                       QUA00094
      DO 30 J=1,8                                                       QUA00095
   30 S(I,J)=0.                                                         QUA00096
      IST=3                                                             QUA00097
      IF (ITYPE.EQ.0) IST=4                                             QUA00098
      DO 80 LX=1,NINT                                                   QUA00099
      RI=XG(LX,NINT)                                                    QUA00100
      DO 80 LY=1,NINT                                                   QUA00101
      SI=XG(LY,NINT)                                                    QUA00102
C                                                                       QUA00103
C     EVALUATE DERIVATIVE OPERATOR B AND THE JACOBIAN DETERMINANT DET   QUA00104
C                                                                       QUA00105
      CALL STDM (XX,B,DET,RI,SI,XBAR,NEL,ITYPE,IOUT)                    QUA00106
C                                                                       QUA00107
C     ADD CONTRIBUTION TO ELEMENT STIFFNESS                             QUA00108
C                                                                       QUA00109
      IF (ITYPE.GT.0) XBAR=THIC                                         QUA00110
      WT=WGT(LX,NINT)*WGT(LY,NINT)*XBAR*DET                             QUA00111
      DO 70 J=1,8                                                       QUA00112
      DO 40 K=1,IST                                                     QUA00113
      DB(K)=0.0                                                         QUA00114
      DO 40 L=1,IST                                                     QUA00115
   40 DB(K)=DB(K) + D(K,L)*B(L,J)                                       QUA00116
      DO 60 I=J,8                                                       QUA00117
      STIFF=0.0                                                         QUA00118
      DO 50 L=1,IST                                                     QUA00119
   50 STIFF=STIFF + B(L,I)*DB(L)                                        QUA00120
   60 S(I,J)=S(I,J) + STIFF*WT                                          QUA00121
   70 CONTINUE                                                          QUA00122
   80 CONTINUE                                                          QUA00123
C                                                                       QUA00124
      DO 90 J=1,8                                                       QUA00125
      DO 90 I=J,8                                                       QUA00126
   90 S(J,I)=S(I,J)                                                     QUA00127
C                                                                       QUA00128
      RETURN                                                            QUA00129
C                                                                       QUA00130
      END                                                               QUA00131
      SUBROUTINE STDM (XX,B,DET,R,S,XBAR,NEL,ITYPE,IOUT)                QUA00132
C                                                                       QUA00133
C . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . QUA00134
C .                                                                   . QUA00135
C .   P R O G R A M                                                   . QUA00136
C .     TO EVALUATE THE STRAIN-DISPLACEMENT TRANSFORMATION MATRIX B   . QUA00137
C .     AT POINT (R,S) FOR A QUADRILATERAL ELEMENT                    . QUA00138
C .                                                                   . QUA00139
C . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . QUA00140
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)                               QUA00141
      DIMENSION XX(2,4),B(4,8),H(4),P(2,4),XJ(2,2),XJI(2,2)             QUA00142
C                                                                       QUA00143
      RP = 1.0 + R                                                      QUA00144
      SP = 1.0 + S                                                      QUA00145
      RM = 1.0 - R                                                      QUA00146
      SM = 1.0 - S                                                      QUA00147
C                                                                       QUA00148
C     INTERPOLATION FUNCTIONS                                           QUA00149
C                                                                       QUA00150
      H(1) = 0.25* RP* SP                                               QUA00151
      H(2) = 0.25* RM* SP                                               QUA00152
      H(3) = 0.25* RM* SM                                               QUA00153
      H(4) = 0.25* RP* SM                                               QUA00154
C                                                                       QUA00155
C     NATURAL COORDINATE DERIVATIVES OF THE INTERPOLATION FUNCTIONS     QUA00156
C                                                                       QUA00157
C        1. WITH RESPECT TO R                                           QUA00158
C                                                                       QUA00159
      P(1,1) = 0.25* SP                                                 QUA00160
      P(1,2) = - P(1,1)                                                 QUA00161
      P(1,3) = - 0.25* SM                                               QUA00162
      P(1,4) = - P(1,3)                                                 QUA00163
C                                                                       QUA00164
C        2. WITH RESPECT TO S                                           QUA00165
C                                                                       QUA00166
      P(2,1) = 0.25* RP                                                 QUA00167
      P(2,2) = 0.25* RM                                                 QUA00168
      P(2,3) = - P(2,2)                                                 QUA00169
      P(2,4) = - P(2,1)                                                 QUA00170
C                                                                       QUA00171
C     EVALUATE THE JACOBIAN MATRIX AT POINT (R,S)                       QUA00172
C                                                                       QUA00173
   10 DO 30 I=1,2                                                       QUA00174
      DO 30 J=1,2                                                       QUA00175
      DUM = 0.0                                                         QUA00176
      DO 20 K=1,4                                                       QUA00177
   20 DUM=DUM + P(I,K)*XX(J,K)                                          QUA00178
   30 XJ(I,J)=DUM                                                       QUA00179
C                                                                       QUA00180
C     COMPUTE THE DETERMINANT OF THE JACOBIAN MATRIX AT POINT (R,S)     QUA00181
C                                                                       QUA00182
      DET = XJ(1,1)* XJ(2,2) - XJ(2,1)* XJ(1,2)                         QUA00183
      IF (DET.GT.0.00000001) GO TO 40                                   QUA00184
      WRITE (IOUT,2000) NEL                                             QUA00185
      GO TO 800                                                         QUA00186
C                                                                       QUA00187
C     COMPUTE INVERSE OF THE JACOBIAN MATRIX                            QUA00188
C                                                                       QUA00189
   40 DUM=1./DET                                                        QUA00190
      XJI(1,1) = XJ(2,2)* DUM                                           QUA00191
      XJI(1,2) =-XJ(1,2)* DUM                                           QUA00192
      XJI(2,1) =-XJ(2,1)* DUM                                           QUA00193
      XJI(2,2) = XJ(1,1)* DUM                                           QUA00194
C                                                                       QUA00195
C     EVALUATE GLOBAL DERIVATIVE OPERATOR B                             QUA00196
C                                                                       QUA00197
      K2=0                                                              QUA00198
      DO 60 K=1,4                                                       QUA00199
      K2=K2 + 2                                                         QUA00200
      B(1,K2-1) = 0.                                                    QUA00201
      B(1,K2  ) = 0.                                                    QUA00202
      B(2,K2-1) = 0.                                                    QUA00203
      B(2,K2  ) = 0.                                                    QUA00204
      DO 50 I=1,2                                                       QUA00205
      B(1,K2-1) = B(1,K2-1) + XJI(1,I) * P(I,K)                         QUA00206
   50 B(2,K2  ) = B(2,K2  ) + XJI(2,I) * P(I,K)                         QUA00207
      B(3,K2  ) = B(1,K2-1)                                             QUA00208
   60 B(3,K2-1) = B(2,K2  )                                             QUA00209
C                                                                       QUA00210
C     IN CASE OF PLANE STRAIN OR PLANE STRESS ANALYSIS DO NOT INCLUDE   QUA00211
C     THE NORMAL STRAIN COMPONENT                                       QUA00212
C                                                                       QUA00213
      IF (ITYPE.GT.0) GO TO 900                                         QUA00214
C                                                                       QUA00215
C     COMPUTE THE RADIUS AT POINT (R,S)                                 QUA00216
C                                                                       QUA00217
      XBAR=0.0                                                          QUA00218
      DO 70 K=1,4                                                       QUA00219
   70 XBAR=XBAR + H(K)*XX(1,K)                                          QUA00220
C                                                                       QUA00221
C     EVALUATE THE HOOP STRAIN-DISPLACEMENT RELATION                    QUA00222
C                                                                       QUA00223
      IF (XBAR.GT.0.00000001) GO TO 90                                  QUA00224
C                                                                       QUA00225
C     FOR THE CASE OF ZERO RADIUS EQUATE RADIAL TO HOOP STRAIN          QUA00226
C                                                                       QUA00227
      DO 80 K=1,8                                                       QUA00228
   80 B(4,K)=B(1,K)                                                     QUA00229
      GO TO 900                                                         QUA00230
C                                                                       QUA00231
C     NON-ZERO RADIUS                                                   QUA00232
C                                                                       QUA00233
   90 DUM=1./XBAR                                                       QUA00234
      K2=0                                                              QUA00235
      DO 100 K=1,4                                                      QUA00236
      K2=K2 + 2                                                         QUA00237
      B(4,K2  ) = 0.                                                    QUA00238
  100 B(4,K2-1) = H(K)*DUM                                              QUA00239
      GO TO 900                                                         QUA00240
C                                                                       QUA00241
  800 STOP                                                              QUA00242
  900 RETURN                                                            QUA00243
C                                                                       QUA00244
 2000 FORMAT (//,' *** ERROR *** ',                                     QUA00245
     1    ' ZERO OR NEGATIVE JACOBIAN DETERMINANT FOR ELEMENT (',I8,')')QUA00246
C                                                                       QUA00247
      END                                                               QUA00248