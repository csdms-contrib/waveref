************************************************************************
*       WAVE REFRACTION PROGRAM                                        *
************************************************************************
      Implicit None
      Character*80 formt
      Integer*2 LI,LJ,i,j,n,im,jm,nm,sm,dm,st,p,q,k,iread
      PARAMETER (LI=50,LJ=50)
C-----------------------------------------------------------------------
C  LI and LJ are the dimensions of the grid; N=im*jm
C-----------------------------------------------------------------------
      Real*8 c(LI,LJ),d(LI,LJ),dt,xo(LI),yo(LI),x(500),y(500),
     #ds,dx(LI),dy(LJ),cx(LI,LJ),cy(LI,LJ),g,pi,thetao(500),theta,
     #sumx,sumy,Lo,L,mn,thetabar,a1,a2,a3,a4,t,dtheta,co,tm,A,xx,yy,
     #cc,ks,cxx,cyy,xhold(500,200),yhold(500,200),deltax,deltay
      PARAMETER(g=9.81,pi=3.14159)
C-----------------------------------------------------------------------
C    Open the input and output files
C-----------------------------------------------------------------------
      OPEN(10, FILE='input.dat')
      OPEN(11, FILE='bathy.dat')
      OPEN(12, FILE='COORD.DAT')
C-----------------------------------------------------------------------
C  Read the Input Variables
C-----------------------------------------------------------------------
      READ (10, 103) IM, JM, NM
      READ (10, 104) DT, IREAD
      READ (10, 102) DM
      READ (10, 102) DELTAX, DELTAY
      READ (10, 104) T, SM
      DO N = 1, NM
         READ (10, *) XO(N), YO(N), THETAO(N)
      END DO
      DO I = 1, IM - 1
         DX(I) = DELTAX
      END DO
      DO J = 1, JM - 1
         DY(J) = DELTAY
      END DO
      FORMT = '(25x,2i6)'
C-----------------------------------------------------------------------
C  Either read the bathymetry or generate it
C------------------------------------------------------------------------
      IF (IREAD .EQ. 1) THEN
         P = 1
         Q = 7
         DO WHILE(P .LE. JM)
            Q = MIN0(JM,Q)
            READ (11, FORMT) K
            DO I = 1, IM
               READ (11, 100) (D(I,J), J = P, Q)
  100 format(4x,7(1x,1e9.2))
            END DO
            P = Q + 1
            Q = Q + 7
         END DO
      ELSE
         DO J = 1, JM
            DO I = 1, IM
               D(I,J) = (JM-J)/2.0 - 5.0*SIN(2.0*3.14159*I/33)
               IF (D(I,J) .LE. 0.0) D(I,J) = 0.0
            END DO
         END DO
      ENDIF
  101   format(65x,f10.3)
  102   format(55x,f10.3,3x,f10.3)
  103   format(/60x,I5,2x,I5,2x,I5)
  104   format(55x,f10.3,3x,I5)
C-----------------------------------------------------------------------
C  Calculate the x- and y-directed total distances
C-----------------------------------------------------------------------
      SUMX = 0.0
      DO I = 1, IM - 1
         SUMX = SUMX + DX(I)
      END DO
      SUMY = 0.0
      DO J = 1, JM - 1
         SUMY = SUMY + DY(J)
      END DO
C-----------------------------------------------------------------------
C  Calculate L and C at each node, based on the water depth there and T
C-----------------------------------------------------------------------
      LO = 9.81*T**2/(2.0*3.14159)
      CO = 9.81*T/(2.0*3.14159)
      DO I = 1, IM
         DO J = 1, JM
            L = LO
            TM = L
   10       CONTINUE
            MN = (L+TM)/2.0
            TM = L
            A = 2.0*3.14159*D(I,J)/MN
            L = LO*(EXP(A)-EXP((-A)))/(EXP(A)+EXP((-A)))
            IF (ABS(L-TM) .GT. 0.1) GO TO 10
            C(I,J) = CO*(EXP(A)-EXP((-A)))/(EXP(A)+EXP((-A)))
         END DO
      END DO
C-----------------------------------------------------------------------
C  Calculate the local x- and y-directed velocity gradients
C-----------------------------------------------------------------------
      DO I = 1, IM - 1
         DO J = 1, JM - 1
            CX(I,J) = (C(I+1,J)+C(I+1,J+1)-C(I,J)-C(I,J+1))/(2.0*DX(I))
            CY(I,J) = (C(I,J+1)+C(I+1,J+1)-C(I,J)-C(I+1,J))/(2.0*DY(J))
         END DO
      END DO
C-----------------------------------------------------------------------
C  Start the main loop by setting the initial wave ray locations
C  and angles
C-----------------------------------------------------------------------
      DO 22 N = 1, NM
         THETA = 2.0*3.14159*THETAO(N)/360.0
         X(1) = XO(N)
         Y(1) = YO(N)
         DO ST = 2, SM
            XX = 0.0
            DO I = 1, IM
               XX = XX + DX(I)
               IF (XX .GT. X(ST-1)) GO TO 16
            END DO
   16       CONTINUE
            XX = (X(ST-1)-XX+DX(I))/DX(I)
            YY = 0.0
            DO J = 1, JM
               YY = YY + DY(J)
               IF (YY .GT. Y(ST-1)) GO TO 18
            END DO
   18       CONTINUE
            YY = (Y(ST-1)-YY+DY(J))/DY(J)
C-----------------------------------------------------------------------
C  Interpolate the celerity and celerity gradients at point n
C  from adjacent nodes
C-----------------------------------------------------------------------
            A1 = (XX-1.0)*(YY-1.0)
            A2 = -XX*(YY-1.0)
            A3 = XX*YY
            A4 = -YY*(XX-1.0)
            CC = C(I,J)*A1 + C(I+1,J)*A2 + C(I+1,J+1)*A3 + C(I,J+1)*A4
            IF (CC .LT. 1.0E-05) GO TO 20
            CXX=CX(I,J)*A1+CX(I+1,J)*A2+CX(I+1,J+1)*A3+CX(I,J+1)*A4
            CYY=CY(I,J)*A1+CY(I+1,J)*A2+CY(I+1,J+1)*A3+CY(I,J+1)*A4
C-----------------------------------------------------------------------
C  Now calculate the theta at the next node and its x,y location
C-----------------------------------------------------------------------
            KS = SIN(THETA)/CC*CXX - COS(THETA)/CC*CYY
            DS = DT*CC
            DTHETA = KS*DS
            THETABAR = THETA + DTHETA/2.0
            X(ST) = X(ST-1) + DS*COS(THETABAR)
            Y(ST) = Y(ST-1) + DS*SIN(THETABAR)
            THETA = THETA + DTHETA
            IF (X(ST).LT.0.0 .OR. X(ST).GT.SUMX .OR. Y(ST).LT.0.0 .OR. Y
     #         (ST).GT.SUMY) GO TO 20
         END DO
C-----------------------------------------------------------------------
C  Write the Results to a Regular File
C-----------------------------------------------------------------------
   20    CONTINUE
         ST = ST - 1
         WRITE (12, FORMT) N, ST
         DO I = 1, ST
            XHOLD(I,N) = X(I)
            YHOLD(I,N) = Y(I)
            WRITE (12, 105) XHOLD(I,N), YHOLD(I,N)
  105 format(4x,2e12.5)
         END DO
   22 CONTINUE
      P = 1
      Q = 7
      DO WHILE(P .LE. JM)
         Q = MIN0(JM,Q)
         DO I = 1, IM
            WRITE (11, 106) (D(I,J), J = P, Q)
  106 format(4x,7(1x,1e9.2))
         END DO
         P = Q + 1
         Q = Q + 7
      END DO
      STOP 
      END
