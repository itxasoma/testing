C===============================================================
C  SIS  Gillespie
C===============================================================
      PROGRAM SIS
      IMPLICIT NONE
C----- PARAMETROS ------------------------------------------------
      INTEGER          MAXN , MAXE
      PARAMETER       (MAXN = 1133,
     &                 MAXE = 5451)

      REAL             LAMBDA , DELTA , TMAX
      PARAMETER       (LAMBDA = 0.1 ,
     &                 DELTA  = 0.1 ,
     &                 TMAX   = 200.0)

      INTEGER          SEEDINF
      PARAMETER       (SEEDINF = 57)
C----- ARRAYS ----------------------------------------------------
      INTEGER          LINKS (MAXE ,2)
      INTEGER          STATE (MAXN)
      INTEGER          VINF  (MAXN)
      INTEGER          WACT  (2 ,MAXE)
C----- VARIABLES -------------------------------------------------
      INTEGER          NUMN , NUML , NINF , EACT
      INTEGER          I , NODE1 , NODE2 , INFN , SUSN , POS
      REAL             RND , TAU , LAMBDA_TOT , T, FRAC
      INTEGER          OUTU , IOS
C-------------------------------------------------
C  1) LEER RED
C-------------------------------------------------
      CALL READNET (LINKS , NUML , NUMN)
C-------------------------------------------------
C  2) INICIALIZAR RNG Y ESTADOS
C-------------------------------------------------
      CALL SRAND (123456)
      DO  10 I = 1 , MAXN
         STATE(I) = 0
  10  CONTINUE
      NINF = 0
      DO  20 I = 1 , SEEDINF
         RND   = RAND()
         NODE1 = INT(RND*NUMN) + 1
         IF (STATE(NODE1) .EQ. 0) THEN
            NINF         = NINF + 1
            VINF(NINF)   = NODE1
            STATE(NODE1) = 1
         END IF
  20  CONTINUE
      CALL UPDEDS (LINKS , NUML , STATE , WACT , EACT)
C-------------------------------------------------
C  3) ABRIR SALIDA
C-------------------------------------------------
      OUTU = 20
      OPEN (UNIT=OUTU,
     & FILE='C:\Users\miria\Documents\Python Scripts\Redes complejas\'//
     &      '5.dat',
     & STATUS='UNKNOWN', IOSTAT=IOS)
      IF (IOS .NE. 0) STOP 'No se ha podido abrir el archivo'

C-------------------------------------------------
C  4) BUCLE DE GILLESPIE----------------------------
      T = 0.0
 100  CONTINUE
      IF (T .GE. TMAX) GO TO 200
      IF ( (NINF .LE. 0) .AND. (EACT .LE. 0) ) GO TO 200

      LAMBDA_TOT = DELTA*NINF + LAMBDA*EACT
      IF (LAMBDA_TOT .LE. 0.0) GO TO 200

      RND = RAND()
      TAU = -LOG(RND)/LAMBDA_TOT
      T   = T + TAU

      RND = RAND()
      IF (RND .LT. (LAMBDA*EACT)/LAMBDA_TOT) THEN
C----- evento de infeccion -------------------------------------
         RND = RAND()
         POS = INT(RND*EACT) + 1
         INFN = WACT(1,POS)
         SUSN = WACT(2,POS)
         IF (STATE(SUSN) .EQ. 0) THEN
            NINF         = NINF + 1
            VINF(NINF)   = SUSN
            STATE(SUSN)  = 1
         END IF
      ELSE
C----- evento de recuperacion ----------------------------------
         RND = RAND()
         POS = INT(RND*NINF) + 1
         INFN = VINF(POS)
         STATE(INFN) = 0
         VINF(POS)  = VINF(NINF)
         NINF       = NINF - 1
      END IF

      CALL UPDEDS (LINKS , NUML , STATE , WACT , EACT)
      FRAC = FLOAT(NINF)/FLOAT(NUMN)    ! RO(t)
      WRITE (6,*) 't=',T,'  rho=',FRAC
      WRITE (OUTU,*) T , FRAC
      GO TO 100

 200  CONTINUE
      CLOSE (OUTU)
      END
C-------------------------------------------------
C  SUBRUTINA: LECTURA DE LA RED
C-------------------------------------------------
      SUBROUTINE READNET (LINKS , NUML , NUMN)
      IMPLICIT NONE
      INTEGER          MAXN , MAXE
      PARAMETER       (MAXN = 1133,
     &                 MAXE = 5451)
      INTEGER          LINKS(MAXE,2)
      INTEGER          NUML , NUMN
      INTEGER          U , V , IU , IOS , I
      LOGICAL          PRES(MAXN)

      DO  5 I = 1 , MAXN
         PRES(I) = .FALSE.
   5  CONTINUE
      NUML  = 0
      IU    = 10
      OPEN (UNIT=IU, FILE='URV.txt', STATUS='OLD', IOSTAT=IOS)
      IF (IOS .NE. 0) STOP 'No se puede abrir URV.txt'

   20 READ (IU,*,IOSTAT=IOS) U , V
      IF (IOS .NE. 0) GO TO 30
      NUML           = NUML + 1
      LINKS(NUML,1)  = U
      LINKS(NUML,2)  = V
      PRES(U) = .TRUE.
      PRES(V) = .TRUE.
      GO TO 20
   30 CONTINUE
      CLOSE (IU)
C----- calcula numero de nodos ---------------------------------
      NUMN = 0
      DO  40 I = 1 , MAXN
         IF (PRES(I)) NUMN = I
  40  CONTINUE
      RETURN
      END
C-------------------------------------------------
C  SUBRUTINA: ACTUALIZA ARISTAS ACTIVAS S-I
C-------------------------------------------------
      SUBROUTINE UPDEDS (LINKS , NUML , STATE , WACT , EACT)
      IMPLICIT NONE
      INTEGER          MAXN , MAXE
      PARAMETER       (MAXN = 1133,
     &                 MAXE = 5451)
      INTEGER          LINKS(MAXE,2)
      INTEGER          STATE(MAXN)
      INTEGER          WACT(2,MAXE)
      INTEGER          NUML , EACT
      INTEGER          I , U , V
      EACT = 0
      DO  60 I = 1 , NUML
         U = LINKS(I,1)
         V = LINKS(I,2)
         IF (STATE(U) .EQ. 1 .AND. STATE(V) .EQ. 0) THEN
            EACT = EACT + 1
            WACT(1,EACT) = U
            WACT(2,EACT) = V
         ELSE IF (STATE(V) .EQ. 1 .AND. STATE(U) .EQ. 0) THEN
            EACT = EACT + 1
            WACT(1,EACT) = V
            WACT(2,EACT) = U
         END IF
   60 CONTINUE
      RETURN
      END

