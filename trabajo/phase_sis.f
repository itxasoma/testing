C=================================================================
C  PHASE DIAGRAM  (SIS – Gillespie, red URV)   g77 / F77
C=================================================================
      PROGRAM PHASE
      IMPLICIT NONE
C----- constantes de la red -------------------------------------
      INTEGER         MAXN , MAXE
      PARAMETER      (MAXN=1133 , MAXE=5451)
C----- parámetros epidemiológicos fijos --------------------------
      REAL            DELTA , TMAX
      PARAMETER      (DELTA = 0.1 , TMAX = 200.0)
      INTEGER         SEEDINF
      PARAMETER      (SEEDINF = 57)
C----- barrido y replicas ---------------------------------------
      INTEGER         NVAL , NSIM
      PARAMETER      (NVAL = 30 , NSIM = 100)
      
      

      INTEGER         LINKS(MAXE,2)
      INTEGER         NUML , NUMN
C----- variables de la fase -------------------------------------
      INTEGER         IV , IS , ISEED
      REAL            LAM , SUMRHO , AVGRHO , RATIO
C----- fichero de salida ----------------------------------------
      INTEGER         FOUT , IOS
      CHARACTER*20    FNAME
      FNAME = 'ph_diag_ins.dat'
      FOUT  = 40
C----------------------------------------
C  1) LEER LA RED UNA SOLA VEZ
C----------------------------------------
      CALL READNET (LINKS , NUML , NUMN)
C----------------------------------------
C  2) ABRIR ARCHIVO DE RESULTADOS
C=================================================================
      OPEN (UNIT=FOUT, FILE=FNAME, STATUS='UNKNOWN', IOSTAT=IOS)
      IF (IOS .NE. 0) STOP 'No se pudo abrir phase_diagram.dat'
C----------------------------------------
C  3) BUCLE SOBRE LOS 100 VALORES
C----------------------------------------
      DO 300 IV = 1 , NVAL
C----- equiespaciado en [0 , 0.5]
         LAM = (FLOAT(IV-1)*0.011) / FLOAT(NVAL-1)
         SUMRHO = 0.0
C----- 100 replicas para este
         DO 200 IS = 1 , NSIM
            ISEED = 12345 + IV*1000 + IS
            CALL SRAND (ISEED)
            CALL SIMUL (LAM , DELTA , TMAX , LINKS , NUML , NUMN ,
     &                   SEEDINF , SUMRHO)
 200     CONTINUE
C----- promedio y escritura
         AVGRHO = SUMRHO / FLOAT(NSIM)
         RATIO  = LAM / DELTA
         WRITE (FOUT,'(F9.5,1X,F9.5)') RATIO , AVGRHO
         WRITE (6   ,'(A,F6.3,A,F6.3)') 'lambda=',LAM,
     &      '  <rho>=',AVGRHO
 300  CONTINUE
      CLOSE (FOUT)
      END
C----------------------------------------
C  SUBRUTINA SIMUL: una trayectoria Gillespie, devuelve ? final
C----------------------------------------
      SUBROUTINE SIMUL (LAM, DELT , TMAX , LINKS , NUML , NUMN ,
     &                   SEEDINF , ACCUM)
      IMPLICIT NONE
      INTEGER        MAXN , MAXE
      PARAMETER     (MAXN=1133 , MAXE=5451)
      REAL           LAM , DELT , TMAX , ACCUM
      INTEGER        LINKS(MAXE,2) , NUML , NUMN , SEEDINF
C----- arrays locales -------------------------------------------
      INTEGER        STATE(MAXN) , VINF(MAXN)
      INTEGER        WACT(2,MAXE)
C----- variables locales ----------------------------------------
      INTEGER        NINF , EACT , I , NODE1 , POS , INFN , SUSN
      REAL           RND , TAU , LAMT , T
C----- 1) inicializar estados -----------------------------------
      DO 10 I = 1 , NUMN
         STATE(I) = 0
  10  CONTINUE
      NINF = 0
      DO 20 I = 1 , SEEDINF
         RND   = RAND()
         NODE1 = INT(RND*NUMN) + 1
         IF (STATE(NODE1) .EQ. 0) THEN
            NINF = NINF + 1
            VINF(NINF) = NODE1
            STATE(NODE1) = 1
         END IF
  20  CONTINUE
      CALL UPDEDS (LINKS , NUML , STATE , WACT , EACT)
C----- 2) bucle Gillespie ---------------------------------------
      T = 0.0
 100  CONTINUE
      IF (T .GE. TMAX) GO TO 150
      IF ( (NINF .LE. 0) .AND. (EACT .LE. 0) ) GO TO 150
      LAMT = DELT*FLOAT(NINF) + LAM*FLOAT(EACT)
      IF (LAMT .LE. 0.0) GO TO 150
      RND = RAND()
      TAU = -LOG(RND)/LAMT
      T  = T + TAU
      RND = RAND()
      IF (RND .LT. (LAM*FLOAT(EACT))/LAMT) THEN
C----- infecci¢n -------------------------------------------------
         RND = RAND()
         POS = INT(RND*EACT) + 1
         INFN = WACT(1,POS)
         SUSN = WACT(2,POS)
         IF (STATE(SUSN) .EQ. 0) THEN
            NINF = NINF + 1
            VINF(NINF) = SUSN
            STATE(SUSN) = 1
         END IF
      ELSE
C----- recuperaci¢n ---------------------------------------------
         RND = RAND()
         POS  = INT(RND*NINF) + 1
         INFN = VINF(POS)
         STATE(INFN) = 0
         VINF(POS)  = VINF(NINF)
         NINF       = NINF - 1
      END IF
      CALL UPDEDS (LINKS , NUML , STATE , WACT , EACT)
      GO TO 100
C----- 3) devolver final --------------------------------------
 150  CONTINUE
      ACCUM = ACCUM + FLOAT(NINF)/FLOAT(NUMN)
      RETURN
      END
C----------------------------------------
C  SUBRUTINA: LECTURA DE LA RED  (URV.txt)
C----------------------------------------
      SUBROUTINE READNET (LINKS , NUML , NUMN)
      IMPLICIT NONE
      INTEGER        MAXN , MAXE
      PARAMETER     (MAXN=1133 , MAXE=5451)
      INTEGER        LINKS(MAXE,2) , NUML , NUMN
      INTEGER        U , V , IU , IOS , I
      LOGICAL        PRES(MAXN)
      DO 5 I = 1 , MAXN
         PRES(I) = .FALSE.
  5   CONTINUE
      NUML = 0
      IU   = 10
      OPEN (UNIT=IU, FILE='URV.txt', STATUS='OLD', IOSTAT=IOS)
      IF (IOS .NE. 0) STOP 'No se puede abrir URV.txt'
 20   READ (IU,*,IOSTAT=IOS) U , V
      IF (IOS .NE. 0) GO TO 30
      NUML = NUML + 1
      LINKS(NUML,1) = U
      LINKS(NUML,2) = V
      PRES(U) = .TRUE.
      PRES(V) = .TRUE.
      GO TO 20
 30   CONTINUE
      CLOSE (IU)
      NUMN = 0
      DO 40 I = 1 , MAXN
         IF (PRES(I)) NUMN = I
 40   CONTINUE
      RETURN
      END
C----------------------------------------
C  SUBRUTINA: ACTUALIZA ARISTAS ACTIVAS S-I
C----------------------------------------
      SUBROUTINE UPDEDS (LINKS , NUML , STATE , WACT , EACT)
      IMPLICIT NONE
      INTEGER        MAXN , MAXE
      PARAMETER     (MAXN=1133 , MAXE=5451)
      INTEGER        LINKS(MAXE,2) , STATE(MAXN)
      INTEGER        WACT(2,MAXE) , NUML , EACT
      INTEGER        I , U , V
      EACT = 0
      DO 60 I = 1 , NUML
         U = LINKS(I,1)
         V = LINKS(I,2)
         IF (STATE(U).EQ.1 .AND. STATE(V).EQ.0) THEN
            EACT = EACT + 1
            WACT(1,EACT) = U
            WACT(2,EACT) = V
         ELSE IF (STATE(V).EQ.1 .AND. STATE(U).EQ.0) THEN
            EACT = EACT + 1
            WACT(1,EACT) = V
            WACT(2,EACT) = U
         END IF
 60   CONTINUE
      RETURN
      END

