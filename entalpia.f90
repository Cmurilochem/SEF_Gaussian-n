 PROGRAM entalpia
 IMPLICIT NONE

 INTEGER :: atomtypetot
 INTEGER :: err, i
 CHARACTER (LEN = 2), ALLOCATABLE, DIMENSION( : ) :: ATOMS  ! atomos da molecula
 INTEGER, ALLOCATABLE, DIMENSION( : ) :: NATOMS  ! o numero de atomos de cada tipo na molecula
 DOUBLE PRECISION, ALLOCATABLE, DIMENSION( : ) :: ATOMENERGIES
 DOUBLE PRECISION :: MOLENERGY ! energia da molecula
 DOUBLE PRECISION :: ATOMIZATION, sum  ! atomizacao
 DOUBLE PRECISION :: SEF_0K, sum2  ! entalpia padrao de formacao a 0K
 DOUBLE PRECISION :: SEF_298K      ! entalpia padrao de formacao a 298,15K
 DOUBLE PRECISION, PARAMETER :: kcal = 627.5095D+00 ! conversao de hartrees para kcal mol-1
 REAL, ALLOCATABLE, DIMENSION( : ) :: ATOMSENTHALPY  ! vetor de entalpia de formacao experimentais do atomos
 REAL, ALLOCATABLE, DIMENSION( : ) :: ATOMSENTHALPYTHERMALCORR ! vetor da correcao a 298.15 das entalpia de formacao experimentais do atomos
 INTEGER :: freqsize   ! numero de frequencies vibracionais harmonicas
 DOUBLE PRECISION, PARAMETER :: c = 2.99792458D+10 ! Velocidade de Luz em cm.s-1   
 DOUBLE PRECISION, PARAMETER :: h =  6.62606896D-034 ! Constante de Planck em J.s
 DOUBLE PRECISION, PARAMETER :: kb = 1.3806504D-023 ! Constante de Boltzmann em J.K-1
 DOUBLE PRECISION, PARAMETER :: jh = 4.3597439D-18 ! 1 Hartree e igual a 4.35974417E-18 Joules 
 DOUBLE PRECISION :: RT   ! TERMO RT 
 DOUBLE PRECISION :: S, FAC
 DOUBLE PRECISION :: EZPE     ! energia do ponto zero total
 DOUBLE PRECISION :: Etrans   ! Correcao termica na energia devido ao movimento translacional
 DOUBLE PRECISION :: Erotlin  ! Correcao termica na energia devido a rotacao para moleculas lineares
 DOUBLE PRECISION :: Erotang  ! Correcao termica na energia devido a rotacao para moleculas nao-lineares
 DOUBLE PRECISION :: Evib     ! Correcao termica na energia devido a vibracao a 298k
 DOUBLE PRECISION :: Evibtotal ! energia do ponto zero total + Correcao termica na energia devido a vibracao a 298k
 DOUBLE PRECISION :: Etermtotal ! energia termica a 298 incluindo a energia do ponto zero 
 DOUBLE PRECISION :: Eterm ! energia termica a 298 excluindo a energia do ponto zero 
 DOUBLE PRECISION :: Htermtotal ! Correcao na Entalpia incluindo a energia do ponto zero
 DOUBLE PRECISION :: Hterm ! Correcao na Entalpia excluindo a energia do ponto zero
 DOUBLE PRECISION, PARAMETER :: T = 298.15D+00        ! Temperatura a 298.15K
 DOUBLE PRECISION :: scale ! fator de escalonamento das frequencias
 DOUBLE PRECISION, ALLOCATABLE, DIMENSION( : ) :: FREQS   ! vetor de frequencies vib. harmonias em cm-1
 DOUBLE PRECISION, ALLOCATABLE, DIMENSION( : ) :: FREQSSEG ! vetor de frequencies vib. harmonias em s-1
 DOUBLE PRECISION, ALLOCATABLE, DIMENSION( : ) :: VIBTERM  ! fracao da correcao vibracional depente da temperatura
 CHARACTER(LEN=6) :: molecula  ! molecula linear ou nao 
 DOUBLE PRECISION :: sum3, sum4, sum5 
 DOUBLE PRECISION :: ATOMCORR ! somatoria das correcoes a 298.15 das entalpia de formacao experimentais do atomos


 !##################################################################################################################
 !RETIRANDO OS DADOS DO ARQUIVO mol.entalpia
 !##################################################################################################################

 OPEN(UNIT=11,FILE='mol.entalpia',IOSTAT=err,FORM='FORMATTED',STATUS='OLD',ACTION='READ')
 IF (err /= 0) THEN
 PRINT*, "ERRO NA ABERTURA DO ARQUIVO"
 END IF

 READ(UNIT=11,FMT='(I4)') atomtypetot ! tipos de atomos que a molecual contem

 ALLOCATE (ATOMS(atomtypetot), STAT=err)
 IF (err /= 0) PRINT*, "FALHOU A ALOCACAO DO VETOR ATOMS"
 
 READ(11,*) ATOMS  ! vetor com os atomos presentes na molecula

 ALLOCATE (NATOMS(atomtypetot), STAT=err)
 IF (err /= 0) PRINT*, "FALHOU A ALOCACAO DO VETOR NATOMS"

 READ(11,*) NATOMS ! vetor com o numero de atomos de cada tipo presente na molecula

 ALLOCATE (ATOMENERGIES(atomtypetot), STAT=err)
 IF (err /= 0) PRINT*, "FALHOU A ALOCACAO DO VETOR ATOMENERGIES"

 DO i=1,atomtypetot
 READ(UNIT=11, FMT='(F30.10)') ATOMENERGIES(i)        ! energia dos atomos
 END DO 

 READ(UNIT=11, FMT='(F30.10)') MOLENERGY              ! energia das moleculas 

 READ(UNIT=11,FMT='(I4)') freqsize     ! Numero de modos normais de vibracao 

 ALLOCATE (FREQS(freqsize), STAT=err)
 IF (err /= 0) PRINT*, "FALHOU A ALOCACAO DO VETOR FREQS" 

 DO i=1,freqsize
 READ(UNIT=11,FMT='(F10.5)') FREQS(i)     ! Vetor de Frequencies de Vibrações Harmonicas em (cm-1), ou seja, numeros de onda v
 END DO

 READ(UNIT=11,FMT='(F10.5)') scale         ! parametro de escalonamento das frequencias

 READ(UNIT=11,FMT='(A)') molecula          ! molecula linear ou nao

 CLOSE(UNIT=11)

 !##################################################################################################################
 !FIM DOS DADOS DO ARQUIVO mol.entalpia
 !##################################################################################################################

 !##################################################################################################################
 !CALCULANDO AS ENERGIAS DE ATOMIZAÇÂO
 !##################################################################################################################
 sum = 0

 DO i=1,atomtypetot                                         ! calculando a energia de atomizacao
   sum = sum + ATOMENERGIES(i)*NATOMS(i)
 END DO 

 ATOMIZATION = (sum - MOLENERGY)*kcal

 PRINT*, "###############################################################################################"
 PRINT*, "A ENERGIA DE ATOMIZAÇÃO EM kcal mol-1 É ...", ATOMIZATION 
 PRINT*, "###############################################################################################" 

 !##################################################################################################################
 !FIM DAS ENERGIAS DE ATOMIZAÇÂO
 !##################################################################################################################

 !##################################################################################################################
 !CONSTRUINDO O VETOR DE ENTALPIAS PADRAO DE FORMACAO EXPERIMENTAL DOS ATOMOS
 !##################################################################################################################

 ALLOCATE (ATOMSENTHALPY(atomtypetot), STAT=err)
 IF (err /= 0) PRINT*, "FALHOU A ALOCACAO DO VETOR ATOMSENTHALPY"

 DO i=1,atomtypetot
  IF (ATOMS(i) == "H") THEN
    ATOMSENTHALPY(i) = 51.63D+00
  ELSE IF (ATOMS(i) == "Li") THEN
    ATOMSENTHALPY(i) = 37.69D+00 
  ELSE IF (ATOMS(i) == "Be") THEN
    ATOMSENTHALPY(i) = 76.48D+00 
  ELSE IF (ATOMS(i) == "B") THEN
    ATOMSENTHALPY(i) = 136.2D+00
  ELSE IF (ATOMS(i) == "C") THEN
    ATOMSENTHALPY(i) = 169.98D+00
  ELSE IF (ATOMS(i) == "N") THEN 
    ATOMSENTHALPY(i) = 112.53D+00
  ELSE IF (ATOMS(i) == "O") THEN 
     ATOMSENTHALPY(i) = 58.99D+00
  ELSE IF (ATOMS(i) == "F") THEN
      ATOMSENTHALPY(i) = 18.47D+00  
  ELSE IF (ATOMS(i) == "Na") THEN
     ATOMSENTHALPY(i) = 25.69D+00 
  ELSE IF (ATOMS(i) == "Mg") THEN 
     ATOMSENTHALPY(i) = 34.87D+00
  ELSE IF (ATOMS(i) == "Al") THEN
     ATOMSENTHALPY(i) = 78.23D+00 
  ELSE IF (ATOMS(i) == "Si") THEN
      ATOMSENTHALPY(i) = 106.6D+00 
  ELSE IF (ATOMS(i) == "P") THEN
       ATOMSENTHALPY(i) = 75.42D+00
  ELSE IF (ATOMS(i) == "S") THEN 
       ATOMSENTHALPY(i) = 65.66D+00 
  ELSE IF (ATOMS(i) == "Cl") THEN 
      ATOMSENTHALPY(i) = 28.59D+00 
  ELSE IF (ATOMS(i) == "Br") THEN
      ATOMSENTHALPY(i) = 28.18D+00
  ENDIF 
 END DO

 !##################################################################################################################
 !FIM DO VETOR DE ENTALPIAS PADRAO DE FORMACAO EXPERIMENTAL DOS ATOMOS
 !##################################################################################################################

 !##################################################################################################################
 !CALCULANDO A ENTALPIA PADRAO DE FORMACAO DO MOLECULA A 0K
 !##################################################################################################################

 sum2 = 0

 DO i=1,atomtypetot                                         ! calculando a energia de atomizacao
   sum2 = sum2 + ATOMSENTHALPY(i)*NATOMS(i)
 END DO

 SEF_0K = sum2 - ATOMIZATION

 PRINT*, "###############################################################################################"
 PRINT*, "A ENTALPIA PADRÃO DE FORMAÇÃO A 0K E EM kcal mol-1 É ...", SEF_0K
 PRINT*, "###############################################################################################"

 !##################################################################################################################
 ! FIM DO CALCULO DA ENTALPIA PADRAO DE FORMACAO DO MOLECULA A 0K
 !##################################################################################################################
 
 !##################################################################################################################
 !CONSTRUINDO O VETOR DAS CORRECOES DAS ENTALPIAS PADRAO DE FORMACAO EXPERIMENTAL DOS ATOMOS a 298K
 !##################################################################################################################

 ALLOCATE (ATOMSENTHALPYTHERMALCORR(atomtypetot), STAT=err)
 IF (err /= 0) PRINT*, "FALHOU A ALOCACAO DO VETOR ATOMSENTHALPYTHERMALCORR"

 DO i=1,atomtypetot
  IF (ATOMS(i) == "H") THEN
    ATOMSENTHALPYTHERMALCORR(i) = 1.01D+00
  ELSE IF (ATOMS(i) == "Li") THEN
    ATOMSENTHALPYTHERMALCORR(i) = 1.10D+00
  ELSE IF (ATOMS(i) == "Be") THEN
    ATOMSENTHALPYTHERMALCORR(i) = 0.46D+00
  ELSE IF (ATOMS(i) == "B") THEN
    ATOMSENTHALPYTHERMALCORR(i) = 0.29D+00
  ELSE IF (ATOMS(i) == "C") THEN
    ATOMSENTHALPYTHERMALCORR(i) = 0.25D+00
  ELSE IF (ATOMS(i) == "N") THEN
    ATOMSENTHALPYTHERMALCORR(i) = 1.04D+00
  ELSE IF (ATOMS(i) == "O") THEN
     ATOMSENTHALPYTHERMALCORR(i) = 1.04D+00
  ELSE IF (ATOMS(i) == "F") THEN
      ATOMSENTHALPYTHERMALCORR(i) = 1.05D+00
  ELSE IF (ATOMS(i) == "Na") THEN
     ATOMSENTHALPYTHERMALCORR(i) = 1.54D+00
  ELSE IF (ATOMS(i) == "Mg") THEN
     ATOMSENTHALPYTHERMALCORR(i) = 1.19D+00
  ELSE IF (ATOMS(i) == "Al") THEN
     ATOMSENTHALPYTHERMALCORR(i) = 1.08D+00
  ELSE IF (ATOMS(i) == "Si") THEN
      ATOMSENTHALPYTHERMALCORR(i) = 0.76D+00
  ELSE IF (ATOMS(i) == "P") THEN
       ATOMSENTHALPYTHERMALCORR(i) = 1.28D+00
  ELSE IF (ATOMS(i) == "S") THEN
       ATOMSENTHALPYTHERMALCORR(i) = 1.05D+00
  ELSE IF (ATOMS(i) == "Cl") THEN
      ATOMSENTHALPYTHERMALCORR(i) = 1.10D+00
  ELSE IF (ATOMS(i) == "Br") THEN
      ATOMSENTHALPYTHERMALCORR(i) = 3.00D+00
  ENDIF
 END DO

 !##################################################################################################################
 !FIM DO VETOR DAS CORRECOES DAS ENTALPIAS PADRAO DE FORMACAO EXPERIMENTAL DOS ATOMOS a 298K
 !################################################################################################################## 

 
 !##################################################################################################################
 !CALCULANDO A CORRECAO TERMICA (SEM O ZPE) NA ENTALPIA A 298.15K
 !##################################################################################################################

 ALLOCATE (FREQSSEG(freqsize), STAT=err)
 IF (err /= 0) PRINT*, "FALHOU A ALOCACAO DO VETOR FREQSSEG" 

 DO i=1,freqsize
 FREQSSEG(i) = (0.5D+00*h*(scale*c*FREQS(i)))/jh     ! Vetor de Energias do Ponto Zero (em hartrees) para cada modo normal de vibracao
 END DO

 sum3 = 0 
 
 DO i=1,freqsize
  sum3 = sum3 + FREQSSEG(i)
 END DO
 
 EZPE = sum3   ! energia do ponto zero total
  
 Etrans = ((1.5D+00*kb*T)/jh)  ! Correcao termica na energia devido a tranlacao (em hartrees)

 Erotlin = ((kb*T)/jh)      ! Correcao termica na energia devido a rotacao linear (em hartrees)

 Erotang = ((1.5D+00*kb*T)/jh)  ! Correcao termica na energia devido a rotacao nao linear (em hartrees)

 ALLOCATE (VIBTERM(freqsize), STAT=err)
 IF (err /= 0) PRINT*, "FALHOU A ALOCACAO DO VETOR VIBTERM"

 S = h*scale*c
 
 RT = (kb*T)

 DO i=1,freqsize
 VIBTERM(i) = ((S*FREQS(i))/(exp((S*FREQS(i))/RT) - 1))/jh   ! Vetor de correcoes da energia interna devido a vibracao
 END DO

 sum4 = 0

 DO i=1,freqsize
  sum4 = sum4 + VIBTERM(i)
 END DO

 Evib = sum4

 IF (molecula == 'LIN') THEN    ! para moleculas linearess.. devemos considerar RT para rotacao
  
 Etermtotal = Etrans + Erotlin + Evibtotal  ! correcao termica na energia interna devido rotacao, vibracao(incluindo ZPE)  e translacao
 Eterm = Etrans + Erotlin + Evib  ! correcao termica na energia interna devido rotacao, vibracao(excluindo ZPE)  e translacao

 ELSEIF (molecula == 'NOLIN') THEN  ! para moleculas nao lineares ... devemos considerar 3/2RT para rotacao

 Etermtotal = Etrans + Erotang + Evibtotal 

 Eterm = Etrans + Erotang + Evib 

 END IF 
 
 Htermtotal = Etermtotal + RT/jh   ! Correcao na Entalpia incluindo a energia do ponto zero

 Hterm = Eterm + RT/jh   ! Correcao na Entalpia excluindo a energia do ponto zero

 !##################################################################################################################
 !FIM DO CALCULO DA CORRECAO TERMICA (SEM O ZPE) NA ENTALPIA A 298.15K
 !##################################################################################################################

 !##################################################################################################################
 !CALCULANDO A ENTALPIA PADRAO DE FORMACAO DO MOLECULA A 298.15 K
 !##################################################################################################################

 sum5 = 0

 DO i=1,atomtypetot                                         ! calculando a correção a 298K para os atomos
   sum5 = sum5 + ATOMSENTHALPYTHERMALCORR(i)*NATOMS(i)
 END DO 

 ATOMCORR = sum5

 SEF_298K = SEF_0K + kcal*(Hterm) - ATOMCORR

 PRINT*, "###############################################################################################"
 PRINT*, "A ENTALPIA PADRÃO DE FORMAÇÃO A 298,15K E EM kcal mol-1 É ...", SEF_298K
 PRINT*, "###############################################################################################"  

 !##################################################################################################################
 !FIM DO CALCULO DA ENTALPIA PADRAO DE FORMACAO DO MOLECULA A 298.15 K
 !##################################################################################################################

 
 IF (ALLOCATED(ATOMS)) DEALLOCATE (ATOMS, STAT=err)
 IF (ALLOCATED(NATOMS)) DEALLOCATE (NATOMS, STAT=err)
 IF (ALLOCATED(ATOMENERGIES)) DEALLOCATE (ATOMENERGIES, STAT=err)
 IF (ALLOCATED(FREQS)) DEALLOCATE (FREQS, STAT=err)
 IF (ALLOCATED(ATOMSENTHALPY)) DEALLOCATE (ATOMSENTHALPY, STAT=err)
 IF (ALLOCATED(ATOMSENTHALPYTHERMALCORR)) DEALLOCATE (ATOMSENTHALPYTHERMALCORR, STAT=err)
 IF (ALLOCATED(FREQSSEG)) DEALLOCATE (FREQSSEG, STAT=err)
 IF (ALLOCATED(VIBTERM)) DEALLOCATE (VIBTERM, STAT=err)

 END PROGRAM entalpia

