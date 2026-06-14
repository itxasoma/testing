! Network dynamics: SIS dynamics
! ARMANDO AGUILAR CAMPOS 
! JUN 2025
PROGRAM SIS_dyn
IMPLICIT NONE
! VALORS
INTEGER :: L,ndim,iseed,i,ios,tmp1,tmp2,tmp3,j,k,n,step,num_chosen,idx_lambda
INTEGER :: N_links,N_nodes,NI,Eact
INTEGER :: index1, old_infected, infected_node, old_pos_new_active, new_infected_node
INTEGER :: neighbor,old_pos, recovered_node
PARAMETER (ndim=2)
INTEGER,DIMENSION(:),ALLOCATABLE :: infected,Pin,Pout
INTEGER,DIMENSION(:,:),ALLOCATABLE :: links,Wact
INTEGER,DIMENSION(:),ALLOCATABLE :: neigh,neigh_pos
INTEGER,DIMENSION(2) :: old_link_value
LOGICAL,DIMENSION(:),ALLOCATABLE :: is_infected
REAL :: r1279,RN1,RN2,RN3
REAL :: t,delta,lambda,lambda_tot,NI_perc
REAL, DIMENSION(:), ALLOCATABLE :: lambda_values
LOGICAL :: already_infected,neighbor_is_infected
CHARACTER(LEN=100) :: data_file,file_name
CHARACTER(LEN=20)  :: lambda_str
COMMON/constants/iseed

iseed = 2025
CALL setr1279(iseed) !initialize RNG seed

data_file = 'web-polblogs.mtx'!'links_exemple.txt'
N_links = 0
! Open the file for number rows
OPEN(25, file=data_file, status='old', action='read', iostat=ios)
IF (ios /= 0) THEN
	WRITE(*,*)'Error opening file.'
	STOP
END IF
! Count lines by reading values
DO
	READ(25, *, iostat=ios) tmp1, tmp2
	IF (ios /= 0) EXIT
	N_links = N_links + 1
END DO
ALLOCATE(links(N_links,2),Wact(N_links,2))
CLOSE(25)

! Open the file for reading edge list
OPEN(25, file=data_file, status='old', action='read', iostat=ios)
! Read values into array
WRITE(*,*)'Number links:', N_links
DO i = 1, N_links
	READ(25, *, iostat=ios) links(i,1), links(i,2)
	Wact(i,:) = links(i,:)
	! WRITE(*,*)links(i,1),links(i,2)
	IF (ios /= 0) EXIT
END DO
CLOSE(25)

N_nodes = MAXVAL(links(:,1))
ALLOCATE(infected(N_nodes),is_infected(N_nodes))
ALLOCATE(neigh(N_links),neigh_pos(N_links))
ALLOCATE(Pin(N_nodes),Pout(N_nodes))

ALLOCATE(lambda_values(12))
lambda_values = [0.001,0.005,0.01,0.03,0.06,0.10,0.50,1.00,2.00,3.00,4.00,5.00]
! ALLOCATE(lambda_values(3))
! lambda_values = [0.001,0.005,0.01]

DO idx_lambda = 1,SIZE(lambda_values)
	! Open the file for reading edge list
	OPEN(25, file=data_file, status='old', action='read', iostat=ios)
	! Read values into array
	WRITE(*,*)'Number links:', N_links
	DO i = 1, N_links
		READ(25, *, iostat=ios) links(i,1), links(i,2)
		Wact(i,:) = links(i,:)
		! WRITE(*,*)links(i,1),links(i,2)
		IF (ios /= 0) EXIT
	END DO
	CLOSE(25)

	! list of all nodes
	is_infected(:) = .FALSE.
	! Print what we read
	WRITE(*,*)'Number nodes:', N_nodes
	DO i = 1, N_nodes
	    infected(i) = i
	    ! WRITE(*,*)infected(i)
	END DO

	! create list of neighbours
	k = 1
	DO i = 1,N_nodes
		Pin(i) = k
		DO j = 1,N_links
			IF (links(j, 1) == i) THEN
				neigh(k) = links(j,2)
				neigh_pos(k) = j !position in Wact 
				k = k+1
			END IF
		Pout(i) = k-1
		END DO
	END DO
	! WRITE(*,*)neigh,neigh_pos
	! WRITE(*,*)Pin,Pout

	! Initialize Infected
	NI = MAX(INT(N_nodes * 0.2), 1)
	num_chosen = 0
	DO WHILE (num_chosen < NI)
	    RN1 = r1279()
	    index1 = INT(RN1 * N_nodes) + 1
	    IF (index1 > N_nodes) index1 = N_nodes

	    IF (.NOT. is_infected(index1)) THEN
	        num_chosen = num_chosen + 1
	        infected(num_chosen) = index1
	        is_infected(index1) = .TRUE.
	    END IF
	END DO

	! Initialize active links
	Wact = links
	Eact = 0
	DO i = 1, NI
	    infected_node = infected(i)

	    DO j = Pin(infected_node), Pout(infected_node)
	        neighbor = neigh(j)

	        ! Check if neighbor is already infected
	        already_infected = .FALSE.
	        DO k = 1, NI
	            IF (infected(k) == neighbor) THEN
	                already_infected = .TRUE.
	                EXIT
	            END IF
	        END DO

	        old_pos = neigh_pos(j)
	        IF (already_infected) THEN
	            ! REMOVE link if it is active
	            IF (old_pos <= Eact) THEN
	                ! Swap with the last active link (Eact)
	                DO n = 1, N_links
	                    IF (neigh_pos(n) == Eact) THEN
	                        neigh_pos(n) = old_pos
	                        EXIT
	                    END IF
	                END DO
	                neigh_pos(j) = Eact

	                ! Swap active link rows
	                old_link_value(:) = Wact(Eact,:)
	                Wact(Eact,:) = Wact(old_pos,:)
	                Wact(old_pos,:) = old_link_value(:)

	                Eact = Eact - 1
	            END IF

	        ELSE
	            ! ADD link if not already active
	            IF (old_pos > Eact) THEN
	                Eact = Eact + 1

	                neigh_pos(j) = Eact
	                DO n = 1, N_links
	                    IF (neigh_pos(n) == Eact .AND. n /= j) THEN
	                        neigh_pos(n) = old_pos
	                        EXIT
	                    END IF
	                END DO

	                ! Swap active link rows
	                old_link_value(:) = Wact(Eact,:)
	                Wact(Eact,:) = Wact(old_pos,:)
	                Wact(old_pos,:) = old_link_value(:)
	            END IF
	        END IF
	        IF (Eact < 1 .OR. Eact > N_links) THEN
			    WRITE(*,*)"Eact out of bounds!", Eact, N_links
			    EXIT
			END IF
	    END DO
	END DO 	
	! DO i = 1, N_links
	!     WRITE(*,*)Wact(i,:)
	! END DO

	! main program
	delta = 1.0
	t = 0.0
	lambda = lambda_values(idx_lambda)
	WRITE(*,*) 'Lambda:', lambda
	WRITE(lambda_str, '(F5.3)') lambda
	file_name = "SISdata_lambda"//TRIM(lambda_str)//".dat"
	OPEN(11,file=file_name,STATUS='UNKNOWN')
	NI_perc = float(NI)/N_nodes
	WRITE(11,'(F10.5, F10.5)')t,NI_perc
	step = 0
	DO WHILE (t < 14.0 .AND. NI > 0)
		step = step+1
		IF (mod(step,10**3).eq.0) THEN
			WRITE(*,*)step
		END IF
		lambda_tot = delta * NI + lambda * Eact
		IF (ABS(lambda_tot) < 1.0E-10) EXIT
		RN1 = r1279() 
		! IF (step==394) THEN 
		! 	WRITE(*,*)RN1,(lambda * Eact / lambda_tot),lambda*Eact,lambda_tot,"NI",NI,"Eact",Eact
		! END IF
		IF (RN1 < (lambda * Eact / lambda_tot)) THEN
			! IF (step==394) THEN  
			! 	WRITE(*,*)"INFECTED"!,infected
			! END IF
			! Infection event
		  	! select infected node
		  	RN2 = r1279()
			index1 = INT(RN2 * Eact) + 1
			IF (index1 > Eact) index1 = Eact
			new_infected_node = Wact(index1,2)
			! change positions to add new infected node

			! maybe look for position before infection and swap nodes

			IF (.NOT. is_infected(new_infected_node)) THEN
			    NI = NI + 1
			    infected(NI) = new_infected_node
			    is_infected(new_infected_node) = .TRUE.
			END IF
			! IF (step==394) THEN 
			! 	WRITE(*,*)"NI",NI,"Eact",Eact, "index1", index1, new_infected_node
			! END IF
			! WRITE(*,*)infected
		ELSE
			! IF (step==394) THEN  
			! 	WRITE(*,*)"RECOVERY"!,infected
			! END IF
			! Recovery event
			RN2 = r1279()
			index1 = INT(RN2 * NI) + 1
			IF (index1 > NI) index1 = NI
			! change positions to delete infected node
			old_infected = infected(NI)
			infected(NI) = infected(index1)
			infected(index1) = old_infected
			NI = NI-1

			! delete from active links all links related to recovery node
			recovered_node = infected(NI+1)
			is_infected(recovered_node) = .FALSE.
			IF (Eact > 0) THEN
				DO j = Eact, 1, -1
					old_pos = neigh_pos(j)
			        IF (Wact(j,1) == recovered_node) THEN
			            ! Swap with last active
			            DO n = 1, N_links
			                IF (neigh_pos(n) == Eact) THEN
			                    neigh_pos(n) = old_pos
			                    EXIT
			                END IF
			            END DO
			            neigh_pos(j) = Eact
			            ! Swap active link positions
			            old_link_value(:) = Wact(Eact,:)
			            Wact(Eact,:) = Wact(old_pos,:)
			            Wact(old_pos,:) = old_link_value(:)
			            Eact = Eact - 1
			        END IF
				END DO
			END IF

			! IF (step==394) THEN  
			! 	WRITE(*,*)"NI",NI,"Eact",Eact, "index1", index1,recovered_node
			! END IF
			! WRITE(*,*)infected
		END IF

		DO i = 1, NI 
		    infected_node = infected(i)
		    ! IF (step==394) THEN  
		    ! 	WRITE(*,*)infected_node,"PINS", Pin(infected_node), Pout(infected_node)
		    ! END IF	

		    DO j = Pin(infected_node), Pout(infected_node)
		        neighbor = neigh(j)
		        ! WRITE(*,*)"NEIGH", j, neighbor

		        ! Check if neighbor is already infected
		        already_infected = .FALSE.
		        DO k = 1, NI
		            IF (infected(k) == neighbor) THEN
		                ! WRITE(*,*)"infected", infected(k), neighbor
		                already_infected = .TRUE.
		                EXIT
		            END IF
		        END DO

		        old_pos = neigh_pos(j)
		        ! WRITE(*,*)old_pos,Eact
		        IF (already_infected) THEN
		            ! REMOVE link if it is active
		            IF (old_pos <= Eact) THEN
		                ! Swap with the last active link (Eact)
		                DO n = 1, N_links
		                    IF (neigh_pos(n) == Eact) THEN
		                        neigh_pos(n) = old_pos
		                        EXIT
		                    END IF
		                END DO
		                neigh_pos(j) = Eact

		                ! Swap active link rows
		                old_link_value(:) = Wact(Eact,:)
		                Wact(Eact,:) = Wact(old_pos,:)
		                Wact(old_pos,:) = old_link_value(:)

		                Eact = Eact - 1
		            END IF

		        ELSE
		            ! ADD link if not already active
		            IF (old_pos > Eact) THEN
		                Eact = Eact + 1

		                neigh_pos(j) = Eact
		                DO n = 1, N_links
		                    IF (neigh_pos(n) == Eact .AND. n /= j) THEN
		                        neigh_pos(n) = old_pos
		                        EXIT
		                    END IF
		                END DO

		                ! Swap active link rows
		                old_link_value(:) = Wact(Eact,:)
		                Wact(Eact,:) = Wact(old_pos,:)
		                Wact(old_pos,:) = old_link_value(:)
		            END IF
		        END IF
		        IF (Eact < 1 .OR. Eact > N_links) THEN
				    WRITE(*,*)"Eact out of bounds!", Eact, N_links
				    STOP
				END IF
		    END DO
		END DO

		! IF (step==394) THEN  
		! 	WRITE(*,*)"Active links:",Eact
		! END IF
		! DO n = 1, N_links
		! 	WRITE(*,*)Wact(n,:)
		! 	IF (ios /= 0) EXIT
		! END DO 

		RN3 = r1279()
		t = t - LOG(RN3)/lambda_tot
		NI_perc = float(NI)/N_nodes
		WRITE(11,'(F10.5, F10.5)')t,NI_perc
	END DO
	CLOSE(11)
	WRITE(*,*)"File produced"
END DO

END PROGRAM SIS_dyn

	  FUNCTION r1279()

	  IMPLICIT NONE
	  INCLUDE "r1279block.h"
	  REAL    r1279, inv_max
	  REAL    INV_MAXINT 
	  PARAMETER (INV_MAXINT = 1.0/2147483647.0)

	  ioffset = iand(ioffset + 1, 2047)
	  irand(ioffset) = (irand(index1(ioffset))*irand(index2(ioffset)))
	  r1279 = ishft(irand(ioffset), -1) * INV_MAXINT

	  END 
	        
	  SUBROUTINE setr1279(iseed)

	  IMPLICIT NONE
	  INCLUDE "r1279block.h"
	  INTEGER ibit, ispoke, one_bit, iseed, localseed, NBITM1
	  REAL    ran2
	  PARAMETER (NBITM1 = 31)
	  !
	  !     Initialize ioffset. This will be increased by (1 mod 2048) for
	  !     each random number which is called. 
	  !
	  ioffset = 0
	  !
	  !     Set up the two arrays which give locations of the two random
	  !     numbers which will be multiplied to get the new random number
	  !
	  do ispoke = 0, 2047
	    index1(ispoke) = iand(ispoke - 1279, 2047)
	    index2(ispoke) = iand(ispoke - 418, 2047)
	  end do
	  !
	  ! set up the initial array of 2048 integer random numbers
	  ! Each bit is separately initialized using ran2 from numerical recipes
	  !
	  localseed = -abs(iseed)


	  ! Matteo:  I have substituted lshift with ishft


	      do ispoke = 0, 2047

	        irand(ispoke) = 0
	        do ibit = 0, NBITM1
	            one_bit = 0
	            if (ran2(localseed) > 0.5) one_bit = 1
	            irand(ispoke) = ior(irand(ispoke), ishft(one_bit, ibit))
	        end do
	        irand(ispoke) = 2 * irand(ispoke) + 1

	      end do

	  END 

	  FUNCTION ran2(idum)
	  INTEGER idum,IM1,IM2,IMM1,IA1,IA2,IQ1,IQ2,IR1,IR2,NTAB,NDIV
	  REAL ran2,AM,EPS,RNMX
	  PARAMETER (IM1=2147483563,IM2=2147483399,AM=1./IM1,IMM1=IM1-1,&
	  IA1=40014,IA2=40692,IQ1=53668,IQ2=52774,IR1=12211,IR2=3791,&
	  NTAB=32,NDIV=1+IMM1/NTAB,EPS=1.2e-7,RNMX=1.-EPS)
	  INTEGER idum2,j,k,iv(NTAB),iy
	  SAVE iv,iy,idum2
	  DATA idum2/123456789/, iv/NTAB*0/, iy/0/
	  if (idum.le.0) then
	    idum=max(-idum,1)
	    idum2=idum
	    do 11 j=NTAB+8,1,-1
	      k=idum/IQ1
	      idum=IA1*(idum-k*IQ1)-k*IR1
	      if (idum.lt.0) idum=idum+IM1
	      if (j.le.NTAB) iv(j)=idum
11      continue
	    iy=iv(1)
	  endif
	  k=idum/IQ1
	  idum=IA1*(idum-k*IQ1)-k*IR1
	  if (idum.lt.0) idum=idum+IM1
	  k=idum2/IQ2
	  idum2=IA2*(idum2-k*IQ2)-k*IR2
	  if (idum2.lt.0) idum2=idum2+IM2
	  j=1+iy/NDIV
	  iy=iv(j)-idum2
	  iv(j)=idum
	  if(iy.lt.1)iy=iy+IMM1
	  ran2=min(AM*iy,RNMX)
	  return
  	  END