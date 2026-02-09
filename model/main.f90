
	program main

	use run_param
	use particles
	use output
	use constants
!$ use omp_lib

	implicit none
	DOUBLE PRECISION, DIMENSION(3) :: R12, E12
	DOUBLE PRECISION :: D12, RV12
        DOUBLE PRECISION :: yaw, pitch, roll, q0, q1, q2, q3
        INTEGER :: I
	real :: xint
	INTEGER :: num_threads

	write(*,*) 'Reading and initializing'
	call read_input()
	call INITIALIZE()

	! Initialize buffer indices
!$OMP PARALLEL
	buffer_idx = 0
!$OMP END PARALLEL

	! Report thread count
	num_threads = 1
!$ num_threads = omp_get_max_threads()
	write(*,*) 'Using ', num_threads, ' thread(s)'

	write(*,*) 'Beginning collisions'

!$OMP PARALLEL PRIVATE(R12, E12, D12, RV12) &
!$OMP SHARED(NTRY, NHIT, SIM_CONTINUE, NSAMPLES)
	DO WHILE(SIM_CONTINUE)
		CALL INIT_PART()
		DO WHILE(.TRUE.)
			PREV_CONTACT = CONTACT
			CALL INTEGRATE_EOM
			IF(.NOT.PREV_CONTACT.AND.CONTACT) NPHIT = NPHIT + 1 

			! Exit Conditions
			R12 = POS(2,:) - POS(1,:)
			D12 = SQRT(DOT_PRODUCT(R12,R12))
			E12 = VEL(2,:) - VEL(1,:)
			E12 = E12/SQRT(DOT_PRODUCT(E12,E12))
			RV12 = DOT_PRODUCT(E12,R12)

			!write(2925, '(I1,/)') 2
			!DO I = 1, 2
			!     ! Calculate Euler angles from orientation vector
	 		!     ! Ovito's body frame x-axis is in the negative direction of z-x-y
	                !     yaw = atan2(U(I,2), U(I,1))
	                !     pitch = atan2(-U(I,3), sqrt(U(I,1)**2 + U(I,2)**2))
	                !     roll = 0.0d0

	                !     ! Convert Euler angles to quaternion
	                !     q0 = cos(yaw/2) * cos(pitch/2) * cos(roll/2) + sin(yaw/2) * sin(pitch/2) * sin(roll/2)
    	                !     q1 = sin(yaw/2) * cos(pitch/2) * cos(roll/2) - cos(yaw/2) * sin(pitch/2) * sin(roll/2)
    	                !     q2 = cos(yaw/2) * sin(pitch/2) * cos(roll/2) + sin(yaw/2) * cos(pitch/2) * sin(roll/2)
    	                !     q3 = cos(yaw/2) * cos(pitch/2) * sin(roll/2) - sin(yaw/2) * sin(pitch/2) * cos(roll/2)
                        !     xint = real(i)
	 
                        !     write(2925,'(18(E14.8,2X))') xint, POS(I,1), POS(I,2), POS(I,3), VEL(I,1), VEL(I,2), VEL(I,3), &
			!		DIA*0.5D0, DIA*0.5D0, 0.0D0, lcyl, q1, q2, q3, -q0, OMEGA(I,1), OMEGA(I,2), OMEGA(I,3)        
	                !END DO
		        !flush(2925)
			IF(RV12.GE.0.D0.AND.D12.GT.BMAX&
				.AND.SUM(F).LT.SMALL_NUM.AND.&
				SUM(TAU).LT.SMALL_NUM) EXIT
		END DO
		IF(HIT) CALL MEASURE_DEM

!$OMP ATOMIC
		NTRY = NTRY + 1
!$OMP END ATOMIC

		! Progress update (only from master thread)
!$OMP MASTER
		IF(MOD(NTRY, 1000) == 0) write(*,*) 'Number of samples: ', NTRY
!$OMP END MASTER
	END DO
!$OMP END PARALLEL

	write(*,*) 'Number of hits: ', NHIT-1
	
		
	write(*,*) 'Recording Final State'
	call measure_final()

	return
	end program main
