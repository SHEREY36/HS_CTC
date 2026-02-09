
	SUBROUTINE INIT_PART

! Global variables
!------------------------------------------------------------------//
	use particleS
	USE CONSTANTS
	use run_param
	use output
!$ use omp_lib

	implicit none
	INTEGER :: I, J
	DOUBLE PRECISION :: RR
	DOUBLE PRECISION :: PHI, THETA
	DOUBLE PRECISION :: VS, pV
	DOUBLE PRECISION :: WS, pW
	DOUBLE PRECISION :: spdmax

	double precision, dimension(3) :: v1com, v2com, vcom
	double precision, dimension(3) :: RVEL, RPOS

	! Thread-safe random number generation
	INTEGER :: thread_id, seed_size
	INTEGER, ALLOCATABLE :: seed(:)

	! Get thread ID (0 in serial mode)
	thread_id = 0
!$ thread_id = omp_get_thread_num()

	! Initialize thread-specific random seed
	CALL RANDOM_SEED(size=seed_size)
	ALLOCATE(seed(seed_size))
	seed = 12345 + thread_id * 100000 + NTRY  ! Unique per thread/iteration
	CALL RANDOM_SEED(put=seed)
	DEALLOCATE(seed)

	! Sample Position
	POS(1,:) = (/0.D0, 0.D0, 0.D0/)
	POS(2,1) = BMAX
	DO I = 2,3
		CALL RANDOM_NUMBER(RR)
		POS(2,I) = BMAX*RR
	END DO
	
	! Sample Orientation Vectors
	DO I = 1,2
		CALL RANDOM_NUMBER(RR); PHI = ACOS(RR*2.D0 - 1.D0)
		CALL RANDOM_NUMBER(RR); THETA = 2.D0*RR*PI
		U(I,:) = (/SIN(PHI)*COS(THETA), SIN(PHI)*SIN(THETA), COS(PHI)/)
		U(I,:) = U(I,:)/SQRT(DOT_PRODUCT(U(I,:),U(I,:)))
		CALL RANDOM_NUMBER(RR); PHI = ACOS(RR*2.D0 - 1.D0)
		CALL RANDOM_NUMBER(RR); THETA = 2.D0*RR*PI
		UX(I,:) = (/SIN(PHI)*COS(THETA), SIN(PHI)*SIN(THETA), COS(PHI)/)
		UX(I,:) = UX(I,:) - DOT_PRODUCT(UX(I,:),U(I,:))*U(I,:)
		UX(I,:) = UX(I,:)/SQRT(DOT_PRODUCT(UX(I,:),UX(I,:)))
		UY(I,:) = CROSSPRDCT(UX(I,:),U(I,:))
		UY(I,:) = UY(I,:)/SQRT(DOT_PRODUCT(UY(I,:),UY(I,:)))
	END DO

! Constant Temperature	
	! Sample Velocity from distribution of relative velocity of colliding particles
	VEL(2,:) = (/0.D0, 0.D0, 0.D0/)
	VEL(1,2) = 0.D0; VEL(1,3) = 0.D0
	OMEGA(1,1) = 0.D0; OMEGA(2,1) = 0.D0
	IF(ToE.EQ.'T') THEN
		DO WHILE(.TRUE.)
			CALL RANDOM_NUMBER(RR)
			VS = RR*VMAX; pV = (VS**3.D0)*EXP(-(VS**2.D0)*0.25D0)
			CALL RANDOM_NUMBER(RR)
			IF(pV.GT.VMAX) THEN
				write(*,*) '>vmax'
				stop
			END IF
			IF(pV.GT.RR*VMAX) EXIT
		END DO
		VEL(1,1) = VS*SQRT(kTm)
		! Sample Angular Velocities
		DO I = 1,2
			J = 1
			DO WHILE(J.LE.2)
				CALL RANDOM_NUMBER(RR)
				WS = 2.D0*RR*WMAX-WMAX; pW = exp(-(ws**2.D0)*0.5D0)
				CALL RANDOM_NUMBER(RR)
				IF(pW.GT.RR*WMAX) THEN
					OMEGA(I,J+1) = WS*SQRT(kTI)
					J = J + 1
				END IF
			END DO
		END DO
	ELSE IF(ToE.EQ.'E') THEN
		CALL RANDOM_NUMBER(RR)
		OMEGA(1,2) = SQEr*RR; OMEGA(1,3) = SQEr*(1.D0-RR)
		CALL RANDOM_NUMBER(RR)
		OMEGA(2,2) = SQEr*RR; OMEGA(2,3) = SQEr*(1.D0-RR)
		VEL(1,1) = SQEk
	ELSE
		write(*,*) 'Particle init not defined'
		stop
	END IF
	
	! Initialize forces
	F(:,:) = 0.D0
	TAU(:,:) = 0.D0
	
	! Record initial conditions
	UX(I,:) = UX(I,:) - DOT_PRODUCT(UX(I,:),U(I,:))*U(I,:)
	RPOS = POS(2,:) - POS(1,:); RVEL = VEL(2,:) - VEL(1,:)
	RPOS = RPOS - DOT_PRODUCT(RVEL,RPOS)&
		*RVEL/DOT_PRODUCT(RVEL,RVEL)
	b_impact = SQRT(DOT_PRODUCT(RPOS,RPOS))
	!TMEAN = TMEAN + MASS*(VEL(1,1)*VEL(1,1))
	!RMEAN = RMEAN + MOI(2)*(OMEGA(1,2)*OMEGA(1,2)+&
	!	OMEGA(1,3)*OMEGA(1,3)+&
	!	OMEGA(2,2)*OMEGA(2,2)+&
	!	OMEGA(2,3)*OMEGA(2,3))
	!write(100,'(2(D10.4,2X))') MASS*(VEL(1,1)*VEL(1,1)),&
	!	MOI(2)*(OMEGA(1,2)*OMEGA(1,2)+OMEGA(1,3)*OMEGA(1,3)+&
	!	OMEGA(2,2)*OMEGA(2,2)+OMEGA(2,3)*OMEGA(2,3))
	!flush(100)
	!write(*,*) TMEAN/(3.D0*DBLE(NTRY+1)), RMEAN/(2.D0*DBLE(NTRY+1))

	VCOM = (VEL(1,:) + VEL(2,:))*0.5D0
	V1COM = VEL(1,:) - VCOM; V2COM = VEL(2,:) - VCOM
	E0 = MASS*(DOT_PRODUCT(V1COM,V1COM) + DOT_PRODUCT(V2COM,V2COM))
	E0 = E0 + &
		moI(2)*(DOT_PRODUCT(OMEGA(1,:),OMEGA(1,:)) +&
		DOT_PRODUCT(OMEGA(2,:),OMEGA(2,:)))
        Er_00 = moI(2)*(DOT_PRODUCT(OMEGA(1,:),OMEGA(1,:)) +&
		DOT_PRODUCT(OMEGA(2,:),OMEGA(2,:)))
	Er_1 = moI(2)*DOT_PRODUCT(OMEGA(1,:),OMEGA(1,:))
	Er_2 = moI(2)*DOT_PRODUCT(OMEGA(2,:),OMEGA(2,:))
        Et_00 = MASS*(DOT_PRODUCT(V1COM,V1COM) + DOT_PRODUCT(V2COM,V2COM))
	
	
	! Max relative speed at impact
	! Additional contributions from rotation
	spdmax = SQRT(DOT_PRODUCT(VEL(1,:)-VEL(2,:),VEL(1,:)-VEL(2,:)))
	spdmax = spdmax + &
		(SQRT(DOT_PRODUCT(OMEGA(1,:),OMEGA(1,:))) +&
		SQRT(DOT_PRODUCT(OMEGA(2,:),OMEGA(2,:))))*(LCYL+DIA)*0.5D0
	! Factor of 100 arbitrary
	dt = TCOLL/50.D0
        !*(SPDMAX**(-0.2D0))
	! For calls to outputs
	VREL0 = VEL(2,:) - VEL(1,:)
	WREL0 = OMEGA(2,:) - OMEGA(1,:)
	HIT = .FALSE.
	CONTACT = .FALSE.; NPHIT = 0
	RETURN
	end subroutine INIT_PART

