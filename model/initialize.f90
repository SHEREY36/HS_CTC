	subroutine initialize
	
	use particles
	use run_param
	use constants
	use output
	
	implicit none
	INTEGER :: I
	DOUBLE PRECISION :: P

! Output Files
!---------------------------------------------------------
	CHARACTER(LEN=300) :: filepath

	! Create output directory if it does not exist
	CALL EXECUTE_COMMAND_LINE('mkdir -p ' // TRIM(output_dir))

	filepath = TRIM(output_dir) // '/csx.txt'
	open(unit=1000, status='replace', file=filepath)

	filepath = TRIM(output_dir) // '/chi.txt'
	open(unit=1001, status='replace', file=filepath)

	filepath = TRIM(output_dir) // '/Ef.txt'
	open(unit=1002, status='replace', file=filepath)

	filepath = TRIM(output_dir) // '/EnergyCons.txt'
	open(unit=1003, status='replace', file=filepath)

	filepath = TRIM(output_dir) // '/NPhit.txt'
	open(unit=1004, status='replace', file=filepath)

	filepath = TRIM(output_dir) // '/projArea.txt'
	open(unit=2000, status='replace', file=filepath)

	filepath = TRIM(output_dir) // '/init.txt'
	open(unit=100,  status='replace', file=filepath)

	filepath = TRIM(output_dir) // '/PreRotEnergy.txt'
	open(unit=1111, status='replace', file=filepath)

	!open(unit=2925, status='replace', file='ovito.txt')

! Particle Properties
!---------------------------------------------------------
	! Particle Geometry/Material
	RAD = DIA*0.5D0; hLCYL = LCYL*0.5D0
	DIASQ = DIA**2.D0
	BMAX = LCYL + DIA; BMAX = BMAX*1.01D0
	PVOL = PI*(DIA**3.D0)/6.D0 + PI*LCYL*(RAD**2.D0)
	RHO = MASS/PVOL
	
	! Moment of Inertia
	moI(2) = PI/48.D0*RHO*(DIA**2.D0)*(LCYL**3.D0) +&
		3.D0*PI/64.D0*RHO*(DIA**4.D0)*LCYL +&
		PI/60.D0*RHO*(DIA**5.D0) +&
		PI/24.D0*RHO*(DIA**3.D0)*(LCYL**2.D0)
	moI(3) = moI(2)
	moI(1) = PI/32.D0*RHO*(DIA**4.D0)*LCYL + &
		PI/60.D0*RHO*(DIA**5.D0)
	DO I = 1,3
		omoI(I) = 1.D0/moI(I)
	END DO

	! Hertzian Spring
	EStar = EYoung/(2.D0*(1.D0-(GPoisson**2.D0)))
	KN = 4.D0/3.D0*SQRT(RAD*0.5D0)*EStar
	! Damper
	BETA = -LOG(ALPHA_PP)/(SQRT((PI**2.D0)+(LOG(ALPHA_PP))**2.D0))
	CN = 2.D0*BETA*SQRT(MASS*0.5D0*KN)
	! Time Step
	TCOLL = PI/SQRT(KN/(MASS*0.5D0) - ((CN/(MASS*0.5D0))**2.D0)/4.D0)
	!3.21D0*(MASS*0.5D0/KN)**(0.4D0)

! Sampling Temperatures
!---------------------------------------------------------
	! STD = kT/m / kT/I
	IF(ToE.EQ.'E') THEN
		SQEk = SQRT(Ek/MASS); SQEr = SQRT(Er*OMOI(2))
	ELSE
		kTm = kTm/MASS; kTI = kTI*OMOI(2)
		okTm= 1.D0/kTm; okTI = 1.D0/kTI
		VMAX = 2.D0
		DO WHILE(.TRUE.)
			P = (VMAX**3.D0)*EXP(-(VMAX**2.D0)*0.25D0)
			IF(P.LT.1.0D-4) EXIT
			VMAX = VMAX * 1.05D0
		END DO
		WMAX = 2.D0
		DO WHILE(.TRUE.)
			P = EXP(-(WMAX**2.D0)*0.5D0)
			IF(P.LT.1.0D-4) EXIT
			WMAX = WMAX * 1.05D0
		END DO
	END IF

! For Outputs
!---------------------------------------------------------
	NHIT = 1
	TMEAN = 0.D0; RMEAN = 0.D0
	SIM_CONTINUE = .TRUE.
	return
	end subroutine initialize
