	subroutine initialize
	
	use particles
	use run_param
	use constants
	use output
	
	implicit none
	INTEGER :: I

! Output Files
!---------------------------------------------------------
	CHARACTER(LEN=300) :: filepath

	! Create output directory if it does not exist
	CALL EXECUTE_COMMAND_LINE('mkdir -p ' // TRIM(output_dir))

	IF (WRITE_LEGACY) THEN
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

		filepath = TRIM(output_dir) // '/ftr_data.txt'
		open(unit=1005, status='replace', file=filepath)

		filepath = TRIM(output_dir) // '/orient.dat'
		open(unit=1006, status='replace', file=filepath)

		filepath = TRIM(output_dir) // '/uvec.dat'
		open(unit=1007, status='replace', file=filepath)
	END IF

	IF (WRITE_CLOSURE) THEN
		filepath = TRIM(output_dir) // '/closure_events.bin'
		open(unit=1010, status='replace', file=filepath, access='stream', &
			form='unformatted', convert='little_endian')
	END IF

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
	! Use one alpha-independent integration step for common random numbers along
	! an alpha line. The elastic Hertz contact time is also the conservative
	! reference replay time scale; damping only changes the force law.
	TCOLL = PI/SQRT(KN/(MASS*0.5D0))

! Sampling Temperatures
!---------------------------------------------------------
	! STD = kT/m / kT/I
	IF(ToE.EQ.'E') THEN
		SQEk = SQRT(Ek/MASS); SQEr = SQRT(Er*OMOI(2))
	ELSE
		kTm = kTm/MASS; kTI = kTI*OMOI(2)
	END IF

	IF (WRITE_CLOSURE) THEN
		filepath = TRIM(output_dir) // '/metadata.txt'
		open(unit=1011, status='replace', file=filepath)
		write(1011,'(A)') 'schema_version=1'
		write(1011,'(A)') 'dtype=<f8'
		write(1011,'(A,I0)') 'n_columns=', N_CLOSURE_COLS
		write(1011,'(A)') 'columns=event_id,c1_x,c1_y,c1_z,c2_x,c2_y,c2_z,' // &
			'omega1_x,omega1_y,omega1_z,omega2_x,omega2_y,omega2_z,' // &
			'u1_x,u1_y,u1_z,u2_x,u2_y,u2_z,delta_tr,delta_rot,delta_total,' // &
			'et_elastic,er_elastic,et_inelastic,er_inelastic,e_initial,' // &
			'elastic_rel_error,n_contact,impact_parameter,contact_n_x,' // &
			'contact_n_y,contact_n_z,centerline_x,centerline_y,centerline_z,' // &
			'contact_lambda,contact_mu'
		write(1011,'(A,ES24.16)') 'alpha=', ALPHA_PP
		write(1011,'(A,ES24.16)') 'theta=', TTR_INPUT/TROT_INPUT
		write(1011,'(A,ES24.16)') 'aspect_ratio=', AR_INPUT
		write(1011,'(A,ES24.16)') 'temperature_translational=', TTR_INPUT
		write(1011,'(A,ES24.16)') 'temperature_rotational=', TROT_INPUT
		write(1011,'(A,ES24.16)') 'velocity_scale=', SQRT(2.D0*kTm)
		write(1011,'(A,ES24.16)') 'omega_scale=', SQRT(2.D0*kTI)
		write(1011,'(A,ES24.16)') 'mass=', MASS
		write(1011,'(A,ES24.16)') 'moi_perpendicular=', moI(2)
		write(1011,'(A,I0)') 'nsamples=', NSAMPLES
		write(1011,'(A,I0)') 'seed=', RUN_SEED
		write(1011,'(A,I0)') 'event_id_offset=', RUN_SEED*EVENT_ID_STRIDE
		write(1011,'(A)') 'rng_contract=event_id_stream_common_across_alpha'
		write(1011,'(A,A)') 'output_mode=', TRIM(OUTPUT_MODE)
		close(1011)
	END IF

! For Outputs
!---------------------------------------------------------
	NHIT = 0
	NTRY = 0
	TMEAN = 0.D0; RMEAN = 0.D0
	SIM_CONTINUE = .TRUE.
	return
	end subroutine initialize
