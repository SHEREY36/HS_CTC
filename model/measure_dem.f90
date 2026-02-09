
	subroutine measure_dem()
	
	use run_param
	use particles
	use output
	
	implicit none
	integer :: i, j
	double precision :: vrelf(3), wrelf(3)
	double precision, dimension(3) :: v1com, v2com, vcom
	double precision :: Ef(2)
	double precision :: chi, psi, Er10, Er20, Er1_prime, Er2_prime

	! Scattering Angle
	VRELF = VEL(2,:) - VEL(1,:)
	chi = SQRT(DOT_PRODUCT(VRELF,VRELF))
	chi = SQRT(DOT_PRODUCT(VREL0,VREL0))*chi
	chi = ACOS(DOT_PRODUCT(VRELF,VREL0)/chi)

	! Scattering Angle for angular velocity W
	WRELF = OMEGA(2,:) - OMEGA(1,:)
	psi = SQRT(DOT_PRODUCT(WRELF,WRELF))
	psi = SQRT(DOT_PRODUCT(WREL0,WREL0))*psi
	psi = ACOS(DOT_PRODUCT(WRELF,WREL0)/psi)
	
	! Post-Collisinal Energy
	VCOM = (VEL(2,:)+VEL(1,:))*0.5D0
	V1COM = VEL(1,:) - VCOM; V2COM = VEL(2,:) - VCOM
	Ef(1) = MASS*(DOT_PRODUCT(V1COM,V1COM) + DOT_PRODUCT(V2COM,V2COM))
	Ef(2) = MOI(2)*(DOT_PRODUCT(OMEGA(1,2:3),OMEGA(1,2:3)) + &
		DOT_PRODUCT(OMEGA(2,2:3),OMEGA(2,2:3)))
	
        Er10 = Er_1
        Er20 = Er_2

	Er1_prime = MOI(2)*DOT_PRODUCT(OMEGA(1,2:3),OMEGA(1,2:3))
	Er2_prime = MOI(2)*DOT_PRODUCT(OMEGA(2,2:3),OMEGA(2,2:3))

	! Buffer outputs instead of direct write
	buffer_idx = buffer_idx + 1

	chi_buffer(buffer_idx, 1) = b_impact
	chi_buffer(buffer_idx, 2) = chi
	chi_buffer(buffer_idx, 3) = psi

	ef_buffer(buffer_idx, 1) = Et_00
	ef_buffer(buffer_idx, 2) = Er10
	ef_buffer(buffer_idx, 3) = Er20
	ef_buffer(buffer_idx, 4) = Ef(1)
	ef_buffer(buffer_idx, 5) = Er1_prime
	ef_buffer(buffer_idx, 6) = Er2_prime
	ef_buffer(buffer_idx, 7) = b_contact

	econs_buffer(buffer_idx) = SUM(Ef)/E0
	nphit_buffer(buffer_idx) = NPHIT
	prerot_buffer(buffer_idx, 1) = sqrt(DOT_PRODUCT(VREL0,VREL0))
	prerot_buffer(buffer_idx, 2) = sqrt(DOT_PRODUCT(VRELF,VRELF))

	! Flush if buffer full
	IF (buffer_idx >= MAX_BUFFER) THEN
		CALL FLUSH_BUFFERS()
	END IF

	NHIT = NHIT + 1
	IF(NHIT.GT.NSAMPLES) SIM_CONTINUE = .False.
	return	
	end subroutine measure_dem
