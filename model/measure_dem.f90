
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
	double precision :: b_out, vrelf_unit(3), rpos12(3), rposf_out(3), vrelf_norm
	double precision :: delta_Et_inel, delta_Et_el, delta_E_diss, f_tr
	logical :: keep_hit

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

	! Post-Collision Energy
	VCOM = (VEL(2,:)+VEL(1,:))*0.5D0
	V1COM = VEL(1,:) - VCOM; V2COM = VEL(2,:) - VCOM
	Ef(1) = MASS*(DOT_PRODUCT(V1COM,V1COM) + DOT_PRODUCT(V2COM,V2COM))
	Ef(2) = MOI(2)*(DOT_PRODUCT(OMEGA(1,2:3),OMEGA(1,2:3)) + &
		DOT_PRODUCT(OMEGA(2,2:3),OMEGA(2,2:3)))

	! --- Elastic pass: store reference energy and return early ---
	IF (ELASTIC_PASS) THEN
		Et_f_elastic = Ef(1)
		Er_f_elastic = Ef(2)
		RETURN
	END IF

	! --- Inelastic pass below ---

	Er10 = Er_1
	Er20 = Er_2
	Er1_prime = MOI(2)*DOT_PRODUCT(OMEGA(1,2:3),OMEGA(1,2:3))
	Er2_prime = MOI(2)*DOT_PRODUCT(OMEGA(2,2:3),OMEGA(2,2:3))

	! Post-collision orientation vectors
	U1_post = U(1,:)
	U2_post = U(2,:)

	! f_tr: translational dissipation fraction
	! Uses Et_f_elastic set during the preceding elastic pass
	delta_Et_inel = Ef(1) - Et_00
	delta_Et_el   = Et_f_elastic - Et_00
	delta_E_diss  = E0 - (Ef(1) + Ef(2))
	IF (ABS(delta_E_diss) > 1.0D-30) THEN
		f_tr = (delta_Et_inel - delta_Et_el) / delta_E_diss
	ELSE
		f_tr = -999.0D0   ! sentinel: elastic collision, f_tr undefined
	END IF

	! b_out: outgoing impact parameter (corrected b_contact)
	! Computed from post-collision outgoing velocities/positions
	vrelf_norm = SQRT(DOT_PRODUCT(VRELF, VRELF))
	IF (vrelf_norm > 1.0D-30) THEN
		vrelf_unit = VRELF / vrelf_norm
		rpos12     = POS(2,:) - POS(1,:)
		rposf_out  = rpos12 - DOT_PRODUCT(rpos12, vrelf_unit) * vrelf_unit
		b_out      = SQRT(DOT_PRODUCT(rposf_out, rposf_out)) / (LCYL + DIA)
	ELSE
		b_out = 0.0D0
	END IF

	keep_hit = .FALSE.
!$OMP CRITICAL(hit_count)
	IF (NHIT < NSAMPLES) THEN
		NHIT = NHIT + 1
		keep_hit = .TRUE.
		IF (NHIT >= NSAMPLES) SIM_CONTINUE = .FALSE.
	END IF
!$OMP END CRITICAL(hit_count)
	IF (.NOT.keep_hit) RETURN

	! Buffer all outputs
	buffer_idx       = buffer_idx       + 1
	buffer_ftr_idx   = buffer_ftr_idx   + 1
	buffer_orient_idx = buffer_orient_idx + 1
	buffer_uvec_idx  = buffer_uvec_idx  + 1

	chi_buffer(buffer_idx, 1) = b_impact
	chi_buffer(buffer_idx, 2) = chi
	chi_buffer(buffer_idx, 3) = psi
	chi_buffer(buffer_idx, 4) = mu_in

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

	ftr_buffer(buffer_ftr_idx, 1) = f_tr
	ftr_buffer(buffer_ftr_idx, 2) = delta_Et_el
	ftr_buffer(buffer_ftr_idx, 3) = delta_E_diss

	orient_buffer(buffer_orient_idx,  1) = S2_pair
	orient_buffer(buffer_orient_idx,  2) = S2_1n
	orient_buffer(buffer_orient_idx,  3) = S2_2n
	orient_buffer(buffer_orient_idx,  4) = S2_1v
	orient_buffer(buffer_orient_idx,  5) = S2_2v
	orient_buffer(buffer_orient_idx,  6) = cos_u1_n
	orient_buffer(buffer_orient_idx,  7) = cos_u2_n
	orient_buffer(buffer_orient_idx,  8) = cos_u1_v
	orient_buffer(buffer_orient_idx,  9) = cos_u2_v
	orient_buffer(buffer_orient_idx, 10) = u1u2_dot
	orient_buffer(buffer_orient_idx, 11) = contact_lambda
	orient_buffer(buffer_orient_idx, 12) = contact_mu
	orient_buffer(buffer_orient_idx, 13) = E_n_pre
	orient_buffer(buffer_orient_idx, 14) = b_out

	uvec_buffer(buffer_uvec_idx,  1:3)  = U1_pre
	uvec_buffer(buffer_uvec_idx,  4:6)  = U2_pre
	uvec_buffer(buffer_uvec_idx,  7:9)  = U1_post
	uvec_buffer(buffer_uvec_idx, 10:12) = U2_post

	! Flush if buffer full
	IF (buffer_idx >= MAX_BUFFER) THEN
		CALL FLUSH_BUFFERS()
	END IF
	return
	end subroutine measure_dem
