
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
	
        Er10 = Er1
        Er20 = Er2

	Er1_prime = MOI(2)*DOT_PRODUCT(OMEGA(1,2:3),OMEGA(1,2:3))
	Er2_prime = MOI(2)*DOT_PRODUCT(OMEGA(2,2:3),OMEGA(2,2:3))

	write(1001,'(3(E14.8,2X))') b_impact, chi, psi
	write(1002,'(7(E14.8,2X))') Et_00, Er10, Er20, Ef(1), Er1_prime, Er2_prime, b_contact
	write(1003, '(E14.8)') SUM(Ef)/E0
	write(1004, *) NPHIT
        write(1111, '(2(E14.8,2X))') sqrt(DOT_PRODUCT(VREL0,VREL0)), sqrt(DOT_PRODUCT(VRELF,VRELF)) 
	flush(1001); flush(1002); flush(1003); flush(1004); flush(1111)

	!write(1003,'(E14.8)') proj_area
	NHIT = NHIT + 1
	IF(NHIT.GT.NSAMPLES) SIM_CONTINUE = .False.
	return	
	end subroutine measure_dem
