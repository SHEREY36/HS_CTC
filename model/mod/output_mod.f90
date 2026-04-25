
	module output

	implicit none

	INTEGER :: NHIT
	LOGICAL :: HIT

	DOUBLE PRECISION :: CSX
	DOUBLE PRECISION, DIMENSION(3) :: VREL0, WREL0
	DOUBLE PRECISION :: PROJ_AREA

	DOUBLE PRECISION :: E0, Er_00, Et_00, Er_1, Er_2
	DOUBLE PRECISION :: TMEAN, RMEAN
	DOUBLE PRECISION :: b_impact, b_contact

	! Elastic-pass energy storage (for f_tr calculation)
	DOUBLE PRECISION :: Et_f_elastic, Er_f_elastic

	! Contact geometry (saved at first contact, normalised by hLCYL)
	DOUBLE PRECISION :: contact_lambda, contact_mu
	DOUBLE PRECISION :: mu_in

	! Nematic order scalars at first contact
	DOUBLE PRECISION :: S2_pair, S2_1n, S2_2n, S2_1v, S2_2v
	DOUBLE PRECISION :: cos_u1_n, cos_u2_n, cos_u1_v, cos_u2_v
	DOUBLE PRECISION :: u1u2_dot

	! Pre/post orientation vectors
	DOUBLE PRECISION, DIMENSION(3) :: U1_pre, U2_pre, U1_post, U2_post

	! Normal-direction translational energy at contact
	DOUBLE PRECISION :: E_n_pre

	LOGICAL :: SIM_CONTINUE

	! Output directory for files (shared - read-only after initialization)
	CHARACTER(LEN=255) :: output_dir

	! Per-collision output variables: each thread tracks its own collision
	!$OMP THREADPRIVATE(HIT, VREL0, WREL0, PROJ_AREA, E0, Er_00, Et_00,        &
	!$OMP&             Er_1, Er_2, TMEAN, RMEAN, b_impact, b_contact,            &
	!$OMP&             Et_f_elastic, Er_f_elastic,                                &
	!$OMP&             contact_lambda, contact_mu, mu_in,                         &
	!$OMP&             S2_pair, S2_1n, S2_2n, S2_1v, S2_2v,                      &
	!$OMP&             cos_u1_n, cos_u2_n, cos_u1_v, cos_u2_v, u1u2_dot,         &
	!$OMP&             U1_pre, U2_pre, U1_post, U2_post, E_n_pre)

	! Elastic-pass flag: when TRUE, calc_force_dem uses CN=0
	LOGICAL :: ELASTIC_PASS
	!$OMP THREADPRIVATE(ELASTIC_PASS)

	! Output buffers (thread-private: each thread accumulates its own data)
	INTEGER, PARAMETER :: MAX_BUFFER = 1000
	DOUBLE PRECISION, DIMENSION(MAX_BUFFER, 4) :: chi_buffer
	DOUBLE PRECISION, DIMENSION(MAX_BUFFER, 7) :: ef_buffer
	DOUBLE PRECISION, DIMENSION(MAX_BUFFER)    :: econs_buffer
	INTEGER,          DIMENSION(MAX_BUFFER)    :: nphit_buffer
	DOUBLE PRECISION, DIMENSION(MAX_BUFFER, 2) :: prerot_buffer
	INTEGER :: buffer_idx
	!$OMP THREADPRIVATE(chi_buffer, ef_buffer, econs_buffer, nphit_buffer, &
	!$OMP&             prerot_buffer, buffer_idx)

	! ftr buffer: [f_tr, delta_Et_el, delta_E_diss] per collision
	DOUBLE PRECISION, DIMENSION(MAX_BUFFER, 3) :: ftr_buffer
	INTEGER :: buffer_ftr_idx
	!$OMP THREADPRIVATE(ftr_buffer, buffer_ftr_idx)

	! orient_buffer (14 cols):
	!   1: S2_pair  2: S2_1n   3: S2_2n   4: S2_1v   5: S2_2v
	!   6: cos_u1_n 7: cos_u2_n 8: cos_u1_v 9: cos_u2_v 10: u1u2_dot
	!  11: lambda_c 12: mu_c   13: E_n_pre  14: b_out
	DOUBLE PRECISION, DIMENSION(MAX_BUFFER, 14) :: orient_buffer
	INTEGER :: buffer_orient_idx
	!$OMP THREADPRIVATE(orient_buffer, buffer_orient_idx)

	! uvec_buffer (12 cols): U1_pre(3), U2_pre(3), U1_post(3), U2_post(3)
	DOUBLE PRECISION, DIMENSION(MAX_BUFFER, 12) :: uvec_buffer
	INTEGER :: buffer_uvec_idx
	!$OMP THREADPRIVATE(uvec_buffer, buffer_uvec_idx)

	contains

	SUBROUTINE FLUSH_BUFFERS()
		INTEGER :: i
		IF (buffer_idx > 0) THEN
!$OMP CRITICAL(file_write)
			DO i = 1, buffer_idx
				write(1001,'(4(E14.8,2X))') chi_buffer(i,:)
				write(1002,'(7(E14.8,2X))') ef_buffer(i,:)
				write(1003,'(E14.8)')        econs_buffer(i)
				write(1004,*)                nphit_buffer(i)
				write(1111,'(2(E14.8,2X))') prerot_buffer(i,:)
			END DO
			flush(1001); flush(1002); flush(1003); flush(1004); flush(1111)
!$OMP END CRITICAL(file_write)
			buffer_idx = 0
		END IF
		CALL FLUSH_FTR_BUFFER()
		CALL FLUSH_ORIENT_BUFFER()
		CALL FLUSH_UVEC_BUFFER()
	END SUBROUTINE FLUSH_BUFFERS

	SUBROUTINE FLUSH_FTR_BUFFER()
		INTEGER :: i
		IF (buffer_ftr_idx > 0) THEN
!$OMP CRITICAL(ftr_write)
			DO i = 1, buffer_ftr_idx
				write(1005,'(3(E14.8,2X))') ftr_buffer(i,:)
			END DO
			flush(1005)
!$OMP END CRITICAL(ftr_write)
			buffer_ftr_idx = 0
		END IF
	END SUBROUTINE FLUSH_FTR_BUFFER

	SUBROUTINE FLUSH_ORIENT_BUFFER()
		INTEGER :: i
		IF (buffer_orient_idx > 0) THEN
!$OMP CRITICAL(orient_write)
			DO i = 1, buffer_orient_idx
				write(1006,'(14(E14.8,2X))') orient_buffer(i,:)
			END DO
			flush(1006)
!$OMP END CRITICAL(orient_write)
			buffer_orient_idx = 0
		END IF
	END SUBROUTINE FLUSH_ORIENT_BUFFER

	SUBROUTINE FLUSH_UVEC_BUFFER()
		INTEGER :: i
		IF (buffer_uvec_idx > 0) THEN
!$OMP CRITICAL(uvec_write)
			DO i = 1, buffer_uvec_idx
				write(1007,'(12(E14.8,2X))') uvec_buffer(i,:)
			END DO
			flush(1007)
!$OMP END CRITICAL(uvec_write)
			buffer_uvec_idx = 0
		END IF
	END SUBROUTINE FLUSH_UVEC_BUFFER

	end module output
