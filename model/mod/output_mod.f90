
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

	LOGICAL :: SIM_CONTINUE

	! Output directory for files
	CHARACTER(LEN=255) :: output_dir

	! Output buffers (thread-private, will be used with OpenMP)
	INTEGER, PARAMETER :: MAX_BUFFER = 1000
	DOUBLE PRECISION, DIMENSION(MAX_BUFFER, 3) :: chi_buffer
	DOUBLE PRECISION, DIMENSION(MAX_BUFFER, 7) :: ef_buffer
	DOUBLE PRECISION, DIMENSION(MAX_BUFFER) :: econs_buffer
	INTEGER, DIMENSION(MAX_BUFFER) :: nphit_buffer
	DOUBLE PRECISION, DIMENSION(MAX_BUFFER, 2) :: prerot_buffer
	INTEGER :: buffer_idx

	contains

	SUBROUTINE FLUSH_BUFFERS()
		INTEGER :: i
		IF (buffer_idx > 0) THEN
			DO i = 1, buffer_idx
				write(1001,'(3(E14.8,2X))') chi_buffer(i,:)
				write(1002,'(7(E14.8,2X))') ef_buffer(i,:)
				write(1003,'(E14.8)') econs_buffer(i)
				write(1004,*) nphit_buffer(i)
				write(1111,'(2(E14.8,2X))') prerot_buffer(i,:)
			END DO
			flush(1001); flush(1002); flush(1003); flush(1004); flush(1111)
			buffer_idx = 0
		END IF
	END SUBROUTINE FLUSH_BUFFERS

	end module output	
