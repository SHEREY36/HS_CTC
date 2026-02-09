
	subroutine read_input

	use particles
	use run_param
	use output

	implicit none

	integer :: I
	integer :: lunit = 10
	character(len=255) :: filename_input
	double precision :: utemp, wtemp
	integer :: num_args
	character(len=100) :: arg_val

	! Check for command-line arguments
	num_args = COMMAND_ARGUMENT_COUNT()

	IF (num_args >= 3) THEN
		! Command-line mode: read ALPHA_PP, kTm, kTI from arguments
		write(*,*) 'Reading parameters from command line'

		! Read system_input.dat for geometry and material properties
		filename_input = 'system_input.dat'
		OPEN(unit=lunit, FILE=filename_input, status='old')
		READ(lunit,*) NSAMPLES
		READ(lunit,*) DIA, LCYL, MASS
		READ(lunit,*) EYoung, GPoisson
		close(lunit)

		! Override parameters from command line
		CALL GET_COMMAND_ARGUMENT(1, arg_val)
		READ(arg_val, *) ALPHA_PP

		CALL GET_COMMAND_ARGUMENT(2, arg_val)
		READ(arg_val, *) utemp
		kTm = utemp

		CALL GET_COMMAND_ARGUMENT(3, arg_val)
		READ(arg_val, *) wtemp
		kTI = wtemp

		ToE = 'T'  ! Force temperature mode for command-line

		! Optional 4th argument: output directory
		IF (num_args >= 4) THEN
			CALL GET_COMMAND_ARGUMENT(4, output_dir)
		ELSE
			output_dir = './'
		END IF

		write(*,*) 'ALPHA_PP = ', ALPHA_PP
		write(*,*) 'kTm = ', kTm
		write(*,*) 'kTI = ', kTI
		write(*,*) 'Output dir: ', TRIM(output_dir)

	ELSE
		! Legacy mode: read all from system_input.dat
		write(*,*) 'Reading parameters from system_input.dat'
		output_dir = './'

		filename_input	= 'system_input.dat'
		OPEN(unit=lunit, FILE=filename_input, status='old')
		READ(lunit,*) NSAMPLES
		READ(LUNIT,*) DIA, LCYL, MASS
		READ(lunit,*) EYoung, GPoisson
		READ(lunit,*) ALPHA_PP
		READ(lunit,*) ToE, utemp, wtemp
		IF(ToE.EQ.'E') THEN
			Ek = utemp;	Er = wtemp
		ELSE
			! Reading in kT (not kT/m) - will be adjusted in initialize
			kTm = utemp; kTI = wtemp
		END IF
		close(lunit)
	END IF

	return
	end subroutine read_input
