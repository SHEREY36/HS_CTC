
	subroutine read_input
	
	use particles
	use run_param
	use output
	
	implicit none
	
	integer :: I
	integer :: lunit = 10
	character(len=255) :: filename_input
	double precision :: utemp, wtemp
	
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
	return
	end subroutine read_input
