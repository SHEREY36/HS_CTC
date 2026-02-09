
	subroutine measure_final()
	
	use output
	use run_param
	use particles
	use constants
	
	implicit none
	integer :: I,J

	write(1000,'(E14.8)') DBLE(NHIT-1)/DBLE(NTRY)*BMAX*BMAX
	FLUSH(1000)
	close(1000); close(1001); close(1002); close(1003)
	RETURN
	
	end subroutine measure_final
