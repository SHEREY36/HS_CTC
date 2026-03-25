
	subroutine measure_final()

	use output
	use run_param
	use particles
	use constants

	implicit none
	integer :: I,J

	! Flush any remaining buffered data
	CALL FLUSH_BUFFERS()

	write(1000,'(E14.8)') 4.D0*DBLE(NHIT)/DBLE(NTRY)*BMAX*BMAX
	FLUSH(1000)
	close(1000); close(1001); close(1002); close(1003); close(1004); close(1111)
	close(1005); close(1006); close(1007)
	RETURN

	end subroutine measure_final
