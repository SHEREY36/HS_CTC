	module rng_mod

	use, intrinsic :: iso_fortran_env, only: int64

	implicit none

	integer(int64), parameter :: RNG_MODULUS = 2147483647_int64
	integer(int64), parameter :: RNG_MULTIPLIER = 16807_int64
	integer(int64), parameter :: RNG_QUOTIENT = 127773_int64
	integer(int64), parameter :: RNG_REMAINDER = 2836_int64
	integer(int64), parameter :: RNG_BASE_SEED = 12345_int64
	integer(int64), parameter :: RNG_THREAD_STRIDE = 104729_int64

	integer(int64) :: rng_state = 1_int64
	logical :: rng_initialized = .false.

!$OMP THREADPRIVATE(rng_state, rng_initialized)

	contains

	subroutine init_thread_rng(thread_id)
		integer, intent(in) :: thread_id
		integer(int64) :: seed_value
		integer :: i

		seed_value = modulo(RNG_BASE_SEED + RNG_THREAD_STRIDE * int(thread_id, int64), &
			RNG_MODULUS - 1_int64) + 1_int64
		rng_state = seed_value

		! Warm up the stream so adjacent thread seeds diverge further.
		do i = 1, 8
			call advance_rng_state()
		end do

		rng_initialized = .true.
	end subroutine init_thread_rng

	double precision function rng_uniform()
		if (.not. rng_initialized) then
			write(*,*) 'Thread RNG used before initialization'
			stop 1
		end if

		call advance_rng_state()
		rng_uniform = dble(rng_state) / dble(RNG_MODULUS)
	end function rng_uniform

	subroutine advance_rng_state()
		integer(int64) :: hi, lo, test_value

		hi = rng_state / RNG_QUOTIENT
		lo = modulo(rng_state, RNG_QUOTIENT)
		test_value = RNG_MULTIPLIER * lo - RNG_REMAINDER * hi

		if (test_value > 0_int64) then
			rng_state = test_value
		else
			rng_state = test_value + RNG_MODULUS
		end if
	end subroutine advance_rng_state

	end module rng_mod
