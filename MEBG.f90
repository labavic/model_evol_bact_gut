module glob

	
	implicit none
	
	
	
	! global variables
	integer		:: tI, lI, tw, mut_add
	real(8)		:: a, r, k, d, v, Fin, rm 		! model parameters
	real(8)		:: dx, dt, length, time, wt		! algorithm parameters
	real(8)		:: fixed, Bprop, Fc, M0
	character	:: path*200, chD*10, chV*10, chr*10, chk*10, chf*10, chL*10
	logical		:: loglab
	
	
	contains !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	
		subroutine calc_steady_state()
			
			implicit none
			
			real(8), allocatable			:: F(:), B(:), Ft(:), Bt(:)
			integer							:: it, ix, i, j
			
				
			! allocating dimensions of concentration arrays 
			allocate(F(0:lI), B(0:lI), Ft(0:lI), Bt(0:lI))
			
				
			! initial conditions:
			F = 0.9 * Fin; B = a*(Fin-F);
			
			
					
			
			
			do it = 1, tI  ! integration in time:
				
				! boundary at x = 0
				Ft(0) = F(0) + 2.*d*dt/dx**2 * ( F(1) - F(0)*(1.+dx*v/d) + dx*v/d*Fin ) - &
					  & v**2*dt/d * (F(0)-Fin) - r*dt/a*(F(0)/(k + F(0)	)) * B(0)
				
				Bt(0) = B(0) + 2.*d*dt/dx**2 * ( B(1) - B(0)*(1.+dx*v/d)) - &
					  & v**2*dt/d * (B(0)) + r*dt*(F(0)/(k + F(0))) * B(0)
					  
				! middle of the lattice
				do ix = 1, lI - 1
				
					Ft(ix) = F(ix) + d*dt/(dx)**2 * ( F(ix+1) - 2.*F(ix) + F(ix-1) ) - &
						  & v*dt/dx * (  F(ix+1) - F(ix-1) )/2. - r*dt/a * (F(ix)/(k + F(ix))) * B(ix)
							
					Bt(ix) = B(ix) + d*dt/(dx)**2 * ( B(ix+1) - 2.*B(ix) + B(ix-1) ) - &
						  & v*dt/dx * (  B(ix+1) - B(ix-1) )/2. + r*dt *( F(ix)/(k + F(ix))) * B(ix)
				
				end do
				
				! boundary at x = L
				Ft(lI) = F(lI) + 2.*d*dt/dx**2 * ( F(lI-1) - F(lI) ) - r*dt/a*(F(lI)/(k + F(lI))) * B(lI)
				Bt(lI) = B(lI) + 2.*d*dt/dx**2 * ( B(lI-1) - B(lI) ) + r*dt*(F(lI)/(k + F(lI))) * B(lI)
				
				
				! replace old values with new ones
				F = Ft
				B = a*(Fin-F)
			
				
			end do
			

			
			
			! write the final state to save it for initial conditions for adding mutant later
	
			open(unit=1, file=trim(path)//'/FinalSS.dat', form='formatted', status='unknown')
			do i = 0, lI
				write(1, *) B(i), F(i)
			end do
			close(1)
			
			deallocate(F, B, Ft, Bt)
			
			100 return
			
			
			
			
		end subroutine
		
		!***************************************************************
		!***************************************************************
		!***************************************************************
		!***************************************************************
		
		
		
		
		
		subroutine calc_add_mut( )
		
		
		
			implicit none
			real(8), allocatable	:: F(:), B(:), Ft(:), Bt(:), M(:), Mt(:), B0(:), F0(:)
			real(8)					:: M0
			integer					:: it, ix, i, P, fln, j
			
			
			
			
			
			allocate(F(0:lI), B(0:lI), Ft(0:lI), Bt(0:lI), M(0:lI), Mt(0:lI), B0(0:lI), F0(0:lI))
			
			
			
			
			
			! read steady state of B and F from FinalSS.dat and store it to F0, B0 for
			! initial conditions when adding mutant
						
				
			open(unit=1, file=trim(path)//'/FinalSS.dat',  form='formatted', status='old')
            do i = 0, lI
				read(1, *) B0(i), F0(i)
			end do
			close(1)
			
			
			
			mut_add = floor(.2d0/dx) 
			
			
			open(unit=2, file=trim(path)//'/dataMutantEnd.dat', form='formatted')
			
			do i = 0, lI ! different positions for a mutant seeding
			
				if (mod(i, mut_add) == 0) then
					
					P = i ! seed position
										
					! initial condions
					B = B0; F = F0; M = 0.d0
					M(P) = M0/dx *( fixed + Bprop *r*F(P)/(k+F(P)) )
					
					if (P == 0 .or. P == lI) M(P) = M(P)*2. ! needed due to the lattice discretization choice
					
					
					
					do it = 1, tI ! integration in time:
						
						! boundary conditions at x = 0
						Ft(0) = F(0) + 2.*d*dt/dx**2 * ( F(1) - F(0)*(1.+dx*v/d) + dx*v/d*Fin ) - &
							  & v**2*dt/d * (F(0)-Fin) - dt/a*F(0)/(k + F(0)) * (r*B(0)+r*M(0))	
						
						Bt(0) = B(0) + 2.*d*dt/dx**2 * ( B(1) - B(0)*(1.+dx*v/d)) - &
							  & v**2*dt/d * (B(0)) + r*dt*F(0)/(k + F(0)) * B(0)
							  
						Mt(0) = M(0) + 2.*d*dt/dx**2 * ( M(1) - M(0)*(1.+dx*v/d)) - &
							  & v**2*dt/d * (M(0)) + rm*dt*F(0)/(k + F(0)) * M(0)
							  
						! middle of the lattice	  
						do ix = 1, lI - 1
						
							Ft(ix) = F(ix) + d*dt/(dx)**2 * ( F(ix+1) - 2.*F(ix) + F(ix-1) ) - &
								  & v*dt/dx * (  F(ix+1) - F(ix-1) )/2. - dt/a * F(ix)/(k + F(ix)) * (r*B(ix)+r*M(ix))
									
							Bt(ix) = B(ix) + d*dt/(dx)**2 * ( B(ix+1) - 2.*B(ix) + B(ix-1) ) - &
								  & v*dt/dx * (  B(ix+1) - B(ix-1) )/2. + r*dt * F(ix)/(k + F(ix)) * B(ix)
								  
							Mt(ix) = M(ix) + d*dt/(dx)**2 * ( M(ix+1) - 2.*M(ix) + M(ix-1) ) - &
								  & v*dt/dx * (  M(ix+1) - M(ix-1) )/2. + rm*dt * F(ix)/(k + F(ix)) * M(ix)
						
						end do
						
						! boundary conditions at x = L
						Ft(lI) = F(lI) + 2.*d*dt/dx**2 * ( F(lI-1) - F(lI) ) - dt/a*F(lI)/(k + F(lI)) * (r*B(lI)+r*M(lI))	
						Bt(lI) = B(lI) + 2.*d*dt/dx**2 * ( B(lI-1) - B(lI) ) + r*dt*F(lI)/(k + F(lI)) * B(lI)
						Mt(lI) = M(lI) + 2.*d*dt/dx**2 * ( M(lI-1) - M(lI) ) + rm*dt*F(lI)/(k + F(lI)) * M(lI)
						
						! replace old values with new ones
						F = Ft
						B = Bt
						M = Mt
						
						
					
					end do
					
					
					write(2, *)P*dx, B, M ,F
					
					
				end if
				
			end do
			
			!***********************************************************
			
			close(2)
		
			deallocate(F,B, Ft, Bt, M, Mt, B0, F0)
			
			
		
		
		end subroutine calc_add_mut
		
		
		!***************************************************************
		!***************************************************************
		!***************************************************************
		!***************************************************************
		
		
		subroutine var_char(var, ch_var)
		
			real(8), intent(in)				:: var
			character(len=10), intent(out)	:: ch_var
			character						:: ch_num*2, ft*10
			integer							:: i
			
			
			if (var < 1.d0) then
			write(ch_var, '(f10.8)') var
				else
						do i = 1, 9
			
								if (var<dble(10**i)) then
										write(ch_num, '(i2.2)') 9-i
										ft = '(f10.'//ch_num//')'
										write(ch_var, ft) var
										exit
								end if
					   end do
			
				end if
			
		
		end subroutine


end module

!***********************************************************************
!***********************************************************************
!***********************************************************************
!***********************************************************************
!***********************************************************************



program main
	
	
	use glob
	
	implicit none
	
	! define type
	
	type vars
		real(8)				:: vlu
		character(len = 5)	:: nme
	end type
	
	integer		:: i, j, cac, mi, process, ierror, P, proc_tot
	real(8)		:: miR 
	!logical	::
	character	:: fname*100, arg*30, Ch_mi*20
	type(vars)	:: var(1:16)
	
	
	! default values if no input from a command line
	
	d 		= 0.2				! cm^2/h 
	v 		= 0.5		  		! cm/h
	a		= 0.184				! OD/mM (OD is optical density)
	k		= 0.10				! mM
	r		= 0.42				! 1/h 
	rm		= 0.42				! 1/h
	dx	 	= 0.01				! cm
	dt		= dx**2*0.5/d*0.8	! h
	time	= 1.0				! h 
	length	= 6.0				! cm
	Fc		= 1.0
	Fin		= Fc/v				! mM
	
	wt		= 1000000.			! write state every wt integration steps (1/dt), not used here
	P		= 0.				! mutant seed position
	miR		= 0.				! mutant initial concentration, if 0. fixed, if 1. proportional to repr. rate
	M0		= 1.d-20			! initial mutant amount 
	
	var%nme = (/ "d    " , "v    ", "a    ", "k    ", "r    ", "rm   ", "Fin  ", "dx   " &
				& , "dt   ", "time ", "len  ", "wt   ", "P    ", "mi   " , "Fc   ", "M0   "  /)
	var%vlu = (/ d, v, a, k, r, rm, Fin, dx, dt, time, length, wt, 0.d0, 0.d0, Fc, M0 /)
	
	! read the variables from a command line
	
	cac = command_argument_count()
	
	if (cac == 0) goto 100 ! if no input from the command line take default values
	
	! read the parameter values from the command line
	do i = 1, cac-1, 2
		call getarg(i, arg)
		do j = 1, 16
			if (arg == trim(var(j)%nme)) then
				call getarg(i+1,arg)
				read(arg,*) var(j)%vlu
				exit
			end if
		end do
	end do
	
	
	100 continue
	
	d 		= var(1)%vlu	; v 		= var(2)%vlu		
	a		= var(3)%vlu	; k			= var(4)%vlu			
	r		= var(5)%vlu	; rm		= var(6)%vlu
	Fin		= var(7)%vlu	; dx	 	= var(8)%vlu	
	dt		= var(9)%vlu	; time		= var(10)%vlu	
	length	= var(11)%vlu	; wt		= var(12)%vlu
	P		= floor(var(13)%vlu); miR	= var(14)%vlu
	Fc		= var(15)%vlu	;M0			= var(16)%vlu
	
	mi = Ceiling(miR)
	
	
	
	! write all parameters of the program in a param file and save in the same
	! directory as integration results
	
	call var_char(v, chV)
	call var_char(D, chD)
	call var_char(length, chL)
	call var_char(Fc,chf)
	call var_char(r, chr)
	call var_char(k, chk)

	
	
	! create the directory to save the data
	if (mi == 0) then
		Ch_mi="Fixed"
		fixed = 1.d0
		Bprop = 0.d0
	elseif (mi == 1) then
		Ch_mi="B"
		fixed = 0.d0
		Bprop = 1.d0
	else
		print*, "something not right"
		goto 500
	end if
	
	path = '/media/darka/data/Data/'//trim(Ch_mi)//&
					&'/k_'//trim(adjustl(chk))//&
					&'/r_'//trim(adjustl(chr))//&
					&'/F_'//trim(adjustl(chf))//&
					&'/L_'//trim(adjustl(chL))//&
					&'/D_'//trim(adjustl(chD))//&
					&'/v_'//trim(adjustl(chV))

	call system('mkdir -p '//path)
	
	
	
	open(unit = 2, file = trim(path)//'/parameters.dat', status = 'unknown')
	! write parameters 
	write(2, *) 'diffusion [cm^2/h]            ', d
	write(2, *) 'velocity  [cm/h]              ', v
	write(2, *) 'growth rate Bac  [1/h]        ', r
	write(2, *) 'growth rate Mut  [1/h]        ', rm
	write(2, *) 'food inflow  [mM]             ', Fin
	write(2, *) 'Fin * v                       ', Fc
	write(2, *) 'scaling growth food [OD/mM]   ', a
	write(2, *) 'Monod constant   [mM]         ', k
	write(2, *) '  '
	write(2, *) 'space step  [mm]              ', dx
	write(2, *) 'time step   [h]               ', dt
	write(2, *) 'length      [cm]              ', length
	write(2, *) 'integration time   [h]        ', time
	write(2, *) 'initial mutant amount 		   ', M0
	close(2)

	
	
	
	
	tI  = min(abs(ceiling(time/dt)),2147483647)	! number of time iteration steps
	lI  = ceiling(length/dx)					! number of points on a 1D spatial grid
		
	
	
	
	
	
	!call subroutine to integrate the system without a mutant to a steady state
		
	call calc_steady_state( )
	
	call calc_add_mut( )
	
	500 print*,"end"
	
	
	
end program
