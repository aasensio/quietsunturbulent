module atomic_mod
use types_mod, only : spectrum_type, atmosphere_type, config_type, atomic_transition_type
use vars_mod, only : PI, PE, PC, PK, PH, UMA, OPA, PHK, ionization_state, EV_ERG, SQRTPI, conf, LARMOR
use maths_mod, only : extremes, saha, vecfvoigt, planck_vect_lambda, planck_lambda, formal_sol_vect, formal_sol_polarized,&
	strength_zeeman, vecfvoigt_zeeman, fvoigt
use background_opacity_mod, only : background_opacity
implicit none

contains

!-----------------------------------------------------------------
! Calculate the damping constant a for the Voigt profile
! following Gray 1976 "The observation and analysis of stellar photospheres"
!  Natural damping: eq (11.19)
!  van der Waals damping: eq (11.35) & (11.36)
!  Stark damping: eq (11.35) & (11.36)
!  Damping: eq (11.69)
! ion_pot : ionization potential in eV
! exc_pot_lower : excitation potential of the lower level in eV
! T : temperature array
! Pg : gas pressure array
! lambda : wavelength of the transition in A
! doppler_width : Doppler width in cm/s
!-----------------------------------------------------------------
	subroutine calculate_damping(T, nHtot, Pe, Pg, doppler_width, ion_pot, exc_pot_lower, gamma_rad, gamma_stark, gamma_vdw, lambda, &
		damping_out)
	real(kind=8) :: ion_pot, exc_pot_lower, T(:), nHtot(:), Pe(:), Pg(:), lambda, doppler_width(:)
	real(kind=8) :: damping_out(size(T)), gamma_rad, gamma_stark, gamma_vdw
	integer :: n, trans
	real(kind=8) :: c6, chi_lambda, cte, a, b
	real(kind=8), allocatable :: gamma6(:), gamma_nat(:), gamma_s(:)

		n = size(T)
		
		allocate(gamma6(n))
		allocate(gamma_s(n))
		allocate(gamma_nat(n))


! If the damping is not in the linelist, then use the formulas. In the other case, use those
! given by Kurucz

! RADIATIVE DAMPING
! Radiative damping
! 			Reference: Gray 1976, "The observation and analysis of stellar
!  		photospheres", pag 227, just after Eq. (11-19). This approximation
!  		is poor, but the radiative damping is usually negligible anyway.		

		gamma_nat = 0.22233d0 / (lambda*1.d-8)**2	

! STARK DAMPING
! Stark damping
! Formula used: gamma_4=38.8*(v**(1./3.))*(C4**(2./3.))*N , from
!   Unsold 1955, "Physik der Sternatmospharen", 2nd ed.,
!   Springer-Verlag, pp. 326 and 331. According to Gray (ref above),
!   pag 237, this is similar to
!   log (gamma_4) = 19.4 + 2./3*log C4 + log Pe - 5./6.*log T
!   The value of log C4 is approximated by -14. (see Gray, Table 11-3)

		gamma_s = 19.4 + 2.d0/3.d0*(-14.d0) + dlog10(Pe) - 5.d0/6.d0*dlog10(T)
		gamma_s = 10.d0**(gamma_s)

! HYDROGEN DAMPING
! Van der Waals damping:
! Reference: Gray (see above, pag 239, Eqs. (11-35) and (11-36))
! Formula used:
!    log (gamma_vdw) = 19.6 + 2./5.*log (C6(H)) + log Pg - 7./10.*log T

		chi_lambda = (PH*PC/EV_ERG) / (lambda*1.d-8)
		a = ion_pot-exc_pot_lower-chi_lambda
		b = ion_pot-exc_pot_lower
		c6 = 0.3d-30 * (1.d0/a**2 - 1.d0/b**2)
! It seems that the ionization potential and the energy of some levels are inconsistent, and E_exc>I
		if (c6 > 0.d0) then
			gamma6 = 10.d0**(19.6d0 + 0.4d0*dlog10(c6) + dlog10(Pg) - 0.7*dlog10(T))
		else
			gamma6 = 0.d0
		endif

		cte = 1.d-8 / (4.d0 * PI * PC)
		damping_out = cte * (gamma6+gamma_nat+gamma_s) * lambda**2 / (doppler_width * lambda / PC)

		deallocate(gamma6)
		deallocate(gamma_nat)
		deallocate(gamma_s)		

	end subroutine calculate_damping

!-----------------------------------------------------------------
! Fill the strength and splitting variables of the atomic_transition structure
!-----------------------------------------------------------------
	subroutine generate_atomic_zeeman_components(at)
	type(atomic_transition_type) :: at
	integer :: i, n, nlow, nup, iup, ilow, i_pi, i_red, i_blue, cual
	real(kind=8) :: Mup, Mlow, strength


		nup = 2*at%J_up+1
		nlow = 2*at%J_low+1

		i_pi = 0
		i_blue = 0
		i_red = 0

! First count the number of components of each type (
		at%n_components = 0
		do iup = 1, nup
			Mup = at%J_up + 1 - iup  ! Mup=J...-J
			do ilow = 1, 3
				Mlow = Mup-2+ilow			! Mlow=Mup-1,Mup,Mup+1 (the allowed transitions)
				if (abs(Mlow) <= at%J_low) then
					at%n_components(ilow) = at%n_components(ilow)+1
				endif
			enddo
		enddo

		allocate(at%splitting(3,maxval(at%n_components)))
		allocate(at%strength(3,maxval(at%n_components)))

! Now generate all the data
		do iup = 1, nup
			Mup = at%J_up + 1 - iup  ! Mup=J...-J
			do ilow = 1, 3
				Mlow = Mup-2+ilow			! Mlow=Mup-1,Mup,Mup+1 (the allowed transitions)
				if (abs(Mlow) <= at%J_low) then
					if (ilow == 1) then
						i_blue = i_blue + 1
						cual = i_blue
					endif
					if (ilow == 2) then
						i_pi = i_pi + 1
						cual = i_pi
					endif
					if (ilow == 3) then
						i_red = i_red + 1
						cual = i_red
					endif

					at%strength(ilow,cual) = strength_zeeman(at%J_up,at%J_low,Mup,Mlow)
					at%splitting(ilow,cual) = (at%lande_up*Mup - at%lande_low*Mlow)

				endif
			enddo
		enddo

	end subroutine generate_atomic_zeeman_components

	

!-----------------------------------------------------------------
! Return the profiles weighted by the strength of the components for a given frequency
! It returns zeeman_profile(q,n_depths) with
!  q=1  Mlow=Mup-1  (sigma blue)
!  q=2  Mlow=Mup    (sigma pi)
!  q=3  Mlow=Mup+1  (sigma red)
!-----------------------------------------------------------------
	subroutine zeeman_profile(damping,freq,at,vmacros,B_field,dnu,zeeman_voigt,zeeman_faraday)
	real(kind=8) :: freq
	type(atomic_transition_type) :: at
	real(kind=8) :: damping(:), vmacros(:), dnu(:), B_field(:)
	real(kind=8) :: zeeman_voigt(3,size(damping)), zeeman_faraday(3,size(damping))
	integer :: n, nlow, nup, iup, ilow, i_pi, i_blue, i_red, cual
	real(kind=8) :: Mup, Mlow, strength
	real(kind=8), allocatable :: splitting(:), profile(:,:)

		n = size(damping)

		allocate(splitting(n))
		allocate(profile(2,n))

		nup = 2*at%J_up+1
		nlow = 2*at%J_low+1

		zeeman_voigt = 0.d0
		zeeman_faraday = 0.d0

		i_red = 0
		i_pi = 0
		i_blue = 0

		do iup = 1, nup
			Mup = at%J_up + 1 - iup  ! Mup=J...-J
			do ilow = 1, 3
				Mlow = Mup-2+ilow			! Mlow=Mup-1,Mup,Mup+1 (the allowed transitions)
				if (abs(Mlow) <= at%J_low) then

					if (ilow == 1) then
						i_blue = i_blue + 1
						cual = i_blue
					endif
					if (ilow == 2) then
						i_pi = i_pi + 1
						cual = i_pi
					endif
					if (ilow == 3) then
						i_red = i_red + 1
						cual = i_red
					endif

					strength = at%strength(ilow,cual)

					splitting = LARMOR * B_field * at%splitting(ilow,cual)

					profile = vecfvoigt_zeeman(damping,(at%freq - freq - at%freq*vmacros / PC + splitting ) / dnu)

					zeeman_voigt(ilow,:) = zeeman_voigt(ilow,:) + strength * profile(1,:)
					zeeman_faraday(ilow,:) = zeeman_faraday(ilow,:) + strength * profile(2,:)
				endif
			enddo
		enddo

		deallocate(splitting)
		deallocate(profile)

	end subroutine zeeman_profile

!-----------------------------------------------------------------
! Add line opacity to the total opacity
!-----------------------------------------------------------------
	subroutine add_line_opacity(atm, spectrum, mu)
	type(config_type) :: conf
	type(atmosphere_type) :: atm
	type(spectrum_type) :: spectrum
	real(kind=8) :: mu, dnumax

	real(kind=8), allocatable :: u1(:), u2(:), u3(:), n1overn0(:), n2overn1(:), niovern(:), at_opacity(:)
	real(kind=8), allocatable :: dnu(:), doppler_atomic(:), perfil(:), ui(:), damping(:)
	integer, allocatable :: from_to(:,:)
	real(kind=8), allocatable :: zeeman_voigt(:,:), zeeman_faraday(:,:), coefficients(:,:)
	integer :: i, j, k, from_to_line(2)


		allocate(u1(atm%n_depths))
		allocate(u2(atm%n_depths))
		allocate(u3(atm%n_depths))
		allocate(ui(atm%n_depths))
		allocate(n1overn0(atm%n_depths))
		allocate(n2overn1(atm%n_depths))
		allocate(niovern(atm%n_depths))
		allocate(at_opacity(atm%n_depths))
		allocate(damping(atm%n_depths))
		allocate(dnu(atm%n_depths))
		allocate(doppler_atomic(atm%n_depths))
		allocate(perfil(atm%n_depths))

! If magnetic, allocate memory for the Voigt and Faraday profiles
		if (index(atm%magnetic, 'NONMAGNETIC') == 0) then
			allocate(zeeman_voigt(3,atm%n_depths))
			allocate(zeeman_faraday(3,atm%n_depths))
			allocate(coefficients(7,atm%n_depths))		
			zeeman_voigt = 0.d0
			zeeman_faraday = 0.d0

! Additionally, compute the angles between the magnetic field and the LOS
! for the inclination angle of the LOS and the angles of B wrt to the local vertical
			atm%cos_gamma = (1.d0-mu**2) * sin(atm%thetaB) * cos(atm%chiB) + mu * cos(atm%thetaB)
			atm%sin_gamma = sqrt(1.d0 - atm%cos_gamma**2)
			atm%cos_2chi = cos(2.d0*atm%chiB)
			atm%sin_2chi = sin(2.d0*atm%chiB)
		endif

! Calculate the line opacity
		do i = 1, spectrum%nlines

! Compute line integrated opacity
! If atomic line			
			if (spectrum%line(i)%atomic_molecular == 1) then
				u1 = atm%partition_functions(:, spectrum%line(i)%element, 1)
				u2 = atm%partition_functions(:, spectrum%line(i)%element, 2)
				u3 = atm%partition_functions(:, spectrum%line(i)%element, 3)

				n1overn0 = saha(atm%T, PK * atm%ne * atm%T, u1, u2, atm%ionization_pot(1,spectrum%line(i)%element))
				n2overn1 = saha(atm%T, PK * atm%ne * atm%T, u2, u3, atm%ionization_pot(2,spectrum%line(i)%element))

				if (spectrum%line(i)%ionization == 0) then
					niovern = 1.d0 / (1.d0 + n1overn0 + n2overn1 * n1overn0)
					ui = u1
				else if (spectrum%line(i)%ionization == 1) then
					niovern = 1.d0 / (1.d0 + 1.d0/n1overn0 + n2overn1)
					ui = u2
				else
					print *, 'Unknown ionization stage...'
					niovern = 0.d0
				endif

! Compute the line center opacity
				at_opacity = OPA * spectrum%line(i)%gf / ui * dexp(-spectrum%line(i)%expot / (PK * atm%T)) *&
									(1.d0 - dexp(-PHK * spectrum%line(i)%freq / atm%T)) * niovern * &
									(atm%nhtot * atm%abundances(spectrum%line(i)%element))

! Doppler width in velocity and frequency units
				doppler_atomic = dsqrt(atm%microturb**2.d0 + 2.d0 * PK * atm%T / (spectrum%line(i)%mass * UMA))
				dnu = spectrum%line(i)%freq * doppler_atomic / PC

! Compute damping
				call calculate_damping(atm%T, atm%nHtot, atm%Pe, atm%P_total, doppler_atomic, &
					atm%ionization_pot(spectrum%line(i)%ionization+1,spectrum%line(i)%element), spectrum%line(i)%expot/EV_ERG, &
					spectrum%line(i)%gamma_rad, spectrum%line(i)%gamma_stark, spectrum%line(i)%gamma_vdw, spectrum%line(i)%wlength, damping)

! Try to estimate the range of wavelengths in which the line has opacity contribution
				dnumax = estimate_range_line(spectrum%line(i)%wlength, dnu, damping, at_opacity, atm%height) * maxval(dnu)

! If the line is estimated to be very broad, take into account that the result
! might be truncated if the incorrect chunks are used
				from_to_line = extremes(spectrum%freq, spectrum%line(i)%freq, dnumax)

				spectrum%line(i)%lambda_from = from_to_line(1)
				spectrum%line(i)%lambda_to = from_to_line(2)
					
			else

! If molecular line				
				at_opacity = OPA * spectrum%line(i)%gf / atm%partition_functions_molecular(:,spectrum%line(i)%molecule) * &
									dexp(-spectrum%line(i)%expot / (PK * atm%T)) *&
									(1.d0 - dexp(-PHK * spectrum%line(i)%freq / atm%T)) * &
									atm%mol_density(spectrum%line(i)%molecule,:)

				doppler_atomic = dsqrt(atm%microturb**2.d0 + 2.d0 * PK * atm%T / (spectrum%line(i)%mass * UMA))
				dnu = spectrum%line(i)%freq * doppler_atomic / PC

				damping = 0.d0
									
			endif
 


! Compute opacity
			if (index(atm%magnetic, 'NONMAGNETIC') /= 0) then
			
! Non-magnetic			
				do k = spectrum%line(i)%lambda_from, spectrum%line(i)%lambda_to
					perfil = vecfvoigt(damping,(spectrum%line(i)%freq - spectrum%freq(k) - spectrum%line(i)%freq * atm%vmacro / PC) / dnu)
					spectrum%chi(1,:,k) = spectrum%chi(1,:,k) + at_opacity * perfil / (dnu * SQRTPI)
				enddo
				
			else
			
! Magnetic				
				do k = spectrum%line(i)%lambda_from, spectrum%line(i)%lambda_to

					call zeeman_profile(damping, spectrum%freq(k), spectrum%line(i), atm%vmacro, atm%B, dnu, zeeman_voigt, zeeman_faraday)
					call zeeman_opacity(zeeman_voigt, zeeman_faraday, atm, coefficients)
					
					do j = 1, 7
						spectrum%chi(j,:,k) = spectrum%chi(j,:,k) + at_opacity * coefficients(j,:) / (dnu * SQRTPI)
					enddo
					
				enddo
				
			endif

		enddo

		deallocate(u1)
		deallocate(u2)
		deallocate(u3)
		deallocate(ui)
		deallocate(n1overn0)
		deallocate(n2overn1)
		deallocate(niovern)
		deallocate(at_opacity)
		deallocate(dnu)		
		deallocate(doppler_atomic)
		deallocate(perfil)
		deallocate(damping)

		if (index(atm%magnetic, 'NONMAGNETIC') == 0) then
			deallocate(zeeman_voigt)
			deallocate(zeeman_faraday)
			deallocate(coefficients)
		endif
		
	end subroutine add_line_opacity

!-----------------------------------------------------------------
! Add line opacity to the total opacity
!-----------------------------------------------------------------
	subroutine add_background_opacity(atm, spectrum)
	integer :: i, j
	type(atmosphere_type) :: atm
	type(spectrum_type) :: spectrum

		spectrum%chi = 0.d0
		
		do i = 1, atm%n_depths
			do j = 1, spectrum%nlambda				
				spectrum%chi(1,i,j) = background_opacity(atm%T(i), atm%Pe(i), atm%PH(i), atm%PHminus(i), &
					atm%PHplus(i), atm%PH2(i), atm%PH2plus(i), PC / spectrum%freq(j)*1.d8)
			enddo
		enddo

	end subroutine add_background_opacity

!-----------------------------------------------------------------
! Add line opacity to the total opacity
!-----------------------------------------------------------------
	subroutine solve_RT(atm, spectrum, mu, weight)
	integer :: i
	type(atmosphere_type) :: atm
	type(spectrum_type) :: spectrum
	real(kind=8) :: mu, weight

		spectrum%boundary = 0.d0

! Source function and boundary condition		
		do i = 1, spectrum%nlambda
			spectrum%source(:,i) = planck_vect_lambda(PC/spectrum%freq(i),atm%T)
			spectrum%boundary(1,i) = planck_lambda(PC / spectrum%freq(i), atm%T(atm%n_depths))
		enddo

		if (index(atm%magnetic, 'NONMAGNETIC') /= 0) then
! Non-magnetic			
			spectrum%flux(1,:) = spectrum%flux(1,:) + &
				weight * formal_sol_vect(atm%height, spectrum%chi(1,:,:), spectrum%source(:,:), mu, spectrum%boundary(1,:), atm%direction)
			

		else
! Magnetic			
			do i = 1, spectrum%nlambda
				spectrum%flux(:,i) = spectrum%flux(:,i) + &
					weight * formal_sol_polarized(atm%height, spectrum%chi(:,:,i), spectrum%source(:,i), mu, spectrum%boundary(:,i), atm%direction)
	
			enddo

		endif
						
	end subroutine solve_RT

!-----------------------------------------------------------------
! Generate wavelength axis
!-----------------------------------------------------------------
	subroutine setup_spectrum(conf, atm, spectrum)
	type(atmosphere_type) :: atm
	type(spectrum_type) :: spectrum
	type(config_type) :: conf
	integer :: i

		spectrum%nlambda = (spectrum%end_lambda - spectrum%ini_lambda) / conf%step_lambda
		allocate(spectrum%source(atm%n_depths, spectrum%nlambda))

		allocate(spectrum%freq(spectrum%nlambda))

		if (index(atm%magnetic, 'NONMAGNETIC') /= 0) then
! Non-magnetic			
			allocate(spectrum%flux(1,spectrum%nlambda))
			allocate(spectrum%chi(1,atm%n_depths, spectrum%nlambda))			
			allocate(spectrum%boundary(1,spectrum%nlambda))

		else
! Magnetic
			allocate(spectrum%flux(4,spectrum%nlambda))
			allocate(spectrum%chi(7,atm%n_depths, spectrum%nlambda))
			allocate(spectrum%boundary(4,spectrum%nlambda))

		endif		

! Generate the frequency axis
 		do i = 1, spectrum%nlambda
 			spectrum%freq(i) = PC / (1.d-8 * (spectrum%ini_lambda + (i-1.d0) * conf%step_lambda))
 		enddo

 	end subroutine setup_spectrum

end module atomic_mod