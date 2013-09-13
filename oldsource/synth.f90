module synthModule
use globalModule, only : OPA, PH, PC, PK, PHK, UMA, EV_ERG, SQRTPI, atmosphereType, lineListType
use atomicPartitionModule, only : partitionAtomic
use mathsModule, only : saha, vecfvoigt, calculateDamping, planckFrequency, shortCharacteristics, computeHeight, gasPressure, vecfvoigt2
use backgroundOpacityModule, only : backgroundOpacity
use hydrostaticModule, only : hydrostaticEquilibrium
implicit none
contains

!------------------------------------------------
! Synthesize all lines
!------------------------------------------------
	subroutine synthLines(atmosphere,lineList)
	type(atmosphereType) :: atmosphere
	type(lineListType) :: lineList
	real(kind=8) :: n1overn0, n2overn1, ei1, ei2, weight, abundance
	real(kind=8), allocatable :: voigtProfile(:)
	integer :: i, j, loop, fromLambda, toLambda
		
		if (.not.allocated(voigtProfile)) allocate(voigtProfile(atmosphere%nDepths))
	
! Since all lines are Fe I, precompute the partition functions
		call partitionAtomic(atmosphere%T, 26, atmosphere%u1, atmosphere%u2, atmosphere%u3, ei1, ei2, weight, abundance)

! Put the atmosphere in hydrostatic equilibrium to compute the electron pressure		
 		call hydrostaticEquilibrium(atmosphere) 	 		
 		
! Compute the chemical equilibrium
		call gasPressure(atmosphere%Pe, atmosphere%T, atmosphere%PH, atmosphere%PHminus, atmosphere%PHplus, atmosphere%PH2, atmosphere%PH2plus, atmosphere%PTotal)
				
		atmosphere%n1overn0 = saha(atmosphere%T, atmosphere%Pe, atmosphere%u1, atmosphere%u2, ei1)
		atmosphere%n2overn1 = saha(atmosphere%T, atmosphere%Pe, atmosphere%u2, atmosphere%u3, ei2)
		atmosphere%niovern = 1.d0 / (1.d0 + atmosphere%n1overn0 + atmosphere%n2overn1 * atmosphere%n1overn0)
		atmosphere%ui = atmosphere%u1		 		

! Total hydrogen density
		atmosphere%nhtot = (atmosphere%PH + atmosphere%PHminus + atmosphere%PHplus + atmosphere%PH2 + atmosphere%PH2plus) / (PK * atmosphere%T)
		
! Generate height axis
		do j = 1, atmosphere%nDepths
			atmosphere%opacity500(j) = backgroundOpacity(atmosphere%T(j), atmosphere%Pe(j), atmosphere%PH(j), atmosphere%PHminus(j), &
				atmosphere%PHplus(j), atmosphere%PH2(j), atmosphere%PH2plus(j), 5000.d0)		
		enddo
		
! Compute the height depth scale
		call computeHeight(atmosphere%lTau500, atmosphere%opacity500, atmosphere%height)			
		
		fromLambda = 1
		loop = 1
		
		do i = 1, lineList%nLines			
			if (lineList%transition(i)%active) then
			
! Compute line opacity
				lineList%transition(i)%lineOpacity = OPA * lineList%transition(i)%gf / atmosphere%ui * dexp(-lineList%transition(i)%Elow / (PK * atmosphere%T)) *&
					(1.d0 - dexp(-PHK * lineList%transition(i)%frequency0 / atmosphere%T)) * atmosphere%niovern * &
					(atmosphere%nhtot * atmosphere%abundance)													
					
! Doppler width
				lineList%transition(i)%dopplerWidth = dsqrt(atmosphere%microturbulence**2.d0 + 2.d0 * PK * atmosphere%T / (weight * UMA))
				lineList%transition(i)%deltaNu = lineList%transition(i)%dopplerWidth * lineList%transition(i)%frequency0 / PC
				
! Background opacity
				do j = 1, atmosphere%nDepths
					lineList%transition(i)%backOpacity(j) = backgroundOpacity(atmosphere%T(j), atmosphere%Pe(j), atmosphere%PH(j), atmosphere%PHminus(j), &
						atmosphere%PHplus(j), atmosphere%PH2(j), atmosphere%PH2plus(j), lineList%transition(i)%lambdaNew)
				enddo
				
! Compute the damping
				lineList%transition(i)%damping = calculateDamping(atmosphere%T, atmosphere%nHtot, atmosphere%Pe, atmosphere%Pg, &
					lineList%transition(i)%dopplerWidth, ei1, lineList%transition(i)%Elow / EV_ERG, 0.d0, 0.d0, 0.d0, lineList%transition(i)%lambdaNew, &
					lineList%transition(i)%alphaABO, lineList%transition(i)%sigmaABO, weight)
									
! Compute frequency dependency of line opacity, source function and boundary condition
				do j = 1, lineList%transition(i)%nLambda
 					voigtProfile = vecfvoigt2(lineList%transition(i)%damping,(lineList%transition(i)%frequency0 - lineList%transition(i)%frequency(j)) / lineList%transition(i)%deltaNu)
					
					lineList%transition(i)%opacity(:,j) = lineList%transition(i)%lineOpacity * voigtProfile / (lineList%transition(i)%deltaNu * SQRTPI) + &
						lineList%transition(i)%backOpacity
						
					lineList%transition(i)%opacityContinuum(:,j) = lineList%transition(i)%backOpacity					
					
					lineList%transition(i)%source(:,j) = planckFrequency(lineList%transition(i)%frequency(j), atmosphere%T)
					
					lineList%transition(i)%boundary = lineList%transition(i)%source(atmosphere%nDepths,j)
				enddo
				
				toLambda = fromLambda + lineList%transition(i)%nLambda
				
				loop = loop + lineList%transition(i)%nLambda
								
				lineList%opacity(:,fromLambda:toLambda) = lineList%transition(i)%opacity
				lineList%opacityContinuum(:,fromLambda:toLambda) = lineList%transition(i)%opacityContinuum
				lineList%source(:,fromLambda:toLambda) = lineList%transition(i)%source
				lineList%boundary(fromLambda:toLambda) = lineList%transition(i)%boundary
				
				fromLambda = toLambda
								
				
!  				lineList%transition(i)%intensity = shortCharacteristics(atmosphere%height, lineList%transition(i)%opacity, lineList%transition(i)%source, 1.d0, lineList%transition(i)%boundary, -1) 				
!  				lineList%transition(i)%continuumIntensity = shortCharacteristics(atmosphere%height, lineList%transition(i)%opacityContinuum, lineList%transition(i)%source, 1.d0, lineList%transition(i)%boundary, -1)
 				 				
!  				do j = 1, lineList%transition(i)%nLambda
! 					write(12,*) lineList%transition(i)%lambda(j), lineList%transition(i)%intensity(j) / lineList%transition(i)%continuumIntensity(j)
! 				enddo				
								
			endif			
						
		enddo
				
		lineList%intensity = shortCharacteristics(atmosphere%height, lineList%opacity, lineList%source, 1.d0, lineList%boundary, -1) 				
 		lineList%intensityContinuum = shortCharacteristics(atmosphere%height, lineList%opacityContinuum, lineList%source, 1.d0, lineList%boundary, -1)
		
! 		loop = 1
! 		do i = 1, lineList%nLines			
! 			if (lineList%transition(i)%active) then
! 				do j = 1, lineList%transition(i)%nLambda
! 					write(12,*) lineList%transition(i)%lambda(j), lineList%intensity(loop) / lineList%intensityContinuum(loop)
! 					loop = loop + 1
! 				enddo
! 			endif			
! 		enddo				
		
	end subroutine synthLines
end module synthModule