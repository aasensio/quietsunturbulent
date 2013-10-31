module ioModule
use globalModule, only : atlasType, lineListType, PH, PC, nVariables, priorLower, priorUpper
use mathsModule, only: extremes
use dualClass
implicit none
contains

!------------------------------------------------
! Read input data
!------------------------------------------------
	subroutine readInput(lineList, atlas)	
	type(lineListType) :: lineList
	type(atlasType) :: atlas
	integer :: i
	
		open(unit=12,file='turbHIER.conf',action='read',status='old')
		
		read(12,*)
		read(12,*) lineList%file
		
! Read the line list
		call readLineList(lineList)
		
		read(12,*)
		read(12,*)
		read(12,*) atlas%file
		
! Read the atlas
		call readAtlas(atlas)
				
! Set-up line information
 		call setupLines(lineList,atlas)
 		
! Read prior information
		allocate(priorLower(12))
		allocate(priorUpper(12))
		
		read(12,*)		
		read(12,*)
		do i = 1, 12
			read(12,*)
			read(12,*) priorLower(i), priorUpper(i)
			print *, priorLower(i), priorUpper(i)
		enddo
 		close(12)
		
	end subroutine readInput
		
!------------------------------------------------
! Read the model atmosphere
!------------------------------------------------
	subroutine readLineList(lineList)
	integer :: i, active, activeLines
	type(lineListType) :: lineList
	type(dual) :: test
	
		open(unit=13,file=trim(adjustl(lineList%file)),action='read',status='old')
		read(13,*) lineList%nLines
		
		allocate(lineList%transition(lineList%nLines))
		
		lineList%nActiveLines = 0
				
		do i = 1, lineList%nLines
			read(13,*) active, lineList%transition(i)%lambda0, lineList%transition(i)%lambdaNew, lineList%transition(i)%lambdaKurucz, lineList%transition(i)%gf, lineList%transition(i)%Elow, &
				lineList%transition(i)%geff, lineList%transition(i)%Gt, lineList%transition(i)%sigmaABO,&
				lineList%transition(i)%alphaABO, lineList%transition(i)%lambdaLeft, lineList%transition(i)%lambdaRight, lineList%transition(i)%lambdaContLeft, lineList%transition(i)%lambdaContRight
							
			if (active == 1 .and. lineList%transition(i)%lambda0 /= 0.d0 .and. lineList%transition(i)%Gt /= 0.d0) then
				lineList%nActiveLines = lineList%nActiveLines + 1
			else						
				lineList%transition(i)%lambda0 = 0.d0
			endif
			
			lineList%transition(i)%gf = 10.d0**lineList%transition(i)%gf
			lineList%transition(i)%Elow = lineList%transition(i)%Elow * PH * PC         ! Transform to erg
			lineList%transition(i)%frequency0 = PC / (lineList%transition(i)%lambdaNew * 1.d-8)
		enddo
		
		close(13)
				
		write(*,*) 'Number of active lines : ', lineList%nActiveLines
		
! Number of variables (etal, beta0, deltaLambda and B for each line, plus the hyperparameters)
		nVariables = lineList%nActiveLines * 6 + 4 + 3
		
! Test that the code is compiled to cope with this problem using automatic derivatives
		if (size(test%ip) /= nVariables) then
			write(*,FMT='(A)') 'Code not ready for this problem if using automatic differentiation.'
			write(*,FMT='(A,I5)') 'Recompile the code changing the NVARIABLES variable in the makefile to ', nVariables			
		endif
		
		
	end subroutine readLineList
	
!------------------------------------------------
! Read the atlas
!------------------------------------------------
	subroutine readAtlas(atlas)
	type(atlasType) :: atlas	
	integer :: i
	
		open(unit=13,file=trim(adjustl(atlas%file)),action='read',status='old',access='stream')
		read(13,pos=1) atlas%nLambda
		write(*,*) 'Reading atlas with ', atlas%nLambda, ' wavelength points...'
				
		allocate(atlas%lambda(atlas%nLambda),atlas%intensity(atlas%nLambda))
				
		read(13,pos=9) atlas%lambda		
  		read(13,pos=9+8*atlas%nLambda) atlas%intensity
  						
		close(13)
		
	end subroutine readAtlas
	
!------------------------------------------------
! Find indices
!------------------------------------------------
	subroutine setupLines(lineList,atlas)	
	type(lineListType) :: lineList
	type(atlasType) :: atlas
	integer :: i, limits(2), nTotal
	real(kind=8) :: a, b, x1, y1, x2, y2
		
		lineList%nLambdaTotal = 0
		lineList%nLambdaMax = 0
		nTotal = 0
		
		do i = 1, lineList%nLines
			
			lineList%transition(i)%active = .FALSE.
			
			if (lineList%transition(i)%lambda0 /= 0.d0) then
			
				nTotal = nTotal + 1
! Get ranges of the line
				lineList%transition(i)%active = .TRUE.
				limits = extremes(atlas%lambda, lineList%transition(i)%lambdaLeft, lineList%transition(i)%lambdaRight)
				write(*,FMT='(A,I3,A,I3,A,I7,A,I7,A,F9.3,A,F9.3)') 'Line ', i, '(', nTotal, ') - Range : ', limits(1), ' -> ', limits(2), ' - lambda: ', atlas%lambda(limits(1)), ' -> ', atlas%lambda(limits(2))
				
				lineList%transition(i)%nLambda = limits(2) - limits(1) + 1
				allocate(lineList%transition(i)%lambda(lineList%transition(i)%nLambda))
				allocate(lineList%transition(i)%frequency(lineList%transition(i)%nLambda))
				allocate(lineList%transition(i)%intensity(lineList%transition(i)%nLambda))
				allocate(lineList%transition(i)%observed(lineList%transition(i)%nLambda))
								
				lineList%transition(i)%lambda = atlas%lambda(limits(1):limits(2))
				lineList%transition(i)%observed = atlas%intensity(limits(1):limits(2))
				lineList%transition(i)%frequency = PC / (lineList%transition(i)%lambda * 1.d-8)
				
				lineList%nLambdaTotal = lineList%nLambdaTotal + lineList%transition(i)%nLambda				

! Normalize to the continuum
				limits = extremes(atlas%lambda, lineList%transition(i)%lambdaContLeft, lineList%transition(i)%lambdaContRight)
				x1 = atlas%lambda(limits(1))
				x2 = atlas%lambda(limits(2))
				y1 = atlas%intensity(limits(1))
				y2 = atlas%intensity(limits(2))
				
				a = (y1-y2) / (x1-x2)
				b = (y2*x1 - x2*y1) / (x1-x2)
				
				lineList%transition(i)%observed = lineList%transition(i)%observed / (a*lineList%transition(i)%lambda + b)
				
				if (lineList%transition(i)%nLambda > lineList%nLambdaMax) then
					lineList%nLambdaMax = lineList%transition(i)%nLambda
				endif
				
			endif
		
		enddo
		
		print *, 'Total number of wavelengths : ', lineList%nLambdaTotal
						
	end subroutine setupLines

	
end module ioModule
