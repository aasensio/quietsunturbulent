module samplingModule
use globalModule
use likelihoodModule, only : negLogPosterior, writeHMCProcess, writeHMCProcessBurnin
use geodesiclmModule, only : maximumLikelihoodGeodesicLM
use mathsModule, only : sigmoid, invSigmoid, fvoigtScalar
implicit none
contains

!------------------------------------------------
! Carry out the sampling using the HMC
!------------------------------------------------
	subroutine doSampling(action)
	character(len=100) :: action
	real(kind=8), allocatable :: st(:), stepSize(:), logPGradient(:), logPGradientNew(:), st2(:), meanOld(:)
	real(kind=8) :: scaleFactor, logP, logP2
	integer :: seed, fbInt, maxStep, resume, nburn, i, nStepsBurn
	character(len=128) :: flPfx
	
		bestLogPosterior = -1.d100
		scaleFactor = 0.5d0
		seed = 1234
		fbInt = 10
		maxStep = 10
		if (action(1:4) == 'CONT') then
			resume = 1		
		else
			resume = 0
		endif
		nburn = 0
		flPfx = 'test'
		
		allocate(st(nVariables))
		allocate(stepSize(nVariables))
		allocate(parsOld(nVariables))
		allocate(parsInitial(nVariables))
		allocate(parsVariance(nVariables))
		allocate(parsMean(nVariables))
		allocate(bestPars(nVariables))
		
		allocate(beta0(lineList%nActiveLines))
		allocate(etal(lineList%nActiveLines))
		allocate(deltaLambda(lineList%nActiveLines))
		allocate(B(lineList%nActiveLines))
		allocate(damping(lineList%nActiveLines))
		allocate(sigma(lineList%nActiveLines))
		allocate(xB(2))
		
		allocate(beta02(lineList%nActiveLines))
		allocate(etal2(lineList%nActiveLines))
		allocate(deltaLambda2(lineList%nActiveLines))
		allocate(B2(lineList%nActiveLines))
		allocate(damping2(lineList%nActiveLines))
		allocate(sigma2(lineList%nActiveLines))
		allocate(xB2(2))
		allocate(profileH(lineList%nLambdaMax))
		allocate(profileL(lineList%nLambdaMax))
		allocate(v(lineList%nLambdaMax))
		allocate(Psi(lineList%nLambdaMax))
		allocate(synthesis(lineList%nLambdaMax))
		allocate(diffProfile(lineList%nLambdaMax))
		
		nStepsBurn = 250
				
		if (resume == 0) then
			call initialValues(st, stepSize)
			
 			stepSize = stepSize / maxStep
 											
			parsOld = st
			parsInitial = st
			nStep = 1
			parsVariance = 0.d0
			parsMean = 0.d0
					
			open(unit=20,file= flPfx(1:len_trim(flPfx))//".burnin",action='write',status='replace',access='stream')					
			call run_guided_hmc(nVariables,st,scaleFactor,maxStep,stepSize,flPfx(1:len_trim(flPfx)),seed,0,&
				fbInt, negLogPosterior, writeHMCProcessBurnin, nBurn, nStepsBurn)			
			close(20)
		endif
			
! Estimate width of distributions
		open(unit=20,file=flPfx(1:len_trim(flPfx))//".burnin",action='read',status='old',access='stream')
		parsMean = 0.d0
		parsVariance = 0.d0
		nStep = 1
		allocate(meanOld(nVariables))
		do i = 1, nStepsBurn
			read(20) st
			if (i > 100) then
				meanOld = parsMean
				parsMean = meanOld + (st - meanOld) / (nStep + 1.d0)		
				parsVariance = (nStep - 1.d0) / nStep * parsVariance + (st - meanOld)**2 / (nStep+1.d0)**2 + (st - meanOld)**2 / nStep
				nStep = nStep + 1
			endif
		enddo		
		deallocate(meanOld)
		
		close(20)
		
		stepSize = sqrt(parsVariance)
		st = parsMean
		
		open(unit=20,file='variances',action='write',status='replace')
		do i = 1, nVariables
			write(20,FMT='(I4,1X,F12.3,1X,F12.3,1X,F12.3)') i, parsMean(i), parsInitial(i), stepSize(i)
		enddo
		close(20)
				
		maxStep = 10
		scaleFactor = 0.2d0
		   					
  		stepSize = stepSize / maxStep
		
		if (resume == 0) then					
			open(unit=20,file= (trim(flPfx)//".extract"),action='write',status='replace',access='stream')			
		else
			open(unit=20,file='final.parameters',action='read',status='old',access='stream')
			read(20) st
			close(20)
			open(unit=20,file= (trim(flPfx)//".extract"),action='write',position='append',access='stream')			
		endif
		
		call run_guided_hmc(nVariables,st,scaleFactor,maxStep,stepSize,flPfx(1:len_trim(flPfx)),seed,0,&
			fbInt, negLogPosterior, writeHMCProcess, nBurn, nSteps)
		close(20)
		
		if (resume == 0) then
			open(unit=20,file='posterior.sizes',action='write',status='replace')
			write(20,*) nVariables, nSteps
			close(20)
			open(unit=20,file='final.parameters',action='write',status='replace',access='stream')
			write(20) st
			close(20)
			open(unit=20,file='best.parameters',action='write',status='replace')
			open(unit=21,file='best.profiles',action='write',status='replace')
			call writeBestProfiles(20,21)
			close(20)
			close(21)
		else
			open(unit=20,file='posterior.sizes',action='read',status='old')
			read(20,*) i, nStepsBurn
			close(20)
			open(unit=20,file='posterior.sizes',action='write',status='replace')
			write(20,*) nVariables, nStepsBurn+nSteps
			close(20)
			open(unit=20,file='final.parameters',action='write',status='replace',access='stream')
			write(20) st
			close(20)
			open(unit=20,file='best.parameters',action='write',status='replace')
			open(unit=21,file='best.profiles',action='write',status='replace')
			call writeBestProfiles(20,21)
			close(20)
			close(21)
		endif
		
	end subroutine doSampling
	
!------------------------------------------------------------------
! Set initial values for the parameters
!------------------------------------------------------------------
	subroutine initialValues(pars, stepSize)
	real(kind=8) :: pars(:), stepSize(:)
	integer :: loop, i, j
	real(kind=8) :: value, x(4), stddev, xInput(6)
	integer :: fail
					
		loop = 1

		do i = 1, lineList%nLines
			if (lineList%transition(i)%active) then				
				x = (/1.d0,0.2d0,0.05d0,0.1d0/)
				call maximumLikelihoodGeodesicLM(4, x, i, stddev)
				
				xInput(1) = 10.d0
				xInput(2:3) = x(2:3)
				xInput(4) = 100.d0
				xInput(5) = x(4)
				xInput(6) = stddev
				
				fail = 0
				do j = 1, 6
					if (xInput(j) > priorUpper(j)) then
						xInput(j) = 0.9*priorUpper(j) + 0.1*priorLower(j)
					endif
					if (xInput(j) < priorLower(j)) then
						xInput(j) = 0.9*priorLower(j) + 0.1*priorUpper(j)
					endif					
				enddo
				
				write(*,FMT='(A,I3,A,E12.3,4(A,F11.4))') 'Line ', i, ' ->  beta0=', x(1), ' - etal=', x(2),' - deltaL=', x(3), ' -a=', x(4), ' - sigma=', stddev
				
				do j = 1, 4
					
				enddo
																				
! beta0, etal, deltaL, B and damping
				pars(loop) = invSigmoid(xInput(1), priorLower(1), priorUpper(1))				                  ! B0
				pars(loop + lineList%nActiveLines) = invSigmoid(xInput(2), priorLower(2), priorUpper(2))       ! etal
				pars(loop + 2*lineList%nActiveLines) = invSigmoid(xInput(3), priorLower(3), priorUpper(3))     ! deltaL
				pars(loop + 3*lineList%nActiveLines) = invSigmoid(xInput(4), priorLower(4), priorUpper(4))    ! B
				pars(loop + 4*lineList%nActiveLines) = invSigmoid(xInput(5), priorLower(5), priorUpper(5))     ! damping
				pars(loop + 5*lineList%nActiveLines) = invSigmoid(xInput(6), priorLower(6), priorUpper(6))     ! sigma
				loop = loop + 1			
			endif			
		enddo
						
		stepSize = 1.d0
						
	end subroutine initialValues
	
!------------------------------------------------------------------
! Write profiles that produce the best MAP
!------------------------------------------------------------------
	subroutine writeBestProfiles(unitPars, unitProfiles)
	integer :: unitPars, unitProfiles, i, j, loop
	real(kind=8) :: deltaT2, prof, syn, beta0Local, etalLocal, deltaLambdaLocal, BLocal, aLocal
	
		loop = 1
		do i = 1, lineList%nLines
			if (lineList%transition(i)%active) then
				beta0Local = bestPars(loop)
				beta0Local = sigmoid(beta0Local, priorLower(1), priorUpper(1))
				
				etalLocal = bestPars(loop + lineList%nActiveLines)
				etalLocal = sigmoid(etalLocal, priorLower(2), priorUpper(2))
				
				deltaLambdaLocal = bestPars(loop + 2*lineList%nActiveLines)
				deltaLambdaLocal = sigmoid(deltaLambdaLocal, priorLower(3), priorUpper(3))
				
				BLocal = bestPars(loop + 3*lineList%nActiveLines)
				BLocal = sigmoid(BLocal, priorLower(4), priorUpper(4))
				
				aLocal = bestPars(loop + 4*lineList%nActiveLines)
				aLocal = sigmoid(aLocal, priorLower(5), priorUpper(5))
								
				deltaT2 = deltaLambdaLocal**2 + 4.d0 * lineList%transition(i)%Gt * (4.6686d-13 * lineList%transition(i)%lambda0**2 * BLocal )**2
								
				write(unitPars,*) beta0Local, etalLocal, deltaLambdaLocal, BLocal, aLocal
 				do j = 1, lineList%transition(i)%nLambda
 					
 					prof = fvoigtScalar(aLocal, (lineList%transition(i)%lambda(j) - lineList%transition(i)%lambdaNew) / sqrt(deltaT2)) / sqrt(PI)
 					syn = (1.d0 + beta0Local / (1.d0 + etalLocal * prof)) / (1.d0 + beta0Local)
 					
  					write(unitProfiles,*) lineList%transition(i)%lambda(j), syn, lineList%transition(i)%observed(j)
  										
				enddo
				
				loop = loop + 1
			endif
		enddo
	end subroutine writeBestProfiles
	
end module samplingModule