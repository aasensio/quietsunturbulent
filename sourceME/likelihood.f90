module likelihoodModule
use globalModule
use mathsModule, only : sigmoid, diffSigmoid, invSigmoid, voigtZeeman
use dualClass
implicit none

contains

!------------------------------------------------
! Negative log-posterior and its gradient
!------------------------------------------------
 	subroutine negLogPosterior(nVariables,trial,logP,logPGradient)
	integer :: nVariables
	real(kind=8), dimension(nVariables) :: trial
   real(kind=8), dimension(nVariables) :: logPGradient
   real(kind=8) :: logP, logPrior2, deltaT2, ratioDelta2, k
	type(dual) :: logPosterior, logPrior, prof, syn, deltaT, ratioDelta
	integer :: i, j, loop, ierr
										

		loop = 1		 		
		do i = 1, lineList%nActiveLines
			beta02(i) = trial(loop)
			loop = loop + 1						
		enddo		
		beta02 = sigmoid(beta02, priorLower(1), priorUpper(1))
		   					
		do i = 1, lineList%nActiveLines
			etal2(i) = trial(loop)
			loop = loop + 1			
		enddo
		etal2 = sigmoid(etal2, priorLower(2), priorUpper(2))
		
		do i = 1, lineList%nActiveLines
			deltaV2(i) = trial(loop)
			loop = loop + 1			
		enddo
		deltaV2 = sigmoid(deltaV2, priorLower(3), priorUpper(3))		
						
		do i = 1, lineList%nActiveLines
			B2(i) = trial(loop)
			loop = loop + 1			
		enddo
		B2 = sigmoid(B2, priorLower(4), priorUpper(4))
		
		do i = 1, lineList%nActiveLines
			damping2(i) = trial(loop)
			loop = loop + 1			
		enddo
		damping2 = sigmoid(damping2, priorLower(5), priorUpper(5))
		
		do i = 1, lineList%nActiveLines
			sigma2(i) = trial(loop)
			loop = loop + 1			
		enddo
		sigma2 = sigmoid(sigma2, priorLower(6), priorUpper(6))
							
		do i = 1, 4
			xB2(i) = trial(loop)
			loop = loop + 1
		enddo
		xB2(1) = sigmoid(xB2(1), priorLower(7), priorUpper(7))
		xB2(2) = sigmoid(xB2(2), priorLower(8), priorUpper(8))
		xB2(3) = sigmoid(xB2(3), priorLower(9), priorUpper(9))
		xB2(4) = sigmoid(xB2(4), priorLower(10), priorUpper(10))
															
		logPGradient = 0.d0
		logP = 0.d0
				
		loop = 1

!---------------------
! Log-normal prior for B
!---------------------
 		logP = logP - sum( (log(B2) - xB2(1))**2 / (2.d0*xB2(2)**2) ) - sum(log(B2)) - lineList%nActiveLines * log(xB2(2))
		
! dIG/dB_i
 		logPGradient(3*lineList%nActiveLines+1:4*lineList%nActiveLines) = logPGradient(3*lineList%nActiveLines+1:4*lineList%nActiveLines) - 1.d0 / B2 - (log(B2) - xB2(1)) / (B2*xB2(2)**2) 
! dIG/dalpha
 		logPGradient(6*lineList%nActiveLines+1) = sum( (log(B2) - xB2(1)) / xB2(2)**2 )
! dIG/dbeta
 		logPGradient(6*lineList%nActiveLines+2) = -lineList%nActiveLines / xB2(2) + sum( (log(B2) - xB2(1))**2 / xB2(2)**3)
 		
! Jeffreys prior for sigma
  		logP = logP - 3.d0*log(xB2(2))

! dP/dbeta
  		logPGradient(6*lineList%nActiveLines+2) = logPGradient(6*lineList%nActiveLines+2) - 3.d0 / xB2(2)
  		

!---------------------
! Log-normal prior for deltaV
!---------------------		
 		logP = logP - sum( (log(deltaV2) - xB2(3))**2 / (2.d0*xB2(4)**2) ) - sum(log(deltaV2)) - lineList%nActiveLines * log(xB2(4))
		
! dIG/dB_i
 		logPGradient(2*lineList%nActiveLines+1:3*lineList%nActiveLines) = logPGradient(2*lineList%nActiveLines+1:3*lineList%nActiveLines) - 1.d0 / deltaV2 - (log(deltaV2) - xB2(3)) / (deltaV2*xB2(4)**2) 
! dIG/dalpha
 		logPGradient(6*lineList%nActiveLines+3) = sum( (log(deltaV2) - xB2(3)) / xB2(4)**2 )
! dIG/dbeta
 		logPGradient(6*lineList%nActiveLines+4) = -lineList%nActiveLines / xB2(4) + sum( (log(deltaV2) - xB2(3))**2 / xB2(4)**3)
 		
! Jeffreys prior for sigma
  		logP = logP - 3.d0*log(xB2(4))

! dP/dbeta
  		logPGradient(6*lineList%nActiveLines+4) = logPGradient(6*lineList%nActiveLines+4) - 3.d0 / xB2(4)
 		
!---------------------
! Inverse Gamma prior for B
!---------------------
!  		logP = logP + sum( -(xB2(1)+1.d0) * log(B2) - xB2(2) / B2 ) + lineList%nActiveLines * xB2(1) * log(xB2(2)) - lineList%nActiveLines * alngam(xB2(1), ierr)
! 		
! ! dIG/dB_i
!  		logPGradient(3*lineList%nActiveLines+1:4*lineList%nActiveLines) = logPGradient(3*lineList%nActiveLines+1:4*lineList%nActiveLines) + xB2(2) / B2**2 - (1.d0+xB2(1)) / B2
! ! dIG/dalpha
!  		logPGradient(6*lineList%nActiveLines+1) = lineList%nActiveLines * log(xB2(2)) - sum(log(B2)) - lineList%nActiveLines * digama(xB2(1), ierr)
! ! dIG/dbeta
!  		logPGradient(6*lineList%nActiveLines+2) = sum(-1.d0 / B2) + lineList%nActiveLines * xB2(1) / xB2(2)
!  		
!  		logP = logP - (xB2(1)-2.d0) - (xB2(1)-1.d0)
!  		logPGradient(6*lineList%nActiveLines+1) = logPGradient(6*lineList%nActiveLines+1) - 2.d0 		

!---------------------
! Gaussian prior for B
!---------------------
!  		logP = logP - 0.5d0 * sum( (B2 - xB2(1))**2 / xB2(2)**2) - lineList%nActiveLines * log(xB2(2))
!  				
! ! dIG/dB_i
!  		logPGradient(3*lineList%nActiveLines+1:4*lineList%nActiveLines) = logPGradient(3*lineList%nActiveLines+1:4*lineList%nActiveLines) - (B2 - xB2(1)) / xB2(2)**2
! ! dIG/dalpha
!  		logPGradient(6*lineList%nActiveLines+1) = sum( (-xB2(1) + B2) / xB2(2)**2 )
! ! dIG/dbeta
!  		logPGradient(6*lineList%nActiveLines+2) = -lineList%nActiveLines / xB2(2) + sum((B2 - xB2(1))**2 / xB2(2)**3)
! ! Jeffreys prior for sigma
!   		logP = logP - log(xB2(2))
! 
! ! dP/dbeta
!   		logPGradient(6*lineList%nActiveLines+2) = logPGradient(6*lineList%nActiveLines+2) - 1.d0 / xB2(2)
		
		
!---------------------
! Modified Jeffreys for beta0
!---------------------
 		logP = logP - sum(log(beta02 + priorLower(1)))

! dP/dbeta0
 		logPGradient(1:lineList%nActiveLines) = logPGradient(1:lineList%nActiveLines) - 1.d0 / (beta02 + priorLower(1))
		
!---------------------
! Modified Jeffreys for sigma
!---------------------
 		logP = logP - sum(log(sigma2 + priorLower(6)))

! dP/dsigma
 		logPGradient(5*lineList%nActiveLines+1:6*lineList%nActiveLines) = logPGradient(5*lineList%nActiveLines+1:6*lineList%nActiveLines) - 1.d0 / (sigma2 + priorLower(6))
 		
!---------------------
! Modified Jeffreys for damping
!---------------------
 		logP = logP - sum(log(damping2 + priorLower(5)))

! dP/ddamping
 		logPGradient(4*lineList%nActiveLines+1:5*lineList%nActiveLines) = logPGradient(4*lineList%nActiveLines+1:5*lineList%nActiveLines) - 1.d0 / (damping2 + priorLower(5))
		
		
!---------------------
! Data likelihood
!---------------------
 		do i = 1, lineList%nLines
 			if (lineList%transition(i)%active) then
 			
				deltaLambda2(loop) = deltaV2(loop) * lineList%transition(i)%lambda0 / 3.d5
				
				k = 4.d0 * lineList%transition(i)%Gt * (4.6686d-13 * lineList%transition(i)%lambda0**2)**2
				
				deltaT2 = sqrt(deltaLambda2(loop)**2 + k * B2(loop)**2)
 				
 				ratioDelta2 = deltaLambda2(loop) / deltaT2
 				
! Voigt and dispersion profiles
				v(1:lineList%transition(i)%nLambda) = (lineList%transition(i)%lambda - lineList%transition(i)%lambdaNew) / deltaT2
				aDamp = ratioDelta2 * damping2(loop)
 				call voigtZeeman(aDamp, v(1:lineList%transition(i)%nLambda), profileH(1:lineList%transition(i)%nLambda), profileL(1:lineList%transition(i)%nLambda))
 				Psi(1:lineList%transition(i)%nLambda) = profileH(1:lineList%transition(i)%nLambda) * ratioDelta2 / sqrt(PI)
 				
 				synthesis(1:lineList%transition(i)%nLambda) = (1.d0 + beta02(loop) / (1.d0 + etal2(loop) * Psi(1:lineList%transition(i)%nLambda))) / (1.d0 + beta02(loop))
 				
 				logP = logP - 0.5d0 * sum( (lineList%transition(i)%observed(1:lineList%transition(i)%nLambda) - synthesis(1:lineList%transition(i)%nLambda))**2 / sigma2(loop)**2 ) - &
					lineList%transition(i)%nLambda * log(sigma2(loop))
					
! dL/dbeta0
				diffProfile(1:lineList%transition(i)%nLambda) = 1.d0 / ( (1.d0+beta02(loop))*(1.d0+etal2(loop)*Psi(1:lineList%transition(i)%nLambda)) ) - &
					synthesis(1:lineList%transition(i)%nLambda) / (1.d0 + beta02(loop))
					
				logPGradient(loop) = logPGradient(loop) + &
					sum( (lineList%transition(i)%observed(1:lineList%transition(i)%nLambda) - synthesis(1:lineList%transition(i)%nLambda)) / sigma2(loop)**2 * diffProfile(1:lineList%transition(i)%nLambda) )
					
! dL/detal
				diffProfile(1:lineList%transition(i)%nLambda) = - beta02(loop) * Psi(1:lineList%transition(i)%nLambda) / ( (1.d0+beta02(loop)) * (1.d0+etal2(loop)*Psi(1:lineList%transition(i)%nLambda))**2 )
					
				logPGradient(lineList%nActiveLines+loop) = logPGradient(lineList%nActiveLines+loop) + &
					sum( (lineList%transition(i)%observed(1:lineList%transition(i)%nLambda) - synthesis(1:lineList%transition(i)%nLambda)) / sigma2(loop)**2 * diffProfile(1:lineList%transition(i)%nLambda) )
					
! dL/ddeltaL
				diffProfile(1:lineList%transition(i)%nLambda) = (-2.d0*v(1:lineList%transition(i)%nLambda)*profileH(1:lineList%transition(i)%nLambda) + 2.d0*aDamp*profileL(1:lineList%transition(i)%nLambda)) * &
					(-deltaLambda2(loop) * v(1:lineList%transition(i)%nLambda) / deltaT2**2) + &
					(-2.d0/sqrt(PI) + 2.d0*aDamp*profileH(1:lineList%transition(i)%nLambda) + 2.d0*v(1:lineList%transition(i)%nLambda)*profileL(1:lineList%transition(i)%nLambda)) * &
					damping2(loop) * (1.d0/deltaT2 - deltaLambda2(loop)**2 / deltaT2**3)
					
				diffProfile(1:lineList%transition(i)%nLambda) = (1.d0/deltaT2 - deltaLambda2(loop)**2 / deltaT2**3) * profileH(1:lineList%transition(i)%nLambda) / sqrt(PI) + &
					diffProfile(1:lineList%transition(i)%nLambda) * ratioDelta2 / sqrt(PI)
				
				diffProfile(1:lineList%transition(i)%nLambda) = - beta02(loop) * etal2(loop) / ( (1.d0+beta02(loop)) * (1.d0+etal2(loop)*Psi(1:lineList%transition(i)%nLambda))**2 ) * &
					diffProfile(1:lineList%transition(i)%nLambda)
					
				logPGradient(2*lineList%nActiveLines+loop) = logPGradient(2*lineList%nActiveLines+loop) + &
					sum( (lineList%transition(i)%observed(1:lineList%transition(i)%nLambda) - synthesis(1:lineList%transition(i)%nLambda)) / sigma2(loop)**2 * diffProfile(1:lineList%transition(i)%nLambda) ) * &
					lineList%transition(i)%lambda0 / 3.d5
					
! dL/dB
				diffProfile(1:lineList%transition(i)%nLambda) = (-2.d0*v(1:lineList%transition(i)%nLambda)*profileH(1:lineList%transition(i)%nLambda) + 2.d0*aDamp*profileL(1:lineList%transition(i)%nLambda)) * &
					(-k * v(1:lineList%transition(i)%nLambda) * B2(loop) / deltaT2**2) + &
					(-2.d0/sqrt(PI) + 2.d0*aDamp*profileH(1:lineList%transition(i)%nLambda) + 2.d0*v(1:lineList%transition(i)%nLambda)*profileL(1:lineList%transition(i)%nLambda)) * &
					(-damping2(loop) * k * B2(loop) * deltaLambda2(loop) / deltaT2**3)
					
				diffProfile(1:lineList%transition(i)%nLambda) = (-k * B2(loop) * deltaLambda2(loop) / deltaT2**3) * profileH(1:lineList%transition(i)%nLambda) / sqrt(PI) + &
					diffProfile(1:lineList%transition(i)%nLambda) * ratioDelta2 / sqrt(PI)
				
				diffProfile(1:lineList%transition(i)%nLambda) = - beta02(loop) * etal2(loop) / ( (1.d0+beta02(loop)) * (1.d0+etal2(loop)*Psi(1:lineList%transition(i)%nLambda))**2 ) * &
					diffProfile(1:lineList%transition(i)%nLambda)
					
				logPGradient(3*lineList%nActiveLines+loop) = logPGradient(3*lineList%nActiveLines+loop) + &
					sum( (lineList%transition(i)%observed(1:lineList%transition(i)%nLambda) - synthesis(1:lineList%transition(i)%nLambda)) / sigma2(loop)**2 * diffProfile(1:lineList%transition(i)%nLambda) )
					
! dL/da
				diffProfile(1:lineList%transition(i)%nLambda) = (-2.d0/sqrt(PI) + 2.d0*aDamp*profileH(1:lineList%transition(i)%nLambda) + &
					2.d0*v(1:lineList%transition(i)%nLambda)*profileL(1:lineList%transition(i)%nLambda)) * (deltaLambda2(loop) / deltaT2)**2 / sqrt(PI)
									
				diffProfile(1:lineList%transition(i)%nLambda) = - beta02(loop) * etal2(loop) / ( (1.d0+beta02(loop)) * (1.d0+etal2(loop)*Psi(1:lineList%transition(i)%nLambda))**2 ) * &
					diffProfile(1:lineList%transition(i)%nLambda)
					
				logPGradient(4*lineList%nActiveLines+loop) = logPGradient(4*lineList%nActiveLines+loop) + &
					sum( (lineList%transition(i)%observed(1:lineList%transition(i)%nLambda) - synthesis(1:lineList%transition(i)%nLambda)) / sigma2(loop)**2 * diffProfile(1:lineList%transition(i)%nLambda) )
					
! dL/dsigma									
				logPGradient(5*lineList%nActiveLines+loop) = logPGradient(5*lineList%nActiveLines+loop) + &
					sum( (lineList%transition(i)%observed(1:lineList%transition(i)%nLambda) - synthesis(1:lineList%transition(i)%nLambda))**2 / sigma2(loop)**3 - 1.d0 / sigma2(loop) )

					
 				loop = loop + 1
 			endif
 		enddo

! Use chain rule to take into account the sigmoid transformation
		do i = 1, 6
			logPGradient((i-1)*lineList%nActiveLines+1:i*lineList%nActiveLines) = logPGradient((i-1)*lineList%nActiveLines+1:i*lineList%nActiveLines) * &
				diffSigmoid(trial((i-1)*lineList%nActiveLines+1:i*lineList%nActiveLines), priorLower(i), priorUpper(i))
		enddo
		do i = 1, 4
			logPGradient(6*lineList%nActiveLines+i) = logPGradient(6*lineList%nActiveLines+i) * diffSigmoid(trial(6*lineList%nActiveLines+i), priorLower(7+i-1), priorUpper(7+i-1))
		enddo
				
		if (logP > bestLogPosterior) then
			bestPars = trial
			bestLogPosterior = logP
		endif
		
		logP = -logP
		logPGradient = -logPGradient				
								
	end subroutine negLogPosterior
	
!------------------------------------------------
! Negative log-posterior and its gradient
!------------------------------------------------
 	subroutine negLogPosteriorAutomatic(nVariables,trial,logP,logPGradient)
	integer :: nVariables
	real(kind=8), dimension(nVariables) :: trial
   real(kind=8), dimension(nVariables) :: logPGradient
   real(kind=8) :: logP, logPrior2, deltaT2, ratioDelta2, k
	type(dual) :: logPosterior, logPrior, prof, syn, deltaT, ratioDelta
	integer :: i, j, loop, ierr
										
! Insert parameters into the dual variables to compute automatic derivatives
 		loop = 1		 		
		do i = 1, lineList%nActiveLines
			call beta0(i)%init(trial(loop), loop)
			loop = loop + 1						
		enddo		
		beta0 = priorLower(1) + (priorUpper(1)-priorLower(1)) / (1.d0 + exp(-beta0))
		   					
		do i = 1, lineList%nActiveLines
			call etal(i)%init(trial(loop), loop)
			loop = loop + 1			
		enddo
		etal = priorLower(2) + (priorUpper(2)-priorLower(2)) / (1.d0 + exp(-etal))
		
		do i = 1, lineList%nActiveLines
			call deltaLambda(i)%init(trial(loop), loop)
			loop = loop + 1			
		enddo
		deltaLambda = priorLower(3) + (priorUpper(3)-priorLower(3)) / (1.d0 + exp(-deltaLambda))
		
		do i = 1, lineList%nActiveLines
			call B(i)%init(trial(loop), loop)			
			loop = loop + 1			
		enddo
		B = priorLower(4) + (priorUpper(4)-priorLower(4)) / (1.d0 + exp(-B))
		
		do i = 1, lineList%nActiveLines
			call damping(i)%init(trial(loop), loop)
			loop = loop + 1			
		enddo
		damping = priorLower(5) + (priorUpper(5)-priorLower(5)) / (1.d0 + exp(-damping))
		
		do i = 1, lineList%nActiveLines
			call sigma(i)%init(trial(loop), loop)
			loop = loop + 1			
		enddo
		sigma = priorLower(6) + (priorUpper(6)-priorLower(6)) / (1.d0 + exp(-sigma))
							
		do i = 1, 2			
			call xB(i)%init(trial(loop), loop)
			loop = loop + 1
		enddo
		xB(1) = priorLower(7) + (priorUpper(7)-priorLower(7)) / (1.d0 + exp(-xB(1)))
		xB(2) = priorLower(8) + (priorUpper(8)-priorLower(8)) / (1.d0 + exp(-xB(2)))
		
		
										
		call logPosterior%init()
		call logPrior%init()
		call deltaT%init()
		call ratioDelta%init()
		call prof%init()
		call syn%init()
		
		logPGradient = 0.d0
		logP = 0.d0
				
!----------------
! Data likelihood + prior
!----------------
		loop = 1
		do i = 1, lineList%nLines
			if (lineList%transition(i)%active) then
			
! Inverse Gamma prior for B
   				logPrior = logPrior + xB(1) * log(xB(2)) - (xB(1)+1.d0) * log(B(loop)) - xB(2) / B(loop) - lnGamma(xB(1))
  				
! Modified Jeffreys for beta0
   				logPrior = logPrior - log(beta0(loop) + priorLower(1))
  				  				
! Modified Jeffreys for sigma
   				logPrior = logPrior - log(sigma(loop) + priorLower(6))
  				
! Modified Jeffreys for damping
   				logPrior = logPrior - log(damping(loop) + priorLower(5))
				
! Data likelihood
 				deltaT = sqrt(deltaLambda(loop)**2 + 4.d0 * lineList%transition(i)%Gt * (4.6686d-13 * lineList%transition(i)%lambda0**2 * B(loop) )**2)
 				
				ratioDelta = deltaLambda(loop) / deltaT
				
 				do j = 1, lineList%transition(i)%nLambda
 					
 					prof = voigt(ratioDelta * damping(loop), (lineList%transition(i)%lambda(j) - lineList%transition(i)%lambdaNew) / deltaT) * ratioDelta / sqrt(PI)
 					syn = (1.d0 + beta0(loop) / (1.d0 + etal(loop) * prof)) / (1.d0 + beta0(loop))
 					
   					logPosterior = logPosterior - 0.5d0 * (lineList%transition(i)%observed(j) - syn)**2 / sigma(loop)**2 - log(sigma(loop))
  										
				enddo				

				loop = loop + 1		
			endif
		enddo
		
		if (logPosterior%rp > bestLogPosterior) then
			bestPars = trial
			bestLogPosterior = logPosterior%rp
		endif
		
		logPosterior = logPrior + logPosterior
		
		logPosterior = -logPosterior
		logP = logPosterior%rp				
		logPGradient = logPosterior%ip
											
	end subroutine negLogPosteriorAutomatic


!------------------------------------------------
! A subroutine to write the extract file
! I have assumed that the unit=20 is opened for
! writing (append) earlier. In general only write
! those parametes which are estimated. The files
! can be really big depending on the dimensionality
!------------------------------------------------
	subroutine writeHMCProcess(nVariables,x,v,g)
	integer :: nVariables
	real(kind=8), dimension(nVariables) :: x, x2
   real(kind=8), dimension(nVariables) :: g, meanOld
   real(kind=8) :: v, xwrite(4)
	integer i
				
		do i = 1, 6
			x2((i-1)*lineList%nActiveLines+1:i*lineList%nActiveLines) = sigmoid(x((i-1)*lineList%nActiveLines+1:i*lineList%nActiveLines), priorLower(i), priorUpper(i))
		enddo
		x2(6*lineList%nActiveLines+1) = sigmoid(x(6*lineList%nActiveLines+1), priorLower(7), priorUpper(7))
		x2(6*lineList%nActiveLines+2) = sigmoid(x(6*lineList%nActiveLines+2), priorLower(8), priorUpper(8))
		x2(6*lineList%nActiveLines+3) = sigmoid(x(6*lineList%nActiveLines+3), priorLower(9), priorUpper(9))
		x2(6*lineList%nActiveLines+4) = sigmoid(x(6*lineList%nActiveLines+4), priorLower(10), priorUpper(10))
			
 		write(20) x2
				
	end subroutine writeHMCProcess
	
!------------------------------------------------
! A subroutine to write the extract file
! I have assumed that the unit=20 is opened for
! writing (append) earlier. In general only write
! those parametes which are estimated. The files
! can be really big depending on the dimensionality
!------------------------------------------------
	subroutine writeHMCProcessBurnin(nVariables,x,v,g)
	integer :: nVariables
	real(kind=8), dimension(nVariables) :: x
   real(kind=8), dimension(nVariables) :: g, meanOld
   real(kind=8) :: v, xwrite(8)
	integer i
							
 		write(20) x		
						
	end subroutine writeHMCProcessBurnin
		
end module likelihoodModule