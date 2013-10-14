module likelihoodModule
use globalModule
use mathsModule, only : sigmoid, diffSigmoid, invSigmoid, voigtZeeman, logSum
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
   real(kind=8) :: logP, logPrior2, deltaT2, ratioDelta2, k, logP1, logP2, der1, der2, dlogPdlogP1, dlogPdlogP2
   real(kind=8) :: dlogP1dmu1, dlogP1dmu2, dlogP1dsigma1, dlogP1dsigma2, dlogP1dp, dlogP1dB
   real(kind=8) :: dlogP2dmu1, dlogP2dmu2, dlogP2dsigma1, dlogP2dsigma2, dlogP2dp, dlogP2dB
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
			damping2(i) = trial(loop)
			loop = loop + 1			
		enddo
		damping2 = sigmoid(damping2, priorLower(4), priorUpper(4))
		
		do i = 1, lineList%nActiveLines
			sigma2(i) = trial(loop)
			loop = loop + 1			
		enddo
		sigma2 = sigmoid(sigma2, priorLower(5), priorUpper(5))
		
		B2 = trial(loop)
		loop = loop + 1		
		B2 = sigmoid(B2, priorLower(6), priorUpper(6))
									
		do i = 1, 2
			xB2(i) = trial(loop)			
			loop = loop + 1
		enddo
		xB2(1) = sigmoid(xB2(1), priorLower(7), priorUpper(7))   ! mu1_V
		xB2(2) = sigmoid(xB2(2), priorLower(8), priorUpper(8))   ! sigma1_V
				
		logPGradient = 0.d0
		logP = 0.d0
						
		loop = 1

!---------------------
! Log-normal prior for deltaV
!---------------------		
 		logP = logP - sum( (log(deltaV2) - xB2(1))**2 / (2.d0*xB2(2)**2) ) - sum(log(deltaV2)) - lineList%nActiveLines * log(xB2(2))
 				
! dIG/dB_i
 		logPGradient(2*lineList%nActiveLines+1:3*lineList%nActiveLines) = logPGradient(2*lineList%nActiveLines+1:3*lineList%nActiveLines) - 1.d0 / deltaV2 - (log(deltaV2) - xB2(1)) / (deltaV2*xB2(2)**2)
! dIG/dalpha
 		logPGradient(5*lineList%nActiveLines+2) = sum( (log(deltaV2) - xB2(1)) / xB2(2)**2 )
! dIG/dbeta
 		logPGradient(5*lineList%nActiveLines+3) = -lineList%nActiveLines / xB2(2) + sum( (log(deltaV2) - xB2(1))**2 / xB2(2)**3)
 		
! Jeffreys prior for sigma
  		logP = logP - 3.d0*log(xB2(2))  		
  		
  		
! dP/dbeta
  		logPGradient(5*lineList%nActiveLines+3) = logPGradient(5*lineList%nActiveLines+3) - 3.d0 / xB2(2)
 				
!---------------------
! Modified Jeffreys for beta0
!---------------------
 		logP = logP - sum(log(beta02 + priorLower(1))) 		
 		
! dP/dbeta0
 		logPGradient(1:lineList%nActiveLines) = logPGradient(1:lineList%nActiveLines) - 1.d0 / (beta02 + priorLower(1))
		
!---------------------
! Modified Jeffreys for sigma
!---------------------
 		logP = logP - sum(log(sigma2 + priorLower(5))) 		
 		
! dP/dsigma
 		logPGradient(4*lineList%nActiveLines+1:5*lineList%nActiveLines) = logPGradient(4*lineList%nActiveLines+1:5*lineList%nActiveLines) - 1.d0 / (sigma2 + priorLower(5))
 		
!---------------------
! Modified Jeffreys for damping
!---------------------
 		logP = logP - sum(log(damping2 + priorLower(4)))
 		
! dP/ddamping
 		logPGradient(3*lineList%nActiveLines+1:4*lineList%nActiveLines) = logPGradient(3*lineList%nActiveLines+1:4*lineList%nActiveLines) - 1.d0 / (damping2 + priorLower(4))
 		
!---------------------
! Data likelihood
!---------------------
 		do i = 1, lineList%nLines
 			if (lineList%transition(i)%active) then
 			
				deltaLambda2(loop) = deltaV2(loop) * lineList%transition(i)%lambda0 / 3.d5
				
				k = 4.d0 * lineList%transition(i)%Gt * (4.6686d-13 * lineList%transition(i)%lambda0**2)**2
				
				deltaT2 = sqrt(deltaLambda2(loop)**2 + k * B2**2)
 				
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
										
! dL/da
				diffProfile(1:lineList%transition(i)%nLambda) = (-2.d0/sqrt(PI) + 2.d0*aDamp*profileH(1:lineList%transition(i)%nLambda) + &
					2.d0*v(1:lineList%transition(i)%nLambda)*profileL(1:lineList%transition(i)%nLambda)) * (deltaLambda2(loop) / deltaT2)**2 / sqrt(PI)
									
				diffProfile(1:lineList%transition(i)%nLambda) = - beta02(loop) * etal2(loop) / ( (1.d0+beta02(loop)) * (1.d0+etal2(loop)*Psi(1:lineList%transition(i)%nLambda))**2 ) * &
					diffProfile(1:lineList%transition(i)%nLambda)
					
				logPGradient(3*lineList%nActiveLines+loop) = logPGradient(3*lineList%nActiveLines+loop) + &
					sum( (lineList%transition(i)%observed(1:lineList%transition(i)%nLambda) - synthesis(1:lineList%transition(i)%nLambda)) / sigma2(loop)**2 * diffProfile(1:lineList%transition(i)%nLambda) )
					
! dL/dsigma									
				logPGradient(4*lineList%nActiveLines+loop) = logPGradient(4*lineList%nActiveLines+loop) + &
					sum( (lineList%transition(i)%observed(1:lineList%transition(i)%nLambda) - synthesis(1:lineList%transition(i)%nLambda))**2 / sigma2(loop)**3 - 1.d0 / sigma2(loop) )
					
					
! dL/dB
				diffProfile(1:lineList%transition(i)%nLambda) = (-2.d0*v(1:lineList%transition(i)%nLambda)*profileH(1:lineList%transition(i)%nLambda) + 2.d0*aDamp*profileL(1:lineList%transition(i)%nLambda)) * &
					(-k * v(1:lineList%transition(i)%nLambda) * B2 / deltaT2**2) + &
					(-2.d0/sqrt(PI) + 2.d0*aDamp*profileH(1:lineList%transition(i)%nLambda) + 2.d0*v(1:lineList%transition(i)%nLambda)*profileL(1:lineList%transition(i)%nLambda)) * &
					(-damping2(loop) * k * B2 * deltaLambda2(loop) / deltaT2**3)
					
				diffProfile(1:lineList%transition(i)%nLambda) = (-k * B2 * deltaLambda2(loop) / deltaT2**3) * profileH(1:lineList%transition(i)%nLambda) / sqrt(PI) + &
					diffProfile(1:lineList%transition(i)%nLambda) * ratioDelta2 / sqrt(PI)
				
				diffProfile(1:lineList%transition(i)%nLambda) = - beta02(loop) * etal2(loop) / ( (1.d0+beta02(loop)) * (1.d0+etal2(loop)*Psi(1:lineList%transition(i)%nLambda))**2 ) * &
					diffProfile(1:lineList%transition(i)%nLambda)
					
				logPGradient(5*lineList%nActiveLines+1) = logPGradient(5*lineList%nActiveLines+1) + &
					sum( (lineList%transition(i)%observed(1:lineList%transition(i)%nLambda) - synthesis(1:lineList%transition(i)%nLambda)) / sigma2(loop)**2 * diffProfile(1:lineList%transition(i)%nLambda) )

					
 				loop = loop + 1
 			endif
 		enddo

! Use chain rule to take into account the sigmoid transformation
		do i = 1, 5
			logPGradient((i-1)*lineList%nActiveLines+1:i*lineList%nActiveLines) = logPGradient((i-1)*lineList%nActiveLines+1:i*lineList%nActiveLines) * &
				diffSigmoid(trial((i-1)*lineList%nActiveLines+1:i*lineList%nActiveLines), priorLower(i), priorUpper(i))
		enddo
		do i = 1, 3
			logPGradient(5*lineList%nActiveLines+i) = logPGradient(5*lineList%nActiveLines+i) * diffSigmoid(trial(5*lineList%nActiveLines+i), priorLower(6+i-1), priorUpper(6+i-1))
		enddo
						
		if (logP > bestLogPosterior) then
			bestPars = trial
			bestLogPosterior = logP
		endif
		
		logP = -logP
		logPGradient = -logPGradient				
										
	end subroutine negLogPosterior
	
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
				
		do i = 1, 5
			x2((i-1)*lineList%nActiveLines+1:i*lineList%nActiveLines) = sigmoid(x((i-1)*lineList%nActiveLines+1:i*lineList%nActiveLines), priorLower(i), priorUpper(i))
		enddo
		do i = 1, 3
			x2(5*lineList%nActiveLines+i) = sigmoid(x(5*lineList%nActiveLines+i), priorLower(6+i-1), priorUpper(6+i-1))
		enddo				
			
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