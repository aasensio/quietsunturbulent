program turbulent
use globalModule, only : atlas, lineList, nSteps, nVariables
use ioModule, only : readInput
use likelihoodModule, only : negLogPosterior
use samplingModule, only : doSampling
use dualClass

implicit none

	integer :: i, j, nargs
	character(len=100) :: action, arg
	
! Find the parameters passed to the program
	nargs = iargc()
	
	if (nargs <= 1) then
		print *, 'Parameters missing. The options are:'
		print *, '  - START nsteps -> start a MCMC with nsteps'
		print *, '  - CONT nsteps -> continue a previous MCMC and proceed another nsteps steps'		
		stop
	endif
	
	if (nargs == 2) then
		call getarg(1,action)
		call getarg(2,arg)
		read(arg,*) nSteps
	endif

	call readInput(lineList, atlas)
		
	call doSampling(action)
end program turbulent