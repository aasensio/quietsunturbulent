program turbulent
use globalModule, only : atmosphere, atlas, lineList
use ioModule, only : readInput
use synthModule, only : synthLines

implicit none
integer :: i

	call readInput(atmosphere, lineList, atlas)

	do i = 1, 20
		print *, i
		call synthLines(atmosphere,lineList)
	enddo

end program turbulent