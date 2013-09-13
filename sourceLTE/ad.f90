module dualClass
implicit none

	type, public :: dual	
		real(kind=8) :: rp
		real(kind=8), dimension(dualSize) :: ip   ! Real and infinitesimal part
	contains
		procedure, pass :: init
	end type dual
	
! Operator overload
	interface operator(+)
		module procedure additionDualDual, additionDualReal, additionRealDual
	end interface
	
	interface operator(-)
		module procedure substractionDualDual, substractionDualReal, substractionRealDual, signChangeDual
	end interface
	
	interface operator(*)
		module procedure productDualDual, productDualReal, productRealDual
	end interface
	
	interface operator(/)
		module procedure divisionDualDual, divisionDualReal, divisionRealDual
	end interface
	
	interface operator(**)
		module procedure powerDual, powerIntDual
	end interface
	
	interface exp
		module procedure expDual
	end interface
		
	interface sin
		module procedure sinDual
	end interface
	
	interface cos
		module procedure cosDual
	end interface
	
	interface log
		module procedure logDual
	end interface
	
	interface tanh
		module procedure tanhDual
	end interface
	
	
	contains

!---------------------------
! Creation of the dual number
! value and whichIndependent are optional arguments
! value is used to initialize the value of the number
! whichIndependent is used to set independent variables. If not present, the variable is dependent
!---------------------------
	subroutine init(this, value, whichIndependent)
	class(dual), intent(inout) :: this
	real(kind=8), optional :: value
	integer, optional :: whichIndependent
				
		this%ip = 0.d0
		
		if (present(value)) then
			this%rp = value
		else
			this%rp = 0.d0
		endif
		
		this%ip = 0.d0
		if (present(whichIndependent)) then
			this%ip(whichIndependent) = 1.d0
		endif
		
	end subroutine init
	
! Functions for defining the basic operations

! ---------------
! MULTIPLICATION
!----------------
	elemental function productDualDual(a,b) result(c)
		type (dual), intent(in) :: a, b
		type (dual) :: c
		c%rp = a%rp * b%rp
		c%ip = a%rp * b%ip + b%rp * a%ip
	end function productDualDual
	
	elemental function productDualReal(a,b) result(c)
		type (dual), intent(in) :: a
		real(kind=8), intent(in) :: b
		type (dual) :: c
		c%rp = a%rp * b
		c%ip = a%ip * b
	end function productDualReal
	
	elemental function productRealDual(a,b) result(c)
		real(kind=8), intent(in) :: a
		type (dual), intent(in) :: b		
		type (dual) :: c
		c%rp = b%rp * a
		c%ip = b%ip * a
	end function productRealDual

! ---------------
! DIVISION
!----------------
	elemental function divisionDualDual(a,b) result(c)
		type (dual), intent(in) :: a, b
		type (dual) :: c
		c%rp = a%rp / b%rp
		c%ip = (b%rp * a%ip - a%rp * b%ip) / b%rp**2
	end function divisionDualDual
	
	elemental function divisionDualReal(a,b) result(c)
		type (dual), intent(in) :: a
		real(kind=8), intent(in) :: b
		type (dual) :: c
		c%rp = a%rp / b
		c%ip = a%ip / b
	end function divisionDualReal
	
	elemental function divisionRealDual(a,b) result(c)
		real(kind=8), intent(in) :: a
		type (dual), intent(in) :: b		
		type (dual) :: c
		c%rp = a / b%rp
		c%ip = -a * b%ip / b%rp**2
	end function divisionRealDual

!---------------
! ADDITION
!---------------
	elemental function additionDualDual(a,b) result(c)
		type (dual), intent(in) :: a, b
		type (dual) :: c
		c%rp = a%rp + b%rp
		c%ip = a%ip + b%ip
	end function additionDualDual
		
	elemental function additionDualReal(a,b) result(c)
		type (dual), intent(in) :: a
		real(kind=8), intent(in) :: b
		type (dual) :: c
		c%rp = a%rp + b
		c%ip = a%ip
	end function additionDualReal
	
	elemental function additionRealDual(a,b) result(c)
		real(kind=8), intent(in) :: a
		type (dual), intent(in) :: b		
		type (dual) :: c
		c%rp = b%rp + a
		c%ip = b%ip
	end function additionRealDual

!---------------
! SUBSTRACTION
!---------------
	elemental function substractionDualDual(a,b) result(c)
		type (dual), intent(in) :: a, b
		type (dual) :: c
		c%rp = a%rp - b%rp
		c%ip = a%ip - b%ip
	end function substractionDualDual
	
	elemental function substractionDualReal(a,b) result(c)
		type (dual), intent(in) :: a
		real(kind=8), intent(in) :: b
		type (dual) :: c
		c%rp = a%rp - b
		c%ip = a%ip
	end function substractionDualReal
	
	elemental function substractionRealDual(a,b) result(c)
		real(kind=8), intent(in) :: a
		type (dual), intent(in) :: b		
		type (dual) :: c
		c%rp = b%rp - a
		c%ip = b%ip
	end function substractionRealDual
	
	elemental function signChangeDual(a) result(c)
		type (dual), intent(in) :: a
		type (dual) :: c
		c%rp = -a%rp
		c%ip = -a%ip
	end function signChangeDual

!---------------
! POWER
!---------------	
	elemental function powerDual(a,b) result(c)
		type (dual), intent(in) :: a
		real(kind=8), intent(in) :: b
		type (dual) :: c		
		c%rp = a%rp**b
		c%ip = b * a%rp**(b-1.d0) * a%ip		
	end function powerDual
	
	elemental function powerIntDual(a,b) result(c)
		type (dual), intent(in) :: a
		integer, intent(in) :: b
		type (dual) :: c		
		c%rp = a%rp**b
		c%ip = b * a%rp**(b-1.d0) * a%ip		
	end function powerIntDual
	
!---------------
! POWER
!---------------
	elemental function expDual(a) result(c)
		type (dual), intent(in) :: a
		type (dual) :: c				
		c%rp = exp(a%rp)
		c%ip = a%ip * exp(a%rp)
	end function expDual

!---------------
! SIN
!---------------
	elemental function sinDual(a) result(c)
		type (dual), intent(in) :: a
		type (dual) :: c
		c%rp = sin(a%rp)
		c%ip = a%ip * cos(a%rp)
	end function sinDual

!---------------
! COS
!---------------
	elemental function cosDual(a) result(c)
		type (dual), intent(in) :: a
		type (dual) :: c
		c%rp = cos(a%rp)
		c%ip = -a%ip * sin(a%rp)
	end function cosDual

!---------------
! LN
!---------------
	elemental function logDual(a) result(c)
		type (dual), intent(in) :: a
		type (dual) :: c
		c%rp = log(a%rp)
		c%ip = a%ip / a%rp
	end function logDual
		
!---------------
! TANH
!---------------
	elemental function tanhDual(a) result(c)
		type (dual), intent(in) :: a
		type (dual) :: c
		real(kind=8) :: t
		c%rp = tanh(a%rp)
		c%ip = 1.d0/cosh(a%rp)**2 * a%ip
	end function tanhDual
	
end module dualClass