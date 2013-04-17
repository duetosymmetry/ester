
c******************************************************************

	SUBROUTINE shell(n,x)

c routine public du module mod_numerique	
c tri d'un tableau

c entr�e
c	n : nombre de points

c entr�e/sortie
c	x : tableau

c	B.Pichon, P.Morel, D�partement Cassiop�e, O.C.A.

c----------------------------------------------------------------------

	USE mod_kind

	IMPLICIT NONE

	INTEGER, INTENT(in) :: n	
	REAL (kind=dp), INTENT(inout), DIMENSION(n) :: x
	
	INTEGER :: i, j
	
c----------------------------------------------------------------------	

2000	FORMAT(8es10.3)

	DO i=1,n-1
	 j=i-1+MINLOC(x(i:),dim=1)
	 x((/i,j/))=x((/j,i/))	
	ENDDO
		
	RETURN
	
	END SUBROUTINE shell
 
