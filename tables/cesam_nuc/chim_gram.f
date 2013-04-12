 
c****************************************************************

	SUBROUTINE chim_gram(xchim,dxchim)
	
c routine public du module mod_variables

c	transforme les abondances par mole ---> par gramme
	
c	Auteur: P.Morel, D�partement J.D. Cassini, O.C.A.
c	CESAM2k

c entr�es/sorties
c	xchim,dxchim : comp. chim et d�riv�e par mole ---> par gramme

c-----------------------------------------------------------------------

	USE mod_donnees, ONLY : nchim, nucleo
	USE mod_kind	

	IMPLICIT NONE
	
	REAL (kind=dp), INTENT(inout), DIMENSION(:) :: dxchim, xchim	

c--------------------------------------------------------------

	xchim=xchim*nucleo ; dxchim=dxchim*nucleo

	RETURN
	
	END SUBROUTINE chim_gram
