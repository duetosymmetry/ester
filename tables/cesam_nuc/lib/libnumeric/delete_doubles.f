
c*********************************************************************

	SUBROUTINE delete_doubles(n,a)
	
c routine public du module mod_numerique

c suppression des r�p�titions dans un tableau
c adaptation de la dimension

c entr�e/sortie
c	a : tableau ALLOCATABLE
c	n : nombre de points

c	P.Morel, D�partement Cassiop�e, O.C.A.

c----------------------------------------------------------------

	USE mod_kind

	IMPLICIT NONE
		 
	REAL (kind=dp), INTENT(inout), ALLOCATABLE, DIMENSION(:) :: a
	INTEGER, INTENT(inout) :: n
	 
	REAL (kind=dp), ALLOCATABLE, DIMENSION(:) :: temp
	INTEGER :: i, newn

c-----------------------------------------------------------------------

2000	FORMAT(8es10.3)

c on ordonne le tableau pour que les multiplicit�s soient cons�cutives
	CALL shell(n,a)
	
c s'il y a r�p�tition
	newn=n ; i=1
	B1: DO
	 IF(i >= n)EXIT B1
	 IF(a(i) == a(i+1))THEN	 
	  a(i)=HUGE(1.d0) ; newn=newn-1
	 ENDIF
	 i=i+1 
	ENDDO B1
	
c	PRINT*,newn ; WRITE(*,2000)a

c �limination des mutiplicit�s, on r�ordonne et on coupe
	IF(newn /= n)THEN	
	 CALL shell(n,a)	 
	 ALLOCATE(temp(newn)) ; temp(1:newn)=a(1:newn) ; DEALLOCATE(a)
	 n=newn ; ALLOCATE(a(n)) ; a=temp ; DEALLOCATE(temp)
	ENDIF	 
	 
	RETURN
	 
	END SUBROUTINE delete_doubles
