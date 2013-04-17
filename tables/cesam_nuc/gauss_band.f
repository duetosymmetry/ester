
c**************************************************************

	SUBROUTINE gauss_band(a,b,indpc,nl,n,nemr,ns,inversible)

c subroutine public du module mod_numerique

c r�solution par la methode du pivot d'un syst�me lin�aire bande
c en ne gardant que les coefficients non identiquement nuls,
c avec un nombre d'�quations sup�rieur au rang comme dans la m�thode de 
c Galerkin (convient �galement, et sans surco�t, si nl=n)
c il y a ns seconds membres

c m�thode: pivot partiel avec equilibrage

c Les indpc (indices de premi�re colonne) ne sont pas n�cessairement croissants.

c A l'issue du pivotage, le rang �tant n, les nl-n derni�res
c �quations sont 0=0 aux arrondis pr�s

c Auteur: P.Morel, Departement J.D. Cassini, O.C.A.

c entr�es
c	n: nombre de colonnes = rang
c	nl : nombre de lignes >= au rang = n
c	nemr : longueur max. d'une ligne
c	ns: nombre de seconds membres

c entr�es/sorties
c	a(nl,nemr) : matrice des coefficients non identiquement nuls
c	indpc(nl) : indice de la premi�re colonne de chaque ligne
c	b(ns,nl) : seconds membres de 1 � nl / solution de 1 � n

c sortie
c	inversible=.true. : solution obtenue

c------------------------------------------------------------------

	USE mod_kind
	
	IMPLICIT NONE
	
	INTEGER, INTENT(in) :: n, nemr, nl, ns	
	REAL (kind=dp), INTENT(inout), DIMENSION (:,:) :: a, b	
	INTEGER, INTENT(inout), DIMENSION(:) :: indpc
		
	LOGICAL, INTENT(out) :: inversible	
	
	REAL (kind=dp) :: ai, pivot, sto
	INTEGER :: i, ipivot, isto, j, k, ligne

c----------------------------------------------------------------------

2000	FORMAT(8es10.3)

c	PRINT*,'entr�e gaus_band' ; PRINT*,indpc
c	DO i=1,nl
c	 WRITE(*,2000)(a(i,j),j=1,nemr),(b(j,i),j=1,ns)
c	ENDDO
	 
c �quilibrage de la matrice
	DO ligne=1,nl
	 sto=MAXVAL(ABS(a(ligne,:)))
	 IF(sto /= 0.d0)THEN
	  a(ligne,:)=a(ligne,:)/sto ; b(:,ligne)=b(:,ligne)/sto
	 ENDIF
	ENDDO	!ligne

c mise en ordre des lignes pour que indpc soit croissant
	DO ligne=1,nl-1
	 IF(indpc(ligne) > indpc(ligne+1))THEN
	  i=ligne
	  DO WHILE(i >= 1)
	   IF(indpc(i) > indpc(i+1))THEN
	    DO j=1,ns
	     sto=b(j,i) ; b(j,i)=b(j,i+1) ; b(j,i+1)=sto
	    ENDDO	!j
	    isto=indpc(i) ; indpc(i)=indpc(i+1) ; indpc(i+1)=isto
	    DO j=1,nemr
	     sto=a(i,j) ; a(i,j)=a(i+1,j) ; a(i+1,j)=sto
	    ENDDO	!j
	   ENDIF
	   i=i-1
	  ENDDO		!while
	 ENDIF
	ENDDO		!ligne
	   
c �limination de gauss avec pivots
	B1: DO ligne=1,n		!�limination jusqu'� ligne n=rang
c	 PRINT* ; PRINT*,ligne ; PRINT*
	 pivot=0.d0
	 B2: DO i=ligne,nl		!recherche du pivot jusqu'� nl
	  IF(indpc(i) > ligne)EXIT B2
c	  PRINT*,i,a(i,1)
	  IF(pivot < abs(a(i,1)))THEN
	   pivot=abs(a(i,1)) ; ipivot=i
c	   PRINT*,'pivot',pivot,ipivot
	  ENDIF	   
	 ENDDO B2
	 inversible=pivot > 0.d0
	 IF(.NOT. inversible)THEN
	  PRINT*,'pivot nul dans gauss_band, ligne',ligne,', pivot',pivot,
	1 ', indpc',indpc(ligne)
	
	DO i=MAX(1,ligne-8),MIN(nl,ligne+8)
	 PRINT*,i, indpc(i)
	 WRITE(*,2000)a(i,:),b(1,i)
	ENDDO
	
	  RETURN
	 ENDIF
c	 PRINT*,'pivot pour ligne',ligne,pivot,ipivot

c permutation de ligne et de ipivot, division par pivot
	 DO j=1,ns
	  sto=b(j,ligne) ; b(j,ligne)=b(j,ipivot)/a(ipivot,1)
	  IF(ligne /= ipivot)b(j,ipivot)=sto
	 ENDDO
	 DO j=nemr,1,-1
	  sto=a(ligne,j) ; a(ligne,j)=a(ipivot,j)/a(ipivot,1)
	  IF(ligne /= ipivot)a(ipivot,j)=sto
	 ENDDO	!j

c d�calage pour avoir la diagonale en colonne 1
	 DO i=ligne+1,nl
	  IF(indpc(i) > ligne)CYCLE B1
c	  PRINT*,'i de perm',i
	  ai=a(i,1)
	  DO j=1,ns
	   b(j,i)=b(j,i)-b(j,ligne)*ai
	  ENDDO 
	  DO j=2,nemr
	   a(i,j-1)=a(i,j)-a(ligne,j)*ai
	  ENDDO	!j
	  a(i,nemr)=0.d0
	 ENDDO

c	 DO i=1,nl
c	  WRITE(*,1000)(a(i,j),j=1,nemr),(b(j,i),j=1,ns)
c1000	  FORMAT((1x,1p16e8.1))
c	 ENDDO	!i

	ENDDO B1

c on remonte de n (seulement) a 1, pour les ns seconds membres
	DO k=1,ns
	 DO i=n-1,1,-1
	  DO j=2,nemr
	   IF(i+j-1 <= n)b(k,i)=b(k,i)-a(i,j)*b(k,i+j-1)
	  ENDDO	!j
	 ENDDO	!i
	ENDDO	!k

	RETURN

	END SUBROUTINE gauss_band
