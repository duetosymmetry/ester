
c******************************************************************

	SUBROUTINE bsp_gal(n,f,m,xt,knot,coeff)

c routine PUBLIC du module mod_numerique

c d�termination des coefficients splines pour interpolation de n fonctions
c connues en tout point. Le but est de d�terminer une approximation des
c d�riv�es qui ne sont pas accessibles, la base de spline est quelconque
c en particulier elle peut comporter des discontinuit�s de tous ordre.
c Exemple d'application d�termination de d ln ro / d lnR
c On utilise le formalisme des �l�ments finis pour former le syst�me lin�aire
c qui a n = nombre de fonctions seconds membres
c dont la solution donne les coefficients d'interpolation spline.
c Les coefficients sont exploitables par bsp1dn, bsp1ddn.
c l'int�grale de Gauss �tant d'ordre 2m-1 > m il n'y a pas de perte de pr�cision

c entr�es : 
c	n : nombre de fonctions
c	x : abscisses
c	nx : nombre d'abscisses
c	m : ordre des splines
c	xt : vecteur nodal
c	knot : dimension du vecteur nodal

c sorties :
c	f(n,knot-m) : coefficients des splines (knot-m dimension de la base)

c external :
c	coeff : routine du calcul des n valeurs des fonctions 

c Auteur: P.Morel, D�partement Cassiop�e, O.C.A., CESAM2k

c----------------------------------------------------------------

	USE mod_kind
		
	IMPLICIT NONE
	
	INTERFACE
	 SUBROUTINE coeff(n,xx,y)
	  USE mod_kind
	  IMPLICIT NONE	  
	  REAL(kind=dp), INTENT(in) :: xx	   
	  INTEGER, INTENT(in) :: n
	  REAL(kind=dp), INTENT(out), DIMENSION(n) :: y  
	 END SUBROUTINE coeff
	END INTERFACE
	
c mg : nb. points pour l'int�gration de Gauss d'ordre 2mg-1
	INTEGER, PARAMETER :: mg=2
	
	REAL(kind=dp), INTENT(in) , DIMENSION(:) :: xt
	INTEGER, INTENT(in) :: knot, m, n
	REAL(kind=dp), INTENT(out) , DIMENSION(:,:) :: f
			
	REAL(kind=dp), ALLOCATABLE, DIMENSION(:,:) :: a
	REAL(kind=dp), DIMENSION(m+1) :: q
	REAL(kind=dp), DIMENSION(mg) :: qg, wg
	REAL(kind=dp), DIMENSION(n) :: y
	
	INTEGER, ALLOCATABLE, DIMENSION(:) :: indpc
	INTEGER :: bloc, col, colonne, dim, i, ig, j, k, l=1, ligne, spi

	LOGICAL :: inversible
	
c----------------------------------------------------------------

2000	FORMAT(8es10.3)
2003	FORMAT(i3,9es8.1)

c initialisations
c bloc: longueur du bloc des coeff. non idt. nuls
	bloc=2*m-1

c dim : dimension de l'espace des splines 
c dim = nl = rang : nombre de lignes du syst�me lin�aire
	dim=knot-m

c allocations et initialisations		
	ALLOCATE(a(dim,bloc),indpc(dim))
	a=0.d0 ; f=0.d0
	 
c indices de premi�re colonne pour les produits scalaires
	DO i=1,dim      	!pour chaque spline de la base
	 indpc(i)=MAX(1,i-m+1)
	ENDDO 
	 
c---------------d�but de la construction du syt�me lin�aire------------	
	    
c les produits scalaires, pour chaque intervalle de m � dim
	B1: DO k=m,dim

c on �vite les points multiples
	 IF(xt(k) >= xt(k+1))CYCLE B1	 
	  	 
c spi : indice-1 de la premi�re spline non id.nulle  
c comme xt(l) <= qg(ig) < xt(l+1), l est le m�me pour toutes
c les abscisses de Gauss dans l'intervalle [xt(k),xt(k+1)]
	 CALL linf(xt(k),xt,knot,l) ; spi=l-m	  

c poids, wg et abscisses, qg pour int�gration Gauss d'ordre 2mg-1
	 CALL intgauss(xt(k),xt(k+1),qg,wg,mg)

c pour chaque abscisse de GAUSS
	 DO ig=1,mg

c variables
	  CALL coeff(n,qg(ig),y)
	   		   
c les B-splines en qg(ig), d�riv�es 0 � 1
	  CALL bval0(qg(ig),xt,m,l,q)   
	   	   
c contribution au syst�me lin�aire, avec n second membres
	  DO i=1,m		!pour chaque spline d'indice i
	   ligne=spi+i
	   f(:,ligne)=f(:,ligne)+wg(ig)*q(i)*y(:)	   
		
c la matrice compress�e est le jacobien 'diagonal' ie. sans les
c �l�ments 'non diagonaux' identiquement nuls
	   DO j=1,m		!pour chaque spline j
	    colonne=spi+j 
	    col=colonne-indpc(ligne)+1   !colonne matrice compress�e	     
	    a(ligne,col)=a(ligne,col)+wg(ig)*q(i)*q(j)
	   ENDDO	!j
	  ENDDO		!i	   
	 ENDDO 		!ig
	ENDDO B1
	 
c----------------fin de la construction du syst�me lin�aire-------------- 

c	 DO i=1,dim
c	  PRINT*,i,indpc(i)
c	  WRITE(*,2000)a(i,:),b(1,i)
c	 ENDDO
c	 PAUSE'a'
	 	
c r�solution du syst�me lin�aire
	CALL gauss_band(a,f,indpc,dim,dim,bloc,n,inversible)	
	IF(.NOT.inversible)THEN
	 PRINT*,'ARRRET, matrice singuli�re dans bsp_gal' ; STOP
	ENDIF	 
c	WRITE(*,2000)f(1,:) ; PAUSE'solution'

c sortie	
	DEALLOCATE(a,indpc)
	
	RETURN
	
	END SUBROUTINE bsp_gal
