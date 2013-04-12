
c*******************************************************************

	SUBROUTINE newspl(n,x,x1t,kno1,m1,x2t,kno2,m2,s1,s2,duale)

c subroutine public du module mod_numerique

c pour n fonctions
c transforme la spline s1 d'ordre m1 sur le r�seau x1t \ kno1,
c en la spline s2 d'ordre m2 sur le r�seau x2t \ kno2
c on n'utilise pas la base duale

c le cas ou la spline s2 a des discontinuit�s est envisag�:
c la discontinuit� ne peut �tre localis�e avec la pr�cision
c machine, on se d�place l�g�rement � gauche et � droite

c Auteur: P. Morel, D�partement J.D. Cassini, O.C.A.
c CESAM2k (version f95 de newspl1.f)

c entr�es
c	n: nombre de fonctions
c	x, x1t, x2t, kno1, kno2: abcisses et vecteurs nodaux
c	m1, m2: ordres des splines
c	s1: premi�re spline
c	duale=.TRUE. : on utilise la base duale

c sorties:
c	s2: nouvelle spline

c--------------------------------------------------------------------------

	USE mod_kind
		
	IMPLICIT NONE

	REAL (kind=dp), INTENT(in), DIMENSION (:) :: x, x2t
	INTEGER, INTENT(in) :: kno2, m1, m2, n
	LOGICAL, INTENT(in), OPTIONAL :: duale
	REAL (kind=dp), INTENT(inout), DIMENSION (:,:) :: s1
	REAL (kind=dp), INTENT(inout), DIMENSION (:) :: x1t
	INTEGER, INTENT(inout) :: kno1
	REAL (kind=dp), INTENT(out), DIMENSION (:,:) :: s2
	
	REAL (kind=dp), DIMENSION (kno2-m2,m2) :: a
	REAL (kind=dp), DIMENSION (MAX(m1,m2)+1) :: q
	REAL (kind=dp), DIMENSION (n) :: dfxdx, fx	
	
	REAL (kind=dp), PARAMETER :: dx=1.d-3 , un_dx=1.d0-dx
	REAL (kind=dp) :: p
			
	INTEGER, DIMENSION(kno2-m2) :: indpc
		
	INTEGER :: i, l=1, n1=0, n2
	LOGICAL :: inversible

c------------------------------------------------------------------------

2000	FORMAT(8es10.3)
	
c dimension de la spline 2
	n2=kno2-m2
		
c calcul des tau selon de-Boor p.123

	a=0.d0 ; i=0
	B1: DO
	 i=i+1 ; IF(i > n2)EXIT B1
	 p=SUM(x2t(i+1:i+m2-1))/REAL(m2-1,dp)	!abscisses de de Boor
	 
c d�tection d'une discontinuit� pour s2	 
	 IF(i > 1 .AND. i < n2 .AND. p == x2t(i+1)) THEN
	 	  
c on se place � gauche de la discontinuit� 	  
	  p=x2t(i+1)*un_dx+x2t(i)*dx	  	  
	  CALL bsp1dn(n,s1,x,x1t,n1,m1,kno1,.TRUE.,p,l,fx,dfxdx)
	  IF(no_croiss)PRINT*,'Pb. en 1 dans newspl'	  	  
	  s2(:,i)=fx
	  CALL linf(p,x2t,kno2,l) ; CALL bval0(p,x2t,m2,l,q)
	  a(i,1:m2)=q(1:m2) ; indpc(i)=l-m2+1
c	  PRINT*,'� gauche, i=',i
c	  PRINT*,x2t(i),x2t(i+1),p
	  
c on se place � droite de la discontinuit�
	  i=i+1
	  p=x2t(i)*un_dx+x2t(i+m2)*dx	  	  
	  CALL bsp1dn(n,s1,x,x1t,n1,m1,kno1,.TRUE.,p,l,fx,dfxdx)
	  IF(no_croiss)PRINT*,'Pb. en 2 dans newspl'	  	  
	  s2(:,i)=fx
	  CALL linf(p,x2t,kno2,l) ; CALL bval0(p,x2t,m2,l,q)
	  a(i,1:m2)=q(1:m2) ; indpc(i)=l-m2+1
c	  PRINT*,'� droite, i=',i
c	  PRINT*,x2t(i),x2t(i+m2),p
	  	  
c sans discontinuit�	  	  
	 ELSE	  
c	  PRINT*,i,p,x1t(1),x2t(1)
	  CALL bsp1dn(n,s1,x,x1t,n1,m1,kno1,.TRUE.,
	1 MAX(x1t(1),MIN(p,x1t(kno1))),l,fx,dfxdx)
	  IF(no_croiss)PRINT*,'Pb. en 3 dans newspl'
	  s2(:,i)=fx
	  CALL linf(p,x2t,kno2,l) ; CALL bval0(p,x2t,m2,l,q)
	  a(i,1:m2)=q(1:m2) ; indpc(i)=l-m2+1	  
	 ENDIF 
	ENDDO	B1

c utilisation de la base duale	
	IF(PRESENT(duale))THEN
	 IF(duale)RETURN
	ENDIF

c coefficients de s2
	CALL gauss_band(a,s2,indpc,n2,n2,m2,n,inversible)
	
c si s2 est pas inversible on utilise la base duale	
	IF(.NOT.inversible)THEN
	 PRINT*,'matrice non inversible dans newspl'
	 PRINT*,'tentative d''utilisation de la base duale'	
	 DO i=1,n2
	  CALL bsp1dn(n,s1,x,x1t,n1,m1,kno1,.TRUE.,x2t(i),l,fx,dfxdx)
	  s2(:,i)=fx	 
	 ENDDO	 
	ENDIF
		
	RETURN

	END SUBROUTINE newspl
