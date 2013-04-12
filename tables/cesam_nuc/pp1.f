
c******************************************************************

	SUBROUTINE pp1(t,ro,comp,dcomp,jac,deriv,fait,
	1 epsilon,et,ero,ex,hhe,be7e,b8e,n13e,o15e,f17e)

c	routine private du module mod_nuc

c	r�actions thermonucl�aires selon GONG

c	Auteur: P. Morel, D�partement J.D. Cassini, O.C.A.
c	CESAM2k

c entree :
c	t : temp�rature cgs
c	ro : densit� cgs
c	comp : abondances
c	deriv=.TRUE. : on calcule les d�riv�es dont le jacobien
c	fait=1 : initialisation de la composition chimique
c	    =2 : d�riv�e premiere et jacobien si deriv
c	    =3 : �nergie nucl�aire et d�riv�es / t
c	    =4 : production de neutrinos

c sorties
c	dcomp : d�riv�e temporelle (unit� de temps : 10**6 ans)
c	jac : jacobien (unit� de temps : 10**6 ans)
c	e, et, ero, ex : �nergie thermonucl�aire (unit� de temps : s)
c	et d�riv�es /t, ro ,X
c	hhe, be7e, b8e, n13e, o15e, f17e : nombre de neutrinos g/s
c	pour les r�actions designees par le symbole (e pour electron +)

c----------------------------------------------------------------

	USE mod_donnees, ONLY : ab_ini, ab_min, ah, amu, fmin_abon, ihe4, i_ex,
	1 nchim, nom_elem, nucleo, secon6, t_inf, t_sup,
	2 x0, zi
	USE mod_kind

	IMPLICIT NONE
	
	REAL (kind=dp), INTENT(in):: t, ro		
	integer, INTENT(in) :: fait
	logical, INTENT(in) :: deriv	  
	REAL (kind=dp), INTENT(inout), DIMENSION(:) :: comp
	REAL (kind=dp), INTENT(out), DIMENSION(:,:) :: jac		 	
	REAL (kind=dp), INTENT(out), DIMENSION(:) :: dcomp, ex, epsilon
	REAL (kind=dp), INTENT(out) :: et, ero, hhe, be7e, b8e, n13e,
	1 o15e, f17e

	REAL (kind=dp), PARAMETER :: re=6.5d0, a11=4.21d-15, b=3.6d0,
	1 qe=6.5d-5, cte3=1.d-9, cte5=1.d0/3.d0
	REAL (kind=dp), SAVE :: cte1, cte6
	REAL (kind=dp) :: t913, t9
	
	INTEGER :: i

	LOGICAL, SAVE :: init=.TRUE.
	
c-----------------------------------------------------------------

2000	FORMAT(8es10.3)

	IF(init)THEN
	 init=.FALSE. ; cte1=-re*a11*secon6/2.d0*ah
	 cte6=qe*a11/2.d0/amu*ah**2 ; l_vent=.FALSE. ; l_planet=.FALSE.
	ENDIF

	SELECT CASE(fait)
	CASE(0)		!initialisation, nombre d'�l�ments
	
c d�finition de nchim: nombre d'�l�ments chimiques dont on
c calcule l'abondance: H1
	 nchim=1 ; nreac=1 ; t_inf=1.d6 ; t_sup=20.d6 ; ihe4=-100 ; i_ex=-100

	 ALLOCATE(ab_ini(nchim),ab_min(nchim),nom_elem(nchim),
	1 nucleo(nchim),zi(nchim))

	 nom_elem(1)=' H1 ' ; nucleo(1)=ah ; zi(1)=1.d0	 

	CASE(1)		!initialisation des abondances initiales
	 ab_ini(1)=x0 ; ab_min=ab_ini*fmin_abon ; comp(1)=x0/nucleo(1)
	
	 WRITE(2,*)'R�action thermonucl�aire du cycle PP simplifie (GONG)'
	 WRITE(2,*)' ' ; WRITE(2,*)'nombre de r�actions: ',nreac
	 WRITE(2,*)'nombre d''�l�ments chimiques: ',nchim
	 WRITE(2,*)' ' ; WRITE(2,3)comp(1:nchim)
3	 FORMAT(1x,'Abondance initiale H: ',es10.3)
2	 FORMAT(1x,'Abondance n�gligeable H:', es10.3)
	 WRITE(2,2)ab_min(1:nchim) ; WRITE(2,*)' '
	 WRITE(2,*)'pour l''�volution temporelle, test de pr�cision sur H'
	 WRITE(2,*)' '
	
	 WRITE(*,*)'R�action thermonucl�aire du cycle PP simplifi� (GONG)'
	 WRITE(*,*)' '
	 WRITE(*,*)'nombre de r�actions: ',nreac
	 WRITE(*,*)'nombre d''�l�ments chimiques: ',nchim
	 WRITE(*,*)' ' ; WRITE(*,3)comp(1:nchim)
	 WRITE(*,2)ab_min(1:nchim) ; WRITE(*,*)' '
	 WRITE(*,*)'pour l''�volution temporelle, test de pr�cision sur H'
	 WRITE(*,*)' '
	 
c par mole	 
	 ab_min=ab_min/nucleo
	 	 
	CASE(2)
	 IF(t < t_inf)THEN		!si t<t_inf
	  dcomp=0.d0 ; jac=0.d0 ; epsilon(1:4)=0.d0
	 ELSE
	  t913=(t*1.d-9)**(1.d0/3.d0)		!evolution chimique
	  dcomp(1)=cte1*comp(1)**2*ro/t913**2*exp(-b/t913)	!f=aX**2	
	  IF(deriv)THEN
	   jac=0.d0		!pour Mw jac=0
	   jac(1,1)=dcomp(1)*2.d0/comp(1)	!jac=d y'/dX=2*aX
	  ENDIF	 
	 ENDIF
	 
	CASE(3)
c	 PRINT*,'pp1 case3'
	 ex=0.d0
c	 WRITE(*,2000)t,ex,comp,ro
	 IF(t <= t_inf)THEN
	  epsilon(1:4)=0.d0 ; et=0.d0 ; ero=0.d0
	 ELSE		!�nergie nucl�aire
	  t9=t*cte3 ; t913=t9**cte5 ; epsilon=0.d0
	  epsilon(2)=cte6*comp(1)**2*ro/t913/t913*exp(-b/t913)	!4.1
	  DO i=2,4
	   epsilon(1)=epsilon(1)+epsilon(i)
	  ENDDO
	  IF(deriv)THEN
	   ero=epsilon(1)/ro ; et=epsilon(1)*((b/t913-2.d0)/3.d0/t9*cte3)
	   ex(1)=epsilon(1)*2.d0/comp(1)
	  ENDIF
	 ENDIF
	 
	CASE(4)			!neutrino
	 hhe=0.d0 ; be7e=0.d0 ; b8e=0.d0 ; n13e=0.d0
	 o15e=0.d0  ; f17e=0.d0 
	CASE DEFAULT
	 PRINT*,'erreur dans pp1, fait <=4, fait=', fait
	 PRINT*,'ARRET' ; STOP	 
	END SELECT
 
	END SUBROUTINE pp1
