
c**************************************************************************

	SUBROUTINE pp3(t,ro,comp,dcomp,jac,deriv,fait,
	1 epsilon,et,ero,ex,hhe,be7e,b8e,n13e,o15e,f17e)

c	routine private du module mod_nuc

c	CYCLE PP cf. Clayton p. 380, 392 et 430
c	pp3 a pour vocation de servir pour les tests de mise au point

c	�l�ments pris en compte:
c	H1, He3, He4, Ex
c	Ex est l'�l�ment fictif compl�ment, il n'int�resse que la diffusion
c	H2, Li7, Be7 � l'�quilibre

c	au premier appel une table de r�actions est cr�e automatiquement
c	par rq_reac --> tabul_nuc

c	Auteur: P.Morel, D�partement J.D. Cassiop�e, O.C.A.
c	CESAM2k

c entr�es :
c	t : temp�rature cgs
c	ro : densit� cgs
c	comp : abondances
c	deriv=.true. : on calcule le jacobien
c	fait=1 : initialisation de la composition chimique
c	    =2 : calcul de dcomp et jacobien si deriv
c	    =3 : �nergie nucl�aire et d�riv�es / t et ro
c	    =4 : production de neutrinos

c sorties
c	dcomp : d�riv�e temporelle (unit� de temps : 10**6 ans)
c	jac : jacobien (unit� de temps : 10**6 ans)
c	epsilon, et, ero, ex : �nergie thermonucl�aire (unit� de temps : s)
c			   : et d�riv�es /t, ro ,X
c	hhe, be7e, b8e, n13e, o15e, f17e : nombre de neutrinos g/s
c	hhe r�action : H1(p,e+ nu)H2
c	be7e r�action : Be7(e-,nu g)Li7
c	b8e r�action : B8(,e+ nu)Be8
c	n13e r�action : N13(,e+ nu)C13 mis � 0
c	o15e r�action : O15(e+,nu)N15 mis � 0 
c	f17e r�action : F17(,e+ nu)O17 mis � 0

c initialisation de
c	ab_min : abondances n�gligeables
c	ab_ini : abondances initiales

c	r(1) : r�action H1(p,e+ nu)H2			PP
c	r(2) : r�action H2(p,g)H3
c	r(3) : r�action He3(He3,2p)He4
c	r(4) : r�action He3(a,g)Be7
c	r(5) : r�action Li7(p,a)He4
c	r(6) : r�action Be7(e-,nu g)Li7
c	r(7) : r�action Be7(p,g)B8(,e+ nu)Be8(a)He4

c	indices des �l�ments

c	H1 : 1
c	He3 : 2
c	He4 : 3
c	Ex : 4

c----------------------------------------------------------------------

	USE mod_donnees, ONLY : ab_ini, ab_min, ah, amu, fmin_abon, ihe4,
	1 i_ex, langue, nchim, nom_elem, nom_xheavy,
	2 nucleo, secon6, t_inf, x0, y0, zi, z0
	USE mod_kind
		
	IMPLICIT NONE
	
	INTEGER, INTENT(in) :: fait
	LOGICAL, INTENT(in) :: deriv
	REAL (kind=dp), INTENT(in):: t, ro
	REAL (kind=dp), INTENT(inout), DIMENSION(:) :: comp
	REAL (kind=dp), INTENT(out), DIMENSION(:,:) :: jac	
	REAL (kind=dp), INTENT(out), DIMENSION(:) :: dcomp, ex, epsilon
	REAL (kind=dp), INTENT(out) :: et, ero, hhe, be7e, b8e, n13e,
	1 o15e, f17e
	
	REAL (kind=dp), ALLOCATABLE, SAVE, DIMENSION(:,:) :: drx, dqx
	REAL (kind=dp), ALLOCATABLE, SAVE, DIMENSION(:) :: anuc, comp_dex,
	1 dmuex, dh2x, denx, dbe7x, dli7x, drt, dro, r, q, dqt, dqo		
	REAL (kind=dp) :: mue, nbz, den, be7, h2, li7, dh2t, dh2ro,
	1 dent, denro, dbe7t, dbe7ro, dli7t, dli7ro,
	2 mass_ex, charge_ex

	INTEGER :: i, j
	
	CHARACTER (len=2) :: text
		
c--------------------------------------------------------------------------

2000	FORMAT(8es10.3)
2001	FORMAT(5es15.8)
2002	FORMAT(11es8.1)
	
c	initialisations

	SELECT CASE(fait)
	CASE(0)
	 
c	 d�finition de nchim: nombre d'�l�ments chimiques dont on
c	 calcule l'abondance H1, He3, He4, C13, C13, N14, N15, O16, O17, Ex

	 nchim=3+1

c	 appel d'initialisation pour tabulation des r�actions nucl�aires
c	 allocations fictives

	 ALLOCATE(drx(1,1),dqx(1,1),r(1),drt(1),dro(1),q(1),
	1 dqt(1),dqo(1),dmuex(1))
	 CALL rq_reac(comp,1.d7,1.d0,r,drt,dro,drx,q,dqt,dqo,dqx,mue,dmuex) 
c	 PAUSE'apr�s case0'
	 
	 DEALLOCATE(dqx,drx) ; ALLOCATE(dqx(nreac,nchim),drx(nreac,nchim))
	 
	CASE(1)

c	 d�termination des abondances initiales
c	 He3+He4=Y0
c	 Z0 = somme des �l�ments plus lourds que helium
c	 dans Z rapports en nombre

	 CALL abon_ini
	 
c	 Ex : �l�ment fictif compl�mentaire
c	 pour pp3 il correspond � Z, on utilise O16

	 charge_ex=8.d0 ; mass_ex=16.d0
	 WRITE(text,10)NINT(mass_ex)
10	 FORMAT(i2)

	 nucleo(nchim)=mass_ex	!nucleo de l'�l�ment chimique reliquat
	 zi(nchim)=charge_ex	!charge de l'�l�ment chimique reliquat
	 i=NINT(charge_ex)
	 nom_elem(nchim)=elem(i)//text	!nom du reliquat
	 nom_xheavy=nom_elem(nchim)	!idem
	 i_ex=nchim 	!indice de l'�l�ment chimique reliquat	 
	 SELECT CASE(langue)	  
	 CASE('english')	
	  WRITE(*,1023)TRIM(nom_elem(nchim)),NINT(mass_ex),NINT(charge_ex)
	  WRITE(2,1023)TRIM(nom_elem(nchim)),NINT(mass_ex),NINT(charge_ex)	 
1023	  FORMAT(a,': fictitious species /= CNO, of mass : ',i3,/,
	1 'and charge :',i3)	 
	 CASE DEFAULT	 
	  WRITE(*,23)TRIM(nom_elem(nchim)),NINT(mass_ex),NINT(charge_ex)
	  WRITE(2,23)TRIM(nom_elem(nchim)),NINT(mass_ex),NINT(charge_ex)	 
23	  FORMAT(a,': �l�ment fictif de masse : ',i3,' et de charge :',i3)
	 END SELECT
	 	 
c	 PRINT*,nchim
c	 WRITE(*,2000)nucleo(1:nchim) 
	 
c d�termination des abondances initiales

	 comp(1)=x0/nucleo(1)				!H1
	 comp(3)=y0/(he3she4z*nucleo(2)+nucleo(3))	!He4
	 comp(2)=comp(3)*he3she4z			!He3
	 comp(4)=z0/nucleo(4)				!Ex

c allocations diverses

	 DEALLOCATE(drt,dro,r,q,dqt,dqo,dmuex)
	 ALLOCATE(ab_ini(nchim),ab_min(nchim),drt(nreac),dro(nreac),
	1 r(nreac),q(nreac),dqt(nreac),dqo(nreac),anuc(nchim),
	2 dmuex(nchim),dh2x(nchim),denx(nchim),dbe7x(nchim),dli7x(nchim))

c abondances initiales et abondances n�gligeables
	 
	 ab_ini(1:nchim)=comp(1:nchim)*nucleo(1:nchim)
	 ab_min=ab_ini*fmin_abon	 
	 
c nombre/volume des m�taux dans Z
		
	 nbz=SUM(comp(ihe4+1:nchim))	 
	 
c abondances en DeX, H=12

	 ALLOCATE(comp_dex(nchim))
	 comp_dex=12.d0+LOG10(comp/comp(1))
	 
c �critures
 
	 SELECT CASE(langue)	  
	 CASE('english')	
	  WRITE(2,1002) ; WRITE(*,1002) 
1002	  FORMAT(/,'PP thermonuclear reactions',/)
	  WRITE(2,1003)nreac ; WRITE(*,1003)nreac 
1003	  FORMAT('number of reaction : ',i3)
	  WRITE(2,1004)nchim ; WRITE(*,1004)nchim
1004	  FORMAT('number of species : ',i3)
	  WRITE(2,1020)x0,y0,z0,z0/x0 ; WRITE(*,1020)x0,y0,z0,z0/x0
1020	  FORMAT(/,'Initial abundances/mass computed with :',/,
	1 'X0=',es10.3,', Y0=',es10.3,', Z0=',es10.3,/,'Z0/X0=',es10.3,/,
	2 'H1=X0, H2+He3+He4=Y0, with H2 in He3',/,
	3 'Z0 = 1-X0-Y0 = Ex',/)	
	  WRITE(2,1)ab_ini(1:nchim) ; WRITE(*,1)ab_ini(1:nchim)
1	  FORMAT('H1 :',es10.3,', He3 :',es10.3,', He4 :',es10.3,
	1 ', Ex :',es10.3)
	  WRITE(2,1009)comp_dex ; WRITE(*,1009)comp_dex
1009	  FORMAT(/,'Initial abundances/number: 12+Log10(Ni/Nh)',/,
	1 'H1:',es10.3,', He3:',es10.3,', He4:',es10.3,', Ex :',es10.3)
	  WRITE(2,1021)comp(4)/nbz
1021	  FORMAT(/,'mass ratio by number within Z : Ex/Z :',es10.3)	
	  WRITE(2,1022)ab_ini(4)/z0
1022	  FORMAT(/,'mass ratio by mass within Z : Ex/Z :',es10.3)	
	  WRITE(2,1014)he3she4z,c13sc12,n15sn14,o17so16
1014	  FORMAT(/,'Isotopic ratios by number : He3/He4=',es10.3)	
	  WRITE(2,1005)ab_min(1:nchim) ; WRITE(*,1005)ab_min(1:nchim)
1005	  FORMAT(/,'threhold for neglectable abundances/mass :',/,
	1 'H1:',es10.3,', He3:',es10.3,', He4:',es10.3,', Ex:',es10.3)
	  WRITE(2,1006) ; WRITE(*,1006)
1006	  FORMAT(/,'H2, Li7, Be7 at equilibrium')
	  WRITE(2,1007) ; WRITE(*,1007)
1007	  FORMAT('Use of a table')
	  WRITE(2,1008) ; WRITE(*,1008)
1008	  FORMAT('Temporal evolution, test of precision on H1 He4')
	 CASE DEFAULT
	  WRITE(2,2) ; WRITE(*,2) 
2	  FORMAT(/,'R�actions thermonucl�aires des cycles PP',/)
	  WRITE(2,3)nreac ; WRITE(*,3)nreac 
3	  FORMAT('nombre de r�actions : ',i3)
	  WRITE(2,4)nchim ; WRITE(*,4)nchim
4	  FORMAT('nombre d''�l�ments chimiques : ',i3)
	  WRITE(2,20)x0,y0,z0,z0/x0 ; WRITE(*,20)x0,y0,z0,z0/x0
20	  FORMAT(/,'abondances initiales/gramme d�duites de:',/,
	1 'X0=',es10.3,', Y0=',es10.3,', Z0=',es10.3,/,'Z0/X0=',es10.3,/,
	2 'H1=X0, H2+He3+He4=Y0, avec H2 dans He3',/,
	3 'Z0 = 1-X0-Y0 = Ex',/)	
	  WRITE(2,1)ab_ini(1:nchim) ; WRITE(*,1)ab_ini(1:nchim)
	  WRITE(2,9)comp_dex ; WRITE(*,9)comp_dex
9	  FORMAT(/,'Abondances initiales en nombre: 12+Log10(Ni/Nh)',/,
	1 'H1:',es10.3,', He3:',es10.3,', He4:',es10.3,', Ex:',es10.3)
	  WRITE(2,21)comp(4)/nbz
	  WRITE(*,21)comp(4)/nbz
21	  FORMAT(/,'rapports en nombre dans Z: Ex/Z:',es10.3)	
	  WRITE(2,22)ab_ini(4)/z0
	  WRITE(*,22)ab_ini(4)/z0
22	  FORMAT(/,'rapports en masse dans Z: Ex/Z:',es10.3)	
	  WRITE(2,14)he3she4z
	  WRITE(*,14)he3she4z
14	  FORMAT(/,'Rapports isotopiques en nombre: He3/He4=',es10.3)	
	  WRITE(2,5)ab_min(1:nchim) ; WRITE(*,5)ab_min(1:nchim)
5	  FORMAT(/,'abondances/gramme n�gligeables:',/,
	1 'H1:',es10.3,', He3:',es10.3,', He4:',es10.3,', Ex:',es10.3)
	  WRITE(2,6) ; WRITE(*,6)
6	  FORMAT(/,'H2, Li7, Be7 � l''�quilibre')
	  WRITE(2,7) ; WRITE(*,7)
7	  FORMAT(/,'on utilise une table')
	  WRITE(2,8) ; WRITE(*,8)
8	  FORMAT(/,'�vol. temporelle, test de pr�cision sur H1 et He4')
	 END SELECT

c d�finitions diverses

	 ab_min(1:nchim)=ab_min(1:nchim)/nucleo(1:nchim)
	 anuc(1:nchim)=ANINT(nucleo(1:nchim))		!nombre atomique

c nettoyage

	 DEALLOCATE(comp_dex)
	 	 
c les r�actions	 
	 
	CASE(2)
	 dcomp=0.d0 ; jac=0.d0

	 IF(t < t_inf)RETURN	
	
	 CALL rq_reac(comp,t,ro,r,drt,dro,drx,q,dqt,dqo,dqx,mue,dmuex)	

c	 WRITE(*,*)'comp' ; WRITE(*,2000)comp(1:nchim)
c	 WRITE(*,*)'r�actions' ; WRITE(*,2000)r(1:nreac)

c	 �quations d'�volution

	 dcomp(1)=-3.d0*r(1)*comp(1)**2
	1 +(2.d0*r(3)*comp(2)-r(4)*comp(3))*comp(2)		!H1
	 dcomp(2)=r(1)*comp(1)**2-(2.d0*r(3)*comp(2)
	1 +r(4)*comp(3))*comp(2) 				!He3
	 dcomp(3)=(r(3)*comp(2)+r(4)*comp(3))*comp(2)	!He4

c	   Pour v�rifications SUM dcomp*nucleo=0

c	 PRINT*,'pp3, v�rifications SUM dcomp*nucleo=0'
c	 WRITE(*,2000)DOT_PRODUCT(dcomp,anuc) ; PAUSE'v�rif' 

c	conservation des baryons
	 
	 dcomp(nchim)=-DOT_PRODUCT(dcomp,anuc)/anuc(nchim)
	
c	 WRITE(*,*)'somme << dcomp / dcomp'
c	 WRITE(*,2000)SUM(nucleo(1:nchim)*dcomp(1:nchim))
c	 WRITE(*,2000)dcomp(1:nchim) ; WRITE(*,*)'reac,rx'
c	 WRITE(*,2000)r(1:nreac) ; WRITE(*,2000)rx(1:nreac)
c	 PAUSE'apres case2 avant deriv'
	 
c	 calcul du jacobien

	 IF(deriv)THEN	!jac(i,j) : �quation, i : �l�ment j

c	 dcomp(1)=-3.d0*r(1)*comp(1)**2
c	1 +(2.d0*r(3)*comp(2)-r(4)*comp(3))*comp(2)		!H1				!He3

	  jac(1,1)=-6.d0*r(1)*comp(1)			!d /H1
	  jac(1,2)=4.d0*r(3)*comp(2)-r(4)*comp(3)	!d /He3
	  jac(1,3)=-r(4)*comp(2)			!d /He4
	 
	  DO i=1,nchim		!d�pendances dues � l'effet d'�cran
	   jac(1,i)=jac(1,i)-3.d0*drx(1,i)*comp(1)**2
	2  +(2.d0*drx(3,i)*comp(2)-drx(4,i)*comp(3))*comp(2)
	  ENDDO
			 
c	  �quation dcomp(2)
c	  dcomp(2)=r(1)*comp(1)**2-(2.*r(3)*comp(2)
c	  +r(4)*comp(3))*comp(2)!He3

	  jac(2,1)=2.d0*r(1)*comp(1)			!d /H1
	  jac(2,2)=-4.d0*r(3)*comp(2)-r(4)*comp(3)	!d /He3
	  jac(2,3)=-r(4)*comp(2)			!d /He4

	  DO i=1,nchim		!d�pendances dues � l'effet d'�cran
	   jac(2,i)=jac(2,i)
	1  +drx(1,i)*comp(1)**2-(2.d0*drx(3,i)*comp(2)
	2  +drx(4,i)*comp(3))*comp(2)
	  ENDDO
	 
c	  �quation dcomp(3)
c	  dcomp(3)=(r(3)*comp(2)+r(4)*comp(3))*comp(2)

	  jac(3,2)=2.d0*r(3)*comp(2)+r(4)*comp(3)	!d /He3
	  jac(3,3)=r(4)*comp(2)				!d /He4
	 
	  DO i=1,nchim		!d�pendances dues � l'effet d'�cran
	   jac(3,i)=jac(3,i)
	1  +(drx(3,i)*comp(2)+drx(4,i)*comp(3))*comp(2)
	  ENDDO
	 
c 	  �quation dcomp(4)	 
c	  DO i=1,nchim-1				!Ex
c	   dcomp(nchim)=dcomp(nchim)+anuc(i)*dcomp(i)
c	  ENDDO
c	  dcomp(nchim)=-dcomp(nchim)/anuc(nchim) !conservation des baryons

	  DO j=1,nchim
	   DO i=1,nchim-1
	    jac(nchim,j)=jac(nchim,j)+anuc(i)*jac(i,j)
	   ENDDO
	   jac(nchim,j)=-jac(nchim,j)/anuc(nchim)
	  ENDDO

c	  unit�s de temps pour int�gration temporelle

	  jac=jac*secon6
	  
c	  PAUSE'apr�s case2 deriv'
	 ENDIF		!deriv

	 dcomp=dcomp*secon6

	CASE(3)

c	 calcul de la production d'�nergie nucl�aire et d�riv�es
c	 pour H2(H,g)He3, q(2)H**2=q(2)*r(1)/r(2)
	 
	 epsilon(1:4)=0.d0 ; et=0.d0 ; ero=0.d0 ; ex=0.d0
	 IF(t <= t_inf)RETURN
	
	 CALL rq_reac(comp,t,ro,r,drt,dro,drx,q,dqt,dqo,dqx,mue,dmuex)

c	 mue : nombre d'�lectrons / mole /g = 1/poids mol. moy. par e-

	 IF(comp(1) > 0.d0)THEN
	  h2=r(1)/r(2)*comp(1) ; den=r(6)*mue+r(7)*comp(1)
	  be7=r(4)*comp(2)*comp(3)/den ; li7=r(6)*be7*mue/r(5)/comp(1)
	 ELSE
	  h2=0.d0 ; be7=0.d0 ; li7=0.d0
	 ENDIF
	
c	 PRINT*,'h2,li7,be7' ; WRITE(*,2000)h2,li7,be7

	 epsilon(2)=(q(1)*comp(1)+q(2)*h2+q(5)*li7+q(7)*be7)*comp(1)
	1 +(q(3)*comp(2)+q(4)*comp(3))*comp(2)+q(6)*mue*be7
	 DO i=2,4
	  epsilon(1)=epsilon(1)+epsilon(i)
	 ENDDO
c	 PAUSE'apres case3 avant deriv'
	 
	 IF(deriv)THEN	
	  IF(h2 > 0.d0)THEN
	   dh2t=h2*(drt(1)/r(1)-drt(2)/r(2))
	   dh2ro=h2*(dro(1)/r(1)-dro(2)/r(2))
	   DO i=1,nchim
	    dh2x(i)=h2*(drx(1,i)/r(1)-drx(2,i)/r(2))
	   ENDDO
	   dh2x(1)=dh2x(1)-h2/comp(1)	  
	  ELSE
	   dh2t=0.d0 ; dh2ro=0.d0 ; dh2x=0.d0
	  ENDIF

	  IF(be7 > 0.d0)THEN
	   dent= drt(6)*mue+drt(7)*comp(1)
	   denro=dro(6)*mue+dro(7)*comp(1)
	   DO i=1,nchim	  
	    denx(i)=drx(6,i)*mue+r(6)*dmuex(i)+drx(7,i)*comp(1)
	   ENDDO
	   denx(1)=denx(1)+r(7)
	   	   
	   dbe7t= be7*(drt(4)/r(4)- dent/den)
	   dbe7ro=be7*(dro(4)/r(4)-denro/den)
	   DO i=1,nchim
	    dbe7x(i)=be7*(drx(4,i)/r(4)-denx(i)/den)
	   ENDDO
	   dbe7x(2)=dbe7x(2)+be7/comp(2)
	   dbe7x(3)=dbe7x(3)+be7/comp(3)	   
	  ELSE
	   dbe7t=0.d0 ; dbe7ro=0.d0 ; dbe7x=0.d0
	  ENDIF

	  IF(li7 > 0.d0)THEN
	   dli7t= li7*(drt(6)/r(6) +dbe7t/be7-drt(5)/r(5))
	   dli7ro=li7*(dro(6)/r(6)+dbe7ro/be7-dro(5)/r(5))
	   DO i=1,nchim	   
	    dli7x(i)=li7*(drx(6,i)/r(6)+dbe7x(i)/be7
	1	+dmuex(i)/mue-drx(5,i)/r(5))
	   ENDDO
	   dli7x(1)=dli7x(1)-li7/comp(1)	   
	  ELSE
	   dli7t=0.d0 ; dli7ro=0.d0 ; dli7x=0.d0
	  ENDIF
	  	
 	  et=(dqt(1)*comp(1)+dqt(2)*h2+dqt(5)*li7+dqt(7)*be7)*comp(1)
	1  +(dqt(3)*comp(2)+dqt(4)*comp(3))*comp(2)
	2  +dqt(6)*mue*be7
	5  +(q(2)*dh2t+q(5)*dli7t+q(7)*dbe7t)*comp(1)+q(6)*mue*dbe7t
	
 	  ero=(dqo(1)*comp(1)+dqo(2)*h2+dqo(5)*li7+dqo(7)*be7)*comp(1)
	1  +(dqo(3)*comp(2)+dqo(4)*comp(3))*comp(2)
	2  +dqo(6)*mue*be7
	5  +(q(2)*dh2ro+q(5)*dli7ro+q(7)*dbe7ro)*comp(1)+q(6)*mue*dbe7ro

	  ex(1)=2.d0*q(1)*comp(1)+q(2)*h2+q(5)*li7+q(7)*be7
	  ex(2)=2.d0*q(3)*comp(2)+q(4)*comp(3) ; ex(3)=q(4)*comp(2)
	 
	  DO i=1,nchim	!contributions des �crans
	   ex(i)=ex(i)+(dqx(1,i)*comp(1)+dqx(2,i)*h2
	1   +q(2)*dh2x(i)+dqx(5,i)*li7+q(5)*dli7x(i)
	2   +dqx(7,i)*be7+q(7)*dbe7x(i))*comp(1)
	3   +(dqx(3,i)*comp(2)+dqx(4,i)*comp(3))*comp(2)
	4   +dqx(6,i)*mue*be7+q(6)*dmuex(i)*be7+q(6)*mue*dbe7x(i)
	  ENDDO
c	  PAUSE'apres case3 deriv'
	 
	 ENDIF	!deriv
	 
	CASE(4)		!taux de production des neutrinos

	 IF(t >= t_inf)THEN
	  CALL rq_reac(comp,t,ro,r,drt,dro,drx,q,dqt,dqo,dqx,mue,dmuex)	 
	  be7=r(4)*comp(2)*comp(3)/(r(6)*mue+r(7)*comp(1))
	  hhe=r(1)*comp(1)**2/amu ; be7e=r(6)*mue*be7/amu
	  b8e=r(7)*comp(1)*be7/amu ; n13e=0.d0 ; o15e=0.d0 ; f17e=0.d0
	 ELSE
	  hhe=0.d0 ; be7e=0.d0 ; b8e=0.d0
	  n13e=0.d0 ; o15e=0.d0 ; f17e=0.d0
	 ENDIF
c	 PAUSE'apres case4'

	CASE DEFAULT
	 PRINT*,'ppcno0, fait ne peut prendre que les valeurs 1, 2, 3 ou 4'
	 PRINT*,'ERREUR: fait a la valeur:',fait
	 PRINT*,'ARRET' ; PRINT* ; STOP
	 
	END SELECT
	
	RETURN

	END SUBROUTINE pp3
