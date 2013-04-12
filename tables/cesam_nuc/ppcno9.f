
!**************************************************************************

	SUBROUTINE ppcno9(t,ro,comp,dcomp,jac,deriv,fait,
	1 epsilon,et,ero,ex,hhe,be7e,b8e,n13e,o15e,f17e)

!	routine private du module mod_nuc

!	CYCLEs PP, CNO 
!	cf. Clayton p. 380, 392 et 430

!	�l�ments pris en compte:
!	H1, He3, He4, C12, C13, N14, N15, O16, O17, Ex
!	Ex est l'�l�ment fictif compl�ment, il n'int�resse que la diffusion
!	H2, Li7, Be7 � l'�quilibre

!	au premier appel une table ppcno est cr�e automatiquement par
!	rq_reac --> tabul_nuc

!	Auteur: P. Morel, D�partement J.D. Cassini, O.C.A.
!	CESAM2k

! entr�es :
!	t : temp�rature cgs
!	ro : densit� cgs
!	deriv=.true. : on calcule le jacobien
!	fait=1 : initialisation de la composition chimique
!	    =2 : calcul de dcomp et jacobien si deriv
!	    =3 : �nergie nucl�aire et d�riv�es / t et ro
!	    =4 : production de neutrinos

! entr�es/sorties :
!	comp : abondances par mole

! sorties
!	dcomp : d�riv�e temporelle (unit� de temps : 10**6 ans)
!	jac : jacobien (unit� de temps : 10**6 ans)
!	epsilon, et, ero, ex : �nergie thermonucl�aire (unit� de temps : s)
!			   : et d�riv�es /t, ro ,X
!	hhe, be7e, b8e, n13e, o15e, f17e : nombre de neutrinos g/s
!	hhe r�action : H1(p,e+ nu)H2
!	be7e r�action : Be7(e-,nu g)Li7
!	b8e r�action : B8(,e+ nu)Be8
!	n13e r�action : N13(,e+ nu)C13
!	o15e r�action : O15(e+,nu)N15 
!	f17e r�action : F17(,e+ nu)O17

! initialisation de
!	ab_min : abondances n�gligeables
!	ab_ini : abondances initiales

!	r(1) : r�action H1(p,e+ nu)H2			PP
!	r(2) : r�action H2(p,g)H3
!	r(3) : r�action He3(He3,2p)He4
!	r(4) : r�action He3(a,g)Be7
!	r(5) : r�action Li7(p,a)He4
!	r(6) : r�action Be7(e-,nu g)Li7
!	r(7) : r�action Be7(p,g)B8(,e+ nu)Be8(a)He4

!	r(8) : r�action C12(p,g)N13(,e+ nu)C13		CNO
!	r(9) : r�action C13(p,g)N14
!	r(10) : r�action N14(p,g)O15(e+,nu)N15
!	r(11) : r�action N15(p,g)O16
!	r(12) : r�action N15(p,a)C12
!	r(13) : r�action O16(p,g)F17(,e+ nu)O17
!	r(14) : r�action O17(p,a)N14

!	indices des �l�ments

!	H1 : 1
!	He3 : 2
!	He4 : 3
!	C12 : 4
!	C13 : 5
!	N14 : 6
!	N15 : 7
!	O16 : 8
!	O17 : 9
!	Ex : 10

!----------------------------------------------------------------------

	USE mod_donnees, ONLY : ab_ini, ab_min, ah, amu, fmin_abon, ihe4,
	1 i_ex, langue, nchim, nom_elem, nom_xheavy,
	2 nucleo, secon6, t_inf, x0, y0, zi, z0
	USE mod_kind
	USE mod_numerique, ONLY : gauss_band
		
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
	REAL (kind=dp), ALLOCATABLE, DIMENSION(:,:) :: a, b
	REAL (kind=dp), ALLOCATABLE, SAVE, DIMENSION(:) :: anuc, comp_dex,
	1 dmuex, dh2x, denx, dbe7x, dli7x, drt, dro, r, q, dqt, dqo		
	REAL (kind=dp) :: mue, nbz, den, be7, h2, li7, dh2t, dh2ro,
	1 dent, denro, dbe7t, dbe7ro, dli7t, dli7ro,
	2 mass_ex, charge_ex, sum_a
		
	INTEGER, ALLOCATABLE, DIMENSION(:) :: indpc
	INTEGER :: i, j
	
	LOGICAL :: inversible
	
	CHARACTER (len=2) :: text
		
!--------------------------------------------------------------------------

2000	FORMAT(8es10.3)
2001	FORMAT(5es15.8)
2002	FORMAT(11es8.1)
	
! initialisations
	SELECT CASE(fait)
	CASE(0)
	 
! d�finition de nchim: nombre d'�l�ments chimiques dont on
! calcule l'abondance H1, He3, He4, C13, C13, N14, N15, O16, O17, Ex
	 nchim=9+1

! appel d'initialisation pour tabulation des r�actions nucl�aires
! allocations fictives
	 ALLOCATE(drx(1,1),dqx(1,1),r(1),drt(1),dro(1),q(1),
	1 dqt(1),dqo(1),dmuex(1))
	 CALL rq_reac(comp,1.d7,1.d0,r,drt,dro,drx,q,dqt,dqo,dqx,mue,dmuex) 
!	 PAUSE'apr�s case0'
	 
	 DEALLOCATE(dqx,drx) ; ALLOCATE(dqx(nreac,nchim),drx(nreac,nchim))
	 
	CASE(1)

! d�termination des abondances initiales, He3+He4=Y0
! Z0 = somme des �l�ments plus lourds que helium, dans Z rapports en nombre
	 CALL abon_ini
	 
! Ex : �l�ment fictif moyenne des �l�ments # CNO
	 charge_ex=0.d0 ; mass_ex=0.d0 ; sum_a=0.d0
	 B1: DO i=3,nelem_ini		!� partir de Li=3
	  IF(elem(i) == ' C')CYCLE B1
	  IF(elem(i) == ' N')CYCLE B1
	  IF(elem(i) == ' O')CYCLE B1
	  charge_ex=charge_ex+c(i)*ab(i) ; mass_ex=mass_ex+m(i)*ab(i)
	  sum_a=sum_a+ab(i)
	 ENDDO B1
	 charge_ex=NINT(charge_ex/sum_a) ; mass_ex=mass_ex/sum_a
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
23	  FORMAT(a,': �l�ment fictif /= CNO, de masse : ',i3,/,
	1 'et de charge :',i3)
	 END SELECT
	 	 
!	 PRINT*,nchim
!	 WRITE(*,2000)nucleo(1:nchim) 
	 
! d�termination des abondances initiales, a(�quation,�l�ment)
	 ALLOCATE(a(nchim,nchim),indpc(nchim),b(1,nchim))
	 a=0.d0 ; b=0.d0 ; indpc=1	
		
	 a(1,1)=nucleo(1)  	!H1
	 b(1,1)=x0
	
	 a(2,2)=nucleo(2)  	!He3
	 a(2,3)=nucleo(3)	!He4
	 b(1,2)=y0

	 DO j=4,10
	  a(3,j)=nucleo(j)	!somme i > 4 comp(i)*nucleo(i)=Z0
	  a(4,j)=-abon_rela(6)	!somme comp(i) C, C/Z
	  a(5,j)=-abon_rela(7)	!somme comp(i) N, N/Z	 
	  a(6,j)=-abon_rela(8)	!somme comp(i) O, O/Z
	 ENDDO
	 b(1,3)=z0			!Z
	
	 a(4,4)=a(4,4)+1.d0	!C12	 		
	 a(4,5)=a(4,5)+1.d0	!C13
	
	 a(5,6)=a(5,6)+1.d0	!N14	 		
	 a(5,7)=a(5,7)+1.d0	!N15
	
	 a(6,8)=a(6,8)+1.d0	!O16	 		
	 a(6,9)=a(6,9)+1.d0	!O17
	
! rapports isotopiques		
!	 a(7,1)=1.d0		!H1		
	 a(7,2)=1.d0		!He3
	 a(7,3)=-he3she4z	!He3/He4 avec H2 dans He3		

	 a(8,5)=1.d0		!C13
	 a(8,4)=-c13sc12	!C13/C12
	
	 a(9,7)=1.d0		!N15
	 a(9,6)=-n15sn14	!N15/N14
	
	 a(10,9)=1.d0		!O17
	 a(10,8)=-o17so16	!O17/O16
	
!	 PRINT*,nchim
!	 DO i=1,nchim
!	  WRITE(*,2002)a(i,1:nchim),b(1,i)
!	 ENDDO

	 CALL gauss_band(a,b,indpc,nchim,nchim,nchim,1,inversible)
	 IF(.not.inversible)THEN
	  PRINT*,'ppcno9, matrice du calcul des abondances non inversible'
	  PRINT*,'ARRET'
	  STOP
	 ENDIF

! allocations diverses
	 DEALLOCATE(drt,dro,r,q,dqt,dqo,dmuex)
	 ALLOCATE(ab_ini(nchim),ab_min(nchim),drt(nreac),dro(nreac),
	1 r(nreac),q(nreac),dqt(nreac),dqo(nreac),anuc(nchim),
	2 dmuex(nchim),dh2x(nchim),denx(nchim),dbe7x(nchim),dli7x(nchim))

! abondances initiales et abondances n�gligeables	 
	 comp(1:nchim)=MAX(1.d-29,b(1,1:nchim))
	 ab_ini(1:nchim)=comp(1:nchim)*nucleo(1:nchim)
	 ab_min=ab_ini*fmin_abon	 
	 
! nombre/volume des m�taux dans Z		
	 nbz=SUM(comp(ihe4+1:nchim))	 
	 
! abondances en DeX, H=12
	 ALLOCATE(comp_dex(nchim))
	 comp_dex=12.d0+LOG10(comp/comp(1))
	 
! �critures
	 SELECT CASE(langue)	  
	 CASE('english')	
	  WRITE(2,1002) ; WRITE(*,1002) 
1002	  FORMAT(/,'PP + CNO thermonuclear reactions',/)
	  WRITE(2,1003)nreac ; WRITE(*,1003)nreac 
1003	  FORMAT('number of reaction : ',i3)
	  WRITE(2,1004)nchim ; WRITE(*,1004)nchim
1004	  FORMAT('number of species : ',i3)
	  WRITE(2,1020)x0,y0,z0,z0/x0 ; WRITE(*,1020)x0,y0,z0,z0/x0
1020	  FORMAT(/,'Initial abundances/mass computed with :',/,
	1 'X0=',es10.3,', Y0=',es10.3,', Z0=',es10.3,/,'Z0/X0=',es10.3,/,
	2 'H1=X0, H2+He3+He4=Y0, with H2 in He3',/,
	3 'Z0 = 1-X0-Y0 = C12+C13+N14+N15+O16+O17+Ex',/)	
	  WRITE(2,1)ab_ini(1:nchim) ; WRITE(*,1)ab_ini(1:nchim)
1	  FORMAT('H1 :',es10.3,', He3 :',es10.3,', He4 :',es10.3,
	1 ', C12 :',es10.3,', C13 :',es10.3,/,'N14 :',es10.3,
	2 ', N15 :',es10.3,', O16 :',es10.3,', O17 :',es10.3,
	3 ', Ex :',es10.3)
	  WRITE(2,1009)comp_dex ; WRITE(*,1009)comp_dex
1009	  FORMAT(/,'Initial abundances/number: 12+Log10(Ni/Nh)',/,
	1 'H1:',es10.3,', He3:',es10.3,', He4:',es10.3,
	2 ', C12:',es10.3,', C13:',es10.3,/,'N14:',es10.3,
	3 ', N15 :',es10.3,', O16 :',es10.3,', O17 :',es10.3,
	4 ', Ex :',es10.3)
	  WRITE(2,1021)(comp(4)+comp(5))/nbz,(comp(6)+comp(7))/nbz,
	1 (comp(8)+comp(9))/nbz,comp(10)/nbz
	  WRITE(*,1021)(comp(4)+comp(5))/nbz,(comp(6)+comp(7))/nbz,
	1 (comp(8)+comp(9))/nbz,comp(10)/nbz
1021	  FORMAT(/,'mass ratio by number within Z :',/,'C/Z :',es10.3,
	1 ', N/Z :',es10.3,', O/Z :',es10.3,', Ex/Z :',es10.3)	
	  WRITE(2,1022)(ab_ini(4)+ab_ini(5))/z0,(ab_ini(6)+ab_ini(7))/z0,
	1 (ab_ini(8)+ab_ini(9))/z0,ab_ini(10)/z0
	  WRITE(*,1022)(ab_ini(4)+ab_ini(5))/z0,(ab_ini(6)+ab_ini(7))/z0,
	1 (ab_ini(8)+ab_ini(9))/z0,ab_ini(10)/z0
1022	  FORMAT(/,'mass ratio by mass within Z :',/,'C/Z :',es10.3,
	1 ', N/Z :',es10.3,', O/Z :',es10.3,', Ex/Z :',es10.3)	
	  WRITE(2,1014)he3she4z,c13sc12,n15sn14,o17so16
	  WRITE(*,1014)he3she4z,c13sc12,n15sn14,o17so16
1014	  FORMAT(/,'Isotopic ratios by nomber :',/,
	1 'He3/He4=',es10.3,', C13/C12=',es10.3,
	2 ', N15/N14=',es10.3,', O17/O16=',es10.3)	
	  WRITE(2,1005)ab_min(1:nchim) ; WRITE(*,1005)ab_min(1:nchim)
1005	  FORMAT(/,'threhold for neglectable abundances/mass :',/,
	1 'H1:',es10.3,', He3:',es10.3,', He4:',es10.3,
	2 ', C12:',es10.3,', C13:',es10.3,/,'N14:',es10.3,
	3 ', N15:',es10.3,', O16:',es10.3,', O17:',es10.3,', Ex:',es10.3)
	  WRITE(2,1006) ; WRITE(*,1006)
1006	  FORMAT(/,'H2, Li7, Be7 at equilibrium')
	  WRITE(2,1007) ; WRITE(*,1007)
1007	  FORMAT('Use of a table')
	  WRITE(2,1008) ; WRITE(*,1008)
1008	  FORMAT('Temporal evolution, test of precision on H1 He4')
	 CASE DEFAULT
	  WRITE(2,2) ; WRITE(*,2) 
2	  FORMAT(/,'R�actions thermonucl�aires des cycles PP, CNO',/)
	  WRITE(2,3)nreac ; WRITE(*,3)nreac 
3	  FORMAT('nombre de r�actions : ',i3)
	  WRITE(2,4)nchim ; WRITE(*,4)nchim
4	  FORMAT('nombre d''�l�ments chimiques : ',i3)
	  WRITE(2,20)x0,y0,z0,z0/x0 ; WRITE(*,20)x0,y0,z0,z0/x0
20	  FORMAT(/,'abondances initiales/gramme d�duites de:',/,
	1 'X0=',es10.3,', Y0=',es10.3,', Z0=',es10.3,/,'Z0/X0=',es10.3,/,
	2 'H1=X0, H2+He3+He4=Y0, avec H2 dans He3',/,
	3 'Z0 = 1-X0-Y0 = C12+C13+N14+N15+O16+O17+Ex',/)	
	  WRITE(2,1)ab_ini(1:nchim) ; WRITE(*,1)ab_ini(1:nchim)
	  WRITE(2,9)comp_dex ; WRITE(*,9)comp_dex
9	  FORMAT(/,'Abondances initiales en nombre: 12+Log10(Ni/Nh)',/,
	1 'H1:',es10.3,', He3:',es10.3,', He4:',es10.3,
	2 ', C12:',es10.3,', C13:',es10.3,/,'N14:',es10.3,
	3 ', N15:',es10.3,', O16:',es10.3,', O17:',es10.3,', Ex:',es10.3)
	  WRITE(2,21)(comp(4)+comp(5))/nbz,(comp(6)+comp(7))/nbz,
	1 (comp(8)+comp(9))/nbz,comp(10)/nbz
	  WRITE(*,21)(comp(4)+comp(5))/nbz,(comp(6)+comp(7))/nbz,
	1 (comp(8)+comp(9))/nbz,comp(10)/nbz
21	  FORMAT(/,'rapports en nombre dans Z:',/,'C/Z:',es10.3,', N/Z:',
	1 es10.3,', O/Z:',es10.3,', Ex/Z:',es10.3)	
	  WRITE(2,22)(ab_ini(4)+ab_ini(5))/z0,(ab_ini(6)+ab_ini(7))/z0,
	1 (ab_ini(8)+ab_ini(9))/z0,ab_ini(10)/z0
	  WRITE(*,22)(ab_ini(4)+ab_ini(5))/z0,(ab_ini(6)+ab_ini(7))/z0,
	1 (ab_ini(8)+ab_ini(9))/z0,ab_ini(10)/z0
22	  FORMAT(/,'rapports en masse dans Z:',/,'C/Z:',es10.3,', N/Z:',
	1 es10.3,', O/Z:',es10.3,', Ex/Z:',es10.3)	
	  WRITE(2,14)he3she4z,c13sc12,n15sn14,o17so16
	  WRITE(*,14)he3she4z,c13sc12,n15sn14,o17so16
14	  FORMAT(/,'Rapports isotopiques en nombre:',/,
	1 'He3/He4=',es10.3,', C13/C12=',es10.3,
	2 ', N15/N14=',es10.3,', O17/O16=',es10.3)	
	  WRITE(2,5)ab_min(1:nchim) ; WRITE(*,5)ab_min(1:nchim)
5	  FORMAT(/,'abondances/gramme n�gligeables:',/,
	1 'H1:',es10.3,', He3:',es10.3,', He4:',es10.3,
	2 ', C12:',es10.3,', C13:',es10.3,/,'N14:',es10.3,
	3 ', N15:',es10.3,', O16:',es10.3,', O17:',es10.3,', Ex:',es10.3)
	  WRITE(2,6) ; WRITE(*,6)
6	  FORMAT(/,'H2, Li7, Be7 � l''�quilibre')
	  WRITE(2,7) ; WRITE(*,7)
7	  FORMAT(/,'on utilise une table')
	  WRITE(2,8) ; WRITE(*,8)
8	  FORMAT(/,'�vol. temporelle, test de pr�cision sur H1 et He4')
	 END SELECT

! d�finitions diverses
	 ab_min=ab_min/nucleo ; anuc=ANINT(nucleo)	!nombre atomique

! nettoyage
	 DEALLOCATE(a,b,comp_dex,indpc)
	 	 
! les r�actions	 	 
	CASE(2)
	 dcomp=0.d0 ; jac=0.d0
	
	 IF(t < t_inf)RETURN	
	
	 CALL rq_reac(comp,t,ro,r,drt,dro,drx,q,dqt,dqo,dqx,mue,dmuex)	

!	 WRITE(*,*)'comp' ; WRITE(*,2000)comp(1:nchim)
!	 WRITE(*,*)'r�actions' ; WRITE(*,2000)r(1:nreac)

! �quations d'�volution
	 dcomp(1)=-(3.d0*r(1)*comp(1)+r(8)*comp(4)+r(9)*comp(5)+r(10)*comp(6)
	1 +(r(11)+r(12))*comp(7)+r(13)*comp(8)+r(14)*comp(9))*comp(1)
	2 +(2.d0*r(3)*comp(2)-r(4)*comp(3))*comp(2)			!H1
	 dcomp(2)=r(1)*comp(1)**2-(2.d0*r(3)*comp(2)+r(4)*comp(3))*comp(2)	!He3
	 dcomp(3)=(r(3)*comp(2)+r(4)*comp(3))*comp(2)
	1 +(r(12)*comp(7)+r(14)*comp(9))*comp(1)			!He4
	 dcomp(4)=(-r(8)*comp(4)+r(12)*comp(7))*comp(1)			!C12
	 dcomp(5)=(r(8)*comp(4)-r(9)*comp(5))*comp(1)			!C13
	 dcomp(6)=(r(9)*comp(5)-r(10)*comp(6)+r(14)*comp(9))*comp(1)	!N14
	 dcomp(7)=(r(10)*comp(6)-(r(11)+r(12))*comp(7))*comp(1)		!N15
	 dcomp(8)=(r(11)*comp(7)-r(13)*comp(8))*comp(1)			!O16
	 dcomp(9)=(r(13)*comp(8)-r(14)*comp(9))*comp(1)			!O17

! Pour v�rifications SUM dcomp*nucleo=0
!	 PRINT*,'ppcno9, v�rifications SUM dcomp*nucleo=0'
!	 WRITE(*,2000)DOT_PRODUCT(dcomp,anuc) ; PAUSE'v�rif' 

	 dcomp(10)=-DOT_PRODUCT(dcomp,anuc)/anuc(10) 	!cons. des baryons	
	 	 
! calcul du jacobien
	 IF(deriv)THEN	!jac(i,j) : �quation, i : �l�ment j
	
! �quation dcomp(1)
!	  dcomp(1)=-(3.*r(1)*comp(1)+r(8)*comp(4)+r(9)*comp(5)+r(10)*comp(6)
!	1 +(r(11)+r(12))*comp(7)+r(13)*comp(8)+r(14)*comp(9))*comp(1)
!	2 +(2.*r(3)*comp(2)-r(4)*comp(3))*comp(2)			!H1

	  jac(1,1)=-6.d0*r(1)*comp(1)-r(8)*comp(4)-r(9)*comp(5)
	1 -r(10)*comp(6)-(r(11)+r(12))*comp(7)-r(13)*comp(8)
	2 -r(14)*comp(9)				!d /H1
	  jac(1,2)=4.d0*r(3)*comp(2)-r(4)*comp(3)	!d /He3
	  jac(1,3)=-r(4)*comp(2)			!d /He4
	  jac(1,4)=-r(8)*comp(1)			!d /C12
	  jac(1,5)=-r(9)*comp(1)			!d /C13
	  jac(1,6)=-r(10)*comp(1)			!d /N14
	  jac(1,7)=-(r(11)+r(12))*comp(1)		!d /N15
	  jac(1,8)=-r(13)*comp(1)			!d /O16
	  jac(1,9)=-r(14)*comp(1)			!d /O17
	 
	  DO i=1,nchim		!d�pendances dues � l'effet d'�cran
	   jac(1,i)=jac(1,i) 
	1  -(3.d0*drx(1,i)*comp(1)+drx(8,i)*comp(4)
	2  +drx(9,i)*comp(5)+drx(10,i)*comp(6)
	3  +(drx(11,i)+drx(12,i))*comp(7)
	4  +drx(13,i)*comp(8)+drx(14,i)*comp(9))*comp(1)
	5  +(2.*drx(3,i)*comp(2)-drx(4,i)*comp(3))*comp(2)
	  ENDDO
			 
! �quation dcomp(2)
!	  dcomp(2)=r(1)*comp(1)**2-(2.*r(3)*comp(2)+r(4)*comp(3))*comp(2)!He3

	  jac(2,1)=2.d0*r(1)*comp(1)			!d /H1
	  jac(2,2)=-4.d0*r(3)*comp(2)-r(4)*comp(3)	!d /He3
	  jac(2,3)=-r(4)*comp(2)			!d /He4

	  DO i=1,nchim		!d�pendances dues � l'effet d'�cran
	   jac(2,i)=jac(2,i)
	1  +drx(1,i)*comp(1)**2-(2.d0*drx(3,i)*comp(2)
	2  +drx(4,i)*comp(3))*comp(2)
	  ENDDO
	 
! �quation dcomp(3)
!	  dcomp(3)=(r(3)*comp(2)+r(4)*comp(3))*comp(2)
!	1	+(r(12)*comp(7)+r(14)*comp(9))*comp(1)	!He4

	  jac(3,1)=r(12)*comp(7)+r(14)*comp(9)		!d /H1
	  jac(3,2)=2.d0*r(3)*comp(2)+r(4)*comp(3)		!d /He3
	  jac(3,3)=r(4)*comp(2)				!d /He4
	  jac(3,7)=r(12)*comp(1)			!d /N15
	  jac(3,9)=r(14)*comp(1)			!d /O17
	 
	  DO i=1,nchim		!d�pendances dues � l'effet d'�cran
	   jac(3,i)=jac(3,i)
	1  +(drx(3,i)*comp(2)+drx(4,i)*comp(3))*comp(2)
	2  +(drx(12,i)*comp(7)+drx(14,i)*comp(9))*comp(1)
	  ENDDO
	 
! �quation dcomp(4)
!	  dcomp(4)=(-r(8)*comp(4)+r(12)*comp(7))*comp(1)	!C12

	  jac(4,1)=-r(8)*comp(4)+r(12)*comp(7)			!d /H1
	  jac(4,4)=-r(8)*comp(1)				!d /C12
	  jac(4,7)=r(12)*comp(1)				!d /N15
	 
	  DO i=1,nchim		!d�pendances dues � l'effet d'�cran
	   jac(4,i)=jac(4,i)
	1  +(-drx(8,i)*comp(4)+drx(12,i)*comp(7))*comp(1)
	  ENDDO
	 	 
! �quation dcomp(5)
!	  dcomp(5)=(r(8)*comp(4)-r(9)*comp(5))*comp(1)	!C13
	
	  jac(5,1)=r(8)*comp(4)-r(9)*comp(5)		!d /H1
	  jac(5,4)=r(8)*comp(1)				!d /C12
	  jac(5,5)=-r(9)*comp(1)				!d /C13

	  DO i=1,nchim		!d�pendances dues � l'effet d'�cran
	   jac(5,i)=jac(5,i)+(drx(8,i)*comp(4)-drx(9,i)*comp(5))*comp(1)
	  ENDDO
	
! �quation dcomp(6)
!	  dcomp(6)=(r(9)*comp(5)-r(10)*comp(6)+r(14)*comp(9))*comp(1)	!N14

	  jac(6,1)=r(9)*comp(5)-r(10)*comp(6)+r(14)*comp(9)	!d /H1
	  jac(6,5)=r(9)*comp(1)				!d /C13
	  jac(6,6)=-r(10)*comp(1)			!d /N14
	  jac(6,9)=r(14)*comp(1)			!d /O17

	  DO i=1,nchim		!d�pendances dues � l'effet d'�cran
	   jac(6,i)=jac(6,i)
	1  +(drx(9,i)*comp(5)-drx(10,i)*comp(6)
	2  +drx(14,i)*comp(9))*comp(1)
	  ENDDO
	 	 
! �quation dcomp(7)
!	  dcomp(7)=(r(10)*comp(6)-(r(11)+r(12))*comp(7))*comp(1)	!N15
	
	  jac(7,1)=r(10)*comp(6)-(r(11)+r(12))*comp(7)	!d /H1
	  jac(7,6)=r(10)*comp(1)			            !d /N14
	  jac(7,7)=-(r(11)+r(12))*comp(1)		        !d /N15
	 
	  DO i=1,nchim		!d�pendances dues � l'effet d'�cran
	   jac(7,i)=jac(7,i)
	1  +(drx(10,i)*comp(6)-(drx(11,i)+drx(12,i))*comp(7))*comp(1)
	  ENDDO
	 	 
! �quation dcomp(8)
!	  dcomp(8)=(r(11)*comp(7)-r(13)*comp(8))*comp(1)			!O16

	  jac(8,1)=r(11)*comp(7)-r(13)*comp(8)		!d /H1
	  jac(8,7)=r(11)*comp(1)			!d /N15
	  jac(8,8)=-r(13)*comp(1)			!d /O16
	 
	  DO i=1,nchim		!d�pendances dues � l'effet d'�cran
	   jac(8,i)=jac(8,i)+(drx(11,i)*comp(7)-drx(13,i)*comp(8))*comp(1)
	  ENDDO
	 
! �quation dcomp(9)
!	  dcomp(9)=(r(13)*comp(8)-r(14)*comp(9))*comp(1)			!O17

	  jac(9,1)=r(13)*comp(8)-r(14)*comp(9)		!d /H1
	  jac(9,8)=r(13)*comp(1)			!d /O16
	  jac(9,9)=-r(14)*comp(1)			!d /O17
	 
	  DO i=1,nchim		!d�pendances dues � l'effet d'�cran
	   jac(9,i)=jac(9,i)+(drx(13,i)*comp(8)-drx(14,i)*comp(9))*comp(1)
	  ENDDO		 

! �quation dcomp(10)	 
!	  dcomp(10)=-SUM(anuc*dcomp)/anuc(10)!conservation des baryons

	  DO j=1,10
	   DO i=1,9
	    jac(10,j)=jac(10,j)+anuc(i)*jac(i,j)
	   ENDDO
	   jac(10,j)=-jac(10,j)/anuc(10)	   
	  ENDDO

! unit�s de temps pour int�gration temporelle
	  jac=jac*secon6
	  
!	  PAUSE'apr�s case2 deriv'
	 ENDIF		!deriv

	 dcomp=dcomp*secon6

	CASE(3)

! calcul de la production d'�nergie nucl�aire et d�riv�es
! pour H2(H,g)He3, q(2)H**2=q(2)*r(1)/r(2)	 
	 epsilon(1:4)=0.d0 ; et=0.d0 ; ero=0.d0 ; ex=0.d0
	 IF(t <= t_inf)RETURN
	
	 CALL rq_reac(comp,t,ro,r,drt,dro,drx,q,dqt,dqo,dqx,mue,dmuex)

! mue : nombre d'�lectrons / mole /g = 1/poids mol. moy. par e-
	 IF(comp(1) > 0.d0)THEN
	  h2=r(1)/r(2)*comp(1) ; den=r(6)*mue+r(7)*comp(1)
	  be7=r(4)*comp(2)*comp(3)/den ; li7=r(6)*be7*mue/r(5)/comp(1)
	 ELSE
	  h2=0.d0 ; be7=0.d0 ; li7=0.d0
	 ENDIF
	
!	 PRINT*,'h2,li7,be7' ; WRITE(*,2000)h2,li7,be7

! �nergie PP
	 epsilon(2)=(q(1)*comp(1)+q(2)*h2+q(5)*li7+q(7)*be7)*comp(1)
	1 +(q(3)*comp(2)+q(4)*comp(3))*comp(2)+q(6)*mue*be7
         print*,'In ppcno9 q(1)=',q(1)
         print*,'In ppcno9 q(2)=',q(2)
         print*,'In ppcno9 q(3)=',q(3)
         print*,'In ppcno9 q(4)=',q(4)
         print*,'In ppcno9 q(5)=',q(5)
         print*,'In ppcno9 q(6)=',q(6)
         print*,'In ppcno9 q(7)=',q(7)
         print*,'In ppcno9 eps_ppI=',comp(1)**2*q(1)+q(2)*h2*comp(1)+q(3)*comp(2)**2
! �nergie CNO
	 epsilon(3)=(q(8)*comp(4)+q(9)*comp(5)+q(10)*comp(6)
	1 +(q(11)+q(12))*comp(7)+q(13)*comp(8)+q(14)*comp(9))*comp(1)
! �nergie totale
	 epsilon(1)=epsilon(2)+epsilon(3)
!	 PAUSE'apres case3 avant deriv'
	 
	 IF(deriv)THEN	
	  IF(h2 > 0.d0)THEN
	   dh2t=h2*(drt(1)/r(1)-drt(2)/r(2))
	   dh2ro=h2*(dro(1)/r(1)-dro(2)/r(2))
	   DO i=1,nchim
	    dh2x(i)=h2*(drx(1,i)/r(1)-drx(2,i)/r(2))
	   ENDDO
	   dh2x(1)=dh2x(1)+h2/comp(1)	  
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
	1 +(dqt(3)*comp(2)+dqt(4)*comp(3))*comp(2)
	2 +dqt(6)*mue*be7+(dqt(8)*comp(4)+dqt(9)*comp(5)
	3 +dqt(10)*comp(6)+(dqt(11)+dqt(12))*comp(7)+dqt(13)*comp(8)
	4 +dqt(14)*comp(9))*comp(1)
	5 +(q(2)*dh2t+q(5)*dli7t+q(7)*dbe7t)*comp(1)+q(6)*mue*dbe7t
	
 	  ero=(dqo(1)*comp(1)+dqo(2)*h2+dqo(5)*li7+dqo(7)*be7)*comp(1)
	1 +(dqo(3)*comp(2)+dqo(4)*comp(3))*comp(2)
	2 +dqo(6)*mue*be7+(dqo(8)*comp(4)+dqo(9)*comp(5)
	3 +dqo(10)*comp(6)+(dqo(11)+dqo(12))*comp(7)+dqo(13)*comp(8)
	4 +dqo(14)*comp(9))*comp(1)
	5 +(q(2)*dh2ro+q(5)*dli7ro+q(7)*dbe7ro)*comp(1)+q(6)*mue*dbe7ro

	  ex(1)=2.d0*q(1)*comp(1)+q(2)*h2+q(5)*li7+q(7)*be7
	1 +q(8)*comp(4)+q(9)*comp(5)+q(10)*comp(6)+(q(11)+q(12))*comp(7)
	2 +q(13)*comp(8)+q(14)*comp(9)
	  ex(2)=2.d0*q(3)*comp(2)+q(4)*comp(3)
	  ex(3)=q(4)*comp(2) ; ex(4)=q(8)*comp(1) ; ex(5)=q(9)*comp(1)
	  ex(6)=q(10)*comp(1) ; ex(7)=(q(11)+q(12))*comp(1)
	  ex(8)=q(13)*comp(1) ; ex(9)=q(14)*comp(1)
	 
	  DO i=1,nchim	!contributions des �crans
	   ex(i)=ex(i)+(dqx(1,i)*comp(1)+dqx(2,i)*h2
	1  +q(2)*dh2x(i)+dqx(5,i)*li7+q(5)*dli7x(i)
	2  +dqx(7,i)*be7+q(7)*dbe7x(i))*comp(1)
	3  +(dqx(3,i)*comp(2)+dqx(4,i)*comp(3))*comp(2)
	4  +dqx(6,i)*mue*be7+q(6)*dmuex(i)*be7+q(6)*mue*dbe7x(i)
	5  +(dqx(8,i)*comp(4)+dqx(9,i)*comp(5)
	6  +dqx(10,i)*comp(6)+(dqx(11,i)+dqx(12,i))*comp(7)
	7  +dqx(13,i)*comp(8)+dqx(14,i)*comp(9))*comp(1)
	  ENDDO
!	  PAUSE'apres case3 deriv'
	 
	 ENDIF	!deriv
	 
	CASE(4)		!taux de production des neutrinos

	 IF(t >= t_inf)THEN
	  CALL rq_reac(comp,t,ro,r,drt,dro,drx,q,dqt,dqo,dqx,mue,dmuex)	 
	  be7=r(4)*comp(2)*comp(3)/(r(6)*mue+r(7)*comp(1))
	  hhe=r(1)*comp(1)**2/amu ; be7e=r(6)*mue*be7/amu
	  b8e=r(7)*comp(1)*be7/amu ; n13e=r(8)*comp(1)*comp(4)/amu
	  o15e=r(10)*comp(1)*comp(6)/amu ; f17e=r(13)*comp(1)*comp(8)/amu
	 ELSE
	  hhe=0.d0 ; be7e=0.d0 ; b8e=0.d0
	  n13e=0.d0 ; o15e=0.d0 ; f17e=0.d0
	 ENDIF
!	 PAUSE'apres case4'

	CASE DEFAULT
	 PRINT*,'ppcno0, fait ne peut prendre que les valeurs 1, 2, 3 ou 4'
	 PRINT*,'ERREUR: fait a la valeur:',fait
	 PRINT*,'ARRET' ; PRINT* ; STOP
	 
	END SELECT
	
	RETURN

	END SUBROUTINE ppcno9
