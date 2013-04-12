	CHARACTER (len=7), PARAMETER, PUBLIC :: version='V3.2.12'

!Signification du num�ro de version : Va.b.c
!a augmente si les fichiers binaires de reprise *.pms, *.rep, *.dat... changent
!m�me si certains restent compatibles
!b augmente si l'un des fichiers, donn�es (*.don), r�glage, sortie (*.osc)...
!change
!c augmente s'il y a des modifications d'algorithmes (bug, nouvelles
!implantations, ...)

!Journal des am�nagements de CESAM2k

!15/11/08 CESAM2k.V3.2.12
! 1 - bsp_gal.f, newspl_gal.f : optimisation, introduction de mg
! 2 - eq_ini_rota4.f, initialise_rota4 : initialisation de U par eq. diff. lin.
! 3 - eq_diff_rota4.f : condition limite externe sur Tau au lieu de U=0
! 4 - coeff_rota4.f, eq_diff_rota4.f : multiplication par nu des coefficients de
!     l'�quation de U
! 5 - remplacement des routines conv_cm_reza.f et conv_cgm_reza.f, (corrections
!     communiqu�es par Reza Samadi)
! 6 - fcmax.f: le cas precision = 'mx' n'�tait pas envisag�
! 7 - eq_diff_chim.f: mcz_ext --> nuzc_ext
! 8 - evol.f: d�placement des d�finitions de nuzc_ext et mzc_ext,
!     suppression des "bricolages techniques" pour la base de la rotation
! 9 - mod_atm: extension des tables de Marcs introduction/modifications de
!     trho_4000.f, marcs.f, mod_atm.f (B.Pichon), correction d'une coquille
! 10 - pertm*.f: reformulation des pertes/gains de masse, pertm_waldrom.f
!      correction d'un (petit) bug dans la d�termination de la luminosit�.
! 11 - newspl_gal.f: appels � bsp1dn limit�s � l'intervalle [x(1),x(nx)],
!      augmentation du nombre d'abscisses pour int�gration de gauss
! 12 - ppcno10.f : Li7 est hors �quilibre
! 13 - diffus.f: suppression des "bricolages"
! 14 - eq_diff_rota4.f : optimisation de la condition limite de Phi en surface
! 15 - cesam.f : am�nagement des commentaires de la logique de baratine,
!      grille_fixe est impos� avec diffusion du moment angulaire
!      Il est possible de poursuivre un calcul avec des r�saux nucl�aires
!      diff�rents si les �l�ments chimiques pris en compte sont identiques Ex
!      Exemple: ppcno9 et ppcno3a9. Am�lioration du calcul des pourcentages
!      d'�nergie PP, CNO..GRAV, incidences pour les �critures et dans list.f
! 16 - lim_zc.f : adaptation du changement de grille fixe
! 17 - abon_ini.f : introduction des mixtures photosph�rique et m�t�oritique
!      de Asplund Grevesse Sauval 05, table cr�e par Yveline Lebreton
! 18 - calib2k_zams.f, calib2k_pms.f introduction de la possibilit� de modifier
!      le Z/X de calibration
! 19 - implantation de diffw_toul.f : coefficients de diffusion simplifi�s pour
!      la diffusion du moment cin�tique
! 20 - marcs.f, roger.f : �criture de la valeur utilis�e de [Fe/H]
! 21 - Am�nagements dans la notice
! 22 - module mod_exploit : des2k_abon.f introduction du dessin de grad_mu,
!      cr�ation de des2k_dvson.f
! 23 - Suppression des routines ctes*.f, les initialisations sont effectu�es
!      directement dans ini_ctes. Les noms des ensembles de constantes
!      i.e. ctes_94 restent inchang�s. Cr�ation d'ensembles de constantes
!      ctes_31 (de Toulouse), ctes_ba (Basu & Antia) et ctes_gaia.
! 24 - lit_nl.f : possibilit� de faire le calcul avec rot_0 "ie.sans rotation"
!      bien que la vitesse angulaire initiale w_rot soit (par erreur) non nulle
! 25 - Introduction de la correction de Jorgen (astro-ph 0811.100v1) coder
!      NOM_OPA=opa_yveline_jorgen
! 26 - write_nl.f : utilisation de grad_ovi/s pour le signe de ovshti/s
! 27 - Atmosph�res roger & marcs : utilisation des valeurs limites en cas de
!      sorties de table en teff et/ou logg, adaptations dans tdetau.f,
!      trho.f, trho_4000.f
! 28 - Introduction de tacho_ext/full.f: diffusion sous/sur et sous les ZC 
!      selon les Toulousains, pour Li7 solaire calibration, utilisation avec
!      difft_nut_ext/full 

!20/09/08 CESAM2k.V3.2.11
! 1 - d�placement de coeff_vth de mod_etat --> mod_evol
! 2 - diffw_mps.f prise en compte du statut OPTIONAL de dhv
! 3 - coeff_rota4.f correction sur chi_t
! 4 - coeff_rota.f, coeff_rota4.f, eq_diff_rota4.f suppression des d�riv�es
!     dhv et dfrl inutilisables
! 5 - lit_nl.f : suppression de l'utilisation de la condition limite JPZ
!     pour la rotation (il conviendrait de passer � la version V3.3...)
! 6 - mise�� jour des SCRIPTS_BASH

!09/08/08 CESAM2k.V3.2.10
! 1 - ecrit_ascii.f zbar --> z_bar ligne 366
! 2 - bsp_dis.f n --> nx ligne 81
! 3 - cesam.f psi0 : 0.07 -->0.08 par d�faut, 
!    Krot=3,4 : l_demi=.FALSE., ini0=3
! 4 - coeff_rota.f : dfrl(Krot,.. --> dfrl(3,..
! 5 - diffw.f dhv en OPTIONAL, coeff_rota3 suppression de dhv
! 6 - eq_diff_rota4: equation de psi dans les ZC
! 7 - evol : simplification de la gestion des ZC/ZR

!14/07/08 CESAM2k.V3.2.9
! 1 - cesam.f : am�nagement des r�glages pour la pr�cision 'av'
! proportions d'�nergies en % de l'�nergie nucl�aire totale
! ajustement des r�glages avec Krot=3,4
! 2 - evol.f : am�lioration de l'approche de la fin de la MS, 
!  introduction de l_demi dans les pr�cisions (cesam.f, mod_donnees.f)
!  am�lioration de l'algorithme imposant le m�lange radiatif de l'enveloppe
! 3 - lim_zc.f : am�nagement de la suppression des limites ZR/ZC trop proches
! 4 - diffw.f : suppression de deff, dv, dh >= 50

!06/05/08 CESAM2k.V3.2.8
!1 - list : attribution de SAVE � fes, bug signal� par B.Pichon
!2 - cesam : introduction de la possibilit� de poursuivre avec diffusion
! un mod�le initialement sans diffusion et r�ciproquement (souhait de
! F.Th�venin). Limitation du pas temporel maximal suivant la valeur de la masse
! initiale. Introduction de l_demi. Cr�ation du fichier ASCII 4d-2.pms pour
! initialisation PMS de mod�les M>10Msol, et ZAMS m=15Msol, le fichier ASCII
! 4d-2.pms et m150.zams sont mis dans le sous directory EXPLOIT
!3 - routines de r�actions thermonucl�aires : suppression de rot_solid inutile
!4 - eq_diff_chim : r�introduction de Deff avec rotation, voir commentaires
!5 - lim_atm : adjonction d'une loi P(tau) pour Teff > 12000
!6 - mod_donnees : augmentation de pnzc=10. Introduction de l_demi
!7 - evol : introduction de demi, simplification de la logique de suppression
! d'une ZR externe, d�placement de la logique des mvt_dis
!8 - opa.f : la limite pour les opacit�s compton est ramen�e � de 7 � 6.D7
!9 - upddate : suppression de la limitation du pas temporel qui existe dans evol

!01/05/08 CESAM2k.V3.2.7
!1 - etat_opalX/Z : correction de dut=duro*drot+eos(4)*1.e6 au lieu de eos(8)
!    erreur signal�e par Yveline, Jos�phina, Laurent
!2 - ecrit_ascii correction d'un bug signal� par Yveline. Plantage lors de la
! formation d'un fichier ASCII pour un mod�le totalemnt convectif avec ajout
! d'un point au voisinage du centre
!3 - coeff_vth :  adaptation au cas en_masse=.FALSE.

!25/04/08 CESAM2k.V3.2.6
!1 - cesam.f : prise en compte de blabla
!              am�nagement des param�tres pour mod�les en eul�rien
!2 - static_r : remise en fonction
!3 - lim_zc : correction d'un bug affectant la suppression des noyaux radiatifs
!             ou convectifs d'�tendue insignifiante

!20/03/08 CESAM2k.V3.2.5
!1 - ascii.f : correction d'un bug signal� par B.Pichon :
! ligne 144 lire evar(2,:)=EXP(evar(2,:)) et non evar(2,:)=10.d0**evar(2,:)
!2 - le fichier journal change de nom, il devient journal.f

!10/02/08 CESAM2k.V3.2.4
!1 - Rotation : rassemblement des routines resout_rota3/4 dans resout_rota.f.
! Cr�ation de la routine jacob_rota.f; conditions limites identiques pour les
! formalismes tz97 et mz04; mise en place de l'option lim_jpz pour les deux
! formalismes. Introduction de ord_rot=m_rot+r_qs pour diffusion du moment
! cin�tique, ord_rot=m_rot sinon; adaptation de lit_binaire.f
!2 - Suppression de la possibilit� d'utiliser les anciens fichiers de donn�es.
!3 - Mise � jour du MAKEFILE du sous-directory EXPLOIT
!4 - Correction d'un commentaire dans marcs.f (cf. B.Pichon).

!01/01/08 CESAM2k.V3.2.3
!1 - cesam.f : Cr�ation de la pr�cision 'aj' (ajuste) pr�cision r�aliste ('pr')
! avec ajustement du dernier pas temporel pour obtenir  Tc, He core, Xc �
! la fin d'une �volution.
!2 - mod_numerique : cr�ation de la routine bsp_gal, permettant la
! d�termination des coefficients splines d'interpolation de n fonctions connues
! en tout point. Le but �tant d'obtenir une estimation des d�riv�es.
! Cr�ation de la routine de changement de base newspl_gal.
! Ces routines utilisent le formalisme des �l�ments finis Galerkin. 
!3 - mod_evol : cr�ation de la routine tab_vth, tabulation de grandeurs
! thermodynamiques sur le base de la composition chimique. Le but �tant le
! calcul des d�riv�es spatiales de ro, mu.. d'o� vaissala. Cr�ation de la
! routine coeff_vth permettant le calcul des quantit�s � tabuler.
! Les tableaux vth et vth_t sont des �l�m�nts de mod_variables.
! Allocation, exploitation de ces tableaux � divers endroits : cesam,
! ecrit_ascii, evol, coeff_rota4.
! 4 - Am�nagements pour l'atmosph�re. Introduction des lois t(tau) hsra et marcs
! ces derni�re sont disponibles sur demande aupr�s de B.Pichon). Correction
! de bugs dans trho et roger, am�liorations de fesh+00.data et fesh+02.data
!(B.Pichon)
! lim_atm : r�tablissement de la limite tau_max=20 avec loi T(tau) non
! purement radiative

!01/10/07 CESAM2k.V3.2.2
!1 - Implantation de la routine mu_mol (Module mod_etat) de calcul de diverses
! quantit�s li�es � �. Introduction dans la routine dgrad pour le crit�re de
! Ledoux, et implantation dans thermo, lim_zc, for037. Introduction de mu_mol
! dans le MAKEFILE du directory EXPLOIT.
!2 - Am�lioration de la d�finition de phi, incidences dans coeff_rota3/4,
! eq_diff_chim, cesam. et ecrit_ascii sur la fr�quence de BV avec new_bv
! Utilisation de new_bv=TRUE pour toutes les pr�cisions sauf 'sa'.
! Dans les fichiers de sortie ASCII, s'il n'est pas tenu compte de la pression
! turbulente, replacement de Ptot/Pgaz=1 [var(23,i)], par grad_mu=dln mu/dln P.
!3 - Le crit�re de Ledoux peut �tre utilis� avec et sans convection.
!4 - Adjonction de g, msol, rsol et lsol dans les glob des fichiers ASCII de
! sortie.
!5 - Dans lim_zc : r�duction du nombre d'allocations, d�placement de la
! d�allocation de new toujours allou� (ligne 1157,remarque de B.Pichon).
! D�placement du calcul et des allocations de v_son et r_son 
! (remarque de B.Pichon).
!6 - Impl�mentation de la semi-convection : routine difft_smc du type difft.
! Introduction dans le MAKEFILE du directory EXPLOIT.
!7 - Dans dnunl suppression des allocations (remarque de B.Pichon).
!8 - Dans cesam, update, (reglages) introduction de la variable de "pr�cision"
! dlntc limitant la variation de la temp�rature centrale sur un pas temporel,
! Correction d'une erreur pour le calcul BV dans la cas new_bv=.TRUE..
!9 - Dans lit_nl : suppression de la possibilit� de l'utilisation des anciens
! fichiers de donn�es (suggestion de B. Pichon).
!10 - Sous-directory EXPLOIT cr�ation du programme de dessin des2k_grad,
! des2k_diff_spl, am�lioration de f037_2k, et de divers programmes de dessin.
!11 - Mise � jour de la NOTICE.
!12 - Sous-directory SCRIPTS, cr�ation de scripts en BASH.

!14/07/07 CESAM2k.V3.2.1
!Dans cesam :
!1- Avec diffusion du moment cin�tique (Krot=3, 4) afin de laisser deff jouer
! pleinement son r�le, re_nu est fix� � 0
!2- Introduction du facteur fmin_abon, d�fini par la pr�cison ou les r�glages,
! fmin_abon permet de r�gler ab_min dans le but d'une am�lioration
! du pas temporel lorsque H br�le en couche. Pr�c�demment
! fmin_abon=0.01 �tait fix�.
!3- Dans ecrit_ascii, ajout de points dans les couches centrales, 
! localisation suivant les param�tres q0 et l0
!4- Dans mod_evol, evol : introduction de coll_rot et colpnt_rot pour la
! d�termination des points de collocation pour la rotation. Un d�centrage des
! points au milieu des l'intervalles permet de stabiliser la solution avec
! m_rot=2, 4,..pair
!5- Dans EXPLOIT : introduction de coll_rot et colpnt dans le MAKEFILE
!6- Dans evol : am�lioration de l'algorithme de l'adaptation des la ZC externe.
!7- Dans coeff_rota3/4 phi_e ne peut �tre nul, + sch�ma implicite

!01/07/07 CESAM2k.V3.2.0
!Dans le fichier reglages : ajout de l0 points � disposer autour des
! discontinuit�s dans les fichiers de sortie ASCII, ajout de la variable
! new_bv calcul de la fr�quence BV utilisant Phi (d ln ro / d ln mu),
! avec, faute de mieux, Phi=1
!Dans ecrit_ascii : adaptation pour ajout des points
!Dans bsp_dis et noeud_dis, adaptation au formalisme de cesam2k 
! du traitement des discontinuit�s dans les routines
!Dans evol : s'il y a retrait du coeur convectif
! 1- sans diffusion, interpolation lin�aire de la comp.chim. dans la
! la zone de retrait, 
! 2- avec et sans diffusion suppression de la discontinuit� au niveau du
! raccord ZC/ZR.
! Mise en r�serve du traitement de l'augmentation de l'abscisse inf�rieure
! des ZC externes.

!Dans cesam : m_ch=2 (interpolation lin�aire) de la composition chimique
! pour toutes les pr�cisions
!Dans mod_numerique implantation de max_local : d�termination des maxima locaux
! pour normalisation dans les dessins
!Dans EXPLOIT, modification des programmes des2k_vaiss et des2k_abon, mise en
! place de max_local dans le MAKEFILE

!14/06/07 CESAM2k.V3.1.3
!Ajout de la routine ctes_94m, identique � ctes_94 avec des valeurs enti�res
!pour les masses atomiques.
!Dans la routine ecrit_ascii, ajout �ventuel d'un point � une distance q0 du
!centre. Ajout de q0 dans le fichier reglage.
!Dans le sous-directory EXPLOIT :
! Installation du zoom dans le programme des2k_abon.
! Modification du zoom dans le programme des2k_vaiss.

!25/05/07 CESAM2k.V3.0.2
!Dans eq_diff_rota3/4, resout_rota3/4 :
! am�nagement des conditions limites
!Dans sous-directory EXPLOIT :
! update des jacobiens pour les calibrations en tenant compte de la diffusion du
! moment cin�tique, programmes calib2k_pms, calib2k_zams
!Dans le sous-directory NOTICE :
! correction de quelques points, des conditions limites pour la rotation

!03/05/07 CESAM2k.V3.0.1
!Dans cesam.f :
! correction d'un bug concernant les proportions pp, CNO, 3a etc..
! m_rot=2 pour la pr�cision 'np'
!Dans coeff_rota4 : �limination du cas nui <= 0
!Dans coeff_rota3 : introduction de phi/delta*gad_mu dans dgrad (remarque de Andy)
!Dans eq_diff_rota3 et eq_diff_rota4 :
! correction d'un bug concernant la d�riv�e de bs(2)
! calcul direct de Lambda, reprenant (id�e de Phi)
!Mise en place de la variable logique baratine permettant de
! d�tourner une grande partie des commentaires 'on line' vers les fichiers
! mon_modele_static, *_atmos, *_evol pour, respectivement, l'�quilibre
! quasi-statique, l'atmosph�re, l'�volution de la composition chimique et
! la rotation (id�e de B.Pichon).
!Dans lit_nl introduction de la possibilit� de lire un fichier .don
! ancien et un fichier .don sans diffusion ni rotation (id�e de B.Pichon).
!Cr�ation du fichier pgplot_factice.f permettant d'�viter l'utilisation du
! PGPLOT pour les dessins "on line"

!03/05/07 CESAM2k.V3.0.0
!Pour la rotation, remplacement de la m�thode des �l�ments finis par celle
!de collocation
!Pour des initialisations, les fichiers binaires *.hom, *.pms sont utilisables

!03/02/07 CESAM2k.V2.5.2
!Pour la rotation : mise en place de formalismes identiques dans ZR et ZC
!avec, dans ZC, des coefficients de diffusion >> 1

!04/02/07 CESAM2k.V2.5.1
!Modification du fichier reglages
!Dans evol.f :
! Suppression de demi
!Dans cesam.f :
! Avec Krot=3,4 on ignore les discontinuit�s de composition
! chimique et on les lisse mvt_dis=.FALSE., lisse = .TRUE.
! pour la diffusion du moment cin�tique utilisation de la formulation approch�e
! de mu(Krot=3,4), mu_saha=.FALSE. (sauf avec precision='rg')
! Pour pr�cision r�aliste & corot (pr & co),
! fonction d'espacement limit�e � Ln P, ajustement de psi0 pour avoir un nb de
! couches de l'ordre de 1000 pour une 1.5Msol sur la ZAMS
! addition de ajuste et lisse dans le fichier reglages
!Dans base_chim.f :
! Utilisation d'une base continue non d�rivable
!Dans base_rota.f :
! Utilisation d'une base continue non d�rivable
!Dans print_ctes.f :
! permutation de Li6 et Li7 dans la liste d'�criture
!Dans coeff_rota3.f, coeff_rota4.f, eq_diff_poisson.f :
! allocation du tableau ion, (cf. B.Pichon)

!19/01/07 CESAM2k.V2.4.6
!Dans coeff_rota4
! Addition de pi dans cte5_0
! Changement de signe de cte_0(4) Bug signal� par Yveline
! Corrections des dln_mu Bug signal� par ANDY
! Addition de mu pour les trac�s :
!  Augmentation � 30 de nb_coeff dans mod_evol
!  Modification de des2k_rot dans EXPLOIT
!Dans resout_rota4, ligne 345 2 au lieu de iv

!11/12/06 CESAM2k.V2.4.5
!Dans coeff_rota3 et 4
!Introduction d'une approximation de dln ro / dln mu
!Utilisation de l'approximation num�rique de d ln P/d nu
!Dans cesam.f, �nergies pp, cno, 3a+C+O, grav en Lsol
!Dans evol permutation de l'ordre des int�grations de diffusion
!composition chimique puis moment cin�tique
!Dans coeff_rota4, utilisation de chi_T
!Dans resout_rota4 suppression du controle des corrections NR

!06/12/06 CESAM2k.V2.4.4
!Am�nagement des tests de d�rivation dans coeff_rota4, eqdiff_rota4
!Correction d'Yveline : dans cesam.f, utilisation de -TdS au lieu de ABS(TdS)
!pour l'estimation de l'�nergie graviphique.

!25/11/06 CESAM2k.V2.4.3
!Am�lioration de la logique du dtmax dans cesam.f
!Suppression des discontinuit�s dans la base de la rotation
!les variables deviennent continues non d�rivable ==> suppression du
!r�tablissement des continuit�s dans resout_rota3/4.
!Augmentation du nombre de variables dans le fichier ASCII *_coeff-rota.dat
!Am�nagement du programme des2k_rot du sous directory EXPLOIT pour le trac� des
!des variables de la rotation
!Limitation � 50My du pas temporel s'il y a diffusion du moment cin�tique
!et limitation � 10% la variation du pas temporel
!Augmentation/diminution du nombre de couches limit� � 5% (au lieu de 10%)

!27/10/06 CESAM2k.V2.4.2
!rectification lnt923" --> lnt932  r�action 5 de NACRE dans taux_nuc.f
!rectification des commentaires concernant mu_e, var(15) est mu_e
!Mise � jour de la notice et de l'aide m�moire

!12/09/06 CESAM2k.V2.4.1
!Introduction de la variable logique nc_max permettant d'imposer le nombre
!maximum de couches n_max pour le calcul du dernier mod�le : nc_max=n_max < 0
!Dans cesam.f, introduction de la variable coox et du fichier *.coox
!combustion de l'oxyg�ne
!Dans opa_houdek9 appel � opa_opal2_co en cas de sortie de table (pour
!l'atmosph�re opa_opal2_co, utilise opa_yveline)

!22/08/06 CESAM2k.V2.4.0
!Introduction des r�actions nucl�aires de la combustion du carbone
!C12(C12,..).. C12(O16,..) et de l'oxyg�ne O16(O16,..).
!Intervalle de tabulation des taux des r�actions fix� � 1MK
!Cr�ation de la pr�cision 'av' avec des am�nagement permettant d'atteindre
!les stades avanc�s
!Cr�ation des routines ppcno3aco  ppcno3acos de la combustion de H � O
!Introduction du vecteur iter_qs permettant d'adapter les variables contr�l�es
!pour la r�solution de l'�quilibre quasi-statique
!Cr�ation de la routine opa_compton et utilisation d�s que T > 0.007T9
!Ecriture de Teff dans le dessin
!Calcul des poids statistiques dans saha
!Dans des_m abondances centrales avec 'av'

!05/07/06 CESAM2k.V2.3.3
!Introduction des r�actions nucl�aires de NACRE pour le 3alpha et C12(a,g)O16
!D�but de la combustion du carbone mise � 6d8K, tabulation pour ppcno3ac10
!repous�e � T9
!Cr�ation de la pr�cision 'av' : long runs pour les longues �volutions
!R�am�nagement du module mod_nuc
!Mise en place de SAVE dans etat_opal* (remarque de JP Marques)

!28/06/06 CESAM2k.V2.3.2
!Modification de la localisation des limites ZR/ZC dans lim_zc
!Pour le calcul des coefficients de diffusion du moment cin�tique, utilisation
!de Omega et U au temps t

!13/06/06 CESAM2k.V2.3.1
!Adaptation de l'expression approch�e de \phi
!Possibilit�s de dessin on line et off line des coefficients de la diffusion
!du moment cin�tique

!24/05/06 CESAM2k.V2.3.0
!Introduction du formalisme de diffusion du moment cin�tique
!selon Mathis & Zahn 2004 occasionnant diverses adaptations du fichier de
!donn�es et de r�glages.
!Cr�ation d'une routine de dessin des variables de la diffusion du moment
!cin�tique
!Utilisation de l'ancienne formule de la fr�quence de Brunt Vaissala
!Correction d'un bug dans les routines PPCNO12, PPCNO12Be, PPCNO12BeBFe,
!PPCNO12Li
!Acc�l�ration de la convergence de la diffusion du moment cin�tique

!24/04/06 CESAM2k.V2.2.0
!Restructuration permettant diff�rents calculs de la vitesse angulaire
!Modification du fichier de donn�es
!Introduction du formalisme de Matis & Zahn 2004 (d�but)
!Introduction de la conservation locale du moment cin�tique

!11/04/06 CESAM2k.V2.1.0
!V�rification du Jacobien de ppcno3a12Ne (bug)
!Red�finition de scale dans rkimps
!Suppression, dans cesam, de Kipp=.TRUE. pour les mod�les apr�s la s�quence
!principale
!Suppression, dans cesam, de n_max de la NAMELIST nl_rlg

!20/03/06 CESAM2k.V2.0.8
!Am�nagements dans etat_opal, ZFSinterppeos, opal_ascii_bin, calib2k
!Suppression des SUM dans opa_opal2 (bug)
!Compl�ments de formules de r�actions nucl�aires NACRE
!Am�nagements mineurs dans des_m, des_r, resout, z14xcotrin21

!20/03/06 CESAM2k.V2.0.7
!Cr�ation du programme d'exploitation des2k_opa.f
!Cr�ation de la routine g�n�rique coeff_rota.f et des routines coeff_rota_saha.f,
!coeff_rota_ioni.f, coeff_rota_z16.f
!Dans ces routines mise � 0 de chi_T mal calcul� avec les donn�es dont on dispose
!Cr�ation de la routine difft_sun.f
!Addition de la viscosit� cin�matique au coefficient Dv
!Cr�ation de la routine coeff_rota_ioni.f
!Introduction de mini, valeur minimale de Dv, dans le module mod_evol
!Dans lit_nl, avec diffusion du moment cin�tique, on impose D_turb >= 50
!Introduction du nom de la routine de calcul des coefficients de rotation dans
!le type de pr�cision et dans le fichier r�glages
!Introduction dans mod_donnees de la variable logique ecrit_rot conditionnant
!l'�criture du fichier mon_modele_coeff_rota.dat pour dessin des coefficients
!de rotation, addition de ecrit_rot dans le fichier reglages
!Cr�ation de la routine difft_gab.f
!Mise � jour des param�tres de pr�cision et de l'aide m�moire
!Am�nagements du test de s�curit� MODIF_CHIM dans abon_ini.f vent.f,
!planetoides.f

!02/02/06 CESAM2k.V2.0.6
!Suppression de w_form du fichier de donn�es, mis dans le fichier reglage
!Cr�ation de la routine g�n�rique coeff_rota.f appelant l'une des deux routines
!coeff_rota_saha.f et coeff_rota_z16.f
!Mise en place des chutes de plan�to�des : modification de la composition
!chimique de la ZC externe et apport/retrait de moment cin�tique, adaptation du
!fichier planet

!01/02/06 CESAM2k.V2.0.5
!Limitation de la source du vent � la ZC externe, suppression du
!param�tre p_vent du fichier vent, simplification du traitement du vent,
!am�nagement du programme fichier_vent.f et des fichiers exemple.vent et vent
!Mise en place des chutes de plan�to�des, cr�ation de la routine plan�to�des et
!des fichiers exemple.planet et planet
!Ajonction de la masse terrestre dans les fichiers ctes85 et ctes94
!Ajonction de W_FORM dans le fichier "reglages"
!En cours :
!Cr�ation, suppression, am�nagement et mise au point des divers routines et
!programmes concern�s par la diffusion du moment cin�tique, principalement :
!cesam, ecrit_rota, diffus,
!resout_rota, resout_chim, coeff_rota, eq_diff_rota, eq_diff_chim, des2k_rot

!15/12/05 CESAM2k.V2.0.4
!Introduction de l'argument optionnel duale de la routine newspl
!Construction et am�nagements de diverses routines pour la diffusion du moment
!cin�tique

!20/10/05 CESAM2k.V2.0.3
!Am�nagements mineurs dans inter
!Correction du calcul de d ln l / d ln m, et d ln ro / d m^2/3
!Possibilit� de calcul direct des coefficients de diffusion du moment cin�tique
!Introduction de tab_coeff_rota
!Introduction dans resout de la variable et fonction logique cmax et fcmax pour
!l'utilisation du nombre maximum de couches avant de sortir
!Suppression du dessin de ro sans diffusion (discontinuit�s)
!Dans evol, sans diffusion,l�ger lissage par contour de la composition chimique
!pour lisser le retrait des ZC
!Ajout de lim et de model_num � la fin des fichiesr binaires

!13/10/05 CESAM2k.V2.0.2
!Calcul direct de Deff dans eq_diff_chim
!Permutation de l'ordre diffusion du moment cin�tique <==> diffusion des
!�l�ments chimiques
!Cr�ation du programme de dessin des2k_coeff_rota
!Ecriture du num�ro du mod�le dans des_m et des_r
!Dans les fichiers de sortie ASCII, ajout de 20 points de grille de part et
!d'autre des limites ZR/ZC pour affiner le profil de la fr�q.BV
!En abscence de diffusion am�lioration de la formulation de la fr�q.BV

!05/10/05 CESAM2k.V2.0.1
!Am�lioration du choix de no_croiss dans noein
!Am�lioration d'�critures et introduction de no_croiss dans linf
!SAVE pour les quantit�s ***0 des conditions limites de static_m

!01/10/05 CESAM2k.V2.0.0
!SAVE dans opa_yveline, etat_opalX, etat_opalZ
!Mise en place des num�ros des mod�les, sorties de tous les mod�les en ASCII
!et .rep avec leur num�ro, conservation du num�ro dans les fichiers binaires
!Choix de grad_ad ou grad_rad dans les zones overshoot�es

!20/09/05 CESAM2k.V1.1.15
!Adjonction de v dans l'expression de teta dans les tests de d�rivation de
!static_m et static_r
!Cr�ation de la routine ppcno3a13Ne22

!01/09:05 CESAM2k.V1.1.14
!Correction de dgravr dans thermo et thermo_atm
!Facteur 2/3 sur l'acc�l�ration centrifuge dans coll_atm et eq_atm
!Am�nagements des �quations relatives � la diffusion du moment cin�tique
!Cr�ation des programmes de dessin des2k_dhve, des2k_rot, des2k_bin

!30/08/05 CESAM2k.V1.1.13
!Suppression de commentaires dans les modules mis dans la notice
!Suppression de nom_elem en dp des d�finition de mod_nuc
!Suppression de la variable pmw du module mod_donnees
!Suppression des tableaux xlim_rot, xcin et xcint du module mod_evol

!03/08/05 CESAM2k.V1.1.12
!Inversion de la chronologie du journal
!Corrections de bugs signal�s par A.Moya dans coeff_rota :
!C12=1 et signes - pour C16,17,18
!R�tablissement de d2U/dnu2=0, et �quation de diffusion de Omega dans ZC
!D�placement de l'allocation de frot, dfrot dans lim_zc
!Annulation de C15*, C8 et C9 dans coeff_rota

!27/06/05 CESAM2k.V1.1.11
! coeff_rota, utilisation de rho ie. sans passer par l'�quation d'�tat,
! pour coh�rence avec dln ro
! des_m, augmentation du nb. de chiffres significatifs pour les abondances max
! diffm_mp changement de signe de l'acc�l�ration centrifuge
! diffm_mp & diffm_br coefficient de l'acc�l�ration centrifuge
! am�nagements dans resout_rota, eqdiff_rota, diffus, coef_rota
! cr�ation du programme de dessin des2k_dhve du sous-directory EXPLOIT

!16/06/05 CESAM2k.V1.1.10
! Suppression du fichier *.atm pour initialiser ZAMS ou PMS
! Trac� de ro dans des_m et des_r
! Abondances des �l�ments au centre dans list
! Augmentation de m_rot --> 4 dans cesam.f pour tous les r�glages
! Apr�s la ZAMS on impose l'approximation de Kippenhahan

!13/06/05 CESAM2k.V1.1.9
! Correction C12(a,g)O16 dans ppcno3ac10
! Elargissement des dessins du HR jusqu'� nb nmax mod�les
! etat_opal SAVE pour la variable iri
! Lxchim(nchim) dans etat_opal et etat_ceff
! lit_nl, nb_max_models pour lit_nl_2korg
! R�duction du pas temporel � la fin de la ZAMS
! Augmentation du nombre de couches apr�s la TAMS, He et C burning
! D�finition de z_table=z0 dans opa_opal2
! diffw routine PUBLIC de mod_evol
! Cr�ation de lit_binaire dans mod_exploit
! Cr�ation du programme des2k_dhve dans EXPLOIT

!03/06/05 CESAM2k.V1.1.8
! Rectification d'une virgule dans mod_exploit
! Suppression de la r�f�rence � compg(ihe4,1) dans list
! Save de cte1 dans colatm
! Mise du num�ro de version dans journal

!31/05/05 CESAM2k.V1.1.7
! Mise en service du programme des2k_rot du sous-directory EXPLOIT
! SAVE et allocation des tables de donn�es dans opa_yveline

!26/05/05 CESAM2k.V1.1.6
! Restriction de l'utilisation de ln ro = bp(7,:) au cas avec diffusion
! et ord_qs  2 � cause de la discontinuit� de ro aux limites ZR/ZC
! Introduction du num�ro de version par un include dans mod_donnees

!20/05/05 CESAM2k.V1.1.5
! Corrections dans lim_gong1, lim_tau1, add_ascii, opa_yveline lisse,
! pour r�implantation des mod�les de GONG

!14/05/05 CESAM2k.V1.1.4
!Sous directory SOURCE:
! Correction du bug dans la formule de normalisation dans evol
! lire	chim(1:nchim,i)=chim(1:nchim,i)/norm
! et non  chim(1:nchim,:)=chim(1:nchim,:)/norm
!Sous directory EXPLOIT:
! Cr�ation des programmes de dessin des2k_abonts et des2k_abontc
! Suppression du programme des2k_abont
! Dessin de X, Y, Z dans des2k_abon

!05/05/05 CESAM2k.V1.1.3
! Correction d'un bug li� au calcul de Teff si n_atm=1 ie. dans les
! cas GONG1 et GONG2
! Calcul ("exact") de la f�quence de Brunt-Vaissala en utilisant ln ro
! Suppression de contour

!02/05/05 CESAM2k.V1.1.2
! Am�lioration des algorithmes g�rant les arr�ts sur t_stop et x_stop
! Implantation de l'arr�t sur He_core

!27/04/05 CESAM2k.V1.1.1
!  Correction de bugs engendr�s par ihe4=-100 et lvent=.TRUE.avec PP1

!23/04/05 CESAM2k.V1.1.0
! Introduction dans resout de l'arr�t sur x_stop
! Introduction dans resout de l'arr�t sur t_stop

