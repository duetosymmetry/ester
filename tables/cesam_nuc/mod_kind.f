
c******************************************************************

	MODULE mod_kind
	
c module de d�finition des types
c sp: simple pr�cision, dp: double pr�cision

c Auteurs: P.Morel + B. Pichon Laboratoire JD. Cassini, OCA
c CESAM2k

c-----------------------------------------------------------------
	
	INTEGER, PARAMETER, public :: dp=kind(1.d0), sp=kind(1.)
	
	END module mod_kind
