/* Jean Philippe EIMER - 1999 *
 *  phil.eimer@9online.fr     */

#include <math.h>

#include "include.h"


Vecteur	Vout;


void	AstreInit(Astre *a, Ephemerides *pEphe, unsigned short int id)
{
	// coordonnées géocentriques pour le Soleil et la Lune
	// éléments orbitaux des planètes

	const wchar_t *Vnom[11] = { \
		MER_NOM, VEN_NOM, MAR_NOM, \
		JUP_NOM, SAT_NOM, URA_NOM, \
		NEP_NOM, PLU_NOM, TER_NOM, \
		SOL_NOM, LUN_NOM	};


	a->pE = pEphe;
	a->astreid = id;
	wcscpy(a->nom, Vnom[id]);

	OrbiteInit(&a->ElemOrbit, a);

	if(id == TER)
		VecteurCopy(&a->Coord[HELIO], AstreElemOrbit2EcliptHelio(&a->ElemOrbit));

	AstreParallaxe(a);
}


Vecteur *AstreEcliptGeo(Astre *a)
{
	OrbiteElements(&a->ElemOrbit, &a->pE->Observ.Inst);

	if(a->astreid == LUN)
	{
		VecteurCopy(&a->Coord[GEO], AstreElemOrbit2EcliptGeo(a, &a->ElemOrbit, a->pE->Observ.Inst.t, &a->parallaxe));
		VecteurCopy(&a->Coord[HELIO], AstreEcliptGeo2EcliptHelio(&a->Coord[GEO], &a->pE->Soleil.Coord[GEO]));
	}
	else
	{
		OrbitePerturbations(&a->ElemOrbit, a->pE->Observ.Inst.t);
		OrbiteAnomalies(&a->ElemOrbit);
		VecteurCopy(&Vout, AstreElemOrbit2EcliptHelio(&a->ElemOrbit));

		if(a->astreid != SOL)
		{
			VecteurCopy(&a->Coord[HELIO], &Vout);
			VecteurCopy(&a->Coord[GEO], AstreEcliptHelio2EcliptGeo(&a->Coord[HELIO], &a->pE->Soleil.Coord[GEO]));
		}
		else
			VecteurCopy(&a->Coord[GEO], &Vout);
	}


	return(&Vout);
}


Vecteur	*AstreEcliptGeo2Horizon(Astre *a, double t)
{
	VecteurCopy(&a->Coord[GEOAPP], AstreEcliptGeo2EcliptGeoApp(a, &a->Coord[GEO], a->pE->Soleil.Coord[GEO].Spher[COORD_TETA], t));
	VecteurCopy(&a->Coord[EQUA], AstreEcliptGeoApp2Equa(a, &a->Coord[GEOAPP]));
	VecteurCopy(&a->Coord[HORAIRE], AstreEqua2Horaire(a, &a->Coord[EQUA]));
	VecteurCopy(&a->Coord[HORIZON], AstreHoraire2Horizon(a, &a->Coord[HORAIRE]));

	return(&a->Coord[HORIZON]);
}


Vecteur *AstreElemOrbit2EcliptHelio(const Orbite *Orb)
{
	Vout.Spher[COORD_R] = Orb->r;
	Vout.Spher[COORD_TETA] = Orb->omega + Orb->v;
	Vout.Spher[COORD_PHI] = 0;

	VecteurRectangulaire(&Vout);

	VecteurRotationX(&Vout, -Orb->i);

	Vout.Spher[COORD_TETA] += Orb->Omega;

	Vout.Spher[COORD_R] += Orb->rper;		// Corrections
	Vout.Spher[COORD_TETA] += Orb->lper;
	Vout.Spher[COORD_PHI] += Orb->bper;

	VecteurRectangulaire(&Vout);

	return(&Vout);
}


Vecteur *AstreEcliptHelio2EcliptGeo(const Vecteur *pVin, const Vecteur *pGeoSol)
{
	Vecteur V1, V2;


	VecteurCopy(&V1, pVin);
	VecteurCopy(&V2, pGeoSol);

	V1.Spher[COORD_TETA] -= pGeoSol->Spher[COORD_TETA];		// l-lamdaS
	VecteurRectangulaire(&V1);

	V2.Spher[COORD_TETA]=0;
	V2.Spher[COORD_PHI]=0;
	VecteurRectangulaire(&V2);

	VecteurPlus(&Vout, &V1, &V2);

	Vout.Spher[COORD_TETA] += pGeoSol->Spher[COORD_TETA];	// l+lamdaS = lamda
	VecteurRectangulaire(&Vout);

	return(&Vout);		// d, lamda et beta
}


Vecteur *AstreEcliptGeo2EcliptHelio(const Vecteur *pVin, const Vecteur *pGeoSol)
{
	Vecteur V1, V2;


	VecteurCopy(&V1, pVin);
	VecteurCopy(&V2, pGeoSol);

	V1.Spher[COORD_TETA] -= pGeoSol->Spher[COORD_TETA];		// lamda-lamdaS
	VecteurRectangulaire(&V1);


	V2.Spher[COORD_TETA] = 0;
	V2.Spher[COORD_PHI] = 0;
	VecteurRectangulaire(&V2);

	VecteurMinus(&Vout, &V1, &V2);

	Vout.Spher[COORD_TETA] += pGeoSol->Spher[COORD_TETA];	// lamda+lamdaS = l
	VecteurRectangulaire(&Vout);

	return(&Vout);		// r, l et b
}


Vecteur *AstreEcliptGeo2EcliptGeoApp(Astre *a, const Vecteur *pVin, double lamdaS, double t)
{
	Orbite	Orb;


	OrbiteCopy(&Orb, &a->ElemOrbit);

	if(a->astreid < a->pE->maxastr)
	{
		// trajet de la lumière
		double	tl = AstreTempsLumiere(pVin->Spher[COORD_R]);
		double	dL = tl * Orb.deltaL;
		Orb.L -= dL;				// modification de la longitude moyenne
		Orb.MRsp -= DEG2RAD(dL);	// modification de l'anomalie moyenne

		if(a->astreid < JUP)
		{
			Orb.L -= Orb.Lper;
			OrbitePerturbations(&Orb, t-tl);
		}
		else
			Orb.Lper = Orb.vkper = Orb.eper = Orb.aper = 0;

		OrbiteAnomalies(&Orb);

		VecteurCopy(&Vout, AstreEcliptHelio2EcliptGeo(AstreElemOrbit2EcliptHelio(&Orb), &a->pE->Soleil.Coord[GEO]));
	}
	else
		VecteurCopy(&Vout, pVin);

	Vout.Spher[COORD_TETA] += a->pE->Nutat.lamda;

	if(a->astreid != LUN)
	{
		double	lsl=DEG2RAD(lamdaS - pVin->Spher[COORD_TETA]);
		Vout.Spher[COORD_TETA] -= COEF_APP * cos(lsl)/cos(DEG2RAD(pVin->Spher[COORD_PHI]));
		Vout.Spher[COORD_PHI] -= COEF_APP * sin(lsl)*sin(DEG2RAD(pVin->Spher[COORD_PHI]));
	}

	VecteurRectangulaire(&Vout);

	return(&Vout);			// lamda App et beta App
}

#if 0
Vecteur *AstreEcliptGeoApp2EcliptGeo(Astre *a, const Vecteur *pVin, double lamdaS)
{
	double	lsl = DEG2RAD(lamdaS - pVin->Spher[COORD_TETA]);


	VecteurCopy(&Vout, pVin);
	Vout.Spher[COORD_TETA] -= a->pE->Nutat.lamda;

	if(a->astreid != LUN)
	{
		Vout.Spher[COORD_TETA] += COEF_APP * cos(lsl)/cos(pVin->Spher[COORD_PHI]);
		Vout.Spher[COORD_PHI] += COEF_APP * sin(lsl)*sin(pVin->Spher[COORD_PHI]);
	}

	VecteurRectangulaire(&Vout);

	return(&Vout);		// lamda et beta
}
#endif

Vecteur *AstreEcliptGeoApp2Equa(Astre *a, const Vecteur *pVin)
{
	VecteurCopy(&Vout, pVin);
	VecteurRotationX(&Vout, -(a->pE->epsilon + a->pE->Nutat.epsilon));

	return(&Vout);		// alpha et delta
}

#if 0
Vecteur *AstreEqua2EcliptGeoApp(Astre *a, const Vecteur *pVin)
{
	VecteurCopy(&Vout, pVin);
	VecteurRotationX(&Vout, a->pE->epsilon + a->pE->Nutat.epsilon);

	return(&Vout);		// lamda et beta
}
#endif

Vecteur *AstreEqua2Horaire(Astre *a, const Vecteur *pVin)
{
	VecteurCopy(&Vout, pVin);

	Vout.Spher[COORD_TETA] = a->pE->Observ.Inst.ts * KH2D - pVin->Spher[COORD_TETA];	// TS-alpha
	VecteurRectangulaire(&Vout);

	return(&Vout);		// H et delta
}

#if 0
Vecteur *AstreHoraire2Equa(Astre *a, const Vecteur *pVin)
{
	VecteurCopy(&Vout, pVin);

	Vout.Spher[COORD_TETA] = a->pE->Observ.Inst.ts * KH2D - pVin->Spher[COORD_TETA];		// TS-H
	VecteurRectangulaire(&Vout);

	return(&Vout);			// alpha et delta
}
#endif

Vecteur *AstreHoraire2Horizon(Astre *a, const Vecteur *pVin)
{
	VecteurCopy(&Vout, pVin);
	VecteurRotationY(&Vout, a->pE->Observ.latitude - 90.0);

	return(&Vout);			// A et h, azimuth et hauteur
}

#if 0
Vecteur *AstreHorizon2Horaire(Astre *a, const Vecteur *pVin)
{
	VecteurCopy(&Vout, pVin);
	VecteurRotationY(&Vout, 90.0 - a->pE->Observ.latitude);

	return(&Vout);			// H et delta
}
#endif

// Calcul du temps de trajet de la lumière
double	AstreTempsLumiere(double dist)
{
	return(dist*COEF_LUM);
}


// calcul des levers et couchers
void	AstreLeverCoucher(Astre *a)
{
	Vecteur	Eq12;
	double	h0,tu,tu0,trig[3];


	h0 = a->parallaxe - a->diametre/7200.0 - a->pE->Observ.refraction - a->pE->Observ.eta1;

	trig[2]=DEG2RAD(a->pE->Observ.latitude);
	trig[1]=sin(trig[2]);
	trig[2]=cos(trig[2]);

	// lever
	tu0 = 12.0;	// départ à 12 h TU
	VecteurCopy(&Eq12, Astre_tu2Eq(a, tu0));
	trig[0] = sin(DEG2RAD(h0 + a->pE->Observ.eta2L));
	tu = AstreEq2tu(a, &Eq12, trig, LEVER);

	while(AstreTUcond(&tu,tu0,PRECISION_LC))
	{
		tu0 = tu;
		tu = AstreEq2tu(a, Astre_tu2Eq(a, tu), trig, LEVER);
	}

	Instant_tu2t(&a->Lever, &a->pE->Observ.Inst, tu, a->pE->Observ.zone);

	// coucher
	tu0 = 12.0;	// départ à 12 h TU
	trig[0] = sin(DEG2RAD(h0 + a->pE->Observ.eta2C));
	tu = AstreEq2tu(a, &Eq12, trig, COUCHER);

	while(AstreTUcond(&tu,tu0,PRECISION_LC))
	{
		tu0 = tu;
		tu = AstreEq2tu(a, Astre_tu2Eq(a, tu), trig, COUCHER);
	}

	Instant_tu2t(&a->Coucher, &a->pE->Observ.Inst, tu, a->pE->Observ.zone);

	// passage au méridien
	tu0 = 12.0;	// départ à 12 h TU
	tu = AstreEq2tu(a, &Eq12, trig, MERIDIEN);

	while(AstreTUcond(&tu, tu0, PRECISION_M))
	{
		tu0 = tu;
		tu = AstreEq2tu(a, Astre_tu2Eq(a, tu), trig, MERIDIEN);
	}

	Instant_tu2t(&a->Meridien, &a->pE->Observ.Inst, tu, a->pE->Observ.zone);
}


// calcul du temps universel à partir de la hauteur et des coordonnées équatoriales
double	AstreEq2tu(Astre *a, const Vecteur *pEq,double trig[3],unsigned short int lcm)
{
	double	H, tu=0;


	if(lcm == MERIDIEN)
		H = 0;
	else
	{
		H = DEG2RAD(pEq->Spher[COORD_PHI]);
		H = (trig[0] - trig[1] * sin(H)) / trig[2] / cos(H);

		if(H >= -1.0 && H <= 1.0)
			H = RAD2DEG(acos(H));
		else
			tu = TUSING;

		if(lcm == LEVER)
			H = -H;
	}

	if(tu != TUSING)
	{
		tu = (H + a->pE->Observ.longitude + pEq->Spher[COORD_TETA]) / KH2D - a->pE->Observ.Inst.tsg0;
		tu = K24(tu);
		tu = tu / SOL2SID;
	}

	return(tu);
}


// calcul des coordonnées équatoriales à partir du temps universel
Vecteur	*Astre_tu2Eq(Astre *a, double tu)
{
	Orbite	Orb;
	Vecteur	Geo, GeoSol;
	Instant	I;


	OrbiteCopy(&Orb, &a->ElemOrbit);

	Instant_tu2t(&I, &a->pE->Observ.Inst, tu, NOZONE);

	OrbiteElements(&Orb, &I);

	if(a->astreid == LUN)
		VecteurCopy(&Geo, AstreElemOrbit2EcliptGeo(a, &Orb, I.t, NULL));
	else
	{
		OrbitePerturbations(&Orb, I.t);
		OrbiteAnomalies(&Orb);

		VecteurCopy(&Geo, AstreElemOrbit2EcliptHelio(&Orb));

		if(a->astreid != SOL)
		{
			OrbiteCopy(&Orb, &a->pE->Soleil.ElemOrbit);
			OrbiteElements(&Orb, &I);
			OrbitePerturbations(&Orb, I.t);
			OrbiteAnomalies(&Orb);
			VecteurCopy(&GeoSol, AstreElemOrbit2EcliptHelio(&Orb));
			VecteurCopy(&Geo, AstreEcliptHelio2EcliptGeo(&Geo, &GeoSol));
		}
	}

	return(AstreEcliptGeoApp2Equa(a, &Geo));
}


unsigned char	AstreTUcond(double *ptu, double tu0, double prec)
{
	unsigned char	cond;


	if(*ptu == TUSING)
		cond = FALSE;
	else
	{
		double abs = fabs(*ptu - tu0);

		if(abs > 12.1)
		{
			*ptu = TUSING;
			cond = FALSE;
		}
		else
			cond = abs > prec;
	}

	return(cond);
}

// calcul de la phase
double	AstrePhase(Astre *a)
{
	a->phase = acos((a->Coord[HELIO].Spher[COORD_R] * a->Coord[HELIO].Spher[COORD_R] + \
				a->Coord[GEOAPP].Spher[COORD_R] * a->Coord[GEOAPP].Spher[COORD_R] - \
				a->pE->Soleil.Coord[GEO].Spher[COORD_R] * a->pE->Soleil.Coord[GEO].Spher[COORD_R]) \
                / 2 / a->Coord[HELIO].Spher[COORD_R] / a->Coord[GEOAPP].Spher[COORD_R]);

	a->phase100 = 100.0 * ( 1.0 + cos(a->phase) ) / 2.0;
	a->phase = RAD2DEG(a->phase);

	return(a->phase);
}


// calcul de l'élongation
double	AstreElongation(Astre *a)
{
	a->elongation = acos((a->pE->Soleil.Coord[GEO].Spher[COORD_R] * a->pE->Soleil.Coord[GEO].Spher[COORD_R] + \
                a->Coord[GEOAPP].Spher[COORD_R] * a->Coord[GEOAPP].Spher[COORD_R] - \
				a->Coord[HELIO].Spher[COORD_R] * a->Coord[HELIO].Spher[COORD_R]) \
                / 2 / a->pE->Soleil.Coord[GEO].Spher[COORD_R] / a->Coord[GEOAPP].Spher[COORD_R]);

	a->elongation = RAD2DEG(a->elongation);

	return(a->elongation);
}


// calcul du diamètre apparent
double	AstreDiametre(Astre *a)
{
	static const double	Vdia[11]={	MER_DIA, VEN_DIA, \
									MAR_DIA, JUP_DIA, SAT_DIA, \
									URA_DIA, NEP_DIA, PLU_DIA, 0, SOL_DIA, LUN_DIA};

	if(a->astreid == LUN)
	{
		a->diametre = 2 * asin(Vdia[LUN] * sin(DEG2RAD(a->parallaxe)));
		a->diametre = RAD2DEG(a->diametre) * 3600.0;				// en secondes d'arc
	}
	else
		a->diametre = Vdia[a->astreid] / a->Coord[GEOAPP].Spher[COORD_R];	// en secondes d'arc

	return(a->diametre);
}


// calcul de la parallaxe
double	AstreParallaxe(Astre *a)
{
	a->parallaxe=0;

	return(a->parallaxe);
}


// calcul de la magnitude
double	AstreMagnitude(Astre *a)
{
	static const VectData	Vmag[11]={MER_VMAG, VEN_VMAG, \
									MAR_VMAG, JUP_VMAG, SAT_VMAG, \
									URA_VMAG, NEP_VMAG, PLU_VMAG, \
									{{0,0,0,0}}, SOL_VMAG, LUN_VMAG};
	VectData	Vphi;

	Vphi.d[0] = 1.0;
	Vphi.d[1] = fabs(a->phase) / 100.0;
	Vphi.d[2] = Vphi.d[1] * Vphi.d[1];
	Vphi.d[3] = Vphi.d[2] * Vphi.d[1];

	a->magnitude = VectDataDot(&Vphi, &Vmag[a->astreid]);

	if(a->astreid != SOL)
		a->magnitude += 5.0 * log10(a->Coord[GEOAPP].Spher[COORD_R] * a->Coord[HELIO].Spher[COORD_R]);

	return(a->magnitude);
}

// Spécifique à la Lune

Vecteur *AstreElemOrbit2EcliptGeo(Astre *a, const Orbite *pOrb, double t, double *pparalx)
{
	static const VectData	Vlamda[50]=LUN_L_VV,Vbeta[45]=LUN_B_VV;
	static const VectData	Vpara[31]=LUN_P_VV;
	static const double	Clamda[50]=LUN_L_VC,Cbeta[45]=LUN_B_VC;
	static const double	Cpara[31]=LUN_P_VC;
	VectData	V,S;
	unsigned short int	i;
	double	paralx;


	S.d[0]=sin(DEG2RAD(VectDataDotT(&pOrb->Va, t)));
	S.d[1]=sin(pOrb->OmegaR);
	S.d[2]=LUN_PER_CB * sin(DEG2RAD(VectDataDotT(&pOrb->Vb, t)));
	S.d[3]=sin(pOrb->CCR);
	S.d[4]=0.0;

	V.d[0]=DEG2RAD(pOrb->M + VectDataDot(&pOrb->VMper, &S));
	V.d[1]=DEG2RAD(pOrb->F + VectDataDot(&pOrb->VFper, &S));
	V.d[2]=DEG2RAD(pOrb->D + VectDataDot(&pOrb->VDper, &S));
	V.d[3]=a->pE->Soleil.ElemOrbit.MRsp + DEG2RAD(LUN_PER_CMS * S.d[0]);
	V.d[4]=0.0;

	//************************************
	// calcul de la longitude géocentrique
	Vout.Spher[COORD_TETA] = 0;
	for(i=45;i<50;i++)
		Vout.Spher[COORD_TETA] += Clamda[i]*sin(VectDataDot(&V, &Vlamda[i]));

	Vout.Spher[COORD_TETA] *= pOrb->EE;

	for(i=27;i<45;i++)
		Vout.Spher[COORD_TETA] += Clamda[i]*sin(VectDataDot(&V, &Vlamda[i]));

	Vout.Spher[COORD_TETA] *= pOrb->EE;

	for(i=0;i<27;i++)
		Vout.Spher[COORD_TETA] += Clamda[i]*sin(VectDataDot(&V, &Vlamda[i]));

	Vout.Spher[COORD_TETA] += pOrb->L + VectDataDot(&pOrb->VLper, &S);

	//***********************************
	// calcul de la latitude géocentrique
	Vout.Spher[COORD_PHI] = pOrb->EE * Cbeta[44] * sin(VectDataDot(&V, &Vbeta[44]));

	for(i=29;i<44;i++)
		Vout.Spher[COORD_PHI] += Cbeta[i] * sin(VectDataDot(&V, &Vbeta[i]));

	Vout.Spher[COORD_PHI] *= pOrb->EE;

	for(i=0;i<29;i++)
		Vout.Spher[COORD_PHI] += Cbeta[i] * sin(VectDataDot(&V, &Vbeta[i]));

	Vout.Spher[COORD_PHI] *= (1.0 - LUN_PER_CW1 * cos(pOrb->OmegaR) - LUN_PER_CW2 * cos(pOrb->CCR));

	//**********************
	// calcul de la parallaxe
	paralx = pOrb->EE * Cpara[30] * cos(VectDataDot(&V, &Vpara[30]));

	for(i=18;i<30;i++)
		paralx += Cpara[i] * cos(VectDataDot(&V, &Vpara[i]));

	paralx *= pOrb->EE;
	for(i=0;i<18;i++)
		paralx += Cpara[i] * cos(VectDataDot(&V, &Vpara[i]));

	Vout.Spher[COORD_R] = TER_RT / sin(DEG2RAD(paralx)) / KUA2KM;
	VecteurRectangulaire(&Vout);

	if(pparalx!=NULL)
		*pparalx=paralx;

	return(&Vout);
}


double	AstreAge(Astre *a)
{
	a->age = fabs(((a->Coord[GEO].Spher[COORD_TETA] > a->pE->Soleil.Coord[GEO].Spher[COORD_TETA] ? a->phase : -a->phase) - 180.0) / 360.0 * LUN_MOIS)/* +1.0*/;

	return(a->age);
}

// Gestion de l'affichage des coordonnées
#if 0
void	AstreAffiche(Astre *a, unsigned short int coord, unsigned char hms)
{
	wprintf(L"\n%S\t: ", a->nom);

	switch(coord)
	{
		case	HELIO:
			wprintf(_(L"Coordonnées Ecliptiques Héliocentriques"));
			VecteurAffiche(&a->Coord[HELIO], HELIO, FALSE, hms);
			break;

		case	GEO:
			wprintf(_(L"Coordonnées Ecliptiques Géocentriques"));
			VecteurAffiche(&a->Coord[GEO], GEO, FALSE, hms);
			break;

		case	GEOAPP:
			wprintf(_(L"Coordonnées Ecliptiques Géocentriques Apparentes"));
			VecteurAffiche(&a->Coord[GEOAPP], GEOAPP, FALSE, hms);
			break;

		case	EQUA:
			wprintf(_(L"Coordonnées Equatoriales"));
			VecteurAffiche(&a->Coord[EQUA], EQUA, TRUE, hms);
			break;

		case	HORAIRE:
			wprintf(_(L"Coordonnées Horaires"));
			VecteurAffiche(&a->Coord[HORAIRE], HORAIRE, TRUE, hms);
			break;

		case	HORIZON:
			wprintf(_(L"Coordonnées Horizontales"));
			VecteurAffiche(&a->Coord[HORIZON], HORIZON, FALSE, hms);
			break;
	}
}
#endif

//////////////////////////////////////////////////////////////////////////////
void	NutationCalcule(Nutation *n, Ephemerides *pEphe)
{
	static const double	Vlam[13][5]=NUT_L_VV,Veps[9][5]=NUT_E_VV;
	static const VectData	Vlc1[4]=NUT_L_VC1,Vec1[2]=NUT_E_VC1;
	static const double Vlc2[9]=NUT_L_VC2,Vec2[7]=NUT_E_VC2;
	VectData	V;
	unsigned short int	i;


	n->pE = pEphe;

	V.d[0] = DEG2RAD(n->pE->Soleil.ElemOrbit.L);
	V.d[1] = DEG2RAD(n->pE->Soleil.ElemOrbit.M);
	V.d[2] = DEG2RAD(n->pE->Lune.ElemOrbit.L);
	V.d[3] = DEG2RAD(n->pE->Lune.ElemOrbit.M);
	V.d[4] = n->pE->Lune.ElemOrbit.OmegaR;

	// Nutation en longitude
	n->lamda = 0;

	for(i=0;i<4;i++)
		n->lamda += VectDataDotT(&Vlc1[i], n->pE->Observ.Inst.t) * sin(VectDataDot5(&V, Vlam[i]));

	for(i=4;i<13;i++)
		n->lamda += Vlc2[i-4] * sin(VectDataDot5(&V, Vlam[i]));

	// Nutation en obliquité
	n->epsilon = VectDataDotT(&Vec1[0], n->pE->Observ.Inst.t) * cos(VectDataDot5(&V, Veps[0]));
	n->epsilon += VectDataDotT(&Vec1[1], n->pE->Observ.Inst.t) * cos(VectDataDot5(&V, Veps[1]));

	for(i=2;i<9;i++)
		n->epsilon += Vec2[i-2] * cos(VectDataDot5(&V, Veps[i]));

	// passage des secondes d'arc en degrés
	n->lamda /= 3600.0;
	n->epsilon /= 3600.0;
}

#if 0
void	NutationAffiche(Nutation *n)
{
	wprintf(L"\n%S : %.8f, %S : %.8f", _(L"Nutation en longitude"), n->lamda, _(L"en obliquité"), n->epsilon);
}
#endif
