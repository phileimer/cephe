/* Jean Philippe EIMER - 1999 *
 *  phil.eimer@9online.fr     */

#include <math.h>

#include "include.h"

//#ifndef NO_CONFIGURE
//	#include "../config.h"
//#endif


void	EphemeridesInitAstres(Ephemerides *Eph)
{
	unsigned char i;


	Eph->maxastr=PLU;
	Eph->maxastr++;

	EphemeridesObliquite(Eph);		// calcul de l'obliquité (epsilon)

	AstreInit(&Eph->Soleil, Eph, SOL);
	AstreEcliptGeo(&Eph->Soleil);

	AstreInit(&Eph->Lune, Eph, LUN);
	AstreEcliptGeo(&Eph->Lune);

	NutationCalcule(&Eph->Nutat, Eph);	// calcul de la nutation

	AstreEcliptGeo2Horizon(&Eph->Soleil, Eph->Observ.Inst.t);
	AstreDiametre(&Eph->Soleil);
	AstreMagnitude(&Eph->Soleil);
	AstreLeverCoucher(&Eph->Soleil);


	AstreEcliptGeo2Horizon(&Eph->Lune, Eph->Observ.Inst.t);
	AstreDiametre(&Eph->Lune);
	AstrePhase(&Eph->Lune);
	AstreElongation(&Eph->Lune);
	AstreAge(&Eph->Lune);
	AstreLeverCoucher(&Eph->Lune);


	for(i = 0 ; i < Eph->maxastr ; i++)
		AstreInit(&Eph->Planete[i], Eph, i);

	for(i = 0 ; i < Eph->maxastr ; i++)
	{
		AstreEcliptGeo(&Eph->Planete[i]);

		AstreEcliptGeo2Horizon(&Eph->Planete[i], Eph->Observ.Inst.t);
		AstreDiametre(&Eph->Planete[i]);
		AstrePhase(&Eph->Planete[i]);
		AstreElongation(&Eph->Planete[i]);
		AstreMagnitude(&Eph->Planete[i]);
		AstreLeverCoucher(&Eph->Planete[i]);
	}

	AstreInit(&Eph->Planete[TER], Eph, TER);
}


void	EphemeridesObliquite(Ephemerides *Eph)
{
	VectData	Veps = TER_EPS;


	Eph->epsilon = VectDataDotT(&Veps, Eph->Observ.Inst.t - 1.0);
}


void	EphemeridesAffiche(Ephemerides *Eph, unsigned char hms)
{
	double	eps = Eph->epsilon;
	wchar_t	dunit[3] = L"  ";


	if(hms)
	{
		eps = Deci2Sexa(eps, NULL, NULL, NULL);
		wcscpy(dunit, L"\'\"");
	}


	wprintf(_(L"Ephémérides Astronomiques"));
    wprintf(L" ("PACKAGE" "VERSION")\n");
	wprintf(L"%S\n", _(L"-------------------------"));
	ObservateurAffiche(&Eph->Observ);
	wprintf(L"\n%S : %7.4f °%2S\n",_(L"Obliquité"), eps, dunit);
	EphemeridesAfficheCoord(Eph, HELIO,hms);
	EphemeridesAfficheCoord(Eph, GEO,hms);
	EphemeridesAfficheCoord(Eph, GEOAPP,hms);
	EphemeridesAfficheCoord(Eph, EQUA,hms);
	EphemeridesAfficheCoord(Eph, HORAIRE,hms);
	EphemeridesAfficheCoord(Eph, HORIZON,hms);
	EphemeridesAfficheAPP(Eph);
	EphemeridesAfficheLMC(Eph);
	wprintf(L"\n");
}


void	EphemeridesAfficheCoord(Ephemerides *Eph, unsigned char coord, unsigned char hms)
{
	wchar_t	symbol1[5], symbol2[5], symbol3[5];
	wchar_t	dunit[] = L"°'\"", hunit[] = L"hms", *unit[3];	// ° is 2 bytes in UTF-8
	unsigned char	i;


	if(!hms)
	{
		wcscpy(dunit, L"°  ");
		wcscpy(hunit, L"h  ");
	}

	unit[COORD_R] = _(L"ua");
	unit[COORD_TETA] = dunit;
	unit[COORD_PHI] = dunit;
	wcscpy(symbol1, L"d   ");


	// Séparation
	wprintf(L"\n");
	for(i = 0 ; i < Eph->maxastr + 3 ; i++)
		wprintf(L"---------");

	// Type de coordonnées
	wprintf(L"\n");

	switch(coord)
	{
		case	HELIO:
			wprintf(_(L"Coordonnées Ecliptiques Héliocentriques"));
			wcscpy(symbol1, L"r   ");
			wcscpy(symbol2, L"l   ");
			wcscpy(symbol3, L"b   ");
			break;

		case	GEO:
			wprintf(_(L"Coordonnées Ecliptiques Géocentriques"));
			wcscpy(symbol2, L"lamb");
			wcscpy(symbol3, L"beta");
			break;

		case	GEOAPP:
			wprintf(_(L"Coordonnées Ecliptiques Géocentriques Apparentes"));
			wcscpy(symbol2, L"lamA");
			wcscpy(symbol3, L"betA");
			break;

		case	EQUA:
			wprintf(_(L"Coordonnées Equatoriales"));
			unit[COORD_TETA] = hunit;
			wcscpy(symbol2, L"alph");
			wcscpy(symbol3, L"delt");
			break;

		case	HORAIRE:
			wprintf(_(L"Coordonnées Horaires"));
			unit[COORD_TETA] = hunit;
			wcscpy(symbol2, L"H   ");
			wcscpy(symbol3, L"delt");
			break;

		case	HORIZON:
			wprintf(_(L"Coordonnées Horizontales"));
			wcscpy(symbol2, L"A   ");
			wcscpy(symbol3, L"h   ");
			break;
	}


	// Separation
	wprintf(L"\n");
	for(i = 0 ; i < Eph->maxastr + 3 ; i++)
		wprintf(L"---------");

	// Noms
	if(coord==HELIO)
		wprintf(L"\n%8S|%8S|%8S|", L"", Eph->Planete[TER].nom, Eph->Lune.nom);
	else
		wprintf(L"\n%8S|%8S|%8S|", L"", Eph->Soleil.nom, Eph->Lune.nom);

	for(i = 0 ; i < Eph->maxastr ; i++)
		wprintf(L"%8S|", Eph->Planete[i].nom);


	if(coord<EQUA)
		EphemeridesAfficheLigneCoord(Eph, symbol1, unit[COORD_R], L"%8.5f|", coord, COORD_R, FALSE);
	EphemeridesAfficheLigneCoord(Eph, symbol2, unit[COORD_TETA], L"%8.4f|", coord, COORD_TETA, hms);
	EphemeridesAfficheLigneCoord(Eph, symbol3, unit[COORD_PHI], L"% 8.4f|", coord, COORD_PHI, hms);
}


void	EphemeridesAfficheLigneCoord(Ephemerides *Eph, wchar_t *symbol, wchar_t *unit, wchar_t *form, unsigned char coord, unsigned char coorditem, unsigned char hms)
{
	unsigned char i;
	double	co,dint;
	unsigned char	heure=(coord==EQUA || coord==HORAIRE) && coorditem==COORD_TETA;


	wprintf(L"\n%4S %3S|",symbol,unit);


	// Soleil ou Terre
	if(coord==HELIO)
		co = Eph->Planete[TER].Coord[coord].Spher[coorditem];
	else
		co = Eph->Soleil.Coord[coord].Spher[coorditem];

	if(heure)
		co /= KH2D;
	if(hms)
		co=Deci2Sexa(co, NULL, NULL, NULL);

	wprintf(form, co);


	// Lune
	co = Eph->Lune.Coord[coord].Spher[coorditem];

	if(heure)
		co/=KH2D;
	if(hms)
		co=Deci2Sexa(co, NULL, NULL, NULL);

	if(coord==GEO && coorditem==COORD_R)
	{
		co = Eph->Lune.parallaxe*60.0;
		wprintf(L"P %2u\'%2u\"|",(unsigned int)co,(unsigned char)(modf(co,&dint)*60.0));
	}
	else
		wprintf(form, co);


	// Planètes
	for(i = 0 ; i < Eph->maxastr ; i++)
	{
		co = Eph->Planete[i].Coord[coord].Spher[coorditem];

		if(heure)
			co/=KH2D;

		if(hms)
			co = Deci2Sexa(co, NULL, NULL, NULL);

		wprintf(form, co);
	}
}


void	EphemeridesAfficheAPP(Ephemerides * Eph)
{
	unsigned char	i;


	// Séparation
	wprintf(L"\n");
	for(i = 0 ; i < Eph->maxastr + 3 ; i++)
		wprintf(L"---------");

	wprintf(L"\n%S", _(L"Phase, Elongation, Diamètre Apparent, Magnitude"));

	// Séparation
	wprintf(L"\n");
	for(i = 0 ; i < Eph->maxastr + 3 ; i++)
		wprintf(L"---------");

	// Noms
	wprintf(L"\n%8S|%8S|%8S|", L"", Eph->Soleil.nom, Eph->Lune.nom);

	for(i = 0 ; i < Eph->maxastr ; i++)
		wprintf(L"%8S|", Eph->Planete[i].nom);

	// phase
	wprintf(L"\n%S °|%8S|%8.2f|", _(L"Phase "), L"-", Eph->Lune.phase);

	for(i = 0 ; i < Eph->maxastr ; i++)
		wprintf(L"%8.2f|", Eph->Planete[i].phase);

	// phase en %
	wprintf(L"\n%S %%|%8S|%8.2f|", _(L"Phase "), L"-", Eph->Lune.phase100);

	for(i = 0 ; i < Eph->maxastr ; i++)
		wprintf(L"%8.2f|", Eph->Planete[i].phase100);

	// élongation
	wprintf(L"\n%S °|%8S|%8.2f|",_(L"Elong "), L"-", Eph->Lune.elongation);

	for(i = 0 ; i < Eph->maxastr ; i++)
		wprintf(L"%8.2f|", Eph->Planete[i].elongation);

	// diametre apparent
	wprintf(L"\n%S \"|%4u\'%2u\"|%4u\'%2u\"|", _(L"Diam  "), (unsigned int)Eph->Soleil.diametre/60, (unsigned int)Eph->Soleil.diametre%60, (unsigned int)Eph->Lune.diametre/60, (unsigned int)Eph->Lune.diametre%60);

	for(i = 0 ; i < Eph->maxastr ; i++)
		wprintf(L"%8.2f|", Eph->Planete[i].diametre);

	// magnitude, age pour la Lune
	i=(unsigned char)Eph->Lune.age;
	wprintf(L"\n%S |%8.2f|A %2uj%2uh|",_(L"Magnit."), Eph->Soleil.magnitude, i, (unsigned short int)((Eph->Lune.age-i)*24.0));

	for(i = 0 ; i < Eph->maxastr ; i++)
		wprintf(L"%8.2f|", Eph->Planete[i].magnitude);
}



void	EphemeridesAfficheLMC(Ephemerides *Eph)
{
	unsigned char	i;


	// Séparation
	wprintf(L"\n");
	for(i = 0 ; i < Eph->maxastr + 3 ; i++)
		wprintf(L"---------");

	wprintf(L"\n%S", _(L"Heures des Lever, Passage au Méridien, Coucher (heures locales)"));

	// Séparation
	wprintf(L"\n");
	for(i = 0 ; i < Eph->maxastr + 3 ; i++)
		wprintf(L"---------");

	// Noms
	wprintf(L"\n%8S|%8S|%8S|", L"", Eph->Soleil.nom, Eph->Lune.nom);

	for(i = 0 ; i < Eph->maxastr ; i++)
		wprintf(L"%8S|", Eph->Planete[i].nom);

	// levers
	wprintf(L"\n%S hms|",_(L"Lev "));
	EphemeridesAfficheHeureLMC(&Eph->Soleil.Lever.tmlocal);
	EphemeridesAfficheHeureLMC(&Eph->Lune.Lever.tmlocal);

	for(i = 0 ; i < Eph->maxastr ; i++)
		EphemeridesAfficheHeureLMC(&Eph->Planete[i].Lever.tmlocal);

	// passages au méridien
	wprintf(L"\n%S hms|",_(L"Mér "));
	EphemeridesAfficheHeureLMC(&Eph->Soleil.Meridien.tmlocal);
	EphemeridesAfficheHeureLMC(&Eph->Lune.Meridien.tmlocal);

	for(i = 0 ; i < Eph->maxastr ; i++)
		EphemeridesAfficheHeureLMC(&Eph->Planete[i].Meridien.tmlocal);


	// couchers
	wprintf(L"\n%S hms|",_(L"Cou "));
	EphemeridesAfficheHeureLMC(&Eph->Soleil.Coucher.tmlocal);
	EphemeridesAfficheHeureLMC(&Eph->Lune.Coucher.tmlocal);

	for(i = 0 ; i < Eph->maxastr ; i++)
		EphemeridesAfficheHeureLMC(&Eph->Planete[i].Coucher.tmlocal);
}


void	EphemeridesAfficheHeureLMC(struct tm *tmt)
{
	if(tmt->tm_hour==TUSING)
		wprintf(L"%8S|", L"-");
	else
		wprintf(L"%2u:%02u:%02u|", tmt->tm_hour, tmt->tm_min, tmt->tm_sec);
}
