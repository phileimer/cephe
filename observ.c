/* Jean Philippe EIMER - 1999 *
 *  phil.eimer@9online.fr     */

#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include "include.h"


#define LINE_LENGTH 80
#define	TER_RT_M	(TER_RT*1000.0)


float	LitSexa(wchar_t *line)
{
	float l,result;


	l = wcstof(line, NULL);
	result = fabs(l);
	result += wcstof(wcschr(line, L' '), NULL)/60.0;
	result += wcstof(wcsrchr(line, L' '), NULL)/3600.0;

	result = copysign(result, l);

	return(result);
}


unsigned char	ObservateurChargeLieu(Observateur *Obs)
{
	FILE	*file;
	char *home;
	char filename[LINE_LENGTH];
	unsigned char	ret = FALSE;


    if( !(home = getenv("HOME")) || strlen(home) >= LINE_LENGTH )
		strcpy(filename, ".");
	else
		strcpy(filename, home);

	strcat(filename, "/.epherc");


	file = fopen(filename, "r");

	if( file == NULL )
		wprintf(L"\r\nFichier de configuration ~/.epherc non trouvé.\r\nUtilisation des coordonnées par défaut.\r\n\r\n");
	else
	{
		wchar_t word[LINE_LENGTH];
		wchar_t line[LINE_LENGTH];

		do
		{
			fgetws(line, LINE_LENGTH, file);
//			fwscanf(file, L"%S\n", &line);
//			wprintf(line);
		}
		while( *line == L'#' );

		unsigned char	index2 = wcscspn(line, L"\t");
		wcsncpy(Obs->nom, line, (index2 > 15 ? 16 : index2));

		unsigned char	index1 = index2 + 1;
		index2 = wcscspn(line + index1, L"\t");
		wcsncpy(word, line + index1, index2);
		word[index2]='\0';
		Obs->longitude = LitSexa(word);

		index1 += index2 + 1;
		index2 = wcscspn(line + index1, L"\t");
		wcsncpy(word, line + index1, index2);
		word[index2]='\0';
		Obs->latitude = LitSexa(word);

		index1 += index2 + 1;
		index2 = wcscspn(line+index1, L"\t");
		wcsncpy(word, line + index1, index2);
		word[index2]='\0';
		Obs->altitude = wcstol(word, NULL, 10);

		index1 += index2 + 1;
		wcscpy(word, line + index1);
		Obs->zone = wcstol(word, NULL, 10);
	//	Obs->zone++;		// dst
		Obs->zone *= 3600.0;

		Obs->refraction = 0.61;
		Obs->eta1 = RAD2DEG(acos(TER_RT_M / (TER_RT_M + Obs->altitude)));
		Obs->eta2L = RAD2DEG(0);	// atan(h/d);
		Obs->eta2C = RAD2DEG(0);	// atan(h/d);

		fclose(file);

		ret = TRUE;
	}

	return( ret );
}


void	ObservateurInitLieu(Observateur *Obs, wchar_t *n, double lng, double lat, int alt, long zone)
{
	wcscpy(Obs->nom, n);
	Obs->longitude = lng;
	Obs->latitude = lat;
	Obs->altitude = alt;
	Obs->zone = zone;

	Obs->refraction = 0.61;
	Obs->eta1 = RAD2DEG(acos(TER_RT_M / (TER_RT_M + Obs->altitude)));
	Obs->eta2L = RAD2DEG(0);	// atan(h/d);
	Obs->eta2C = RAD2DEG(0);	// atan(h/d);
}


double	ObservateurInitTempsSyst(Observateur *Obs)
{
	time_t	timet;
	struct tm	*tmtl;

	time(&timet);
	tmtl = localtime(&timet);
	Obs->Inst.tmlocal = *tmtl;

	tmtl = gmtime(&timet);
	Obs->Inst.tmtu = *tmtl;

#ifndef _WIN32
	Obs->zone = Obs->Inst.tmlocal.tm_gmtoff;
#endif	// _WIN32

	InstantInit(&Obs->Inst, Obs->longitude, Obs->zone, TPSLOC);

	return(Obs->Inst.tu);
}

#if 0
double	ObservateurInitTemps(Observateur *Obs, unsigned char jour, unsigned char mois, int annee, unsigned char heure, unsigned char min, unsigned char sec, unsigned char btu)
{
/*	time_t	timet;
	struct tm	*tmtl;
*/

	if(btu == TPSTU)		// heure entrée en TU
	{
		Obs->Inst.tmtu.tm_mday = jour;
		Obs->Inst.tmtu.tm_mon = mois-1;
		Obs->Inst.tmtu.tm_year = annee-1900;

		Obs->Inst.tmtu.tm_hour = heure;
		Obs->Inst.tmtu.tm_min = min;
		Obs->Inst.tmtu.tm_sec = sec;
/*
		timet = timegm(&Obs->Inst.tmtu);
		tmtl = localtime(&timet);
		Obs->Inst.tmlocal = *tmtl;

		InstantInit(&Obs->Inst, Obs->longitude);
*/	}
	else		// heure entrée comme heure locale
	{
		Obs->Inst.tmlocal.tm_mday = jour;
		Obs->Inst.tmlocal.tm_mon = mois-1;
		Obs->Inst.tmlocal.tm_year = annee-1900;

		Obs->Inst.tmlocal.tm_hour = heure;
		Obs->Inst.tmlocal.tm_min = min;
		Obs->Inst.tmlocal.tm_sec = sec;

/*		Obs->Inst.tmlocal.tm_isdst = 0;

		timet = mktime(&Obs->Inst.tmlocal);
		tmtl = localtime(&timet);
		Obs->Inst.tmlocal = *tmtl;		// pour fixer le DST

		Obs->Inst.tmlocal.tm_mday = jour;
		Obs->Inst.tmlocal.tm_mon = mois-1;
		Obs->Inst.tmlocal.tm_year = annee-1900;

		Obs->Inst.tmlocal.tm_hour = heure;
		Obs->Inst.tmlocal.tm_min = min;
		Obs->Inst.tmlocal.tm_sec = sec;
		timet = mktime(&Obs->Inst.tmlocal);
		tmtl = gmtime(&timet);
		Obs->Inst.tmtu = *tmtl;

		InstantInit(&Obs->Inst, Obs->longitude);
*/	}

//	Obs->zone = Obs->Inst.tmlocal.tm_gmtoff;
//	Obs->Inst.tmlocal.tm_gmtoff = Obs->zone / 3600;
	InstantInit(&Obs->Inst, Obs->longitude, Obs->zone, btu);

	return(Obs->Inst.tu);
}
#endif

void	ObservateurAffiche(Observateur *Obs)
{
	InstantAffiche(&Obs->Inst);

	wprintf(L"\n%S : %S",_(L"Lieu de l'Observation"), Obs->nom);
	wprintf(L"\n%S : %S %8.4f °\'\"",_(L"Longitude"), (Obs->longitude>0?_(L"O"):_(L"E")), fabs(Deci2Sexa(Obs->longitude, NULL, NULL, NULL)));
	wprintf(L"\n%S : %S %8.4f °\'\"",_(L"Latitude "), (Obs->latitude>0?_(L"N"):_(L"S")), fabs(Deci2Sexa(Obs->latitude, NULL, NULL, NULL)));
	wprintf(L"\n%S :  %4d m",_(L"Altitude "), Obs->altitude);
}
