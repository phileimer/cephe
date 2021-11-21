/* Jean Philippe EIMER - 1999 *
 *  phil.eimer@9online.fr     */

#include <stdio.h>

#include "include.h"
#include <locale.h>

//#ifndef NO_CONFIGURE
//	#include "../config.h"
//#endif



int	main(int argc,char* argv[])
{
	Ephemerides	E;


	setlocale (LC_ALL, "");
	bindtextdomain (PACKAGE, LOCALEDIR);
	textdomain (PACKAGE);

	if(argc > 1)
		wprintf(L" %S %S - %S\n", PACKAGE, VERSION, L"Jean Philippe EIMER (phil.eimer@9online.fr)");
	else
	{
//		ObservateurInitLieu(&E.Observ, L"Paris", -2.3375, 48.8363888889, 67, 3600L);

		if( !ObservateurChargeLieu(&E.Observ) )
			ObservateurInitLieu(&E.Observ, L"Lyon", -4.8333333, 45.76666667, 248, 3600L /*+ 3600.0*/);


//		ObservateurInitTemps(&E.Observ, 4, 12, 2001, 10, 23, 36, TPSTU);
		ObservateurInitTempsSyst(&E.Observ);

		EphemeridesInitAstres(&E);

//		EphemeridesAffiche(&E, DECI);
		EphemeridesAffiche(&E, HMS);
	}

	return(0);
}
