/* Jean Philippe EIMER - 1999 *
 *  phil.eimer@9online.fr     */

#include <math.h>

#include "include.h"



double	InstantInit(Instant *inst, double longitude, long zone, unsigned char btu)
{
	double	dj;


	if (btu == TPSTU)
	{
		dj = InstantJulien(&inst->tmtu, NULL);
		dj += (double)zone / 86400.0;
		InstantInvJulien(&inst->tmlocal, dj);
	}
	else
	{
		dj = InstantJulien(&inst->tmlocal, NULL);
		dj -= (double)zone / 86400.0;
		InstantInvJulien(&inst->tmtu, dj);
	}

	inst->datejulienne = InstantJulien(&inst->tmtu, &inst->tu);
	Instant_dj2t(inst);
	InstantSideral(inst, longitude);

	return(inst->datejulienne);
}


double	Instant_tu2t(Instant *inst, const Instant *Idate,double tuext, long zone)
{
	if(tuext == TUSING)
	{
		inst->tmtu.tm_mday=0;
		inst->tmtu.tm_mon=0;
		inst->tmtu.tm_year=0;
		inst->tmtu.tm_hour=0;
		inst->tmtu.tm_min=0;
		inst->tmtu.tm_sec=0;
		inst->tmlocal.tm_mday=0;
		inst->tmlocal.tm_mon=0;
		inst->tmlocal.tm_year=0;
		inst->tmlocal.tm_hour=TUSING;
		inst->tmlocal.tm_min=0;
		inst->tmlocal.tm_sec=0;
		inst->datejulienne=0;
		inst->t=0;
	}
	else
	{
		inst->tmtu.tm_mday = Idate->tmtu.tm_mday;
		inst->tmtu.tm_mon = Idate->tmtu.tm_mon;
		inst->tmtu.tm_year = Idate->tmtu.tm_year;
		inst->tmtu.tm_hour=(unsigned char)tuext;
		tuext-=inst->tmtu.tm_hour;
		tuext*=60.0;
		inst->tmtu.tm_min=(unsigned char)tuext;
		tuext-=inst->tmtu.tm_min;
		inst->tmtu.tm_sec=(unsigned char)(tuext*60.0);

		inst->datejulienne=InstantJulien(&inst->tmtu, &inst->tu);
		Instant_dj2t(inst);

		if(zone != NOZONE)
			InstantDiffHeure(inst, TRUE, zone);
	}

	return(inst->t);
}


double	InstantJulien(const struct tm *tmt, double *phdec)
{
	double	dj,hdec;
	unsigned int	a=tmt->tm_year+1900;
	unsigned int	a1=a;
	unsigned short int	m=tmt->tm_mon+1;
	unsigned short int	aa;
	unsigned short int	m1=m+1;


	if(m<3)
	{
		a1--;
		m1+=12;
	}

	aa=a1/100;

	dj=1720994.5+floor(365.25*a1)+floor(30.6*m1)+tmt->tm_mday;

	if((a*100+m)*100+tmt->tm_mday > 15821004)
		dj+=(double)(2-aa+aa/4);

	hdec=Instant_hDeci(tmt);
	dj+=hdec/24.0;

	if(phdec!=NULL)
		*phdec=hdec;

	return(dj);
}


double	InstantInvJulien(struct tm *tmt, double dj)
{
	double	hdec,dint;
	long int a;
	int	b,c,d,e,g;


	a=(long int)floor(dj+0.5);


	if(a>2299160)
	{
		int f=(int)(((double)a-1867216.25)/36524.25);
		a+=1+f-f/4;
	}

	b=(int)(a+1524);
	c=(int)floor(((double)b-122.1)/365.25);
	d=(int)floor((double)c*365.25);
	e=(int)((double)(b-d)/30.6001);
	g=(int)((double)e*30.6);

	tmt->tm_mday = b-d-g;
	tmt->tm_mon = (e<14?e-1:e-13);
	tmt->tm_year = (tmt->tm_mon>2?c-4716:c-4715);
	tmt->tm_mon -= 1;
	tmt->tm_year -= 1900;

	hdec = modf(dj-0.5,&dint)*24.0;

	Deci2Sexa(hdec,&g,&tmt->tm_min,&tmt->tm_sec);
	tmt->tm_hour = (unsigned short int)g;

	return(hdec);
}


void	InstantDiffHeure(Instant *inst, unsigned char tudefined, long int zone)
{
	if(!tudefined)
	{
		InstantInvJulien(&inst->tmtu, InstantJulien(&inst->tmlocal, NULL)-zone/86400.0);
		inst->tmtu.tm_sec = inst->tmlocal.tm_sec;	// conservation des secondes
	}
	else
	{
		InstantInvJulien(&inst->tmlocal, inst->datejulienne+zone/86400.0);
		inst->tmlocal.tm_sec = inst->tmtu.tm_sec;	// conservation des secondes
	}
}


double	InstantSideral(Instant *inst, double longitude)
{
	double	vtsg0[4]=TSG0,t1;
	VectData	V;

	VectDataCopyInit(&V, vtsg0);

	t1=(floor(inst->datejulienne-0.5)+0.5-DJ0)/36525.0-1.0;	// t-1 à 0hTU

	inst->tsg0=VectDataDotT(&V, t1);
	inst->tsg0+=(t1-((double)inst->tmtu.tm_year+1900-2000.0)/100.0)*2400.0;
	inst->tsg0=K24(inst->tsg0);

	inst->tsg=inst->tsg0+inst->tu*SOL2SID;
	inst->tsg=K24(inst->tsg);

	inst->ts=inst->tsg-longitude/KH2D;
	inst->ts=K24(inst->ts);

	return(inst->ts);
}



double	Instant_hDeci(const struct tm *tmt)
{
	double	hdec;

	hdec=tmt->tm_hour+(tmt->tm_min+tmt->tm_sec/60.0)/60.0;
	hdec=K24(hdec);

	return(hdec);
}


double	Instant_dj2t(Instant *inst)
{
	inst->t=(inst->datejulienne-DJ0)/36525.0;

	return(inst->t);
}


void	InstantAffiche(Instant * inst)
{
	wchar_t	str[60];


	wcsftime(str, 60, L"%A %d %B %Y, %T", &inst->tmtu);
	wprintf(L"\n%S : %S",_(L"Date et Heure TU     "), str);
	wcsftime(str, 60, L"%A %d %B %Y, %T, %z (%Z)", &inst->tmlocal);
	wprintf(L"\n%S : %S",_(L"Date et Heure Locales"), str);

	wprintf(L"\n%S : %12.4f",_(L"Date Julienne        "), inst->datejulienne);

	wprintf(L"\n%S (h,ms) : %S0=%7.4f\t%S=%7.4f\t%S=%7.4f\n",_(L"Temps Sidéral "),_(L"TSG"),Deci2Sexa(inst->tsg0, NULL, NULL, NULL),_(L"TSG"),Deci2Sexa(inst->tsg, NULL, NULL, NULL),_(L"TS"),Deci2Sexa(inst->ts, NULL, NULL, NULL));
}
