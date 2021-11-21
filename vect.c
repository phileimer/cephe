/* Jean Philippe EIMER - 1999 *
 *  phil.eimer@9online.fr     */

#include <math.h>
#include <string.h>
#include <stdio.h>

#include "include.h"

#if 0
void VecteurInit(Vecteur *v, double a, double b, double c, unsigned char rect)
{
	if(rect)
	{
		v->Rect[COORD_X]=a;
		v->Rect[COORD_Y]=b;
		v->Rect[COORD_Z]=c;
		VecteurSpherique(v);	// mise à jour des coordonnées sphériques
	}
	else
	{
		v->Spher[COORD_R]=a;
		v->Spher[COORD_TETA]=b;
		v->Spher[COORD_PHI]=c;
		VecteurRectangulaire(v);	// mise à jour des coordonnées rectangulaires
	}
}
#endif

// conversion sphérique --> rectangulaire
Vecteur	*VecteurRectangulaire(Vecteur *v)
{
	double t=DEG2RAD(v->Spher[COORD_TETA]),p=DEG2RAD(v->Spher[COORD_PHI]);
	double cp=cos(p);
	
	v->Rect[COORD_X]=v->Spher[COORD_R]*cos(t)*cp;		// formules
	v->Rect[COORD_Y]=v->Spher[COORD_R]*sin(t)*cp;
	v->Rect[COORD_Z]=v->Spher[COORD_R]*sin(p);

	v->Spher[COORD_TETA]=K360(v->Spher[COORD_TETA]);	// centrage des angles
	v->Spher[COORD_PHI]=K360(v->Spher[COORD_PHI]);
	v->Spher[COORD_PHI]=KNEG(v->Spher[COORD_PHI]);

	return(v);
}


// conversion rectangulaire --> sphérique
Vecteur	*VecteurSpherique(Vecteur *v)
{
	double	p=v->Rect[COORD_X]*v->Rect[COORD_X]+v->Rect[COORD_Y]*v->Rect[COORD_Y];

	v->Spher[COORD_R]=sqrt(p+v->Rect[COORD_Z]*v->Rect[COORD_Z]);		// formule rayon vecteur

	v->Spher[COORD_TETA]=RAD2DEG(atan(v->Rect[COORD_Y]/v->Rect[COORD_X]));	// formule angle 1

	if(v->Rect[COORD_X]<0)
		v->Spher[COORD_TETA]+=180;	// lever d'indétermination suite à atan

	v->Spher[COORD_PHI]=RAD2DEG(atan(v->Rect[COORD_Z]/sqrt(p)));	// formule angle 2

	v->Spher[COORD_TETA]=K360(v->Spher[COORD_TETA]);	// centrage des angles
	v->Spher[COORD_PHI]=K360(v->Spher[COORD_PHI]);
	v->Spher[COORD_PHI]=KNEG(v->Spher[COORD_PHI]);

	return(v);
}


// rotation autour de Ox
Vecteur	*VecteurRotationX(Vecteur *v, double angle)
{
	double	t,ar=DEG2RAD(angle),ca=cos(ar),sa=sin(ar);
	
	t=v->Rect[COORD_Y]*ca+v->Rect[COORD_Z]*sa;		// formules
	v->Rect[COORD_Z]=v->Rect[COORD_Z]*ca-v->Rect[COORD_Y]*sa;
	v->Rect[COORD_Y]=t;

	VecteurSpherique(v);		// mise à jour des coordonnées sphériques

	return(v);
}


// rotation autour de Oy
Vecteur	*VecteurRotationY(Vecteur *v, double angle)
{
	double	t,ar=DEG2RAD(angle),ca=cos(ar),sa=sin(ar);
	
	t=v->Rect[COORD_X]*ca+v->Rect[COORD_Z]*sa;		// formules
	v->Rect[COORD_Z]=v->Rect[COORD_Z]*ca-v->Rect[COORD_X]*sa;
	v->Rect[COORD_X]=t;

	VecteurSpherique(v);		//mise à jour des coordonnées sphériques

	return(v);
}


// surcharge +
Vecteur	*VecteurPlus(Vecteur *vs, const Vecteur *v1, const Vecteur *v2)
{
	vs->Rect[COORD_X]=v1->Rect[COORD_X]+v2->Rect[COORD_X];		// somme vectorielle
	vs->Rect[COORD_Y]=v1->Rect[COORD_Y]+v2->Rect[COORD_Y];
	vs->Rect[COORD_Z]=v1->Rect[COORD_Z]+v2->Rect[COORD_Z];

	VecteurSpherique(vs);		// mise à jour des coordonnées sphériques
								// de la somme
	return(vs);
}


// surcharge -
Vecteur	*VecteurMinus(Vecteur *vs, const Vecteur *v1, const Vecteur *v2)
{
	vs->Rect[COORD_X]=v1->Rect[COORD_X]-v2->Rect[COORD_X];		// soustraction vectorielle
	vs->Rect[COORD_Y]=v1->Rect[COORD_Y]-v2->Rect[COORD_Y];
	vs->Rect[COORD_Z]=v1->Rect[COORD_Z]-v2->Rect[COORD_Z];

	VecteurSpherique(vs);		// mise à jour des coordonnées sphériques

	return(vs);
}


// surcharge =
Vecteur	*VecteurCopy(Vecteur *vo, const Vecteur *vi)
{
	vo->Spher[COORD_R]=vi->Spher[COORD_R];			// affectation coordonnées sphériques
	vo->Spher[COORD_TETA]=vi->Spher[COORD_TETA];
	vo->Spher[COORD_PHI]=vi->Spher[COORD_PHI];

	vo->Rect[COORD_X]=vi->Rect[COORD_X];			// affectation coordonnées rectangulaires
	vo->Rect[COORD_Y]=vi->Rect[COORD_Y];
	vo->Rect[COORD_Z]=vi->Rect[COORD_Z];

	return(vo);
}


// formatage des résultats
#if 0
void	VecteurAffiche(Vecteur *v, unsigned char coord, unsigned char heure, unsigned char hms)
{
	char u0[3];
	char u1[3]={' ',' ','\0'}, u2[3]={' ',' ','\0'};
	char *u3="°";
	char d[2]={'d','\0'};
	char n1[7],n2[7];
	double th=v->Spher[COORD_TETA],ph=v->Spher[COORD_PHI];


	strcpy(u0, "ua");


	if(heure)
	{
		strcpy(u3, "h");
		th/=KH2D;
	}


	if(hms)
	{
		th=Deci2Sexa(th, NULL, NULL, NULL);
		ph=Deci2Sexa(ph, NULL, NULL, NULL);

		if(heure)
		{
			u1[0]='m';		
			u1[1]='s';
		}
		else
		{
			u1[0]='\'';	// '
			u1[1]='"';
		}

		u2[1]='\'';		// '
		u2[2]='"';
	}
	
	switch(coord)
	{
		case	HELIO:
			d[0]='r';
			strcpy(n1,"     l");
			strcpy(n2,"     b");
			break;

		case	GEO:
			strcpy(n1,"lambda");
			strcpy(n2,"  beta");
			break;

		case	GEOAPP:
			strcpy(n1,"la App");
			strcpy(n2,"be App");
			break;

		case	EQUA:
			strcpy(n1," alpha");
			strcpy(n2," delta");
			break;

		case	HORAIRE:
			strcpy(n1,"     H");
			strcpy(n2," delta");
			break;

		case	HORIZON:
			strcpy(n1,"     A");
			strcpy(n2,"     h");
			break;

		default:
			d[0]='r';
			strcpy(n1,"  teta");
			strcpy(n2,"   phi");
	}

	printf("\n     %s  %s =  %8.5f\tx = % 8.4f",d,u0,v->Spher[COORD_R],v->Rect[COORD_X]);
	printf("\n%s %s%s = %9.5f\ty = % 8.4f", n1, u3, u1, th, v->Rect[COORD_Y]);
	printf("\n%s %s%s = % 9.5f\tz = % 8.4f\n", n2, u3, u2, ph, v->Rect[COORD_Z]);
}
#endif

//////////////////////////////////////////////////////////////////////////////
double	VectDataDot(const VectData *v1, const VectData *v2)
{
	return(v1->d[0]*v2->d[0]+v1->d[1]*v2->d[1]+v1->d[2]*v2->d[2]+v1->d[3]*v2->d[3]);
}


double	VectDataDot5(const VectData *v1, const double v2[5])
{
	return(v1->d[0]*v2[0]+v1->d[1]*v2[1]+v1->d[2]*v2[2]+v1->d[3]*v2[3]+v1->d[4]*v2[4]);
}


double	VectDataDotT(const VectData *v, double t)
{
	return(v->d[0]+(v->d[1]+(v->d[2]+v->d[3]*t)*t)*t);
}


double	VectDataDotTprecis(const VectData *v, double t)
{
	double	dot,dint;

	dot=modf(v->d[1]*t,&dint)*360.0;
	dot+=v->d[0]+(v->d[2]+v->d[3]*t)*t*t;

	return(dot);
}


double	VectDataDotTprecis2(const VectData *v, double t,double dj)
{
	double	dot,dint;

	dot=modf((dj-DJ0)/v->d[1],&dint)*360.0;
	dot+=v->d[0]+(v->d[2]+v->d[3]*t)*t*t;

	return(dot);
}


VectData	*VectDataCopyInit(VectData *v, const double Vd[4])
{
	v->d[0] = Vd[0];
	v->d[1] = Vd[1];
	v->d[2] = Vd[2];
	v->d[3] = Vd[3];
	v->d[4] = 0.0;

	return(v);
}


VectData	*VectDataCopy(VectData *vo, const VectData *vi)
{
	vo->d[0] = vi->d[0];
	vo->d[1] = vi->d[1];
	vo->d[2] = vi->d[2];
	vo->d[3] = vi->d[3];
	vo->d[4] = 0.0;

	return(vo);
}

#if 0
void	VectDataAffiche(VectData *v)
{
	printf("\n{%.5f,%.5f,%.5f,%.5f}",v->d[0],v->d[1],v->d[2],v->d[3]);

	if(v->d[4]!=0.0)
		printf("\t%.5f",v->d[4]);
}
#endif

//////////////////////////////////////////////////////////////////////////////
// conversion Décimal --> Sexagésimal
//double	Deci2Sexa(double angle, int *ph=NULL, int *pm=NULL, int *ps=NULL)
double	Deci2Sexa(double angle, int *ph, int *pm, int *ps)
{
	double	fh,fm,fs;

	fm=modf(fabs(angle),&fh)*60.0;	// extraction des minutes
	fs=modf(fm,&fm)*60.0;		// extraction des secondes (/10000)
								// minutes en nombre entier

	if(fs>59.9)			// gestion de l'arrondi des secondes
	{
		fs=0.0;
		fm+=1.0;
	}

	if(fm==60.0)			// gestion de l'arrondi des minutes
	{
		fm=0.0;
		fh++;
	}

	if(ph!=NULL)
	{
		*ph=(int)fh;

		if(angle < 0.0)
			*ph=-*ph;
	}

	if(pm!=NULL)
		*pm=(unsigned short int)fm;

	if(ps!=NULL)
		*ps=(unsigned short int)fs;

	fm*=0.01;			// minutes /100
	fs*=0.0001;			// secondes /10000

	return(copysign(fh+fm+fs,angle));	// signe rétabli
}
