struct	StructAstre
{
	Ephemerides	*pE;	// pour l'accès à tous les éléments

	wchar_t	nom[8];		// nom de l'astre
	unsigned short int	astreid;	// numéro d'identification


	Orbite	ElemOrbit;	// éléments de l'orbite

	Vecteur	Coord[6];	// coordonnées

	Instant	Lever;		// instant du lever
	Instant	Coucher;	// instant du coucher
	Instant Meridien;	// instant du passage au méridien

	double	magnitude;	// magnitude relative
	double	diametre;	// diametre apparent
	double	phase;		// phase en degrés
	double	phase100;	// phase en %
	double	elongation;	// élonogation
	double	parallaxe;	// paralaxe

	// spécifique à la Lune
	double	age;
};
//typedef	struct StructAstre	Astre;

void	AstreInit(Astre *a, Ephemerides *pEphe, unsigned short int id);
Vecteur *AstreEcliptGeo(Astre *a);
Vecteur	*AstreEcliptGeo2Horizon(Astre *a, double t);
Vecteur *AstreElemOrbit2EcliptHelio(const Orbite *pOrb);
Vecteur *AstreEcliptHelio2EcliptGeo(const Vecteur *pVin, const Vecteur *GeoSol);
Vecteur *AstreEcliptGeo2EcliptHelio(const Vecteur *pVin, const Vecteur *GeoSol);
Vecteur *AstreEcliptGeo2EcliptGeoApp(Astre *a, const Vecteur *pVin, double lamdaS, double t);
Vecteur *AstreEcliptGeoApp2EcliptGeo(Astre *a, const Vecteur *pVin, double lamdaS);
Vecteur *AstreEcliptGeoApp2Equa(Astre *a, const Vecteur *pVin);
Vecteur *AstreEqua2EcliptGeoApp(Astre *a, const Vecteur *pVin);
Vecteur *AstreEqua2Horaire(Astre *a, const Vecteur *pVin);
Vecteur *AstreHoraire2Equa(Astre *a, const Vecteur *pVin);
Vecteur *AstreHoraire2Horizon(Astre *a, const Vecteur *pVin);
Vecteur *AstreHorizon2Horaire(Astre *a, const Vecteur *pVin);
double	AstreTempsLumiere(double dist);
void	AstreLeverCoucher(Astre *a);
double	AstreEq2tu(Astre *a, const Vecteur *pEq,double trig[3],unsigned short int lcm);
Vecteur	*Astre_tu2Eq(Astre *a, double tu);
unsigned char	AstreTUcond(double *ptu, double tu0, double prec);
double	AstrePhase(Astre *a);
double	AstreElongation(Astre *a);
double	AstreDiametre(Astre *a);
double	AstreParallaxe(Astre *a);
double	AstreMagnitude(Astre *a);
double	AstreAge(Astre *a);
// spécifique à la Lune
Vecteur *AstreElemOrbit2EcliptGeo(Astre *a, const Orbite *pOrb, double t, double *pparalx);
void	AstreAffiche(Astre *a, unsigned short int coord, unsigned char hms);


//////////////////////////////////////////////////////////////////////////////////////
struct 	StructNutation
{
	Ephemerides	*pE;

	double	lamda;		// nutation en longitude g�ocentrique
	double	epsilon;	// nutation en obliquit�

};
typedef	struct StructNutation Nutation;


void	NutationCalcule(Nutation *n, Ephemerides *pEphe);
void	NutationAffiche(Nutation *n);

