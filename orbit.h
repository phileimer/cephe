struct	StructOrbite
{
	Astre	*pA;	// pour l'accès à tous les éléments de l'astre correspondant

	// Vecteurs d'initialisation des éléments orbitaux
	VectData	Ve;
	VectData	Vi;
	VectData	VL;
	VectData	VO;
	VectData	Vob;

	// Eléments calculés
	double	a;	// demi grand axe
	double	e;	// excentricité
	double	i;	// inclinaison sur l'écliptique
	double	L;	// longitude moyenne
	double	Omega;	// longitude du noeud ascendant
	double	omegab;	// longitude du périhélie
	double	deltaL;	// variation élémentaire de L, pour la prise en compte du temps de trajet de la lumière

	// Perturbations calculées
	double	aper;	// demi grand axe
	double	eper;	// excentricité
	double	Lper;	// longitude moyenne
	double	vkper;	// longitude du périhélie (associée avec e)
	double	rper;	// rayon vecteur
	double	lper;	// longitude
	double	bper;	// latitude

	// Eléments calculés d'après les précédents
	double	M;	// anomalie moyenne = L-omegab
	double	MRsp;	// anomalie moyenne sans perturbations, en radians
	double	omega;	// argument de latitude du périhélie = omegab-Omega

	double	E, ER;	// anomalie excentrique (obtenue par Kepler) en degrés et radians

	double	r;	// rayon vecteur
	double	v;	// anomalie vraie

	// spécifique à la Lune
	VectData	VM, VF, VD;				// éléments de l'orbite
	VectData	VLper, VMper, VFper, VDper;	// pour les perturbations
	VectData	Va, Vb, Vc;

	double		F, D;
	double		OmegaR, CCR, EE;

};
typedef	struct StructOrbite Orbite;

void	OrbiteInit(Orbite *o, Astre *pAst);
Orbite	*OrbiteCopyInit(Orbite *o, const double Vinit[6][4]);
Orbite	*OrbiteCopy(Orbite *o, const Orbite *O2);
void	OrbiteElements(Orbite *o, const Instant *pInst);
void	OrbiteAnomalies(Orbite *o);
double	OrbiteAnomalieMoyenne(Orbite *o);
double	OrbiteKepler(Orbite *o);
double	OrbiteAnomalieVraie(Orbite *o);
void	OrbitePerturbations(Orbite *o, double t);
void	OrbitePerturbationsSoleil(Orbite *o, double t);
void	OrbitePerturbationsTellur(Orbite *o, double t);
void	OrbitePerturbationsGeantes(Orbite *o, double t);
void	OrbiteAffiche(Orbite *o);


struct	StructOrbiteLune
{
	AstreLune	*pA;

};
typedef	struct StructOrbiteLune OrbiteLune;

