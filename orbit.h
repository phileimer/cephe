struct	StructOrbite
{
	Astre	*pA;	// pour l'acc�s � tous les �l�ments de l'astre correspondant

	// Vecteurs d'initialisation des �l�ments orbitaux
	VectData	Ve;
	VectData	Vi;
	VectData	VL;
	VectData	VO;
	VectData	Vob;

	// El�ments calcul�s
	double	a;	// demi grand axe
	double	e;	// excentricit�
	double	i;	// inclinaison sur l'�cliptique
	double	L;	// longitude moyenne
	double	Omega;	// longitude du noeud ascendant
	double	omegab;	// longitude du p�rih�lie
	double	deltaL;	// variation �l�mentaire de L, pour la prise en compte du temps de trajet de la lumi�re

	// Perturbations calcul�es
	double	aper;	// demi grand axe
	double	eper;	// excentricit�
	double	Lper;	// longitude moyenne
	double	vkper;	// longitude du p�rih�lie (associ�e avec e)
	double	rper;	// rayon vecteur
	double	lper;	// longitude
	double	bper;	// latitude

	// El�ments calcul�s d'apr�s les pr�c�dents
	double	M;	// anomalie moyenne = L-omegab
	double	MRsp;	// anomalie moyenne sans perturbations, en radians
	double	omega;	// argument de latitude du p�rih�lie = omegab-Omega

	double	E, ER;	// anomalie excentrique (obtenue par Kepler) en degr�s et radians

	double	r;	// rayon vecteur
	double	v;	// anomalie vraie

	// sp�cifique � la Lune
	VectData	VM, VF, VD;				// �l�ments de l'orbite
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

