struct StructVecteur
{
	double	Spher[3];	// coordonnées sphériques
	double	Rect[3];	// coordonnées rectangulaires
};
typedef struct StructVecteur Vecteur;

void VecteurInit(Vecteur *v, double a, double b, double c, unsigned char rect);
Vecteur	*VecteurRectangulaire(Vecteur *v);
Vecteur	*VecteurSpherique(Vecteur *v);
Vecteur	*VecteurRotationX(Vecteur *v, double angle);
Vecteur	*VecteurRotationY(Vecteur *v, double angle);
Vecteur	*VecteurPlus(Vecteur *vs, const Vecteur *v1, const Vecteur *v2);
Vecteur	*VecteurMinus(Vecteur *vs, const Vecteur *v1, const Vecteur *v2);
Vecteur	*VecteurCopy(Vecteur *vo, const Vecteur *vi);
void	VecteurAffiche(Vecteur *v, unsigned char coord, unsigned char heure, unsigned char hms);


struct StructVectData
{
	double	d[5];	// 4 éléments en standard, 1 supplémentaire
};
typedef struct StructVectData VectData;

double	VectDataDot(const VectData *v1, const VectData *v2);
double	VectDataDot5(const VectData *v1, const double v2[5]);
double	VectDataDotT(const VectData *v, double t);
double	VectDataDotTprecis(const VectData *v, double t);
double	VectDataDotTprecis2(const VectData *v, double t,double dj);
VectData	*VectDataCopyInit(VectData *v, const double Vd[4]);
VectData	*VectDataCopy(VectData *vo, const VectData *vi);
void	VectDataAffiche(VectData *v);


// conversion Décimal --> Sexagésimal
double	Deci2Sexa(double,int *,int *,int *);
