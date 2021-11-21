struct StructInstant
{
	struct tm	tmlocal;	// date et heure locale
	struct tm	tmtu;		// date et heure TU

	double	datejulienne;	// date julienne (avec TU)
	double	t;	 	// nombre de jours julien depuis 31/12/1899 à 12h

	double	tu;		// temps universel		
	double	tsg0;	// temps sidéral à Greenwich à 0hTU
	double	tsg;   	// temps sidéral à Greenwich
	double	ts;		// temps sidéral du lieu
};
typedef	struct StructInstant Instant;

double	InstantInit(Instant *inst, double longitude, long zone, unsigned char btu);
double	Instant_tu2t(Instant *inst, const Instant *Idate, double tuext, long int zone);
double	InstantJulien(const struct tm *tmt, double *phdec);
double	InstantInvJulien(struct tm *tmt, double dj);
void	InstantDiffHeure(Instant *inst, unsigned char tudefined, long int zone);
double	InstantSideral(Instant *inst, double longitude);
double	Instant_hDeci(const struct tm *tmt);
double	Instant_dj2t(Instant *inst);
void	InstantAffiche(Instant * inst);
