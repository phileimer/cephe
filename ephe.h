struct	StructEphemerides
{
	Observateur	Observ;

	Astre	Soleil;
	Astre	Planete[9];
	Astre	Lune;

	Nutation	Nutat;

	double	epsilon;	// obliquité : inclinaison de l'équateur terrestre par rapport à l'écliptique

	unsigned short int	maxastr;	// arrêt des calculs à cet astreid
};
//typedef	struct StructEphemerides	Ephemerides;

void	EphemeridesInitAstres(Ephemerides *Eph);
void	EphemeridesObliquite(Ephemerides *Eph);
void	EphemeridesAffiche(Ephemerides *Eph, unsigned char hms);
void	EphemeridesAfficheCoord(Ephemerides *Eph, unsigned char coord, unsigned char hms);
void	EphemeridesAfficheLigneCoord(Ephemerides *Eph, wchar_t *symbol, wchar_t *unit, wchar_t *form, unsigned char coord, unsigned char coorditem, unsigned char hms);
void	EphemeridesAfficheAPP(Ephemerides * Eph);
void	EphemeridesAfficheLMC(Ephemerides *Eph);
void	EphemeridesAfficheHeureLMC(struct tm *tmt);
