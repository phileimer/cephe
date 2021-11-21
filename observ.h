struct	StructObservateur
{
	// données
	wchar_t	nom[17];		// nom du lieu d'observation
	double	longitude;
	double	latitude;
	int		altitude;
	long 	zone;		// fuseau horaire, en secondes - Paris=+1 en hiver -> +3600

	Instant	Inst;		// définition des dates et temps

	double	refraction;	// angle de réfraction

	double	temperature;	// temperature en °C
	double	pression;		// pression en hPa


	// données calculées
	double	eta1,eta2L,eta2C;	// pour le calcul des lever et coucher des astres
};
typedef	struct StructObservateur	Observateur;


unsigned char	ObservateurChargeLieu(Observateur *Obs);
void	ObservateurInitLieu(Observateur *Obs, wchar_t *n, double lng, double lat, int alt, long zone);
double	ObservateurInitTempsSyst(Observateur *Obs);
double	ObservateurInitTemps(Observateur *Obs, unsigned char jour, unsigned char mois, int annee, unsigned char heure, unsigned char min, unsigned char sec, unsigned char btu);
void	ObservateurAffiche(Observateur *Obs);
