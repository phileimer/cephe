#if TIME_WITH_SYS_TIME
# include <sys/time.h>
# include <time.h>
#else
# if HAVE_SYS_TIME_H
#  include <sys/time.h>
# else
#  include <time.h>
# endif
#endif

/* Pour i18n */
#define	__USE_UNIX98
//#define	__USE_ISOC99
#include <wchar.h>
#ifdef	HAVE_LIBINTL
	#include <libintl.h>
#endif	// HAVE_LIBINTL
//#define _(String) gettext(String)
#define _(String) String

#include "defs.h"
#include "data.h"

typedef	struct StructEphemerides	Ephemerides;
typedef	struct StructAstre	Astre;
typedef	struct StructAstreLune	AstreLune;

#include "vect.h"
#include "instant.h"
#include "orbit.h"
#include "body.h"
#include "observ.h"
#include "ephe.h"
