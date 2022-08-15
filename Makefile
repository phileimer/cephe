prog = $(CROSS_COMPILE)gcc
compflags = -c -O2 -Wall -DNO_CONFIG_H -DHAVE_LIBINTL
linkflags = -lm
debugflag =

exe = ephe
debugexe = ephe

objets = main.o vect.o body.o orbit.o observ.o instant.o ephe.o
inc = include.h data.h defs.h vect.h ephe.h orbit.h instant.h observ.h body.h


default: $(exe)
	$(CROSS_COMPILE)strip $(exe)

$(exe) : $(objets)
	$(prog) $^ $(linkflags) -o $@

$(objets) : %.o : %.c $(inc) 
	$(prog) $(compflags) $(debugflag) $< -o $@

distclean: clean
	rm -f $(exe) *~

clean:
	rm -f *.o
