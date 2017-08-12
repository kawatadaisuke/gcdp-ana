FORT = gfortran
EXEC = plotrvrot
OBJS = main.o 
MLIBS = -lm
OPTIONS = -g 
FFLAGD = -fno-backslash
FFLAGS = ${OPTIONS} ${FFLAGD}
# mac
LIBS = ${MLIBS} -L/usr/local/pgplot -lpgplot -L/usr/X11R6/lib -lX11
# kepler
# LIBS = ${MLIBS} /khome/dka/share/pgplot/libpgplot.a -L/usr/lib -lX11
.SUFFIXES : .o .f90

${EXEC} : ${OBJS}
	${FORT} ${FFLAGD} -o ${EXEC} ${OBJS} ${LIBS}
	mv ${EXEC} ..

.f90.o :
	${FORT} ${FFLAGS} -c $<

clean :
	rm *.o
