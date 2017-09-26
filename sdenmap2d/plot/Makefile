FORT = gfortran
EXEC = 2dp 
OBJS = main.o 
EXECVM = pvmet
OBJSVM = main-vmetmap.o
EXECVMN = pvmetn6754
OBJSVMN = main-vmetmapn6754.o
MLIBS = -lm
OPTIONS = -g 
FFLAGD = -fno-backslash
FFLAGS = ${OPTIONS} ${FFLAGD}
# mac
LIBS = ${MLIBS} -L/usr/local/pgplot -lpgplot -L/usr/X11/lib -lX11
# kepler
# LIBS = ${MLIBS} /khome/dka/share/pgplot/libpgplot.a -L/usr/lib -lX11
.SUFFIXES : .o .f90

${EXEC} : ${OBJS}
	${FORT} ${FFLAGD} -o ${EXEC} ${OBJS} ${LIBS}
	\mv ${EXEC} ..

${EXECVM} : ${OBJSVM}
	${FORT} ${FFLAGD} -o ${EXECVM} ${OBJSVM} ${LIBS}
	\mv ${EXECVM} ..

${EXECVMN} : ${OBJSVMN}
	${FORT} ${FFLAGD} -o ${EXECVMN} ${OBJSVMN} ${LIBS}
	\mv ${EXECVMN} ..


.f90.o :
	${FORT} ${FFLAGS} -c $<

clean :
	rm *.o
