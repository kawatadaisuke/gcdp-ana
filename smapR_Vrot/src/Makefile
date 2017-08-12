EXEC=smrvrot
OBJS = module.o main.o rdata.o para.o kernel.o allocate.o setsval.o
DEFFILE=module.f90
OPTIONS= -frecord-marker=4 -g
### for MPI
# mac
MPI_DIR = /usr/local
#  kepler
# MPI_DIR=/cluster/mpich2
FORT = ${MPI_DIR}/bin/mpif90
MPI_INC = ${MPI_DIR}/include
MPI_LIB = ${MPI_DIR}/lib
FFLAGS = ${OPTIONS} -I${MPI_INC}
LIBS = ${MLIBS} -L${MPI_LIB} ${OPTIONS}
# xt4
#FORT=ftn
#OPTIONS= -fastsse -O4 -Mvect=sse
FFLAGS = ${OPTIONS}
LIBS=${LINK_MPI}
.SUFFIXES : .o .f90

${EXEC} : ${OBJS}
	${FORT} -o ${EXEC} ${OBJS} ${LIBS}
	mv ${EXEC} ..

.f90.o :
	${FORT} ${FFLAGS} -c $<

${OBJS} : ${DEFFILE}
%.o : %.mod

clean :
	rm *.o ${EXEC} *.mod

