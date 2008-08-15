
GFORTRAN=gfortran
FFLAGS=-Wall
PROG=waveref
VERSION=0.1
SOURCES=wave_refraction.f
OBJS=${SOURCES:.f=.o}

all: ${PROG}

${PROG}: ${SOURCES}
	${GFORTRAN} ${FFLAGS} -o ${PROG} ${SOURCES}

dist:
	@mkdir ${PROG}-${VERSION}
	@cp ${SOURCES} ${PROG}-${VERSION}
	@cp Makefile ${PROG}-${VERSION}
	@tar cvfz ${PROG}-${VERSION}.tar.gz ${PROG}-${VERSION} 
	@rm -rf ${PROG}-${VERSION}

clean:
	@rm -f ${PROG} ${OBJS} core

