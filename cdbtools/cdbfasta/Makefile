# Useful directories

TGICLASSDIR := ../tgi_cl
# Directories to search for header files
SEARCHDIRS := -I${TGICLASSDIR}


SYSTYPE :=     $(shell uname)

# C compiler

CC      := g++
CFLAGS  = -Wall ${SEARCHDIRS} -O2 -D_FILE_OFFSET_BITS=64 \
-D_LARGEFILE_SOURCE -fno-exceptions -fno-rtti -D_REENTRANT 

%.o : %.c
	${CC} ${CFLAGS} -c $< -o $@

%.o : %.cc
	${CC} ${CFLAGS} -c $< -o $@

%.o : %.C
	${CC} ${CFLAGS} -c $< -o $@

%.o : %.cpp
	${CC} ${CFLAGS} -c $< -o $@

%.o : %.cxx
	${CC} ${CFLAGS} -c $< -o $@

# C/C++ linker

LINKER    := g++
LDFLAGS    = -lz

.PHONY : all
all:    cdbfasta cdbyank 

cdbfasta:  ./cdbfasta.o ./gcdbz.o $(LIBS) ${TGICLASSDIR}/gcl/gcdb.o ${TGICLASSDIR}/gcl/GReadBuf.o ${TGICLASSDIR}/gcl/GBase.o ${TGICLASSDIR}/gcl/GString.o ${TGICLASSDIR}/gcl/GArgs.o
	${LINKER} -o $@ ${filter-out %.a %.so, $^} $(LDFLAGS)

cdbyank :  ./cdbyank.o ./gcdbz.o ${TGICLASSDIR}/gcl/gcdb.o ${TGICLASSDIR}/gcl/GReadBuf.o ${TGICLASSDIR}/gcl/GBase.o ${TGICLASSDIR}/gcl/GString.o ${TGICLASSDIR}/gcl/GArgs.o
	${LINKER} -o $@ ${filter-out %.a %.so, $^} $(LDFLAGS)

# target for removing all object files

.PHONY : tidy
tidy::
	@${RM} core cdbfasta cdbyank *.o ${TGICLASSDIR}/gcl/gcdb.o ${TGICLASSDIR}/gcl/GReadBuf.o ${TGICLASSDIR}/gcl/GBase.o ${TGICLASSDIR}/gcl/GString.o ${TGICLASSDIR}/gcl/GArgs.o

# target for removing all object files

.PHONY : clean
clean:: tidy
	@${RM} core cdbfasta cdbyank *.o ${TGICLASSDIR}/gcl/gcdb.o ${TGICLASSDIR}/gcl/GReadBuf.o ${TGICLASSDIR}/gcl/GBase.o ${TGICLASSDIR}/gcl/GString.o ${TGICLASSDIR}/gcl/GArgs.o


