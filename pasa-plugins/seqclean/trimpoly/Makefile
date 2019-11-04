# Useful directories

#CC      := gcc
CFLAGS  = -O2 -Wall

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

LINKER    := gcc
LDFLAGS    = 
LOADLIBES =
SYSTYPE :=     $(shell uname)
#ifeq ($(SYSTYPE), Linux)
#  LOADLIBES = -lpvm3
#else #some rigid Unix implementation, like SunOS
#  LOADLIBES = -lpvm3 -lsocket -lnsl
#endif

.PHONY : all
all:    trimpoly

trimpoly:  ./trimpoly.o color_c.o
	${CC} ${LDFLAGS} -o $@ ${filter-out %.a %.so, $^} ${LOADLIBES}

# target for removing all object files

.PHONY : tidy
tidy::
	@${RM} core *.o 

# target for removing all object files

.PHONY : clean
clean:: tidy
	@${RM} core *.o trimpoly


