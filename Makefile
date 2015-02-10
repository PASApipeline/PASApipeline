
all: 
	if [ ! -d bin ]; then mkdir bin; fi
	cd pasa_cpp && $(MAKE) && cp pasa ../bin/.
	cd slclust && $(MAKE) && cp src/slclust ../bin/.
	cd cdbtools/cdbfasta && $(MAKE) && cp cdbfasta ../../bin/. && cp cdbyank ../../bin/.

clean:
	cd pasa_cpp && $(MAKE) clean
	cd slclust && $(MAKE) clean
	cd cdbtools/cdbfasta && $(MAKE) clean
	rm -f bin/* 

###################################################################


