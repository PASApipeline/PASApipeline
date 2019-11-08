
## If wanting to build on mac osx using gcc instead of clang:
## run make like so:
##
##   make CC=gcc CXX=g++

all: 
	if [ ! -d bin ]; then mkdir bin; fi
	cd pasa_cpp && $(MAKE) && cp pasa ../bin/.
	cd pasa-plugins/slclust && $(MAKE) && cp src/slclust ../../bin/.
	cd pasa-plugins/cdbtools/cdbfasta && $(MAKE) && cp cdbfasta ../../../bin/. && cp cdbyank ../../../bin/.
	cd pasa-plugins/seqclean/mdust && $(MAKE) && cp mdust ../../../bin
	cd pasa-plugins/seqclean/psx && $(MAKE) && cp psx ../../../bin
	cd pasa-plugins/seqclean/trimpoly && $(MAKE) && cp trimpoly ../../../bin
	cp pasa-plugins/seqclean/seqclean/seqclean pasa-plugins/seqclean/seqclean/cln2qual pasa-plugins/seqclean/seqclean/bin/seqclean.psx ./bin

clean:
	cd pasa_cpp && $(MAKE) clean
	cd pasa-plugins/slclust && $(MAKE) clean
	cd pasa-plugins/cdbtools/cdbfasta && $(MAKE) clean
	cd pasa-plugins/seqclean/mdust && $(MAKE) clean
	cd pasa-plugins/seqclean/psx && $(MAKE) clean
	cd pasa-plugins/seqclean/trimpoly && $(MAKE) clean
	rm -f bin/* 

###################################################################


