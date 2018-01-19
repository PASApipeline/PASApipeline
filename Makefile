
all: 
	if [ ! -d bin ]; then mkdir bin; fi
	cd pasa_cpp && $(MAKE) && cp pasa ../bin/.
	cd slclust && $(MAKE) && cp src/slclust ../bin/.
	cd cdbtools/cdbfasta && $(MAKE) && cp cdbfasta ../../bin/. && cp cdbyank ../../bin/.
	cd seqclean/mdust && $(MAKE) && cp mdust ../../bin
	cd seqclean/psx && $(MAKE) && cp psx ../../bin
	cd seqclean/trimpoly && $(MAKE) && cp trimpoly ../../bin
	cp seqclean/seqclean/seqclean seqclean/seqclean/cln2qual seqclean/seqclean/bin/seqclean.psx ./bin

clean:
	cd pasa_cpp && $(MAKE) clean
	cd slclust && $(MAKE) clean
	cd cdbtools/cdbfasta && $(MAKE) clean
	cd seqclean/mdust && $(MAKE) clean
	cd seqclean/psx && $(MAKE) clean
	cd seqclean/trimpoly && $(MAKE) clean
	rm -f bin/* 

###################################################################


