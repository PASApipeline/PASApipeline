FROM ubuntu:20.04
LABEL maintainer="bhaas@broadinstitute.org"

ARG DEBIAN_FRONTEND=noninteractive
ENV LANG C.UTF-8
ENV TZ 'America/Indiana/Indianapolis'
RUN ln -snf /usr/share/zoneinfo/$TZ /etc/localtime && echo $TZ > /etc/timezone  && \
    apt-get update && \
    apt-get install --yes --no-install-recommends ca-certificates netbase \
    samtools libdb-dev wget curl git \
    build-essential gcc-multilib \
    pkg-config \
    libgd-gd2-perl \
    libgd-dev \
    gcc g++ make cpanminus \
    sqlite3 mysql-server mysql-client libdbd-mysql-perl libssl-dev \
    software-properties-common && \
    rm -rf /var/lib/apt/lists/*

## set up tool config and deployment area:
ENV SRC /usr/local/src
ENV BIN /usr/local/bin

RUN cpanm install DB_File \
                  URI::Escape \
		              GD \
    	  	        DBI \
		              DBD::SQLite

## GMAP installation
WORKDIR $SRC
ENV GSNAP_VER 2021-07-23
RUN wget http://research-pub.gene.com/gmap/src/gmap-gsnap-${GSNAP_VER}.tar.gz && \
        tar xvf gmap-gsnap-${GSNAP_VER}.tar.gz && \
        cd gmap-${GSNAP_VER} && \
        ./configure && \
        make && \
        make install && \
	rm $SRC/gmap-gsnap-${GSNAP_VER}.tar.gz

## BLAT
WORKDIR $BIN
RUN wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/blat/blat && \
        chmod 755 ./blat

## Fasta3
WORKDIR $SRC
RUN wget https://github.com/wrpearson/fasta36/archive/refs/tags/v36.3.8h_04-May-2020.tar.gz && \
        tar zxvf v36.3.8h_04-May-2020.tar.gz && \
        cd ./fasta36-36.3.8h_04-May-2020/src && \
        make -f ../make/Makefile.linux_sse2 all && \
        cp ../bin/fasta36 /usr/local/bin/fasta && \
	rm $SRC/v36.3.8h_04-May-2020.tar.gz

## PASA installation
WORKDIR $SRC

## pasapipeline v2.5.0-dev
ENV PASA_CO e68b3d7f036f8d01b10bc239cde199b0c9ad8c25

RUN git clone https://github.com/PASApipeline/PASApipeline.git && \
    cd PASApipeline && \
    git checkout $PASA_CO && \
    git submodule init && git submodule update && \
    make


ENV PASAHOME /usr/local/src/PASApipeline

RUN cp $PASAHOME/bin/* /usr/local/bin/
COPY conf.txt /usr/local/src/PASApipeline/pasa_conf/

