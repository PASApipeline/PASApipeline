FROM ubuntu:16.04
MAINTAINER bhaas@broadinstitute.org

RUN apt-get update && apt-get install -y gcc g++ perl python automake make \
                                       wget git curl libdb-dev \
                                       zlib1g-dev bzip2 libncurses5-dev \
                                       texlive-latex-base \
                                       default-jre \
                                       python-pip python-dev \
                                       gfortran \
                                       build-essential libghc-zlib-dev libncurses-dev libbz2-dev liblzma-dev libpcre3-dev libxml2-dev \
                                       libblas-dev gfortran git unzip ftp libzmq3-dev nano ftp fort77 libreadline-dev \
                                       libcurl4-openssl-dev libx11-dev libxt-dev \
                                       x11-common libcairo2-dev libpng12-dev libreadline6-dev libjpeg8-dev pkg-config libtbb-dev \
                   && apt-get clean

RUN curl -L https://cpanmin.us | perl - App::cpanminus

RUN cpanm install DB_File

RUN cpanm install URI::Escape


## set up tool config and deployment area:

ENV SRC /usr/local/src
ENV BIN /usr/local/bin

RUN apt-get install -y sqlite lighttpd libgd-tools libgd2-xpm-dev

RUN cpanm install GD
RUN cpanm install DBI
RUN cpanm install DBD::SQLite

RUN apt-get install -y software-properties-common

RUN add-apt-repository 'deb http://archive.ubuntu.com/ubuntu trusty universe' && \
    apt-get update && \
    apt install -y mysql-server-5.6 && \
    apt install -y mysql-client-5.6 && \
    apt-get install -y libdbd-mysql-perl



## GMAP installation
WORKDIR $SRC
RUN wget http://research-pub.gene.com/gmap/src/gmap-gsnap-2017-11-15.tar.gz && \
        tar xvf gmap-gsnap-2017-11-15.tar.gz && \
        cd gmap-2017-11-15 && \
        ./configure && \
        make && \
        make install

## BLAT
WORKDIR $BIN
RUN wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/blat/blat && \
        chmod 755 ./blat


## Fasta3
WORKDIR $SRC
RUN wget http://faculty.virginia.edu/wrpearson/fasta/fasta36/fasta-36.3.8g.tar.gz && \
        tar zxvf fasta-36.3.8g.tar.gz && \
        cd ./fasta-36.3.8g/src && \
        make -f ../make/Makefile.linux_sse2 all && \
        cp ../bin/fasta36 /usr/local/bin/fasta

       
## PASA installation
WORKDIR $SRC

ENV PASA_CO a6965a3

RUN git clone https://github.com/PASApipeline/PASApipeline.git && \
    cd PASApipeline && \
    git checkout $PASA_CO && \
    git submodule init && git submodule update && \
    make

ENV PASAHOME /usr/local/src/PASApipeline



