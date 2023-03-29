#!/bin/bash

VERSION=`cat VERSION.txt`

singularity build pasapipeline.v${VERSION}.simg docker://pasapipeline/pasapipeline:${VERSION}

#singularity exec -e trinityrnaseq.v${VERSION}.simg Trinity --version

