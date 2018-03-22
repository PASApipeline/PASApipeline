#!/bin/bash

VERSION=`cat VERSION.txt`

#docker build -t pasapipelinpasapipeline:${VERSION} .
docker build -t pasapipeline/pasapipeline:latest .


