#!/bin/bash

VERSION=`cat VERSION.txt`

docker build -t pasapipeline/pasapipeline:${VERSION} .
docker build -t pasapipeline/pasapipeline:latest .


