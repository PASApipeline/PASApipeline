#!/bin/bash

VERSION=`cat VERSION.txt`

docker push pasapipeline/pasapipeline:${VERSION} 
docker push pasapipeline/pasapipeline:latest 


