
# Testing the Docker run of PASA:



# PasaWeb using Docker

First, cd to the sample_data/ directory and then run:

docker run --rm -it -v /tmp:/tmp -v `pwd`:`pwd` -p 9000:9000 pasapipeline/pasapipeline:latest /usr/local/src/PASApipeline/run_PasaWeb.pl 9000


Then visit your http://${ip}:9000

to explore pasaweb.  You can choose a different port than 9000 above, such as 8080 or other.



    
