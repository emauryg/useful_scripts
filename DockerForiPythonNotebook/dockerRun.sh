
workDirectory=$(dirname $PWD)

docker run -it -p 8888:8888 --cpus=3 --memory=8GB -w /home/leelab -v ${workDirectory}:/home/leelab leelab_dev
