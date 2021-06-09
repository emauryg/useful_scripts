FROM ubuntu:18.04
LABEL junseok.park@childrens.harvard.edu

RUN apt-get update
RUN apt-get install -y git
RUN apt-get install -y wget apt-utils sudo
RUN apt-get install -y build-essential libcurl4-gnutls-dev libxml2-dev libssl-dev
RUN apt-get install -y zlib1g-dev vim

RUN wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
RUN bash Miniconda3-latest-Linux-x86_64.sh -b -p /miniconda
ENV PATH=$PATH:/miniconda/pcondabin:/miniconda/bin
COPY packages packages

# Tensorflow is installed in the base environment
# By default the python3 kernelspec points to a miniconda location
# /miniconda/share/jupyter/kernels/python3 and an ipykernel install installs the kernel in
# /usr/local/share/jupyter/kernels/python3. This can be solved by first deleting the kernelspec.
# This makes it point to /usr/local/share/jupyter/kernels/python3. A subsequent install then
# overwrites this one.
RUN conda env update -f packages/environment_tensorflow.yml
RUN jupyter kernelspec remove -f python3
RUN python -m ipykernel install --name python3 --display-name "Tensorflow"


RUN conda env update -f packages/environment_stea.yml
SHELL ["conda","run","-n","stea","/bin/bash","-c"]
RUN python -m ipykernel install --name python3 --display-name "stea"

RUN rm -rf packages
RUN rm Miniconda3-latest-Linux-x86_64.sh

RUN mkdir -pv /etc/ipython/
COPY ipython/ipython_config.py /etc/ipython/ipython_config.py

#ENV NB_USER leelab
#ENV NB_UID 1000

#RUN useradd -m -s /bin/bash -N -u $NB_UID $NB_USER
#RUN passwd -d $NB_USER
#RUN usermod -a -G sudo $NB_USER
#WORKDIR /home/${NB_USER}
#USER $NB_USER

RUN mkdir /workspace
WORKDIR /workspace


SHELL ["/bin/bash","-c"]
RUN conda init
RUN echo 'conda activate one' >> ~/.bashrc

EXPOSE 8888
#COPY README.md README.md

ENTRYPOINT ["jupyter", "notebook", "--no-browser","--ip=0.0.0.0","--NotebookApp.token=''","--NotebookApp.password=''", "--allow-root"]
