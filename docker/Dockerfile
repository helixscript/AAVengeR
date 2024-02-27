
#Docker version 20.10.17, build 100c701
# build on the r-base docker image 
FROM r-base

USER root
# create and set the main directory structure 
WORKDIR /home

# copy files from the local directory to the image WORKDIR
COPY . .

ENV CONDA_DIR /opt/conda

RUN wget --quiet https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O ~/miniconda.sh && \
	 /bin/bash ~/miniconda.sh -b -p /opt/conda

# # Put conda in path so we can use conda activate
ENV PATH=$CONDA_DIR/bin:$PATH

RUN conda env create -f python2_env.yml

RUN conda install -n python2_env numpy=1.16.6

RUN conda env create -f bioinformatic_dependencies.yml

RUN conda env create -f blat_dependency.yml

RUN apt-get update

RUN apt-get -y install libcurl4-openssl-dev

RUN apt-get -y install libssl-dev

RUN apt-get -y install libxml2-dev

RUN apt-get -y install libgmp3-dev

RUN apt-get -y install libmariadb-dev

RUN apt-get -y install libmysqlclient-dev

RUN Rscript ./R_dependencies.R

RUN echo "export PATH=$PATH:/opt/conda/envs/bioinformatic_dependencies/bin:/opt/conda/envs/python2_env/bin:/opt/conda/envs/blat_dependency/bin" >> /root/.bashrc

RUN echo "alias python2=/opt/conda/envs/python2_env/bin/python" >> /root/.bashrc 

RUN echo "ln -s /usr/lib/x86_64-linux-gnu/libmariadb.so.3 /usr/lib/x86_64-linux-gnu/libmysqlclient.so.18" >> /root/.bashrc

#https://stackoverflow.com/questions/64323539/docker-run-interactive-with-conda-environment-already-activated
