
# build on the r-base docker image 
FROM r-base

USER root

# create and set the main directory structure 
WORKDIR /home

ENV CONDA_DIR /opt/conda

RUN wget --quiet https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O ~/miniconda.sh && \
	 /bin/bash ~/miniconda.sh -b -p /opt/conda

# # Put conda in path so we can use conda activate
ENV PATH=$CONDA_DIR/bin:$PATH

# copy files from the local directory to the image WORKDIR
COPY . .

RUN conda env create -f python2_env.yml

RUN conda env create -f bioinformatic_dependencies.yml

RUN apt-get update

RUN apt-get -y install libcurl4-openssl-dev

RUN apt-get -y install libssl-dev

RUN apt-get -y install libxml2-dev

RUN apt-get -y install libgmp3-dev

RUN apt-get -y install libmariadb-dev

# Verify R dependencies were installed successfully

RUN Rscript ./R_dependencies.R

# Add to paths

RUN echo "export PATH=$PATH:/opt/conda/envs/bioinformatic_dependencies/bin:/opt/conda/envs/python2_env/bin" >> /root/.bashrc

RUN echo "alias python2=/opt/conda/envs/python2_env/bin/python" >> /root/.bashrc 

RUN source /root/.bashrc

RUN source activate bioinformatic_dependencies


# CMD ["Rscript", "myscript.R"]

# allows for variables to be supplied to main program
# ENTRYPOINT ["Rscript", "myscript.R"]