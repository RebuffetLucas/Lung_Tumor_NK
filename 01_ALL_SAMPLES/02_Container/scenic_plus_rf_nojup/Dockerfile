FROM jupyter/minimal-notebook

# Init mamba and create environment
RUN mamba init bash && \
    /bin/bash -c "source ~/.bashrc" && \
    mamba create --name scenicplus python=3.11 -y

# Switch user to be able to install system dependencies
USER root
RUN apt-get update && apt-get install -y gcc g++
#RUN apt-get update && apt install -y default-jdk #Updated with the following line
RUN apt-get update && apt-get install -y gcc g++ make curl default-jdk libz-dev libbz2-dev libclang-dev liblzma-dev environment-modules flatbuffers-compiler

# Initialize the module system by adding it to bashrc
RUN echo 'source /etc/profile.d/modules.sh' >> /etc/bash.bashrc

# Install BEDTools
RUN curl -L https://github.com/arq5x/bedtools2/releases/download/v2.30.0/bedtools-2.30.0.tar.gz > BEDTools.tar.gz && \
    tar -zxvf BEDTools.tar.gz && \
    cd bedtools2 && \
    make && \
    cp bin/* /usr/local/bin/ && \
    cd .. && \
    rm -rf bedtools2 BEDTools.tar.gz

# Install flatbuffers Python package using conda run within the scenicplus environment
#RUN /opt/conda/bin/conda run -n scenicplus pip install flatbuffers

# Install flatbuffers from conda-forge
RUN pip install flatbuffers



# Allow starting jupyterlab when starting container as another user
RUN chmod -R 777 /home/jovyan

#Switch back to Jovyan
USER jovyan

# Download and install scenicplus
RUN git clone https://github.com/aertslab/scenicplus && \
    cd scenicplus && \
    pip install .

# Change location of caching for biomart database to temp directory to prevent database locked when running in container with limited permissions
RUN sed -i 's#.pybiomart#/tmp/.pybiomart#' /opt/conda/lib/python3.11/site-packages/pybiomart/base.py

#here we removed the start of the Jupyter for steps directly in cmd 
CMD ["/init"]
