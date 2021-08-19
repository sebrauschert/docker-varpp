FROM rocker/tidyverse:latest

#=================
# Set conda variables: This is required to make sure when we call 'conda',
# the container knows where to look.
ENV PATH="/root/miniconda3/bin:${PATH}"
ARG PATH="/root/miniconda3/bin:${PATH}"

#=================
# Setup: jsut to make sure we hace the tools available
RUN apt-get update && \
      apt-get -y install sudo && \
      apt-get -y install make

# wget
RUN sudo apt-get -y install wget

#=================
# Create directories
RUN sudo mkdir /root/dependencies
RUN sudo mkdir /root/varpp-prediction

#=================
# Copy over the dependencies
ADD header.txt /root/dependencies/header.txt
ADD predict.R /root/varpp-prediction/predict.R
ADD varppRule.R /root/varpp-prediction/varppRule.R
ADD varppRuleAndPredict.R /root/varpp-prediction/varppRuleAndPredict.R
RUN sudo chmod +x /root/varpp-prediction/predict.R
RUN sudo chmod +x /root/varpp-prediction/varppRule.R
RUN sudo chmod +x /root/varpp-prediction/varppRuleAndPredict.R

#=================
# Add predict.R to the path
ENV PATH=/root/varpp-prediction:$PATH
ENV PATH=/root/varpp-prediction:$PATH

#=================
# Install R packages
RUN mkdir -p /opt/software/setup/R
ADD install_packages.R /opt/software/setup/R
RUN Rscript /opt/software/setup/R/install_packages.R

#=================
# Install miniconda

RUN wget \
    https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh \
    && mkdir /root/.conda \
    && bash Miniconda3-latest-Linux-x86_64.sh -b \
    && rm -f Miniconda3-latest-Linux-x86_64.sh 

#============================================
# INSTALL BEDTOOLS, BCFTOOLS AND PARALLEL
#============================================

#=================
# Install bedtools
RUN conda install -c bioconda bedtools

#=================
# Install bcftools
RUN conda install -c bioconda bcftools

#=================
# Install GNU parallel
RUN wget http://ftp.gnu.org/gnu/parallel/parallel-latest.tar.bz2
RUN sudo tar xjf parallel-latest.tar.bz2
RUN cd parallel-*; sudo ./configure && make; sudo make install

