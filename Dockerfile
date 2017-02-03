# This is the Dockefile to build PiReT (mshakya/PiReT)
# Base Docker Image
FROM ubuntu:16.04

# Maintainer
MAINTAINER Migun Shakya, migun@lanl.gov

# Update the system
RUN apt-get -y update
RUN apt-get -y install build-essential
RUN apt-get -y install git-all
RUN apt-get -y install wget
RUN apt-get clean

# copy local file
CMD cd /home/PiReT/
COPY scripts/ /home/PiReT/scripts/
COPY thirdParty/String-Approx-3.27.tar.gz /home/PiReT/thirdParty/
COPY bioconda_INSTALL.sh /home/PiReT/
COPY README.md /home/PiReT/
COPY runPiReT.pl /home/PiReT
COPY testdata/fastqs /home/PiReT/testdata/fastqs
COPY testdata/data /home/PiReT/testdata/data
COPY testdata/test_pipeline_linux.sh /home/PiReT/testdata/test_pipeline_linux.sh

CMD sh ./bioconda_INSTALL.sh
