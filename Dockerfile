# Docker image for a cawlign development environment
FROM oraclelinux:8

# Set up environment and install dependencies
RUN yum -y update && \
    #yum install -y cmake gcc-c++ gcc-toolset-10 git make oracle-epel-release-el8 && \
    #echo 'source /opt/rh/gcc-toolset-10/enable' > ~/.bashrc && \
    #source ~/.bashrc
    yum install -y cmake gcc-c++ git make

# To compile cawlign within the development environment:
#   cmake .
#   make
