# Docker image for a cawlign development environment
FROM oraclelinux:8

# Set up environment and install dependencies
RUN yum -y update && \
    yum install -y cmake gcc-c++ git make

# To compile cawlign within the development environment:
#   cmake .
#   make
