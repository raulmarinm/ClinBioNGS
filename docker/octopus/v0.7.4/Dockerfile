ARG BOOST_VERSION=1.76.0
ARG BOOST_VERSION_=1_76_0
ARG HTSLIB_VERION=1.12
ARG GMPLIB_VERSION=6.3.0
ARG OCTOPUS_VERION=0.7.4

FROM ubuntu:22.04 AS build
ARG BOOST_VERSION
ARG BOOST_VERSION_
ARG HTSLIB_VERION
ARG GMPLIB_VERSION
ARG OCTOPUS_VERION

RUN apt-get update && \
    apt-get install --no-install-recommends -y \
        g++ \
        gcc \
        make \
        cmake \
        ca-certificates \
        wget \
        bzip2 \
        libbz2-dev \
        libzip-dev \
        liblzma-dev \
        libcurl4-openssl-dev \
        libssl-dev \
        xz-utils \
        m4 \
        git && \
    rm -rf /var/lib/apt/lists/*

# Compile HTSlib
WORKDIR /tmp
RUN wget https://github.com/samtools/htslib/releases/download/${HTSLIB_VERION}/htslib-${HTSLIB_VERION}.tar.bz2 && \
    tar jxf htslib-${HTSLIB_VERION}.tar.bz2
WORKDIR /tmp/htslib-${HTSLIB_VERION}
RUN ./configure && \
    make && \
    make install && \
    cd .. && \
    rm -fr htslib-${HTSLIB_VERION}*

# Compile GMPlib
WORKDIR /tmp
RUN wget https://gmplib.org/download/gmp/gmp-${GMPLIB_VERSION}.tar.xz && \
    tar -xf gmp-${GMPLIB_VERSION}.tar.xz
WORKDIR /tmp/gmp-${GMPLIB_VERSION}/
RUN ./configure && \
    make && \
    make check && \
    make install && \
    cd .. && \
    rm -fr gmp-${GMPLIB_VERSION}*

# Compile Boost
WORKDIR /tmp
RUN wget https://archives.boost.io/release/${BOOST_VERSION}/source/boost_${BOOST_VERSION_}.tar.gz && \
    tar xzf boost_${BOOST_VERSION_}.tar.gz
WORKDIR /tmp/boost_${BOOST_VERSION_}
RUN ./bootstrap.sh --with-libraries=all --prefix=/tmp/boost/ && \
    ./b2 --prefix=/tmp/boost/ install && \
    cd .. && \
    rm -fr boost_${BOOST_VERSION_}*

# Compile Octopus
WORKDIR /tmp
RUN git clone -b master https://github.com/luntergroup/octopus.git && \
    cd octopus && \
    git checkout v${OCTOPUS_VERION}
WORKDIR /tmp/octopus/build
RUN cmake \
        -DBOOST_ROOT=/tmp/boost \
        -DBoost_INCLUDE_DIR=/tmp/boost/include \
        -DBOOST_LIBRARYDIR=/tmp/boost/lib \
        -DGMP_LIBRARIES=/usr/local/lib/libgmp.so \
        -DGMP_INCLUDES=/usr/local/include \
        -DHTSlib_INCLUDE_DIR=/usr/local/include \
        -DHTSlib_LIBRARY=/usr/local/lib/libhts.so \
        .. && \
    make install

FROM ubuntu:22.04
LABEL maintainer="ODAP Administrators <odap_admin@iconcologia.net>"
RUN apt-get update && \
    apt-get install --no-install-recommends -y libcurl4-openssl-dev && \
    rm -rf /var/lib/apt/lists/*

COPY --from=build /tmp/octopus/bin/octopus /usr/bin/octopus
COPY --from=build /usr/local/lib/libgmp.so* /usr/local/lib/libhts.so* /tmp/boost/lib/lib* /usr/local/lib/
ENV LD_LIBRARY_PATH=/usr/local/lib/
