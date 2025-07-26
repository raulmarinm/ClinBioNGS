ARG R_VERSION=4.0.4

FROM ubuntu:22.04 as build
ARG R_VERSION

RUN DEBIAN_FRONTEND=noninteractive apt-get update && \
    apt-get install --no-install-recommends -y \
        wget \
        ca-certificates \
        cmake \
        libssl-dev \
        libxml2-dev

WORKDIR /tmp

RUN wget https://cran.r-project.org/src/base/R-4/R-${R_VERSION}.tar.gz && \
    tar xzf R-${R_VERSION}.tar.gz

RUN sed -i "/^#.*deb-src.*universe$/s/^# //g" /etc/apt/sources.list

RUN apt-get update && \
    DEBIAN_FRONTEND=noninteractive TZ="Europe/Madrid" \
    apt-get -y --no-install-recommends install tzdata && \
    apt-get -y build-dep r-base

WORKDIR /tmp/R-${R_VERSION}

RUN ./configure \
    --prefix=/opt/R \
    --enable-R-shlib \
    --enable-memory-profiling \
    --with-jpeglib \
    --with-readline \
    --with-pcre1 \
    --with-blas \
    --with-lapack \
    --with-cairo \
    --without-x
RUN make && make install

WORKDIR /tmp
RUN \
    /opt/R/bin/R -e 'install.packages(c("cli", "farver", "glue", "labeling", "lifecycle", "munsell", "R6", "RColorBrewer", "rlang", "viridisLite"), repos ="https://cloud.r-project.org")' && \
    wget https://cran.r-project.org/src/contrib/Archive/scales/scales_1.3.0.tar.gz && \
    /opt/R/bin/R CMD INSTALL scales_1.3.0.tar.gz && \
    wget https://cran.r-project.org/src/contrib/Archive/Matrix/Matrix_1.6-5.tar.gz && \
    /opt/R/bin/R CMD INSTALL Matrix_1.6-5.tar.gz
RUN /opt/R/bin/R -e 'install.packages("tidyverse", repos ="https://cloud.r-project.org")'

RUN \
    /opt/R/bin/R -e 'install.packages(c("checkmate", "Formula", "gridExtra", "htmlwidgets", "latticeExtra", "viridis"), repos ="https://cloud.r-project.org")' && \
    wget https://cran.r-project.org/src/contrib/Archive/htmlTable/htmlTable_2.4.1.tar.gz && \
    /opt/R/bin/R CMD INSTALL htmlTable_2.4.1.tar.gz && \
    wget https://cran.r-project.org/src/contrib/Archive/Hmisc/Hmisc_4.6-0.tar.gz && \
    /opt/R/bin/R CMD INSTALL Hmisc_4.6-0.tar.gz
RUN /opt/R/bin/R -e 'install.packages("BiocManager", repos ="https://cloud.r-project.org")'
RUN /opt/R/bin/R -e 'BiocManager::install("CopyNumberPlots")'

RUN \
    /opt/R/bin/R -e 'install.packages("lme4", repos ="https://cloud.r-project.org")' && \
    wget https://cran.r-project.org/src/contrib/Archive/numDeriv/numDeriv_2016.8-1.tar.gz && \
    wget https://cran.r-project.org/src/contrib/Archive/pbkrtest/pbkrtest_0.5.1.tar.gz && \
    /opt/R/bin/R CMD INSTALL numDeriv_2016.8-1.tar.gz && \
    /opt/R/bin/R CMD INSTALL pbkrtest_0.5.1.tar.gz && \
    /opt/R/bin/R -e 'install.packages("ggpubr", repos ="https://cloud.r-project.org")'

RUN /opt/R/bin/R -e 'install.packages("cowplot", repos ="https://cloud.r-project.org")'
RUN /opt/R/bin/R -e 'install.packages("crosstalk", repos ="https://cloud.r-project.org")'
RUN /opt/R/bin/R -e 'install.packages("DT", repos ="https://cloud.r-project.org")'
RUN /opt/R/bin/R -e 'install.packages("flexdashboard", repos ="https://cloud.r-project.org")'

RUN \
    wget https://cran.r-project.org/src/contrib/Archive/yulab.utils/yulab.utils_0.1.3.tar.gz && \
    /opt/R/bin/R CMD INSTALL yulab.utils_0.1.3.tar.gz && \
    /opt/R/bin/R -e 'install.packages("ggplotify", repos ="https://cloud.r-project.org")'

RUN \
    wget https://cran.r-project.org/src/contrib/Archive/HGNChelper/HGNChelper_0.8.1.tar.gz && \
    /opt/R/bin/R CMD INSTALL HGNChelper_0.8.1.tar.gz 

RUN \
    wget https://cran.r-project.org/src/contrib/Archive/svglite/svglite_2.1.3.tar.gz && \
    /opt/R/bin/R CMD INSTALL svglite_2.1.3.tar.gz && \
    /opt/R/bin/R -e 'install.packages("kableExtra", repos ="https://cloud.r-project.org")'

RUN /opt/R/bin/R -e 'install.packages("openxlsx", repos ="https://cloud.r-project.org")'
RUN /opt/R/bin/R -e 'install.packages("R.utils", repos ="https://cloud.r-project.org")'
RUN /opt/R/bin/R -e 'install.packages("plyr", repos ="https://cloud.r-project.org")'


FROM ubuntu:22.04

ENV DEBIAN_FRONTEND="noninteractive"
ENV TZ="Europe/Madrid"

RUN apt-get update && \
    apt-get install --no-install-recommends -y \
        pandoc \
        libblas-dev \
        liblapack-dev \
        libreadline-dev \
        libicu-dev \
        libgsl-dev \
        libudunits2-dev \
        libbz2-dev \
        libpcre2-dev \
        libgomp1 \
        libxml2 \
        libcurl4-openssl-dev \
        libfontconfig1-dev \
        libjpeg-dev \
        libtiff-dev \
        libpangocairo-1.0.0 \
        xz-utils \
        tzdata \
        locales && \
    locale-gen en_US.UTF-8 && \
    rm -fr /var/lib/apt/lists/*

COPY --from=build /opt/R/ /opt/R/

RUN ln -s /opt/R/bin/R /usr/bin/R && \
    ln -s /opt/R/bin/Rscript /usr/bin/Rscript

ENV LANG="en_US.UTF-8"
ENV LC_ALL="en_US.UTF-8"
ENV LANGUAGE="en_US.UTF-8"

CMD ["R"]
