FROM biocontainers/samtools:1.3 as build

WORKDIR /app

ENV PISCES_VERSION=5.3.0.0

RUN wget https://github.com/tamsen/Pisces/releases/download/v5.3.0.0/pisces_all_${PISCES_VERSION}.tar.gz

RUN tar --no-same-owner -xzvf pisces_all_${PISCES_VERSION}.tar.gz

RUN find pisces_all -type f -perm 764 -exec chmod 755 {} \;

RUN mv pisces_all /app/Pisces_v${PISCES_VERSION}

WORKDIR /app/Pisces_v${PISCES_VERSION}

ENV PATH="/app/Pisces_v${PISCES_VERSION}:$PATH"

CMD [ "Pisces" ]