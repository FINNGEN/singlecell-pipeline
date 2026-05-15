FROM europe-docker.pkg.dev/finngen-refinery-dev/eu.gcr.io/singlecell-suite:1.84
LABEL maintainer="mkanai@broadinstitute.org"

ENV SNAPATAC2_VERSION=2.8.0
RUN pip3 install --break-system-packages snapatac2==${SNAPATAC2_VERSION}
