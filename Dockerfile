# syntax=docker/dockerfile:1

FROM debian:trixie-slim AS builder

RUN apt-get update && apt-get install -y --no-install-recommends \
        ca-certificates \
        git \
        build-essential \
        cmake \
        pkg-config \
        libhts-dev \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /src
COPY . .

RUN cmake -S . -B build \
        -DCMAKE_BUILD_TYPE=RelWithDebInfo \
        -DMAKE_MAIN=ON \
        -DMAKE_DAUGHTER=ON \
    && cmake --build build --parallel "$(nproc)" \
    && install -m 0755 build/expos /usr/local/bin/expos \
    && install -m 0755 build/estimate-entropy /usr/local/bin/estimate-entropy


FROM debian:trixie-slim AS runtime

RUN apt-get update && apt-get install -y --no-install-recommends \
        ca-certificates \
        libhts3t64 \
    && update-ca-certificates \
    && rm -rf /var/lib/apt/lists/*

COPY --from=builder /usr/local/bin/expos /usr/local/bin/expos
COPY --from=builder /usr/local/bin/estimate-entropy /usr/local/bin/estimate-entropy

ENTRYPOINT ["/usr/local/bin/expos"]
