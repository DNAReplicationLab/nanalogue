FROM rust:1.99.0-slim-trixie AS build-env
ARG RUST_TARGET=x86_64-unknown-linux-musl
RUN rustup target add ${RUST_TARGET}
# Use matching Clang/libclang 19 for hts-sys/bindgen compatibility.
RUN apt-get update \
    && apt-get install -y --no-install-recommends \
        musl musl-tools musl-dev build-essential clang-19 libclang-19-dev \
        ca-certificates pkg-config perl \
    && rm -rf /var/lib/apt/lists/*
ENV CLANG_PATH=/usr/bin/clang-19
ENV LIBCLANG_PATH=/usr/lib/llvm-19/lib
# Create symlink for ARM64 musl compiler if building for aarch64
RUN if [ "${RUST_TARGET}" = "aarch64-unknown-linux-musl" ]; then \
        ln -s /usr/bin/musl-gcc /usr/local/bin/aarch64-linux-musl-gcc; \
    fi
WORKDIR /app
COPY . /app
RUN RUSTFLAGS='-C link-arg=-s' cargo build --release --target ${RUST_TARGET}

FROM gcr.io/distroless/static-debian12
ARG RUST_TARGET=x86_64-unknown-linux-musl
COPY --from=build-env /app/target/${RUST_TARGET}/release/nanalogue /
COPY --from=build-env /app/target/${RUST_TARGET}/release/nanalogue_sim_bam /
ENV PATH="$PATH:/"
CMD ["nanalogue"]
