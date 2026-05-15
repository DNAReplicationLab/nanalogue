FROM rust:1.92.0-bullseye as build-env
ARG RUST_TARGET=x86_64-unknown-linux-musl
RUN rustup component add rustfmt
RUN rustup target add ${RUST_TARGET}
RUN apt-get -y update
# Pin libclang/clang to 16 because older libclang releases can be too old
# for the clang-sys/bindgen APIs needed by hts-sys during release builds.
RUN apt-get -y install musl musl-tools musl-dev build-essential clang-16 libclang-16-dev
ENV CLANG_PATH=/usr/bin/clang-16
ENV LIBCLANG_PATH=/usr/lib/llvm-16/lib
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
