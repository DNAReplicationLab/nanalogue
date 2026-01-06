FROM rust:1.92.0-bullseye as build-env
RUN rustup component add rustfmt
RUN rustup target add x86_64-unknown-linux-musl
RUN apt-get -y update
RUN apt-get -y install musl musl-tools
WORKDIR /app
COPY . /app
RUN RUSTFLAGS='-C link-arg=-s' cargo build --release --target x86_64-unknown-linux-musl
RUN cargo test --release --target x86_64-unknown-linux-musl

FROM gcr.io/distroless/static-debian12
COPY --from=build-env /app/target/x86_64-unknown-linux-musl/release/nanalogue /
COPY --from=build-env /app/target/x86_64-unknown-linux-musl/release/nanalogue_sim_bam /
ENV PATH="$PATH:/"
CMD ["nanalogue"]
