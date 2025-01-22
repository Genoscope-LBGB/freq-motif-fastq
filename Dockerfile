FROM rust:1.73

# copy the git repo to the container
COPY . /gitlab/repo

# change directory for the remaining run commands
WORKDIR /gitlab/repo

# install the project module
RUN cd freq-motif-fastq && cargo build --release
RUN cargo install --path freq-motif-fastq --root .

