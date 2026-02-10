ARG MOSDEPTH_VERSION=0.3.12

FROM ubuntu:22.04

# Install dependencies
RUN apt-get update && \
    apt-get install -y --no-install-recommends \
    procps \
    ca-certificates \
    curl \
    && rm -rf /var/lib/apt/lists/*

# Download and install mosdepth
RUN curl -L -o /usr/local/bin/mosdepth https://github.com/brentp/mosdepth/releases/download/v0.3.12/mosdepth_d4 \
    && chmod +x /usr/local/bin/mosdepth \
    && file /usr/local/bin/mosdepth
