FROM python:3.11-slim AS builder

ENV CYTHON_VERSION 0.29.34

# Install build dependencies
RUN apt-get update && \
    apt-get upgrade -y  && \
    apt-get install -yq --no-install-recommends build-essential && \
    pip install --no-cache-dir cython==${CYTHON_VERSION}

WORKDIR /tmp

# Build Python package wheels
WORKDIR /wheels
COPY requirements.txt .
RUN pip wheel --no-cache-dir -r requirements.txt && \
    rm requirements.txt

FROM python:3.11-slim

RUN apt-get update && \
    apt-get upgrade -y  && \
    apt-get install -yq --no-install-recommends \
        libmagic1 \
        procps

WORKDIR /tmp

# Install wheels
COPY --from=builder /wheels /wheels
RUN pip install --upgrade pip \
    pip install --no-warn-script-location /wheels/*.whl && \
    rm -rf /wheels && \
    rm -rf /root/.cache/pip/*

# Install pyquest package
COPY pyproject.toml setup.cfg ./
COPY src/ ./src
RUN pip install --no-cache-dir . && \
    rm -rf src/

WORKDIR /home

CMD ["/bin/bash"]
