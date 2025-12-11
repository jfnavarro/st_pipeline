FROM python:3.11.14-slim-bookworm AS builder

# Environment for deterministic, non-interactive builds
ENV PYTHONUNBUFFERED=1 \
    VENV_PATH=/opt/venv

# System deps for building Python wheels and bio tools
RUN apt-get update && apt-get install -y --no-install-recommends \
    build-essential \
    curl \
    gcc \
    git \
    wget \
    unzip \
    zlib1g-dev \
    libbz2-dev \
    liblzma-dev \
    libcurl4-gnutls-dev \
    libncurses5-dev \
    libffi-dev \
    libssl-dev \
    libpq-dev \
    && rm -rf /var/lib/apt/lists/*

# Create and use a dedicated venv
RUN python -m venv "${VENV_PATH}"
ENV PATH="${VENV_PATH}/bin:${PATH}"

# Working directory
WORKDIR /app/stpipeline

# Copy project
COPY pyproject.toml README.md ./
COPY ./stpipeline/ ./stpipeline

# Install project
RUN pip install .

# Build STAR
WORKDIR /tmp
RUN wget -O STAR-2.7.10b.zip "https://github.com/alexdobin/STAR/archive/refs/tags/2.7.10b.zip" \
    && unzip STAR-2.7.10b.zip \
    && make -C STAR-2.7.10b/source STAR \
    && install -m 0755 STAR-2.7.10b/source/STAR /usr/local/bin/STAR \
    && rm -rf STAR-2.7.10b STAR-2.7.10b.zip

# Build Samtools
RUN wget -O samtools-1.17.tar.bz2 "https://github.com/samtools/samtools/releases/download/1.17/samtools-1.17.tar.bz2" \
    && tar -xjf samtools-1.17.tar.bz2 \
    && cd samtools-1.17 \
    && ./configure \
    && make -j"$(nproc)" \
    && make install \
    && cd /tmp \
    && rm -rf samtools-1.17 samtools-1.17.tar.bz2

# Strip binaries to reduce size
RUN strip /usr/local/bin/STAR || true
RUN find /usr/local/bin -maxdepth 1 -type f -name 'samtools' -exec strip {} \; || true

FROM python:3.11.14-slim-bookworm

# Minimal runtime deps
RUN apt-get update && apt-get install -y --no-install-recommends \
    libbz2-dev \
    liblzma-dev \
    libcurl4-gnutls-dev \
    libncurses5-dev \
    && rm -rf /var/lib/apt/lists/*

# Copy venv and installed binaries from builder
ENV VENV_PATH=/opt/venv
COPY --from=builder ${VENV_PATH} ${VENV_PATH}
COPY --from=builder /usr/local/bin/STAR /usr/local/bin/STAR
COPY --from=builder /usr/local/bin/samtools /usr/local/bin/samtools

# Add scripts to PATH
ENV PATH="${VENV_PATH}/bin:${PATH}"
ENV PYTHONUNBUFFERED=1

# App workdir
WORKDIR /app/stpipeline

# Create a non-root user and ensure permissions
RUN adduser --system --group app \
    && chown -R app:app /app
USER app

CMD ["st_pipeline_run", "--help"]
