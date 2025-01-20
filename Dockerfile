FROM python:3.11-slim

# Set environment variables
ENV POETRY_VERSION=1.7.2 \
    PYTHONUNBUFFERED=1 \
    POETRY_NO_INTERACTION=1 \
    PATH="/root/.local/bin:$PATH"

# Install system dependencies, Poetry, STAR, and Samtools
RUN apt-get update \
    && apt-get install -y --no-install-recommends \
       build-essential \
       curl \
       libpq-dev \
       libffi-dev \
       libssl-dev \
       git \
       gcc \
       wget \
       unzip \
       zlib1g-dev \
       libbz2-dev \
       liblzma-dev \
       libcurl4-gnutls-dev \
    && wget https://github.com/alexdobin/STAR/archive/refs/tags/2.7.10b.zip \
    && unzip 2.7.10b.zip \
    && cd STAR-2.7.10b/source \
    && make STAR \
    && mv STAR /usr/local/bin/ \
    && cd /app \
    && wget https://github.com/samtools/samtools/releases/download/1.17/samtools-1.17.tar.bz2 \
    && tar -xjf samtools-1.17.tar.bz2 \
    && cd samtools-1.17 \
    && ./configure \
    && make \
    && make install \
    && cd /app \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/* 2.7.10b.zip STAR-2.7.10b samtools-1.17 samtools-1.17.tar.bz2

# Set working directory
WORKDIR /app

# Copy project files
COPY pyproject.toml poetry.lock README.md /app/

# Install dependencies using Poetry
RUN poetry install --no-root --only main

# Copy the entire project
COPY . /app

# Ensure scripts are executable
RUN chmod +x /app/stpipeline/scripts/*.py

# Set entrypoint for the container
ENTRYPOINT ["poetry", "run"]
