FROM python:3.11-slim

# Set environment variables
ENV POETRY_VERSION=1.7.2 \
    PYTHONUNBUFFERED=1 \
    POETRY_NO_INTERACTION=1 \
    PATH="/root/.local/bin:$PATH"

# Install system dependencies, Poetry, and STAR
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
    && wget https://github.com/alexdobin/STAR/archive/refs/tags/2.7.10b.zip \
    && unzip 2.7.10b.zip \
    && cd STAR-2.7.10b/source \
    && make STAR \
    && mv STAR /usr/local/bin/ \
    && cd /app \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/* 2.7.10b.zip STAR-2.7.10b

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
