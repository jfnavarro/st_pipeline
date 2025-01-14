FROM python:3.11-slim

# Set environment variables
ENV POETRY_VERSION=1.7.2 \
    PYTHONUNBUFFERED=1 \
    POETRY_NO_INTERACTION=1 \
    PATH="/root/.local/bin:$PATH"

# Install system dependencies and Poetry
RUN apt-get update \
    && apt-get install -y --no-install-recommends \
       build-essential \
       curl \
       libpq-dev \
       libffi-dev \
       libssl-dev \
       git \
       gcc \
    && curl -sSL https://install.python-poetry.org | python3 - \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

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
