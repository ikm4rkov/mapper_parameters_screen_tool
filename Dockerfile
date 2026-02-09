FROM python:3.9.6-slim

# Prevent interactive prompts
ENV DEBIAN_FRONTEND=noninteractive

# Install basic Linux tools needed by bash pipelines
RUN apt-get update && apt-get install -y \
    bash \
    coreutils \
    procps \
    && rm -rf /var/lib/apt/lists/*

# Set working directory
WORKDIR /app

# Copy requirements first (better caching)
COPY requirements.txt .

RUN pip install --no-cache-dir -r requirements.txt

# Copy repository
COPY . /app

# Default shell (important for bash-heavy repos)
SHELL ["/bin/bash", "-c"]

# No ENTRYPOINT on purpose — this is a toolbox container
CMD ["bash"]

