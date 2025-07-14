# Use an official Python runtime as the base image (Python 3.10 is used here)
FROM python:3.10-slim

# Set environment variables to prevent Python from writing pyc files and to buffer outputs
ENV PYTHONDONTWRITEBYTECODE=1
ENV PYTHONUNBUFFERED=1

# Set the working directory inside the Docker container
WORKDIR /usr/src/app

# Install system dependencies required for building Python packages
RUN apt-get update && apt-get install -y \
    build-essential \
    && apt-get clean && rm -rf /var/lib/apt/lists/*

# Upgrade pip to the latest version
RUN pip install --upgrade pip

# Install the latest stable version of liana with common dependencies
RUN pip install liana[common]

# (Optional) If you need all dependencies (including development dependencies), uncomment the following line:
# RUN pip install liana[full]

# Install pydeseq2
RUN pip install pydeseq2

RUN pip install corneto cvxpy-base scipy