FROM python:3.11-slim

# Copy the requirements file
COPY requirements_full.txt requirements.txt

# Install the project dependencies
RUN pip install -r requirements.txt
# Adding ps to the container (required by nextflow)
RUN apt-get update && \
    apt-get install -y procps=2:4.0.2-3 && \
    rm -rf /var/lib/apt/lists/*

# Run your application
CMD [ "python"]
