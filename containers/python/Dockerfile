FROM python:3.11-slim

# Copy the requirements file
COPY requirements_full.txt requirements.txt

# Install the project dependencies
RUN pip install -r requirements.txt

# Run your application
CMD [ "python"]