# set base image (host OS)
FROM python:3.9-slim-buster

# set the working directory in the container
WORKDIR /code

# COPY OilGame.py .
# COPY helper_functions.py .

# COPY requirements.txt .
# RUN pip install --no-cache-dir -r requirements.txt
RUN pip install --no-cache-dir numpy
RUN pip install --no-cache-dir casadi
RUN pip install --no-cache-dir netCDF4