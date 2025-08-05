FROM python:3.13

WORKDIR /app

COPY . .

RUN pip install --no-cache-dir --upgrade pip && pip install --no-cache-dir scanpy

RUN pip install -e .

