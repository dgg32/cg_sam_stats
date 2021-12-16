# syntax=docker/dockerfile:1
FROM continuumio/miniconda

WORKDIR /app

COPY requirements.txt requirements.txt

RUN pip install -r requirements.txt

COPY . .


ENV PATH="/app:${PATH}"