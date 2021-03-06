#!/usr/bin/env bash
# python3.9 -m venv env
# source env/bin/activate
# pip install numpy casadi
# pip freeze > requirements.txt
docker login -u nathanlazarus
docker build -t nathanlazarus/oilgame:v3 .
docker push nathanlazarus/oilgame:v3
