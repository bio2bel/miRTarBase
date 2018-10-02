# 1. Build with: docker build -t bio2bel/mirtarbasebase .
# 2. Run with docker run -p 5000:5000 bio2bel/mirtarbase (-p maps the port)
FROM python:3.7.0
MAINTAINER Charles Tapley Hoyt "cthoyt@gmail.com"

COPY . /app
WORKDIR /app

RUN pip install .[web]

RUN bio2bel_mirtarbase populate

EXPOSE 5000
ENTRYPOINT bio2bel_mirtarbase web --port 5000 --host 0.0.0.0
