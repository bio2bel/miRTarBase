FROM python:3.6.5
MAINTAINER Charles Tapley Hoyt "cthoyt@gmail.com"

RUN pip3 install --upgrade pip
RUN pip3 install gunicorn mysqlclient

ADD requirements.txt /
RUN pip3 install -r requirements.txt

COPY . /app
WORKDIR /app

RUN pip3 install .
