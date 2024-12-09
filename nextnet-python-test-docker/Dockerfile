FROM python:3.9

RUN apt-get update && \
    apt-get install -y git build-essential

RUN pip install matplotlib

COPY test.sh test.py /
RUN chmod 755 /test.sh

ENTRYPOINT ["/test.sh"]
