FROM python:3.9

RUN apt-get update && \
    apt-get install -y git build-essential

COPY id_rsa id_rsa.pub known_hosts /root/.ssh/
RUN chmod 700 /root/.ssh/
RUN chmod 600 /root/.ssh/known_hosts
RUN chmod 600 /root/.ssh/id_rsa
RUN chmod 644 /root/.ssh/id_rsa.pub

RUN pip install matplotlib

COPY test.sh test.py /
RUN chmod 755 /test.sh

ENTRYPOINT ["/test.sh"]
