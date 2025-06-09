#/bin/bash
set -e

docker build -t nextnet-python-test .

rm -f /tmp/test.pdf
docker run --mount type=bind,source=/tmp,target=/out nextnet-python-test
