#/bin/bash
docker build --label nextnet-python-test .
docker run --mount type=bind,source=/tmp,target=/out nextnet-python-test
