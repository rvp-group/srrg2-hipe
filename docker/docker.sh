#!/bin/bash
# This file is covered by the LICENSE file in the root of this project.
docker build -t hipe .
docker run --privileged \
       -ti --rm \
       --net=host \
       -v $1:/root/data/ \
       hipe
