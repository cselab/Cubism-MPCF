#!/bin/bash
# File       : make.core.sh
# Date       : Fri 02 Sep 2016 09:37:06 AM CEST
# Author     : Fabian Wermelinger
# Description: make core-layer
# Copyright 2016 ETH Zurich. All Rights Reserved.
BUILD="$(pwd -P)"
make -C $BUILD clean;
time make -C $BUILD \
    -j 2 CC=gcc \
    bs=32 \
    config=release \
    ap=float \
    hdf=0 \
    fftw=0 \
    qpxemu=1 \
    preclevel=0 \
    microfusion=2 \
    accurateweno=0 \
    extra='-g'
