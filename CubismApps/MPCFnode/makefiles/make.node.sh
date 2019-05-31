#!/bin/bash
# File       : make.node.sh
# Date       : Thu 01 Sep 2016 10:12:50 AM CEST
# Author     : Fabian Wermelinger
# Description: make node-layer
# Copyright 2016 ETH Zurich. All Rights Reserved.
BUILD="$(pwd -P)"
make -C $BUILD cleanall;
time make -C $BUILD \
    -j 2 CC=mpic++ \
    bs=32 \
    config=release \
    ap=float \
    hdf=1 \
    fftw=0 \
    nonuniform=0 \
    qpxemu=0 \
    preclevel=2 \
    microfusion=2 \
    accurateweno=0 \
    hdf-inc="${HDF5_ROOT}/include" \
    hdf-lib="${HDF5_ROOT}/lib" \
    extra='-g'
    # extra='-g -D_WENO3_ -D_BCLABCLOUD1DCHARNREF_ -D_CONVERTCLIP_'
