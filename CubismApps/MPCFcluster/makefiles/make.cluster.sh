#!/bin/bash
# File       : make.cluster.sh
# Date       : Thu 01 Sep 2016 10:21:07 AM CEST
# Author     : Fabian Wermelinger
# Description: make cluster-layer
# Copyright 2016 ETH Zurich. All Rights Reserved.
BUILD="$(pwd -P)"
make -C $BUILD cleanall;
time make -C $BUILD \
    -j 2 CC=mpic++ \
    bs=16 \
    config=release \
    ap=double \
    hdf=1 \
    fftw=0 \
    nonuniform=0 \
    loa1r1=700.0 \
    loa2r2=0.7 \
    use_fpzip=1 \
    use_wavz=1 \
    qpxemu=0 \
    preclevel=2 \
    microfusion=2 \
    accurateweno=0 \
    hpm=0 \
    hdf-inc="${HDF5_ROOT}/include" \
    hdf-lib="${HDF5_ROOT}/lib" \
    extra='-g'
