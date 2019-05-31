#!/bin/bash
set -x #echo on

##################
h5file=/home/chatzidp/gitlab/fabdata/medium/data_010000-p-med.h5
./runone.sh $h5file

##################
h5file=/home/chatzidp/gitlab/fabdata/medium/data_010000-rho-med.h5
./runone.sh $h5file

##################
h5file=/home/chatzidp/gitlab/fabdata/medium/data_010000-E-med.h5
./runone.sh $h5file

##################
h5file=/home/chatzidp/gitlab/fabdata/medium/data_010000-a2-med.h5
./runone.sh $h5file

##################
#h5file=/home/chatzidp/gitlab/fabdata/data_010000-divU-med.h5
#./runone.sh $h5file

##################
#h5file=/home/chatzidp/gitlab/fabdata/data_010000-Omegax-med.h5
#./runone.sh $h5file

##################
#h5file=/home/chatzidp/gitlab/fabdata/data_010000-Ux-med.h5
#./runone.sh $h5file

