#! /bin/bash

GSL_LOCATION="$HOME/gsl/lib"
LD_LIBRARY_PATH="$LD_LIBRARY_PATH:$GSL_LOCATION"
make parallel > /dev/null &2>1
if [ $? == 0 ]
then
    mpirun -np $1 -x LD_LIBRARY_PATH --hostfile my_hosts parallel
else
    echo "Error making"
fi
