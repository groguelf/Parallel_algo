#! /usr/bin/env bash

if [ "$#" -ne 2 ]; then
  echo ""
  echo "Usage: ./script.sh number_of_processes sub_matrix_size"
  echo ""
else
  echo ""
  make clean
  echo ""
  echo "-------------------------------------------------------------------------"
  echo ""
  make
  echo ""
  echo "-------------------------------------------------------------------------"
  echo ""
  mpirun -n $1 ./exercise_1 $2
  echo ""
fi
