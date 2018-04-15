#! /usr/bin/env bash

if [ "$#" -ne 3 ]; then
  echo ""
  echo "Usage: ./script.sh program number_of_processes sub_matrix_size"
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
  mpirun -n $2 ./$1 $3
  echo ""
fi
