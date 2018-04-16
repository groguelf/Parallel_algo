#! /usr/bin/env bash

if [ "$#" -lt 2 ] || [ "$#" -gt 3 ]; then
  echo ""
  echo -e "Usage: \n\t for exercise_1 and exercise_2: \n\t\t./script.sh program number_of_processes sub_matrix_size"
  echo -e "\n\t for fox_2: \n\t\t./script.sh program number_of_processes"
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
