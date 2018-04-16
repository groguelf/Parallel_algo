#! /usr/bin/env bash

if [ "$#" -ne 3 ]; then
  echo ""
  echo -e "Usage: \n\t For exercise_1 and exercise_2: \n\t\t./script.sh program number_of_processes sub_matrix_size"
  echo -e "\n\t For fox_2: \n\t\t./script.sh program number_of_processes matrix_dimension"
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
