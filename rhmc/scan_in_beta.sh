#!/bin/bash

Ls=8
for b in 0.25 0.275 0.3 0.325 0.35 0.375 0.4 0.425 ; do
  dir="L_16_Ls_${Ls}_b_${b}"
  rm -rf $dir
  mkdir -p $dir
  #cp _con condp $dir
  cp remez* batch.sh $dir
  cp cheisdwf "${dir}/xcheis"
  cd $dir

  #cp condp con

  echo "0.125   ${b} 1.0  0.00   0.    1    ${Ls}     2500" > midout
  echo "dt     beta am3    am    h0    imass iterl  iter2" >> midout

  sed "s/heis_40/heis_${Ls}_${b}/" -i batch.sh

  sbatch batch.sh

  cd ..
done
