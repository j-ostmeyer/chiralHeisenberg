#!/bin/bash

for length in 16 20 24 ; do
  echo $length

  touch mycpu.c
  sed "s/ksize=12,ksizet=12/ksize=${length},ksizet=${length}/g" -i cheisdwf.f
  sleep 1
  make
  sed "s/ksize=${length},ksizet=${length}/ksize=12,ksizet=12/g" -i cheisdwf.f

  if [ "$length" = "8" ]; then
    cores=4
  elif [ "$length" = "12" ]; then
    cores=4
  elif [ "$length" = "16" ]; then
    cores=6
  elif [ "$length" = "20" ]; then
    cores=19
  elif [ "$length" = "24" ]; then
    cores=19
  fi

  for Ls in 8 16 24 ; do
    #for b in 0.25 0.275 0.3 0.325 0.35 0.375 0.4 0.425 ; do
    for b in 0.45 0.475 0.5 ; do
      dir="L_${length}_Ls_${Ls}_b_${b}"
      #rm -rf $dir
      mkdir -p $dir
      rm -f $dir/*con*
      cp _con condp $dir
      cp remez* batch.sh $dir
      cp cheisdwf "${dir}/xcheis"
      cd $dir

      #cp condp con

      echo "0.0625   ${b} 1.0  0.00   0.    1    ${Ls}     5000" > midout
      echo "dt     beta am3    am    h0    imass iterl  iter2" >> midout

      sed "s/heis_40/heis_${length}_${Ls}_${b}/; s/-c 4/-c ${cores}/; s/NUM_THREADS=4/NUM_THREADS=${cores}/" -i batch.sh

      for job in {1..5}; do
	sbatch batch.sh
      done

      cd ..
    done
  done

done
