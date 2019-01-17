#/bin/bash
for exe in dkd kdk kdkdk2 kdkdk4
do
  for tick in 512 1024 2048 4096 8192 16384
  do
    ./$exe pl1k.txt $tick | grep err >> eps256i/$exe.dat
  done
done
