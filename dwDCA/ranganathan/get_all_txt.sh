#!/bin/bash

END=100

for i in {1..$END}
do
    arma2ascii -s stat_MC_1p_$i.bin
    arma2ascii -s stat_MC_2p_$i.bin
done
