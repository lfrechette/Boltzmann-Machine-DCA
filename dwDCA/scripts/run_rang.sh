#!/bin/bash

mydir=$1
msa=$2
conf=$3

bmdca -i $msa -d $mydir -r -c $conf -t 0.8
