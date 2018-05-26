#!/bin/csh
#  Dont forget to make this file executable!

cd source
make clean
make BUILDFLAG=OPT
cd ..
