#!/bin/bash

gfortran main.f90 > /dev/null
./a.out
python3 plot.py & disown
