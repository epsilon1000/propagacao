#!/bin/bash
gcc -lm prop.c
rm datos.dat
./a.out > datos.dat
