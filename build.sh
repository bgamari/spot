#!/bin/bash

ghc Spot.hs -threaded -O2 -fno-liberate-case -funfolding-use-threshold1000 -funfolding-keeness-factor1000 -optlo-O3
time ./Spot +RTS -N4 -s >| h
python plot.py < h
