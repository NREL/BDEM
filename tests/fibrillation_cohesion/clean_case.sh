#!/bin/bash

rm -r plt*
rm -r rst*
rm -r tri*
rm -r ebpl*
rm incline_*
rm block_*
rm Backtrace*
rm particle_input.dat
python3 generate_particles.py
