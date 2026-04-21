#!/bin/bash

FLUKA_BIN="{rel_fluka_bin}"

"$FLUKA_BIN/fff" neutron_scoring.f
"$FLUKA_BIN/lfluka" -o fluka_neutron neutron_scoring.o
"$FLUKA_BIN/rfluka" -N0 -M1 -e ./fluka_neutron fluka_input.inp
