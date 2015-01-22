#!/bin/sh

cat ssatest.c | sed "s/_SWITCH_/CMVFEED/" > ssatest_cmvfeed.c
gcc ssatest_cmvfeed.c -o ssatest_cmvfeed -lm -I models -I weights
./ssatest_cmvfeed
gcc ssatest_cmvfeed.c -o ssatest_cmvfeed -lm -I models -I weights

cat ssatest.c | sed "s/_SWITCH_/MINFEED/" > ssatest_minfeed.c
gcc ssatest_minfeed.c -o ssatest_minfeed -lm -I models -I weights
./ssatest_minfeed
gcc ssatest_minfeed.c -o ssatest_minfeed -lm -I models -I weights

cat ssatest.c | sed "s/_SWITCH_/TOGGLE/" > ssatest_toggle.c
gcc ssatest_toggle.c -o ssatest_toggle -lm -I models -I weights
./ssatest_toggle
gcc ssatest_toggle.c -o ssatest_toggle -lm -I models -I weights

cat ssatest.c | sed "s/_SWITCH_/NOLOOP/" > ssatest_noloop.c
gcc ssatest_noloop.c -o ssatest_noloop -lm -I models -I weights
./ssatest_noloop
gcc ssatest_noloop.c -o ssatest_noloop -lm -I models -I weights

cat ssatest.c | sed "s/_SWITCH_/NOCOMP/" > ssatest_nocomp.c
gcc ssatest_nocomp.c -o ssatest_nocomp -lm -I models -I weights
./ssatest_nocomp
gcc ssatest_nocomp.c -o ssatest_nocomp -lm -I models -I weights
