#!/bin/bash


# (Bistable switch with positive feedback loop w/ minimal promoters)
# BFP ->        (test stability)
./ssatest_minfeed output/minfeed_1_0.txt 1 0
# mCITRINE ->       (test stability)
./ssatest_minfeed output/minfeed_2_0.txt 2 0
# BFP -> mCITRINE            (test switching)
./ssatest_minfeed output/minfeed_1_2.txt 1 2
# mCITRINE -> BFP            (test switching)
./ssatest_minfeed output/minfeed_2_1.txt 2 1


# (Bistable switch with positive feedback loop w/ constitutive promoters)
# BFP ->        (test stability)
./ssatest_cmvfeed output/cmvfeed_1_0.txt 1 0
# mCITRINE ->       (test stability)
./ssatest_cmvfeed output/cmvfeed_2_0.txt 2 0
# BFP -> mCITRINE            (test switching)
./ssatest_cmvfeed output/cmvfeed_1_2.txt 1 2
# mCITRINE -> BFP            (test switching)
./ssatest_cmvfeed output/cmvfeed_2_1.txt 2 1


# (Classic toggle switch with no positive feedback) 
# BFP ->        (test stability)
./ssatest_toggle output/toggle_1_0.txt 1 0
# mCITRINE ->       (test stability)
./ssatest_toggle output/toggle_2_0.txt 2 0
# BFP -> mCITRINE            (test switching)
./ssatest_toggle output/toggle_1_2.txt 1 2
# mCITRINE -> BFP            (test switching)
./ssatest_toggle output/toggle_2_1.txt 2 1


#  (Bistable switch with positive feedback loop w/ constitutive promoters, no positive feedback)
# BFP ->        (test stability)
./ssatest_noloop output/noloop_1_0.txt 1 0
# mCITRINE ->       (test stability)
./ssatest_noloop output/noloop_2_0.txt 2 0
# BFP -> mCITRINE            (test switching)
./ssatest_noloop output/noloop_1_2.txt 1 2
# mCITRINE -> BFP            (test switching)
./ssatest_noloop output/noloop_2_1.txt 2 1

# (Bistable switch with positive feedback loop w/ minimal promoters, no competition for binding sites)
# BFP ->        (test stability)
./ssatest_nocomp output/nocomp_1_0.txt 1 0
# mCITRINE ->       (test stability)
./ssatest_nocomp output/nocomp_2_0.txt 2 0
# BFP -> mCITRINE            (test switching)
./ssatest_nocomp output/nocomp_1_2.txt 1 2
# mCITRINE -> BFP            (test switching)
./ssatest_nocomp output/nocomp_2_1.txt 2 1

