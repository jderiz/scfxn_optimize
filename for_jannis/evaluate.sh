#!/bin/sh

mkdir opt_beta16

#parallel -S 20/dig4,20/dig6,20/dig7,20/dig8,20/dig11,20/dig16,20/dig21,20/dig23 --workdir . :::: alljobs ::: beta16 ::: beta_nov16.wts
parallel -j 6 --workdir . :::: alljobs ::: beta16 ::: beta_nov16.wts

./run_cleanup.sh opt_beta16
