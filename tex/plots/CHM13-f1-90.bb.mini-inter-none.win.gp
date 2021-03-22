set t po eps co so enh "Helvetica,18"
set out "chr-plot.eps"
set size 2,1.52
set multiplot layout 23,1
set lmargin screen 0.095
set border 0; unset xtics; unset ytics; set bmargin 0; set tmargin 0.02; set rmargin 0.02
set style line 1 lc rgb "#377eb8" lw 1
set style line 2 lc rgb "#e41a1c" lw 1
set style line 3 lc rgb "#4daf4a" lw 1
set yran [0:164]

set style fill solid 0.8

set origin 0,1.4447826086956521
set xran [0:248.387497]
set size 2,0.06521739130434782
set style rect fc lt -1 fs solid 0.15 noborder
unset obj; unset label
set obj rect from 116.796216, graph 0 to 147.241828, graph 1
set label "chr1" at screen 0.01, graph 0.5
set key at screen 1.95,1.32
plot \
     "<awk '$1==\"chr1\"' CHM13-f1-90.bb.mini-inter-none.win" u ($2*1e-6):(0):($3) t "VNTR" w filledcu ls 1, \
     "<awk '$1==\"chr1\"' CHM13-f1-90.bb.mini-inter-none.win" u ($2*1e-6):($3):($3+$4) not w filledcu ls 2, \
     "<awk '$1==\"chr1\"' CHM13-f1-90.bb.mini-inter-none.win" u ($2*1e-6):($3+$4):($3+$4+$5) not w filledcu ls 3
set origin 0,1.3795652173913044
set xran [0:242.696747]
set size 1.95417845045558,0.06521739130434782
set style rect fc lt -1 fs solid 0.15 noborder
unset obj; unset label
set obj rect from 85.991672, graph 0 to 99.67301599999999, graph 1
set label "chr2" at screen 0.01, graph 0.5
set key at screen 1.95,1.28
plot \
     "<awk '$1==\"chr2\"' CHM13-f1-90.bb.mini-inter-none.win" u ($2*1e-6):(0):($3) not w filledcu ls 1, \
     "<awk '$1==\"chr2\"' CHM13-f1-90.bb.mini-inter-none.win" u ($2*1e-6):($3):($3+$4) t "Intersperse" w filledcu ls 2, \
     "<awk '$1==\"chr2\"' CHM13-f1-90.bb.mini-inter-none.win" u ($2*1e-6):($3+$4):($3+$4+$5) not w filledcu ls 3
set origin 0,1.3143478260869565
set xran [0:201.106605]
set size 1.619297327191956,0.06521739130434782
set style rect fc lt -1 fs solid 0.15 noborder
unset obj; unset label
set obj rect from 85.80519199999999, graph 0 to 101.415517, graph 1
set label "chr3" at screen 0.01, graph 0.5
set key at screen 1.95,1.24
plot \
     "<awk '$1==\"chr3\"' CHM13-f1-90.bb.mini-inter-none.win" u ($2*1e-6):(0):($3) not w filledcu ls 1, \
     "<awk '$1==\"chr3\"' CHM13-f1-90.bb.mini-inter-none.win" u ($2*1e-6):($3):($3+$4) not w filledcu ls 2, \
     "<awk '$1==\"chr3\"' CHM13-f1-90.bb.mini-inter-none.win" u ($2*1e-6):($3+$4):($3+$4+$5) t "Partial/non-repeat" w filledcu ls 3
set origin 0,1.2491304347826087
set xran [0:193.57542999999998]
set size 1.5586567950318369,0.06521739130434782
set style rect fc lt -1 fs solid 0.15 noborder
unset obj; unset label
set obj rect from 44.705247, graph 0 to 59.870604, graph 1
set label "chr4" at screen 0.01, graph 0.5
plot \
     "<awk '$1==\"chr4\"' CHM13-f1-90.bb.mini-inter-none.win" u ($2*1e-6):(0):($3) not w filledcu ls 1, \
     "<awk '$1==\"chr4\"' CHM13-f1-90.bb.mini-inter-none.win" u ($2*1e-6):($3):($3+$4) not w filledcu ls 2, \
     "<awk '$1==\"chr4\"' CHM13-f1-90.bb.mini-inter-none.win" u ($2*1e-6):($3+$4):($3+$4+$5) not w filledcu ls 3
set origin 0,1.1839130434782608
set xran [0:182.045437]
set size 1.4658180399474776,0.06521739130434782
set style rect fc lt -1 fs solid 0.15 noborder
unset obj; unset label
set obj rect from 42.077197, graph 0 to 54.596619, graph 1
set label "chr5" at screen 0.01, graph 0.5
plot \
     "<awk '$1==\"chr5\"' CHM13-f1-90.bb.mini-inter-none.win" u ($2*1e-6):(0):($3) not w filledcu ls 1, \
     "<awk '$1==\"chr5\"' CHM13-f1-90.bb.mini-inter-none.win" u ($2*1e-6):($3):($3+$4) not w filledcu ls 2, \
     "<awk '$1==\"chr5\"' CHM13-f1-90.bb.mini-inter-none.win" u ($2*1e-6):($3+$4):($3+$4+$5) not w filledcu ls 3
set origin 0,1.118695652173913
set xran [0:172.12687]
set size 1.3859543823979192,0.06521739130434782
set style rect fc lt -1 fs solid 0.15 noborder
unset obj; unset label
set obj rect from 53.286919999999995, graph 0 to 66.058622, graph 1
set label "chr6" at screen 0.01, graph 0.5
plot \
     "<awk '$1==\"chr6\"' CHM13-f1-90.bb.mini-inter-none.win" u ($2*1e-6):(0):($3) not w filledcu ls 1, \
     "<awk '$1==\"chr6\"' CHM13-f1-90.bb.mini-inter-none.win" u ($2*1e-6):($3):($3+$4) not w filledcu ls 2, \
     "<awk '$1==\"chr6\"' CHM13-f1-90.bb.mini-inter-none.win" u ($2*1e-6):($3+$4):($3+$4+$5) not w filledcu ls 3
set origin 0,1.0534782608695652
set xran [0:160.567423]
set size 1.2928784656177763,0.06521739130434782
set style rect fc lt -1 fs solid 0.15 noborder
unset obj; unset label
set obj rect from 55.414367999999996, graph 0 to 68.714496, graph 1
set label "chr7" at screen 0.01, graph 0.5
plot \
     "<awk '$1==\"chr7\"' CHM13-f1-90.bb.mini-inter-none.win" u ($2*1e-6):(0):($3) not w filledcu ls 1, \
     "<awk '$1==\"chr7\"' CHM13-f1-90.bb.mini-inter-none.win" u ($2*1e-6):($3):($3+$4) not w filledcu ls 2, \
     "<awk '$1==\"chr7\"' CHM13-f1-90.bb.mini-inter-none.win" u ($2*1e-6):($3+$4):($3+$4+$5) not w filledcu ls 3
set origin 0,0.9882608695652174
set xran [0:146.259322]
set size 1.1776705652780906,0.06521739130434782
set style rect fc lt -1 fs solid 0.15 noborder
unset obj; unset label
set obj rect from 39.243541, graph 0 to 51.325075999999996, graph 1
set label "chr8" at screen 0.01, graph 0.5
plot \
     "<awk '$1==\"chr8\"' CHM13-f1-90.bb.mini-inter-none.win" u ($2*1e-6):(0):($3) not w filledcu ls 1, \
     "<awk '$1==\"chr8\"' CHM13-f1-90.bb.mini-inter-none.win" u ($2*1e-6):($3):($3+$4) not w filledcu ls 2, \
     "<awk '$1==\"chr8\"' CHM13-f1-90.bb.mini-inter-none.win" u ($2*1e-6):($3+$4):($3+$4+$5) not w filledcu ls 3
set origin 0,0.9230434782608696
set xran [0:150.61727399999998]
set size 1.2127605118545883,0.06521739130434782
set style rect fc lt -1 fs solid 0.15 noborder
unset obj; unset label
set obj rect from 39.952788999999996, graph 0 to 81.69403299999999, graph 1
set label "chr9" at screen 0.01, graph 0.5
plot \
     "<awk '$1==\"chr9\"' CHM13-f1-90.bb.mini-inter-none.win" u ($2*1e-6):(0):($3) not w filledcu ls 1, \
     "<awk '$1==\"chr9\"' CHM13-f1-90.bb.mini-inter-none.win" u ($2*1e-6):($3):($3+$4) not w filledcu ls 2, \
     "<awk '$1==\"chr9\"' CHM13-f1-90.bb.mini-inter-none.win" u ($2*1e-6):($3+$4):($3+$4+$5) not w filledcu ls 3
set origin 0,0.8578260869565217
set xran [0:134.758122]
set size 1.0850636495604287,0.06521739130434782
set style rect fc lt -1 fs solid 0.15 noborder
unset obj; unset label
set obj rect from 34.633784, graph 0 to 46.66458, graph 1
set label "chr10" at screen 0.01, graph 0.5
plot \
     "<awk '$1==\"chr10\"' CHM13-f1-90.bb.mini-inter-none.win" u ($2*1e-6):(0):($3) not w filledcu ls 1, \
     "<awk '$1==\"chr10\"' CHM13-f1-90.bb.mini-inter-none.win" u ($2*1e-6):($3):($3+$4) not w filledcu ls 2, \
     "<awk '$1==\"chr10\"' CHM13-f1-90.bb.mini-inter-none.win" u ($2*1e-6):($3+$4):($3+$4+$5) not w filledcu ls 3
set origin 0,0.792608695652174
set xran [0:135.127772]
set size 1.0880400473619651,0.06521739130434782
set style rect fc lt -1 fs solid 0.15 noborder
unset obj; unset label
set obj rect from 46.061948, graph 0 to 59.413484999999994, graph 1
set label "chr11" at screen 0.01, graph 0.5
plot \
     "<awk '$1==\"chr11\"' CHM13-f1-90.bb.mini-inter-none.win" u ($2*1e-6):(0):($3) not w filledcu ls 1, \
     "<awk '$1==\"chr11\"' CHM13-f1-90.bb.mini-inter-none.win" u ($2*1e-6):($3):($3+$4) not w filledcu ls 2, \
     "<awk '$1==\"chr11\"' CHM13-f1-90.bb.mini-inter-none.win" u ($2*1e-6):($3+$4):($3+$4+$5) not w filledcu ls 3
set origin 0,0.7273913043478262
set xran [0:133.324781]
set size 1.0735224808839714,0.06521739130434782
set style rect fc lt -1 fs solid 0.15 noborder
unset obj; unset label
set obj rect from 29.62049, graph 0 to 42.202481999999996, graph 1
set label "chr12" at screen 0.01, graph 0.5
plot \
     "<awk '$1==\"chr12\"' CHM13-f1-90.bb.mini-inter-none.win" u ($2*1e-6):(0):($3) not w filledcu ls 1, \
     "<awk '$1==\"chr12\"' CHM13-f1-90.bb.mini-inter-none.win" u ($2*1e-6):($3):($3+$4) not w filledcu ls 2, \
     "<awk '$1==\"chr12\"' CHM13-f1-90.bb.mini-inter-none.win" u ($2*1e-6):($3+$4):($3+$4+$5) not w filledcu ls 3
set origin 0,0.6621739130434783
set xran [0:114.240146]
set size 0.9198542388790205,0.06521739130434782
set style rect fc lt -1 fs solid 0.15 noborder
unset obj; unset label
set obj rect from 0, graph 0 to 23.171058, graph 1
set label "chr13" at screen 0.01, graph 0.5
plot \
     "<awk '$1==\"chr13\"' CHM13-f1-90.bb.mini-inter-none.win" u ($2*1e-6):(0):($3) not w filledcu ls 1, \
     "<awk '$1==\"chr13\"' CHM13-f1-90.bb.mini-inter-none.win" u ($2*1e-6):($3):($3+$4) not w filledcu ls 2, \
     "<awk '$1==\"chr13\"' CHM13-f1-90.bb.mini-inter-none.win" u ($2*1e-6):($3+$4):($3+$4+$5) not w filledcu ls 3
set origin 0,0.5969565217391305
set xran [0:101.219177]
set size 0.8150102418399908,0.06521739130434782
set style rect fc lt -1 fs solid 0.15 noborder
unset obj; unset label
set obj rect from 0, graph 0 to 17.765925, graph 1
set label "chr14" at screen 0.01, graph 0.5
plot \
     "<awk '$1==\"chr14\"' CHM13-f1-90.bb.mini-inter-none.win" u ($2*1e-6):(0):($3) not w filledcu ls 1, \
     "<awk '$1==\"chr14\"' CHM13-f1-90.bb.mini-inter-none.win" u ($2*1e-6):($3):($3+$4) not w filledcu ls 2, \
     "<awk '$1==\"chr14\"' CHM13-f1-90.bb.mini-inter-none.win" u ($2*1e-6):($3+$4):($3+$4+$5) not w filledcu ls 3
set origin 0,0.5317391304347826
set xran [0:100.338308]
set size 0.8079175418398777,0.06521739130434782
set style rect fc lt -1 fs solid 0.15 noborder
unset obj; unset label
set obj rect from 0, graph 0 to 23.279251, graph 1
set label "chr15" at screen 0.01, graph 0.5
plot \
     "<awk '$1==\"chr15\"' CHM13-f1-90.bb.mini-inter-none.win" u ($2*1e-6):(0):($3) not w filledcu ls 1, \
     "<awk '$1==\"chr15\"' CHM13-f1-90.bb.mini-inter-none.win" u ($2*1e-6):($3):($3+$4) not w filledcu ls 2, \
     "<awk '$1==\"chr15\"' CHM13-f1-90.bb.mini-inter-none.win" u ($2*1e-6):($3+$4):($3+$4+$5) not w filledcu ls 3
set origin 0,0.4665217391304348
set xran [0:96.33049299999999]
set size 0.7756468756557421,0.06521739130434782
set style rect fc lt -1 fs solid 0.15 noborder
unset obj; unset label
set obj rect from 30.848291, graph 0 to 57.219476, graph 1
set label "chr16" at screen 0.01, graph 0.5
plot \
     "<awk '$1==\"chr16\"' CHM13-f1-90.bb.mini-inter-none.win" u ($2*1e-6):(0):($3) not w filledcu ls 1, \
     "<awk '$1==\"chr16\"' CHM13-f1-90.bb.mini-inter-none.win" u ($2*1e-6):($3):($3+$4) not w filledcu ls 2, \
     "<awk '$1==\"chr16\"' CHM13-f1-90.bb.mini-inter-none.win" u ($2*1e-6):($3+$4):($3+$4+$5) not w filledcu ls 3
set origin 0,0.4013043478260869
set xran [0:84.277185]
set size 0.6785944221661044,0.06521739130434782
set style rect fc lt -1 fs solid 0.15 noborder
unset obj; unset label
set obj rect from 18.892709999999997, graph 0 to 32.48723, graph 1
set label "chr17" at screen 0.01, graph 0.5
plot \
     "<awk '$1==\"chr17\"' CHM13-f1-90.bb.mini-inter-none.win" u ($2*1e-6):(0):($3) not w filledcu ls 1, \
     "<awk '$1==\"chr17\"' CHM13-f1-90.bb.mini-inter-none.win" u ($2*1e-6):($3):($3+$4) not w filledcu ls 2, \
     "<awk '$1==\"chr17\"' CHM13-f1-90.bb.mini-inter-none.win" u ($2*1e-6):($3+$4):($3+$4+$5) not w filledcu ls 3
set origin 0,0.33608695652173926
set xran [0:80.542536]
set size 0.6485232708794517,0.06521739130434782
set style rect fc lt -1 fs solid 0.15 noborder
unset obj; unset label
set obj rect from 10.965698, graph 0 to 25.93355, graph 1
set label "chr18" at screen 0.01, graph 0.5
plot \
     "<awk '$1==\"chr18\"' CHM13-f1-90.bb.mini-inter-none.win" u ($2*1e-6):(0):($3) not w filledcu ls 1, \
     "<awk '$1==\"chr18\"' CHM13-f1-90.bb.mini-inter-none.win" u ($2*1e-6):($3):($3+$4) not w filledcu ls 2, \
     "<awk '$1==\"chr18\"' CHM13-f1-90.bb.mini-inter-none.win" u ($2*1e-6):($3+$4):($3+$4+$5) not w filledcu ls 3
set origin 0,0.27086956521739136
set xran [0:61.707359]
set size 0.4968636484951576,0.06521739130434782
set style rect fc lt -1 fs solid 0.15 noborder
unset obj; unset label
set obj rect from 19.655572, graph 0 to 34.768167999999996, graph 1
set label "chr19" at screen 0.01, graph 0.5
plot \
     "<awk '$1==\"chr19\"' CHM13-f1-90.bb.mini-inter-none.win" u ($2*1e-6):(0):($3) not w filledcu ls 1, \
     "<awk '$1==\"chr19\"' CHM13-f1-90.bb.mini-inter-none.win" u ($2*1e-6):($3):($3+$4) not w filledcu ls 2, \
     "<awk '$1==\"chr19\"' CHM13-f1-90.bb.mini-inter-none.win" u ($2*1e-6):($3+$4):($3+$4+$5) not w filledcu ls 3
set origin 0,0.20565217391304347
set xran [0:66.210247]
set size 0.5331206103341023,0.06521739130434782
set style rect fc lt -1 fs solid 0.15 noborder
unset obj; unset label
set obj rect from 21.383653, graph 0 to 37.969530999999996, graph 1
set label "chr20" at screen 0.01, graph 0.5
plot \
     "<awk '$1==\"chr20\"' CHM13-f1-90.bb.mini-inter-none.win" u ($2*1e-6):(0):($3) not w filledcu ls 1, \
     "<awk '$1==\"chr20\"' CHM13-f1-90.bb.mini-inter-none.win" u ($2*1e-6):($3):($3+$4) not w filledcu ls 2, \
     "<awk '$1==\"chr20\"' CHM13-f1-90.bb.mini-inter-none.win" u ($2*1e-6):($3+$4):($3+$4+$5) not w filledcu ls 3
set origin 0,0.1404347826086958
set xran [0:45.827690999999994]
set size 0.3690015927009402,0.06521739130434782
set style rect fc lt -1 fs solid 0.15 noborder
unset obj; unset label
set obj rect from 0, graph 0 to 17.078862, graph 1
set label "chr21" at screen 0.01, graph 0.5
plot \
     "<awk '$1==\"chr21\"' CHM13-f1-90.bb.mini-inter-none.win" u ($2*1e-6):(0):($3) not w filledcu ls 1, \
     "<awk '$1==\"chr21\"' CHM13-f1-90.bb.mini-inter-none.win" u ($2*1e-6):($3):($3+$4) not w filledcu ls 2, \
     "<awk '$1==\"chr21\"' CHM13-f1-90.bb.mini-inter-none.win" u ($2*1e-6):($3+$4):($3+$4+$5) not w filledcu ls 3
set origin 0,0.07521739130434789
set xran [0:51.353905999999995]
set size 0.41349831710732204,0.06521739130434782
set style rect fc lt -1 fs solid 0.15 noborder
unset obj; unset label
set obj rect from 0, graph 0 to 20.739832999999997, graph 1
set label "chr22" at screen 0.01, graph 0.5
plot \
     "<awk '$1==\"chr22\"' CHM13-f1-90.bb.mini-inter-none.win" u ($2*1e-6):(0):($3) not w filledcu ls 1, \
     "<awk '$1==\"chr22\"' CHM13-f1-90.bb.mini-inter-none.win" u ($2*1e-6):($3):($3+$4) not w filledcu ls 2, \
     "<awk '$1==\"chr22\"' CHM13-f1-90.bb.mini-inter-none.win" u ($2*1e-6):($3+$4):($3+$4+$5) not w filledcu ls 3
set origin 0,0.01
set xran [0:154.259625]
set size 1.2420884856374232,0.06521739130434782
set style rect fc lt -1 fs solid 0.15 noborder
unset obj; unset label
set obj rect from 52.820107, graph 0 to 65.927026, graph 1
set label "chrX" at screen 0.01, graph 0.5
plot \
     "<awk '$1==\"chrX\"' CHM13-f1-90.bb.mini-inter-none.win" u ($2*1e-6):(0):($3) not w filledcu ls 1, \
     "<awk '$1==\"chrX\"' CHM13-f1-90.bb.mini-inter-none.win" u ($2*1e-6):($3):($3+$4) not w filledcu ls 2, \
     "<awk '$1==\"chrX\"' CHM13-f1-90.bb.mini-inter-none.win" u ($2*1e-6):($3+$4):($3+$4+$5) not w filledcu ls 3
