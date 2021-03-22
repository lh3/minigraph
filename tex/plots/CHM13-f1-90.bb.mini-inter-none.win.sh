gzip -dc CHM13-f1-90.bb.anno.gz | awk '$12~/mini/' | ./bedutils.js window -l CHM13v1.size -w500000 -s100000 -c /dev/stdin > CHM13-f1-90.bb.mini-win
gzip -dc CHM13-f1-90.bb.anno.gz | awk '$12~/inter|SINE|LINE|SVA|DNA|ERV/' | ./bedutils.js window -l CHM13v1.size -w500000 -s100000 -c /dev/stdin > CHM13-f1-90.bb.inter-win
gzip -dc CHM13-f1-90.bb.anno.gz | awk '$12~/none|partial|self/' | ./bedutils.js window -l CHM13v1.size -w500000 -s100000 -c /dev/stdin > CHM13-f1-90.bb.none-win

paste CHM13-f1-90.bb.mini-win CHM13-f1-90.bb.inter-win CHM13-f1-90.bb.none-win | awk '$1~/^chr([0-9]+|X)$/' | cut -f1-3,6,9 > CHM13-f1-90.bb.mini-inter-none.win

./chr-plot.js -n3 CHM13v1.cen.bed CHM13-f1-90.bb.mini-inter-none.win|gnuplot
