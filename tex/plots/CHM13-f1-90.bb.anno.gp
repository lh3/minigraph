#set t pdfcairo transparent enh font "Helvetica,15"
set t po eps co so enh "Helvetica,18"

set style line 1 lt 1 lc rgb "#FF0000" lw 1;
set style line 2 lt 1 lc rgb "#00C000" lw 1;
set style line 3 lt 1 lc rgb "#0080FF" lw 1;
set style line 4 lt 1 lc rgb "#C000FF" lw 1;
set style line 5 lt 1 lc rgb "#00EEEE" lw 1;
set style line 6 lt 1 lc rgb "#FF80FF" lw 1;

set style line 1 lt 1 lc rgb "#fbb4ae" lw 1;
set style line 2 lt 1 lc rgb "#b3cde3" lw 1;
set style line 3 lt 1 lc rgb "#ccebc5" lw 1;

set out "CHM13-f1-90.bb.anno.cnt.eps"

set size 1,0.9

set style histogram rowstacked
set xtics rotate by 40 right nomirror font "Helvetica,18"
set boxwidth 0.8 relative
set style data histograms
set style fill solid 1.0 border lt -1
#set style fill pattern 7 border lt -1
set ylab "Count ({/Symbol \264}10^3)" off +0.0,0
set bmargin 5
set lmargin 8

set title "CHM13 minigraph (CHM13 +GRCh38 +44 samples)"
plot \
	"<cat CHM13-f1-90.bb.anno.tbl" u ($3*1e-3):xtic(2) t '2 alleles' ls 1, \
	"" u ($4*1e-3) t '3 alleles' ls 3, \
	"" u ($5*1e-3) t '>3 alleles' ls 2

set out "CHM13-f1-90.bb.anno.len.eps"

set ylab "Sum of length on reference (Mbp)" off +0.0,0
set key top left
plot \
	"<cat CHM13-f1-90.bb.anno.tbl" u ($6*1e-6):xtic(2) t '2 alleles' ls 1, \
	"" u ($7*1e-6) t '3 alleles' ls 3, \
	"" u ($8*1e-6) t '>3 alleles' ls 2
