./mgutils.js anno -m -r GRCh38-464.bb.fa.out.gz -g GRCh38-464.bb.gap -d GRCh38-464.bb.sdust.gz -e GRCh38-464.bb.etrf.gz -p GRCh38-464.bb.paf.gz -c GRCh38-464.bb.brnn.gz -b GRCh38-464.bb.bed.gz -x ../GRCh38-464.seg-dust19.bed.gz -s ../GRCh38-464.seg-SD.bed.gz GRCh38-464.bb.base |gzip > GRCh38-464.bb.v2.anno.gz
./mgutils.js anno -m -r CHM13-464.bb.fa.out.gz -g CHM13-464.bb.gap -d CHM13-464.bb.sdust.gz -e CHM13-464.bb.etrf.gz -p CHM13-464.bb.paf.gz -c CHM13-464.bb.brnn.gz -b CHM13-464.bb.bed.gz -x ../CHM13-464.seg-dust19.bed.gz -s ../CHM13-464.seg-SD.bed.gz CHM13-464.bb.base |gzip > CHM13-464.bb.v2.anno.gz

paste call-GRCh38/*.bed|./mgutils.js merge -a anno/GRCh38-464.bb.v2.anno.gz -s 31sample.txt -|bgzip > GRCh38-464.call.bed.gz
paste call-CHM13/*.bed|./mgutils.js merge -a anno/CHM13-464.bb.v2.anno.gz -s 31sample.txt -|bgzip > CHM13-464.call.bed.gz
