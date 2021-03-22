prefix=CHM13-f1-90.bb

all:$(prefix).base $(prefix).gap $(prefix).brnn.gz $(prefix).etrf.gz $(prefix).sdust.gz

$(prefix).base:$(prefix).fa
	seqtk comp $< | cut -f1,2 | tr "_" "\t" > $@

$(prefix).gap:$(prefix).fa
	seqtk gap $< > $@

$(prefix).brnn.gz:$(prefix).fa
	~/dna-nn/dna-brnn -Ai ~/dna-nn/attcc-alpha.knm -t16 $< | htsbox bgzip > $@

$(prefix).etrf.gz:$(prefix).fa
	~/src/etrf/etrf $< | htsbox bgzip > $@

$(prefix).sdust.gz:$(prefix).fa
	~/minimap2/sdust $< | htsbox bgzip > $@

CHM13-f1-90.bb.paf.gz:CHM13-f1-90.bb.fa
	minimap2 -cxasm20 -r2k --cs -t16 ~/ref/CHM13v1Y.fa $< 2> CHM13-f1-90.bb.paf.log | gzip > $@

GRCh38-f1-90.bb.paf.gz:GRCh38-f1-90.bb.fa
	minimap2 -cxasm20 -r2k --cs -t16 ~/ref/hs38.fa $< 2> GRCh38-f1-90.bb.paf.log | gzip > $@
