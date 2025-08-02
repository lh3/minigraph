prefix=CHM13-464.bb

all:$(prefix).base $(prefix).gap $(prefix).brnn.gz $(prefix).etrf.gz $(prefix).sdust.gz

$(prefix).fa:$(prefix).bed.gz
	gzip -dc $< | awk '{print ">"$$1"_"$$2"_"$$3;print $$14}' > $@

$(prefix).base:$(prefix).fa
	seqtk comp $< | cut -f1,2 | perl -ane 'print "$$1\t$$2\t$$3\t$$F[1]\n" if $$F[0]=~/(\S+)_(\d+)_(\d+)/' > $@

$(prefix).gap:$(prefix).fa
	seqtk gap $< > $@

$(prefix).brnn.gz:$(prefix).fa
	~/src/dna-nn/dna-brnn -Ai ~/src/dna-nn/attcc-alpha.knm -t16 $< | htsbox bgzip > $@

$(prefix).etrf.gz:$(prefix).fa
	~/src/etrf/etrf $< | htsbox bgzip > $@

$(prefix).sdust.gz:$(prefix).fa
	~/src/sdust/sdust $< | htsbox bgzip > $@

CHM13-464.bb.paf.gz:CHM13-464.bb.fa
	minimap2 -cxasm20 --cs -t16 ~/ref/chm13v2.fa $< 2> CHM13-464.bb.paf.log | gzip > $@

GRCh38-464.bb.paf.gz:GRCh38-464.bb.fa
	minimap2 -cxasm20 --cs -t16 ~/ref/hs38.fa $< 2> GRCh38-464.bb.paf.log | gzip > $@
