CC=			gcc
CFLAGS=		-g -Wall -Wc++-compat -std=c99 -O2
CPPFLAGS=
INCLUDES=
OBJS=		kalloc.o kthread.o algo.o sys.o gfa-base.o gfa-io.o gfa-aug.o \
            sketch.o misc.o bseq.o options.o fastcmp.o shortk.o \
			index.o lchain.o gchain1.o gcmisc.o map-algo.o cal_cov.o \
			format.o gmap.o ksw2_extd2_sse.o ggsimple.o ggen.o
PROG=		minigraph
LIBS=		-lz -lpthread -lm

ifneq ($(asan),)
	CFLAGS+=-fsanitize=address
	LIBS+=-fsanitize=address
endif

.SUFFIXES:.c .o
.PHONY:all clean depend

.c.o:
		$(CC) -c $(CFLAGS) $(CPPFLAGS) $(INCLUDES) $< -o $@

all:$(PROG)

minigraph:$(OBJS) main.o
		$(CC) $(CFLAGS) $^ -o $@ $(LIBS)

clean:
		rm -fr gmon.out *.o a.out $(PROG) *~ *.a *.dSYM

depend:
		(LC_ALL=C; export LC_ALL; makedepend -Y -- $(CFLAGS) $(DFLAGS) -- *.c)

# DO NOT DELETE

algo.o: kalloc.h algo.h kvec.h ksort.h
bseq.o: bseq.h kvec.h kalloc.h kseq.h
cal_cov.o: mgpriv.h minigraph.h gfa.h gfa-priv.h algo.h kalloc.h
fastcmp.o: mgpriv.h minigraph.h gfa.h algo.h kalloc.h ksw2.h
format.o: kalloc.h mgpriv.h minigraph.h gfa.h
gchain1.o: mgpriv.h minigraph.h gfa.h ksort.h khash.h kalloc.h gfa-priv.h
gcmisc.o: mgpriv.h minigraph.h gfa.h kalloc.h
gfa-aug.o: gfa-priv.h gfa.h ksort.h
gfa-base.o: gfa-priv.h gfa.h kstring.h khash.h kalloc.h ksort.h
gfa-io.o: kstring.h gfa-priv.h gfa.h kseq.h
ggen.o: kthread.h kalloc.h sys.h bseq.h ggen.h minigraph.h gfa.h mgpriv.h
ggen.o: gfa-priv.h
ggsimple.o: mgpriv.h minigraph.h gfa.h gfa-priv.h kalloc.h bseq.h algo.h
ggsimple.o: sys.h ggen.h
gmap.o: kthread.h kalloc.h bseq.h sys.h mgpriv.h minigraph.h gfa.h gfa-priv.h
index.o: mgpriv.h minigraph.h gfa.h khash.h kalloc.h kthread.h kvec.h sys.h
kalloc.o: kalloc.h
ksw2_extd2_sse.o: ksw2.h kalloc.h
kthread.o: kthread.h
lchain.o: mgpriv.h minigraph.h gfa.h kalloc.h
main.o: mgpriv.h minigraph.h gfa.h gfa-priv.h sys.h ketopt.h
map-algo.o: kalloc.h mgpriv.h minigraph.h gfa.h khash.h ksort.h
misc.o: mgpriv.h minigraph.h gfa.h ksort.h
options.o: mgpriv.h minigraph.h gfa.h sys.h
shortk.o: mgpriv.h minigraph.h gfa.h ksort.h kavl.h algo.h khash.h kalloc.h
sketch.o: kvec.h kalloc.h mgpriv.h minigraph.h gfa.h
sys.o: sys.h
