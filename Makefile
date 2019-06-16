CC=			gcc
CFLAGS=		-g -Wall -Wc++-compat -O2
CPPFLAGS=
INCLUDES=	-I.
OBJS=		kalloc.o kthread.o gfa-base.o gfa-io.o gfa-sub.o gfa-aug.o \
            sketch.o misc.o mss.o bseq.o options.o \
			index.o lchain.o gchain1.o gcmisc.o map.o ggen.o ggsimple.o format.o
PROG=		minigraph
LIBS=		-lz -lm

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

bseq.o: bseq.h kvec.h kalloc.h kseq.h
format.o: kalloc.h mgpriv.h bseq.h minigraph.h gfa.h
gchain1.o: mgpriv.h kalloc.h bseq.h minigraph.h gfa.h ksort.h khash.h
gcmisc.o: mgpriv.h kalloc.h bseq.h minigraph.h gfa.h
gfa-aug.o: gfa-priv.h gfa.h ksort.h
gfa-base.o: gfa.h khash.h kalloc.h ksort.h
gfa-io.o: kstring.h gfa.h kseq.h
gfa-sub.o: gfa.h kalloc.h kavl.h khash.h ksort.h
ggen.o: kthread.h kalloc.h mgpriv.h bseq.h minigraph.h gfa.h
ggsimple.o: mgpriv.h kalloc.h bseq.h minigraph.h gfa.h ksort.h mss.h
index.o: mgpriv.h kalloc.h bseq.h minigraph.h gfa.h khash.h kthread.h kvec.h
kalloc.o: kalloc.h
kthread.o: kthread.h
lchain.o: mgpriv.h kalloc.h bseq.h minigraph.h gfa.h
main.o: bseq.h mgpriv.h kalloc.h minigraph.h gfa.h ketopt.h
map.o: kthread.h kvec.h kalloc.h mgpriv.h bseq.h minigraph.h gfa.h khash.h
map.o: ksort.h
misc.o: mgpriv.h kalloc.h bseq.h minigraph.h gfa.h ksort.h
mss.o: kvec.h kalloc.h mss.h
options.o: minigraph.h gfa.h
sketch.o: kvec.h kalloc.h mgpriv.h bseq.h minigraph.h gfa.h
