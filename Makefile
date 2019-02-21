CC=			gcc
CFLAGS=		-g -Wall -Wc++-compat -O2
CPPFLAGS=
INCLUDES=	-I.
OBJS=		kalloc.o kthread.o gfa-base.o gfa-io.o sketch.o misc.o options.o bseq.o index.o lchain.o map.o
PROG=		minigraph
LIBS=		-lz

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

gfa-base.o: gfa.h khash.h kalloc.h ksort.h
gfa-io.o: kstring.h gfa.h kseq.h
index.o: mgpriv.h kalloc.h minigraph.h gfa.h khash.h kthread.h kvec.h
kalloc.o: kalloc.h
kthread.o: kthread.h
main.o: bseq.h minigraph.h gfa.h mgpriv.h kalloc.h ketopt.h
misc.o: mgpriv.h kalloc.h minigraph.h gfa.h ksort.h
options.o: minigraph.h gfa.h
sketch.o: kvec.h kalloc.h mgpriv.h minigraph.h gfa.h
