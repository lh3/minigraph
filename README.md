## Getting Started

```sh
git clone https://github.com/lh3/minigraph
cd minigraph && make
# Map sequence to sequence, similar to minimap2 without base alignment
./minigraph test/MT-human.fa test/MT-orangA.fa > out.paf
# Map sequence to graph
./minigraph test/MT.gfa test/MT-orangA.fa > out.gaf
# Incremental graph generation (-l10k necessary for this toy example)
./minigraph -xggs -l10k test/MT.gfa test/MT-chimp.fa test/MT-orangA.fa > out.gfa
```

## Introduction

<img align="right" width="278" src="doc/example1.png"/>

Minigraph is a *proof-of-concept* sequence-to-graph mapper and graph
constructor. It finds approximate locations of a query sequence in a sequence
graph and incrementally augments an existing graph with long query subsequences
diverged from the graph. It can construct a graph from 15 human assemblies in
an hour using 24 CPU cores.

Minigraph is at an early development stage. It lacks important features and may
produce suboptimal mappings. Please read the [Limitations](#limit) section of
this README before using minigraph.

## User's Guide

To install minigraph, type `make` in the source code directory. The only
non-standard dependency is [zlib][zlib].

## Algorithm Overview

In the following, minigraph command line options have a dash ahead and are
highlighted in bold. The description may help to tune minigraph parameters.

1. Read all reference bases, extract (**-k**,**-w**)-minimizers and index them
   in a hash table.

2. Read **-K** [=*500M*] query bases in the mapping mode, or read all query
   bases in the graph construction mode. For each query sequence, do step 3
   through 5:

3. Find colinear minimizer chains using the [minimap2][minimap2] algorithm,
   assuming segments in the graph are disconnected. These are called *linear
   chains*.

4. Perform another round of chaining, taking each linear chain as an anchor.
   For a pair of linear chains, minigraph finds up to 15 shortest paths between
   them and chooses the path of length closest to the distance on the query
   sequence. Importantly, sequences are ignored in this round of chaining.
   Chains found at this step are called *graph chains*.

5. Identify primary chains and estimate mapping quality with a method similar
   to the one used in minimap2.

6. In the graph construction mode, collect all mappings longer than **-d**
   [=*10kb*] and keep their query and graph segment intervals in two lists,
   respectively.

7. For each mapping longer than **-l** [=*50k*], finds poorly aligned regions.
   A region is filtered if it overlaps two or more intervals collected at step
   6.

8. Insert the remaining poorly aligned regions into the input graph. This
   constructs a new graph.

## <a name="limit"></a>Limitations

* Minigraph needs to find strong colinear chains first. For a graph consisting
  of many short segments (e.g. one generated from rare SNPs in large
  populations), minigraph will fail to map query sequences.

* When connecting colinear chains on graphs, minigraph ignores sequences, and
  only considers the distances of top 15 shortest paths between colinear
  chains. It may miss the optimal alignments.

* Minigraph doesn't give base-level alignment.

[zlib]: http://zlib.net/
[minimap2]: https://github.com/lh3/minimap2
