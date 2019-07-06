## Getting Started

```sh
git clone https://github.com/lh3/minigraph
cd minigraph && make
# Map sequence to sequence
./minigraph test/MT-human.fa test/MT-orangA.fa > out.paf
# Map sequence to graph
./minigraph test/MT.gfa test/MT-orangA.fa > out.gaf
# Incremental graph generation (-l10k necessary for this toy example)
./minigraph -xggs -l10k test/MT.gfa test/MT-chimp.fa test/MT-orangA.fa > out.gfa
```

## Introduction
