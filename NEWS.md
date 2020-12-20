Release 0.14-r415 (19 December 2020)
------------------------------------

Notable changes:

 * Added the `--call` option to find the allele/walk in each bubble.

 * Reduced the default minimum variant length (option `-L`) from 100 to 50 for
   the consistency with the SV community.

(19 December 2020, r415)



Release 0.13-r397 (3 December 2020)
-----------------------------------

Notable change:

 * Fixed incorrect anchors in linear chains. In older versions, a linear chain
   may contain two anchors with identical reference or query coordinates.

(3 December 2020, r397)



Release 0.12-r389 (26 October 2020)
-----------------------------------

Notable changes:

 * Improve alignments towards ends of graph segments. If there is an SV close to
   the ends but not at the ends, older versions may produce an excessively
   large bubble including high-identity matches.

 * Heuristically accelerates alignment in complex subgraphs by skipping
   many unnecessary sequence-aware graph traversals. This speeds up graph
   generation for CHM13 by three folds without obviously affecting accuracy.

 * Added option --inv to optionally disable inversions. Graph traversal is hard
   with inversions.

 * Fixed the bug that prevents large -K.

 * Apply option -K4g to the asm preset.

 * Added option --write-mz to output the positions of minimizer anchors.

(26 October 2020, r389)



Release 0.11-r371 (13 September 2020)
-------------------------------------

Notable changes:

 * Added option --max-rmq-size to limit the max RMQ size, which is set 100k by
   default. This heuristic reduces the long running time for aligning long
   centromeric sequences. The accuracy might be affected in rare cases.

 * Cap the max k-mer occurrence to 250 by default. For maize genomes, the
   current heuristic may choose an occurrence cutoff larger than 1000. This
   makes minigraph too slow to be practical.

 * Added option -S to output more detailed information about linear chains.

 * Added option -D to ignore diagonal minimizer anchors. This is useful to
   mapping a sequence against itself.

(13 September 2020, r371)



Release 0.10-r356 (14 February 2020)
------------------------------------

Notable changes:

 * Older releases miss a small fraction of INDELs involving repeats. This
   release fixes this issue.

 * Added the "stableGaf" command to mgutils.js to convert unstable GAF (e.g. by
   GraphAligner) to stable GAF.

(14 February 2020, r356)



Release 0.9-r343 (31 December 2019)
-----------------------------------

Notable changes:

 * RMQ based linear chaining. The chaining accuracy should be higher for large
   events. The speed remains similar.

 * Use ksw2 to check the sequence divergence of events to be inserted.

 * Treat inversions as special events. Don't insert them as long substitutions.

(31 December 2019, r343)



Release 0.8-r316 (11 December 2019)
-----------------------------------

This release reduces suboptimal chains caused by the chaining heuristics. It
generates slightly simpler human graphs.

(11 December 2019, r316)



Release 0.7-r310 (21 November 2019)
-----------------------------------

Notable changes:

 * Increased the default maximum INDEL/event length from 10kb to 100kb for
   assembly mapping and graph generation.

 * Decreased the default minimum INDEL/event length from 250bp to 100bp.

 * Accelerated graph mapping by pre-filtering isolated anchors and disconnected
   linear chains. This triples the performance when long gaps are desired.

Due to the change of default parameters, this release generates graphs
different from the previous versions.

(21 November 2019, r310)



Release 0.6-r302 (17 November 2019)
-----------------------------------

Notable changes:

 * Assign weight to seeds based on their repetitiveness. This helps chaining in
   repetitive regions a little bit.

 * For short-read mapping, prefer the reference path if the alternate path is
   not much better.

Major changes may be coming in the next release.

(17 November 2019, r302)



Release 0.5-r285 (8 September 2019)
-----------------------------------

Notable changes:

 * Fixed a bug that leads to wrong mapping positions in GAF.

 * Fixed two bugs related to graph chaining.

 * Added option `-j` to set expected sequence divergence and to adjust other
   chaining parameters accordingly.

 * Increased the k-mer thresholds for fast divergence estimate. This improves
   the alignment around low-complexity regions.

 * Tuned the default parameters to add highly divergent events only.

 * Warn about duplicated sequence names in graph construction (#3).

This version generates graphs different from the previous versions. The mapping
accuracy is improved due to the bug fixes and parameter tuning.

(8 September 2019, r285)



Release 0.4-r267 (22 August 2019)
---------------------------------

Notable changes:

 * Support paired-end mapping for short reads.

 * Remap and calculate coverage (see the new --cov option in the manpage).

 * Fixed multiple edges in the generated graphs. The v0.3 14-genome graph
   contains one multiple edge.

 * Use dynamic minimizer occurrence cutoff. For human data, the dynamic cutoff
   is around 137, higher than the default cutoff 100 used in earlier versions.
   As a result, graph generations will become a little slower.

Due to the last two changes, graphs generated with this version are different
from the previous versions.

(22 August 2019, r267)



Release 0.3-r243 (7 August 2019)
--------------------------------

This release generates graphs with SR tags on L-lines. The topology of the
graph is identical to the one generated with v0.2.

(7 August 2019, r243)



Release 0.2-r235 (19 July 2019)
-------------------------------

This release fixes multiple minor bugs. It also considers k-mer matches and
improves the accuracy of graph chaining. Nonetheless, the old chaining
algorithm, albeit simple, works quite well. The improvement is marginal.

(19 July 2019, r235)



Release 0.1-r191 (6 July 2019)
------------------------------

Initial proof-of-concept release.

(6 July 2019, r191)
