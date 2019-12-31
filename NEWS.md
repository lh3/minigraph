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
