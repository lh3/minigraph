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
