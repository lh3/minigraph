Release 0.4-r264 (21 August 2019)
---------------------------------

Notable changes:

 * Support paired-end mapping for short reads.

 * Fixed multiple edges in the generated graphs.

 * Remap and calculate coverage (see the new --cov option in the manpage).

On the 14-genome dataset, this release produces a graph identical to the one
produced by v0.3, except the removal of one multiple edge.

(21 August 2019, r264)



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
