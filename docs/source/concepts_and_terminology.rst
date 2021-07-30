Concepts and Terminology
========================

Terminology
-----------

We refer to a set of genotypes with a given phenotype as a genotype set, which
is typically a very small subset of genotype space – the set of all possible
genotypes. Each genotype set therefore corresponds to a single phenotype. A
genotype set may comprise one or more genotype networks. In such networks,
vertices represent genotypes and edges connect vertices if their corresponding
genotypes are separated by a single small mutation. Vertices that share an edge
are referred to as neighbors. Since individual genotypes may belong to multiple
genotype sets, genotype networks may overlap.

We also construct and visualize phenotype networks. In such networks, vertices
represent genotype sets, and edges connect vertices if any genotypes in the
corresponding genotype sets can be interconverted via a single small mutation.
Vertices that share an edge are referred to as adjacent. We refer to mutations
that lead to a change in phenotype as non-neutral. For further information on
these and related concepts, the reader is referred to [1]_.

Analyses
--------

Please review the following important information that applies to analyses in
general:

* All analyses are performed on the dominant genotype network. That is, if the
  genotype set is fragmented into several genotype networks, analyses will be
  performed on the largest of these networks. All other components are ignored.
  The only exception is :ref:`evolvability <analysis_evolvability>`, for which
  if the :ref:`--use-all-components <use_all_components>` command-line
  argument is set, all connected components are considered.

* Genotypes with a score below the score threshold :ref:`tau <tau>` are not
  considered in any analysis.

* The noise threshold :ref:`delta <input_format_delta>` is only used in
  landscape related analyses, i.e., Paths, Peaks, and Epistasis (see below).

.. _analysis_evolvability:

Evolvability
^^^^^^^^^^^^

We define evolvability as the ability of mutation to bring forth a novel
phenotype. Evolvability can be measured at the scale of an individual genotype
or of a genotype set. In the Genonets Server, evolvability analyses are always
enabled, resulting in the following:

* Calculation of genotype evolvability
* Calculation of phenotype evolvability
* Construction of the phenotype network

These computational routines are based on [2]_, and are described in detail
below.

.. _analysis_genotype_evolvability:

Genotype evolvability
"""""""""""""""""""""

For a genotype g in a genotype set :math:`S`, evolvability is the ratio of the
number of genotype sets to which :math:`g` can evolve via a single mutation, to
the total number of genotype sets in the input data.

The higher the evolvability of a genotype :math:`g`, the higher the number of
genotype sets that can be reached by a single mutation from :math:`g`.

.. _analysis_phenotype_evolvability:

Phenotype evolvability
""""""""""""""""""""""

For a genotype set :math:`S`, phenotype evolvability is the ratio of the number
of unique genotype sets to which genotypes in :math:`S` can evolve via single
mutations, to the total number of genotype sets available in the input data.

.. _analysis_phenotype_network:

Phenotype network
"""""""""""""""""

* Each vertex in the phenotype network represents a genotype set.
* The higher the out-degree of a vertex, the higher the number of genotype sets
  to which genotypes in this genotype set can evolve.
* The higher the in-degree of a vertex, the higher the number of genotype sets
  from which this genotype set is accessible.
* The network is a directed graph, which captures the possibly asymmetric
  relation between vertices. This means that it is possible for a genotype set
  to be accessible from other genotype sets, despite having zero phenotype
  evolvability itself.

.. _analysis_robustness:

Robustness
^^^^^^^^^^

We define robustness as the invariance of a phenotype in the face of genetic
perturbation. Like evolvability, robustness can be measured at the scale of an
individual genotype or of a genotype set. In Genonets, robustness
analysis results in the following computations:

* Genotype robustness
* Average robustness of the genotype set

These computational routines are based on [2]_, and are described in detail
below.

.. _analysis_genotype_robustness:

Genotype robustness
"""""""""""""""""""

For a genotype :math:`g` in a genotype set :math:`S`, robustness is the
fraction of all possible mutational neighbors that are also in :math:`S`. Thus,
:math:`g` is maximally robust if all possible neighbors are members of
:math:`S`.

.. _analysis_phenotype_robustness:

Phenotype robustness
""""""""""""""""""""

Phenotype robustness is the arithmetic mean of the genotype robustness values
for all genotypes in the genotype set :math:`S`.

.. _analysis_accessibility:

Accessibility
^^^^^^^^^^^^^

The accessibility of a genotype set :math:`S` measures the potential for
mutation to generate a genotype in :math:`S` from genotypes in different
genotype sets. In Genonets, accessibility is measured for all genotype sets in
the input data. Specifically, for a genotype set :math:`S`, accessibility is
computed as follows:

#. For each pair of genotype sets :math:`(S, T)`, calculate the ratio of the
   number of genotypes in :math:`S` that are separated by a single mutation from
   any genotype in :math:`T`, to the total number of genotypes that are
   separated by a single mutation from genotypes in :math:`T`.
#. Then calculate the accessibility of :math:`S` as the sum of these ratios for
   all pairs :math:`(S, T)`.

Computational routines for accessibility are based on [3]_.

.. _analysis_neighbor_abundance:

Neighbor Abundance
^^^^^^^^^^^^^^^^^^

The neighbor abundance of a genotype set :math:`S` measures the size of adjacent
genotype sets, in proportion to the probability that a mutation will generate a
genotype in these adjacent genotype sets. In Genonets, neighbor abundance is
measured for all genotype sets in the input data. Specifically, for a genotype
set :math:`S`, neighbor abundance is computed as follows:

#. Calculate the ratio of the number of genotypes in :math:`T` that are
   accessible from :math:`S`, to the total number of genotypes that are
   accessible from :math:`S`.
#. Multiply this ratio by the number of genotypes in :math:`S`.
#. Repeat this process for all genotype set pairs :math:`(S, T)`, taking the sum
   as the neighbor abundance of :math:`S`.

Computational routines for neighbor abundance are based on [3]_.

.. _analysis_diversity_index:

Diversity Index
^^^^^^^^^^^^^^^

The diversity index of a genotype set :math:`S` gives the probability that two
randomly chosen non-neutral mutations to genotypes in :math:`S` yield genotypes
that belong to the same genotype set :math:`T`. In Genonets, diversity index is
measured for all genotype sets in the input data. Specifically, the diversity
index of a genotype set :math:`S` is computed as follows:

#. Calculate the ratio of the number of genotypes in :math:`T` that are
   accessible from :math:`S`, to the total number of genotypes that are
   accessible from :math:`S`.
#. Square this ratio.
#. Repeat this process for all genotype set pairs :math:`(S, T)`, summing up
   along the way.
#. The diversity index of :math:`S` is one minus this sum.

Computational routines for diversity index are based on [3]_.

.. _analysis_structure:

Structure
^^^^^^^^^

The last two decades of research in network science have produced a wealth of
measures for describing the structure of networks. Genonets includes many of
these analyses. These result in measures at the level of individual genotypes
and genotype sets.

Computations performed at the level of the genotype set are:

* Number of connected components, i.e., number of genotype networks within a
  single genotype set
* Sizes of all connected components
* Size of the giant component, i.e., size of the dominant genotype network
* Proportional size of the dominant genotype network
* Diameter of the dominant genotype network
* Edge density of the dominant genotype network
* Average clustering coefficient for the dominant genotype network

Computations performed at the level of genotypes are:

* Coreness
* Clustering coefficient
* Computational routines for structural analysis are described in [4]_.

.. _analysis_overlap:

Overlap
^^^^^^^

Since some genotypes belong to more than one genotype set, genotype networks
sometimes overlap. Using the overlap analysis, Genonets will characterize these
regions of overlap for all pairs of genotype sets. Specifically, for each pair
of genotype sets :math:`(S, T)` available in the input data, this analysis
calculates the number of genotypes that are common to both genotype sets
:math:`S` and :math:`T`.

.. _analysis_epistasis:

Epistasis
^^^^^^^^^

Epistasis – non-additive interactions between mutations – can impose severe
constraints on molecular evolution because the mutations that are beneficial in
one genetic background may be deleterious in another. Epistasis can be
classified as *magnitude*, *simple sign*, or *reciprocal sign* epistasis
depending on the sign (i.e., positive or negative) of the individual mutations
and of the mutations in combination (please see [5]_ for details). In Genonets,
epistasis analysis results in the following calculations:

#. Identify all squares in the dominant genotype network, as these represent
   pairs of mutations.
#. For each square, determine the class of epistasis (magnitude, simple sign,
   reciprocal sign).
#. For each epistasis class, calculate the proportion of all squares in the
   dominant genotype network that belong to this class.

Computational routines for epistasis are based on [5]_.

.. _analysis_peaks:

Peaks
^^^^^

In the input data, the user is required to provide a :ref:`input_format_score`
for each genotype. Since these scores may reflect a quantitative phenotype that
is related to organismal fitness, and because these scores vary amongst the
genotypes in a genotype network, one may think of a genotype network as an
adaptive landscape [6]_. This opens the door to a slew of analyses that
characterize the potential for mutation and selection to explore these
landscapes. One of these analyses comprises determination of peaks in the
landscape.

Peaks analysis in Genonets results in the determination of the global and all
local peaks in the landscape. We refer to the genotype with the highest score
in the genotype network as the *summit*. Please note that even though there can
be multiple genotypes within a peak, when referring to the global peak within
the Genonets documentation, we are in fact referring to the summit.

Computational routines for peaks are based on [5]_.

.. _analysis_paths:

Paths
^^^^^

Another analysis where the genotype network is considered an adaptive
landscape [6]_ (see the introduction to Peaks analysis above) is the computation
of accessible mutational paths.

The paths analysis involves computing all accessible mutational paths from each
genotype in the network, to the summit. A path is accessible, if and only if the
scores for the genotypes on the path increase monotonically (plus or minus the
user-supplied parameter delta), from the source genotype to the target genotype.

Computational routines for paths are based on [5]_.

.. rubric:: References

.. [1] Andreas Wagner. Neutralism and selectionism: a network-based
   reconciliation. Nature Reviews Genetics 9, 965-974 (December 2008).

.. [2] Andreas Wagner. Robustness and evolvability: a paradox resolved. Proc.
   R. Soc. B 2008 275 91-100; DOI: 10.1098/rspb.2007.1137. Published 7 January
   2008.

.. [3] Cowperthwaite MC, Economo EP, Harcombe WR, Miller EL, Meyers LA (2008)
   The Ascent of the Abundant: How Mutational Networks Constrain Evolution.
   PLoS Comput Biol 4(7): e1000110. doi:10.1371/journal.pcbi.1000110.

.. [4] Mark Newman. Networks: An Introduction. Oxford University Press, Inc.,
   New York, NY, USA. (2010).

.. [5] Jose Aguilar Rodriguez, Joshua L. Payne, Andreas Wagner One thousand
   adaptive landscapes and their navigability. In review.

.. [6] Sewall Wright. The roles of mutation, inbreeding, crossbreeding and
   selection in evolution. In Proc. Sixth Int. Congr. Genet. 356–366 (1932).