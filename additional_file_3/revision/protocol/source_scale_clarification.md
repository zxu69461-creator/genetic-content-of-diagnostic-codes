# Source-effect scale clarification

Added during source verification on 6 October 2026 after the first targeted results were computed. No loci, effect values, estimator, test direction or decision threshold were changed.

The revision plan called IGAS effects log-odds effects. The original GWAS Catalog file supplies beta, standard error, OR and OR confidence limits, and its visible records satisfy beta = log(OR). However, the primary IGAS study used FaST-LMM for its principal association analysis. A catalog column label alone does not establish comparability with adjusted logistic-regression effects from another cohort. The manuscript therefore uses the more precise description "released IGAS effects" and defines the alignment coefficient as target adjusted log-odds effect per unit of this fixed source effect. No conversion to clinical odds ratios or patient fractions is asserted. The source header and public acquisition record are retained.

The historical AS-equivalent ratio divides two loadings on this same source vector and is invariant to a common multiplicative rescaling of all source effects and their standard errors. This invariance does not remove variant-specific scale differences, source selection, ascertainment or population differences. It is not evidence for a clinical mixture interpretation.

Sources: the GWAS Catalog GCST005529 original summary file; Cortes and IGAS, Nature Genetics 2013, doi:10.1038/ng.2667; the original ImmunoBase study description https://www.immunobase.org/page/overview/display/study_id/gdxhss00026/.
