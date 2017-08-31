# Multigroup analysis

Abstract: Idea is to provide an easy to use API to process and analyze Quantitative mass spectrometry data from different sources.

## Aims

- Single configuration script allowing to specify imputation, normalization methods, models used, visualizations
- Generate output Report, including figure captions, in HTML and PDF format. Place all figures with reasonable naming in high res png and pdf in folders.
- Allow for variable column naming
     - change column names depending on input software and allow for meaningful column names for factors.
- Ensure that naming is consistent through analysis.
- Report p-values and adjusted p-values but more importantly confidence intervals.
- Provide statistics of fold changes in case of multiple comparisons:
     - Barplots, correlation plots, Venn diagrams.
- Support data on protein, peptide/precursor and transition level (likely peptide/precursor level only).
     - use filtering and aggregation on precursor level if statistics backend does not support those levels (e.g. Limma)
- Support algorithms working on the precursor (only MSStats and MapDIA), peptide (msqrob) and protein level (limma).
- Visualize up to 3 explanatory variables e.g. condition, gender, patient. - This also means that the plots need to be annotated.
- define comparisons using contrasts i.e. condition1 vs condition2
- define hierarchy of conditions i.e. test >> patient >> gender / MainCondition variable
- define set of QC visualizations for each level - transition, precursor and protein.
- model fractionated experiments
- create analysis compatible summaries for the analysis so they can be than compared.
- Aggregate first then fold changes vs fold changes first than statistics on fold changes.
     - $log2(a1+b1+c1) - log(a2+b2+c2)$ vs. $log(a1/a2), log(b1/b2), log(c1/c2)$
- Compare results. There might be different results do to change in parameters.
    - Study parameter sensitivity - which step in the analysis changes the result most.
    - Benchmark datasets.

## Functionality:

* Data processing summaries
   * How many conditions, How many samples per condition (2x2 tables). for implementation see [Pass generic column names to xtabs function in R](https://stackoverflow.com/questions/31992301/pass-generic-column-names-to-xtabs-function-in-r)
```
ftable(xtabs( ~ Time + Treatment, data= annotation))
ftable(xtabs( ~ Strain + Time + Treatment + Plant, data= annotation))
```

   * How many transitions, proteins, peptides, precursors at each filtering step.

* generate QC plots for the data:
   * Scatter plot within condition
       - on which level - Precursor, protein
       - on which condition if up to 3?
   * CV for each condition if N > 2 and overall.
       -  How to interpret it in a paired experiment? Maybe disable it? Or show 2 CV plots based on each condition variable.
   * Distribution plot for each file - Violin plot or density?
       * Densities within condition - before and after normalization
       * compare condition based on common density - before and after normalization
   * Correlation plot for samples - on which level - precursor protein or on all?
   * Clustering on not normalized and normalized data with condition labels
       - specify which variables to encode - work out how to encode - color scheme.
   * Summaries for NAs
   * QC for Retention time only possible on transition and precursor level data.
   * QC of data filtering
        - use sample correlation
        - use transition and peptide correlation
    * Pairsplot of group averages (output of modeling actually)

* Data preprocessing:
     * data filtering
         - based on intensity : top X transitions per precursor, top X precursors per protein
         - correlation
              - can be performed on normalized and not normalized data. What is the consequence?
              - Measure which combines correlation and scale. i.e. it does not matter if transitions not correlated if observed fold change are very small.
         - NA counting and removal, globally and per condition - how to specify it?
         - min X precursors per protein
         - for DIA data - Q value filtering using 2 thresholds.
         - remove / keep not proteotypic peptides
     * data aggregation (protein from precursors, and precursors from transitions)
         - using sum, mean median
     * data normalization
         - vsn
         - median
         - median and variance
         - quantile
         - Tuckey median polish, what for: https://www.youtube.com/watch?v=RtC9ZMOYgk8
     * missing value imputation - must be used if aggregating is enabled
         - row mean and column mean imputation
         - linear regression
         - hot deck encoding.

* Data modeling:
     * ANOVA
     * linear models
     * mixed linear models
     * external packages like _limma_ and _msqrob_ (eventually MSStats) on the data.

* Visualization of modeling results.
     * Global for all proteins:
         * ma plot for each comparison - color code p-values?
         * volcano plot for each comparison (p-values, adjusted p-values)
         * scatterplot of fold changes for 2 comparisons
         * distribution of p-values and distribution of adjusted p-values (FDR)
     * Local (for each protein)
         * boxplot for each protein (unpaired)
         * line plot showing all transitions.
         * Show confidence intervals.


## Protein inference problem

We measure precursor but the subject is proteins. Precursors can be grouped in proteins. Either, unambiguously (N:1) no conflicting peptide protein assignments. Optionally we can allow for ambiguity (N:M), which means that a peptide can be assigned to more than one protein.
- Visualize profiles of proteins which have shared peptides.

## Implementation

- Use long format (tibble) as long as possible and implement functions to switch from wide to long.
- methods work on _tibble's_ and should support _magrittr_ operator
- Use tibbles to represent transition, precursor and protein/peptide.
     - those tibbles should have types since some methods are only applicable to transitions, precursors, peptides and proteins.
- Optional tibble to represent annotation if annotation largish.
- allow for switching between long and wide format (if wide store long format in an tibble attribute).


## Configuration

### What variables are needed for filtering?

### What configuration is needed for modeling?

We want to learn about proteins, therefore the input needs to be:

      ProteinID, fixed effects, random effects, Response, Obsolete columns.

The ProteinID is here the grouping variable. Obsolete columns can be e.g. __Fasta.headers__ containing additional information. But obsolete can be also original intensities if transformed intensities are used. The transformed intensity will be the response. We limit the number of fixed effects to 3.

The specification needs to be

* subject = "ProteinID"
* FixedEffects = c("a","b","c")
* RandomEffects = c("x","y","z")
* Response = "Intensity"

* Formula for fixed and random effects or design matrix if only fixed effects are considered.
* A list of contrasts based on the 'main' Condition.

# What configuration is needed for plotting?

Most of the plots will work with the long format (ggplot2) but some of the plots need the data in a wide format (e.g. pairs plot).

## Diagnostic plots

* Pairsplot within condition - works only with 2 to 4 samples. Group labels required.
* Plot - Matrix (line) plot of transitions, peptide intensities per sample and protein.  No group labels required

* Correlation map samples.
* Violin plot showing CV for each condition.

Visualization of modeling results.

* Altman Bland plot - mean vs. log2 fold change ideally with the labeling of significant proteins.
* Distribution of p-values with description.
* Volcano plot

Is modeling diagnostics for single proteins needed?

* boxplot for single protein when intensity of protein level only
* box plot matrix plot for each peptide and precursor if available.

## Precursors are aggregated to proteins passing several levels

* Precursor
* Modified peptide sequence
* Peptide
* Protein

## Output

* modelling results
       Subject, Contrast, fold change, p-value, adjusted-pvalue, obsolete columns (i.e fasta.headers).

* transformation results:
       Subject, Sample, Intensity

## Evaluation

* How filtering is influenced by the number of samples/conditions?
* How one filtering step can influence a different filtering step? E.g. stronger qValue filtering might in the increase number of proteins left because more are kept after correlation filtering.

# Example Projects

## Example project 2147

factors - levels

* Raw.file
* Measurement.Date
* Measurement.Order
* Strain
* Timepoint
* Treatment
* Plant
* Leave
* BATCH

## fixed effects

interesting to the analysis

* Condition (Strain + Timepoint + Treatment)
* Since the experiment is paired we also need plant as an additional factor

## Random effects

Of no interest to the analysis, in this case, this would be the batch. If we move to protein level it will be peptide or precursor.

## Contrasts

Are computed based on a single factor - condition. Which can be a composition of several factors (see above).

## Protein Level Annotation

We are interested in inferences on the protein or protein groups level usually.

* Nr.Peptides
* Majority.protein.IDs (all proteins belonging to a protein group. See protein grouping  <http://aggrivet.blogspot.ch/2016/12/protein-clusters.html>.
* Top.Protein.Name - label protein of the protein group
* Fasta.headers (additional protein annotation)


## in proteomics we are usually measuring precursors:

* FileName
* TransitionGroupID
* StrippedSequence
* ModifiedSequence
* PrecursorCharge
* Decoy
* LabelType
* Measured.PrecursorMZ
* Measured.PrecursorRT
* Measured.PrecursorScore (some arbitrary quality score).
* Measured.MS1Intensity
* Measured.MS2IntensityAggregated

# Comparable concepts:

![ms QRob Benchmarking](inst/figures/msQRobBenchmark.gif)
![TAPR Package Overview](inst/figures/TAPROverviewImage.jpg)



# References

[TRAPR: R Package for Statistical Analysis and Visualization of RNA-Seq Data](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5389949/)
[Visualization of proteomics data using R and Bioconductor](http://onlinelibrary.wiley.com/doi/10.1002/pmic.201400392/full)
[limma](http://bioconductor.org/packages/release/bioc/html/limma.html)
[MSStats](http://msstats.org/)

[Summarization vs Peptide-Based Models in Label-Free Quantitative Proteomics: Performance, Pitfalls, and Data Analysis Guidelines](http://pubs.acs.org/doi/abs/10.1021/pr501223t)
[MSqRob github](https://github.com/statOmics/MSqRob)
[Benchmark data for MSqRob](https://github.com/statOmics/MSqRobData/)
