# MMI_Calc
Calculation of metrics for multi-metric index (MMI) for benthic macroinvertebrates or fish.  Inputs are 3 files (master taxa list (with phylogeny, autecological information, and any operation taxonomic units (OTUs), stations (with location and any regional information), and a taxa samples file (with StationID, SampleID, TaxonID, Count, and Exclude code).

Installation
------------
library(devtools)  #install if needed
install_github("leppott/MMIcalc")

Purpose
------------
Generate metric values for each provided sample.  Then calculate metric scores and combine into index score.  Scoring is based on a user defined table of IndexName, IndexRegion, MetricName, Direction in response to increasing perturbation, ScoringRegime (1/3/5 or 0:100), and High and Low values.

Status
------------
Development stage.  File input is assumed to already be aggregated.  Benthos and Fish are working.  Algae is a placeholder.  Not all metrics are included.  By default all metrics are calculated and if user doesn't specify metric names then all will be returned.

