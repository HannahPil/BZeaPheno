# BZeaPheno
For plotting and analyzing phenotype data taken on the BZea population. Phenotypes were taken at Central Crops Research Station in Clayton, NC in years 2023 and 2025.

## Project structure

```
data/           # input data (gitignored)
output/         # generated plots and tables (gitignored)
scripts/        # all R analysis scripts
```

## Trait descriptions

DOS:	Date of silking; on which 50% of plants had visible silks.

DOA:	Date of anthesis; on which 50% of plants began shedding pollen on both central and lateral tassel spikes.

PH:	Plant height; measured from base of plant to tip of highest tassel.

EH:	Ear height; measured from base of plant to primary ear bearing node.

EN:	Ear number; number of nodes with ears per plant.

Prolif:	Prolificacy; total number of ears per plant.

NBR:	Number of nodes with brace roots per plant.

SPAD:	Leaf greenness index taken by SPAD meter on one leaf above primary ear bearing node.

LAE:	Number of leaves above primary ear.

## Scripts

| Script | Description |
|--------|-------------|
| `CLY23_phenotype_analysis.Rmd` | Main phenotype analysis notebook |
| `SpATS_analyze.R` | Spatial analysis of field trials using SpATS |
| `broad_sense_heritability.R` | Broad-sense heritability estimation |
| `upd_heritability.R` | Updated heritability analysis |
| `height_sister_jitter.R` | Sister-line jitter plots for plant height by introgression status at a focus gene, plus height-vs-expression correlations |
| `drone_vs_fieldbook_height.R` | Correlation of drone-derived heights (8/24) vs hand-measured fieldbook PH |
| `introgression_frequency_heatmap.R` | Heatmap of introgression frequencies across the genome |
| `inv4m_analysis.R` | Analysis of the Inv4m inversion region |
| `octexp_analysis.R` | Octoploid expression analysis |
| `jitter.R` | General jitter plots for trait comparisons |
| `N_jitter.R` | Jitter plots for ear number (N) |
| `beeswarms.R` | Beeswarm plots for trait distributions |
| `beeswarm_N.R` | Beeswarm plots for ear number |
| `N_actual_v_predicted.R` | Actual vs predicted ear number comparison |
| `violin.R` | Violin plots for trait distributions |
| `ridges.R` | Ridge plots for trait distributions |
| `histograms.R` | Histogram plots |
| `SPAD_density.R` | Density plots for SPAD values |
| `correlation_curves.R` | Trait correlation curves |
| `percentile_test.R` | Percentile-based statistical tests |
