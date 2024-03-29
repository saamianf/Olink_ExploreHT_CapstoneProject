# Olink_ExploreHT_CapstoneProject

Part 1 - LOD Evaluation
1. Read in HT csv file.
2. Determine LOD cutoff value for each assay.
3. Determine whether each NPX value for the samples in each assay fall below the LOD cutoff --> output is 'PASS' or 'WARNING'.
4. Evaluate rates of warnings bsaed on assay, sample, and block groupings.
5. Output a HT LOD check csv file.

Part 2 - General Figures
1. Read in HT LOD check csv file using NPX read-in.
2. Read in manifest xlsx file.
3. Apply pediatric + peripheral blood filter.
4. Generate heatmaps, NPX distribution plots, QC plots, and UMAP plots.
5. Run full t-test + generate volcano plot.
6. Generate LOD-colored volcano plot for full t-test.
7. Run t-tests per time-point pairs + generate volcano plots.

Part 3 - LOD Figures
1. Read in pediatric filtered dataset from part 2.
2. Generate summary graphs indicating LOD cutoff + NPX distribution for each assay.

Part 4 - Olink Cross-Panel Comparisons
1. Read in Target96 NPX file + associated manifest.
2. Read in Flex48 NPX file + associated manifest.
3. Apply pediatric + peripheral blood filter.
4. AdjustedNPX calculations.
5. Iterate through each assay.
6. Iterate through each subject to generate line plots for NPX.
7. Generate histograms + bar graphs by assay.

![image](https://github.com/saamianf/Olink_ExploreHT_CapstoneProject/assets/120279968/4bbe53e0-7281-41d2-a91c-8ec8bdc191ab)
