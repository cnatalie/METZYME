# METZYME
R code and data files used to create heatmaps, ordinations, normalizations, the 18S phyloseq and WCGNA analyses used in the study:

## Dinoflagellates alter their carbon and nutrient metabolic strategies across environmental gradients in the central Pacific Ocean

by Natalie R. Cohen, Matthew R. McIlvin, Dawn M. Moran, Noelle A. Held, Jaclyn Saunders, Nicholas Hawco, Michael Brosnahan, Giacomo DiTullio, Carl Lamborg, John P. McCrow, Chris Dupont, Andrew E. Allen, Mak Saito

In this study we explored the metabolism of unicellular eukaryotic organisms (protists) across a 4,600 km meridional transect in the central Pacific Ocean. The region contains a natural biogeochemical gradient spanning from low nitrogen, oligotrophic waters to a productive equatorial upwelling system. We used a combined geochemical and 'omic approach to characterize the metabolic strategies these organisms rely upon to adapt to changes in their chemical environment.

Example input files are included here, except the master dataframes containing contigs, transcript raw reads/protein spectral counts, and annotation information which is available as supplemental material with the manuscript. Data availability, assembly, annotation and normalization methods are described in Cohen et al. 

## Input data files & generated tables

### Normalizations
* ```annotation_all.filtered.grps.go_TRANSCRIPTS_0.8lpi_Dino_only_sansAnnotations.csv``` - Raw transcript counts belonging to dinoflagellates (Lineage Probability Index (LPI) > 0.8)
* ```orflength.csv``` - open reading frame lengths from PhyloDB output, for normalizing to open reading length (ORF) length during transcripts per million (TPM) transformation 
* ```exclusive_counts_annotated_withscaffoldgroups_dinolpi0.8_NSAFpre.csv``` - Raw exclusive spectral counts belonging to dinoflagellates (LPI > 0.8)
* ```genus2family_dinos.csv``` - Manually curated list of dinoflagellate genera and their associated family - for relating genera to higher classification in the PhyloDB database
### 18S rRNA analysis with phyloseq
* ```OTU.csv``` - 18S V9 OTU table, OTUs were clustered using swarm (see Allen Lab rRNA pipeline here: https://github.com/allenlab/rRNA_pipeline)
* ```TAXA_newPR2.csv``` - OTU taxonomy table. OTUs were searched against Protist Ribosomal Reference (PR2) https://pr2-database.org/
* ```sampledata.csv``` - environmental metadata for ordination
### Heatmaps with pheatmap
* ```kodef.tab``` - KEGG KO to KEGG description key
* ```TPM_TRANSCRIPTS_Dino.lpi0.8_only_annotations_orf_KOpost.csv``` - normalized dinoflagellate transcript counts, retaining contigs with a lineage probability index (LPI) > 0.8. Annotations were subsetted to the KEGG (KO) annotation level
* ```TPM_transcripts_Dino.lpi0.8_KOdef_top100m.csv``` - normalized dinoflagellate transcript counts, retaining contigs with a lineage probability index (LPI) > 0.8. Annotations were subsetted to the KEGG (KO) annotation level, retaining samples from > 100m depth.
* ```annotation_all.filtered.grps.go.lpi_0.8_dino_diatom_hapto_TPM_KOpre.csv``` - Normalized transcript counts from dinoflagellates, diatoms, and haptophytes, using LPI > 0.8. Annotated at the KEGG functional level
* ```exclusive_counts_annotations_dino_lpi0.8_post_NSAF_orf_Pfampost_annotated.csv``` - NSAF-normalized exclusive spectral counts, summed to the PFAM annoation level, using contigs with LPI > 0.8
### Ordinations in vegan
* ```TPM_TRANSCRIPTS_Dino.lpi0.8_only_annotations_orf_allcontigs.csv``` - normalized dinoflagellate transcript counts, retaining contigs with a lineage probability index (LPI) > 0.8. 
* ```CCA_meta.csv``` - environmental metadata from sites/depths for ordination
### Weighted Correlation Network Analysis (WGCNA) and KEGG enrichment with clusterProfiler
* ```TPM_WGCNA_080519_orf_allsamples.csv``` - WGCNA module designations for each KEGG KO (including select PFAM IDs)
* CCA_meta_WCGNA_short.csv - environmental metadata from sites/depths
* ```bg.csv``` - KEGG KOs from all WGCNA modules
* ```blue.csv``` - KEGG KOs in the blue (deep: >200m) module (note: color was changed to black for clarity)
* ```turquoise.csv``` - KEGG KOs in the turquoise (shallow: <200m) module (note: color was changed to white for clarity)
* ```blue_generatios_080519_bars.csv``` - Gene ratio percentages in the blue module compared to total KOs identified, from KEGG enrichment 
* ```turquoise_generatios_080519_bars.csv``` - Gene ratio percentages in the turquoise module compared to total KOs identified, from KEGG enrichment 

## R code
* ```TPMnormalizations_heatmaps.R``` - code for NSAF and TPM normalizations and heatmaps
* ```phyloseq.R``` - OTU transformation, ordination, and plotting with ggplot2
* ```Dinoflagellate_KOG_transcripts_CA_ordination.R``` - code for generating ordinations
* ```WGCNA.R``` - dendrogram generation, network construction, traits to modules relationships, module expression and kegg enrichment analysis. Code was modified from Sarah Davies at Boston University.


## R packages used
* ggplot2
* vegan
* phyloseq
* WGCNA
* flashClust
* clusterProfiler
* ape
* caroline
* tibble
* dplyr
* pheatmap
* genefilter
* viridis
* RColorBrewer
* stringr


### Comments/questions/feedback welcome! ncohen@whoi.edu
