ScalingMicroBiodiversity
=======

Code, data, and text documents related to the study of abundance scaling relationships focused on taxonomic aspects of microbial biodiversity (evenness, dominance, richness, rarity), using datasets of plants and animals for comparative perspective and signals of potentially unifying patterns.

##Recreating analyses of Locey and Lennon (2016)

Download or clone the public GitHub repository https://github.com/klocey/DiversityTools into the home directory and then running the Cython script:https://github.com/klocey/DiversityTools/blob/master/radMetrics/radMetrics.pyx
###Quantifying dominance, evenness, rarity, and richness
For each site in our data compilation, all aspects of diversity (dominance, evenness, rarity, richness) were calculated or estimated from the species abundance distribution (SAD). These and all other analyses can be reproduced by following the directions given here: https://github.com/klocey/ScalingMicroBiodiversity/blob/master/README.md. Which primarily requires the user to download the combined code/data repository and run the code for the respective figure or analysis.
###Global-scale analysis of dominance-abundance (Nmax vs. N) relationshipThe Python script to recreate this analysis and figure is available here, and can be modified to use open or closed reference data as well as to exclude or include singletons: https://github.com/LennonLab/ScalingMicroBiodiversity/blob/master/fig-scripts/Fig2/Fig2.py


### Supplemental figures & analyses:

*1.) Richness* 

2.) Evenness

3.) Dominance

4.) Rarity

5.) Effect of sample size

6.) Effect of the microbe/macrobe categorical variable

7-15.) Per dataset results

Table 1.) Comparison of fits to power-law, linear, semi-log, and exponential functions