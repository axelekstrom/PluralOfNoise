# PluralOfNoise

**Introduction**

This repository contains three scripts that simulate the impact of random data on statistical significance in a given analysis pipeline, derived from Berthet et al. (2025). Both scripts repeat the core analysis (up to and including the generation of Figure 4) 1,000 times, each time replacing the FoCs variables with randomly generated data.

**Scripts**

_1 Random_

This script replaces all FoCs values with random "yes/no" data drawn from a binomial distribution with a 50% probability of "yes". Missing values (NA:s) in the original dataset are not preserved. The analysis is executed 1,000 times, and results from each run are stored. Finally, the script calculates the proportion of iterations that yield statistically significant results.

_2 Random with NA:s retained_

This version uses the same randomization method as Random', drawing from a binomial distribution with a 50% "yes" probability. However, it retains the original dataset's NA:s from Berthet et al., replacing only the non-missing values with randomized entries. As with the first script, all results are saved and the proportion of significant outcomes is computed.

_3 MCA_factors.Rmd_

This script 

(1) runs the MCA described in the paper with 3, 4, and 5 factors and compares the effects of the pairwise contrasts reported in Figure 4 in Berthet et al., 

(2) plots the number of occurrences of each compositional call identified in the paper by individual and group, 

(3) extracts the FOCs associated with each factor in the MCA, and for each FOC calculates the proportion of "no"'s and NAs in the data, and 

(4) performs the likelihood ratio tests described in the section "(iii) The meaning of AB is derived from the meaning of A and B" (Berthet et al., 2025) while directly testing the significance of the interaction effect.


**Contributors**

Andreas Wartel (Centre for Cultural Evolution, Department of Psychology, Stockholm University) wrote Scripts 1-2.

Johan L. Kleberg (Department of Psychology, Stockholm University) wrote Script 3. 

**References:**

Berthet, M., Surbeck, M., & Townsend, S. (2025). Science, 388(6742), 104-108. doi: 10.1126/science.adv1170
