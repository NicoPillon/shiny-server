---
title: "Tutorial"
output: html_document
---

####  <a id="displayText" href="javascript:toggle(1);">I cannot find my favorite gene!</a>
  <div class="div_help" id="toggleText1" style="display: none">

Start typing the gene name and suggestions will appear in the scroll menu. MetaMEx works with official gene symbols, for instance the official gene name of PGC1α is PPARGC1A.

<img src="tutorial/tutorial_select_gene.svg" width="60%"/>

  </div>  

  
####  <a id="displayText" href="javascript:toggle(2);">My favorite gene is not detected in many studies!</a>
  <div class="div_help" id="toggleText2" style="display: none">

Older studies, or custom arrays often have a limited number of probes and therefore fewer detected genes. On the other hand, the more recent RNA sequencing datasets often have more depth and detect non-coding RNAs which are not present in gene arrays.

  </div>  
  
  
####  <a id="displayText" href="javascript:toggle(3);">How do I read a forest plot?</a>
  <div class="div_help" id="toggleText3" style="display: none">

A forest plot is a graphical representation of results from several scientific studies and is typically used to plot meta-analyses. The left-hand columns list the names of the studies, followed by the fold-change (log2), false discovery rate (FDR) and sample size (n) for each individual study. The right-hand column is a plot of the fold-change (log2) represented by a square and the 95% confidence intervals represented by horizontal lines. The area of each square is proportional to the study's weight (sample size) in the meta-analysis. The overall meta-analysed score is represented by a diamond on the bottom line, the lateral points of which indicate confidence intervals. 

<img src="tutorial/tutorial_forestplot.svg" width="80%"/>

  </div> 
  
  
####  <a id="displayText" href="javascript:toggle(4);">How do I select my population of interest?</a>
  <div class="div_help" id="toggleText4" style="display: none;">

MetaMEx compiles more than 90 studies which include volunteers of different age, sex, weight, fitness, weight and health. Studies can be included or excluded from the analysis by scrolling at the bottom of the page and checking the boxes. For instance, select males or females by checking the corresponding tick boxes.

<img src="tutorial/tutorial_select_population.svg" width="70%"/>

*	**Sex.** Choose whether you want males (M) or females (F). Some studies have pooled males and females or did not provide sex information and are labelled as undefined (U).
*	**Age.** Studies in MetaMEx are split into three age groups: young (<35), middle age (35-60) and elderly (>60).
*	**Fitness**. Activity levels were determined based on the description of the cohorts available in the publications. Sedentary is defined as no formal exercise training. Individuals performing exercise for more than 150 min per week and/or having an average VO2max are considered active. Athletes are individuals engaged in formal and regular exercise training and exhibit good to excellent VO2max.
*	**Weight.** Body composition is based on body mass index provided in the publications and the actual definition of lean (BMI<25), overweight (25≤BMI<30), obese (30≤BMI<40) and morbidly obese (BMI≥40).
*	**Muscle.** Most studies do cycling exercise and therefore collect vastus lateralis (quadriceps) biopsies. However, a handful of studies used soleus or biceps biopsies. Sometimes the muscle biopsy is unknow and is therefore annotated as N.A. 
*	**Health.** MetaMEx includes studies from healthy individuals with no history of disease as well as people diagnosed with metabolic diseases or other chronic conditions such as chronic kidney disease or frailty.

  </div>  
  
  
####  <a id="displayText" href="javascript:toggle(5);">What filters can I apply?</a>
  <div class="div_help" id="toggleText5" style="display: none">

After selecting  either acute exercise, exercise training or inactivity, a specific menu will appear on the right of the page. This menu includes parameters such as exercise duration or time of biopsy collection after exercise cessation. Another list will appear under the forest plot to select or unselect specific datasets.

<img src="tutorial/tutorial_filters.svg" width="80%"/>

* **Acute exercise studies.** For acute exercise protocols, it is possible to customize the type of exercise (concentric, eccentric or mixed) and the time of the biopsy collection after exercise cessation.
* **Exercise training studies.** For exercise training protocols, it is possibly to customize the duration of the training (from 1 week to lifelong) and the time of the biopsy collection after exercise cessation. It is also possible to include/exclude specific studies based on their GEO accession number.
* **Inactivity studies.** MetaMEx includes two inactivity protocols: bed rest or limb immobilization. It is also possible to customize the duration of the inactivity and include/exclude specific studies based on their GEO accession number.

  </div>  


####  <a id="displayText" href="javascript:toggle(7);">What analysis and statistical methods were used?</a>
  <div class="div_help" id="toggleText7" style="display: none">

The meta-analysis was created by collecting publicly available studies on mRNA expression levels in human skeletal muscle after exercise or inactivity. Statistics were first performed individually for each array. 

* Robust multiarray averaging was used for affymetrix arrays (oligo package)
* Quantile normalization was used for other microarrays (limma package)
* RNA sequencing datasets were processed using the edgeR package following the standard pipeline (https://bioconductor.org/packages/devel/workflows/vignettes/RnaSeqGeneEdgeRQL/inst/doc/edgeRQL.html). 

When a gene is selected, the corresponding data is collected and t-tests are performed comparing each group to its baseline. Using statistics from the t-tests, the meta-analysis summary is calculated using restricted maximum likelihood (metafor package). The analysis is weighted using sample size (n) to adjust for studies with small number of volunteers. 

The timeline is calculated by collecting all data available in the database. Moderated t-statistics are calculated with empirical Bayes statistics after blocking for other confounding parameters.

Since MetaMEx v3.2205, p values are adjusted for multiple-comparison with the more conservative Bonferroni method. This compensates for the increasing number of studies added to the database leading to more statistical power. The use of Bonferroni correction further reduces the risk of false positives.

  </div>


####  <a id="displayText" href="javascript:toggle(8);">Why are statistics on the website different from what is reported in the original publications?</a>
  <div class="div_help" id="toggleText8" style="display: none">

Whenever possible, we downloaded the raw data and re-processed studies from raw data files (CEL, fastq...). That means that the processing and normalization methods that we used might differ from the ones used by the original authors. In addition, samples were often insufficiently annotated to allow us to run paired statistics comparing pre/post interventions. We therefore have to use unpaired statistics and lose power in the process. Finally, many studies pooled individuals of different age and BMI to have higher sample size. To allow proper comparison in the meta-analysis, we split these studies into sub groups and analyzed them separately, therefore reducing the sample size and statistical power.

  </div>
  

####  <a id="displayText" href="javascript:toggle(9);">What is the version history of the database?</a>
  <div class="div_help" id="toggleText9" style="display: none">

* MetaMEx 4.2208 - September 2024. Major update. Use of machine learning to replace missing metadata (sex, age, BMI). The app now uses individual data instead of summaries by studies.  
* MetaMEx 3.2208 - September 2022. Minor update with recently published datasets. Fix of minor bugs.
* MetaMEx 3.2207 - July 2022. Addition of mouse datasets. Addition of timeline plots for acute and inactivity studies. All RNA sequencing datasets realigned to the most recent annotations (Human GRCh38.102 and mouse GRCm39.104). Addition of the most recently published human and mouse datasets.
* MetaMEx 2.2101 - Jan 2021. Addition of recently published studies.
* MetaMEx 2.2008 - Aug 2020. Major update of the database and style of the app. Addition of the most recent published datasets. Adjustment of colors and style for accessibility. 
* MetaMEx 1.1912 - Dec 19, 2019. Added published HIIT studies.
* MetaMEx 1.1902 - Feb 12, 2019. Updated clinical characteristics of studies with obesity status.
* MetaMEx 1.1809 - Sep 20, 2018. Added feature for correlations, improved speed and added progress bars.
* MetaMEx 1.1805 - Initial release on May 22, 2018.

  </div>
  

####  <a id="displayText" href="javascript:toggle(10);">Where any studies excluded from the database?</a>
  <div class="div_help" id="toggleText10" style="display: none">
  
The following studies were excluded from the database:

* GSE156247 was excluded because it contains the exact same data as GSE53598.
* GSE163434 was excluded because it is a re-analysis of GSE157585.
* GSE230002 and GSE137631 were not included because the intervention was a combination of exercise with weight loss (diet/surgery).
* GSE44818, GSE28998, GSE24235 were excluded because they did not include pre-exercise or baseline biopsies.
* GSE165630. This study compared trained athletes to untrained individuals, but the data clustered with acute exercise studies, suggesting that the biopsies were not taken long enough after the last exercise session.
* GSE1718. This study compared individuals with/without changes in insulin sensitivity following a training intervention. The only data available is already processed in a way that makes it impossible to separate conditions. In addition, the samples are a mix of biopsies coming from males and females.


In addition, a few subgroups not relevant for the database were removed in specific studies: Detraining, bed rest with exercise/nutrition intervention, reloading after inactivity, insulin infusion and functional electrical stimulation in spinal chord injury.


  </div>
