# MSc_Project_2026_Cancer_Parameters
Collection of code used for my MSc Data Science for Biology dissertation project entitled: "Comparing Gene Activity Patterns in Pre and Post Immunotherapy Cells".

Different parts of the project are saved in different folders:

* **1_NB_Methods:** Tests on different methods of estimating negative binomial (NB) parameters from simulated count data in R to choose the most robust method.
* **2_Simulations:** Simulating count data and using the chosen method of estimation to determine errors on estimations
* **3_Real_Data:** Using everything learned during simulations to estimate parameters of real biological data. 
* **estimator_function:** key functions developed during the project for estimating that can easily be used on any gene count data. Use '?function_name' to read documentation.

There are html documents available for all markdown files so outputs can be viewed without downloading and running everything, as well as separate saved figures. 
While immunotherapy's effect on cancer cells is the focus of this study, the methods were designed to be reproducible on any gene count data.

The RNAseq data was obtained from this link:
[https://zenodo.org/records/14511579](https://zenodo.org/records/14511579)

 And the corresponding GitHub for this data source found here:  
[https://github.com/MahnoorNGondal/scRNA-seq-ICB-cohorts](https://github.com/MahnoorNGondal/scRNA-seq-ICB-cohorts)

