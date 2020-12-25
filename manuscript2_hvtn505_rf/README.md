# manuscript2_hvtn505_rf

This document aims to instruct how to reproduce results in the manuscript "Improving Random Forest Predictions in Small Datasets from Two-phase Sampling Designs". The source codes can be generally split into two main R files "helper_rf_hvtn.R" and "ExpeRF505". The detailed explanations are as follows.

## helper_rf_hvtn.R
This R file contains several helper functions that help to easly get results.

### Function list
- *screen_lasso*: implements lasso variable screening.
- *screen.dat.index*: generates a dataset and a screened variable index. Users can specify a candidate marker set and whether lasso variable screening is applied.
- *get_nms_group_all_antigens*: helps to specifify variables depending on assay type.
- *get.rf.cvauc*: fits different types of random forests(standard random forest (RF), RF with under-sampling (RF_under), RF with over-sampling (RF_over), and tuned RF (tRF)). Users can specify inver sampling probability weights (ipw) in model training.
- *get.glm.cvauc*: fits generalized linear models. Users can specify inver sampling probability weights (ipw) in model training.
- *get.st.cvauc*: fits stacking. 

## ExpeRF505.R
This R file conducts all experiments in the manuscript.

Step 1) Import HVTN 505 dataset
- i) First download the HVTN 505 dataset at https://atlas.scharp.org/cpas/project/HVTN%20Public%20Data/HVTN%20505/begin.view
- ii) install R package HVTN505 using the code "devtools::install_local("HVTN505_2019-4-25.tar.gz")".

Step 2) Experiments
- Codes in "2-1. Random forest (RF)" generate all random forest-based results in Table 1, 2, and 3 in the manuscript. You can specify whether clinical covariates are included or whether lasso variable screening is applied.
- Codes in "2-2. Generalized linear models (GLM)" generates all generalized linear models-based results in Table 3 and supplemental table. You can specify whether clinical covariates are included or whether lasso variable screening is applied.
- Codes in "2-3. Stacking (RF + GLM)" generate all stacking results in Table 4. To fit stacking, you should first import five R files that are "helper_functions.R", "caretList.R", "caretEnsemble.R", "caretStack.R", and "method.R". This five files come from the *caretEnsemble* and the *SuperLearner* R packages (Deane-Mayer and Knowles, 2016), (Polley et al. 2019).
