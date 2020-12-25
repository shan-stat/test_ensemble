# Improving Random Forest Predictions in Small Datasets from Two-phase Sampling Designs

This document aims to instruct how to reproduce results in the manuscript "Improving Random Forest Predictions in Small Datasets from Two-phase Sampling Designs". The source codes can be generally split into two main R files "helper_rf_hvtn.R" and "ExpeRF505". The detailed explanations are as follows.

## helper_rf_hvtn.R
This R file contains several helper functions that help to easly get results.

### Function list
- *screen_lasso*: implements lasso variable screening.
- *screen.dat.index*: generates a dataset and a screened variable index. Users can choose a specific marker set and lasso variable screening.
- *get_nms_group_all_antigens*: helps to specifify variables depending on assay type.
- *get.rf.cvauc*: fits different types of random forests(standard random forest (RF), RF with under-sampling (RF_under), RF with over-sampling (RF_over), and tuned RF (tRF)). Users can specify inver sampling probability weights (ipw) in model training.
- *get.glm.cvauc*: fits generalized linear models. Users can specify inver sampling probability weights (ipw) in model training.
- *get.st.cvauc*: fits stacking. 

