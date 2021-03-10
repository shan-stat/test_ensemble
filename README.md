# test_cvauc
This document aims to instruct how to reproduce results in the manuscript "Testing a Global Null Hypothesis Using Cross Validated Area Under the ROC Curve". The source codes are split into four main R files "helper_cvauc.R", "experiments_sec2.R", "experiments_sec3.R", and "experiments_real.R". The details are as follows.

## helper_cvauc.R
The R file contains several helper functions that help reproduce results easily.

### Function list
- *sim.1*: generates simulated dataset for Section 2.
- *sim.rv144*: generates rv144 simulated dataset for Section 3.
- *kfold.split*: performs k-fold cross validation.
- *ran.kfold.split*: performs random 4:1 cross validation.
- *lpo.split*: performs leave-pair-out cross validation
- *get.splits*: produces subject indexes for a variety of cross validation schemes.
- *get.cv.auc*: fits multiple logistic regression and estimates cv-auc.
- *get.cv.auc.LeDell*: fits LeDell's confidence interval.

## experiments_sec2.R
The R file contains codes for experiments in Section 2.

Step 1) Setting process arguments
- You can set conditions of the experiment using the arguments of sim.setting, basic.cv.setting, and proj. The detailed description is in R file.
- To make Figure 1 in the manuscript, you can set the make.plot argument with TRUE.

Step 2) Experiments
- Codes in "2-1. Comparison of CV schemes" are used for results in Table 2 in the manuscript.
- Codes in "2-2. Fitting tests" are used for results in Table 3 in the manuscript.
- Codes in "2-3. Figure" are used for Figure 2 in the manuscript.

## experiments_sec3.R
The R file contains codes for experiments in Section 3.

Step 1) Setting process arguments
- You can set conditions of the experiment using the arguments of sim.setting, fit.setting, proj, and ipw. The detailed description is in R file.

Step 2) Experiments
- Codes in "2-1. Fitting tests" are used for fitting all of tests listed in Section 2 in the manuscript. These codes are also used for results in Table 4 in the manuscript.

Step 3) Figure B.1 in supplementary material
- This codes are to draw heatmap plot in Figure B.1 in supplementary material.

## experiments_real.R
The R file contains codes for experiments in Section 4.

...
