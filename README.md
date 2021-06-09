# test_cvauc
This document aims to instruct how to reproduce results in the manuscript "Testing a Global Null Hypothesis Using Ensemble Machine Learning Methods". The source codes are split into four main R files "helper_cvauc.R", "experiments_sec2.R", "experiments_sec3.R", and "experiments_real.R". The details are as follows.

## helper_cvauc.R
The R file contains several helper functions that help reproduce results easily.

### Function list
- *sim.1*: generates simulated dataset for Section 2.
- *sim.rv144*: generates simulated dataset for Section 3.
- *kfold.split*: performs k-fold cross validation.
- *ran.kfold.split*: performs random 4:1.
- *lpo.split*: performs Leave-Pair-Out cross validation
- *get.splits*: produces subject indexes for a variety of cross validation schemes.
- *get.cv.auc*: fits multiple logistic regression and estimates cvAUC.
- *screen_lasso*: performs lasso variable screening.
- *get.st.auc*: performs stacking minimum P value and random forest.

## experiments_sec2.R
The R file contains codes for simulation studies in Section 2.

Step 0) Importing helper function
- You should import "helper_cvauc.R" file.

Step 1) Setting arguments
- You can set experimental conditions using arguments, such as sim.setting, basic.cv.setting, and proj. The detailed comments for each argument are in the R file.
- To make Figure A.1 in the Supplementary Materials, you can set "make.plot" argument with TRUE.

Step 2) Experiments
- Codes in "2-1. Comparison CV schemes" are for Table 1 in the manuscript.
- Codes in "2-2. Comparison permutation test and bootstrap CI-based test" are for Table 2 in the manuscript.
- Codes in "2-3. Plotting permutation distribution and bootstrap distribution" are for Figure A.1 in the Supplementary Materials in the manuscript.

## experiments_sec3.R
The R file contains codes for simulation studies in Section 3.

Step 0) Importing helper function
- You should import "helper_cvauc.R" file.
- To perform stacking method, you should import four R files under the folder caretEnsembleR, where there are "helper_functions.R", "caretList.R", "caretEnsemble.R", and "caretStack.R".

Step 1) Setting arguments
- You can set experimental conditions using arguments, such as sim.setting, sim.linear, fit.setting, proj, ipw, and scr. The detailed comments for each argument are in the R file.

Step 2) Experiments
- Codes in "2-1. Fitting testing methods" are for fitting all tests, including perm_MP, perm_RR, perm_BG, perm_RF, perm_AB, and perm_stacking, in the manuscript. These codes are used for Table 3 and 4 in the manuscript.

Step 3) Figure C.1 in supplementary material
- This codes are to draw Figure C.1 in the Supplementary Materials.

## experiments_real.R
The R file contains codes for real data experiments in Section 4.

Step 0) Importing helper function
- You should import "helper_cvauc.R" file.
- To perform stacking method, you should import four R files under the folder caretEnsembleR, where there are "helper_functions.R", "caretList.R", "caretEnsemble.R", and "caretStack.R".

Step 1) Real RV144 dataset
- To access RV144 dataset, you can install R package RV144cc using the code "devtools::install_local("RV144cc_2019.08-14.tar")".
- You can generate the full dataset or the driven data, which are all markers without IgAC1.

Step 2) Setting arguments
- You can set experimental conditions using arguments, such as setting, fit.setting, and proj. The detailed comments for each argument are in the R file.

Step 3) Experiments
- Codes in "2-1. Fitting testing methods" are for fitting all tests, including perm_MP, perm_RR, perm_BG, perm_RF, perm_AB, and perm_stacking, in the manuscript. These codes are used for Table 5 in the manuscript.
