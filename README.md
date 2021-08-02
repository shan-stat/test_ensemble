# test_ensemble
This document aims to instruct how to reproduce results in the manuscript "Testing a Global Null Hypothesis Using Ensemble Machine Learning Methods". The source codes are split into three main R files "helper_cvauc.R", "experiments_simul.R", and "experiments_real.R". The details are as follows.

## helper_cvauc.R
The R file contains several helper functions that help reproduce results easily.

### Function list
- *sim.1*: generates simulated dataset for Section 2.1.
- *sim.rv144*: generates simulated dataset for Section 2.2 and 2.3.
- *kfold.split*: performs k-fold cross validation.
- *ran.kfold.split*: performs random 4:1.
- *lpo.split*: performs Leave-Pair-Out cross validation
- *get.splits*: produces subject indexes for a variety of cross validation schemes.
- *get.cv.auc*: fits multiple logistic regression and estimates cvAUC.
- *screen_lasso*: performs lasso variable screening.
- *get.st.auc*: performs stacking Westfall and Young’s method and random forest.

## experiments_simul.R
The R file contains codes for simulation studies in Section 2.

### 1) Section 2.1
Step 0) Importing helper function
- You should import "helper_cvauc.R" file.

Step 1) Setting arguments
- You can set experimental conditions using arguments, such as sim.setting and proj. The detailed comments for each argument are in the R file.

Step 2) Experiments
- Codes in "2-1. Comparison CV schemes" are for Table 1 in the manuscript.

### 2) Section 2.2 and 2.3

Step 0) Importing helper function
- You should import "helper_cvauc.R" file.
- To perform stacking method, you should import four R files under the folder caretEnsembleR, where there are "helper_functions.R", "caretList.R", "caretEnsemble.R", and "caretStack.R".

Step 1) Setting arguments
- You can set experimental conditions using arguments, such as sim.setting, sim.linear, fit.setting, proj, ipw, and scr. The detailed comments for each argument are in the R file.

Step 2) Experiments
- Codes in "2-1. Fitting testing methods" are for fitting all testing methods, including WY, perm_RR, perm_BG, perm_RF, perm_AB, and perm_stacking, in the manuscript. These codes are used for Table 2 and 3 in the manuscript.

Step 3) Figure B.1 in the Supporting Information
- This codes are to draw Figure B.1 in the Supporting Information.

## experiments_real.R
The R file contains codes for real data experiments in Section 3.

Step 0) Importing helper function
- You should import "helper_cvauc.R" file.
- To perform stacking method, you should import four R files under the folder caretEnsembleR, where there are "helper_functions.R", "caretList.R", "caretEnsemble.R", and "caretStack.R".

Step 1) Real RV144 dataset
- While the RV144 dataset cannot be shared because they are proprietary to the US Military HIV Research Program, interested researchers may request access to anonymized patient-level data at PubRequest@hivresearch.org.

Step 2) Setting arguments
- You can set experimental conditions using arguments, such as setting, fit.setting, and proj. The detailed comments for each argument are in the R file.

Step 3) Experiments
- Codes in "2-1. Fitting testing methods" are for fitting all tests, including WY, perm_RR, perm_BG, perm_RF, perm_AB, and perm_stacking, in the manuscript. These codes are used for Table 4 in the manuscript.
