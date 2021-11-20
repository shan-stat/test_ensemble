# test_ensemble
This document aims to instruct how to reproduce results in the manuscript "Testing a Global Null Hypothesis Using Ensemble Machine Learning Methods". The source codes are split into four R files "helper_cvauc.R", "experiments_simul_1.R", "experiments_simul_2.R", and "experiments_real.R". The details are as follows.

## helper_cvauc.R
This file contains several helper functions that help reproduce results easily.

### Function list
- *sim.1*: generates simulated dataset for Section 2.1.
- *sim.rv144*: generates simulated dataset for Section 2.2 and 2.3.
- *sim.rv144.lognorm*: generates simulated dataset for Supplementary Materials Section D
- *get.cv.auc*: fits multiple logistic regression and estimates cvAUC.
- *screen_lasso*: performs lasso variable screening.
- *get.st.auc*: performs stacking Westfall and Young’s method and random forest.
- *get.cv.auc.LeDell*: permform LeDell CI-based test
- *get.cvtmle*: permform Benkeser CI-based test

## experiments_simul_1.R
Theis file contains codes for simulation studies in Section 2.1 and Supplementary Materials Section C.

Step 0) Importing helper function
- You should import "helper_cvauc.R" file.

Step 1) Setting arguments
- You can set experimental conditions using arguments, such as sim.setting and proj. The detailed comments for each argument are in the R file.

Step 2) Experiments
- Codes in "Section 2.1 (Choice of CV schemes)" are for Table 1 of the manuscript and Table B.1 of the Supplementary Materials.
- Codes in "Supplementary Materials Section C (Size of confidence interval-based hypothesis testing approaches)" are for Table C.1 of the Supplementary Materials.

## experiments_simul_2.R
Theis file contains codes for simulation studies in Section 2.2, 2.3, and Supplementary Materials Section D.

Step 0) Importing helper function
- You should import "helper_cvauc.R" file.
- To perform stacking method, you should import four R files under the folder caretEnsembleR, where there are "helper_functions.R", "caretList.R", "caretEnsemble.R", and "caretStack.R".

Step 1) Setting arguments
- You can set experimental conditions using arguments, such as sim.setting, sim.linear, fit.setting, proj, ipw, and scr. The detailed comments for each argument are in the R file.

Step 2) Experiments
- Codes in "2-1. Fitting testing methods" are for fitting all testing methods, including WY, perm_RR, perm_BG, perm_RF, perm_AB, and perm_stacking, in the manuscript and the Supplementary Materials. These codes are used for Table 2 and 3 in the manuscript, and Table D.1 in the Supplementary Materials.

Step 3) Figure A.2 in the Supplementary Materials
- This codes are to draw Figure A.2 in the Supplementary Materials.

## experiments_real.R
The file contains codes for real data experiments in Section 3.

Step 0) Importing helper function
- You should import "helper_cvauc.R" file.
- To perform stacking method, you should import four R files under the folder caretEnsembleR, where there are "helper_functions.R", "caretList.R", "caretEnsemble.R", and "caretStack.R".

Step 1) Real RV144 dataset
- While the RV144 dataset cannot be shared because they are proprietary to the US Military HIV Research Program, interested researchers may request access to anonymized patient-level data at PubRequest@hivresearch.org.

Step 2) Setting arguments
- You can set experimental conditions using arguments, such as setting, fit.setting, and proj. The detailed comments for each argument are in the R file.

Step 3) Experiments
- Codes in "2-1. Fitting testing methods" are for fitting all tests, including WY, perm_RR, perm_BG, perm_RF, perm_AB, and perm_stacking, in the manuscript. These codes are used for Table 4 in the manuscript.
