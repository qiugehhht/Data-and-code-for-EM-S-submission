# Data and code for EM&S submission
Please decompress "./dataset/Stack_Shanghai.gdb.zip" and sequentially run "prepossesing.R", "stacked_model.R" and "inla_prediction.R" to reproduce the findings in the manuscript.

The script "prepossesing.R" performs data preprocessing and variable selection.

The script "stacked_model.R" constructs base models - including RF, SVM, and ERT - and the stacked ensemble model using OLS as the meta-model.

The script "inla_prediction.R" constructs the spatial Bayesian stacked ensemble model.
