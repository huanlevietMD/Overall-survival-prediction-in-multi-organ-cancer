# Overall-survival-prediction-in-multi-organ-cancer
Development and validation of CT-based radiomics signature for overall survival prediction in multi-organ cancer
# 2. Code presentation
# 2.1 Patient cohort
Four TCIA public data sets of cancer patients in three distinct organs were enrolled in the retrospective analysis (420, 157, 137, and 191 patients for Lung 1 training, Lung 2 testing, and two external validation set: Kidney and Head & Neck respectively)
# 2.2 Radiomics feature extraction
The data set's author segmented the tumors in the datasets containing CTscan images from several institutions using various image preprocessing techniques. The 3D Slicer program was used in this study to extract radiomics features from CT scans. The Pyradiomics package in Python was used to extract radiomics features that met the IBSI requirements. For each CT scan from one patient, we received 851 radiomics features, which we divided into four categories: tumor intensity, shape, texture, and wavelet filters. Customizing Bin Width is set at 25.
# 2.3.	Feature selection and best model construction in Training set (Lung 1 data set)
# 2.3.1 Removing Highly Correlated radiomics Features in the Training Set
Pairwise correlations were used in the training set to get rid of redundant radiomics features and prevent over-fitting or bias in the analysis. We employed Pearson's correlation analysis to measure the correlation of every radiomics features pairs in the training set. If the correlation coefficient between two features is greater than 0.7, which illustrated a high correlation, we determined to exclude one with lower correlation coefficient to the target variable (survival outcome) and retain the more significant one
# 2.3.2 LASSO regression model
To identify the most important 15 features for predicting survival, the LASSO regression model is utilized with the glmnet package in the R language, all of the radiomics features were fed into the LASSO regression model. We applied a 10-fold cross-validation in the training set to find optimal lamda to avoiding model simplification and overfitting.
# 2.3.3 Univariable Cox proportional hazard models were constructed to find the 10 radiomics signatures
From the 15 radiomics features identified by LASSO in the training set, we developed and compared univariate Cox proportional hazard models to identify the 10 radiomics signatures associated with overall survival time in the training set using the survival package of the R language (statistical significance was determined when the p-value was less than 0.05).
# 2.3.4 Risk score calculation by multivariable Cox proportional hazard model
Calculating risk score for each patient using the 10 radiomics signatures by the Est.PH function in the “survC1” package
# 2.3.5 Kaplan-Meier plot
Following the median risk score, the Kaplan-Meier curve was used to show the categorization of patients as high- or low-risk groups.
# 2.3.6 Assess the radiomics signature and clinical parameters integration's effectiveness
The iAUC was obtained for each predictive model using the "risksetROC" and “SurvC1” package. Bootstrapped resampling with 1,000 repetitions was used to calculate the differences in iAUC between the multivariable predictive models (radiomics model, clinical model and combine model), if the 95% CI of the iAUC difference did not contain a zero value, the difference was judged statistically significant. Risk prediction abilities were assessed by graphing the Brier score prediction error curves across survival times by the “pec” package in R
# 2.4 Finally, in the testing and 2 validation sets, we locked and independently assessed the risk score created by the best radiomics signatures model for overall survival prediction.
Based on the 10 radiomics signatures that we found in the Training set we locked and chose the similar 10 radiomics signature in Testing set (Lung 2) and two validation set (Head & Neck validation set and Kidney validation set). We perform the same statistical analysis steps to find the corresponding risk score for each data set. From the risk score scale, we draw the Kaplan Meier curve and calculate and compare the effectiveness of each model in survival prognosis.
