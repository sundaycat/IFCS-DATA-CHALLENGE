# IFCS-Data-Challenge

A data challenge competition held by International Federation of Classification Societies(IFCS). Aiming to find a (semi)automatic classification of the lower back pain patients based on the baseline variables in order to find clinically applicable and useful groups.

My teammate and I are submitting our cluster analysis contribution in response to the data challenge. You'll find the following files:

1) "LBP_Project_Report.pdf" — this is our written report which describes the entire process and methods used in our analysis.
2) "cluster membership.xls" — this Excel file contains the membership based on our clustering result.
3) "consolFullOrdinal.Rdata" — this is our imputed dataset which we use to perform cluster analysis.

4) "LBP_ProjectCode.R" — this is our code for the data imputation, clustering, and other analysis.  Running this code will generate 4 csv files saved to current working directory (two files for variable contributions, one file for the selected important variables table, and one file for cluster membership).

5) "Alg_Comparison.R" — we used this file to compare different clustering algorithms only and don’t need to run it since the final selected algorithm is already used in the "LBP_ProjectCode.R".

The Excel files "variables_data_challenge.xlsx" and "data_challenge.xlsx" are the raw data that was given by IFCS. Please put them in the current working directory in R, along with “consolFullOrdinal.Rdata” file, when the project code file is run.
