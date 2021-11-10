library(readxl)
library(switchBox)
library(pracma)
set.seed(123)

# read excel file which has binarized z-scores of metabolites and training histone data
# mention the appropriate sheet name (for reference, I am using threshold 0)
HM_Data <- read_excel("CCLE_metab_prediction_from_histone.xlsx", sheet = "Hist_metabZ0")
HM_Data_LeRoy <- read_excel("LeRoy_metab_prediction_from_histone.xlsx", sheet = "Hist_metabZ0")
metab_names <- colnames(HM_Data[1:136])

# replacing metabolite names with common name - "Target"
colnames(HM_Data)[1:136] <- "Target"
colnames(HM_Data_LeRoy)[1:136] <- "Target"

# Histone Markers which are used as training dataset are present in excel sheet from columns 138:188
# X contains training data
X <- HM_Data[138:188]

# create vectors to store final results
final_acc <- c()
final_mcc <- c()
final_prec <- c()
final_rec <- c()
final_f1 <- c()
final_guess_acc <- c()

# The classification is repeated for 136 metabolites which are present in columns 1:136 in excel file.
# (exclude the columns which have all 0 or 1's in them for Xtest, as they will not represent binary classification.)
for (metab in 1:136)
{
  # Y contains predictor data  
  Y <- HM_Data[metab]
  # Training is done on the entire CCLE dataset
  Xtrain <- X
  Ytrain <- Y
  
  # Testing is done on the LeRoy data
  Xtest <- HM_Data_LeRoy[138:188]
  Ytest <- HM_Data_LeRoy[metab]
  
  # Taking transpose to bring data to correct format
  Xtrain <- t(Xtrain)
  Xtest <- t(Xtest)
  
  df_Ytrain <- data.frame(Ytrain)
  df_Ytrain$Target <- as.factor(df_Ytrain$Target)
  Ytrain_fac <- df_Ytrain$Target
  
  df_Ytest <- data.frame(Ytest)
  df_Ytest$Target <- as.factor(df_Ytest$Target)
  Ytest_fac <- df_Ytest$Target
  
  # Train a classifier using default filtering function based on the Wilcoxon test
  classifier <- SWAP.Train.KTSP(Xtrain, Ytrain_fac, krange=c(3:15))
  
  testPrediction <- SWAP.KTSP.Classify(Xtest, classifier)
  
  # Resubstitution performance in the TEST set
  tab <- table(testPrediction, Ytest_fac)
  acc <- sum(Ytest_fac == testPrediction)/length(Ytest_fac)
  prec <- sum(Ytest_fac == 1 & testPrediction == 1)/sum(testPrediction == 1)
  rec <- sum(Ytest_fac == 1 & testPrediction == 1)/sum(Ytest_fac == 1)
  F1 <- 2*rec*prec/(rec+prec)
  
  # MCC Calculation
  TP <- tab[1]
  TN <- tab[4]
  FP <- tab[3]
  FN <- tab[2]
  N <- TN + TP + FN + FP
  S <- (TP + FN)/N
  P <- (TP + FP)/N
  MCC <- ((TP/N) - (S*P))/(sqrt(P*S*(1-S)*(1-P)))
  
  # Random guess accuracy
  guess_acc <- rep(0,100)
  for (j in 1:100)
  {
    guess <- Ytest_fac[randperm(length(Ytest_fac))]
    guess_acc[j] =sum(Ytest_fac == guess)/length(Ytest_fac)
  }
  
  # Final performance parameter values for each metabolite
  final_acc[metab] <- acc
  final_mcc[metab] <- MCC
  final_prec[metab] <- prec
  final_rec[metab] <- rec
  final_f1[metab] <- F1
  final_guess_acc[metab] <- guess_acc
  
  cat("Finished processing metabolite - ", metab, "\n")
}

# Convert performance parameters into data frames
final_acc_df <- data.frame(final_acc)
final_mcc_df <- data.frame(final_mcc)
final_prec_df <- data.frame(final_prec)
final_rec_df <- data.frame(final_rec)
final_f1_df <- data.frame(final_f1)
final_guess_acc_df <- data.frame(final_guess_acc)

# Assign original metabolite names to result files
rownames(final_acc_df) <- metab_names
rownames(final_mcc_df) <- metab_names
rownames(final_prec_df) <- metab_names
rownames(final_rec_df) <- metab_names
rownames(final_f1_df) <- metab_names
rownames(final_guess_acc_df) <- metab_names

# write the final performance parameter values for each metabolite in text file for storage.
# Path can be added to store it at specified location
write.table(final_acc_df, "acc_LeRoy.txt", sep="\t")
write.table(final_mcc_df, "mcc_LeRoy.txt", sep="\t")
write.table(final_prec_df, " prec_LeRoy.txt", sep="\t")
write.table(final_rec_df, "rec_LeRoy.txt", sep="\t")
write.table(final_f1_df, "f1_LeRoy.txt", sep="\t")
write.table(final_guess_acc_df, "guess_acc_LeRoy.txt", sep="\t")