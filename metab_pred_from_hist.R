library(readxl)
library(switchBox)
library(pracma)
set.seed(123)

# read excel file which has binarized z-scores of metabolites and training histone data
# mention the appropriate sheet name (for reference, I am using threshold 0)
HM_Data <- read_excel("CCLE_metab_prediction_from_histone.xlsx", sheet = "Hist_metabZ0")

# Histone Markers which are used as training dataset are present in excel sheet from columns 138:188
# X contains training data
X <- HM_Data[138:188]
metab_names <- colnames(HM_Data[1:136])

# replacing metabolite names with common name - "Target"
colnames(HM_Data)[1:136] <- "Target"

# randomly divide the data into 10 sets for 10-fold cross validation
indices <- sample(10,nrow(HM_Data), replace = TRUE, prob = c(0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1))

# create vectors to store final results
final_acc <- c()
final_mcc <- c()
final_prec <- c()
final_rec <- c()
final_f1 <- c()
final_guess_acc <- c()
p_val <- c()

# The classification is repeated for 136 metabolites which are present in columns 1:136 in excel file.
for (metab in 1:136)
{
  # Y contains predictor data  
  Y <- HM_Data[metab]
  df <- data.frame(Y)
  Y_fac <- as.factor(df$Target)
  
  # vectors to store the 10 performance paramters values generated from cross validation for every metabolite
  acc <- rep(10,1)
  prec <- rep(10,1)
  rec <- rep(10,1)
  F1 <- rep(10,1)
  MCC <- rep(10,1)
  
  for (k in 1:10)  #10-fold Cross Validation
  {
    test <- indices == k 
    train <- indices != k
    Xtrain <- X[train,]
    Ytrain <- Y[train,]
    Xtest <- X[test,]
    Ytest <- Y[test,]
    
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
    acc[k] <- sum(Ytest_fac == testPrediction)/length(Ytest_fac)
    prec[k] <- sum(Ytest_fac == 1 & testPrediction == 1)/sum(testPrediction == 1)
    rec[k] <- sum(Ytest_fac == 1 & testPrediction == 1)/sum(Ytest_fac == 1)
    F1[k] <- 2*rec[k]*prec[k]/(rec[k]+prec[k])
    
    # MCC Calculation
    TP <- tab[1]
    TN <- tab[4]
    FP <- tab[3]
    FN <- tab[2]
    N <- TN + TP + FN + FP
    S <- (TP + FN)/N
    P <- (TP + FP)/N
    MCC[k] <- ((TP/N) - (S*P))/(sqrt(P*S*(1-S)*(1-P)))
  }
  
  # Random guess accuracy
  guess_acc <- rep(0,100)
  for (j in 1:100)
  {
    guess <- Y_fac[randperm(length(Y_fac))]
    guess_acc[j] =sum(Y_fac == guess)/length(Y_fac)
  }
  
  # Final performance parameter values for each metabolite
  final_acc[metab] <- mean(acc, na.rm = TRUE)
  final_mcc[metab] <- mean(MCC, na.rm = TRUE)
  final_prec[metab] <- mean(prec, na.rm = TRUE)
  final_rec[metab] <- mean(rec, na.rm = TRUE)
  final_f1[metab] <- mean(F1, na.rm = TRUE)
  final_guess_acc[metab] <- mean(guess_acc, na.rm = TRUE)
  p_val[metab] <- t.test(acc, guess_acc)$p.value
  
  cat("Finished processing metabolite - ", metab, "\n")
  
}

# Convert performance parameters into data frames
final_acc_df <- data.frame(final_acc)
final_mcc_df <- data.frame(final_mcc)
final_prec_df <- data.frame(final_prec)
final_rec_df <- data.frame(final_rec)
final_f1_df <- data.frame(final_f1)
final_guess_acc_df <- data.frame(final_guess_acc)
final_pval_df <- data.frame(p_val)

# Assign original metabolite names to result files
rownames(final_acc_df) <- metab_names
rownames(final_mcc_df) <- metab_names
rownames(final_prec_df) <- metab_names
rownames(final_rec_df) <- metab_names
rownames(final_f1_df) <- metab_names
rownames(final_guess_acc_df) <- metab_names
rownames(final_pval_df) <- metab_names

# write the final performance parameter values for each metabolite in text file for storage.
# Path can be added to store it at specified location
write.table(final_acc_df, "acc.txt", sep="\t")
write.table(final_mcc_df, "mcc.txt", sep="\t")
write.table(final_prec_df, " prec.txt", sep="\t")
write.table(final_rec_df, "rec.txt", sep="\t")
write.table(final_f1_df, "f1.txt", sep="\t")
write.table(final_guess_acc_df, "guess_acc.txt", sep="\t")
write.table(final_pval_df, "pval.txt", sep="\t")