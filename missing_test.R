
library("dplyr")
library("mice")
library("mi")
library("ggplot2")
# library("ggthemes")
library("caret")


###### GETTING THE DATASET
setwd("C://Users//matteo.fiorani//Desktop//pjs//heart-disease")
df <- read.csv("data//heart.csv")


###### FUNCTIONS
# ----------------------------------------------------------------
# -------------------- F U N C T I O N S -------------------------
# ----------------------------------------------------------------

generate_data <- function(model_name, df, perc, rr, target_var, positive){
  cat("generating data", fill = T)
  index <- which(names(df) == eval(target_var))
  pattern <- rep(1, length(names(df)))
  pattern[index] <- 0
  amp <- ampute(df, prop = perc, mech = "MAR", patterns = pattern)
  miss <- ifelse(is.na(amp$amp[, eval(target_var)]), 1, 0)

  # MULTIPLE IMPUTATIONS WITH MI PACKAGE
  if(model_name == "mi"){
    df_covariates_mi <- df
    df_covariates_mi[, eval(target_var)][miss == 1] <- NA

    mdf <- mi::missing_data.frame(df_covariates_mi)
    options(mc.cores = 2)
    cat(model_name, " imputations", fill = T)
    imputations <- mi::mi(mdf, n.iter = 5, n.chains = 4, max.minutes = 20)

    dfs <- mi::complete(imputations)
    cat(model_name, " imputations completed", fill = T)

    mi <- dfs$`chain:1`
    mi <- mi[mi[, ncol(mi)] == T, index]
  }
   
  # MULTIPLE IMPUTATIONS WITH MICE PACKAGE
  if(model_name == "mice"){
    df_covariates_mi <- df
    df_covariates_mi[, eval(target_var)][miss == 1] <- NA

    cat(model_name, " imputations", fill = T)
    imp <- mice(df_covariates_mi, method = "logreg", m = 1, maxit = 1)
    imputations <- mice::complete(imp)
    cat(model_name, " imputations completed", fill = T)
    
    mi <- imputations[miss == 1, index]
  }

  # MULTIPLE IMPUTATIONS WITH MICE PACKAGE USING RANDOMFOREST
  if(model_name == "missForest"){
    df_covariates_mi <- df
    df_covariates_mi[, eval(target_var)][miss == 1] <- NA
    
    cat(model_name, " imputations", fill = T)
    
    imp <- missForest::missForest(df_covariates_mi, verbose = FALSE)
    imputations <- imp$ximp
    cat(model_name, " imputations completed", fill = T)
    
    mi <- imputations[miss == 1, index]
  }
  
  label <- df[miss == 1, index]
  try(cm_mi <- confusionMatrix(as.factor(mi), as.factor(label), positive = positive))
  
  cat(model_name, "saving data", fill = T)
  try(saveRDS(cm_mi, file = paste0("test/CM_r_", rr, "_", model_name, "_", perc, ".rds")))
  cat(model_name, "ALL DONE", fill = T)
}

test_increasing <- function(perc, rr, models, df, target_var, positive){
  cat("============================================================", fill = T)
  cat("Round N", rr, " \t Missing  = " , perc * 100, "%", fill = T)
  df <- df
  for(model_name in models){
    cat(model_name, fill = T)
    generate_data(model_name, df = df, perc = perc, rr = rr, target_var = target_var, positive = positive)
  }
}

summarise_SE <- function(df, conf_interval, groupings, statistic){
  qstatistic <- enquo(statistic)
  tmp <- df %>% group_by(!!!groupings) %>% 
    summarise(N = n(), mean = mean(!!qstatistic), sd = sd(!!qstatistic))
  tmp$se <- tmp$sd / sqrt(tmp$N)
  tmp$ci <- tmp$se * qt(conf_interval / 2 + .5, tmp$N - 1)
  names(tmp)
  tmp
}

# ----------------------------------------------------------------
# ----------------------------------------------------------------
# ----------------------------------------------------------------

###### DATA PREPROCESSING

df <- data.frame(df)
names(df)[1] <- "age"

# CONVERTING COLUMNS TO FACTORS
cols <- c("sex", "cp", "fbs","restecg", "exang", "slope", "ca", "thal", "target")
for(col in names(df)){
  if(col %in% cols){
    df[, col] <- as.factor(df[, col])
  }
}

# LEVELS FOR FACTOR VARIABLE
levels(df$target) <- c("No", "Yes")
levels(df$sex) <- c("Female", "Male")
levels(df$fbs) <- c("False", "True")
levels(df$exang) <- c("No", "Yes")

###### IMPUTATIG MISSING DATA
models <- c("mi", "mice", "missForest")
repeats <- 5
probs <- c(1:5)/10
target <- "exang"
positive <- levels(df[, eval(target)])[2]
for(i in 1:repeats){
  for(p in probs){
    test_increasing(p, rr = i, models, df, target_var = target, positive = positive)
  }
}

###### COLLECTING RESULTS
res <- data.frame()
for(i in 1:repeats){
  for(p in probs){
    for(model_name in models){
      file <- paste0("test/CM_r_", i, "_", model_name, "_", p, ".rds")
      name <- paste0(model_name, "_CM")
      try(assign(name, readRDS(file)))
      tmp <- data.frame(t(get(name)$overall))
      tmp <- cbind(package = model_name, missing = p, rep = i, tmp)
      tmp <- cbind(tmp, data.frame(t(get(name)$byClass)) )
      res <- rbind(res, tmp)
    }
  }
}

grouping <- quos(package, missing)
acc <- summarise_SE(df = res, .95, statistic = Accuracy, grouping = grouping)
spec <- summarise_SE(df = res, .95, statistic = Specificity, grouping = grouping)
sens <- summarise_SE(df = res, .95, statistic = Sensitivity, grouping = grouping)
F1 <- summarise_SE(df = res, .95, statistic = F1, grouping = grouping)


######  GENERATING PLOTS
theme_set(theme_fivethirtyeight())
pd <- 0.04

p_acc <- ggplot(acc, aes(x=missing, y=mean, color=package)) + 
  geom_line(size = 1, position = position_dodge(pd)) + 
  geom_point(position = position_dodge(pd)) + 
  scale_fill_brewer(palette = "Spectral") +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se),
                size=.75,
                width=.02, position = position_dodge(pd)) +
  guides(fill = guide_legend(title =  element_blank(), ncol = 8)) +
  theme(legend.text = element_text(size = 10)) +
  xlab("Missing data (%)") +
  ylab("Accuracy (avg + se)") +
  labs(title="Accuracy") +
  expand_limits(y=c(0.4, 0.7)) +
  scale_y_continuous(breaks=0:10*0.1) +
  scale_x_continuous(breaks=0:10*0.1) +
  theme(legend.position="bottom") +
  facet_wrap( ~ package, ncol=4)

p_spec <- ggplot(spec, aes(x=missing, y=mean, color=package)) + 
  geom_line(size = 1, position = position_dodge(pd)) + 
  geom_point(position = position_dodge(pd)) + 
  scale_fill_brewer(palette = "Spectral") +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se),
                size=.75,
                width=.02, position = position_dodge(pd)) +
  guides(fill = guide_legend(title =  element_blank(), ncol = 8)) +
  theme(legend.text = element_text(size = 10)) +
  xlab("Missing data (%)") +
  ylab("Specificity (avg + se)") +
  labs(title="Specificity") +
  expand_limits(y=c(0.1, 0.8)) +
  scale_y_continuous(breaks=0:10*0.1) +
  scale_x_continuous(breaks=0:10*0.1) +
  theme(legend.position="bottom") +
  facet_wrap( ~ package, ncol=4)

p_sens <- ggplot(sens, aes(x=missing, y=mean, color=package)) + 
  geom_line(size = 1, position = position_dodge(pd)) + 
  geom_point(position = position_dodge(pd)) + 
  scale_fill_brewer(palette = "Spectral") +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se),
                size=.75,
                width=.02, position = position_dodge(pd)) +
  guides(fill = guide_legend(title =  element_blank(), ncol = 8)) +
  theme(legend.text = element_text(size = 10)) +
  xlab("Missing data (%)") +
  ylab("Sensitivity (avg + se)") +
  labs(title="Sensitivity") +
  expand_limits(y=c(0.2, 0.9)) +
  scale_y_continuous(breaks=0:10*0.1) +
  scale_x_continuous(breaks=0:10*0.1) + 
  theme(legend.position="bottom") +
  facet_wrap( ~ package, ncol=4)

p_F1 <- ggplot(F1, aes(x=missing, y=mean, color=package)) + 
  geom_line(size = 1, position = position_dodge(pd)) + 
  geom_point(position = position_dodge(pd)) + 
  scale_fill_brewer(palette = "Spectral") +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se),
                size=.75,
                width=.02, position = position_dodge(pd)) +
  guides(fill = guide_legend(title =  element_blank(), ncol = 8)) +
  theme(legend.text = element_text(size = 10)) +
  xlab("Missing data (%)") +
  ylab("F1 (avg + se)") +
  labs(title="F1") +
  expand_limits(y=c(0.2, 0.9)) +
  scale_y_continuous(breaks=0:10*0.1) +
  scale_x_continuous(breaks=0:10*0.1) + 
  theme(legend.position="bottom") +
  facet_wrap( ~ package, ncol=4)

library("gridExtra")
plots <- list(p_acc, p_spec, p_sens, p_F1)
# do.call(grid.arrange, plots)

grid.arrange(arrangeGrob(grobs = plots, top = "Statistics for cross-validated imputation techniques at increasing amounts of missing values"))
