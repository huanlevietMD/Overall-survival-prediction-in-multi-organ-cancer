setwd("D:/Huan/Project/Project multi organ/Upload Github")

#----------------------------------------------------------------
# Training set 
#----------------------------------------------------------------
# We used Python for removing Highly Correlated radiomics Features in the Training Set

# LASSO model was created to determine the most important 15 features 

data1 <- read.csv("Training set 95 radiomics features after Pearson correlation pairwise selection.csv",
                  header = TRUE)
attach(data1)
scale(data1)

glimpse(data1)
set.seed(100) 

library(plyr)
library(readr)
library(dplyr)
library(caret)
library(ggplot2)
library(repr)
library(glmnet)

ctrlspecs <- trainControl(method="cv", number = 10, savePredictions = "all")
lambda_vector <- 10^seq(5, -5, length=5000)
set.seed(100)
model1 <- train(Survival.time ~ .,
                data = data1,
                preProcess=c("center","scale"),
                method="glmnet",
                tuneGrid=expand.grid(alpha=1, lambda = lambda_vector),
                trControl=ctrlspecs,
                na.action=na.omit)
model1$bestTune
model1$bestTune$lambda
round(coef(model1$finalModel, model1$bestTune$lambda), 3)
varImp(model1)
h <- varImp(model1)
  ggplot(data = h, top = 15)+
  ggtitle("Top 15 Features Importance") +
  geom_bar(stat="identity", width=0.1)+
  geom_bar(stat="identity", color="steelblue", fill="steelblue", position="dodge")+
  theme_light(base_size = 15, base_family = "") +
  theme_minimal()



# Univariable Cox proportional hazard models

library(RegParallel)
library(survival)
data5 <- read.csv("Training set - top important 15 features.csv", header = TRUE)
attach(data5)
res5 <- RegParallel(
  data = data5,
  formula = 'Surv(Survival.time, deadstatus.event) ~ [*]',
  FUN = function(formula, data5)
    coxph(formula = formula,
          data = data5,
          ties = 'breslow',
          singular.ok = TRUE),
  FUNtype = 'coxph',
  variables = colnames(data5)[3:ncol(data5)],
  blocksize = 10,
  p.adjust = "BH")
res5 <- res5[!is.na(res5$P),]
res5
data6 = data.frame(res5)
write.csv(data6,'result after Univariable Cox in training set.csv', row.names = FALSE)


# Risk score calculation

library(survival)
library(survC1)
library(survival)
data6 <- read.csv("Training set 10 signature radiomics features.csv",
                  header = TRUE)
attach(data6)
data7 = Est.PH(data6)

data7

data8 = data.frame(data7["rs"])
median(data8[,1])

write.csv(data8,'riskscore Training set.csv', row.names = FALSE)

data10 <- read.csv("riskscore Training set.csv", header = TRUE)  
attach(data10)
riskscore = rs
riskscore[rs < 0.4308121] <- "low risk"
riskscore[rs > 0.4308121] <- "high risk"
data11 = data.frame(rs,riskscore)
data11
write.csv(data11,'Riskgroup Training set.csv', row.names = FALSE)

# Risk Score plot 

library(ggrisk)
library(ggplot2)
library(survival)
data9 <- read.csv("Training set 10 signature radiomics features.csv",
                  header = TRUE)
attach(data9)
Cox1 = coxph(Surv(Survival.time, deadstatus.event)~ .,
             data = data9)   
Cox1
plot8 <- ggrisk(
  Cox1,
  heatmap.genes = NULL,
  new.data = NULL,
  code.0 = "Alive",
  code.1 = "Dead",
  code.highrisk = "High",
  code.lowrisk = "Low",
  cutoff.show = FALSE,
  cutoff.value = "median",
  cutoff.x = NULL,
  cutoff.y = NULL,
  cutoff.label = NULL,
  title.A.ylab = "Risk Score",
  title.B.ylab = "Survival Time",
  title.A.legend = "Risk Group",
  title.B.legend = "Status",
  title.C.legend = "Expression",
  size.ABC = 1.5,
  size.ylab.title = 20,
  size.Atext = 18,
  size.Btext = 18,
  size.Ctext = 18,
  size.yticks = 0.5,
  size.yline = 0.5,
  size.points = 2,
  size.dashline = 1.2,
  size.cutoff = 5,
  size.legendtitle = 20,
  size.legendtext = 15,
  color.A = c(low = "blue", high = "red"),
  color.B = c(code.0 = "blue", code.1 = "red"),
  color.C = c(low = "blue", median = "white", high = "red"),
  vjust.A.ylab = 1,
  vjust.B.ylab = 2,
  family = "sans",
  expand.x = 10,
  relative_heights = c(0.5, 0.5, 0.06, 0.4)
)
plot8
ggsave("Heatmap Lung1 full.png", plot8, width = 17, height = 10,dpi = 300)


# Kaplan-Meier plot 

data35 <- read.csv("Training set Kaplan + Riskset.csv", header = TRUE)  
attach(data35)
library(survminer)
library(survival)
library(survival)

fit <- survfit(Surv(os, status) ~ riskscore,
               data = data35)
ggsurvplot(
  fit,                     # survfit object with calculated statistics.
  data = data35,  # data used to fit survival curves. 
  risk.table = TRUE,       # show risk table.
  pval = TRUE,
  pval.method=TRUE,# show p-value of log-rank test.
  pval.size = 8,
  conf.int = FALSE,         # show confidence intervals for 
  # point estimaes of survival curves.
  xlim = c(0,2000),        # present narrower X axis, but not affect
  # survival estimates.
  break.time.by = 500,     # break X axis in time intervals by 500.
  ggtheme = theme_minimal(), # customize plot and risk table with a theme.
  risk.table.y.text.col = T, # colour risk table text annotations.
  risk.table.y.text = FALSE, # show bars instead of names in text annotations
  font.title    = c(14, "bold", "black"),
  font.subtitle = c(25, "bold", "purple"),
  font.caption  = c(20, "plain", "orange"),
  font.x        = c(18, "plain","black"),
  font.y        = c(18,"plain", "black"),
  font.xtickslab = c(12, "plain", "black"),
  font.ytickslab = c(12, "plain", "black"),
  font.legend    = c(18,"plain", "black"),
  xlab = "Time in days",
  title="Lung 1 training set",
  legend.labs=c("High risk","Low risk"),
  legend.title="Risk Score",
) 


# Assess the radiomics signature and clinical 
# parameters integration's effectiveness (iAUC)

data1 <- read.csv("Training set iAUC compare.csv",
                  header = TRUE)
attach(data1)
library(survC1)
library(survival)
data1$female <- as.numeric(data1$Gender == "female")
data1$stageI <- as.numeric(data1$Stage == "I")
data1$stageII <- as.numeric(data1$Stage == "II")
data1$stageIIIa <- as.numeric(data1$Stage == "IIIa")

tau=2000
# iAUC
# Radiomics model
C=Inf.Cval(mydata = data1[,c("os","status","Riskscore")], tau, itr=1000)
round(c(C$Dhat, C$se, C$low95, C$upp95), digits=3)

# Clinical model
C=Inf.Cval(mydata = data1[,c("os","status","Age","female","stageI","stageII","stageIIIa")], tau, itr=1000)
round(c(C$Dhat, C$se, C$low95, C$upp95), digits=3)
# Combine model
C=Inf.Cval(mydata = data1[,c("os","status","Riskscore","Age","female","stageI","stageII","stageIIIa")], tau, itr=1000)
round(c(C$Dhat, C$se, C$low95, C$upp95), digits=3)

# Combine model vs Radiomics model
model0<-data1[,c(1:3)] ; 
model1<-data1[,c(1:4,7:10)]
covs1<-as.matrix(model1[,c(-1,-2)])
covs0<-as.matrix(model0[,c(-1,-2)])
Delta=Inf.Cval.Delta(model0[,1:2], covs0, covs1, tau, itr=1000)
round(Delta, digits=3)

# Radiomics model vs Clinical model
model0<-data1 [,c(1:2,4,7:10)]; 
model1<-data1[,c(1:3)]
covs1<-as.matrix(model1[,c(-1,-2)])
covs0<-as.matrix(model0[,c(-1,-2)])
tau=2000
Delta=Inf.Cval.Delta(model0[,1:2], covs0, covs1, tau, itr=1000)
round(Delta, digits=3)

# Combine model & Clinical model
model0<-data1 [,c(1:2,4,7:10)]; 
model1<-data1[,c(1:4,7:10)]
covs1<-as.matrix(model1[,c(-1,-2)])
covs0<-as.matrix(model0[,c(-1,-2)])
tau=2000
Delta=Inf.Cval.Delta(model0[,1:2], covs0, covs1, tau, itr=1000)
round(Delta, digits=3)


# iAUC plot

data10 <- read.csv("Training set Kaplan + Riskset.csv", header = TRUE)
attach(data10)
library(survivalROC)
library(ggplot2)
library(pROC)
library(survival)
library(MASS)
library(risksetROC)
library(ggthemes)
library(timeROC)
survival.time=data10$os
survival.status=data10$status
fit1 <- coxph( Surv(survival.time,survival.status)
               ~ Riskscore
               , data=data10)
eta1 <- fit1$linear.predictor

fit2 <- coxph( Surv(survival.time,survival.status)
               ~ Gender + Age + 
                 Stage 
               , data=data10)
eta2 <- fit2$linear.predictor

fit3 <- coxph( Surv(survival.time,survival.status)
               ~ Riskscore + Gender + Age + 
                 Stage 
               , data=data10)
eta3 <- fit3$linear.predictor

tmax=2000
AUC1 = risksetAUC(Stime=survival.time,
                  status=survival.status, marker=eta1, method="Cox", span = 0.4,
                  tmax =tmax, plot = FALSE);
AUC1

AUC2=risksetAUC(Stime=survival.time,
                status=survival.status, marker=eta2, method="Cox", tmax=tmax,
                span = 0.4, plot = FALSE);
AUC3=risksetAUC(Stime=survival.time,
                status=survival.status, marker=eta3, method="Cox", tmax=tmax,
                span = 0.4, plot = FALSE);
auc <- rbind(data.frame(AUC1, Data = 'Radiomics'),
             data.frame(AUC2, Data = 'Clinical'),
             data.frame(AUC3, Data = 'Radiomics + Clinical'))
auc$Data <- factor(auc$Data, levels = unique(auc$Data))
My_Theme = theme(
  axis.title.x = element_text(size = 16),
  axis.text.x = element_text(size = 14),
  axis.title.y = element_text(size = 16),
  legend.text = element_text(size = 16),
  legend.title = element_text(size = 16))

ggplot(auc, aes(x = utimes, y = AUC, color = Data)) +
  geom_line(size = 1.5) + labs(title="Lung 1 training set") +
  scale_size_continuous(guide = F) +
  scale_color_discrete(name = "Models",
                       labels = c(paste0('Radiomics'),
                                  paste0('Clinical'),
                                  paste0('Radiomics + Clinical'))) +
  labs(x ='time', y = 'AUC') +
  scale_y_continuous(limits = c(.5, 1)) + theme_minimal() + My_Theme +  xlim(0, 2000)


# Pred error
data11 <- read.csv("Training set Kaplan + Riskset.csv", header = TRUE)
attach(data11)
library(pec)
time=data11$os
status=data11$status

Models <- list("Clinical"=coxph(Surv(time,status)~Gender + Age + 
                                  Stage,data=data11,na.action = na.omit, x=TRUE,y=TRUE),
               "Radiomics"=coxph(Surv(time,status)~Riskscore,data=data11,na.action = na.omit,x=TRUE,y=TRUE),
               "Radiomics + Clinical"=coxph(Surv(time,status)~Riskscore+Gender + Age + 
                                              Stage ,na.action = na.omit, data=data11,x=TRUE,y=TRUE))
# compute the .632+ estimate of the generalization error

PredError <- pec(object=Models,
                 formula=Surv(time,status)~Riskscore+Gender + Age + 
                   Stage,
                 data=data11,
                 exact=TRUE,
                 cens.model="marginal",
                 splitMethod="none",
                 B=0,
                 verbose=TRUE,
                 na.action = na.fail,
                 maxtime=2000
)

print(PredError)
library(ggplot2)
# packages for the pipe and pivot_wider, you can do it with base functions, I just prefer these
library(tidyr)
library(dplyr)

df <- do.call(cbind, PredError[["AppErr"]]) 
df <- cbind(PredError[["time"]], df) 
colnames(df)[1] <- "time"             
df <- as.data.frame(df) %>% pivot_longer(cols = 2:last_col(), names_to = "Models", values_to = "PredError") 
My_Theme = theme(
  axis.title.x = element_text(size = 16),
  axis.text.x = element_text(size = 14),
  axis.title.y = element_text(size = 16),
  legend.text = element_text(size = 16),
  legend.title = element_text(size = 16))
ggplot(data = df, aes(x = time, y = PredError, color=Models),xlim=c(0,2000)) +
  geom_line(size = 1.5) + labs(title="Lung 1 training set") + theme_minimal()+
  xlim(0, 2000) + My_Theme + scale_color_manual(values=c('Green','tomato','Steelblue','Black'))




#--------------------------------------------------------------------
# In the testing and 2 validation sets
# we locked and independently assessed the 
# risk score created by the best radiomics 
# signatures model for overall survival prediction
#--------------------------------------------------------------------


# Testing set
# Lung 2

library(survival)
library(survC1)
library(survival)
data6 <- read.csv("Testing set 10 signature radiomics features.csv",
                  header = TRUE)
attach(data6)
data7 = Est.PH(data6)

data7

data8 = data.frame(data7["rs"])
median(data8[,1])

write.csv(data8,'riskscore testing set.csv', row.names = FALSE)

data10 <- read.csv("riskscore testing set.csv", header = TRUE)  
attach(data10)
riskscore = rs
riskscore[rs < -0.7948127] <- "low risk"
riskscore[rs > -0.7948127] <- "high risk"
data11 = data.frame(rs,riskscore)
data11
write.csv(data11,'Riskgroup Testing set.csv', row.names = FALSE)

# Kaplan Meier plot

data36 <- read.csv("Testing set Kaplan + Riskset.csv", header = TRUE)  
attach(data36)
library(survminer)
library(survival)
library(survival)


fit <- survfit(Surv(os, status) ~ riskscore,
               data = data36)
ggsurvplot(
  fit,                     # survfit object with calculated statistics.
  data = data36,  # data used to fit survival curves. 
  risk.table = TRUE,       # show risk table.
  pval = TRUE,
  pval.method=TRUE,# show p-value of log-rank test.
  pval.size = 8,
  conf.int = FALSE,         # show confidence intervals for 
  # point estimaes of survival curves.
  xlim = c(0,2000),        # present narrower X axis, but not affect
  # survival estimates.
  break.time.by = 500,     # break X axis in time intervals by 500.
  ggtheme = theme_minimal(), # customize plot and risk table with a theme.
  risk.table.y.text.col = T, # colour risk table text annotations.
  risk.table.y.text = FALSE, # show bars instead of names in text annotations
  font.title    = c(14, "bold", "black"),
  font.subtitle = c(25, "bold", "purple"),
  font.caption  = c(20, "plain", "orange"),
  font.x        = c(18, "plain","black"),
  font.y        = c(18,"plain", "black"),
  font.xtickslab = c(12, "plain", "black"),
  font.ytickslab = c(12, "plain", "black"),
  font.legend    = c(18,"plain", "black"),
  xlab = "Time in days",
  title="Lung 2 testing set",
  legend.labs=c("High risk","Low risk"),
  legend.title="Risk Score",
) 


# iAUC plot 

data10 <- read.csv("Testing set Kaplan + Riskset.csv", header = TRUE)
attach(data10)
library(survivalROC)
library(ggplot2)
library(pROC)
library(survival)
library(MASS)
library(risksetROC)
library(ggthemes)
survival.time=data10$os
survival.status=data10$status
fit1 <- coxph( Surv(survival.time,survival.status)
               ~ Riskscore
               , data=data10)
eta1 <- fit1$linear.predictor

fit2 <- coxph( Surv(survival.time,survival.status)
               ~ Gender + Age + 
                 Stage 
               , data=data10)
eta2 <- fit2$linear.predictor

fit3 <- coxph( Surv(survival.time,survival.status)
               ~ Riskscore + Gender + Age + 
                 Stage 
               , data=data10)
eta3 <- fit3$linear.predictor

tmax=2000
AUC1 = risksetAUC(Stime=survival.time,
                  status=survival.status, marker=eta1, method="Cox", span = 0.4,
                  tmax =tmax, plot = FALSE);
AUC2=risksetAUC(Stime=survival.time,
                status=survival.status, marker=eta2, method="Cox", tmax=tmax,
                span = 0.4, plot = FALSE);
AUC3=risksetAUC(Stime=survival.time,
                status=survival.status, marker=eta3, method="Cox", tmax=tmax,
                span = 0.4, plot = FALSE);
auc <- rbind(data.frame(AUC1, Data = 'Radiomics'),
             data.frame(AUC2, Data = 'Clinical'),
             data.frame(AUC3, Data = 'Radiomics + Clinical'))
auc$Data <- factor(auc$Data, levels = unique(auc$Data))
My_Theme = theme(
  axis.title.x = element_text(size = 16),
  axis.text.x = element_text(size = 14),
  axis.title.y = element_text(size = 16),
  legend.text = element_text(size = 16),
  legend.title = element_text(size = 16))

ggplot(auc, aes(x = utimes, y = AUC, color = Data)) +
  geom_line(size = 1.5) + labs(title="Lung 2 testing set") +
  scale_size_continuous(guide = F) +
  scale_color_discrete(name = "Models",
                       labels = c(paste0('Radiomics'),
                                  paste0('Clinical'),
                                  paste0('Radiomics + Clinical'))) +
  labs(x ='time', y = 'AUC') +
  scale_y_continuous(limits = c(.5, 1)) + theme_minimal() + My_Theme +  xlim(0, 2000)



# Pred error
data11 <- read.csv("Testing set Kaplan + Riskset.csv", header = TRUE)
attach(data11)
library(pec)
time=data11$os
status=data11$status

Models <- list("Clinical"=coxph(Surv(time,status)~Gender + Age + 
                                  Stage,data=data11,na.action = na.omit, x=TRUE,y=TRUE),
               "Radiomics"=coxph(Surv(time,status)~Riskscore,data=data11,na.action = na.omit,x=TRUE,y=TRUE),
               "Radiomics + Clinical"=coxph(Surv(time,status)~Riskscore+Gender + Age + 
                                              Stage ,na.action = na.omit, data=data11,x=TRUE,y=TRUE))
# compute the .632+ estimate of the generalization error

PredError <- pec(object=Models,
                 formula=Surv(time,status)~Riskscore+Gender + Age + 
                   Stage,
                 data=data11,
                 exact=TRUE,
                 cens.model="marginal",
                 splitMethod="none",
                 B=0,
                 verbose=TRUE,
                 na.action = na.fail,
                 maxtime=2000
)

print(PredError)
library(ggplot2)
# packages for the pipe and pivot_wider, you can do it with base functions, I just prefer these
library(tidyr)
library(dplyr)

df <- do.call(cbind, PredError[["AppErr"]]) # contains y values for each model
df <- cbind(PredError[["time"]], df) # values of the x axis
colnames(df)[1] <- "time"             
df <- as.data.frame(df) %>% pivot_longer(cols = 2:last_col(), names_to = "Models", values_to = "PredError") # pivot table to long format makes it easier to use ggplot
My_Theme = theme(
  axis.title.x = element_text(size = 16),
  axis.text.x = element_text(size = 14),
  axis.title.y = element_text(size = 16),
  legend.text = element_text(size = 16),
  legend.title = element_text(size = 16))
ggplot(data = df, aes(x = time, y = PredError, color=Models),xlim=c(0,2000)) +
  geom_line(size = 1.5)+ labs(title="Lung 2 testing set") + theme_minimal()+
  xlim(0, 2000) + My_Theme + scale_color_manual(values=c('Green','tomato','Steelblue','Black'))

# Assess the radiomics signature and clinical 
# parameters integration's effectiveness (iAUC) 

data1 <- read.csv("Testing set iAUC compare.csv",
                  header = TRUE)
attach(data1)
library(survC1)
library(survival)
data1$female <- as.numeric(data1$Gender == "Female")
data1$tis <- as.numeric(data1$Stage == "Tis")
data1$stageI <- as.numeric(data1$Stage == "I")
data1$stageII <- as.numeric(data1$Stage == "II")
data1$stageIIIa <- as.numeric(data1$Stage == "IIIA")
data1$stageIIIb <- as.numeric(data1$Stage == "IIIB")
tau=2000
# iAUC
# Radiomics model
C=Inf.Cval(mydata = data1[,c("os","status","Riskscore")], tau, itr=1000)
round(c(C$Dhat, C$se, C$low95, C$upp95), digits=3)

# Clinical model
C=Inf.Cval(mydata = data1[,c("os","status","Age","female","tis","stageI","stageII","stageIIIa"
                             ,"stageIIIb")], tau, itr=1000)
round(c(C$Dhat, C$se, C$low95, C$upp95), digits=3)

# Combine model
C=Inf.Cval(mydata = data1[,c("os","status","Riskscore","Age","female","tis","stageI","stageII","stageIIIa",
                             "stageIIIb")], tau, itr=1000)
round(c(C$Dhat, C$se, C$low95, C$upp95), digits=3)


# Radiomics model vs Combine model
model0<-data1[,c(1:3)] ; 
model1<-data1[,c(1:3,6:12)]
covs1<-as.matrix(model1[,c(-1,-2)])
covs0<-as.matrix(model0[,c(-1,-2)])
Delta=Inf.Cval.Delta(model0[,1:2], covs0, covs1, tau, itr=1000)
round(Delta, digits=3)

# Radiomics model & Clinical model
model0<-data1 [,c(1:2,6:12)]; 
model1<-data1[,c(1:3)]
covs1<-as.matrix(model1[,c(-1,-2)])
covs0<-as.matrix(model0[,c(-1,-2)])
tau=2000
Delta=Inf.Cval.Delta(model0[,1:2], covs0, covs1, tau, itr=1000)
round(Delta, digits=3)

# Combine model & Clinical model
model0<-data1 [,c(1:3,6:12)] ; 
model1<-data1[,c(1:2,6:12)]
covs1<-as.matrix(model1[,c(-1,-2)])
covs0<-as.matrix(model0[,c(-1,-2)])
tau=2000
Delta=Inf.Cval.Delta(model0[,1:2], covs1, covs0, tau, itr=1000)

round(Delta, digits=3)


#----------------------------------------------------------------------------
# Kidney external validation set
# Risk score calculation
library(survival)
library(survC1)
library(survival)
data6 <- read.csv("Kidney 10 signature radiomics features.csv",
                  header = TRUE)
attach(data6)
data7 = Est.PH(data6)
data7
data8 = data.frame(data7["rs"])
median(data8[,1])

write.csv(data8,'riskscore Kidney validation set.csv', row.names = FALSE)

data10 <- read.csv("riskscore Kidney validation set.csv", header = TRUE)  
attach(data10)
riskscore = rs
riskscore[rs < -1.149632] <- "low risk"
riskscore[rs > -1.149632] <- "high risk"
data11 = data.frame(rs,riskscore)
data11
write.csv(data11,'Riskgroup Kidney validation set.csv', row.names = FALSE)


# Kaplan Meier plot
data37 <- read.csv("Kidney Kaplan + Riskset.csv", header = TRUE)  
attach(data37)
library(survminer)
library(survival)
library(survival)


fit <- survfit(Surv(os, status) ~ riskscore,
               data = data37)
ggsurvplot(
  fit,                     # survfit object with calculated statistics.
  data = data37,  # data used to fit survival curves. 
  risk.table = TRUE,       # show risk table.
  pval = TRUE,
  pval.method=TRUE,# show p-value of log-rank test.
  pval.size = 8,
  conf.int = FALSE,         # show confidence intervals for 
  # point estimaes of survival curves.
  xlim = c(0,2000),        # present narrower X axis, but not affect
  # survival estimates.
  break.time.by = 500,     # break X axis in time intervals by 500.
  ggtheme = theme_minimal(), # customize plot and risk table with a theme.
  risk.table.y.text.col = T, # colour risk table text annotations.
  risk.table.y.text = FALSE, # show bars instead of names in text annotations
  font.title    = c(14, "bold", "black"),
  font.subtitle = c(25, "bold", "purple"),
  font.caption  = c(20, "plain", "orange"),
  font.x        = c(18, "plain","black"),
  font.y        = c(18,"plain", "black"),
  font.xtickslab = c(12, "plain", "black"),
  font.ytickslab = c(12, "plain", "black"),
  font.legend    = c(18,"plain", "black"),
  xlab = "Time in days",
  title="Kidney validation set",
  legend.labs=c("High risk","Low risk"),
  legend.title="Risk Score",
) 



# iAUC plot

data10 <- read.csv("Kidney Kaplan + Riskset.csv", header = TRUE)
attach(data10)
library(survivalROC)
library(ggplot2)
library(pROC)
library(survival)
library(MASS)
library(risksetROC)
library(ggthemes)
survival.time=data10$os
survival.status=data10$status
fit1 <- coxph( Surv(survival.time,survival.status)
               ~ Riskscore
               , data=data10)
eta1 <- fit1$linear.predictor

fit2 <- coxph( Surv(survival.time,survival.status)
               ~ Gender + Age + 
                 Stage 
               , data=data10)
eta2 <- fit2$linear.predictor

fit3 <- coxph( Surv(survival.time,survival.status)
               ~ Riskscore + Gender + Age + 
                 Stage 
               , data=data10)
eta3 <- fit3$linear.predictor

tmax=2000
AUC1 = risksetAUC(Stime=survival.time,
                  status=survival.status, marker=eta1, method="Cox", span = 0.4,
                  tmax =tmax, plot = FALSE);
AUC2=risksetAUC(Stime=survival.time,
                status=survival.status, marker=eta2, method="Cox", tmax=tmax,
                span = 0.4, plot = FALSE);
AUC3=risksetAUC(Stime=survival.time,
                status=survival.status, marker=eta3, method="Cox", tmax=tmax,
                span = 0.4, plot = FALSE);
auc <- rbind(data.frame(AUC1, Data = 'Radiomics'),
             data.frame(AUC2, Data = 'Clinical'),
             data.frame(AUC3, Data = 'Radiomics + Clinical'))
auc$Data <- factor(auc$Data, levels = unique(auc$Data))
My_Theme = theme(
  axis.title.x = element_text(size = 16),
  axis.text.x = element_text(size = 14),
  axis.title.y = element_text(size = 16),
  legend.text = element_text(size = 16),
  legend.title = element_text(size = 16))

ggplot(auc, aes(x = utimes, y = AUC, color = Data)) +
  geom_line(size = 1.5) + labs(title="Kidney validation set") +
  scale_size_continuous(guide = F) +
  scale_color_discrete(name = "Models",
                       labels = c(paste0('Radiomics'),
                                  paste0('Clinical'),
                                  paste0('Radiomics + Clinical'))) +
  labs(x ='time', y = 'AUC') +
  scale_y_continuous(limits = c(.5, 1)) + theme_minimal() + My_Theme +  xlim(0, 2000)

# Pred error
data11 <- read.csv("Kidney Kaplan + Riskset.csv", header = TRUE)
attach(data11)
library(pec)
time=data11$os
status=data11$status

Models <- list("Clinical"=coxph(Surv(time,status)~Gender + Age + 
                                  Stage,data=data11,na.action = na.omit, x=TRUE,y=TRUE),
               "Radiomics"=coxph(Surv(time,status)~Riskscore,data=data11,na.action = na.omit,x=TRUE,y=TRUE),
               "Radiomics + Clinical"=coxph(Surv(time,status)~Riskscore+Gender + Age + 
                                              Stage ,na.action = na.omit, data=data11,x=TRUE,y=TRUE))
# compute the .632+ estimate of the generalization error

PredError <- pec(object=Models,
                 formula=Surv(time,status)~Riskscore+Gender + Age + 
                   Stage,
                 data=data11,
                 exact=TRUE,
                 cens.model="marginal",
                 splitMethod="none",
                 B=0,
                 verbose=TRUE,
                 na.action = na.fail,
                 maxtime=2000
)

print(PredError)
library(ggplot2)
# packages for the pipe and pivot_wider, you can do it with base functions, I just prefer these
library(tidyr)
library(dplyr)

df <- do.call(cbind, PredError[["AppErr"]]) # contains y values for each model
df <- cbind(PredError[["time"]], df) # values of the x axis
colnames(df)[1] <- "time"             
df <- as.data.frame(df) %>% pivot_longer(cols = 2:last_col(), names_to = "Models", values_to = "PredError") # pivot table to long format makes it easier to use ggplot
My_Theme = theme(
  axis.title.x = element_text(size = 16),
  axis.text.x = element_text(size = 14),
  axis.title.y = element_text(size = 16),
  legend.text = element_text(size = 16),
  legend.title = element_text(size = 16))
ggplot(data = df, aes(x = time, y = PredError, color=Models),xlim=c(0,2000)) +
  geom_line(size = 1.5)+ labs(title="Kidney validation set") + theme_minimal()+
  xlim(0, 2000) + My_Theme + scale_color_manual(values=c('Green','tomato','Steelblue','Black'))


# Assess the radiomics signature and clinical 
# parameters integration's effectiveness (iAUC) 

data1 <- read.csv("Kidney iAUC compare.csv",
                  header = TRUE)
attach(data1)
library(survC1)
library(survival)
data1$female <- as.numeric(data1$Gender == "female")
data1$stageI <- as.numeric(data1$Stage == "I")
data1$stageII <- as.numeric(data1$Stage == "II")
data1$stageIII <- as.numeric(data1$Stage == "III")

tau=2000
# iAUC
# Radiomics model
C=Inf.Cval(mydata = data1[,c("os","status","Riskscore")], tau, itr=1000)
round(c(C$Dhat, C$se, C$low95, C$upp95), digits=3)

# Clinical model
C=Inf.Cval(mydata = data1[,c("os","status","Age","female","stageI","stageII","stageIII"
)], tau, itr=1000)
round(c(C$Dhat, C$se, C$low95, C$upp95), digits=3)

# Combine model
C=Inf.Cval(mydata = data1[,c("os","status","Riskscore","Age","female","stageI","stageII","stageIII"
)], tau, itr=1000)
round(c(C$Dhat, C$se, C$low95, C$upp95), digits=3)

# Radiomics model vs Combine model
model0<-data1[,c(1:2,6)] ; 
model1<-data1[,c(1:3,6:10)]
covs1<-as.matrix(model1[,c(-1,-2)])
covs0<-as.matrix(model0[,c(-1,-2)])
Delta=Inf.Cval.Delta(model0[,1:2], covs0, covs1, tau, itr=1000)
round(Delta, digits=3)

# Radiomics model vs Clinical model
model0<-data1 [,c(1:3,7:10)]; 
model1<-data1[,c(1:2,6)]
covs1<-as.matrix(model1[,c(-1,-2)])
covs0<-as.matrix(model0[,c(-1,-2)])
tau=2000
Delta=Inf.Cval.Delta(model0[,1:2], covs0, covs1, tau, itr=1000)
round(Delta, digits=3)

# Combine model vs Clinical model
model0<-data1 [,c(1:3,7:10)]; 
model1<-data1[,c(1:3,6:10)]
covs1<-as.matrix(model1[,c(-1,-2)])
covs0<-as.matrix(model0[,c(-1,-2)])
tau=2000
Delta=Inf.Cval.Delta(model0[,1:2], covs0, covs1, tau, itr=1000)
round(Delta, digits=3)


#----------------------------------------------------------------------------
# Head and neck validation set
# Risk score calculation
library(survival)
library(survC1)
library(survival)
data6 <- read.csv("H&N 10 signature radiomics features.csv",
                  header = TRUE)
attach(data6)
data7 = Est.PH(data6)
data7
data8 = data.frame(data7["rs"])
median(data8[,1])

write.csv(data8,'riskscore H&N validation set.csv', row.names = FALSE)

data10 <- read.csv("riskscore H&N validation set.csv", header = TRUE)  
attach(data10)
riskscore = rs
riskscore[rs < 1.10654] <- "low risk"
riskscore[rs > 1.10654] <- "high risk"
data11 = data.frame(rs,riskscore)
data11
write.csv(data11,'Riskgroup H&N validation set.csv', row.names = FALSE)

# Kaplan Meier plot

data38 <- read.csv("H&N Kaplan + Riskset.csv", header = TRUE)  
attach(data38)
library(survminer)
library(survival)
library(survival)


fit <- survfit(Surv(os, status) ~ riskscore,
               data = data38)
ggsurvplot(
  fit,                     # survfit object with calculated statistics.
  data = data38,  # data used to fit survival curves. 
  risk.table = TRUE,       # show risk table.
  pval = TRUE,
  pval.method=TRUE,# show p-value of log-rank test.
  pval.size = 8,
  conf.int = FALSE,         # show confidence intervals for 
  # point estimaes of survival curves.
  xlim = c(0,2000),        # present narrower X axis, but not affect
  # survival estimates.
  break.time.by = 500,     # break X axis in time intervals by 500.
  ggtheme = theme_minimal(), # customize plot and risk table with a theme.
  risk.table.y.text.col = T, # colour risk table text annotations.
  risk.table.y.text = FALSE, # show bars instead of names in text annotations
  font.title    = c(14, "bold", "black"),
  font.subtitle = c(25, "bold", "purple"),
  font.caption  = c(20, "plain", "orange"),
  font.x        = c(18, "plain","black"),
  font.y        = c(18,"plain", "black"),
  font.xtickslab = c(12, "plain", "black"),
  font.ytickslab = c(12, "plain", "black"),
  font.legend    = c(18,"plain", "black"),
  xlab = "Time in days",
  title="Head & Neck validation set",
  legend.labs=c("High risk","Low risk"),
  legend.title="Risk Score",
) 


# iAUC plot
data10 <- read.csv("H&N Kaplan + Riskset.csv", header = TRUE)
attach(data10)
library(survivalROC)
library(ggplot2)
library(pROC)
library(survival)
library(MASS)
library(risksetROC)
library(ggthemes)
survival.time=data10$os
survival.status=data10$status
fit1 <- coxph( Surv(survival.time,survival.status)
               ~ Riskscore
               , data=data10)
eta1 <- fit1$linear.predictor

fit2 <- coxph( Surv(survival.time,survival.status)
               ~ Gender + Age + 
                 Stage 
               , data=data10)
eta2 <- fit2$linear.predictor

fit3 <- coxph( Surv(survival.time,survival.status)
               ~ Riskscore + Gender + Age + 
                 Stage 
               , data=data10)
eta3 <- fit3$linear.predictor

tmax=2000
AUC1 = risksetAUC(Stime=survival.time,
                  status=survival.status, marker=eta1, method="Cox", span = 0.4,
                  tmax =tmax, plot = FALSE);
AUC2=risksetAUC(Stime=survival.time,
                status=survival.status, marker=eta2, method="Cox", tmax=tmax,
                span = 0.4, plot = FALSE);
AUC3=risksetAUC(Stime=survival.time,
                status=survival.status, marker=eta3, method="Cox", tmax=tmax,
                span = 0.4, plot = FALSE);
auc <- rbind(data.frame(AUC1, Data = 'Radiomics'),
             data.frame(AUC2, Data = 'Clinical'),
             data.frame(AUC3, Data = 'Radiomics + Clinical'))
auc$Data <- factor(auc$Data, levels = unique(auc$Data))
My_Theme = theme(
  axis.title.x = element_text(size = 16),
  axis.text.x = element_text(size = 14),
  axis.title.y = element_text(size = 16),
  legend.text = element_text(size = 16),
  legend.title = element_text(size = 16))

ggplot(auc, aes(x = utimes, y = AUC, color = Data)) +
  geom_line(size = 1.5) + labs(title="Head & Neck validation set") +
  scale_size_continuous(guide = F) +
  scale_color_discrete(name = "Models",
                       labels = c(paste0('Radiomics'),
                                  paste0('Clinical'),
                                  paste0('Radiomics + Clinical'))) +
  labs(x ='time', y = 'AUC') +
  scale_y_continuous(limits = c(.5, 1)) + theme_minimal() + My_Theme +  xlim(0, 2000)

# Pred error
data11 <- read.csv("H&N Kaplan + Riskset.csv", header = TRUE)
attach(data11)
library(pec)
time=data11$os
status=data11$status

Models <- list("Clinical"=coxph(Surv(time,status)~Gender + Age + 
                                  Stage,data=data11,na.action = na.omit, x=TRUE,y=TRUE),
               "Radiomics"=coxph(Surv(time,status)~Riskscore,data=data11,na.action = na.omit,x=TRUE,y=TRUE),
               "Radiomics + Clinical"=coxph(Surv(time,status)~Riskscore+Gender + Age + 
                                              Stage ,na.action = na.omit, data=data11,x=TRUE,y=TRUE))
# compute the .632+ estimate of the generalization error

PredError <- pec(object=Models,
                 formula=Surv(time,status)~Riskscore+Gender + Age + 
                   Stage,
                 data=data11,
                 exact=TRUE,
                 cens.model="marginal",
                 splitMethod="none",
                 B=0,
                 verbose=TRUE,
                 na.action = na.fail,
                 maxtime=2000
)

print(PredError)
library(ggplot2)
# packages for the pipe and pivot_wider, you can do it with base functions, I just prefer these
library(tidyr)
library(dplyr)

df <- do.call(cbind, PredError[["AppErr"]]) # contains y values for each model
df <- cbind(PredError[["time"]], df) # values of the x axis
colnames(df)[1] <- "time"             
df <- as.data.frame(df) %>% pivot_longer(cols = 2:last_col(), names_to = "Models", values_to = "PredError") # pivot table to long format makes it easier to use ggplot
My_Theme = theme(
  axis.title.x = element_text(size = 16),
  axis.text.x = element_text(size = 14),
  axis.title.y = element_text(size = 16),
  legend.text = element_text(size = 16),
  legend.title = element_text(size = 16))
ggplot(data = df, aes(x = time, y = PredError, color=Models),xlim=c(0,2000)) +
  geom_line(size = 1.5)+ labs(title="Head & Neck validation set") + theme_minimal()+
  xlim(0, 2000) + My_Theme + scale_color_manual(values=c('Green','tomato','Steelblue','Black'))

# Assess the radiomics signature and clinical 
# parameters integration's effectiveness (iAUC) 

data1 <- read.csv("H&N iAUC compare.csv",
                  header = TRUE)
attach(data1)
library(survC1)
library(survival)
data1$female <- as.numeric(data1$Gender == "female")
data1$stageI <- as.numeric(data1$Stage == "i")
data1$stageII <- as.numeric(data1$Stage == "ii")
data1$stageIII <- as.numeric(data1$Stage == "iii")
data1$stageIVa <- as.numeric(data1$Stage == "iva")
data1$stageIVb <- as.numeric(data1$Stage == "ivb")

tau=2000
# iAUC
# Radiomics model
C=Inf.Cval(mydata = data1[,c("os","status","Riskscore")], tau, itr=1000)
round(c(C$Dhat, C$se, C$low95, C$upp95), digits=3)
# Clinical model
C=Inf.Cval(mydata = data1[,c("os","status","Age","female","stageI","stageII","stageIII",
                             "stageIVa","stageIVb")], tau, itr=1000)
round(c(C$Dhat, C$se, C$low95, C$upp95), digits=3)
# Combine model
C=Inf.Cval(mydata = data1[,c("os","status","Riskscore","Age","female","stageI","stageII","stageIII",
                             "stageIVa","stageIVb")], tau, itr=1000)
round(c(C$Dhat, C$se, C$low95, C$upp95), digits=3)


# Radiomics model vs Combine model

model0<-data1[,c(1:2,6)] ; 
model1<-data1[,c(1:3,6:12)]
covs1<-as.matrix(model1[,c(-1,-2)])
covs0<-as.matrix(model0[,c(-1,-2)])
Delta=Inf.Cval.Delta(model0[,1:2], covs0, covs1, tau, itr=1000)
round(Delta, digits=3)

# Radiomics model vs Clinical model
model0<-data1 [,c(1:3,7:12)]; 
model1<-data1[,c(1:2,6)]
covs1<-as.matrix(model1[,c(-1,-2)])
covs0<-as.matrix(model0[,c(-1,-2)])
tau=2000
Delta=Inf.Cval.Delta(model0[,1:2], covs0, covs1, tau, itr=1000)
round(Delta, digits=3)

# Combine model vs Clinical model
model0<-data1 [,c(1:3,7:12)] ; 
model1<-data1[,c(1:3,6:12)]
covs1<-as.matrix(model1[,c(-1,-2)])
covs0<-as.matrix(model0[,c(-1,-2)])
tau=2000
Delta=Inf.Cval.Delta(model0[,1:2], covs0, covs1, tau, itr=1000)
round(Delta, digits=3)


