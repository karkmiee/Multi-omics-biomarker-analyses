### Title: Lasso Cox proportional hazards model ----

#This script is used to identify potential cancer-predicting biomarkers of Lynch Syndrome and evaluate their cancer-predicting capacity. 

#Data: The data includes next generation sequencing serum microRNA & metabolite content
#and information about whether the study subject developed cancer during a 5.8-year surveillance period. 

#-------------------------------------------------------------------------------
##Step 1. Load libraries and import data

#Install libraries
library(tidyverse)
library(glmnet)
library(survival)
library(mixOmics) #to impute missing values
library(survminer)
library(SurvMetrics) # to get all the metrics
library(pec) # to make predictions based on Cox model
library(survivalROC)

totalx <- read.delim("multi_omics_Filtered2.txt", header = TRUE)

#If needed impute missing values, NIPALS is used to decompose the dataset.
sum(is.na(totalx)) # number of cells with NA
totalx <- impute.nipals(X = totalx, ncomp = 10)
sum(is.na(totalx)) # number of cells with NA

#Rename responce variables
totalx <- totalx %>%
  dplyr::rename("status" = "Cancer_during_surveillance")
totalx <- totalx %>%
  dplyr::rename("time" = "Survival_time_years")


#-------------------------------------------------------------------------------
##Step 2. Feature selection using Lasso cox

#c-miRNAs

# Extract predictors
x0 <- data.matrix(totalx[ ,c(7:21)]) %>%  
  scale(center = T) %>%
  na.omit()

# Extract response variable
y0 <- totalx[,c(5,6)] %>% 
  na.omit() %>% 
  as.matrix()


# Fit the LASSO model (Lasso: Alpha = 1)

# Inspect Lasso  fit lambdas
fit0 <- glmnet(x0,y0, family = "cox", alpha = 1, maxit=1000000)
plot(fit0, xvar="lambda")
print(fit0)

lambda0 <- coef(fit0, s = 0.040680)  #select Lambda value where the number of features are downregularized to 5

#The selected features and their coefficients can be obtained:
nonZeroIdx0<-which(lambda0[,1]!= 0) 
features0<-rownames(lambda0)[nonZeroIdx0] 
features0


#c-Metab

# Extract predictors
x1 <- data.matrix(totalx[ ,c(22:85)]) %>%  
  scale(center = T) %>%
  na.omit()


# Fit the LASSO model (Lasso: Alpha = 1)

# Inspect Lasso  fit lambdas
fit1 <- glmnet(x1,y0, family = "cox", alpha = 1, maxit=1000000)
plot(fit1, xvar="lambda")
print(fit1)

lambda1 <- coef(fit1, s = 0.051400)  #select Lambda value where the number of features are downregularized to 5

#The selected features and their coefficients can be obtained:
nonZeroIdx1<-which(lambda1[,1]!= 0) 
features1<-rownames(lambda1)[nonZeroIdx1] 
features1


#-------------------------------------------------------------------------------
##Step 3. Scale cmiRNAs and metabolites to same distributions and extract df for CRC only

total1 <- totalx[,-c(1:4)] #remove phenodata from the df
total1[, -c(1, 2)] <- scale(total1[, -c(1, 2)], 
                            center = TRUE, 
                            scale = TRUE)
#Table for CRC only
total2 <- total1[!(rownames(total1) %in% c("Sample113", "Sample112", "Sample110", "Sample109", "Sample107", "Sample105", "Sample104", "Sample102", "Sample21")), ]

#-------------------------------------------------------------------------------
##Step 4. Fit the cox model with 10 biomarkers for all LS related cancers

full.cox <- coxph(Surv(time,status) ~ hsa.miR.101.3p + hsa.miR.182.5p + hsa.miR.183.5p + hsa.miR.4732.3p + hsa.miR.148b.3p + HDL_TG + Tyr + Glucose + Acetate + GlycA, data = total1, x=TRUE)
summary(full.cox)
anova(full.cox) #leave only most promising predictive features

ggforest(full.cox, data = total1) #plot the results

#Check model residual assumptions.
#For each covariate, the function cox.zph() correlates the corresponding set of scaled Schoenfeld residuals with time, to test for independence between residuals and time. 
#Additionally, it performs a global test for the model as a whole.
test.ph <- cox.zph(full.cox) 
test.ph 

#Repeat same for the reduced model 
full.cox.1 <- coxph(Surv(time,status) ~ hsa.miR.101.3p + hsa.miR.183.5p + HDL_TG, data = total1, x=TRUE)
summary(full.cox.1)
anova(full.cox.1)

ggforest(full.cox.1, data = total1)

test.ph <- cox.zph(full.cox.1)
test.ph

#Create a custom function to fix an error in residual plots produced by the ggcoxzph function from the survminer package, which visualizes the Schoenfeld residuals for Cox proportional hazards models.
ggcoxzphFixed <- function (fit, resid = TRUE, se = TRUE, df = 4, nsmo = 40, var, 
                           point.col = "red", point.size = 1, point.shape = 19, point.alpha = 1, 
                           caption = NULL, ggtheme = theme_survminer(), ...) 
{
  x <- fit
  if (!methods::is(x, "cox.zph")) 
    stop("Can't handle an object of class ", class(x))
  xx <- x$x
  yy <- x$y
  d <- nrow(yy)
  df <- max(df)
  nvar <- ncol(yy)
  pred.x <- seq(from = min(xx), to = max(xx), length = nsmo)
  temp <- c(pred.x, xx)
  lmat <- splines::ns(temp, df = df, intercept = TRUE)
  pmat <- lmat[1:nsmo, ]
  xmat <- lmat[-(1:nsmo), ]
  qmat <- qr(xmat)
  if (qmat$rank < df) 
    stop("Spline fit is singular, try a smaller degrees of freedom")
  if (se) {
    bk <- backsolve(qmat$qr[1:df, 1:df], diag(df))
    xtx <- bk %*% t(bk)
    seval <- ((pmat %*% xtx) * pmat) %*% rep(1, df)
  }
  ylab <- paste("Beta(t) for", dimnames(yy)[[2]])
  if (missing(var)) 
    var <- 1:nvar
  else {
    if (is.character(var)) 
      var <- match(var, dimnames(yy)[[2]])
    if (any(is.na(var)) || max(var) > nvar || min(var) < 
        1) 
      stop("Invalid variable requested")
  }
  if (x$transform == "log") {
    xx <- exp(xx)
    pred.x <- exp(pred.x)
  }
  else if (x$transform != "identity") {
    xtime <- as.numeric(dimnames(yy)[[1]])
    indx <- !duplicated(xx)
    apr1 <- approx(xx[indx], xtime[indx], seq(min(xx), max(xx), 
                                              length = 17)[2 * (1:8)])
    temp <- signif(apr1$y, 2)
    apr2 <- approx(xtime[indx], xx[indx], temp)
    xaxisval <- apr2$y
    xaxislab <- rep("", 8)
    for (i in 1:8) xaxislab[i] <- format(temp[i])
  }
  plots <- list()
  plots <- lapply(var, function(i) {
    invisible(pval <- round(x$table[i, 3], 4))
    gplot <- ggplot() + labs(title = paste0("Schoenfeld Individual Test p: ", 
                                            pval)) + ggtheme
    y <- yy[, i]
    yhat <- as.vector(pmat %*% qr.coef(qmat, y))
    if (resid) 
      yr <- range(yhat, y)
    else yr <- range(yhat)
    if (se) {
      temp <- as.vector(2 * sqrt(x$var[i, i] * seval))
      yup <- yhat + temp
      ylow <- yhat - temp
      yr <- range(yr, yup, ylow)
    }
    if (x$transform == "identity") {
      gplot <- gplot + geom_line(aes(x = pred.x, y = yhat)) + 
        xlab("Time") + ylab(ylab[i]) + ylim(yr)
    }
    else if (x$transform == "log") {
      gplot <- gplot + geom_line(aes(x = log(pred.x), 
                                     y = yhat)) + xlab("Time") + ylab(ylab[i]) + 
        ylim(yr)
    }
    else {
      gplot <- gplot + geom_line(aes(x = pred.x, y = yhat)) + 
        xlab("Time") + ylab(ylab[i]) + scale_x_continuous(breaks = xaxisval, 
                                                          labels = xaxislab) + ylim(yr)
    }
    if (resid) 
      gplot <- gplot + geom_point(aes(x = xx, y = y), 
                                  col = point.col, shape = point.shape, size = point.size, 
                                  alpha = point.alpha)
    if (se) {
      gplot <- gplot + geom_line(aes(x = pred.x, y = yup), 
                                 lty = "dashed") + geom_line(aes(x = pred.x, 
                                                                 y = ylow), lty = "dashed")
    }
    ggpubr::ggpar(gplot, ...)
  })
  names(plots) <- var
  class(plots) <- c("ggcoxzph", "ggsurv", "list")
  if ("GLOBAL" %in% rownames(x$table)) 
    global_p <- x$table["GLOBAL", 3]
  else global_p <- NULL
  attr(plots, "global_pval") <- global_p
  attr(plots, "caption") <- caption
  plots
}


test.ph <- cox.zph(full.cox.1)
test.ph

# Test the new function with your object
gg1 <- ggcoxzphFixed(test.ph)

gg1 ## provides plots with correct confidence bands and improved y axes
plot(test.ph)

#-------------------------------------------------------------------------------
##Step 5. Fit the cox model with 10 biomarkers for colorectal cancer

full.cox.CRC <- coxph(Surv(time,status) ~ hsa.miR.101.3p + hsa.miR.182.5p + hsa.miR.183.5p + hsa.miR.4732.3p + hsa.miR.148b.3p + HDL_TG + Tyr + Glucose + Acetate + GlycA, data = total2, x=TRUE)
summary(full.cox.CRC)
anova(full.cox.CRC)

test.ph <- cox.zph(full.cox.CRC)
test.ph

ggforest(full.cox.CRC, data = total2)

#Repeat same for the reduced model
full.cox.CRC2 <- coxph(Surv(time,status) ~ hsa.miR.101.3p + hsa.miR.183.5p + HDL_TG, data = total2, x=TRUE)
summary(full.cox.CRC2)
anova(full.cox.CRC2)

test.ph.CRC <- cox.zph(full.cox.CRC2)
test.ph.CRC

ggforest(full.cox.CRC2, data = total2)


#-------------------------------------------------------------------------------

##Step 6. Model validation: Built prediction models using 5 data splits.

##Fit the Cox Proportional Hazards Model on train data and make predictions with test data.


train0 <- read.csv("train_set_0.csv")
train0<-train0[,-c(1:5)]
train1 <- read.csv("train_set_1.csv")
train1<-train1[,-c(1:5)]
train2 <- read.csv("train_set_2.csv")
train2<-train2[,-c(1:5)]
train3 <- read.csv("train_set_3.csv")
train3<-train3[,-c(1:5)]
train4 <- read.csv("train_set_4.csv")
train4<-train4[,-c(1:5)]

#rename columns
train_data0 <- train0 %>%
  dplyr::rename("status" = "Cancer_during_surveillance")
train_data0 <- train_data0 %>%
  dplyr::rename("time" = "Survival_time_years")

train_data1 <- train1 %>%
  dplyr::rename("status" = "Cancer_during_surveillance")
train_data1 <- train_data1 %>%
  dplyr::rename("time" = "Survival_time_years")

train_data2 <- train2 %>%
  dplyr::rename("status" = "Cancer_during_surveillance")
train_data2 <- train_data2 %>%
  dplyr::rename("time" = "Survival_time_years")

train_data3 <- train3 %>%
  dplyr::rename("status" = "Cancer_during_surveillance")
train_data3 <- train_data3 %>%
  dplyr::rename("time" = "Survival_time_years")

train_data4 <- train4 %>%
  dplyr::rename("status" = "Cancer_during_surveillance")
train_data4 <- train_data4 %>%
  dplyr::rename("time" = "Survival_time_years")


# Standardize the data (excluding the response variable)
train_data_scaled0 <- scale(train_data0[, -c(1, 2)], 
                            center = TRUE, 
                            scale = TRUE)
train_means0 <- attr(train_data_scaled0, "scaled:center") #save train means for scaling the test data
train_sds0 <- attr(train_data_scaled0, "scaled:scale") #save train sds for scaling the test data

#Convert the scaled matrix to a data frame
train_data_scaled0 <- as.data.frame(train_data_scaled0)

#Add the unscaled 'time' and 'status' variables back to the scaled data frame
train_data_scaled0 <- cbind(train_data0[, c("time", "status")], train_data_scaled0)

#Train the model with predictive features
full.cox0 <- coxph(Surv(time,status) ~ hsa.miR.101.3p + hsa.miR.183.5p + HDL_TG, data = train_data_scaled0, x=TRUE)
summary(full.cox0)

test.ph0 <- cox.zph(full.cox0)
test.ph0

ggforest(full.cox0, data = train_data_scaled0)


#Make model predictions with test data
test_data0 <- read.csv("test_set_0.csv")
test_data0 <- test_data0[,-(1:5)]

test_data0 <- test_data0 %>%
  dplyr::rename("status" = "Cancer_during_surveillance")
test_data0 <- test_data0 %>%
  dplyr::rename("time" = "Survival_time_years")

test_data0[, -c(1, 2)] <- scale(test_data0[, -c(1, 2)], 
                                center = train_means0, 
                                scale = train_sds0)


# Extract event times and event indicators
event_times <- full.cox0$y[, "time"]
event_indicator <- full.cox0$y[, "status"]  # This indicates if the event occurred (1) or was censored (0)
distime <- event_times[event_indicator == 1]
distime <- sort(unique(distime)) # Make sure distime is sorted
#distime
mat_cox <- predictSurvProb(full.cox0, test_data0, distime) #Get the survival probability matrix
med_index <- median(1:length(distime)) #The index of median survival time of events
vec_cox <- mat_cox[ ,med_index]

vec_cox


times <- test_data0$time
status <- test_data0$status


#CI BS IBS IAE ISE based on Cox model: Non-standard model input methods
Cindex_cox <- Cindex(Surv(times, status), vec_cox)
BS_cox <- Brier(Surv(times, status), vec_cox, distime[med_index])
IBS_cox <- IBS(Surv(times, status), mat_cox, distime)
IAE_cox <- IAEISE(Surv(times, status), mat_cox, distime)[1]
ISE_cox <- IAEISE(Surv(times, status), mat_cox, distime)[2]
Cindex_cox
BS_cox
IBS_cox
IAE_cox
ISE_cox


##------------------------------------------------------------------------------
#Iteration 1:

###Fit the Cox Proportional Hazards Model on Train Data and make predictions with test data

# Standardize the data (excluding the response variable)
train_data_scaled1 <- scale(train_data1[, -c(1, 2)], 
                            center = TRUE, 
                            scale = TRUE)
train_means1 <- attr(train_data_scaled1, "scaled:center") #save train means for scaling the test data
train_sds1 <- attr(train_data_scaled1, "scaled:scale") #save train sds for scaling the test data

#Convert the scaled matrix to a data frame
train_data_scaled1 <- as.data.frame(train_data_scaled1)

#Add the unscaled 'time' and 'status' variables back to the scaled data frame
train_data_scaled1 <- cbind(train_data1[, c("time", "status")], train_data_scaled1)

#Train the model with predictive features
full.cox1 <- coxph(Surv(time,status) ~ hsa.miR.101.3p + hsa.miR.183.5p + HDL_TG, data = train_data_scaled1, x=TRUE)
summary(full.cox1)

test.ph1 <- cox.zph(full.cox1)
test.ph1


#Plot predictive features' impact on event hazard
ggforest(full.cox1, data = train_data_scaled1)


#Make model predictions with test data
test_data1 <- read.csv("test_set_1.csv")
test_data1 <- test_data1[,-(1:5)]

test_data1 <- test_data1 %>%
  dplyr::rename("status" = "Cancer_during_surveillance")
test_data1 <- test_data1 %>%
  dplyr::rename("time" = "Survival_time_years")

test_data1[, -c(1, 2)] <- scale(test_data1[, -c(1, 2)], 
                                center = train_means1, 
                                scale = train_sds1)

# Extract event times and event indicators
event_times1 <- full.cox1$y[, "time"]
event_indicator1 <- full.cox1$y[, "status"]  # This indicates if the event occurred (1) or was censored (0)
distime1 <- event_times1[event_indicator1 == 1]
distime1 <- sort(unique(distime1)) # Make sure distime is sorted
distime1
mat_cox1 <- predictSurvProb(full.cox1, test_data1, distime1) #get the survival probability matrix
med_index1 <- median(1:length(distime1)) #the index of median survival time of events
vec_cox1 <- mat_cox1[ ,med_index1]

times1 <- test_data1$time
status1 <- test_data1$status

#CI BS IBS IAE ISE based on Cox model: Non-standard model input methods
Cindex_cox1 <- Cindex(Surv(times1, status1), vec_cox1)
BS_cox1 <- Brier(Surv(times1, status1), vec_cox1, distime1[med_index1])
IBS_cox1 <- IBS(Surv(times1, status1), mat_cox1, distime1)
IAE_cox1 <- IAEISE(Surv(times1, status1), mat_cox1, distime1)[1]
ISE_cox1 <- IAEISE(Surv(times1, status1), mat_cox1, distime1)[2]
Cindex_cox1
BS_cox1
IBS_cox1
IAE_cox1
ISE_cox1


#-------------------------------------------------------------------------------
#Iteration 2:
###Fit the Cox Proportional Hazards Model on Train Data and make predictions with test data

# Standardize the data (excluding the response variable)
train_data_scaled2 <- scale(train_data2[, -c(1, 2)], 
                            center = TRUE, 
                            scale = TRUE)
train_means2 <- attr(train_data_scaled2, "scaled:center") #save train means for scaling the test data
train_sds2 <- attr(train_data_scaled2, "scaled:scale") #save train sds for scaling the test data

#Convert the scaled matrix to a data frame
train_data_scaled2 <- as.data.frame(train_data_scaled2)

#Add the unscaled 'time' and 'status' variables back to the scaled data frame
train_data_scaled2 <- cbind(train_data2[, c("time", "status")], train_data_scaled2)

#Train the model with predictive features

full.cox2 <- coxph(Surv(time,status) ~ hsa.miR.101.3p + hsa.miR.183.5p + HDL_TG, data = train_data_scaled2, x=TRUE)
summary(full.cox2)

test.ph2 <- cox.zph(full.cox2)
test.ph2

#Plot predictive features' impact on event hazard
ggforest(full.cox2, data = train_data_scaled2)


#Make model predictions with test data
test_data2 <- read.csv("test_set_2.csv")
test_data2 <- test_data2[,-(1:5)]

test_data2 <- test_data2 %>%
  dplyr::rename("status" = "Cancer_during_surveillance")
test_data2 <- test_data2 %>%
  dplyr::rename("time" = "Survival_time_years")

test_data2[, -c(1, 2)] <- scale(test_data2[, -c(1, 2)], 
                                center = train_means2, 
                                scale = train_sds2)


# Extract event times and event indicators
event_times2 <- full.cox2$y[, "time"]
event_indicator2 <- full.cox2$y[, "status"]  # This indicates if the event occurred (1) or was censored (0)
distime2 <- event_times2[event_indicator2 == 1]
distime2 <- sort(unique(distime2)) # Make sure distime is sorted
distime2
mat_cox2 <- predictSurvProb(full.cox2, test_data2, distime2) #get the survival probability matrix
med_index2 <- median(1:length(distime2)) #the index of median survival time of events
vec_cox2 <- mat_cox2[ ,med_index2]

times2 <- test_data2$time
status2 <- test_data2$status

#CI BS IBS IAE ISE based on Cox model: Non-standard model input methods
Cindex_cox2 <- Cindex(Surv(times2, status2), vec_cox2)
BS_cox2 <- Brier(Surv(times2, status2), vec_cox2, distime2[med_index2])
IBS_cox2 <- IBS(Surv(times2, status2), mat_cox2, distime2)
IAE_cox2 <- IAEISE(Surv(times2, status2), mat_cox2, distime2)[1]
ISE_cox2 <- IAEISE(Surv(times2, status2), mat_cox2, distime2)[2]
Cindex_cox2
BS_cox2
IBS_cox2
IAE_cox2
ISE_cox2

#-------------------------------------------------------------------------------
#Iteration 3:
###Fit the Cox Proportional Hazards Model on Train Data and make predictions with test data

# Standardize the data (excluding the response variable)
train_data_scaled3 <- scale(train_data3[, -c(1, 2)], 
                            center = TRUE, 
                            scale = TRUE)
train_means3 <- attr(train_data_scaled3, "scaled:center") #save train means for scaling the test data
train_sds3 <- attr(train_data_scaled3, "scaled:scale") #save train sds for scaling the test data

#Convert the scaled matrix to a data frame
train_data_scaled3 <- as.data.frame(train_data_scaled3)

#Add the unscaled 'time' and 'status' variables back to the scaled data frame
train_data_scaled3 <- cbind(train_data3[, c("time", "status")], train_data_scaled3)

#Train the model with predictive features
full.cox3 <- coxph(Surv(time,status) ~ hsa.miR.101.3p + hsa.miR.183.5p + HDL_TG, data = train_data_scaled3, x=TRUE)
summary(full.cox3)

#For each covariate, the function cox.zph() correlates the corresponding set of scaled Schoenfeld residuals with time, to test for independence between residuals and time. 
#Additionally, it performs a global test for the model as a whole.
test.ph3 <- cox.zph(full.cox3)
test.ph3


#Plot predictive features' impact on event hazard
ggforest(full.cox3, data = train_data_scaled3)


#Make model predictions with test data
test_data3 <- read.csv("test_set_3.csv")
test_data3 <- test_data3[,-(1:5)]

test_data3 <- test_data3 %>%
  dplyr::rename("status" = "Cancer_during_surveillance")
test_data3 <- test_data3 %>%
  dplyr::rename("time" = "Survival_time_years")

test_data3[, -c(1, 2)] <- scale(test_data3[, -c(1, 2)], 
                                center = train_means3, 
                                scale = train_sds3)


# Extract event times and event indicators
event_times3 <- full.cox3$y[, "time"]
event_indicator3 <- full.cox3$y[, "status"]  # This indicates if the event occurred (1) or was censored (0)
distime3 <- event_times3[event_indicator3 == 1]
distime3 <- sort(unique(distime3)) # Make sure distime is sorted
distime3
mat_cox3 <- predictSurvProb(full.cox3, test_data3, distime3) #get the survival probability matrix
med_index3 <- median(1:length(distime3)) #the index of median survival time of events
vec_cox3 <- mat_cox3[ ,med_index3]

times3 <- test_data3$time
status3 <- test_data3$status

#CI BS IBS IAE ISE based on Cox model: Non-standard model input methods
Cindex_cox3 <- Cindex(Surv(times3, status3), vec_cox3)
BS_cox3 <- Brier(Surv(times3, status3), vec_cox3, distime3[med_index3])
IBS_cox3 <- IBS(Surv(times3, status3), mat_cox3, distime3)
IAE_cox3 <- IAEISE(Surv(times3, status3), mat_cox3, distime3)[1]
ISE_cox3 <- IAEISE(Surv(times3, status3), mat_cox3, distime3)[2]
Cindex_cox3
BS_cox3
IBS_cox3
IAE_cox3
ISE_cox3

#-------------------------------------------------------------------------------
#Iteration 4:
###Fit the Cox Proportional Hazards Model on Train Data and make predictions with test data

# Standardize the data (excluding the response variable)
train_data_scaled4 <- scale(train_data4[, -c(1, 2)], 
                            center = TRUE, 
                            scale = TRUE)
train_means4 <- attr(train_data_scaled4, "scaled:center") #save train means for scaling the test data
train_sds4 <- attr(train_data_scaled4, "scaled:scale") #save train sds for scaling the test data

#Convert the scaled matrix to a data frame
train_data_scaled4 <- as.data.frame(train_data_scaled4)

#Add the unscaled 'time' and 'status' variables back to the scaled data frame
train_data_scaled4 <- cbind(train_data4[, c("time", "status")], train_data_scaled4)

#Train the model with predictive features

full.cox4 <- coxph(Surv(time,status) ~ hsa.miR.101.3p + hsa.miR.183.5p + HDL_TG, data = train_data_scaled4, x=TRUE)
summary(full.cox4)


#For each covariate, the function cox.zph() correlates the corresponding set of scaled Schoenfeld residuals with time, to test for independence between residuals and time. 
#Additionally, it performs a global test for the model as a whole.
test.ph4 <- cox.zph(full.cox4)
test.ph4


#Plot predictive features' impact on event hazard
ggforest(full.cox4, data = train_data_scaled4)


#Make model predictions with test data
test_data4 <- read.csv("test_set_4.csv")
test_data4 <- test_data4[,-(1:5)]

test_data4 <- test_data4 %>%
  dplyr::rename("status" = "Cancer_during_surveillance")
test_data4 <- test_data4 %>%
  dplyr::rename("time" = "Survival_time_years")

test_data4[, -c(1, 2)] <- scale(test_data4[, -c(1, 2)], 
                                center = train_means4, 
                                scale = train_sds4)


# Extract event times and event indicators
event_times4 <- full.cox4$y[, "time"]
event_indicator4 <- full.cox4$y[, "status"]  # This indicates if the event occurred (1) or was censored (0)
distime4 <- event_times4[event_indicator4 == 1]
distime4 <- sort(unique(distime4)) # Make sure distime is sorted
distime4
mat_cox4 <- predictSurvProb(full.cox4, test_data4, distime4) #get the survival probability matrix
med_index4 <- median(1:length(distime4)) #the index of median survival time of events
vec_cox4 <- mat_cox4[ ,med_index4]

times4 <- test_data4$time
status4 <- test_data4$status

#CI BS IBS IAE ISE based on Cox model: Non-standard model input methods
Cindex_cox4 <- Cindex(Surv(times4, status4), vec_cox4)
BS_cox4 <- Brier(Surv(times4, status4), vec_cox4, distime4[med_index4])
IBS_cox4 <- IBS(Surv(times4, status4), mat_cox4, distime4)
IAE_cox4 <- IAEISE(Surv(times4, status4), mat_cox4, distime4)[1]
ISE_cox4 <- IAEISE(Surv(times4, status4), mat_cox4, distime4)[2]
Cindex_cox4
BS_cox4
IBS_cox4
IAE_cox4
ISE_cox4

#-------------------------------------------------------------------------------
##Step 7. Plot ROC curves

#install.packages("ROCR")
library(ROCR)

# Convert survival probabilities to event probabilities (1 - survival probability)
event_probs <- 1 - vec_cox  # vec_cox contains survival probabilities, so we want event probabilities
event_probs1 <- 1 - vec_cox1
event_probs2 <- 1 - vec_cox2
event_probs3 <- 1 - vec_cox3
event_probs4 <- 1 - vec_cox4

# Create prediction object with predicted probabilities (event_probs) and actual outcomes (status)
pred <- prediction(event_probs, status)
pred1 <- prediction(event_probs1, status1)
pred2 <- prediction(event_probs2, status2)
pred3 <- prediction(event_probs3, status3)
pred4 <- prediction(event_probs4, status4)

# Calculate performance (True Positive Rate and False Positive Rate)
perf <- performance(pred, "tpr", "fpr")
perf1 <- performance(pred1, "tpr", "fpr")
perf2 <- performance(pred2, "tpr", "fpr")
perf3 <- performance(pred3, "tpr", "fpr")
perf4 <- performance(pred4, "tpr", "fpr")

dev.off()

# Plot the ROC curves
plot(perf, col = "black", lwd = 2, main = "ROC Curve for Each Iteration")
plot(perf1, col = "blue", add = TRUE, lwd = 2)
plot(perf2, col = "red", add = TRUE, lwd = 2)
plot(perf3, col = "green", add = TRUE, lwd = 2)
plot(perf4, col = "purple", add = TRUE, lwd = 2)

# Add a diagonal reference line
abline(a = 0, b = 1, col = "grey", lty = 2)

# Add a legend
legend("bottomright", legend = c("Iteration 0", "Iteration 1", "Iteration 2", "Iteration 3", "Iteration 4"),
       col = c("black", "blue", "red", "green", "purple"), lty = 1, lwd = 2)

#-------------------------------------------------------------------------------
