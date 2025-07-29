##################################################################################
# Cleaning Empirical Data
# Code to accompany Wysocki, van Bork, Cramer, & Rhemtulla 
# "Cross-Lagged Network Models"
##################################################################################

# install and load required package ##############################################
#install.packages("qgraph") #uncomment to install packages for the 1st time
#install.packages("glmnet") 
#install.packages("lavaan")

library(qgraph)
library(glmnet)
library(lavaan)

####################################################################################

# Load and clean data ###########################################
# The example in the paper uses variables from the Iowa Youth and
# Familes Study, which can be accessed at 
# http://www.icpsr.umich.edu/icpsrweb/ICPSR/studies/26721
# To replicate the analyses, download and unzip the files under: 
#   "DS3:  Wave A 1989: Child Survey Data"
#   "DS9:  Wave B 1990: Child Survey Data"

setwd("~/Ecole/M1/Stage/Internship_repo/Network/CLPN Example/Empirical Example Code") #set working directory to location of data files
W1child <- read.table(file = "WaveAChild/26721-0003-Data.tsv", sep = "\t", header = TRUE)
W2child <- read.table(file = "WaveBChild/26721-0009-Data.tsv", sep = "\t", header = TRUE)
W3child <- read.table(file = "WaveCChild/26721-0016-Data.tsv", sep = "\t", header = TRUE)
W4child <- read.table(file = "WaveDChild/26721-0022-Data.tsv", sep = "\t", header = TRUE)

W1W2child<- merge(W1child, W2child, by="CASEID", all.x = TRUE, all.y = TRUE)
W123child<- merge(W1W2child, W3child, by = "CASEID", all.x = TRUE, all.y = TRUE)
W1234child<-merge(W123child, W4child, by = "CASEID", all.x = TRUE, all.y = TRUE)
#######################################################################################

# select variables of interest (school commitment, self-esteem) #######################
#school commitment

mydata <- NULL
mydata$CASEID <- W1234child$CASEID
mydata <- data.frame(mydata)

#wave 1
mydata$SC01.1 <- 6 - W1234child$at202012 #R
mydata$SC02.1 <- W1234child$at202013
mydata$SC03.1 <- W1234child$at202014
mydata$SC04.1 <- W1234child$at202015
mydata$SC05.1 <- W1234child$at202016
mydata$SC06.1 <- 6 - W1234child$at202017 #R 
mydata$SC07.1 <- 6 - W1234child$at202018 #R

mydata$SE01.1 <- 6 - W1234child$at105008 #R
mydata$SE02.1 <- 6 - W1234child$at105009 #R
mydata$SE03.1 <- W1234child$at105010
mydata$SE04.1 <- 6 - W1234child$at105011 #R
mydata$SE05.1 <- W1234child$at105012
mydata$SE06.1 <- 6 - W1234child$at105013 #R
mydata$SE07.1 <- 6 - W1234child$at105014 #R
mydata$SE08.1 <- W1234child$at105015
mydata$SE09.1 <- W1234child$at105016
mydata$SE10.1 <- W1234child$at105017

#wave 2
mydata$SC01.2 <- 6 - W1234child$bt203012 #R
mydata$SC02.2 <- W1234child$bt203013
mydata$SC03.2 <- W1234child$bt203014
mydata$SC04.2 <- W1234child$bt203015
mydata$SC05.2 <- W1234child$bt203016
mydata$SC06.2 <- 6 - W1234child$bt203017#R 
mydata$SC07.2 <- 6 - W1234child$bt203018#R

mydata$SE01.2 <- 6 - W1234child$bt105008 #R
mydata$SE02.2 <- 6 - W1234child$bt105009 #R
mydata$SE03.2 <- W1234child$bt105010
mydata$SE04.2 <- 6 - W1234child$bt105011 #R
mydata$SE05.2 <- W1234child$bt105012
mydata$SE06.2 <- 6 - W1234child$bt105013 #R
mydata$SE07.2 <- 6 - W1234child$bt105014 #R
mydata$SE08.2 <- W1234child$bt105015
mydata$SE09.2 <- W1234child$bt105016
mydata$SE10.2 <- W1234child$bt105017

#wave 3
mydata$SC01.3 <- 6 - W1234child$ct204012 #R
mydata$SC02.3 <- W1234child$ct204013
mydata$SC03.3 <- W1234child$ct204014
mydata$SC04.3 <- W1234child$ct204015
mydata$SC05.3 <- W1234child$ct204016
mydata$SC06.3 <- 6 - W1234child$ct204017 #R 
mydata$SC07.3 <- 6 - W1234child$ct204018 #R

mydata$SE01.3 <- 6 - W1234child$ct105008 #R
mydata$SE02.3 <- 6 - W1234child$ct105009 #R
mydata$SE03.3 <- W1234child$ct105010
mydata$SE04.3 <- 6 - W1234child$ct105011 #R
mydata$SE05.3 <- W1234child$ct105012
mydata$SE06.3 <- 6 - W1234child$ct105013 #R
mydata$SE07.3 <- 6 - W1234child$ct105014 #R
mydata$SE08.3 <- W1234child$ct105015
mydata$SE09.3 <- W1234child$ct105016
mydata$SE10.3 <- W1234child$ct105017

#wave 4
mydata$SC01.4 <- 6 - W1234child$dt204012 #R
mydata$SC02.4 <- W1234child$dt204013
mydata$SC03.4 <- W1234child$dt204014
mydata$SC04.4 <- W1234child$dt204015
mydata$SC05.4 <- W1234child$dt204016
mydata$SC06.4 <- 6 - W1234child$dt204017#R 
mydata$SC07.4 <- 6 - W1234child$dt204018#R

mydata$SE01.4 <- 6 - W1234child$dt105008 #R
mydata$SE02.4 <- 6 - W1234child$dt105009 #R
mydata$SE03.4 <- W1234child$dt105010
mydata$SE04.4 <- 6 - W1234child$dt105011 #R
mydata$SE05.4 <- W1234child$dt105012
mydata$SE06.4 <- 6 - W1234child$dt105013 #R
mydata$SE07.4 <- 6 - W1234child$dt105014 #R
mydata$SE08.4 <- W1234child$dt105015
mydata$SE09.4 <- W1234child$dt105016
mydata$SE10.4 <- W1234child$dt105017

#recode missings to NA:
mydata[mydata == -3] <- NA
mydata[mydata == 9] <- NA

#remove case ID variable: 
mydata <- mydata[,-1]
names(mydata)
#######################################################################################

# standardize each variable based on mean/var across all 4 waves ######################
mydata.matrix <- as.matrix(mydata)

k <- 17
mydata.stacked <- rbind(mydata.matrix[,1:k], mydata.matrix[,(k+1):(2*k)], 
                        mydata.matrix[,(2*k+1):(3*k)], mydata.matrix[,(3*k+1):(4*k)])
sdvec <- apply(mydata.stacked, 2, sd, na.rm = TRUE)
meanvec <- apply(mydata.stacked, 2, mean, na.rm = TRUE)

sdvec.4 <- rep(sdvec, 4)
meanvec.4 <- rep(meanvec, 4)

mydata.std <- scale(mydata, center = meanvec.4, scale = sdvec.4)
# mydata[1:4, c(1, 18, 35, 52)]     #look to see if that worked
# mydata.std[1:4, c(1, 18, 35, 52)] #yes it did
#######################################################################################

# save data and legend of variable descriptions #######################################
n <- nrow(mydata) #number of observations
groups <- c(rep("Commitment to School", 9), rep("Self-Esteem", 10))

mylegend <- c("SC1: I like school a lot (R)", 
              "SC2: School bores me", 
              "SC3: I don't do well at school", 
              "SC4: I don't belong at school", 
              "SC5: Homework is a waste of time", 
              "SC6: I try hard at school (R)", 
              "SC7: I finish my homework (R)",
              "SE1: I'm a person of worth (R)",
              "SE2: I have a number of good qualities (R)",
              "SE3: I'm a failure",
              "SE4: I am able to do things as well as most people (R)",
              "SE5: I do not have much to be proud of",
              "SE6: I take a positive attitude toward myself (R)",
              "SE7: I am satisfied with myself (R)",
              "SE8: I feel useless at times",
              "SE9: I wish I could have more respect for myself",
              "SE10:I think I am no good at all")

labels1 <- c("SC1", "SC2", "SC3", "SC4", "SC5", "SC6",
             "SC7", "SE1", "SE2", "SE3", "SE4", 
             "SE5", "SE6", "SE7", "SE8", "SE9", 
             "SE10")

# save data and legends (not necessary)###########################
# Add file name and location you would like to save the file to the 
# file argument below

save(mydata, mydata.std, mylegend, labels1, k, n, file = "mydata.Rda")

#load(file = "myData.Rda") #to load the data and legends without running all of the above code
#################################################################