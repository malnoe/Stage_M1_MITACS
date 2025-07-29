# Visualize the Constrained Model Prediction Metrics
######################################################

# Load things
load("predictionResConsModel.Rdat")
load("bootstrapConsPredictSummaryStats.Rdat")

library(ggplot2)
library(dplyr)
library(gridExtra)

# Create helpers/labels

groups <- c(rep("Commitment to School", 7), rep("Self-Esteem", 10)) 


# Join prediction and bootstrap estimates
predictionEstimates$Variable <- sub("\\..*", "" ,predictionEstimates$Variable) 
predictionEstimates %>%
  group_by(PredictionIndex, Variable) %>%
  summarise(meanPredictionValue = mean(PredictionValue)) -> meanPredictionEstimates

left_join(x = meanPredictionEstimates, y = bootPredSummaryStat, 
          by = c("PredictionIndex", "Variable")) -> meanBootPredictionEstimates

meanBootPredictionEstimates$lower <- meanBootPredictionEstimates$meanPredictionValue - meanBootPredictionEstimates$s
meanBootPredictionEstimates$lower[meanBootPredictionEstimates$lower < 0] <- 0
meanBootPredictionEstimates$upper <- meanBootPredictionEstimates$meanPredictionValue + meanBootPredictionEstimates$s

# Visualize CL Prediction
inpredictViz <- meanBootPredictionEstimates %>% 
  filter(PredictionIndex == "CrossLagInPred") %>%
  ggplot(aes(x = Variable, y = meanPredictionValue)) + 
  geom_bar(stat = "identity", color = "black", 
           fill = c(rep("#E69F00", 7), rep("#56B4E9", 10))) +
  geom_errorbar(aes(ymin = lower, ymax = upper), 
                width = .5) + 
  coord_flip() + 
  scale_x_discrete(limits = rev) + 
  labs(x = "in-prediction", y = "") + 
  theme_bw() + 
  theme(text = element_text(size=20))

outpredictViz <- meanBootPredictionEstimates %>% 
  filter(PredictionIndex == "CrossLagOutPred") %>%
  ggplot(aes(x = Variable, y = meanPredictionValue)) + 
  geom_bar(stat = "identity", color = "black", 
           fill = c(rep("#E69F00", 7), rep("#56B4E9", 10))) +
  geom_errorbar(aes(ymin =lower, ymax = upper), 
                width = .5) + 
  coord_flip() + 
  scale_x_discrete(limits = rev) +
  labs(x = "out-prediction", y = "") +
  theme_bw()+ 
  theme(text = element_text(size=20))

jpeg("InOutPrediction.jpg", width=8.6, height=8, res=800, units="in")
grid.arrange(inpredictViz, outpredictViz, nrow  = 2 )
dev.off() 

# Within Construt

withinCinpredictViz <- meanBootPredictionEstimates %>% 
  filter(PredictionIndex == "SameConstructInPred") %>%
  ggplot(aes(x = Variable, y = meanPredictionValue)) + 
  geom_bar(stat = "identity", color = "black", 
           fill = c(rep("#E69F00", 7), rep("#56B4E9", 10))) +
  geom_errorbar(aes(ymin =lower, ymax = upper), 
                width = .5) + 
  coord_flip() + 
  scale_x_discrete(limits = rev) +
  labs(x = "in-prediction", y = "") +
  theme_bw()+ 
  theme(text = element_text(size=20))

withinCoutpredictViz <- meanBootPredictionEstimates %>% 
  filter(PredictionIndex == "SameConstructOutPred") %>%
  ggplot(aes(x = Variable, y = meanPredictionValue)) + 
  geom_bar(stat = "identity", color = "black", 
           fill = c(rep("#E69F00", 7), rep("#56B4E9", 10))) +
  geom_errorbar(aes(ymin =lower, ymax = upper), 
                width = .5) + 
  coord_flip() + 
  scale_x_discrete(limits = rev) +
  labs(x = "out-prediction", y = "") +
  theme_bw()+ 
  theme(text = element_text(size=20))

jpeg("WithinConstructPrediction.jpg", width=8.6, height=8, res=800, units="in")
grid.arrange(withinCinpredictViz, withinCoutpredictViz, nrow  = 2 )
dev.off() 

# Between Constructs
otherCinpredictViz <- meanBootPredictionEstimates %>% 
  filter(PredictionIndex == "OtherConstructInPred") %>%
  ggplot(aes(x = Variable, y = meanPredictionValue)) + 
  geom_bar(stat = "identity", color = "black", 
           fill = c(rep("#E69F00", 7), rep("#56B4E9", 10))) +
  geom_errorbar(aes(ymin =lower, ymax = upper), 
                width = .5) + 
  coord_flip() + 
  scale_x_discrete(limits = rev) +
  labs(x = "in-prediction", y = "") +
  theme_bw()+ 
  theme(text = element_text(size=20))

otherCoutpredictViz <- meanBootPredictionEstimates %>% 
  filter(PredictionIndex == "OtherConstructOutPred") %>%
  ggplot(aes(x = Variable, y = meanPredictionValue)) + 
  geom_bar(stat = "identity", color = "black", 
           fill = c(rep("#E69F00", 7), rep("#56B4E9", 10))) +
  geom_errorbar(aes(ymin =lower, ymax = upper), 
                width = .5) + 
  coord_flip() + 
  scale_x_discrete(limits = rev) +
  labs(x = "out-prediction", y = "") +
  theme_bw()+ 
  theme(text = element_text(size=20))

jpeg("OtherConstructPrediction.jpg", width=8.6, height=8, res=800, units="in")
grid.arrange(otherCinpredictViz, otherCoutpredictViz, nrow  = 2 )
dev.off() 
