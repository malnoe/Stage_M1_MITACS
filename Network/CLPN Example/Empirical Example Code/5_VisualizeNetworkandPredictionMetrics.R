##################################################################################
# Visualizing a CLPN Network 
# Code to accompany Wysocki, van Bork, Cramer & Rhemtulla
# "Cross-Lagged Network Models"
##################################################################################


### Set-Up #####################################################################


# Install and load required package(s) ########
#install.packages("lavaan") #uncomment to install packages for the 1st time
library(lavaan)
library(ggplot2)
library(qgraph)

# Load Model ###########################

load(file = "hybridModel.Rdata") 

load(file = "EmpiricalExampleData.Rda") # Loads legend

# Subset hybrid CLPN results
hybridBList <- hybridModelList$hybridBMatrixSep



### Visualizing a CLPN #########################################################

# create a 'group' vector that says which variables belong to which construct:
groups <- c(rep("Commitment to School", 7), rep("Self-Esteem", 10)) 

# create a vector of variable labels for the figure
labels <- c("SC1","SC2","SC3","SC4","SC5", "SC6", "SC7",
            "SE1","SE2","SE3","SE4","SE5","SE6","SE7", "SE8", 
            "SE9", "SE10")

# Make sure edge thickness is consistent across plots
MaxEdge <- max(unlist(hybridBList))


#plot CLPN results for each wave: 
jpeg("CLPN_hybrid.jpg", width=12, height=12, res=800, units="in")  #(uncomment this code to save plot as jpeg)

layout1_hybrid <- qgraph(hybridBList[[1]], legend = FALSE)$layout
par(mar = rep(0, 4))
layout(mat = matrix(c(1, 2, 3, 4, 4, 4), 3, 2), widths = c(1.2, 1)) #to plot multiple graphs in one figure
#T1 -> T2
qgraph(hybridBList[[1]], groups = groups, labels = labels, legend = FALSE,
       color = c("#E69F00", "#56B4E9"), asize = 5, layout = layout1_hybrid,
       maximum = MaxEdge, title = expression(7^{th}~to~8^{th}~grade), title.cex = 2, 
       edge.width = 2)
#T2 -> T3
qgraph(hybridBList[[2]], groups = groups, labels = labels, legend = FALSE,
       color = c("#E69F00", "#56B4E9"), asize = 5, layout = layout1_hybrid,
       maximum = MaxEdge, title = expression(8^{th}~to~9^{th}~grade), title.cex = 2, 
       edge.width = 2)
#T3 -> T4
qgraph(hybridBList[[3]], groups = groups, labels = labels, legend = FALSE,
       color = c("#E69F00", "#56B4E9"), asize = 5, layout = layout1_hybrid,
       maximum = MaxEdge, title = expression(9^{th}~to~10^{th}~grade), title.cex = 2, 
       edge.width = 2)
#legend
plot.new()
par(mar = rep(0, 4))
legend(x = "center", inset = c(0, 0), bty = "n", legend = mylegend,  
       col = c(rep("#E69F00",7),rep("#56B4E9",10)), y.intersp = 1,
       pch = 19, cex = 1.8) #add legend to identify variables
dev.off() #(uncomment this code to save plot as jpeg)
