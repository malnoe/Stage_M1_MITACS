library(ggplot2)
library(dplyr)
library(gridExtra)
set.seed(13)

# Generate synthetic data
n <- 100
adversity <- runif(n, 0, 8)
noise <- rnorm(n, 0, 2)

# Schema A: Higher = better functioning (↓ with adversity)
outcome_A <- -0.8 * adversity + 8 + noise
df_A <- data.frame(adversity, outcome = outcome_A)

# Schema B: Higher = worse functioning (↑ with adversity)
outcome_B <- 0.8 * adversity + noise
df_B <- data.frame(adversity, outcome = outcome_B)

# Fit models
model_A <- lm(outcome ~ adversity, data = df_A)
model_B <- lm(outcome ~ adversity, data = df_B)

# Plot function
make_plot <- function(df, model, title, is_high_better) {
  x_range <- range(df$adversity)
  y_range <- range(df$outcome)
  
  plot <- ggplot(df, aes(x = adversity, y = outcome)) +
    geom_point(color = "grey50", size = 2) +
    geom_smooth(method = "lm", se = FALSE, color = "grey50", size=0.8) +
    labs(title = title, x = "Adversity", y = "Outcome") +
    theme_minimal(base_size = 12)+
    theme(
      axis.text.x = element_blank(),
      axis.text.y = element_blank()
    )

  if(is_high_better){
    plot <- plot + geom_text(x=max(na.omit(df$adversity))-2.5,y=max(na.omit(df$outcome))-2.5,label="Resilient",alpha=0.2,color="grey50",size=5) + geom_text(x=min(na.omit(df$adversity))+2.5,y=min(na.omit(df$outcome))+2.5,label="Vulnerable",alpha=0.2,color="grey50",size=5)
  }
  else{
    plot <-  plot + geom_text(x=min(na.omit(df$adversity))+2.5,y=max(na.omit(df$outcome))-2.5,label="Vulnerable",alpha=0.2,color="grey50",size=5) + geom_text(x=max(na.omit(df$adversity))-2.5,y=min(na.omit(df$outcome))+2.5,label="Resilient",alpha=0.2,color="grey50",size=5)
  }
  
}

# Create both plots
plot_A <- make_plot(df_A, model_A, "Diagram A: higher score = better functioning", is_high_better = TRUE)

plot_B <- make_plot(df_B, model_B,
                    "Diagram B: higher score = worse functioning", is_high_better = FALSE)

grid.arrange(plot_A, plot_B, ncol = 2)