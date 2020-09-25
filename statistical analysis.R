library(tidyverse)
library(conicfit)
library(sf)
rm(list = ls())

setwd("C:/Users/cc/Desktop/Examples/csv")

a0 <- read.csv("liveSinglet.csv",
               header = TRUE,
               stringsAsFactors = FALSE)

##################################Overview of a0
a0 %>%
  group_by(treatment, stain) %>%
  summarise(
    count = n()
  ) %>% view()

################################################
# calculate TPE increment in control and experimental groups
TPEList <- a0 %>%
  group_by(treatment, stain, replicate) %>%
  summarize(
    medianTPE = median(V450.A),
  ) %>%
  pivot_wider(names_from = stain, values_from = medianTPE) %>%
  mutate(TPE = `TPE-NMI` - DMSO) %>%
  select(c(-`TPE-NMI`, -DMSO))

# Statistical analysis of two inputs
t <- function(d0, d1){
  
  s0 <- d0 %>%
    group_by(treatment) %>%
    summarise(
      weight = mean(TPE),
      sd = sd(TPE),
      n = n(),
      se = sd / sqrt(n)
    ) %>%
    mutate(p = c(1.0))
  
  s1 <- d1 %>%
    group_by(treatment) %>%
    summarise(
      weight = mean(TPE),
      sd = sd(TPE),
      n = n(),
      se = sd / sqrt(n)
    ) %>%
    mutate(p = c(1.0))
  
  v <- var.test(d0$TPE, d1$TPE)
  
  if(v$p.value <= .05) {
    p <- t.test(d0$TPE, d1$TPE, var.equal = FALSE)$p.value
  } else {
    p <- t.test(d0$TPE, d1$TPE, var.equal = TRUE)$p.value
  }
  
  if(s0$weight > s1$weight){
    s0$p <- p
  }else{
    s1$p <- p
  }
  
  rbind(s0, s1) %>%
    mutate(sig = cut(p, 
                     breaks = c(0, .0001, .001, .01, .05, Inf),
                     labels = c("****",
                                "***",
                                "**",
                                "*",
                                ""),
                     right = FALSE)
    )
}

####################  compare control and Tunicamycin

result <- t(TPEList[TPEList$treatment == "Untreated",],
            TPEList[TPEList$treatment == "Tunicamycin",])

#################### plot
result$treatment <- factor(result$treatment,
                           levels = c("Untreated","Tunicamycin"))

ggplot(result, aes(treatment, weight)) +
  geom_errorbar(aes(ymin = weight - se,
                    ymax = weight + se), 
                width =.05,
                alpha = .5,
                size = .25,
                colour = "black") +
  geom_point(alpha = .5,
             size = 1.5,
             shape = 16) +
  scale_y_continuous(name = "Normalized TPE intensity (a.u.)",
                     limits = c(0, 800)) +
  scale_x_discrete(name = NULL) +
  theme_classic() +
  theme(
    strip.background = element_rect(colour = "white"),
    axis.text.x = element_text(
      angle = 45,
      hjust =1
    ),
    axis.text.y = element_text(
      margin = margin(t = 0, r = 5, b = 0, l = 5)
    )) +
  geom_text(aes(y = weight + se + 20,
                label = sig),
            colour = "black")
