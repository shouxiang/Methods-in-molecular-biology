##########################  extract files

library(tidyverse)
library(conicfit)
library(sf)

rm(list = ls())

process <- function(df) {
  
  isInGate <- function(df, gatePoints, x, y) {
    
    gatePoints[nrow(gatePoints), ] <- gatePoints[1, ]
    
    gateSfc <- gatePoints %>% 
      list() %>% 
      st_polygon() %>% 
      st_sfc()
    
    xySfc <- cbind(df[[x]], df[[y]]) %>% 
      st_multipoint() %>%
      st_sfc() %>%
      st_cast("POINT")
    
    st_intersects(xySfc, gateSfc, sparse =  FALSE)[, 1]
  }
  
  #p1 the main population
  p1Points <-calculateEllipse(x = 90000, y = 60000,
                              a = 60000, b= 45000,
                              angle = 60, steps = 10)
  p1 <- isInGate(df, p1Points,
                 "FSC.A", "SSC.A")
  
  # to generate a rectangular polygon points 1e+03
  recRange <- function(xmin, xmax, ymin, ymax) {
    
    rbind(c(xmin * 1e+03,ymin * 1e+03),
          c(xmin * 1e+03,ymax * 1e+03),
          c(xmax * 1e+03,ymax * 1e+03),
          c(xmax * 1e+03,ymin * 1e+03),
          c(xmin * 1e+03,ymin * 1e+03))
  }
  
  #Singlet from FSC
  sFPoints <- recRange(75, 150, 10, 70)
  sFPoints <- 
    sF <- isInGate(df, sFPoints, 
                   "FSC.W", "FSC.H")
  
  #Singlet from SSC
  sSPoints <- recRange(95, 200, 1, 55)
  sS <- isInGate(df, sSPoints, "SSC.W", "SSC.H")
  
  #live from Singlet
  livePoints <- rbind(c(4e04, 5),
                      c(4e04, 1000),
                      c(14e04, 1000),
                      c(14e04, 5),
                      c(4e04, 5))
  live <- isInGate(df, livePoints, "FSC.A", "APC.A")
  
  df %>% mutate(p1 = p1,
                sF = sF,
                sS = sS,
                live = live)
}

readFile <- function(fileName) {
  
  file <- read.csv(fileName,
                   header = TRUE, stringsAsFactors = FALSE,
                   nrow = 1000)
  file %>% select(FSC.H, FSC.W, FSC.A, 
                  SSC.H, SSC.W, SSC.A,
                  V450.A, APC.A) %>%
    mutate("sampleName" = str_remove_all(fileName,
                                         "export_Specimen_001_|.csv")) %>%
    separate(sampleName, into = c("treatment", "stain","replicate"),
             sep = "_", convert = TRUE)
}

setwd("C:/Users/cc/Desktop/Examples/csv")

fileName <- list.files(pattern = "^export_Specimen_001_")
samples <- data.frame()

for (i in 1:length(fileName)) {
  samples <- samples %>% rbind(readFile(fileName[i]))
}

samples <- process(samples)

################################P1 in all cells
df <- samples

summary <- df %>% 
  group_by(treatment, stain) %>% 
  summarise(
    
    `P1 #` = sum(p1),
    `% of the parent` = mean(p1) * 100,
    `parent #` = n()
  )

p1Points <-calculateEllipse(x = 90000, y = 60000,
                            a = 60000, b= 45000,
                            angle = 60, steps = 10) %>% as.data.frame()

ggplot(df, aes(FSC.A, SSC.A)) +
  geom_point(alpha = 1/15,
             size = 1.5) +
  scale_x_continuous(limits = c(0 , 17.5* 1e04),
                     name = "FSC.A") +
  scale_y_continuous(limits = c(0 , 15* 1e04),
                     name = "SSC.A") +
  geom_path(data = p1Points, aes(V1, V2)) +
  theme_classic() +
  theme(
    strip.background = element_rect(colour = "white"),
    axis.text.x = element_text(
      angle = 30,
      hjust = 1,
      margin = margin(t = 5, r = 0, b = 5, l = 0)),
    axis.text.y = element_text(
      margin = margin(t = 0, r = 5, b = 0, l = 5))
  )
#########################################################################

##########################  see FSC single in P1
df <- samples[samples$p1 ,]

summary <- df %>% 
  group_by(treatment, stain) %>% 
  summarise(
    
    `sF #` = sum(sF),
    `% of the parent` = mean(sF) * 100,
    `parent #` = n()
  )

recRange <- function(xmin, xmax, ymin, ymax) {
  
  rbind(c(xmin * 1e+03,ymin * 1e+03),
        c(xmin * 1e+03,ymax * 1e+03),
        c(xmax * 1e+03,ymax * 1e+03),
        c(xmax * 1e+03,ymin * 1e+03),
        c(xmin * 1e+03,ymin * 1e+03))
}

sFPoints <- recRange(75, 150, 10, 70) %>% as.data.frame()

ggplot(df, aes(FSC.W, FSC.H)) +
  geom_point(alpha = 1/20, size = 1.5) +
  scale_x_continuous(limits = c(5e+04 , 22e+04)) +
  scale_y_continuous(limits = c(0e+04 , 7.5e+04)) +
  geom_path(data = sFPoints, aes(V1, V2)) +
  theme_classic() +
  theme(
    strip.background = element_rect(colour = "white"),
    axis.text.x = element_text(
      angle = 30,
      hjust = 1,
      margin = margin(t = 5, r = 0, b = 5, l = 0)),
    axis.text.y = element_text(
      margin = margin(t = 0, r = 5, b = 0, l = 5))
  )

######################################################################

################################### see SSC single in FSC
df <- samples[samples$p1 & 
                samples$sF,]

summary <- df %>%
  group_by(treatment, stain) %>% 
  summarise(
    
    `sS #` = sum(sS),
    `% of the parent` = mean(sS) * 100,
    `parent #` = n()
  )

recRange <- function(xmin, xmax, ymin, ymax) {
  
  rbind(c(xmin * 1e+03,ymin * 1e+03),
        c(xmin * 1e+03,ymax * 1e+03),
        c(xmax * 1e+03,ymax * 1e+03),
        c(xmax * 1e+03,ymin * 1e+03),
        c(xmin * 1e+03,ymin * 1e+03))
}

sSPoints <- recRange(95, 200, 1, 55) %>% as.data.frame()

ggplot(df, aes(SSC.W, SSC.H)) +
  geom_point(alpha = 1/20, size = 1.5) +
  scale_x_continuous(limits = c(9e04 , 30e04)) +
  scale_y_continuous(limits = c(0e04 , 6e04)) +
  geom_path(data = sSPoints, aes(V1, V2)) +
  theme_classic() +
  theme(
    strip.background = element_rect(colour = "white"),
    axis.text.x = element_text(
      angle = 30,
      hjust = 1,
      margin = margin(t = 5, r = 0, b = 5, l = 0)),
    axis.text.y = element_text(
      margin = margin(t = 0, r = 5, b = 0, l = 5))
  )

########################################################################

##############################   see live in singlets
df <- samples[samples$p1 & 
                samples$sF &
                samples$sS,]

summary <- df %>% 
  group_by(treatment, stain) %>% 
  summarise(
    
    'live #' = sum(live),
    '% of the parent' = mean(live) * 100,
    'parent #' = n()
  )

livePoints <- rbind(c(4e04, 5),
                    c(4e04, 1000),
                    c(14e04, 1000),
                    c(14e04, 5),
                    c(4e04, 5)) %>% as.data.frame()

ggplot(df, aes(FSC.A, APC.A)) +
  geom_point(alpha = 1/15, size = 1.5) +
  scale_x_continuous(limits = c(3e04 , 15e04)) +
  scale_y_continuous(trans = "log10") +
  geom_path(data = livePoints, aes(V1, V2)) +
  theme_classic() +
  theme(
    strip.background = element_rect(colour = "white"),
    axis.text.x = element_text(
      angle = 30,
      hjust = 1,
      margin = margin(t = 5, r = 0, b = 5, l = 0)),
    axis.text.y = element_text(
      margin = margin(t = 0, r = 5, b = 0, l = 5))
  )
