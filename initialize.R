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
                header = TRUE, stringsAsFactors = FALSE)
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


# write to all cells
write_csv(samples, "all.csv")

# write to single cells
singlet <- samples %>% filter(p1 & sF & sS)
write_csv(singlet, "singlet.csv")

# write live singlet
liveSinglet <- samples %>%
  filter(p1 & sF & sS & live)
write_csv(liveSinglet, "liveSinglet.csv")

