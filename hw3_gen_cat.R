############Call library################################
library(tidyverse)
library(stringi)
library(data.table)
library(skimr)
library(energy)
library(fda.usc)

#############Getting the data###########################
# load the data file
load("my_data.RData")


#Assign data
df <- my_data$avg_dat

#############Function definition########################

#################Pearson method##############
n_top_corr_gene_p <- function(dat_d){
  
  
  # initialize variable
  cor_val <- c()
  first_gene <- c()
  other_gene <- c()
  
  # get the rows
  n <- nrow(dat_d)
  
  for (i in 1:(n-1)){
    
    xk<- i+1
    
    for(j in xk:n){
      
      # take first x variable
      x<- dat_d[i,]
      nam <- as.character(x$gene)
      x<- x[,-1]
      x<- as.numeric(x)
      
      # take second y variable
      y<- dat_d[j,]
      nam_2 <- as.character(y$gene)
      y<- y[,-1]
      y<- as.numeric(y)
      
      # do correlation test
      test <- cor.test(x,y,method = "pearson")
      
      # take the estimate and p.value
      c1 <- test$estimate
      
      # absolute c1
      c1 <- abs(c1)
      
      
      # put the value in the vector
      cor_val <- c(cor_val,c1)
      first_gene <- c(first_gene,nam)
      other_gene <- c(other_gene,nam_2)
      
    }
  }
  
  # selecting data
  new_d <- cbind(first_gene, other_gene, cor_val) %>% 
    as_tibble() %>%
    arrange(.,desc(cor_val)) %>%
    filter(first_gene != other_gene) %>%
    slice(1:20)
  
  # Write the file to csv
  write.csv(new_d,"gene_pearson.csv")
  print("Please open the gene_pearson.csv file in your current working directory")
  
  
}

################Corr Distance method###################

n_top_corr_gene_d <- function(dat_d){
  
  # initialize variable
  cor_val <- c()
  first_gene <- c()
  other_gene <- c()
  
  # get the rows
  n <- nrow(dat_d)
  
  
  for (i in 1:(n-1)){
    
    xk<- i+1
    
    for(j in xk:n){
      
      # take first x variable
      x<- dat_d[i,]
      nam <- as.character(x$gene)
      x<- x[,-1]
      x<- as.numeric(x)
      
      # take second y variable
      y<- dat_d[j,]
      nam_2 <- as.character(y$gene)
      y<- y[,-1]
      y<- as.numeric(y)
      
      # do correlation test
      test <- dcor(x,y,1.5)
      
      # absolute correlation test
      test <- abs(test)
      
      # put the value in the vector
      cor_val <- c(cor_val,test)
      first_gene <- c(first_gene,nam)
      other_gene <- c(other_gene,nam_2)
      
    }
  }
  
  # selecting data
  new_d <- cbind(first_gene, other_gene, cor_val) %>% 
    as_tibble() %>%
    arrange(.,desc(cor_val)) %>%
    filter(first_gene != other_gene) %>%
    slice(1:20)
  
  # Write the file to csv
  write.csv(new_d,"gene_distance.csv")
  print("Please open the gene_distance.csv file in your current working directory")
  
}

#############CALL for testing ##################

dat_d<-my_data$avg_dat %>% slice(1:50)
start_time <- Sys.time()
n_top_corr_gene_p(dat_d)
end_time <- Sys.time()
Total_time <- end_time - start_time
Total_time


dat_d<-my_data$avg_dat %>% slice(1:50)
start_time <- Sys.time()
n_top_corr_gene_d(dat_d)
end_time <- Sys.time()
Total_time <- end_time - start_time
Total_time
