getDistance <- function(matrix, rownames,colnames){
  mat <- matrix(nrow = length(rownames), ncol = length(colnames),
                dimnames = list(rownames,colnames))
  cat("create new matrix..",'\n')
  #c1 = 0
  for(i in rownames){ #row
    #c1=c1+1
    #c2 = c1
    for(j in colnames){ #every col
      #c2 = c2+1
      mat[i,j] <- sum(matrix[,i] & matrix[,j])/sum(matrix[,i] | matrix[,j])
      #mat[i,j] <- fisher.test(tab)
      #cat('\r', c1, ' ', c2)
      #flush.console() 
      
    }
  }
  return(mat)
} 

#######################################
data <- readRDS("../Desktop/stud/7sem/kursinis2/TF_cg_study.rds")

library("taRifx")
names <- colnames(data) # [1:10]
Sys.time()
study_dist <- getDistance(data, names, names)
Sys.time()
#saveRDS(study_dist, "../Desktop/stud/7sem/kursinis2/study_distance.rds")
Sys.time()
