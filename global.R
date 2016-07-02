# Some of the reversals in this function
# will be difficult to follow, but are necessary
# 6/30/2016
subsetData <- function(data,ranges)
{
  CG3 <- data
  nR <- nrow(CG3)
  CG3 <- CG3[rev(rownames(CG3)),]
  CG3$V7<-rev(CG3$V7)
  
  ranges <- strsplit(ranges,",")[[1]]

  genRows <- function(range)
  {
    if (length(strsplit(range,"-")[[1]])>1)
    {
      rows <- strsplit(range,"-")[[1]]
      start <- as.numeric(rows[1])
      end <- as.numeric(rows[2])
      return(CG3[end:start,]) # end:start is strange, but necessary
      #return(CG3[(nR-end):(nR-start),]) 
      
      # BECAUSE THE ROW LABELS ARE REVERSED!!!
      # It really just selects the rows you want.
    }
    else
    {
      row <- as.numeric(range[[1]][1])
      return(CG3[row,])
    }
  }
  CG5 <- data.frame()
  for (r in rev(ranges))
  {
    CG5 <- rbind(CG5,genRows(r))
  }
  #CG5 <- CG5[rev(rownames(CG5)),]
  CG5$V7<-rev(CG5$V7)
  CG3<-CG5
  return(CG3)
}