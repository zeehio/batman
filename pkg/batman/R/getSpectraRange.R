getSpectraRange <- function(myVector)
{
  ranges <- strsplit(myVector[[1]][2],",")
  sno <- NULL
  for (rgs in 1:length(ranges[[1]]))
  {
    ranges2 <- strsplit(ranges[[1]][rgs], "-")
    ranges3 <-as.numeric(ranges2[[1]])
    if (length(ranges3)>1)
    {
      if (ranges3[1]>ranges3[2])
      {
        tmp1 <- ranges3[1]
        ranges3[1] <- ranges3[2]
        ranges3[2] <- tmp1
      }
    } else {
      ranges3[2] <- ranges3[1]
    }
    sno <- c(sno, ranges3[1]:ranges3[2])
  }
  
  return (sno)
}