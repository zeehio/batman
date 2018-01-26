menuA<-function (choices, stInd, showLine) 
{
  ## written by Dr. Jie Hao, Imperial College London
  ## menue function, for internal use
  nc <- length(choices)
  cat(showLine)
  cat(as.character(choices))
  
  ind <- 0
  while (ind < 1 || ind > 3){
   ind <- readline("\nSelection: ")
   #ind<- scan(n=1, quiet = TRUE)
   #ind <- as.integer(ind)
    ind <- ifelse(grepl("\\D",ind),-1,as.integer(ind))
    #if(is.na(ind)){break}  # breaks when hit enter
  }
  return(ind)

}