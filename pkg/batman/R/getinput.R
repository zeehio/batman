getinput<-function(lowlim=0,highlim=1)
{
  ## written by Dr. Jie Hao, Imperial College London
  ## get input from user, for internal use
  repeat {
    inp<- scan(n=1, quiet = TRUE)
    if (length(inp)) {
      if (highlim>0) {
        if (inp >= lowlim && inp <= highlim ) 
          return(inp)
      } else {
        if (inp >= lowlim )
          return(inp)
      }
    }  
  }
}