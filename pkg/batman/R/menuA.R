menuA<-function (choices, stInd, showLine) 
{
  ## menue function, for internal use
  nc <- length(choices)
  cat(showLine)
  if (stInd == 1) {
    op <- paste(format(seq_len(nc)), ": ", choices, sep = "")
    st <-stInd
    ed<-nc+stInd-1
  } else if (stInd == 0) {
    op <- paste(format(seq_len(nc)-1), ": ", choices, sep = "")
    st <-stInd
    ed<-nc+stInd-1
  }
  cat("  ", op, "", sep = "\n")
  repeat {
    ind <- .Internal(menu(as.character(choices)))
    if (ind >= st && ind <= ed) 
      return(ind)
  }
}