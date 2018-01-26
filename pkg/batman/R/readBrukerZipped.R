readBrukerZipped<-function(BrukerDataZipDir)
{
  ## written by Dr. Jie Hao, Imperial College London
  warnDef<-options("warn")$warn
  warnRead<-options(warn = -1)
  datapath<-BrukerDataZipDir
  unzipfile <- NULL
  ## unzip folder
  zipfile <-list.files(path = datapath, pattern = "*.zip$", all.files = FALSE,full.names = TRUE, recursive = TRUE, ignore.case = TRUE)
  if (length(zipfile)>0)
  {
    for (i in 1:length(zipfile))
    {
      unzip(zipfile[i], exdir = datapath)
      unzipfile[i] <- substr(zipfile[i], 1, nchar(zipfile[1])-4)
    }
  } else if (substr(datapath,nchar(datapath)-3,nchar(datapath)) == '.zip')
  {
    zipfile <- datapath
    ns <- strsplit(zipfile, "/")
    ns <- ns[[1]]
    datapath <- substr(zipfile, 1, nchar(zipfile)-nchar(ns[length(ns)])-1)
    unzip(zipfile, exdir = datapath)
    unzipfile <- substr(zipfile, 1, nchar(zipfile)-4)
  }
  
  sa<-readBruker(datapath)
  if (length(unzipfile)>0)
  {
    unlink(unzipfile, recursive = TRUE, force = FALSE)
  }
  warnRead<-options(warn = warnDef)
  return (sa)
}