library(tools)

validateinput <- function(infile){
  isvalid <- "File validation complete."
  
  #a csv file should be uploaded
  if (!file_ext(inFile$name) %in% c('csv')) {
    isvalid <- "Wrong file format! Please try again with a CSV file."
    break
  }
  #read the csv file for other checks
  in.dat <- read.csv(infile)
  #list of expeected column names
  cols.expected = c("plate_id","well_id","site_id", "trt","total","count" )
  #check if the file has six columns
  if(ncol(in.dat) != 6){
    isvalid <- paste0("Input file has ", 
                      ncol(in.dat), 
                      " columns. Input file header row must have 6 columns",
                      " as specified in the instructions.")
    break
  } else if(!((setequal(colnames(in.dat),cols.expected)))){
    isvalid <- paste0("Input file header does not contain column(s)",
                      setdiff(cols.expected,colnames(in.dat)),".", collapse=",")
    break
  }
  
  return(isvalid)
}