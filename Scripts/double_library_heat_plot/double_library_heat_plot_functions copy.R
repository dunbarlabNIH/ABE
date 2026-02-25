######################
# getLibraryBarcodes #
######################
# INPUT: Given curData from combine file (rownames are barcodes)
# OUTPUT: Extract and return all barcodes of a particular library
getLibraryBarcodes <- function(curData, libraryID) {
  is.part.of.library = grepl(x = rownames(curData) , pattern = paste0("^" , libraryID))
  rownames(curData)[is.part.of.library]
}

#########################
# getTopLibraryBarcodes #
#########################
# INPUT: Given dataframe from combine step process
# OUTPUT: Extract top n barcodes from each sample
# @input:
  # curData: data frame from loading combine file
  # libraryID: 5-6 bp ID of given library
  # n: number of barcodes
# @output:
  # text vector of barcode names which make the list
getTopLibraryBarcodes <- function(curData, libraryID, n = 10) {
  libBarcodes <- getLibraryBarcodes(curData, libraryID)
  libSubset = curData[libBarcodes,]
  getTopBarcodes(curData = libSubset, n = n)
}

##################
# getTopBarcodes #
##################
# get top barcodes of a combine data frame
# OUTPUT: character vector of top barcode values
getTopBarcodes <- function(curData, n = 10){
  top.barcodes.by.sample <- lapply(X = colnames(curData), FUN = function(col){topNBarcodesForSample(curData, col, n = n)} )
  top.barcodes.by.sample <- unlist(top.barcodes.by.sample)
  #return non duplicated values
  top.barcodes.by.sample[!duplicated(top.barcodes.by.sample)]
}

###############
# getLibOrder #
###############
# get library order for a given library in a combine data frame
# INPUT: curData (data frame or Prop table from combine file) & libID
# OUTPUT: vector of barcode names ordered from greatest to least
getLibOrder <- function(curProp, libraryID, curMethod = 'euclidean'){
  #get only library barcodes
  lib.barcodes <- getLibraryBarcodes(curData = curProp, libraryID = libraryID)
  libSubset = curProp[lib.barcodes,]
  #order
  lib.dist <- dist(curProp[lib.barcodes,], method = curMethod)
  #replace na's with very large value
  lib.dist <- dist.replace.nas(lib.dist)
  barcode.order <- hclust(lib.dist)$order
  #return actual barcodes in order of clustered values
  rownames(curProp)[barcode.order]
}

"
orderDataByLibraries
Description:
  given data frame from combine step, order the data into clusters associate with each library
@input:
  curProp: data normalized by to proportion of sample
  libraries: list of libraries where each list element includes a library$ID object of the 6bp ID
#output:
  character vector of barcode order
"
getDataOrderByLibraries <- function( curProp, libraries ){
  order.barcodes = character()
  for(library in libraries){
    libraryOrder = getLibOrder(curProp, library$ID)
    #append to end of order.barcodes
    order.barcodes = c(order.barcodes, libraryOrder)
  }
  order.barcodes
}

# Mapping info
# getReadMeData: Read in readme data from combine ReadMe file
getReadMeData <- function(readme.file){
  temp.read.in = read.delim(file = readme.file, header = T, 
                            skip = 3, stringsAsFactors = F) # read in file (skip first 3 lines that are totals)
  last.line.i = length(temp.read.in[,1]) # get last line (also totals)
  temp.read.in = temp.read.in[ 1:(last.line.i - 1),] # skip last line (also totals)
  colnames(temp.read.in) = c("FILENAME" , "MAPPED" , "READS" , "MAP.PERCENT","THRESH", "LIBID") # rename column names ("%" is problematic)
  return(temp.read.in)
}

"
generate mapping data for annotation making
@input:
  curReadMeData: data frame generated from getReadMeData() function
@output:
  data frame of mapping data with columns of mapped reads and non-mapped reads
"
getMappingDataOneLibrary <- function( curReadMeData ){
  #first get appropriate names
  sample.names = curReadMeData$FILENAME
  mapping.data = data.frame(sample.names,curReadMeData$MAPPED, curReadMeData$READS, row.names = 1)
  colnames(mapping.data) = c("MAPPED", "OTHER")
  mapping.data$OTHER = mapping.data$OTHER - mapping.data$MAPPED
  mapping.data
}

getMappingDataOn2Libraries <- function( readmeData1, libraryInfo1 ,readmeData2, libraryInfo2 ){
  assertthat::are_equal(readmeData1$FILENAME, readmeData2$FILENAME) # must have identical rownames (ie samples) between two libraries
  sample.names <- readmeData1$FILENAME
  total.reads <- readmeData1$READS # total reads will also be the same since samples are identical
  lib1.reads <- readmeData1$MAPPED
  lib2.reads <- readmeData2$MAPPED
  other.reads = total.reads - lib1.reads - lib2.reads  # all other unmapped reads
  mapping.data <- data.frame(sample.names, lib1.reads, lib2.reads, other.reads, row.names = 1)
  lib1.name = paste0("LIB", libraryInfo1$name,"_MAPPED")
  lib2.name = paste0("LIB", libraryInfo2$name, "_MAPPED")
  colnames(mapping.data) = c(lib1.name, lib2.name, "OTHER")
  mapping.data
}

"
generateMappingAnnotation: generate the above row for the mapping annotation given mapping data and library information
@input:
  mapping.data: read in README file from combine experiment
  library.info: list of library information including color, name, ID
@output
  heatmap annotation for mapping percentage
"
generate_1_LibraryMappingAnnotation <- function( mapping.data, 
                                                 lib.info,
                                                 show.total.reads = T,
                                                 title = NULL
                                                 ){
  require(ComplexHeatmap)
  require(circlize)
  #now lets try to make a barplot out of it
  reads.based.column_ha = HeatmapAnnotation(Reads = anno_barplot(as.matrix(mapping.data), 
                                                                 axis = TRUE,
                                                                 gp = gpar( fill = c(lib.info$color,"grey") , 
                                                                            width = unit(2, "cm")),
                                                                 show_legend = T #don't show units,
                                                                 
  ), 
  show_annotation_name = T,
  annotation_name_offset = unit(1.2, "cm"),
  which = "column")
}

############################################
# generate mapping data column on 2 values #
############################################

generate_2_LibraryMappingAnnotation <- function( jointMappingData, library1Info, library2Info ){
  column_ha = HeatmapAnnotation(mapping = anno_barplot(as.matrix(jointMappingData), 
                                                       axis = TRUE,
                                                       gp = gpar( fill = c(library1Info$color, library2Info$color, "grey") , 
                                                                  lwd = unit(1, "cm")),
                                                       show_legend = T
  )
  ) 
}

"
get and return top n barcoddes for a number of libraries

Input:
1. curDF: 
the current data frame: barcodes = rows, sample file names = columns
2. libraries:
list of libraries formated:
- each element is a library :
- each library is a list of:
element 1: library Name
element 2: library ID (6 base-pair barcode ID)
element 3 (optional): barcodes of that library
Note: if have already separated barcodes by library then fill this element. otherwise function will call getLibraryBarcodeList
3. n: 
number of barcodes from each sample to grab
Default = 10
5. samples.by.cell.type: 
list result from included funciton: getSamplesByCellType
required if trying to separate by cell type
6. curCellType (optional): 
cell type with which to subset data sets. If not passed --> assumes to grab top n barcodes from all libraries

@Return:
return list of lists:
- each element of list = library
- inside each sample file list:
- numeric vector of barcodes
names = barcodes
values = proportion

"
topNBarcodesByLibrary <- function( curDF, libraries, n = 10 ){
  
  #get proportion table
  #TO DO: maybe: implemenet this to be saved in the curDF object
  curPropTable = getPropTable(curDF)
  topNBarcodes = list()
  #iterate through each library
  for(library in libraries){
    
    #check if library has barcodes already
    if(! checkLibrary(library)){
      stop("bad library format")
    }
    
    #get barcodes
    if( !("barcodes" %in% names(library) )   ){
      print( paste0("Barcodes not in library ", library, ". getting barcodes") )
      library[['barcodes']] = getLibraryBarcodeVector(curDF, library)
    }
    
    #get proportions of specific barcodes
    curLibProp <- subset(curPropTable, subset = rownames(curPropTable) %in% library$barcodes)
    #get top proportions for each barcode
    topPropBarcodes <- apply( curLibProp, MARGIN = 1, max)
    #order
    topPropBarcodes <- topPropBarcodes[order(topPropBarcodes, decreasing = T)]
    topNBarcodes[[library$name]] = topPropBarcodes[1:n]
  }
  
  if(length(topNBarcodes) == 0){
    stop("topNBarcodesByLibrary didn't catch any top barcodes. Check inputs")
  }
  topNBarcodes
}


"
topNBarcodesByLibrary_andCellType: get and return top n barcoddes for a number of libraries separated by cellType

@Input:
1. curDF: 
the current data frame: barcodes = rows, sample file names = columns
2. libraries: 
list of libraries formated:
element 1: library Name
element 2: library ID (6 base-pair barcode ID)
element 3 (optional): barcodes of that library
Note: if have already separated barcodes by library then fill this element. otherwise function will call getLibraryBarcodeList
3. n: 
number of barcodes from each sample to grab
5. samples.by.cell.type: 
list result from included funciton: getSamplesByCellType
6. curCellType (optional): 
cell type with which to subset data sets. If not passed --> assumes to grab top n barcodes from all libraries for every cell type
@Return:
return list of lists:
- each element of list = cell type name
- inside each cell type list element:
- list of libraries
- inside each library element: list of character vectors
1. barcode name
2. barcode ID
3. top n barcodes

@Output Summary
list of cell types -> list of sample names of cell type --> list of libraries in sample --> list of library info including top n barcodes
"
topNBarcodesByLibrary_andCellType <- function(curDF, libraries, n = 10, samples.by.cell.type, cell.type = null){
  #if pass cell type then subset samples to that cell.type
  topNBarcodes = list()
  cellTypes <- names(samples.by.cell.type)
  if(!(is.null(cell.type) ) ){
    cellTypes = cell.type
  }
  for(curCell in cellTypes ){
    cellTypeData <- curDF[ , samples.by.cell.type[[curCell]] ]
    topNBarcodes[[curCell]] = topNBarcodesByLibrary(curDF = cellTypeData, libraries = libraries, n = n)
  }
  if(length(topNBarcodes == 0)){
    stop("topNBarcodesByLibrary_andCellType didn't catch any barcodes. see input data")
  }
  topNBarcodes
}

"
get topNBarcodes for a sample: Get top n barcodes for a particular sample
INPUT:
1. curDF:
current data frame: rows = barcodes, cols = sample names, values = reads for each barcode
2. sampleName:
sample name to subset and look at
3. n:
number of clones

RETURNS:
character vector for top n barcodes
"
topNBarcodesForSample <- function(curDF, sampleName ,n = 10 ){
  if(class(curDF) != "data.frame"){
    stop("need to pass data frame for curDF")
  }
  if( !(sampleName %in% colnames(curDF)) ){
    stop( "can't find given sample name in data frame passed")
  }
  #get specific sample name vector
  sample.of.interest = curDF[,sampleName]
  names(sample.of.interest) = rownames(curDF)
  #get top n
  sample.of.interest = sample.of.interest[order(sample.of.interest, decreasing = T)]
  #return top n
  names(sample.of.interest)[1:n]
}

"
topNBarcodesForSample_byLibrary: Get top barcodes for across libraries for a particular sample
RETURNS:
List with elements:
1. sampleName
2. library A
3. library B
4. etc....
#To DO:
Maybe return library element with sample name and top barcodes
"
topNBarcodesForSample_byLibrary <- function(curDF, sampleName, libraries ,n = 10){
  topNBarcodesForSample_byLibrary = list(sampleName = sampleName)
  for(library in libraries){
    libBarcodes <- getLibraryBarcodeVector(curDF = curDF, library = library)
    if(is.null(libBarcodes)){
      stop("couldn't get baarcodes")
    }
    library[['barcodes']] <- libBarcodes
    #subset data frame
    libraryDF <- curDF[library$barcodes,]
    topBarcodes <-topNBarcodesForSample(curDF = libraryDF,sampleName = sampleName, n = n )
    topNBarcodesForSample_byLibrary[[library$name]] <- topBarcodes
  }
  topNBarcodesForSample_byLibrary
}

"
getLibraryBarcodeVector

@Description:
get all barcodes of a library in a data frame

@Inputs:
1. curDF (optional):
the current data frame: barcodes = rows, sample file names = columns
2. barcodes (optional need to pass this or curDF):
vector of all barcode names
3. library:
list of information about library:
1. library name
2. library ID (6 base pair barcode ID)

@Return:
list element of library where:
1. library name
2. library ID
3. character vector of all barcodes in that library found

@To Do:
implement way so that you can just pass barcode names instead of entire data frame --> can work with just samples
"
getLibraryBarcodeVector <- function( curDF=null, barcodes = NULL, library){
  # if("barcodes" %in% names(library)){
  #   print("already has barcodes")
  #   return(library$barcodes)
  # }
  if(!is.null(barcodes)){
    all.barcodes = barcodes
  }else{
    all.barcodes <- rownames(curDF)
  }
  
  barcodes <- all.barcodes[grep(all.barcodes, pattern = paste0("^",library["ID"]) ) ]
  
  if(length(barcodes) == 0 ){
    warning("didn't catch any barcodes")
  }
  barcodes
}

"
getSamplesByCellType: Get sample file names separated by cell type
INPUT:
1. curDF:
Curente data frame
2. delimiter:
delimiter separating sample name information:
example: file format: monkeyname_date_cellType_month_etc.
delimiter for above = '_'
3. cellTypeIndex : #Will be optional once 4. is ipemented
once name is split by delimiter, if know where cell type resides in name --> pass this index
example: file format: monkeyname_date_cellType_month_etc.
index = 3

Not Implemented:
4. example.cell.type (Optional):
If don't know which index cell type resides in, pass a cell name 
you know resides in one of data file. Make sure to make this a unique element residing in data frame (ex. _T_ or Grans )
rather than something likely to be in file name multiple times (ex. T )

Note: Need to pass either cellTypeIndex OR example.cell.type: if neither passed --> error thrown

RETURNS:
list of character vectors:
1. vector name: cell type
2. vector elements = sample names
"
getSamplesByCellType <- function( sample.names , delimiter = "_", cellTypeIndex){ #, example.cell.type = NULL){
  curSampleNames <- sample.names
  splitNames <- strsplit(curSampleNames, split = delimiter)
  
  #for each sample in colnames --> find out which cell type it is and add to list element for that cell type
  samples.by.cell.type <- list()
  for( sample in 1:length(splitNames) ){
    cell.type <- splitNames[[sample]][cellTypeIndex]
    
    #if cell.type vector not already an element in samples.by.cell.type then start a new one. 
    if( !(cell.type %in% names(samples.by.cell.type) ) ){
      samples.by.cell.type[[cell.type]] <- c(curSampleNames[sample])
    }else{ 
      #otherwise append name to end
      samples.by.cell.type[[cell.type]]= c(samples.by.cell.type[[cell.type]], curSampleNames[sample])
    }
  }
  samples.by.cell.type
}

"
getSmaplesByDate: return list of given samples, separated by their date
INPUT: vector of sample names
OUTPUT: list of dates and samples in each date

To Implement: regex way of figuring out where _8m_ part is in pattern
To Implement: warning if returning null object: warning doesn't work currently
"
getSamplesByDate <- function( sampleNameVector, delimiter = "_" , dateIndex){
  
  #split names
  splitNames <- strsplit(sampleNameVector, split = delimiter)
  
  samples.by.date <- list()
  for( sample in 1:length(splitNames) ){
    curDate <- splitNames[[sample]][dateIndex]
    #if date vector not already an element in samples.by.cell.type then start a new one. 
    if( !(curDate %in% names(samples.by.date) ) ){
      samples.by.date[[curDate]] <- c(sampleNameVector[sample])
    }else{ 
      #otherwise append name to end
      samples.by.date[[curDate]]= c(samples.by.date[[curDate]], sampleNameVector[sample])
    }
  }
  
  
  if(is.null(samples.by.date)){
    warning("returning empty samples.by.date")
  }
  samples.by.date
}

#############################
# change.names.to.cell.date #
#############################
"
@input: 
1. data.sample.names: usually column names of barcode data frame
2. delimiter: consistent barcode name (as with lauren's script) renders consistent naming scheme with delimiter between elements of sample: pass that character
#     Default: `_`
# example: name: M1105188_1m_T_XG4_sampled_i502_S2_L004_R1_001.fastq
#                delimiter: `_` 
3. cell.type.index: once names are divided by delimiter, where does the cell.type (ex. grans, mono,T) reside in that division
#       Default: 3
# example: name: M1105188_1m_T_XG4_sampled_i502_S2_L004_R1_001.fastq
#                cell.type.index = 3 (where T is)
4. date index: once names divided by delimiter, where does the date reside in that division
#       Default: 2
# example: name: M1105188_1m_T_XG4_sampled_i502_S2_L004_R1_001.fastq
#                date.index = 2 (where 1m is)
"
change.names.to.cell.date <- function(data.sample.names, delimiter = "_" , cell.type.index = 3, date.index = 2){
  if(sum(grepl(pattern = delimiter, x = data.sample.names)) < length(data.sample.names) ){
    stop("can't find delimiter in column names. Maybe already changed")
  }
  #let's change the names
  #separate by cell type
  samples.by.cell.type <- getSamplesByCellType(sample.names = data.sample.names, cellTypeIndex = 3)
  #separate by month
  samples.by.month <- lapply(X = samples.by.cell.type, FUN = function(cell.type.name.vector){ getSamplesByDate(cell.type.name.vector, delimiter = "_",dateIndex = 2) } )
  #now combine those values
  new.names = unlist(samples.by.month)
  #return new.names
  names(new.names)
}

"
check data frame:
check library:
test to see if format of library list object is correct:

correct input:
element 1: library name
character vector for the library
ex. '7'

element 2: library ID
library ID (ususally 6 base pair) sequence

To Do: implement data frame holding all library Information so it can be accessed at will
To Do: check elements of name and library correspond to repository of libray information
To Do: include monkey
"
checkLibrary <- function(library){
  if( class(library) == "list"){
    libnames = names(library)
    if(length(libnames) == 3){
      test = all(libnames == c("name", "ID", "barcodes"))
      if( !test ){
        print("bad library list format.")
        print( "current library list format: " )
        print( names(library) )
        print( "requisite library list format: ")
        print( "name, ID, barcodes")
      }
      return(test)
    }
    if(length(libnames) == 2){
      test=  all(libnames == c("name" , "ID") )
      if(!test){
        print("bad library list format.")
        print( "current library list format: " )
        print( names(library) )
        print( "requisite library list format: ")
        print( "name, ID")
      }
      return(test)
    }
  }else{
    print("library element needs to be a list")
    return(FALSE)
  }
}

checkDF <- function(curDF){
  if(is.null(curDF)){
    stop("null object passed to getPropDF()")
  }else if(0 %in% dim(curDF)){
    stop("passing empty data object")
  }
}

# getPropTable: return data frame of proportions for every barcode (as data frame)
getPropDF <- function(curDF){
  checkDF(curDF)
  #get proportion table
  #TO DO: maybe: implement this to be saved in the curDF object
  curPropTable = apply(curDF, MARGIN = 2 , function(col){ col/sum(col) })
  curPropTable[ is.nan(curPropTable) ] = 0
  data.frame(curPropTable)
}

custom_log <- function(x, log_choice, vary_log_min){ # Diego's Custom Log
  x <- log(x, log_choice)
  #sets the -inf values to be the minimum of the possible value in the data frame -1
  if(vary_log_min) {
    x[x == -Inf] <- (min(x[x > -Inf]) - 1)
  } else {
    x[x == -Inf] <- (log(100/4000000, log_choice)-1)
  }
  return(x)
}

# take distance vector and replace na's with replace value or max.amp.value * max(cur.dist.vect)
dist.replace.nas <- function(cur.dist.vect, replace.value = NULL, max.amp.value = 100){
  
  na.values = is.na(cur.dist.vect)
  if(is.null(replace.value)){
    cur.dist.no.na <- na.omit(cur.dist.vect)
    cur.max <- max(cur.dist.no.na)
    replace.value = cur.max * max.amp.value
  }
  cur.dist.vect[na.values] = replace.value
  cur.dist.vect
}
