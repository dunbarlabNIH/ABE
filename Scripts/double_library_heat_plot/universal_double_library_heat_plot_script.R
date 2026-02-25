####################################
## "HEAT PLOTS: Universal Script" ##
####################################

# user input
working_dir = "xyz" #"Put your directory's path here"
library.1.id = 'GAGTTC' #"Put in library 19 Barcode ID"
library.1.combine.file.path = file.path(working_dir, "MonkeyID_combined_date_date_lib19_corrected.txt")
library.1.readme.file.path = file.path(working_dir, "MonkeyID_combined_date_date_lib19.txt_README_corrected.txt")

library.2.id = 'GCAACT'  #"Put in library 21 Barcode ID"
library.2.combine.file.path = file.path(working_dir, "MonkeyID_combined_date_date_lib21_corrected.txt")
library.2.readme.file.path = file.path(working_dir, "MonkeyID_combined_date_date_lib21.txt_README_corrected.txt")

  
# Examples for column inputs: c(1,2,3,4)    or   c(3,5,7,8) 
# Can also dictate order you want columns to be in: c(5,3,2,6)
samples.column.numbers.of.interest.in.order.desired = c(# In vitro
                                                53, # 51,52,53,
                                                      #PB
                                                   #  43,36,15,16,17,18,19,14,21,25,1,4,5,6,7,9,10,11,13)
                                                      #BM
                                                   40, 31,22, 2, 8,12,
                                                      #LN
                                                    41,32,23,3,
                                                      #Liver
                                                   42,33,24,
                                                      #Spleen
                                                     44,26,
                                                      #B cell sites
                                                    49,50,47,48)
n.clones = 10

########################################################



###########################
#Step 0: Source necessary functions and packages:
###########################
double_library_dir = "xyz/double_library_heat_plot" #"Substitute xyz with your directory's path"
source(file.path(double_library_dir, "barcode_ggheatmap.R"))
source(file.path(double_library_dir, "double_library_heat_plot_functions.R"))
#load packages:
require(ComplexHeatmap)
require(circlize)


#Step 2: Define Library Annotations: 
library.1.annotations <- list("name" = 'library_1', "ID" = library.1.id, "color" = "gray")
library.2.annotations <- list("name" = 'library_2' , "ID" = library.2.id, "color" = "red")
#push together
library.anotations = list(library.1.annotations, library.2.annotations)


#Step 3: Read in data
lib.1.data = read.table( file = library.1.combine.file.path , header = T )
lib.1.readme.data = getReadMeData( library.1.readme.file.path )
lib.2.data = read.table( file = library.2.combine.file.path , header = T )
lib.2.readme.data = getReadMeData( library.2.readme.file.path )

#Step 4: Combine Data
comb.data <- rbind(lib.1.data , lib.2.data)

#Step 5: Get Mapping data on 2 libraries
both.mapping.data <- getMappingDataOn2Libraries( readmeData1 = lib.1.readme.data , 
                                                 libraryInfo1 = library.1.annotations,
                                                 readmeData2 = lib.2.readme.data ,
                                                 libraryInfo2 = library.2.annotations)

#Step 6: subset columns and get order by specified 
comb.data <- comb.data[,samples.column.numbers.of.interest.in.order.desired]
both.mapping.data <- both.mapping.data[ samples.column.numbers.of.interest.in.order.desired , ]

#Step 7: Run Diego's hclust function to get the top barcodes
diego.heat.plot.data = barcode_ggheatmap(your_data = comb.data , printtable = T, table_option = "logs", n_clones = n.clones)
#now lets get the propDF
curProp <- getPropDF(comb.data[rownames(diego.heat.plot.data),])

#Step 8: organize data so that libraries are separate
lib.1.barcodes<- grepl(x = rownames(diego.heat.plot.data), pattern = paste0("^",library.1.annotations$ID) )
lib.2.barcodes <- grepl(x = rownames(diego.heat.plot.data), pattern = paste0("^",library.2.annotations$ID) )
lib.1.barcode.order <- hclust( dist(curProp[lib.1.barcodes,] , method = "euclidean") )$order
lib.1.barcode.order <- rownames(subset(curProp, subset = lib.1.barcodes))[lib.1.barcode.order]
lib.2.barcode.order <- hclust( dist(curProp[lib.2.barcodes,] , method = "euclidean") )$order
lib.2.barcode.order <- rownames(subset(curProp, subset = lib.2.barcodes))[lib.2.barcode.order]

both.barcode.order <- c(lib.1.barcode.order, lib.2.barcode.order)
diego.heat.plot.data <- diego.heat.plot.data[both.barcode.order,]

#Step 9: make library barcode data object for side bar:
barcode.colors = sapply(X = rownames(diego.heat.plot.data),
                        FUN = function(barcode){
                          if( grepl(pattern = paste0("^",library.1.annotations$ID), x = barcode) ){
                            "lib1"
                          }else{
                            "lib2"
                          }
                        }
)
#turn into a dataframe
barcode.colors <- data.frame(barcode.colors , stringsAsFactors = F)
colnames(barcode.colors) = "library"

#Step 10 Generate Side Bar Heatmap Annotation
barcode.source <- HeatmapAnnotation( df = barcode.colors, 
                                     col = list(library = c(lib1 = library.1.annotations$color, lib2 = library.2.annotations$color)), 
                                     which = 'row')

#Step 11: Generate mapping data annotation
both.mapping.annotation <- generate_2_LibraryMappingAnnotation(both.mapping.data, 
                                                               library1Info = library.1.annotations, 
                                                               library2Info = library.2.annotations)


#Step 12: Plot
log_choice <- exp(1)
vals = c(log(100/4000000, log_choice) - 1, log(100/4000000, log_choice), log(0.001, log_choice), log(0.01,log_choice), log(0.1, log_choice), 0)
labels = c("Defined 0", "0.0025% ", "0.1%", "1%", "10%", "100%")
Heatmap(as.matrix(diego.heat.plot.data), 
        show_row_dend = FALSE,     # don't show dendrograms
        show_column_dend = FALSE, 
        show_row_names = FALSE,    # don't show names
        row_order = rownames(diego.heat.plot.data), # specify order so goes faster
        cluster_columns = FALSE,   #don't cluster the columns
        cluster_rows = FALSE,
        col = colorRamp2(vals, c("#4575B4", "#4575B4", "lightblue", "#fefeb9", "#D73027", "red4")),  #specify colors
        #col = colorRamp2( vals , c("#4575B4", "#4575B4", "lightblue", "#fefeb9", "#D73027", "red4") ),
        heatmap_legend_param = list(at = vals,
                                    labels = labels, #vals,#, 
                                    color_bar = 'continuous',
                                    title = "%",
                                    labels_gp = gpar(fontsize = 10),
                                    legend_height = unit(5, "cm")
        ),
        rect_gp = gpar(col = "black", lwd = 0.30),  # add borders between heatmap elements,
        top_annotation = both.mapping.annotation
)+ barcode.source



########################################
##Extra:   Generate Plot With Normalized Stacked Bar Plots
########################################
no.junk.mapped = both.mapping.data[,c(1,2)]
for( i in 1: length(no.junk.mapped[,1])){
  cur.row = no.junk.mapped[i,]
  total.mapped.reads = cur.row[1] + cur.row[2]
  no.junk.mapped[i,1] = cur.row[1]/total.mapped.reads
  no.junk.mapped[i,2] = cur.row[2]/total.mapped.reads
}

no.junk.mapping.annotation <- generate_2_LibraryMappingAnnotation(no.junk.mapped, 
                                                                  library1Info = library.1.annotations, 
                                                                  library2Info = library.2.annotations)

#Plot:
vals = c(log(100/4000000, log_choice) - 1, log(100/4000000, log_choice), log(0.001, log_choice), log(0.01,log_choice), log(0.1, log_choice), 0)
labels = c("Defined 0", "0.0025% ", "0.1%", "1%", "10%", "100%")
Heatmap(as.matrix(diego.heat.plot.data), 
        show_row_dend = FALSE,     # don't show dendrograms
        show_column_dend = FALSE, 
        show_row_names = FALSE,    # don't show names
        row_order = rownames(diego.heat.plot.data), # specify order so goes faster
        cluster_columns = FALSE,   #don't cluster the columns
        cluster_rows = FALSE,
        col = colorRamp2(vals, c("#4575B4", "#4575B4", "lightblue", "#fefeb9", "#D73027", "red4")),  #specify colors
    #  show_column_names =FALSE, # don't show names
        #col = colorRamp2( vals , c("#4575B4", "#4575B4", "lightblue", "#fefeb9", "#D73027", "red4") ),
        heatmap_legend_param = list(at = vals,
                                    labels = labels, #vals,#, 
                                    color_bar = 'continuous',
                                    title = "%",
                                    labels_gp = gpar(fontsize = 10),
                                    legend_height = unit(5, "cm")
        ),
        rect_gp = gpar(col = "black", lwd = 0.25),  # add borders between heatmap elements,
        top_annotation = no.junk.mapping.annotation
)+ barcode.source








