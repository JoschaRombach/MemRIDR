# Functions for HD-peptide Array analysis
# Joscha Rombach
# 24/03/2023

# Note - SWaFi stands for Smoothed Weighted-average Fluorescence intensity,
# which is referred to as "binding score" in the publication

#-############################################################################-#
#-########################### INSTALL PACKAGES ###############################-#
#-############################################################################-#

# Install packages from CRAN, Bioconductor and Github
# cran <- c('LearnGeom','tidyverse','foreach','doParallel','berryFunctions',
#           'parallelly','seqinr','gtools','broom','readxl','tools','stringr',
#           'data.table','viridis','ggrepel','stringr','dynamicTreeCut',
#           'plotrix','fields','circlize','gamlss','gamlss.dist','gamlss.add',
#           'ggExtra','ggprism','bio3d','ggnewscale','reshape2','helixvis',
#           'gridExtra','pbapply','zoo','magick','rbioapi','RColorBrewer',
#           'misc3d','eulerr','MASS','svMisc','umap','stringdist','devtools')

cran <- c('tidyverse','berryFunctions','seqinr','stringr','data.table',
            'viridis','bio3d','ggprism','ggExtra','gridExtra','zoo','RColorBrewer',
            'shiny','DT','shiny.fluent','imola','cowplot','shiny.router','shiny.react',
            'shinyWidgets','circlize','msa','ggtree')

bioconductor <- c('msa','ggtree')#,'ComplexHeatmap','fgsea')

github_packages <- c('r3dmol')#,'drawCell')

# Function for installing missing packages from CRAN and Bioconductor
install_CRAN_and_Bioconductor_packages <- function(cran,bioconductor){
  
  cran_missing <- cran[which(!cran %in% rownames(installed.packages()))]
  bioconductor_missing <- bioconductor[which(!bioconductor %in% rownames(installed.packages()))]
  n_missing <- length(cran_missing)+length(bioconductor_missing)
  
  message(paste0('Installing ',n_missing,' packages from CRAN and Bioconductor'))
  
    if( n_missing > 0){
    
    continue <- NA
    while(!continue%in%c('y','Y','ye','YE','Ye','yes','YES','Yes','YeS','YEs','yeS','yES','n','N','no','No','NO')){
      continue <- readline(paste0("There are ",n_missing," packages that need to be installed - Do you want to install them? "))
    }
    
    if(continue%in%c('y','Y','ye','YE','Ye','yes','YES','Yes','YeS','YEs','yeS','yES')){
      
      if(length(cran_missing)>0){
        install.packages(cran_missing)
      }
      
      if(length(bioconductor_missing)>0){
        if (!require("BiocManager", quietly = TRUE))
          install.packages("BiocManager")
        
        BiocManager::install(bioconductor_missing)
      }
    }
    
    cat('\n')

  } else {
    message('All packages already installed')
  }
}

# Install from CRAN and Bioconductor
install_CRAN_and_Bioconductor_packages(cran,bioconductor)

# Install r3dmol from Github
message('Testing whether to install r3dmol from Github')
if(!'r3dmol' %in% rownames(installed.packages())){
  
  continue <- NA
  while(!continue %in% c('yes','YES','Yes','YeS','YEs','yeS','no','No','NO')){
    continue <- readline("Install r3dmol from page swsoyee/r3dmol? (yes/no) : ")
  }
  
  if(continue%in%c('yes','YES','Yes','YeS','YEs','yeS')){
       devtools::install_github("swsoyee/r3dmol", upgrade = 'never')
  }
  
  rm(continue)
  
} else {
  message('r3dmol already installed')
}
    
# if(!'drawCell' %in% rownames(installed.packages())){
#    devtools::install_github("svalvaro/drawCell")
# }

#-############################################################################-#
#-########################### LOAD PACKAGES ##################################-#
#-############################################################################-#

# Function for loading packages
load_packages <- function(packages){
  message("Loading packages ...")
  for (p in packages) {
    if (p %in% rownames(installed.packages())) {
      suppressPackageStartupMessages(library(p, character.only=TRUE))
    } else {
      message(paste0('Error: package ',p,' is not installed ... '))
    }
  }
}

# Load packages
load_packages(c(cran, bioconductor, github_packages))

rm(cran, 
   bioconductor, 
   github_packages,
   install_CRAN_and_Bioconductor_packages, 
   load_packages)

#-############################################################################-#
#-###################### FUNCTION TO IMPORT FUNCTIONS ########################-#
#-############################################################################-#

Import_functions <- function(path='code/Functions.R'){
  if(!require('here')){install.packages('here')}
  library('here')
  source(here(path))
}

#-############################################################################-#
#-############################# AESTHETICS ###################################-#
#-############################################################################-#

# Set ggplot theme for script
message("Setting ggplot theme to prism style ...")
theme_set(new = theme_prism())

# Custom prism theme for ggplot
theme_update(axis.line  = element_line(linewidth = .5),
             axis.ticks = element_line(linewidth = .5),
             legend.position = 'top', 
             aspect.ratio = 1, 
             axis.text  = element_text(face = 'plain'), 
             axis.title = element_text(face = 'plain'))

# Optional axis offset prism style
prism <- function(){
         guides(x = "prism_offset", y = "prism_offset")
}

# Easy functions for removing axes
remove_x_axis <- function(){
                 list(
                 theme(axis.title.x = element_blank(),
                       axis.ticks.x = element_blank(),
                       axis.line.x  = element_blank(),
                       axis.text.x  = element_blank())
                 )}
remove_y_axis <- function(){
  list(
    theme(axis.title.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.line.y  = element_blank(),
          axis.text.y  = element_blank())
  )}

# Add large color bar to the right matplotlib style
matplotlib_colorbar <- function(position = 'right',
                                specify_plot_dim_for_saving = NULL, 
                                adjust   =  8){
  
  if(position %in% c('top','bottom')){
    angle <- 0
  }
  if(position == 'right'){
    angle <- 270
  }
  if(position == 'left'){
    angle <- 90
  }
  
  if(is.null(specify_plot_dim_for_saving)){
    specify_plot_dim_for_saving <- dev.size()
  }
  
  if(position %in% c('left','right')){
    res <-  theme(legend.position = position, 
                  legend.title.position = position, 
                  legend.title = element_text(angle = angle),
                  legend.key.height = unit(specify_plot_dim_for_saving[2]/adjust, "inches"))
  }
  
  if(position %in% c('top','bottom')){
    res <-  theme(legend.position = position, 
                  legend.title.position = position, 
                  legend.title = element_text(angle = angle),
                  legend.key.width = unit(specify_plot_dim_for_saving[1]/adjust, "inches"))
  }
  
  return(res)
  }

# Set progress bar settings
pbapply::pboptions(char='=')

#-############################################################################-#
#-########################### SHINY APP FUNCTIONS ############################-#
#-############################################################################-#
search_for_protein <- function(search, keys){
  
  if(search==''){
    search <- 'vGluT1'
  }
  
  res <- keys[str_detect(keys$Search,tolower(search)),]
  
  if(length(res)==0){
    res <- keys[keys$UniProtKB.AC=='Q9P2U7',]
  }
  
  return(res)
}


#-############################################################################-#
#-################### FUNCTIONS FOR IMPORTING NEW DATASETS ###################-#
#-############################################################################-#

message("Loading functions ...")

# Function for reading HD array files
# read.array <- function(file, 
#                        cols=9, 
#                        input.order=T, 
#                        doCalcs=F, 
#                        legacy=F){
#   
#   if(input.order==F & cols==9){
#     cols <- 8
#   }
#   
#   #Reading in .csv data file         
#   RawData <- read.csv(file, sep="\t")
#   
#   #Number of cols in data file
#   no_cols <- cols
#   
#   #Start of Dataframe
#   Start <- as.numeric(which(grepl(">------------------------------",
#                                   RawData[,1])))[1]
#   #Removing Data Commentary
#   RawData <- as.data.frame(RawData[c((Start+3):length(RawData[,1])),]) 
#   
#   FillRows <- as.numeric(which(grepl(">------------------------------",
#                                      RawData[,1])))
#   
#   #Removing subheaders
#   y <- 0
#   for(i in FillRows){
#     print(paste("Processing Line:", i))
#     i <- i+y
#     RawData <- as.data.frame(RawData[-c(i:(i+no_cols+2)),])
#     y <- y-(no_cols+3)
#   }
#   
#   if(legacy==T){RawData <- RawData[RawData[,1] != "NO GROUP",]}
#   
#   if(legacy==F){
#     #Finding rows with missing data
#     FillRows <- as.numeric(which(grepl("-0001",RawData[,1])))
#     
#     Raw_mtrx <- as.matrix(RawData)
#     Raw_Data <- matrix(NA,nrow=(nrow(RawData)+(length(FillRows)*3)),ncol=1)
#     Empty <- sort(c(FillRows+5,FillRows+6,FillRows+7))+(sort(c(rep(0:(length(FillRows)-1),3)))*3)
#     Raw_Data[-Empty,] <- Raw_mtrx 
#     remove(Empty)
#   }
#   
#   remove(y,i,FillRows,Start)
#   
#   #Creating Data Matrix
#   Organised_Data <- as.data.frame(unique(t(matrix(Raw_Data,no_cols,byrow=F))))
#   
#   column_names <- c("signal","coresequence","sector","row",
#                     "col","signal-µblanks","z-score","p-value")
#   
#   if(input.order==T){column_names <- c("input",column_names)}
#   
#   if(legacy==T){column_names <- c("signal","coresequence","sector","row","col")}
#   
#   colnames(Organised_Data) <- column_names
#   
#   Organised_Data <- Organised_Data[-1,]
#   for(i in which(colnames(Organised_Data) %notin% "coresequence")){
#     Organised_Data[,i]<-as.numeric(Organised_Data[,i])}
#   
#   Organised_Data <- Organised_Data[order(Organised_Data[,"sector"], 
#                                          Organised_Data[,"row"], 
#                                          Organised_Data[,"col"]),]
#   
#   # Removing rawdata
#   remove(RawData,Raw_mtrx,Raw_Data,no_cols)
#   
#   # Creating specific peptide location ID
#   Organised_Data$location <- paste(Organised_Data$sector,
#                                    Organised_Data$row,Organised_Data$col, 
#                                    sep = ".")
#   
#   Organised_Data$signal[Organised_Data$signal==-1] <- NA
#   
#   Blanks <- Organised_Data[nchar(Organised_Data$coresequence)<=5,]
#   Blanks <- Blanks[order(Blanks$coresequence),]
#   rownames(Blanks) <- NULL
#   
#   Organised_Data <- Organised_Data[nchar(Organised_Data$coresequence)>5,]
#   
#   rownames(Organised_Data) <- NULL
#   
#   if(doCalcs == T){
#     # Calculating peptide properties
#     print("Calculating peptide properties - this step takes ~5min")
#     Organised_Data <- PeptideProperties(Organised_Data, "coresequence", pH=7.4)
#     
#     # Unlisting columns
#     Organised_Data[] <- lapply(Organised_Data, function(x){
#       if(is.list(x)) unlist(x) else x})
#     
#     # Making Cluster
#     cl <- makeCluster(detectCores())
#     registerDoParallel(cl)
#     
#     # Statement 
#     print("Calculating avg. charges & average signals surrounding each field")
#     print("This step takes a while depending on the amount of threads utilised")
#     print("System time: ~20min")
#     print("Intel(R) Core(TM) i7-8665U CPU @ 1.90GHz - 2.11 GHz + 16GB RAM")
#     
#     # Calculating the surrounding charges & z-scores
#     Surrounding <- foreach (i = 1:nrow(Organised_Data), .combine = rbind, 
#                             .multicombine = T, .inorder = T) %dopar% {
#                               
#         j <- as.numeric(Organised_Data[i, "row"])
#         c  <- as.numeric(Organised_Data[i, "col"])
#         sec <- as.numeric(Organised_Data[i, "sector"])
#         
#         n <- j-1
#         s <- j+1
#         e <- c+1
#         w <- c-1
#         
#         points <- data.frame(AdjacentRows = c(n,n,n,j,j,s,s,s), 
#                              AdjacentCols = c(e,c,w,e,w,e,c,w), 
#                              AdjacentSector = c(rep(sec, 8)))
#         colnames(points)[3] <- "AdjacentSector"
#         points <-  points[points$AdjacentRows > 0 & points$AdjacentCols > 0,]
#         points <- transform(points,location = paste(AdjacentSector, 
#                                                     AdjacentRows, AdjacentCols, 
#                                                     sep = "."))
#         
#         filtered <-Organised_Data[Organised_Data$location %in% points$location,]
#         
#         Adjacent <- c(mean(filtered$signal), mean(filtered$z.cont))
#         names(Adjacent) <- c("Surrounding.signal", "Surrounding.z")
#         
#         Adjacent
#       }
#     
#     # Stopping Cluster to free up resources 
#     stopCluster(cl)
#     
#     # Binding dataframes
#     Organised_Data <- cbind(Organised_Data, as.data.frame(Surrounding))
#     
#     # removing unneccesary matrix
#     remove(Surrounding)
#     
#     print("Locating charged residues - takes ~ 3min")
#     
#     # Creating and registering cluster for parallel processing
#     cl <- makeCluster(detectCores())
#     registerDoParallel(cl)
#     
#     # Finding positions of charged residues
#     Charge.position <- foreach(i = 1:nrow(Organised_Data), .combine = rbind,
#                                .multicombine = T, .inorder = T, 
#                                .packages = c("stringr", "tidyverse")) %dopar% {
#                                  
#                # Creating characters containing positions
#                if(Organised_Data[i, "No.of.Charges"] >= 1){
#                  temp <- unlist(str_split(Organised_Data[i,"coresequence"], ""))
#                  
#                  pos <- as.character(which(temp %in% c("R","K"))) %>%
#                    str_flatten(collapse = ",")
#                  
#                  neg <- as.character(which(temp %in% c("E","D"))) %>%
#                    str_flatten(collapse = ",")
#                  
#                } else{pos <- NA
#                neg <- NA}
#                
#                # creating tibble with positions of positive and negative charges
#                strloc <- tibble("Locate.K.R" = pos, "Locate.D.E" = neg)
#                strloc[strloc == ""] <- NA
#                strloc
#              }
#     
#     # Stopping cluster to free up resources
#     stopCluster(cl)  
#     remove(cl)
#     
#     # Extracting positions from a single string to a list of characters
#     Charge.position$Locate.K.R <- str_extract_all(Charge.position$Locate.K.R,
#                                                   pattern = "[[:digit:]]+")
#     Charge.position$Locate.D.E <- str_extract_all(Charge.position$Locate.D.E,
#                                                   pattern = "[[:digit:]]+")
#     
#     # Binding dataframes
#     Organised_Data <- cbind(Organised_Data, Charge.position)
#     
#     # removing unnecessary tibble
#     remove(Charge.position)
#   }
#   return(list(Peptides = Organised_Data, Blanks = Blanks))
# }
# 
# # Functions for gradient correction and normalization
# convert_to_MatrixList <- function(data, c=T){
#   
#   if(c==T){
#     data <- do.call(rbind,data)
#   }
#   
#   secs <- 1:max(data$sector)
#   m <- lapply(secs, function(x){
#     sector <- data[data$sector==x,]
#     sector <- sector[order(sector$row,sector$col),]
#     matrix(sector$signal, nrow=max(sector$row), byrow=T)
#   })
#   names(m) <- paste0('sector_',secs)
#   return(m)
# }
# 
# smooth_array <- function(MatrixList, 
#                          aRange=c(20,20,20), 
#                          get_gradient=0){
#   
#   #c <- F
#   #alldata <- do.call(rbind,rawdata)
#   #sectors <- convert_to_MatrixList(data=alldata, c=c)
#   
#   sectors <- lapply(MatrixList, function(x){
#     # first 2d spline 
#     s1 <- fields::image.smooth(x, aRange = aRange[1])$z
#     s1 <- (max(s1)/s1)
#     x <- x*s1
#     
#     # second 2d spline
#     s2 <- fields::image.smooth(x, aRange = aRange[2])$z
#     s2 <- max(s2)/s2
#     y <- x*s2
#     
#     s3 <- fields::image.smooth(y, aRange = aRange[3])$z
#     s3 <- max(s3)/s3
#     z <- y*s3
#     
#     s4 <- fields::image.smooth(z, aRange = aRange[3])$z
#     s4 <- max(s4)/s4
#     
#     g <- list(s1,s2,s3,s4)
#     
#     if(get_gradient==0){
#       return(z)
#     } else {
#       return(g[[get_gradient]])
#     }
#   })
#   return(sectors)
# }
# 
# melt_MatrixList <- function(MatrixList,
#                             rawdata){
# 
#   alldata <- do.call(rbind,rawdata)
#   
#   molten <- lapply(1:length(MatrixList),function(x){
#     y <- reshape2::melt(MatrixList[[x]],value.name = 'normalised')
#     colnames(y) <- c('row','col','normalised')
#     y <- merge(x=alldata[alldata$sector==x,],y=y,by=c('row','col'))
#     return(y)
#   })
#   return(molten)
# }
# 
# normalise_array <- function(MatrixList,
#                             rawdata, 
#                             control='LKKALKKLKKALKKAL', 
#                             blanks='', 
#                             subtract_blank=T){
#   
#   n <- melt_MatrixList(MatrixList,rawdata)
#   
#   n <- lapply(n,function(y){
#     b <- mean(y$normalised[y$coresequence==blanks],na.rm=T)
#     h <- mean(y$normalised[y$coresequence==control],na.rm=T) 
#     return(c(b,h))
#   })
#   
#     m <- lapply(1:length(MatrixList),function(x){
#       if(subtract_blank==T){
#         y <- MatrixList[[x]]-n[[x]][1] #subtract empty field
#       } else {
#         y <- MatrixList[[x]]
#       }
#     y <- y/n[[x]][2] #normalise to hecate
#       if(subtract_blank==T){
#         y[y<0] <- 0
#       }
#     return(y)
#   })
#   
#   return(m)
# }
# 
# # Function for plotting array
# plot.array <- function(MatrixList,
#                        dims=c(2,3),
#                        byRow=T,
#                        col=NULL,
#                        cRange=NULL){
#   
#   MatrixList <- lapply(MatrixList, function(x){
#     x <- rbind(rep(NA,ncol(x)),rep(NA,ncol(x)),x,rep(NA,ncol(x)),rep(NA,ncol(x)))
#     x <- cbind(rep(NA,nrow(x)),rep(NA,nrow(x)),x,rep(NA,nrow(x)),rep(NA,nrow(x)))
#     return(x)
#   })
#   
#   l <- length(MatrixList)/dims[1]
#   m <- do.call(rbind,lapply(c(1,l+1),function(x){
#     unlist(do.call(cbind,MatrixList[x:(x+l-1)]))
#   }))
#   
#   if(is.null(col)){
#     t <- ceiling(max(m,na.rm=T))
#     b <- floor(min(m,na.rm=T))
#     if(!is.null(cRange)){
#       t <- cRange[2]
#       b <- cRange[1]
#     }
#     col = colorRamp2(seq(b,t,t/10),magma(n=length(seq(b,t,t/10))))
#   }
#   
#   m <- rbind(rep(NA,ncol(m)),rep(NA,ncol(m)),m,rep(NA,ncol(m)),rep(NA,ncol(m)))
#   m <- cbind(rep(NA,nrow(m)),rep(NA,nrow(m)),m,rep(NA,nrow(m)),rep(NA,nrow(m)))
#   
#   Heatmap(m, cluster_columns = F, cluster_rows = F, col = col,
#           heatmap_legend_param = list(title='Intensity (a.u.)'), na_col='black')
# }

#-############################################################################-#
#-###################### FUNCTIONS FOR DATASET COMPARISONS ###################-#
#-############################################################################-#

# Function for comparing datasets and extracting overlapping sites
# Should be able to compare any number of data sets at once, 
# but I haven't tried it yet ... 
# Might need to be optimized for comparisons of more than 2 datasets
# Option to set overlap=1 when sites are within each other added on 30/11/2023
# compare_datasets <- function(datasets, 
#                              pthresh=0.05,
#                              overlap_threshold=0,
#                              if_within_set_overlap_to_1=T,
#                              binding=T,
#                              useAdjusted = F,
#                              get_min=F,
#                              proteins=NULL,
#                              array_object=NULL){
#   
#   # Find array_object in environment
#   if(is.null(array_object)){
#     array_object <- Find_array_object()
#   }
#   
#   if(length(pthresh)==1){pthresh <- rep(pthresh,length(datasets))}
#   
#   cat(paste('Extracting SWaFi binding sites ... \n'))
#   # Get SWaFi scores for all datasets and combine into one data.frame
#   s <- do.call(rbind,lapply(1:length(datasets),function(dataset){
#     SWaFi <- get_SWaFi_binding_site_details(dataset=datasets[dataset],
#                                             pthresh=pthresh[dataset],
#                                             useAdjusted = useAdjusted,
#                                             get_min = get_min,
#                                             proteins = proteins,
#                                             array_object = array_object, 
#                                             binding = binding)
#     SWaFi$from <- datasets[dataset]
#     return(SWaFi)
#   }))
#   
#   s <- as.data.frame(s)
#   
#   # Split data.frame into list by Array identifiers ("UniProtKB.AC_Terminal")
#   s <- split(s,s$UniProtKB.AC)
#   
#   cat('Comparing datasets ... \n')
#   # Identify overlaps of bindingsites between datasets for each protein
#   res <- pblapply(s,function(x){ # x is a data.frame
#     
#     # Get protein ID
#     name <- x$UniProtKB.AC[1]
# 
#     # Create identifier for each region of each protein for each dataset
#     id_list <- lapply(datasets,function(dataset){
#       
#       if(length(x$end[x$from==dataset])>0){
#         ids <- 1:length(x$end[x$from==dataset])
#       } else {
#         ids <- NULL
#       }
#       return(ids)
#     })
#     
#     # Add to x
#     x$new_id <- unlist(id_list)
#     
#     # Identify all possible dataset comparisons 
#     combs <- combinations(n = length(datasets), r = 2)
#     combs <- split.data.frame(combs,1:nrow(combs))
#     
#     # Create an IRanges object for all sites in each dataset (stored as list)
#     site_ranges <- lapply(seq_along(datasets),function(n){
#       IRanges(start=x$start[x$from==datasets[n]],end = x$end[x$from==datasets[n]], names = id_list[[n]])
#     })
#     
#     # Find overlaps between all sites in all datasets
#     ol <- lapply(combs,function(y){
#       as.data.frame(IRanges::findOverlaps(site_ranges[[y[1,1]]],site_ranges[[y[1,2]]]))
#     })
#     
#     # Calculate the degree of overlap of all sites in all datasets
#     degree_ol <- lapply(seq_along(combs),function(n){
#       y <- combs[[n]] # get the datasets to compare
#       
#       if(nrow(ol[[n]])>0){ # if statement to test that overlaps exist
#         
#         # calculate intersect between sites
#         intsec <- as.data.frame(IRanges::pintersect(site_ranges[[y[1,1]]][ol[[n]]$queryHits],site_ranges[[y[1,2]]][ol[[n]]$subjectHits]))
#         
#         # Identify the total range covered by both of the sites being compared 
#         combined <- c()
#         for(i in 1:nrow(ol[[n]])){
#           combined <- c(combined,IRanges::reduce(c(site_ranges[[y[1,1]]][ol[[n]]$queryHits[i]],site_ranges[[y[1,2]]][ol[[n]]$subjectHits[i]])))
#         }
#         
#         full_range <- as.data.frame(do.call(c,combined))
#         
#         contained_within <- c()
#         for(i in 1:nrow(ol[[n]])){
#           frw <- IRanges::findOverlaps(site_ranges[[y[1,1]]][ol[[n]]$queryHits[i]],site_ranges[[y[1,2]]][ol[[n]]$subjectHits[i]],type='within')
#           rev <- IRanges::findOverlaps(site_ranges[[y[1,2]]][ol[[n]]$subjectHits[i]],site_ranges[[y[1,1]]][ol[[n]]$queryHits[i]],type='within')
#           contained_within <- c(contained_within, ifelse(nrow(as.data.frame(c(frw,rev)))>0,T,F))
#         }
#         
#         # Create data.frame with results
#         res <- ol[[n]]
#         res$intersect <- intsec$width
#         res$full_width <- full_range$width
#         res$contained_within <- contained_within
#         res$percent_overlap <- res$intersect/res$full_width # percent overlap
#         
#         if(if_within_set_overlap_to_1==T){
#           res$percent_overlap[res$contained_within==T] <- 1
#         }
#         
#         
#       } else {
#         # dummy data.frame if no overlaps exist
#         res <- data.frame(queryHits=NA,subjectHits=NA,intersect=NA,full_width=NA,percent_overlap=NA)[-1,]
#       }
#       return(res)
#     })
#     
#     # Account for overlap of multiple binding sites ---------------------------
#     
#     for(i in seq_along(combs)){
#       
#       # extract datasets to compare
#       col1 <- combs[[i]][1,1]
#       col2 <- combs[[i]][1,2]
#       dataset1 <- datasets[col1]
#       dataset2 <- datasets[col2]
#       
#       # get overlap data
#       o <- degree_ol[[i]]
#       
#       # test if there are sites that overlap with multiple other sites first as dataset1 vs dataset2, 
#       # if so we need to duplicate rows in x to accompany this ...
#       if(sum(duplicated(o$queryHits))>0){
#         
#         dup <- o$queryHits[duplicated(o$queryHits)] # test to see if anything is duplicated, indicate multiple overlaps
#         
#         already_dup <- x$new_id[x$from==dataset1][duplicated(x$new_id[x$from==dataset1])] # test if this has already been corrected for another dataset
#         
#         for(ad in already_dup){
#           dup <- dup[-match(ad, dup)] # remove previously corrected duplicates (important so that we don't add unnecessary rows)
#         }
#         
#         # test if any duplicates are identified
#         if(sum(is.na(dup))==0 & length(dup)>0){
#           # duplicate rows
#           for(d in dup){
#             dup_row <- rownames(x[x$new_id%in%d & x$from==dataset1,])[1]
#             dup_row <- which(rownames(x)==dup_row)
#             x <- x[append(1:nrow(x), dup_row, after = dup_row),]
#           }
#         }
#       } 
#       
#       # and vice versa, dataset2 vs dataset1 (almost same code as above)
#       if(sum(duplicated(o$subjectHits))>0){
#         
#         dup <- o$subjectHits[duplicated(o$subjectHits)]
#         
#         already_dup <- x$new_id[x$from==dataset2][duplicated(x$new_id[x$from==dataset2])]
#         
#         for(ad in already_dup){
#           dup <- dup[-match(ad, dup)]
#         }
#         
#         if(sum(is.na(dup))==0 & length(dup)>0){
#           for(d in dup){
#             dup_row <- rownames(x[x$new_id%in%d & x$from==dataset2,])[1]
#             dup_row <- which(rownames(x)==dup_row)
#             x <- x[append(1:nrow(x), dup_row, after = dup_row),]
#           }
#         }
#       }
#     }
#     
#     # Create empty matrices ---------------------------------------------------
#     
#     # it is faster to fill marices than build them along the way so we pre-compute matrices to fill out later
#     m1 <- matrix(NA,nrow=nrow(x),ncol=length(datasets)) # matrix to store percentage of overlap between regions
#     m2 <- matrix(NA,nrow=nrow(x),ncol=length(datasets)) # matrix to store the IDs of the regions that overlap
#     
#     # add 1 to all diagonal values (because each dataset overlaps 100% with itself)
#     f <- 0
#     for(j in 1:length(datasets)){
#       if(nrow(x[x$from==datasets[j],])>0){
#         m1[(1:nrow(x[x$from==datasets[j],])+f),j] <- 1
#         f <- f+nrow(x[x$from==datasets[j],])
#       }
#     }
#     
#     
#     # Fill out matrices -------------------------------------------------------
#     
#     for(i in 1:length(combs)){
#       
#       # get names of datsets to compare
#       col1 <- combs[[i]][1,1]
#       col2 <- combs[[i]][1,2]
#       dataset1 <- datasets[col1]
#       dataset2 <- datasets[col2]
#       
#       # get overlap data
#       o <- degree_ol[[i]]
#       
#       # Extract fill for fill matrix m1 for dataset1
#       p_overlap1 <- x$new_id[x$from%in%dataset1] # get bindingsite ids of dataset1
#       overlap_with1 <- p_overlap1 # duplicate
#       
#       p_overlap1[p_overlap1%notin%o$queryHits] <- 0 # add 0 to all sites with no overlap
#       p_overlap1[p_overlap1%in%o$queryHits]   <- o$percent_overlap # add percentage overlap
#       
#       # Extract fill for fill matrix m2 for dataset1
#       overlap_with1[overlap_with1%notin%o$queryHits] <- NA # add NA to all sites with no overlap ID 
#       overlap_with1[overlap_with1%in%o$queryHits]   <- x$ID[x$from==dataset2 & x$new_id%in%o$subjectHits] # Add ID of the regions that overlap
#       
#       # Extract fill for fill matrix m1 for dataset2
#       p_overlap2 <- x$new_id[x$from%in%dataset2] # get bindingsite ids of dataset2
#       overlap_with2 <- p_overlap2 # duplicate
#       
#       p_overlap2[p_overlap2%notin%o$subjectHits] <- 0 # add 0 to all sites with no overlap
#       p_overlap2[p_overlap2%in%o$subjectHits] <- o$percent_overlap # add percentage overlap
#       
#       # Extract fill for matrix m2 for dataset2
#       overlap_with2[overlap_with2%notin%o$subjectHits] <- NA # add NA to all sites with no overlap ID
#       overlap_with2[overlap_with2%in%o$subjectHits] <- x$ID[x$from==dataset1 & x$new_id%in%o$queryHits] # Add ID of the regions that overlap
#       
#       # fill m1 matrices
#       m1[which(x$from==dataset1),col2] <- p_overlap1 
#       m1[which(x$from==dataset2),col1] <- p_overlap2
#       
#       # fill m2 matrices
#       m2[which(x$from==dataset1),col2] <- overlap_with1
#       m2[which(x$from==dataset2),col1] <- overlap_with2
#     }
#     
#     # Merge data.frame and return result
#     m <- cbind(m2,m1) # merge matrices
#     res <- as.data.frame(m) # convert to data.frame
#     colnames(res) <- c(paste0(datasets,'_overlap'),paste0(datasets,'_p_overlap')) # rename columns of data.frame
#     res <- cbind(x,res) # bind the overlap data to data.frame x
#     res <- res[,colnames(res)!='new_id'] # remove the unique ids used by the function 
#     
#     return(res) 
#   })
#   
#   res <- do.call(rbind,res)
#   
#   temp <- as.matrix(res[,str_detect(colnames(res),'_p_overlap')])
#   temp[temp>overlap_threshold] <- 1
#   temp[temp<=overlap_threshold]  <- 0
#   res$overlap_group <- apply(temp,1,paste0,collapse='')
#   
#   
#   cat('Creating common site identifiers ... \n')
#   
#   # order dataframe by datasets and split into list by dataset
#   res <- res %>% arrange(factor(from, levels = datasets))
#   res <- split(res, res$from)
# 
#   # Create unique ids
#   res <- pblapply(datasets,function(dataset){
#     
#     # get dataset from list
#     dataset_df <- res[[dataset]]
#     
#     # find columns with site IDs for the other datasets
#     dataset_cols <- which(colnames(dataset_df)%in%paste0(datasets,'_overlap'))
#     other_dataset_ids <- dataset_cols[dataset_cols%notin%which(colnames(dataset_df)==paste0(dataset,'_overlap'))]
#     
#     # find column with site ID for the current dataset
#     current_dataset_ids <- which(colnames(dataset_df)=='ID')
#     
#     # order the dataset columns by the input order in the function input: "datasets"
#     dataset_order <- c(dataset,gsub('_.*','',colnames(dataset_df)[other_dataset_ids]))
#     dataset_cols <- setNames(c(current_dataset_ids,other_dataset_ids),dataset_order)
#     dataset_cols <- dataset_cols[match(datasets,names(dataset_cols))]
#     
#     # create unique id
#     dataset_df$site_id <- gsub(' ','',apply(dataset_df[,c(which(colnames(dataset_df)=='UniProtKB.AC'),dataset_cols)],1,paste,collapse='-'))
#     return(dataset_df)
#   })
# 
#   res <- do.call(rbind,res)
#   
#   cat('Done!\n')
#   return(res)
#   
# }

# sumlog function adapted from metap package 
# because it is not available for current R verison
# see: https://rdrr.io/cran/metap/src/R/sumlog.R
# sumlog <- function(p, log.p = FALSE, log.input = FALSE, only_return_p=T) {
#   if(log.input) {
#     keep <- p <= 0
#   } else {
#     keep <- (p > 0) & (p <= 1)
#   }
#   invalid <- sum(1L * keep) < 2
#   if(invalid) {
#     warning("Must have at least two valid p values")
#     res <- list(chisq = NA_real_, df = NA_integer_,
#                 p = NA_real_, validp = p[keep])
#   } else {
#     if(log.input) {
#       lnp <- p[keep] # already logged
#     } else {
#       lnp <- log(p[keep])
#     }
#     chisq <- (-2) * sum(lnp)
#     df <- 2 * length(lnp)
#     if(length(lnp) != length(p)) {
#       warning("Some studies omitted")
#     }
#     res <- list(chisq = chisq, df = df,
#                 p = pchisq(chisq, df, lower.tail = FALSE,
#                            log.p = log.p), validp = p[keep])
#   }
#   class(res) <- c("sumlog", "metap")
#   if(only_return_p){
#     res <- res$p
#   }
#   res
# }
# print.sumlog <- function(x, ...) {
#   cat("chisq = ", x$chisq, " with df = ", x$df, " p = ", x$p, "\n")
#   invisible(x)
# }
# 
# extract_overlapping_sites <- function(datasets,
#                                       keep_sites_from_dataset=NULL,
#                                       select_overlapping_datasets=NULL,
#                                       proteins=NULL,
#                                       pthresh=0.05,
#                                       overlap_threshold = 0,
#                                       if_within_set_overlap_to_1=T, 
#                                       useAdjusted=F, 
#                                       get_min=F, 
#                                       array_object=NULL){
#   
#   if(is.null(select_overlapping_datasets)){
#     select_overlapping_datasets <- paste(rep(1,length(datasets)), collapse = '')
#   }
#   
#   if(is.numeric(keep_sites_from_dataset)){
#     keep_sites_from_dataset <- datasets[keep_sites_from_dataset]
#   }
#   
#   if(is.null(array_object)){
#     array_object <- Find_array_object()
#   }
#   
#   cat('\nCalculating site overlaps ... \n')
#   sites     <- compare_datasets(datasets, binding = NULL, pthresh = pthresh, overlap_threshold = overlap_threshold, if_within_set_overlap_to_1=if_within_set_overlap_to_1, useAdjusted = useAdjusted, get_min = get_min, proteins = proteins, array_object = array_object)
#   
#   # Old implementation !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#   #cat('Calculating binding site overlaps ... \n') 
#   #binding     <- compare_datasets(datasets, binding = T, pthresh = pthresh, overlap_threshold = overlap_threshold, if_within_set_overlap_to_1=if_within_set_overlap_to_1, useAdjusted = useAdjusted, get_min = get_min, proteins = proteins, array_object = array_object)
#   #cat('Calculating non-binding site overlaps ... \n')
#   #non_binding <- compare_datasets(datasets, binding = F, pthresh = pthresh, overlap_threshold = overlap_threshold, if_within_set_overlap_to_1=if_within_set_overlap_to_1, useAdjusted = useAdjusted, get_min = get_min, proteins = proteins, array_object = array_object)
#   # !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#   
#   if(is.null(keep_sites_from_dataset)){
#     cat('Returning the overlapping sites from the dataset with the highest number of sites ... \n')
#     temp <- table(unique(sites[sites$overlap_group%in%select_overlapping_datasets,][,1:17])$from)
#     keep_sites_from_dataset <- names(temp[temp==max(temp)])
#     
#     if(length(keep_sites_from_dataset)>1){
#       continue <- NA
#       while(continue%notin%keep_sites_from_dataset){
#         continue <- readline('Please specify which dataset to extract ... default parameters not able to select optimal dataset: ')
#       }
#       keep_sites_from_dataset <- continue
#       cat('\n')
#     }
#   }
#   
#   cat('Extracting overlapping sites ... \n')
#   
#   # Old implementation !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#   # combined <- rbind(binding,non_binding) |> 
#   #   filter(overlap_group%in%select_overlapping_datasets) |> 
#   #   filter(from==keep_sites_from_dataset) |> 
#   #   group_by(site_id) |> 
#   #   filter(pval==min(pval)) |> 
#   #   filter(score==max(score)) |> 
#   #   ungroup() |> 
#   #   dplyr::select(c(1:16)) |> 
#   #   unique() 
#   # !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#   
#   sites <- sites |> 
#     filter(overlap_group%in%select_overlapping_datasets)
#   
#   sites_by_dataset <- split(sites[,c('score','pval','site_id')], sites$from)
#   
#   sites_signif <- sites_by_dataset[[1]]
#   for(i in 2:length(sites_by_dataset)){
#     sites_signif <- merge(x=sites_signif, y=sites_by_dataset[[i]], by ='site_id', all=T)
#   }
#   
#   sites_signif$average_score <- apply(sites_signif[,str_detect(colnames(sites_signif),'score')],1,mean)
#   cat("Running Fisher's combined probability test  ... \n")
#   sites_signif$fisher_pval   <- apply(sites_signif[,str_detect(colnames(sites_signif),'pval')],1,sumlog)
#   
#   sites_signif <- unique(sites_signif[,c('site_id','average_score','fisher_pval')])
#   sites <- merge(x=sites, y=sites_signif, by='site_id', all=T)
#   sites <- sites |> filter(from==keep_sites_from_dataset)
#   
#   sites <- sites |> 
#     mutate(highest_overlap = Reduce('+',sites[,str_detect(colnames(sites),'p_overlap')])) |> 
#     group_by(site_id) |> 
#     filter(highest_overlap==max(highest_overlap)) |> 
#     filter(fisher_pval==min(fisher_pval)) |> 
#     filter(average_score==max(average_score)) |> 
#     ungroup() |> 
#     dplyr::select(c('UniProtKB.AC','ID','start','end','average_score',
#                     'fisher_pval','signif','length','full_protein_sequence',
#                     'sequence','terminal_start','terminal_end','terminal',
#                     'terminal_length','distance_from_TM',
#                     'relative_distance_from_TM','from','highest_overlap')) |> 
#     mutate(site_id = paste0(UniProtKB.AC,'_',start,'_',end)) |> 
#     group_by(site_id) |> 
#     filter(highest_overlap==max(highest_overlap)) |> 
#     filter(fisher_pval==min(fisher_pval)) |> 
#     filter(average_score==max(average_score)) |> 
#     ungroup() |>
#     unique() 
#   
#   sites$signif <- sites$fisher_pval<pthresh
#   
#   cat("Final test for potential overlaps  ... \n")
#   
#   sites_list <- split(sites,sites$UniProtKB.AC)
#   
#   sites <- do.call(rbind,lapply(sites_list, function(site){
#     
#     rngs <- IRanges(site$start, site$end, names = site$site_id)
#     ol <- as.data.frame(IRanges::findOverlapPairs(rngs, rngs, type="within"))
#     ol <- ol[ol$first.names!=ol$second.names,]
#     
#     if(nrow(ol)>0){ # if overlaps are found
#       keep <- unique(unlist(lapply(split(ol,ol$second.names), function(id){
#         
#         if(min(site$fisher_pval[site$site_id %in% id$first.names]) <  min(site$fisher_pval[site$site_id %in% id$second.names])){
#           keep_temp <- id$first.names
#         } else {
#           keep_temp <- id$second.names
#         }
#         return(keep_temp)
#       }), use.names = F))
#       
#       site <- site[site$site_id%in%keep,]
#     }
#     
#     return(site)
#     
#   }))
#   
#   mistakes <- sites$UniProtKB.AC[sites$signif==T & sites$length>100]
#   
#   if(length(mistakes)>0){
#     cat('There seems to be a large mismatch between traces for:\n')
#     print(mistakes)
#     cat('Signif has been set to FALSE for these sites!')
#   }
#   
#   sites$signif[sites$signif==T & sites$length>100] <- F
#   
#   cat('Done!\n')
#   cat('-----------------------------------------------------------\n')
#   return(sites)
# }
# 
# # Function for plotting venndiagram of comparison
# plot_comparison_venndiagram <- function(dataset_comparison,
#                                         colors = c("red","steelblue"),
#                                         alpha = 0.5,
#                                         lty = 1,
#                                         plot = T){
#   
#   datasets <- gsub('_p_overlap','',colnames(dataset_comparison)[str_detect(colnames(dataset_comparison),'_p_overlap')])
#   
#   euler_input <- c()
#   for(i in datasets){
#     n <- table(dataset_comparison$overlap_group[dataset_comparison$from==i])
#     
#     m <- do.call(rbind,str_split(names(n), pattern = ''))
#     colnames(m) <- datasets
#     
#     combs <- apply(m,1,function(x){paste(names(x)[x=='1'],collapse='&')})
#     names(n) <- combs
#     euler_input <- c(euler_input,n)
#   }
#   
#   euler_input <- euler_input[unique(names(euler_input))]
#   
#   p <- plot(euler(euler_input),
#             quantities=TRUE,
#             fills =list(fill=colors,alpha=alpha),
#             lty=lty)
#   if(plot==T){
#     print(p)
#   }
#   return(p)
# }


#-############################################################################-#
#-########################### CONVENIENCE FUNCTIONS ##########################-#
#-############################################################################-#

# identify variables in environment
.my_env <- environment() 

# save_array <- function(array_object, 
#                        filepath=NULL, 
#                        filename=NULL){ #remember to save the output of this function to an array_object (best to overwrite loaded array_object) else it will print everything in the console and your object dates will not be updated in the open array_object
#   
#   overwrite <- F
#   if(is.null(filepath)){
#     filepath <- here('data')
#   }
#   
#   if(is.null(filename)){
#     filename <- 'Array_object'
#     overwrite <- T
#   }
#   
#   if(overwrite==T){
#     continue <- NA
#     while(continue%notin%c('y','Y','ye','YE','Ye','yes','YES','Yes','YeS','YEs','yeS','yES','n','N','no','No','NO')){
#       continue <- readline('!!!! You are trying to overwrite the current Array_object !!!! \nAre you sure you want to do this? (options: "yes" or "no") \n')
#     }
#     if(continue%in%c('n','N','no','No','NO')){
#       overwrite <- F
#       stop()
#     }
#     cat('\n')
#   }
#   
#   if(!str_detect(filename,'\\.')){
#     filename <- paste0(filename,'.rds')
#   }
#   
#   if(!str_detect(filepath,'\\.')){
#     fullpath <- paste0(filepath,'/',filename)
#   } else {
#     fullpath <- filepath
#   }
#   
#   updated_array_object <- update_date(array_object,saved = T)
#   
#   rdata_file = file(fullpath, blocking = TRUE)
#   saveRDS(updated_array_object, file=rdata_file)
#   close(rdata_file)
#   
#   about_array_object(updated_array_object)
#   
#   return(updated_array_object)
# }
# 
# create_array_backup <- function(array_object){
#   
#   if(!dir.exists(here('backup'))){
#     dir.create(here('backup'))
#   }
# 
#   name <- deparse(substitute(array_object))
#   
#   now <- gsub(' ','_',date())
#   colons <- str_locate(now,':')[,'start']
#   substr(now, colons, colons) <- 'h'
#   substr(now, colons+3, colons+3) <- 'm'
#   now <- gsub(substr(now,colons-2,colons+6), paste0(substr(now,colons-2,colons+5),'s_'),now) 
#   
#   foldername <- paste0(name,'_',now)
#   
#   dir.create(here('backup',foldername))
#   
#   save_array(array_object, here('backup',foldername),name)
# 
# }

Array_Class_Filter <- function(x){inherits(get(x,envir=.my_env),'array_object')}

Find_array_object <- function(...){
  array_object_names <- Filter(Array_Class_Filter, ls(name = .my_env))
  if(length(array_object_names)==1){
    object_found <- get(array_object_names, envir = .my_env)  
  } else{
    print("Multiple or no objects with class 'array_object' loaded, try to specify which one to use ...")
    stop()
  }
  return(object_found)
}

update_date <- function(array_object=NULL,
                        saved=F){
  if(is.null(array_object)){
    array_object <- Find_array_object()
  }
  if(saved==F){
    array_object$About$Last_updated <- date()
  }
  if(saved==T){
    array_object$About$Last_saved <- date()
  }
  
  return(array_object)
}

strip_terminal_specifier <- function(ArrayID){
  gsub('_.*','',ArrayID)
}

detect_identifier <- function(id){
  logi <- str_detect(id,'_')
  if(sum(logi)==length(id) | sum(logi)==0){
    id_type <- 'ArrayID'
  } else {
    print('Your identifiers are not all of the same type, please revise them...')
    stop()
  } 
  if(sum(logi)==0){
    id_type <- 'UniProtKB.AC'
  } 
  return(id_type)
}

find_associated_terminals <- function(proteins,
                                      array_object=NULL){
  
  if(is.null(array_object)){
    array_object <- Find_array_object()
  }
  
  if(!is.null(proteins)){
    id_type <- detect_identifier(proteins)
    
    if(id_type=="UniProtKB.AC"){
      proteins <- names(get_terminals(array_object,proteins))
    } 
  }
  return(proteins)
}

# sample_dataset <- function(dataset,
#                            sample_size,
#                            bootstraps=100000,
#                            array_object=NULL,
#                            verbose=T){
#   
#   if(is.null(array_object)){
#     array_object <- Find_array_object()
#   }
#   
#   normalised_scores <- na.omit(get_dataset(array_object,dataset)$normalised)
#   
#   if(verbose==T){
#     cat(paste0("\nRunning ",bootstraps," bootstrap samplings from dataset '",dataset,"' ... \n"))
#     sampling <- pblapply(1:bootstraps,function(x){
#       sample(normalised_scores,sample_size,replace = T)
#     })
#   } else {
#     sampling <- lapply(1:bootstraps,function(x){
#       sample(normalised_scores,sample_size,replace = T)
#     })
#   }
#   
#   return(sampling)
# }
# 
# extract_controls <- function(array_object,
#                              normalised_data,
#                              dataset,
#                              column_with_peptide_sequences){
#   
#   protein_peptides <- get_dataset(array_object,dataset)[,column_with_peptide_sequences]
#   control_peptides <- normalised_data[normalised_data[,column_with_peptide_sequences]%notin%protein_peptides,]
#   
#   return(control_peptides)
# }
# 
# density_function <- function(x,
#                              model,
#                              fit,
#                              log=F){
#   den_func <- get(paste0('d',model))
#   
#   params <- fit$parameters
#   
#   if(length(params)==2){
#     d <- den_func(x,fit[[params[1]]],fit[[params[2]]], log=log)
#   }
#   
#   if(length(params)==3){
#     d <- den_func(x,fit[[params[1]]],fit[[params[2]]],fit[[params[3]]], log=log)
#   }
#   
#   if(length(params)==4){
#     d <- den_func(x,fit[[params[1]]],fit[[params[2]]],fit[[params[3]]],fit[[params[4]]], log=log)
#   }
#   
#   if(length(params)>4 | length(params)==1){
#     print('model has more or less parameters than accepted... min parameters = 1 & max parameters = 4')
#   }  
#   return(d)
# }
# 
# random_generation_function <- function(n,
#                                        model,
#                                        fit){
#   rnd_func <- get(paste0('r',model))
#   
#   params <- fit$parameters
#   
#   if(length(params)==2){
#     r <- rnd_func(n,fit[[params[1]]],fit[[params[2]]])
#   }
#   
#   if(length(params)==3){
#     r <- rnd_func(n,fit[[params[1]]],fit[[params[2]]],fit[[params[3]]])
#   }
#   
#   if(length(params)==4){
#     r <- rnd_func(n,fit[[params[1]]],fit[[params[2]]],fit[[params[3]]],fit[[params[4]]])
#   }
#   
#   if(length(params)>4 | length(params)==1){
#     print('model has more or less parameters than accepted... min parameters = 1 & max parameters = 4')
#   }  
#   return(r)
# }
# 
# extract_binding_site_parameters <- function(calculation_id){
#   info <- str_split(calculation_id,pattern = '_')[[1]][c(2,5,8)]
#   info[1] <- ifelse(nchar(info[1]<6),paste0(info[1],paste(rep(' ',6-nchar(info[1])),collapse='')),info[1])
#   info[2] <- ifelse(info[2]=='FALSE','NO  ','YES ')
#   info[3] <- str_to_title(info[3])
#   return(info)
# }
# 
# create_SWaFi_trace_mask <- function(datasets = c('FBBE1','FBBE2'),
#                                     proteins=NULL,
#                                     array_object=NULL,
#                                     pthresh=NULL,
#                                     get_min=F, ...){
#   
#   if(is.null(array_object)){
#     array_object <- Find_array_object()
#   }
#   
#   sites  <- get_binding_sites(proteins = proteins, binding = NULL, confidence = pthresh, datasets = datasets, array_object = array_object, ...)
#   terminals <- split(sites,sites$UniProtKB.AC)
#   
#   cat('Creating masks ... \n ')
#   mask <- pblapply(terminals, function(terminal){
#     terminal   <- as.data.frame(terminal)
#     SWaFi_mask <- data.frame(position=unique(terminal$terminal_start):unique(terminal$terminal_end), mask=rep(0,terminal$terminal_length[1]))
#     
#     terminal <- terminal[terminal$signif==T,]
# 
#     if(nrow(terminal)>0){
#       keep <- unlist(apply(terminal,1,function(x){as.numeric(x[3]):as.numeric(x[4])}))
#       SWaFi_mask$mask[SWaFi_mask$position%in%keep] <- 1
#     }
#     
#     return(SWaFi_mask)
#     
#   })
#   return(mask)
# }

mask_proteins_in_SWaFi_matrix <- function(SWaFi_matrix,
                                          proteins){
  SWaFi_matrix[rownames(SWaFi_matrix)%in%proteins & !is.na(SWaFi_matrix)] <- 0
  return(SWaFi_matrix)
}

species_filter <- function(proteins=NULL,
                           taxid=9606,
                           array_object=NULL,
                           df=NULL,
                           id_column=NULL){
  
  if(!is.null(proteins) & !is.null(df)){
    cat('Be sure that you want to filter the dataframe based on the provided proteins as well as taxid ... \n You should check if you get the desired result, the function was not meant for this\n ')
  }
  
  taxid <- as.character(taxid)
  species <- unlist(get_metadata(array_object,'Organism.ID',proteins))
  ids <- names(species)[species%in%taxid] 
  
  if(!is.null(proteins)){
    type <- detect_identifier(proteins)
    
    if(type=="UniProtKB.AC"){
      ids <- unique(strip_terminal_specifier(ids))
    }
  }
  
  if(!is.null(df)){
    if(is.null(id_column)){cat('No id_column provided ... trying with default: UniProtKB.AC\n')
      id_column <- 'UniProtKB.AC'
      if(!id_column %in% colnames(df)){
        cat('Default did not work ... exiting\n')
        stop()
      }
    }
    
    ids <- df[df[,id_column]%in%ids,]
    
  }
  
  return(ids)
}

find_SWaFi_matrix_TM_segment <- function(SWaFi_matrix){
  columns <- which(str_detect(colnames(SWaFi_matrix),'TM')==T) 
  columns <- list(start=min(columns),end=max(columns),columns=columns)
  return(columns)
}

subset_SWaFi_matrix <- function(SWaFi_matrix,
                                proteins){
  if(detect_identifier(proteins)=='ArrayID'){
    proteins <- strip_terminal_specifier(proteins)
  }
  SWaFi_matrix <- SWaFi_matrix[rownames(SWaFi_matrix)%in%proteins,]
  return(SWaFi_matrix)                                       
}

set_SWaFi_matrix_width <- function(SWaFi_matrix,
                                   N_aa_length=600, 
                                   C_aa_length=600){
  
  N_aa_length <- N_aa_length/2
  C_aa_length <- C_aa_length/2
  
  TM_segment <- find_SWaFi_matrix_TM_segment(SWaFi_matrix)
  range <- c(TM_segment$start-N_aa_length-1,TM_segment$end+C_aa_length+1)
  
  if(range[1]<1){
    add_to_N <- matrix(NA,ncol=abs(range[1])+1,nrow=nrow(SWaFi_matrix))
    new_range <- c(1,range[2]+ncol(add_to_N))
    cat('Note: empty columns have been added to the N-terminus!\n')
  } else {
    add_to_N <- NULL
    new_range <- range
    }
  
  if(range[2]>ncol(SWaFi_matrix)){
    add_to_C <- matrix(NA,ncol=range[2]-ncol(SWaFi_matrix),nrow=nrow(SWaFi_matrix))
    range <- c(new_range[1],new_range[2]+ncol(add_to_C))
    cat('Note: empty columns have been added to the C-terminus!\n')
  } else {add_to_C <- NULL}
  
  if(!is.null(add_to_N)){
    SWaFi_matrix <- cbind(add_to_N,SWaFi_matrix)
  }
  if(!is.null(add_to_C)){
    SWaFi_matrix <- cbind(SWaFi_matrix,add_to_C)
  }
  
  return(SWaFi_matrix[,c(range[1]:range[2])])
  
}

is.outdated <- function(UniProtKB.AC, array_object=NULL){
  
  if(is.null(array_object)){
    array_object <- Find_array_object()
  }
  
  res <- unlist(lapply(UniProtKB.AC, function(x){
    quiet(get_metadata(array_object, 'status', x)[[1]]) == 'outdated'
  }))
  
  return(setNames(res,UniProtKB.AC))
  
}

outdated <- function(UniProtKB.AC, # function changes status of proteins to outdated
                     array_object=NULL){
  if(is.null(array_object)){
    array_object <- Find_array_object()
  }
  
  terminals <- get_terminals(array_object)
  
  index <- which(strip_terminal_specifier(names(terminals))==UniProtKB.AC)
  
  for(i in index){
    terminals[[i]]$status <- 'outdated'
  }
  
  array_object$Protein_terminals <- terminals
  
  updated_array_object <- update_date(array_object)
  
  return(updated_array_object)
}

filter_outdated <- function(proteins,
                            array_object=NULL){
  
  if(is.null(array_object)){
    array_object <- Find_array_object()
  }
  
  terminals <- get_terminals(array_object, proteins)
  
  status <- unlist(lapply(terminals,function(x){
    x$status
  }))
  
  status <- status[status=='outdated']
  status <- names(status)
  
  if(detect_identifier(proteins)=='UniProtKB.AC'){
    status <- unique(strip_terminal_specifier(status))
  }
  
  proteins <- proteins[proteins%notin%status]
  return(proteins)
}

extract_part_of_sequence <- function(proteins,
                                     from,
                                     to,
                                     array_object=NULL){
  
  if(is.null(array_object)){
    array_object <- Find_array_object()
  }
  
  sequence <- get_metadata(array_object, colnames='Full.Sequence', proteins = proteins)
  
  sub_seqs <- lapply(sequence,substr,start=from,stop=to)
  
  return(sub_seqs)
  
}

#-############################################################################-#
#-############## GET FUNCTIONS ('get' items from array_object) ###############-#
#-############################################################################-#
about_array_object <- function(array_object=NULL) {

  if(is.null(array_object)){
    array_object <- Find_array_object()
  }
  
  all_proteins <- get_available_proteins(array_object)
  outdated <- all_proteins[all_proteins %notin% quiet(filter_outdated(all_proteins,Array_object))]
  
  cat("\n-----------------------------------------------------------",
      "\n Authors -", array_object$About$Authors, 
      "\n-----------------------------------------------------------",
      "\n Date created:           ", as.character(array_object$About$Date_created), 
      "\n Last updated:           ", as.character(array_object$About$Last_updated),
      "\n Last saved:             ", as.character(array_object$About$Last_saved),
      "\n # of proteins:          ", length(all_proteins),
      "\n # of outdated proteins: ", paste0(length(outdated),' (', paste(outdated, collapse = ' '),')'),
      "\n # of terminals:         ", length(get_terminals(array_object)), 
      "\n NCBI taxonomic IDs:     ", paste(unique(unlist(get_metadata(array_object,'Organism.ID'))), collapse = ' '),
      "\n # of Experiments:       ", length(array_object$Experiments),
      "\n Experiment IDs:         ", names(array_object$Experiments),
      "\n Affiliation:            ", array_object$About$Affiliation, "\n")

}

about_experiment <- function(dataset=NULL,
                             array_object=NULL){
  
  if(is.null(array_object)){
    array_object <- Find_array_object()
  }
  
  if(is.null(dataset)){
    dataset <- get_dataset_ids(array_object)
  }
  
  for(i in dataset){
    experiment <- get_experiment(array_object,i)

    if('SWaFi_binding_sites'%in%names(experiment)){
      sites <- names(experiment$SWaFi_binding_sites)
      
      sites <- unlist(sapply(sites,function(x){
        info <- extract_binding_site_parameters(x)
        paste0('\n Significance threshold: ',info[1],' | p-values adjusted: ',info[2],' | ',info[3],' sites have been extracted')
      }))
      
      sites <- paste(c(sites,'\n'),collapse = '')
      
    } else {
      sites <- NULL
    }
    
    spacer <- '\n------------------------------------------------------------------------------------------------------------------'
    
    cat(paste0("\n Dataset '",i,"' was added on the ",experiment$Added,
               "\n From the file: ", gsub('.*/','',experiment$From_file), spacer, ifelse('Imputation_fit'%in%names(experiment),
               paste0("\n Missing values were imputed with the function '",get_imputation_information(array_object,i,verbose = F)$name,"', with the following parameters:",
               "\n ",paste(riffle(paste0(names(get_imputation_information(array_object,i,verbose = F)$parameters),':'),get_imputation_information(array_object,i,verbose = F)$parameters), collapse=' ')),'\n No function for imputation of missing values has been provided'), spacer,
               "\n SWaFi traces ", ifelse('SWaFi_traces'%in%names(experiment),'have', 'have NOT')," been calculated", 
               ifelse('SWaFi_traces_parameters'%in%names(experiment),paste0(' keeping ',experiment$SWaFi_traces_parameters$nthresh*100,'% of the frequencies when smoothing'),''),
               "\n Background distributions ", ifelse('Fits'%in%names(experiment),'have', 'have NOT')," been fitted", ifelse('Fits'%in%names(experiment),paste0(' with ',experiment$Fits[[1]]$family[2]),''),
               "\n SWaFi regions ", ifelse('SWaFi_regions'%in%names(experiment),'have', 'have NOT')," been determined",
               ifelse('SWaFi_regions_parameters'%in%names(experiment),paste0(' with a threshold of ',experiment$SWaFi_regions_parameters$lcut*100,'% for merging regions'),''),spacer,
               ifelse('none'==experiment$SWaFi_regions_parameters$pAdjust,paste0('\n p-values have NOT been corrected',spacer), paste0('\n p-values have been corrected using: ',experiment$SWaFi_regions_parameters$pAdjust,spacer)),"\n",
               ifelse('SWaFi_binding_sites'%in%names(experiment),'', '\n SWaFi binding sites have NOT been calculated'),
               ifelse('SWaFi_binding_sites'%in%names(experiment),"\n SWaFi binding sites with the following parameters have been calculated:","\n"),
               ifelse(length(sites)>0,sites,''),spacer,'\n \n'))
    
  }
  
}

count_binding_sites <- function(dataset,
                                array_object=NULL){
  if(is.null(array_object)){
    array_object <- Find_array_object()
  }
  
  site_list <- get_experiment(array_object,dataset)$SWaFi_binding_sites
  spacer <- '------------------------------------------------------------------------------------------------------------------\n'
  cat(paste0("Dataset '",dataset,"' binding site summary statistics\n",spacer))
  for(i in names(site_list)){
    sites <- rbindlist(site_list[[i]])
    
    no_sites <- nrow(sites[sites$signif==T,])
    no_proteins <- length(unique(strip_terminal_specifier(sites[sites$signif==T,]$UniProtKB.AC)))
    info <- extract_binding_site_parameters(i)
    
    cat(paste0('Significance threshold: ',info[1],' | p-values adjusted: ',info[2],' | ',info[3],' sites have been extracted\n',
               'With these settings ',no_sites,' binding sites are detected in ',no_proteins,' proteins. \n',spacer))
    
  }
  
}

get_dataset <- function(array_object=NULL,
                        dataset=NULL,
                        proteins=NULL){
  
  if(is.null(dataset)){
    print('Please specify which dataset you want to get ...')
    stop()
  }
  
  if(is.null(array_object)){
    array_object <- Find_array_object()
  }
  
  terminals <- get_terminals(array_object, proteins = proteins)
  df <- do.call(rbind,lapply(terminals,function(x){
    x[[dataset]]
  }))
  
  df <- rownames_to_column(df)
  df[,1] <- gsub('\\..*','',df[,1])
  return(df)
}

get_metadata <- function(array_object=NULL, 
                         colnames=NULL, 
                         proteins=NULL){
  
  if(is.null(array_object)){
    array_object <- Find_array_object()
  }
  
  if(is.null(colnames)){
    colnames <- get_colnames(array_object)
  }
  
  terminals <- get_terminals(array_object,proteins)
  lapply(terminals,function(x){
    x <- x[colnames]
    if(length(colnames)>1){
      if(sum(unlist(lapply(x,is.data.frame)))==0){
        x <- as.data.frame(matrix(unlist(x), ncol=length(colnames)))
        colnames(x) <- colnames
      }
    } else {
      x <- x[[colnames]]
    }
    return(x)
  })
}

get_colnames <- function(array_object=NULL){
  
  if(is.null(array_object)){
    array_object <- Find_array_object()
  }
  
  terminals <- array_object$Protein_terminals
  sort(unique(unlist(lapply(terminals,function(x){
    names(x)
  }))))
}

get_dataset_ids <- function(array_object=NULL){
  
  if(is.null(array_object)){
    array_object <- Find_array_object()
  }
  
  names(array_object$Experiments)
}

get_available_proteins <- function(array_object=NULL){
  unique(strip_terminal_specifier(names(quiet(get_terminals(array_object)))))
}

get_name_from_uniprot <- function(uniprotID=NULL, 
                                  array_object=NULL){
  
  if(is.null(array_object)){
    array_object <- Find_array_object()
  }
  
  if(is.null(uniprotID)){
    uniprotID <- get_available_proteins(array_object)
  }
  
  Name <- get_metadata(array_object, colnames = 'Protein.name', proteins=uniprotID)

  if(detect_identifier(uniprotID)=='UniProtKB.AC'){
    
    Name <- unlist(Name)
    names(Name) <- strip_terminal_specifier(names(Name))
    Name <- Name[unique(names(Name))]
    Name <- Name[order(factor(names(Name), levels = uniprotID))]
    
    res <- unique(tibble(UniProtKB.AC=uniprotID,ProteinName=Name))
    
  } else {
    Name <- Name[uniprotID]
    res <- unique(tibble(UniProtKB.AC=uniprotID,ProteinName=unlist(Name)))
  }
  return(res)
}

get_uniprot_from_name <- function(search, 
                                  array_object=NULL){
  
  if(is.null(array_object)){
    array_object <- Find_array_object()
  }
  
  keys  <- get_metadata(array_object, colnames=c('Protein.name','Genes'))
  search <- tolower(search)
  cat("\nSearching ... \n")
  res <- unlist(pblapply(keys,function(x){
    x <- unlist(paste(x,collapse=' '))
    str_detect(tolower(x),search)
  }))
  
  res <- names(res[res==T])
  
  if(length(res)==0){
    cat("\nCouldn't find a match please try with another search term, e.g. gene name ...\n")
    stop()
  }
  
  cat("\nResult:\n")
  
  res <- do.call(rbind,get_metadata(array_object, 
                                    colnames=c('Protein.name','Genes'), 
                                    proteins = res))
  res <- as_tibble(rownames_to_column(res,var = 'UniProtKB.AC'))
  
  return(res)
}

get_gene_from_uniprot <- function(uniprotID=NULL, 
                                  only_first=F, 
                                  array_object=NULL){
  
  if(is.null(array_object)){
    array_object <- Find_array_object()
  }
  
  if(is.null(uniprotID)){
    uniprotID <- get_available_proteins(array_object)
  }
  
  Name <- unlist(get_metadata(array_object, colnames = 'Genes', proteins=uniprotID))
  
  if(detect_identifier(uniprotID)=='UniProtKB.AC'){
    names(Name) <- strip_terminal_specifier(names(Name))
    Name <- Name[unique(names(Name))]
  }
  
  if(only_first==T){
    Name <- gsub(' .*','',Name)
  }
  
  return(unique(tibble(UniProtKB.AC=names(Name),Gene=Name)))
}

get_terminals <- function(array_object=NULL, 
                          proteins=NULL){
  if(is.null(array_object)){
    array_object <- Find_array_object()
  }
  
  terminals <- array_object$Protein_terminals
  
  if(!is.null(proteins)){
    id_type <- detect_identifier(proteins)
    
    if(id_type=="UniProtKB.AC"){
      print('Returning all available terminals for the specified proteins...')
      proteins_in_array_object <- strip_terminal_specifier(names(terminals)) 
    } else {
      proteins_in_array_object <- names(terminals)
    }
    terminals <- terminals[which(proteins_in_array_object%in%proteins)]
  }
  
  return(terminals)
}

get_experiment <- function(array_object=NULL, 
                           dataset=NULL){
  
  if(is.null(dataset)){
    print('Please specify which dataset you want  ...')
    stop()
  }
  
  if(is.null(array_object)){
    array_object <- Find_array_object()
  }
  
  return(array_object$Experiments[[dataset]])
  
}

get_random_samplings <- function(array_object=NULL,
                                 dataset=NULL){
  
  if(is.null(dataset)){
    print('Please specify which dataset you want ...')
    stop()
  }
  
  if(is.null(array_object)){
    array_object <- Find_array_object()
  }
  
  return(array_object$Experiments[[dataset]]$Random_Sampling)
  
}

get_path_to_experiment_file <- function(dataset=NULL,
                                        array_object=NULL){
  
  if(is.null(array_object)){
    array_object <- Find_array_object()
  }
  
  if(is.null(dataset)){
    datasets <- get_dataset_ids(array_object)
    paths <- unlist(lapply(datasets,function(x){
      get_experiment(array_object,x)$From_file
    }))
    names(paths) <- datasets
  } else {
    paths <- get_experiment(array_object,dataset)$From_file
  }
  return(paths)
}

# get_imputation_information <- function(array_object=NULL, 
#                                        dataset=NULL, 
#                                        verbose=T){
#   
#   if(is.null(array_object)){
#     array_object <- Find_array_object()
#   }
#   
#   experiment <- get_experiment(array_object, dataset)
#   
#   if('Imputation_fit' %notin% names(experiment)){
#     print(paste0("Error, fit for imputation does not exist in experiment '",dataset,"' ... please add a fit to continue"))
#     stop()
#   }
#   
#   impute_distribution <- experiment$Imputation_fit$family[1]
#   if(verbose==T){
#     cat(paste0('Function used for imputation of missing values: ',paste0('r',impute_distribution),'\n \n'))
#   }
#   impute_distribution_function <- get(paste0('r',impute_distribution))
#   
#   impute_fit <- experiment$Imputation_fit
#   impute_parameters <- sapply(impute_fit$parameters,function(x){impute_fit[[x]]})
# 
#   return(list(Function=impute_distribution_function, parameters=impute_parameters, name=paste0('r',impute_distribution)))
#   
# }

# get_fits <- function(array_object=NULL, 
#                      dataset=NULL, 
#                      verbose=T){
#  
#   if(is.null(array_object)){
#     array_object <- Find_array_object()
#   }
#    
#   experiment <- get_experiment(array_object, dataset)
#   
#   if('Fits' %notin% names(experiment)){
#     print(paste0("Error, fits do not exist in experiment '",dataset,"' ... please add fits to continue"))
#     stop()
#   }
#   
#   distribution <- experiment$Fits[[1]]$family[1]
#   if(verbose==T){
#     cat(paste0('Function used for calculating p-values: ',paste0('p',distribution),'\n'))
#   }
#   distribution_function <- get(paste0('p',distribution))
#   
#   fits <- experiment$Fits
#   
#   return(list(Fits=fits,Function=distribution_function,name=paste0('p',distribution)))
#   
# }

get_SWaFi_traces <- function(array_object=NULL, 
                             dataset=NULL, 
                             proteins=NULL, 
                             parameters=F){
  
  if(is.null(array_object)){
    array_object <- Find_array_object()
  }
  
  experiment <- get_experiment(array_object, dataset)
  
  if('SWaFi_traces' %notin% names(experiment)){
    print(paste0("Error, SWaFi traces do not exist in experiment '",dataset,"' ... please add traces to continue"))
    stop()
  }
  
  traces <- experiment$SWaFi_traces
  if(!is.null(proteins)){
    id_type <- detect_identifier(proteins)
    
    if(id_type=="UniProtKB.AC"){
      print('Returning all available terminals for the specified proteins...')
      proteins_in_array_object <- strip_terminal_specifier(names(traces)) 
    } else {
      proteins_in_array_object <- names(traces)
    }
    traces <- traces[which(proteins_in_array_object%in%proteins)]
  } 
  if(parameters==T){
    parameters <- experiment$SWaFi_traces_parameters
    traces <- list(SWaFi_traces=traces,SWaFi_parameters=parameters)
  }
 return(traces)
}

# get_SWaFi_regions <- function(array_object=NULL, 
#                               dataset=NULL, 
#                               proteins=NULL,
#                               parameters=F){
#   
#   if(is.null(array_object)){
#     array_object <- Find_array_object()
#   }
#   
#   experiment <- get_experiment(array_object, dataset)
#   
#   if('SWaFi_regions' %notin% names(experiment)){
#     print(paste0("Error, SWaFi regions do not exist in experiment '",dataset,"' ... please define regions to continue"))
#     stop()
#   }
#   
#   regions <- experiment$SWaFi_regions
#   if(!is.null(proteins)){
#     id_type <- detect_identifier(proteins)
#     
#     if(id_type=="UniProtKB.AC"){
#       print('Returning all available terminals for the specified proteins...')
#       proteins_in_array_object <- strip_terminal_specifier(names(regions)) 
#     } else {
#       proteins_in_array_object <- names(regions)
#     }
#     regions <- regions[which(proteins_in_array_object%in%proteins)]
#   }
#   if(parameters==T){
#     parameters <- experiment$SWaFi_regions_parameters
#     regions <- list(SWaFi_traces=regions,SWaFi_parameters=parameters)
#   }
#   return(regions)
# }

get_SWaFi_binding_sites <- function(array_object=NULL,
                                    dataset,proteins=NULL,
                                    pthresh=NULL,
                                    useAdjusted=F,
                                    get_min=F){
  
  if(is.null(array_object)){
    array_object <- Find_array_object()
  }
  
  proteins <- find_associated_terminals(proteins,array_object)
  
  SWaFi_binding_sites <- get_experiment(array_object,dataset)$SWaFi_binding_sites
  
  if(is.null(pthresh)){
    cat('Returning sites with all available thresholds ... \n \n')
    sites <- SWaFi_binding_sites
    
    if(!is.null(proteins)){
      sites <- lapply(sites,function(x){x[proteins]})
    }
    
  } else {
    Type <- ifelse(get_min==F,'widest','lowest')
    calculation_id <- paste0('Threshold_',pthresh,'__Adjusted_',useAdjusted,'__Type_',Type)
    
    if(calculation_id%notin%names(SWaFi_binding_sites)){
      print("'SWaFi_binding_sites' with these parameters do not exist, calculate them with 'add_SWaFi_binding_sites' ...")
      stop()
    } 
    sites <- SWaFi_binding_sites[[calculation_id]]
    
    if(!is.null(proteins)){
      sites <- sites[proteins]
    }
  }

  return(sites)
  
}

get_SWaFi_binding_sites_long <- function(array_object=NULL,
                                         dataset,proteins=NULL,
                                         pthresh,useAdjusted=F,
                                         get_min=F){
  
  if(is.null(array_object)){
    array_object <- Find_array_object()
  }
  
  experiment <- get_experiment(array_object,dataset)
  Type <- ifelse(get_min==F,'widest','lowest')
  calculation_id <- paste0('Threshold_',pthresh,'__Adjusted_',useAdjusted,'__Type_',Type)
  
  if(calculation_id%notin%names(experiment$SWaFi_binding_sites)){
    print("'SWaFi_binding_sites' with these parameters do not exist, calculate them with 'add_SWaFi_binding_sites' ...")
    stop()
  } 
  
  sites <- experiment$SWaFi_binding_sites[[calculation_id]]
  
  if(!is.null(proteins)){
    proteins <- find_associated_terminals(proteins,array_object)
    sites <- sites[proteins]
  }
  
  sites <- lapply(sites,pivot_binding_site_longer)
  ids <- names(sites)
  
  sites <- lapply(ids,function(id){
    site <- sites[[id]]
    trace <- experiment$SWaFi_traces[[id]]
    
    merged <- merge(trace,site,by='position',all.x=T)
    return(merged)
  })
  
  names(sites) <- ids
  
  return(sites)
  
}

get_SWaFi_binding_site_details <- function(dataset,
                                           pthresh,
                                           useAdjusted=F,
                                           get_min=F,
                                           proteins=NULL,
                                           array_object=NULL,
                                           binding=T, 
                                           extract_amphipathic_helices=F, 
                                           pH=7.4,
                                           cassette_len=11){
  
  if(is.null(array_object)){
    array_object <- Find_array_object()
  }
  
  proteins <- find_associated_terminals(proteins,array_object)
  
  sites <- do.call(rbind,get_SWaFi_binding_sites(array_object,dataset,proteins,pthresh,useAdjusted,get_min))
  sites$length <- sites$end-sites$start+1
  
  TM_start_end  <- rownames_to_column(do.call(rbind,get_metadata(array_object,colnames = c('sub.start','sub.end','terminal'))))  
  colnames(TM_start_end)[1] <- 'UniProtKB.AC'
  TM_start_end$sub.start <- as.numeric(TM_start_end$sub.start)
  TM_start_end$sub.end   <- as.numeric(TM_start_end$sub.end)
  names(TM_start_end)[2:3] <- c('terminal_start','terminal_end')
  
  sequences <- get_metadata(array_object,colnames = 'Full.Sequence',proteins=proteins)
  sequences <- rownames_to_column(as.data.frame(unlist(sequences)))
  names(sequences) <- c('UniProtKB.AC','full_protein_sequence')  
  sites <- merge(sites,sequences,by='UniProtKB.AC',all.x=T)
  sites$sequence <- substr(sites$full_protein_sequence, start=sites$start, stop=sites$end)
  
  sites <- merge(sites,TM_start_end,by='UniProtKB.AC',all.x=T)
  sites$terminal_length <- sites$terminal_end-sites$terminal_start+1 
  
  sites$distance_from_TM <- ifelse(sites$terminal=='C',sites$start-sites$terminal_start+1,sites$terminal_end-sites$end+1)
  sites$relative_distance_from_TM <- sites$distance_from_TM/sites$terminal_length
  
  if(!is.null(binding)){
    if(binding==T){
      sites <- sites[sites$signif==T,]
    } else {
      sites <- sites[sites$signif==F,]
    }
  }
  
  sites <- na.omit(sites)
  
  if(extract_amphipathic_helices==T){
    cat(paste0('Extracting best Amphipathic Helices using parameters cassette_len=',cassette_len,', and pH=',pH,':\n'))
    sites <- add_helical_confidence_to_SWaFi_binding_sites(sites, cassette_len = cassette_len, pH = pH, array_object = array_object)
    cat('Done!\n')
  }
  return(sites)
}

get_SWaFi_score <- function(array_object=NULL,
                            dataset=NULL,
                            proteins=NULL,
                            pthresh=NULL,
                            cut=NA,
                            get_min=F,
                            lcut=NULL,
                            pAdjust=NULL,
                            Plot=F){
  
  # Find array_object if missing -----------------------------------------------
  if(is.null(array_object)){
    array_object <- Find_array_object()
  } 
  
  # Messages and Errors --------------------------------------------------------
  cat(paste0('Calculating Smoothed Weighted average Fluorescence intensity scores\nDataset: ',dataset,'\n \n'))
  
  if(class(Array_object)!='array_object'){
    cat('Object is not of class "array_object"... \n  Process aborted \n')
    stop()
  }
  
  if(is.null(dataset)){
    print('Please specify which dataset you want to add SWaFi scores to  ...')
    stop()
  }
  
  if(is.null(pthresh)){
    print(" 'pthresh' not specified ... defaulting to significance threshold of 0.05 ... ")
    pthresh <- 0.05
  }
  
  if(is.null(lcut)){
    print(" 'lcut' not specified    ... defaulting to peak separation of peaks with a height difference of 50% (lcut=0.5) ... ")
    lcut <- 0.5
  }
  
  if(is.null(pAdjust)){
    print(" 'pAdjust' not specified ... not adjusting for multiple correction (recommended) ... ")
    pAdjust <- 'none'
  }
  
  cat('\n')
  
  if(Plot==T & length(proteins)!=1){
    continue <- NA
    while(continue%notin%c('y','Y','ye','YE','Ye','yes','YES','Yes','YeS','YEs','yeS','yES','n','N','no','No','NO')){
      continue <- readline('You are trying to plot data for more than 1 protein... Do you want to turn off plotting? \n (options: "y" or "n") ')
    }
    if(continue%in%c('n','N','no','No','NO')){
      Plot <- F
    }
    cat('\n')
  }
  
  if('Fits'%notin%names(array_object$Experiments[[dataset]])){
    cat(paste0('No fits calculated for the specified dataset "',dataset,'" ...\n'))
    stop()
  }
  
  if('Imputation_fit'%notin%names(array_object$Experiments[[dataset]])){
    cat(paste0('No fit for imputation has been determined for the specified dataset "',dataset,'" ...\n'))
    stop()
  }
  
  # Define input to get_score function -----------------------------------------
  selected_proteins <- find_associated_terminals(proteins,array_object)

  # Run get_score function -----------------------------------------------------
  list_of_scores <- lapply(selected_proteins,
                           .get_score, 
                           array_object=array_object,
                           pthresh=pthresh,
                           lcut=lcut,
                           cut=cut,
                           pAdjust=pAdjust,
                           Plot=Plot,
                           get_min=get_min,
                           dataset=dataset)
  
  table_of_scores <- rbindlist(list_of_scores)
  
  # Save parameters used -------------------------------------------------------
  fits   <- get_fits(array_object, dataset, verbose = F)
  impute <- get_imputation_information(array_object, dataset, verbose = F)
  
  call <- list(pthresh=pthresh, 
               lcut=lcut ,
               pAdjust=pAdjust,
               impute=impute$parameters,
               impute_distribution=impute$name,
               distribution=fits$name,
               Plot=Plot,
               get_min=get_min,
               dataset=dataset) 
  
  # Return results as list -----------------------------------------------------
  result <- list(SWaFi_scores=table_of_scores,Parameters=call) 
  return(result)
}

# this is the final function to get bidning sits
get_binding_sites <- function(proteins=NULL,
                              binding=NULL, 
                              datasets=c('FBBE1','FBBE2'),
                              confidence='low',
                              merizo.overlap = 'within', 
                              merizo.pIoU_cutoff = 0.75,
                              merizo.plDDT_cutoff = 70, 
                              keep = c('untested','no','partially'),
                              array_object=NULL){
  
  # Allows for different inputs as confidence numeric and character
  if(is.character(confidence)==T){
    if(str_detect(confidence,'l|low|L|Low')){
      confidence <- 'low_confidence'
    }
    
    if(str_detect(confidence,'m|medium|M|Medium')){
      confidence <- 'medium_confidence'
    }
    
    if(str_detect(confidence,'h|high|H|High|Hi|hi')){
      confidence <- 'high_confidence'
    }
  }
  
  if(is.numeric(confidence)==T){
    if(confidence==0.1){
      confidence <- 'low_confidence'
    }
    
    if(confidence==0.05){
      confidence <- 'medium_confidence'
    }
    
    if(confidence==0.01){
      confidence <- 'high_confidence'
    }
  }
  
  # get array_oject
  if(is.null(array_object)){
    array_object <- Find_array_object()
  }
  
  # get binding site experiment details
  item_name <- paste(datasets,collapse = '_')
  
  if(!item_name%in%names(array_object$Binding_sites)){
    cat('Sites have not been calculated with these parameters ... \nOr check your input ... \n')
    stop()
  }
  
  # get sites 
  sites <- array_object$Binding_sites[[item_name]][[confidence]]
  
  # filter for proteins queried
  if(!is.null(proteins)){
    sites <- sites[strip_terminal_specifier(sites$UniProtKB.AC)%in%strip_terminal_specifier(proteins),]
  }
  
  # keep binding sites only if requested
  if(!is.null(binding)){
    sites <- sites[sites$signif==binding,]
  }
  
  # get merizo domains
  merizo <- get_merizo_domains(overlap = merizo.overlap,
                               pIoU_cutoff = merizo.pIoU_cutoff,
                               plDDT_cutoff = merizo.plDDT_cutoff,
                               proteins = proteins,
                               array_object = array_object)
  
  # add merizo info to binding site df
  sites <- proportion_site_in_domain(binding_sites = sites, 
                                     merizo_domain_regions = merizo)
  
  sites <- sites[sites$in_domain%in%keep,]
  
  return(sites)
} 

# get_theoretical_maximum_SWaFi_score <- function(array_object=NULL, 
#                                                 dataset=NULL){
#   
#   # Find array_object if missing -----------------------------------------------
#   if(is.null(array_object)){
#     array_object <- Find_array_object()
#   } 
#   
#   if(is.null(dataset)){
#     print('Please specify the dataset ...')
#     stop()
#   }
#   
#   normalised_scores <- get_dataset(array_object,dataset)$normalised 
#   maximum_terminal_length <- get_longest_terminal(amino_acid_resolution = 2,array_object)
#   mx <- max(.get_signal_fast(rep(max(normalised_scores,na.rm=T),maximum_terminal_length)))
#   return(mx)
# }
# 
# get_longest_terminal <- function(amino_acid_resolution=1 ,
#                                  array_object=NULL){
#   
#   # amino acid resolution can be used to get the maximum number probes representing a terminal
#   
#   # Find array_object if missing -----------------------------------------------
#   if(is.null(array_object)){
#     array_object <- Find_array_object()
#   }
#   
#   floor(max(unlist(get_metadata(array_object,'seq.length')))/amino_acid_resolution)
# }
# 
# get_qq_plot_of_fit_residuals <- function(dataset, array_object=NULL){
#   
#   if(is.null(array_object)){
#     array_object <- Find_array_object()
#   }
#   
#   fits <- get_fits(array_object,dataset)$Fits
#   
#   qq <- pblapply(1:length(fits),function(x){
#     f <- fits[[x]]
#     q <- qqnorm(f$residuals,plot.it=F)
#     data.frame(Terminal_length=x,Empirical=q$y,Theoretical=q$x)
#   })
#   
#   return(qq)
# }
# 
# get_qq_plot <- function(dataset,
#                         terminal_length,
#                         model=NULL,
#                         fit=NULL,
#                         array_object=NULL){
#   
#   experiment <- get_experiment(array_object,dataset)
#   
#   if(is.null(fit)){
#     fit <- get_fits(array_object,dataset,verbose = F)$Fits[[terminal_length]]
#   }
#   
#   if(is.null(model)){
#     model <- fit$family[1]
#   }
#   
#   sampling <- get_random_samplings(array_object,dataset)[[terminal_length]]
#   
#   r <- random_generation_function(length(sampling),model,fit)
#   
#   qq <- qqplot(sampling,r, plot.it = F)
#   qq <- data.frame(Terminal_length=terminal_length,Empirical=qq$y,Theoretical=qq$x)
#   return(qq)
# }

get_masked_SWaFi_traces <- function(dataset,
                                    proteins=NULL,
                                    array_object=NULL,
                                    pthresh=NULL,
                                    get_min=F, ...){
  
  if(is.null(array_object)){
    array_object <- Find_array_object()
  }
  
  traces    <- get_SWaFi_traces(array_object = array_object,
                                dataset = dataset,
                                proteins = proteins,
                                parameters=F)
  
  all_masks <- create_SWaFi_trace_mask(proteins = proteins,
                                       array_object = array_object,
                                       pthresh = pthresh,
                                       get_min = get_min, ...)
  
  masked_traces <- lapply(names(traces),function(id){
    trace <- traces[[id]]
    trace_mask <- all_masks[[id]]
    trace_mask <- trace_mask[trace_mask$position%in%trace$position,]
    trace$score <- trace$score*trace_mask$mask
    return(trace)
  })
  
  names(masked_traces) <- names(traces)
  return(masked_traces)
}

get_SWaFi_matrix <- function(dataset,
                             proteins=NULL,
                             array_object=NULL,
                             pthresh=NULL,
                             get_min=F, 
                             TM_width=21, ...){
  
  if(is.null(array_object)){
    array_object <- Find_array_object()
  }
  
  traces <- get_SWaFi_traces(array_object,dataset,proteins,parameters=F)
  
  if(!is.null(pthresh)){
    cat(paste0('\nOnly showing significant binding sites (p<',pthresh,') ... \n'))
    traces <- get_masked_SWaFi_traces(dataset,proteins,array_object,pthresh,get_min, ...)
  }
  
  traces <- lapply(traces,function(x){x$score}) # keep scores only, positions will be relative to TM segment from now on

  # divide into N and C termini
  
  if(sum(str_detect(names(traces),'_N'))>0){
    N_term <- traces[str_detect(names(traces),'_N')]
    # Get longest terminal coverage
    N_len <- max(unlist(lapply(N_term,length)))
    # Extend terminals with NA to same length
    N_term <- lapply(N_term,function(x){c(rep(NA,N_len-length(x)),x)})
    # Remove terminal identifiers
    names(N_term) <- strip_terminal_specifier(names(N_term))
    # Convert to data.frames
    N_term <- as.data.frame(do.call(rbind,N_term))
    colnames(N_term) <- paste0('V',1:ncol(N_term))
    N_term <- rownames_to_column(N_term,var='ID')
    # make unique colnames for merging
    colnames(N_term) <- c('ID',paste0('N-',1:(ncol(N_term)-1)))
  } else {
    N_term <- data.frame(ID=strip_terminal_specifier(names(traces)), Empty=NA)
  }
  
  if(sum(str_detect(names(traces),'_C'))>0){
    C_term <- traces[str_detect(names(traces),'_C')]
    # Get longest terminal coverage
    C_len <- max(unlist(lapply(C_term,length)))
    # Extend terminals with NA to same length
    C_term <- lapply(C_term,function(x){c(x,rep(NA,C_len-length(x)))})
    # Remove terminal identifiers
    names(C_term) <- strip_terminal_specifier(names(C_term))
    # Convert to data.frames
    C_term <- as.data.frame(do.call(rbind,C_term))
    colnames(C_term) <- paste0('V',1:ncol(C_term))
    C_term <- rownames_to_column(C_term,var='ID')
    # make unique colnames for merging
    colnames(C_term) <- c('ID',paste0('C-',1:(ncol(C_term)-1)))
  } else {
    C_term <- data.frame(ID=strip_terminal_specifier(names(traces)), Empty=NA)
  }
  
  # merge data.frames
  transmembrane_segment <- as.data.frame(matrix(NA,ncol=TM_width,nrow=nrow(N_term)))
  colnames(transmembrane_segment) <- paste0('TM-',1:(ncol(transmembrane_segment)))
  N_term <- cbind(N_term,transmembrane_segment)
  
  # Merge terminals
  merged <- merge(x=N_term,y=C_term,by='ID',all=T,.name_repair = "minimal")
  merged <- column_to_rownames(merged,var='ID')
  empty  <- which(str_detect(colnames(merged),'Empty')==T)
  if(length(empty)>0){
    merged <- merged[,-empty]
  }
  merged <- as.matrix(merged)
  return(merged)
}

# get_panther_overrepresentation <- function(proteins,
#                                            database='all',
#                                            Background=NULL,
#                                            cutoff=0.05,
#                                            pAdjust='fdr',
#                                            test_type='fisher',
#                                            organism=9606){
#   
#   opts <- c('all','GO','GOslim','PantherPathway','ReactomePathway','PantherProteinClass')
#   
#   if(database%notin%opts){
#     cat("\nERROR: 'database' has to be one of: \n")
#     cat(paste(c('all','GO','GOslim','PantherPathway','ReactomePathway','PantherProteinClass'),collapse = '|'))
#     cat('\n')
#     stop()
#   }
#   
#   proteins <- unique(proteins)
# 
#   id_type <- detect_identifier(proteins)
#   
#   if(id_type!="UniProtKB.AC"){
#     proteins <- unique(strip_terminal_specifier(proteins))
#   }
#   
#   if(!is.null(Background)){
#     Background <- unique(Background)
#     id_type <- detect_identifier(Background)
#     
#     if(id_type!="UniProtKB.AC"){
#       Background <- unique(strip_terminal_specifier(Background))
#     }
#     
#     ref_organism <- organism
#   } else {
#     ref_organism <- NULL
#   }
#   
#   if(is.null(database)){
#     cat('\nPlease specify a database, the options in PANTHER are: \n ')
#     print(c('all',rba_panther_info(what='datasets')[,4]))
#     stop()
#   } 
#   
#   if(database%in%opts){
#     selected <- database
#     database <- rba_panther_info(what='datasets')[,3:4]
#     name <- database$label
#     database <- database$id
#     
#     if(selected=='GO'){
#       name <- name[1:3]
#       database <- database[1:3]
#     }
#     
#     if(selected=='GOslim'){
#       name <- name[4:6]
#       database <- database[4:6]
#     }
#     
#     if(selected=='PantherPathway'){
#       name <- name[8]
#       database <- database[8]
#     }
#     
#     if(selected=='ReactomePathway'){
#       name <- name[9]
#       database <- database[9]
#     }
#     
#     if(selected=='PantherProteinClass'){
#       name <- name[7]
#       database <- database[7]
#     }
#     
#   } else {
#     possibilites <- rba_panther_info(what='datasets')[,3:4]
#     if(is.character(database)){
#       database <- possibilites$id[possibilites$label==database]
#       name <- possibilites$label[possibilites$label==database]
#     } else {
#       database <- possibilites$id[database-1]
#       name <-  possibilites$label[database-1]
#     }
#     
#   }
#   
#   pAdjust <- toupper(pAdjust)
#   test_type <- toupper(test_type)
#   
#   results <- lapply(database,function(x){
#     enriched <- rba_panther_enrich(genes = proteins, 
#                                    annot_dataset = x, 
#                                    organism = organism,
#                                    cutoff = cutoff, 
#                                    test_type = test_type,
#                                    correction = pAdjust,
#                                    ref_genes = Background, 
#                                    ref_organism = ref_organism)
#     
#     enriched <- enriched$result
#     enriched <- enriched[order(enriched$pValue,decreasing = T),]
#     return(enriched)
#   })
#   
#   names(results) <- name
#   
#   results <- do.call(rbind,results)
#   
#   results <- rownames_to_column(results,var = 'database')
#   results$database <- gsub('\\..*','',results$database)
#   results$database <- gsub('_',' ',results$database)
#   
#   colnames(results) <- c('database','#','fold_change',tolower(pAdjust),'expected','#_ref','pval','+/-','term_id','term_label')
#   
#   print(as_tibble(results), )
#   
#   return(results)
#   
# }
# 
# get_panther_annotations <- function(UniProtKB.AC,
#                                     database='all',
#                                     organism=9606){
#   
#   opts <- c('all','GO','GOslim','PantherPathway','ReactomePathway','PantherProteinClass')
#   
#   if(database%notin%opts){
#     cat("\nERROR: 'database' has to be one of: \n")
#     cat(paste(c('all','GO','GOslim','PantherPathway','ReactomePathway','PantherProteinClass'),collapse = '|'))
#     cat('\n')
#     stop()
#   }
#   
#   if(length(UniProtKB.AC)>1000){
#     intervals <- seq(1,length(UniProtKB.AC),1000)
#   } else {
#     intervals <- 1
#   }
#   
#   annotations <- lapply(intervals, function(index){
#     
#     to <- ifelse((index+999)>length(UniProtKB.AC),length(UniProtKB.AC),(index+999))
#     
#     annotations <- rba_panther_mapping(UniProtKB.AC[index:to], organism = organism)
#     
#     unmapped <- unlist(annotations$unmapped_list$unmapped)
#     
#     mapped <- annotations$mapped_genes$gene
#     
#     if('accession' %in% names(mapped)){
#       mapped <- list(mapped)
#     }
#     
#     names(mapped) <- unlist(lapply(mapped,function(x){x$accession}))
#     names(mapped) <- gsub('.*UniProtKB=','',names(mapped))
#     
#     df <- lapply(mapped,function(protein){
#       UniProtid <-  gsub('.*UniProtKB=','',protein$accession)
#       protein <- protein$annotation_type_list$annotation_data_type
#       
#       if('annotation_list' %in% names(protein)){
#         protein <- list(protein)
#       }
#       
#       results <- do.call(rbind,lapply(protein,function(x){
#         db      <- x$content
#         version <- x$release_version
#         if(is.list(x$annotation_list$annotation[[1]])){
#           name    <- unlist(lapply(x$annotation_list$annotation,function(y){y$name}))
#           id      <- unlist(lapply(x$annotation_list$annotation,function(y){y$id}))
#         } else {
#           name    <- x$annotation_list$annotation$name
#           id      <- x$annotation_list$annotation$id
#         }
#         
#         terms   <- as.data.frame(cbind(id,name))
#         terms$database <- db
#         terms$version  <- as.character(version)
#         return(terms)
#       }))
#       
#       results$UniProtKB.AC <- UniProtid
#       
#       return(results)
#       
#     })
#     
#     df <- df[names(unlist(lapply(df,ncol)))]
#     
#     df <- do.call(rbind,df)
#     
#     rownames(df) <- NULL
#     
#     if(database=='GO'){
#       terms <- df[df$database%in%c("GO:0003674","GO:0008150","GO:0005575"),]
#       unmapped <- c(unmapped,df$UniProtKB.AC[df$UniProtKB.AC%notin%terms$UniProtKB.AC])
#       df <- terms
#     }
#     
#     if(database=='GOslim'){
#       terms <- df[df$database%in%c("ANNOT_TYPE_ID_PANTHER_GO_SLIM_CC","ANNOT_TYPE_ID_PANTHER_GO_SLIM_BP","ANNOT_TYPE_ID_PANTHER_GO_SLIM_MF"),]
#       unmapped <- c(unmapped,df$UniProtKB.AC[df$UniProtKB.AC%notin%terms$UniProtKB.AC])
#       df <- terms
#     }
#     
#     if(database=='PantherPathway'){
#       terms <- df[df$database%in%'ANNOT_TYPE_ID_PANTHER_PATHWAY',]
#       unmapped <- c(unmapped,df$UniProtKB.AC[df$UniProtKB.AC%notin%terms$UniProtKB.AC])
#       df <- terms
#     }
#     
#     if(database=='ReactomePathway'){
#       terms <- df[df$database%in%'ANNOT_TYPE_ID_REACTOME_PATHWAY',]
#       unmapped <- c(unmapped,df$UniProtKB.AC[df$UniProtKB.AC%notin%terms$UniProtKB.AC])
#       df <- terms
#     }
#     
#     if(database=='PantherProteinClass'){
#       terms <- df[df$database%in%'ANNOT_TYPE_ID_PANTHER_PC',]
#       unmapped <- c(unmapped,df$UniProtKB.AC[df$UniProtKB.AC%notin%terms$UniProtKB.AC])
#       df <- terms
#     }
#     return(list(unmapped=unmapped,mapped=df))
#   })
#   
#   mapped <- do.call(rbind,lapply(annotations,function(x){x$mapped}))
#   rownames(mapped) <- NULL
#   unmapped <- unlist(lapply(annotations,function(x){x$unmapped})) 
#   annotations <- list(unmapped=unmapped,mapped=mapped)
#   return(annotations)
# }
# 
# map_panther_annotations <- function(proteins,
#                                     panther_term_name=NULL,
#                                     panther_term_id=NULL,
#                                     database='all',
#                                     organism=9606,
#                                     array_object=NULL){
#   
#   if(is.null(array_object)){
#     array_object <- Find_array_object()
#   }
#   
#   annotations <- get_panther_annotations(UniProtKB.AC = proteins,database=database,organism = organism)
#   if(!is.null(panther_term_name)){
#     annotations$mapped <- annotations$mapped[annotations$mapped$name==panther_term_name,]
#   }
#   if(!is.null(panther_term_id)){
#     annotations$mapped <- annotations$mapped[annotations$mapped$id==panther_term_id,]
#   }
#   return(names(get_terminals(array_object = array_object, proteins = annotations$mapped$UniProtKB.AC)))
# }
# 
# get_GO_terms <- function(UniProtKB.AC){
#   
#   GOterms <- get_panther_annotations(UniProtKB.AC,'GO')
#   
#   mapped <- GOterms$mapped
#   
#   labels <- data.frame(database=c("GO:0003674","GO:0008150","GO:0005575"),db_label=c('GOMF','GOBP','GOCC'))
#   
#   mapped <- merge(mapped,labels,by='database',all.x=T)
#   mapped$fgsea_name <- gsub('-',' ',mapped$name)
#   mapped$fgsea_name <- gsub("\\(|\\)|\\[|\\]|'","",mapped$fgsea_name)
#   mapped$fgsea_name <- toupper(gsub(' |/','_',mapped$fgsea_name))
#   mapped$fgsea_name <- paste0(mapped$db_label,'_',mapped$fgsea_name)
#   
#   terms <- split(mapped$UniProtKB.AC,mapped$fgsea_name)
#   terms <- list(Associations=terms, metadata=mapped)
#   return(terms)
# }
# 
# get_proteins_by_keyword <- function(keywords, array_object=NULL){
#   
#   if(length(keywords)>1){
#     cat(paste0("\n                                      >>>>> ATTENTION <<<<<",
#                "\nIf you supply a vector ALL criteria in the vector have to be",
#                " fullfilled for a protein to be returned ... \nIf you want it to return",
#                " a protein if it matches ANY of the criteria,\n   supply the keywords",
#                " with the 'or' operator ('|') like this:\n\n 'olfaction|g-protein coupled receptor' \n\n",
#                "Always double check the output using the function: 'get_name_from_uniprot'\n"))
#   }
#   
#   keys  <- get_metadata(array_object, colnames='Keywords')
#   words <- tolower(paste(paste0('(?=.*', keywords, ')'), collapse = ''))
#   
#   res <- unlist(pblapply(keys,function(x){
#     str_detect(tolower(x),words)
#   }))
# 
#   cat('\n')
#   
#   return(names(res[res==T]))
# }

get_binding_proteins <- function(dataset,
                                 proteins=NULL,
                                 pthresh=0.05,
                                 useAdjusted=F,
                                 get_min=F,
                                 array_object=NULL){
  all_sites <- get_SWaFi_binding_sites(array_object=NULL,dataset=dataset,proteins=proteins,pthresh=pthresh,useAdjusted=useAdjusted,get_min=get_min)
  all_sites <- do.call(rbind,all_sites)
  binding_sites <- all_sites[all_sites$signif==T,]
  binding_proteins <- strip_terminal_specifier(binding_sites$UniProtKB.AC)
  return(unique(na.omit(binding_proteins)))
}

get_non_binding_proteins <- function(dataset, 
                                     proteins=NULL, 
                                     pthresh=0.05,
                                     useAdjusted=F,
                                     get_min=F,
                                     array_object=NULL){
  all_sites <- get_SWaFi_binding_sites(array_object=NULL,dataset=dataset,proteins=proteins,pthresh=pthresh,useAdjusted=useAdjusted,get_min=get_min)
  all_sites <- do.call(rbind,all_sites)
  binding_sites <- all_sites[all_sites$signif==T,]
  non_binding_sites <- all_sites[strip_terminal_specifier(all_sites$UniProtKB.AC)%notin%strip_terminal_specifier(binding_sites$UniProtKB.AC),]
  non_binding_sites <- strip_terminal_specifier(non_binding_sites$UniProtKB.AC)
  return(unique(na.omit(non_binding_sites)))
}

get_sequence_alignment <- function(seq_colname,
                                   proteins,
                                   method='Muscle',
                                   substitution_matrix=NULL,
                                   array_object=NULL, 
                                   gap_opening=NULL,
                                   gap_extension=NULL){
  
  seqs <- unlist(get_metadata(array_object, colnames=seq_colname,proteins = proteins))
  names(seqs) <- strip_terminal_specifier(names(seqs))
  seqs <- seqs[unique(names(seqs))]
  stringset <- AAStringSet(seqs)
  
  if(is.null(substitution_matrix)){
    data("BLOSUM62")
    substitution_matrix <- BLOSUM62
    cat('\nUsing default substitution matrix - BLOSUM62 ... \n')
    if(method%in%c('ClustalOmega','ClustalW')){
      substitution_matrix <- "BLOSUM65"
    }
  }
  
  # Run Multiple Sequence Alignment ----------------------------------------------
  if(is.null(gap_opening)){
    cat("\nGap penalties not specified using method default: \n")
    alignment <- msa(stringset, method = method, substitutionMatrix = substitution_matrix, verbose=F)
  }
  
  if(!is.null(gap_opening)){
    cat('\nUsing provided parameters: \n')
    alignment <- msa(stringset, method = method, substitutionMatrix = substitution_matrix, 
                     gapOpening = gap_opening, gapExtension = gap_extension, verbose = F)
  }
  
  seqinr_alignment <- msaConvert(alignment, type="seqinr::alignment")
  d <- seqinr::dist.alignment(seqinr_alignment) 
  d[is.na(d)] <- 1
  
  return(list(dist_matrix=d,alignment=unmasked(alignment)))
}

get_AA_frequency <- function(sequences=NULL,
                             proteins=NULL,
                             freq=T,
                             remove_AA=c('X','U')){
  if(!is.null(sequences)){
    AAs <- unlist(str_split(sequences,pattern=''))
    AAs <- AAs[AAs %notin% remove_AA]
  } else {
    if(!is.null(proteins)){
      AAs <- unlist(lapply(get_metadata(colnames='sub.sequence',proteins=proteins),str_split,pattern=''))
      AAs <- AAs[AAs %notin% remove_AA]
    } else {
      print('Error ... please provide input to the function ...')
      stop()
    }
  }
  
  if(freq==F){
    freq <- table(AAs)/length(AAs)
  } else {
    freq <- table(AAs)
  }
  return(freq)
}

get_number_of_TM_segments <- function(proteins=NULL,
                                      no.of.segments=NULL,
                                      array_object=NULL){
  
  if(is.null(array_object)){
    array_object <- Find_array_object()
  }
  
  TM_segments <- do.call(rbind,lapply(get_metadata(array_object = array_object, colnames = 'Transmembrane_segments', proteins = proteins),nrow))
  TM_segments <- rownames_to_column(as.data.frame(TM_segments), var = 'UniProtKB.AC')
  TM_segments$UniProtKB.AC <- strip_terminal_specifier(TM_segments$UniProtKB.AC)
  TM_segments <- unique(TM_segments)
  colnames(TM_segments)[2] <- 'TM_Segments'
  
  if(!is.null(no.of.segments)){
    TM_segments <- TM_segments[TM_segments$TM_Segments%in%no.of.segments,]
  }
  
  return(TM_segments)
}

# get_subcellular_localization <- function(keyword,
#                                          proteins=species_filter(),
#                                          confidence=3, 
#                                          array_object=NULL){
#   
#   if(is.null(array_object)){
#     array_object <- Find_array_object()
#   }
#   
#   compartments <- get_metadata(array_object = array_object, colnames = 'COMPARTMENTS', proteins = proteins)
#   
#   compartments <- rownames_to_column(na.omit(do.call(rbind,lapply(compartments, function(x){
#     x <- x[x$Confidence>=confidence,]
#     x[grepl(tolower(keyword), tolower(x$GO_term), fixed = T),]
#   }))), var='UniProtKB.AC')
#   
#   if(nrow(compartments)<1){
#     cat('No good matches found, consider using another keyword \n')
#     stop()
#   }
#   
#   cat('GO terms matched to keyword:\n')
#   print(sort(table(compartments$GO_term)))
#   
#   return(compartments)
# }

get_short_protein_name <- function(proteins, array_object=NULL){
  
  if(is.null(array_object)){
    array_object <- Find_array_object() 
  }
  
  name <- unlist(get_metadata(array_object = array_object,
                              colnames = "Protein.name",
                              proteins = proteins))
  
  name <- gsub(' \\(.*','',name)
  
  return(name)
}

# calculate_protein_sequence_dissimilarity <- function(proteins, 
#                                                      alignment_type = 'local',
#                                                      colname = 'Full.Sequence',
#                                                      array_object=NULL,
#                                                      sequences=NULL){
#   
#   cat('array_object wrapper for protr::parSeqSim() function\nRunning alignment:\n"global" - Needleman-Wunsch global alignment\n"local"  - Smith-Waterman local alignment\nYou can run this for any metadata sequences available with function get_colnames()\nYou can also add a named vector of sequences manually using argument sequences\n')
#   
#   if(is.null(array_object)){
#     array_object <- Find_array_object()
#   }
#   
#   if(!is.null(sequences)){
#     cat('Calculating similarity for provided sequences ... \n')
#     seqs <- sequences[unique(names(sequences))]
#     
#   } else {
#     if(colname%notin%get_colnames(array_object=array_object)){
#       cat('Provided colname not in array_object ... choose one of the following options: \n')
#       print(get_colnames(array_object=array_object))
#       stop()
#     }
#     cat(paste0('Getting sequences from array_object, with colname: ', colname,'\n'))
#     seqs <- get_metadata(colnames = 'Full.Sequence', proteins=proteins)
#     names(seqs) <- strip_terminal_specifier(names(seqs))
#     seqs <- unlist(seqs)[unique(names(unlist(seqs)))]
#   }
#   
#   seq_sim <- protr::parSeqSim(seqs, type = alignment_type)
#   rownames(seq_sim) <- names(seqs)
#   colnames(seq_sim) <- names(seqs)
#   
#   seq_dissim <- as.dist(1-seq_sim)
#   return(seq_dissim)
# }

pdb_structures <- function(SWaFi_score){
  apply(SWaFi_score,1,function(x){
    name <- x[1] #paste0(x[1],'_C')
    #name <- ifelse(name %in% names(Array_object),name,paste0(x[1],'_N'))
    y <- Array_object$Protein_terminals[[name]]$AlphaFold
    y <- y[y$residue %in% c(x[3]:x[4]),]
    y <- y[,c('aa','structure','phi','psi','confidence')]
    y$aa <- a(str_to_title(tolower(y$aa)))
    y
  })
}

map_uniprot_to_genes <- function(uniprot, get_first=T){
  
  genes <- unlist(lapply(uniprot,function(x){unique(unlist(get_metadata(proteins = x,colnames='Genes')))}))
  
  if(get_first==T){
    genes <- gsub(' .*','',genes)
  }
  return(genes)
}

convert_df_to_SWaFi_score_object <- function(df, 
                                             UniProtKB.AC_column='UniProtKB.AC', 
                                             site_start_column='start', 
                                             site_end_column='end',
                                             array_object=NULL){
  
  if(is.null(array_object)){
    array_object <- Find_array_object()
  }
  
  if(sum(c(UniProtKB.AC_column,site_start_column,site_end_column)%notin%colnames(df))>0){
    cat('The provided column names are not in the df ... \n')
    stop()
  }
  
  sites <- df
  colnames(sites)[which(colnames(sites)==UniProtKB.AC_column)] <- 'UniProtKB.AC'
  colnames(sites)[which(colnames(sites)==site_start_column)] <- 'start'
  colnames(sites)[which(colnames(sites)==site_end_column)] <- 'end'
  
  proteins <- sites$UniProtKB.AC
  
  proteins <- find_associated_terminals(proteins,array_object)
  
  sites$length <- sites$start-sites$end+1
  
  TM_start_end  <- rownames_to_column(do.call(rbind,get_metadata(array_object,colnames = c('sub.start','sub.end','terminal'),proteins = proteins)))  
  colnames(TM_start_end)[1] <- 'UniProtKB.AC'
  TM_start_end$sub.start <- as.numeric(TM_start_end$sub.start)
  TM_start_end$sub.end   <- as.numeric(TM_start_end$sub.end)
  names(TM_start_end)[2:3] <- c('terminal_start','terminal_end')
  
  sequences <- get_metadata(array_object,colnames = 'Full.Sequence',proteins=proteins)
  sequences <- rownames_to_column(as.data.frame(unlist(sequences)))
  names(sequences) <- c('UniProtKB.AC','full_protein_sequence')  
  sites <- merge(sites,sequences,by='UniProtKB.AC',all.x=T)
  sites$sequence <- substr(sites$full_protein_sequence, start=sites$start, stop=sites$end)
  
  sites <- merge(sites,TM_start_end,by='UniProtKB.AC',all.x=T)
  sites$terminal_length <- sites$terminal_end-sites$terminal_start+1 
  
  sites$distance_from_TM <- ifelse(sites$terminal=='C',sites$start-sites$terminal_start+1,sites$terminal_end-sites$end+1)
  sites$relative_distance_from_TM <- sites$distance_from_TM/sites$terminal_length
  
  sites <- na.omit(sites)
  
  return(sites)
}

get_offset_between_datasets <- function(dataset1='FBBE1', 
                                        dataset2='FBBE2', 
                                        array_object=NULL){
  
  if(is.null(array_object)){
    array_object <- Find_array_object()
  }
  
  #get average offset between datasets
  slope <- unlist(lapply(names(get_terminals(array_object = array_object)),function(x){
    
    SWaFi_trace_dataset1 <- get_SWaFi_traces(dataset = dataset1, 
                                             proteins = x,
                                             array_object = array_object)[[1]]$score
    
    SWaFi_trace_dataset2 <- get_SWaFi_traces(dataset = dataset2, 
                                             proteins = x,
                                             array_object = array_object)[[1]]$score
    
    slope <- lm(SWaFi_trace_dataset2 ~ SWaFi_trace_dataset1)$coefficients[2]
    return(slope)
    
  }))
  
  quiet(print({
    hist(slope, breaks = 200, xlim=c(-5,10))
    abline(v=median(slope,na.rm=T), col='red', lty=1, lwd=2)
  }))
  
  k <- median(slope,na.rm=T)
  
  return(k)
}

get_binding_site_conservation <- function(binding_sites = NULL,
                                          pthresh = 0.05,
                                          test = t.test,
                                          impute = F,
                                          impute_mean = NULL,
                                          impute_sd = NULL,
                                          impute_n = 10,
                                          one_sample_t.test = F,
                                          array_object = NULL, ...){
  
  # Find array object if missing
  if(is.null(array_object)){
    array_object <- Find_array_object()
  }
  
  # Get binding sites if missing
  if(is.null(binding_sites)){
    binding_sites <- get_binding_sites(confidence = pthresh, 
                                       binding = T, 
                                       array_object = array_object)
  }
  
  if(max(binding_sites$fisher_pval, na.rm=T)>pthresh){
    message(paste0('binding_sites contains one or more sites above the significance threshold of ',pthresh))
    message('exiting ... ')
    stop()
  }
  
  if(one_sample_t.test==T){
    message('Running One-sample t.test: tests site conservation vs protein average (mu=1)')
  }
  
  # Get list of all binding sites in case a subset is provided
  if(one_sample_t.test==F){
    all_sites <- get_binding_sites(confidence = pthresh, 
                                   proteins = binding_sites$UniProtKB.AC,
                                   array_object = array_object)
  }
  
  # Get all conservation data for terminals
  terminal.conservation <- get_metadata(array_object = array_object, 
                                        colnames = 'terminal.conservation')
  
  # If missing estimate good imputation parameters
  if(impute==T){
    
    if(is.null(impute_mean) & one_sample_t.test==F){
      message('No imputation parameters provided: estimating from sites with p-value > 0.95 ... ')
      
      non_binding_sites <- get_binding_sites(confidence = pthresh,
                                             binding = F,
                                             array_object = array_object)
      
      non_binding_sites <- non_binding_sites[non_binding_sites$fisher_pval>0.95,]
      
      non_binding_conservation <- na.omit(unlist(apply(non_binding_sites,1,function(non_binding_site){
        
        cons  <- terminal.conservation[[non_binding_site['UniProtKB.AC']]]
        cons$WRSRscore[cons$residue %in% c(non_binding_site['start']:non_binding_site['end'])]
        
      })))
      
      impute_sd   <- sd(non_binding_conservation,na.rm=T)/2
      impute_mean <- mean(non_binding_conservation,na.rm=T)
      hist(non_binding_conservation,
           breaks = 100,
           main='Distribution used to determine imputation parameters', 
           xlab='Substitution rate')
      abline(v=impute_mean, col='red')
      abline(v=c(impute_mean-impute_sd,impute_mean+impute_sd), col='red', lty=2)
      
      message(paste0('impute_mean set to: ',impute_mean))
      message(paste0('impute_sd   set to: ',impute_sd))
    }
  }
  
  
  # Loop through all sites and calculate fold change and significance
  site_conservation <- apply(binding_sites,1,function(main_site){
    
    # Get site info
    id    <- main_site['UniProtKB.AC'] 
    start <- main_site['start']
    end   <- main_site['end']

    # Get site conservation
    conservation <- do.call(rbind,terminal.conservation[paste0(strip_terminal_specifier(id),c('_C','_N'))])
    
    # Test if conservation data exists and has the correct content
    if(is.null(conservation)){
      conservation <- NA
    } else {
      if(sum(conservation$residue %in% start:end) != length(start:end)){
        conservation <- NA
      }
    }
    
    # If conservation data exists continue
    if(is.data.frame(conservation)==T){
      
      # Extract conservation of main_site
      main_site_conservation <- conservation$WRSRscore[conservation$residue %in% start:end]
      main_site_conservation <- na.omit(main_site_conservation)
      
      if(one_sample_t.test==F){
        # Identify other sites in terminal
        
        # other_binding_sites   <- all_sites[all_sites$signif==T,]
        # all_sites_in_terminal <- other_binding_sites[strip_terminal_specifier(other_binding_sites$UniProtKB.AC)==strip_terminal_specifier(id),]
        # 
        # all_sites_in_terminal <- unlist(apply(all_sites_in_terminal,1,function(site){
        #   site['start']:site['end']
        # }))
        
        # Extract non binding sites
        non_sites   <- all_sites[all_sites$signif==F,]
        all_non_sites_in_terminal <- non_sites[strip_terminal_specifier(non_sites$UniProtKB.AC)==strip_terminal_specifier(id),]
        
        if(nrow(all_non_sites_in_terminal)>0){
          
          all_non_sites_in_terminal <- unlist(apply(all_non_sites_in_terminal,1,function(site){
            site['start']:site['end']
          }))
          
        } else {
          all_non_sites_in_terminal <- NULL
        }
        
        non_site_conservation <- conservation$WRSRscore[conservation$residue %in% all_non_sites_in_terminal]
        
        
        # Impute data to have at least 10 points
        if(length(non_site_conservation)<10 & impute==T){
          
          non_site_conservation <- c(non_site_conservation,
                                     rnorm(n = impute_n-length(non_site_conservation), 
                                           mean = impute_mean, 
                                           sd = impute_sd)) 
          
          non_site_conservation[non_site_conservation<0] <- 0
          
        }
        
        non_site_conservation <- na.omit(non_site_conservation)
        
      } else {
        non_site_conservation <- c(1,1,1)
      }

      if(length(non_site_conservation)>2){
      
        # Add a little noise to data to prevent 0 variance
        main_site_conservation <- main_site_conservation + rnorm(length(main_site_conservation),
                                                                 mean = 0.0000001, 
                                                                 sd = 0.0000001) #add some noise
        
        if(one_sample_t.test==F){
          non_site_conservation <- non_site_conservation + rnorm(length(non_site_conservation),
                                                                 mean = 0.0000001,
                                                                 sd = 0.0000001) #add some noise
        } else {
          non_site_conservation <- 1
        }
        
        # higher score = higher conservation so positive log2 FC means higher conservation of binding vs non-binding
        main_site_conservation_mean <- mean(main_site_conservation)
        non_site_conservation_mean  <- mean(non_site_conservation)
        logFC <- log2(non_site_conservation_mean)-log2(main_site_conservation_mean)
        
        if(one_sample_t.test==F){
          # Welch t.test 
          test.res <- test(main_site_conservation, non_site_conservation, ...)
        }
        
        if(one_sample_t.test==T){
          test.res <- t.test(main_site_conservation, mu=1)
        }
        
        p_val <- test.res$p.value
        statistic <- test.res$statistic
        
      } else {
        main_site_conservation_mean <- NA
        non_site_conservation_mean <- NA
        logFC <- NA
        p_val <- NA
        statistic <- NA
      }
      
    } else {
      main_site_conservation_mean <- NA
      non_site_conservation_mean <- NA
      logFC <- NA
      p_val <- NA
      statistic <- NA
    }
    
    return(data.frame(WRSR_mean = main_site_conservation_mean, 
                      WRSR_mean_non_site = non_site_conservation_mean, 
                      WRSR_logfc = logFC, 
                      WRSR_pval  = p_val,
                      WRSR_statistic = statistic))
  })
  
  # Add conservation data to binding_site dataframe
  combined <- cbind(binding_sites,do.call(rbind,site_conservation))
  
  # return result
  return(combined)
  
}

# get_proteins_located_in <- function(term, 
#                                     confidence = 3, 
#                                     proteins=NULL, 
#                                     array_object=NULL){
#   
#   # Find array_object in environment
#   if(is.null(array_object)){
#     array_object <- Find_array_object()
#   }
#   
#   cat(paste0('\nFetching proteins associated with subcellular localisation ',term,' [Confidence - ',confidence,']',
#              '\n---------------------------------------------------------------------------------\n'))
#   
#   compartments <- get_metadata(colnames = 'COMPARTMENTS',proteins = proteins, array_object = array_object)
#   
#   associated <- do.call(rbind,lapply(compartments, function(x){
#     
#     associated_terms <- x[x$Confidence>=confidence,]
#     
#     associated_terms <- associated_terms[str_detect(tolower(associated_terms$GO_term),tolower(term)),]
#     
#     return(associated_terms)
#   }))
#   
#   associated <- rownames_to_column(associated,var = 'UniProtKB.AC')
#   associated$Array.ID <- gsub('\\..*','',associated$UniProtKB.AC)
#   associated$UniProtKB.AC <- strip_terminal_specifier(associated$Array.ID)
#   
#   message('\nFound: ',length(unique(associated$UniProtKB.AC)),' associated proteins!\n')
#   
#   return(as_tibble(associated))
# }


#-############################################################################-#
#-########## UPDATE FUNCTIONS ('update' metadata from array_object) ##########-#
#-############################################################################-#

# retrieve_uniprot_info <- function(proteins){
#   
#   connected <- rbioapi::rba_connection_test(print_output = F)$UniProt
#   if(connected==T){
#   info <- pblapply(proteins,function(protein){
#     protein <- strip_terminal_specifier(protein)
#     
#     status <- is.error(rbioapi::rba_uniprot_proteins(protein,verbose=F))
#     
#     if(status==F){
#       res <- rbioapi::rba_uniprot_proteins(protein,verbose=F)
#       res$status <- 'up-to-date'
#     } 
#     
#     if(status==T){
#       res <- list(status='outdated')
#     }
#     return(res)
#   })
#   names(info) <- proteins
#   
#   outdated <- which(unlist(lapply(info,length))==1)
#   outdated <- names(info)[outdated]
#   
#   if(length(outdated)>0){
#     outdated <- riffle(outdated, rep("\n",length(outdated)))
#     spacer <- "----------------------------------------------------------------"
#     cat(paste0("\nThe following proteins could not be retrieved from UniProt ... \n",
#                "The entries might be 'outdated' ... \n",spacer,"\n ", paste(outdated,collapse = ' '),"\n",spacer))
#   }
#   } else {
#     cat('\nERROR: UniProt server not responding ...\nQuick fix: please check your internet connection.\n')
#     info <- 'ERROR'
#   }
#   return(info)
# }
# 
# update_metadata <- function(array_object=NULL){
#   
#   if(is.null(array_object)){
#     array_object <- Find_array_object()
#   }
#   
#   terminals <- get_terminals(array_object)
#   info <- retrieve_uniprot_info(names(terminals))
#   
#   updated_terminals <- pblapply(names(terminals),function(name){
#     
#     old <- terminals[[name]]
#     new <- info[[name]]
#     
#     if(new$status=='up-to-date'){
#       
#       old$Keywords <- paste(new$keywords$value, collapse=';')
#       old$Protein.name <- new$protein$recommendedName$fullName$value
#       
#       References_exists <- 'References'%in%names(old)
#       old$References <- new$references$citation[c('type','publicationDate','title','authors','publication','dbReferences')]
#       
#       if(References_exists==F){
#         old <- old[c(1:which(names(old)=='Reviewed'),
#                        which(names(old)=='References'),
#                       (which(names(old)=='Reviewed')+1):(length(old)-1))]
#         }
#       
#       old$Genes <- new$gene$name$value
#       
#       if('FUNCTION'%in%new$comments$type){
#         old$Function <-  new$comments$text[new$comments$type=='FUNCTION'][[1]]$value
#         }
#       
#       if('SUBCELLULAR_LOCATION'%in%new$comments$type){
#         Topology_exists <- 'Topology'%in%names(old)
#         old$Subcellular.location <-  new$comments$locations[new$comments$type=='SUBCELLULAR_LOCATION'][[1]]$location$value
#         old$Topology <- new$comments$locations[new$comments$type=='SUBCELLULAR_LOCATION'][[1]]$topology$value
#         
#         if(Topology_exists==F){
#           old <- old[c(1:which(names(old)=='Subcellular.location'),
#                          which(names(old)=='Topology'),
#                         (which(names(old)=='Subcellular.location')+1):(length(old)-1))]
#           }
#         }
#       
#       } else {
#       
#       References_exists <- 'References'%in%names(old)
#       old$References <- NA
#       if(References_exists==F){
#         old <- old[c(1:which(names(old)=='Reviewed'),
#                        which(names(old)=='References'),
#                       (which(names(old)=='Reviewed')+1):(length(old)-1))]
#         }
#       
#       Topology_exists <- 'Topology'%in%names(old)
#       old$Topology <- NA
#       if(Topology_exists==F){
#         old <- old[c(1:which(names(old)=='Subcellular.location'),
#                        which(names(old)=='Topology'),
#                       (which(names(old)=='Subcellular.location')+1):(length(old)-1))]
#         }
#       
#       old$status <- new$status
#       }
#     
#     return(old)
#     
#   })
#   
#   names(updated_terminals) <- names(terminals)
#   array_object$Protein_terminals <- updated_terminals
#   updated_array_object <- update_date(array_object)
#   return(updated_array_object)
# }

#-############################################################################-#
#-############# PLOT FUNCTIONS ('plot' data from array_object) ###############-#
#-############################################################################-#

# plot_ramachandran <- function(SWaFi_score,
#                               plot=T,
#                               normalize_to_max=T,
#                               color='#6082B6'){
#   
#   b <- do.call(rbind,pdb_structures(transferase_scores))
#   b  <- na.omit(b) #remove NAs (occurs if the terminal contains either the first or last residue of a protein)
#   xrng <- range(b$phi)
#   yrng <- range(b$psi)
#   bw.x <- bandwidth.nrd(b$phi)
#   bw.y <- bandwidth.nrd(b$psi)
#   d1 <- kde2d(b$phi, b$psi,  lims=c(xrng, yrng), n=400, h=c(bw.x,bw.y))
#   d1_norm <- d1
#   if(normalize_to_max==T){
#     d1_norm$z <- d1_norm$z/max(d1_norm$z) 
#   }
#   
#   ## Melt data into long format
#   # First, add row and column names (x and y grid values) to the z-value matrix
#   rownames(d1_norm$z) <- d1_norm$x
#   colnames(d1_norm$z) <- d1_norm$y
#   d1_norm <- melt(d1_norm$z, id.var=rownames(d1_norm))
#   names(d1_norm) <- c("phi","psi","z")
#   
#   p <- ggplot() +
#     geom_tile(aes(x=phi, y=psi, fill=z),data=d1_norm) +
#     stat_contour(aes(x=phi,y=psi,z=z,colour=after_stat(level)),binwidth=0.001,data=d1_norm) +
#     scale_fill_gradient2(low="red",mid="white", high=color, midpoint=0) +
#     scale_colour_gradient2(low=alpha("red"),mid="white",high=color,midpoint=0) +
#     guides(colour='none') + 
#     prism() + 
#     theme(legend.position = 'bottom') +
#     scale_y_continuous(breaks = c(150,100,50,0,-50,-100,-150), limits = c(-180,180)) + 
#     scale_x_continuous(breaks = c(150,100,50,0,-50,-100,-150), limits = c(-180,180)) +
#     ggnewscale::new_scale_colour() +
#     geom_point(aes(x=phi, y=psi), data=b, pch='') + #add raw data to create density plots with ggMarginal
#     scale_colour_manual(values=alpha(color,0.8))
#   
#   p <- ggExtra::ggMarginal(p, type = "density",groupColour = F , size = 3) 
#   
#   if(plot==T){
#     print(p)
#   }
#   return(list(angles=b, kde2d=d1, denisty_matrix=d1_norm, plot=p))
# }
# 
# plot_ramachandran_difference <- function(SWaFi_score_up,
#                                          SWaFi_score_down,
#                                          plot=T,
#                                          up_color="#6082B6",
#                                          down_color='red'){
#   
#   b <- do.call(rbind,pdb_structures(SWaFi_score_up))
#   b  <- na.omit(b)
#   
#   nb <- do.call(rbind,pdb_structures(SWaFi_score_down))
#   nb  <- na.omit(nb)
#   
#   b$identifier  <- 'up'
#   nb$identifier <- 'down'
#   
#   # Calculate the common phi and psi ranges
#   b  <- na.omit(b) #remove NAs (occurs if the terminal contains either the first or last residue of a protein)
#   nb <- na.omit(nb)
#   
#   xrng <- range(c(b$phi, nb$phi))
#   yrng <- range(c(b$psi, nb$psi))
#   
#   # Calculate the 2d density estimate over the common range
#   bw.x <- bandwidth.nrd(c(b$phi,nb$phi))
#   bw.y <- bandwidth.nrd(c(b$psi,nb$psi))
#   
#   d1 <- kde2d(b$phi , b$psi,  lims=c(xrng, yrng), n=400, h=c(bw.x,bw.y))
#   d2 <- kde2d(nb$phi, nb$psi, lims=c(xrng, yrng), n=400, h=c(bw.x,bw.y))
#   
#   # Confirm that the grid points for each density estimate are identical
#   identical(d1$x, d2$x) # TRUE
#   identical(d1$y, d2$y) # TRUE
#   
#   d1_norm <- d1
#   d2_norm <- d2
#   
#   # Normalize density to be from 0-1
#   d1_norm$z <- d1_norm$z/max(d1_norm$z) 
#   d2_norm$z <- d2_norm$z/max(d2_norm$z)
#   
#   # Calculate the difference between the 2d density estimates
#   diff   <- d1_norm 
#   diff$z <- d1_norm$z - d2_norm$z
#   
#   ## Melt data into long format
#   # First, add row and column names (x and y grid values) to the z-value matrix
#   rownames(diff$z) <- diff$x
#   colnames(diff$z) <- diff$y
#   diff <- melt(diff$z, id.var=rownames(diff))
#   names(diff) <- c("phi","psi","z")
#   
#   # Plot density difference between binding and non-binding 
#   p <- ggplot() +
#     geom_tile(aes(x=phi, y=psi, fill=z),data=diff) +
#     stat_contour(aes(x=phi,y=psi,z=z,colour=after_stat(level)),binwidth=0.001,data=diff) +
#     scale_fill_gradient2(low=down_color,mid="white", high=up_color, midpoint=0) +
#     scale_colour_gradient2(low=alpha(down_color),mid="white",high=alpha(up_color),midpoint=0) +
#     guides(colour='none') + 
#     prism() + 
#     theme(legend.position = 'bottom') +
#     scale_y_continuous(breaks = c(150,100,50,0,-50,-100,-150), limits = c(-180,180)) + 
#     scale_x_continuous(breaks = c(150,100,50,0,-50,-100,-150), limits = c(-180,180)) +
#     ggnewscale::new_scale_colour() +
#     geom_point(aes(x=phi, y=psi, col=identifier), data=rbind(b,nb), pch='') + #add raw data to create density plots with ggMarginal
#     scale_colour_manual(values=c(alpha(down_color,0.8),alpha(up_color,0.8)))
#   
#   # add marginal density plots and save 
#   p <- ggExtra::ggMarginal(p, type = "density", groupColour = T, size = 3) 
# 
#   if(plot==T){
#     print(p)
#   }
#   
#   return(list(angles=list(b,nb), kde2d=list(d1,d2), denisty_matrix=diff, plot=p))
#   
# }

plot_SWaFi_trace <- function(array_object=NULL, 
                             dataset=NULL,
                             proteins=NULL,
                             ymax=NULL,
                             estimate_backgroud=F,
                             bootstraps=1000, 
                             Plot=T){
  
  traces <- get_SWaFi_traces(array_object,dataset,proteins)
  ids <- names(traces)
  
  if(is.null(ymax)){
    ymax <- max(unlist(lapply(traces,function(x){max(x$score)})))
  }
  cat('\nCreating plots ------------------------- \n')
  trace_plots <- pblapply(ids,function(id){
    trace <- traces[[id]]
    trace$id <- '0'
    max.x <- ifelse(nrow(trace)>75,round(max(trace$position),-2),round(max(trace$position),-1))
    min.x <- ifelse(nrow(trace)>75,round(min(trace$position),-2),round(min(trace$position),-1))
    max.x <- ifelse(max.x<max(trace$position),max(trace$position),max.x)
    min.x <- ifelse(min.x>min(trace$position),min(trace$position),min.x)
    
    cols <- 'black'
    
    if(estimate_backgroud==T){
      nthresh <- get_SWaFi_traces(array_object,dataset,parameters=T)$SWaFi_parameters$nthresh
      # Bootstrap sample 1000 sequences of same length as protein
      sampling <- sample_dataset(dataset,nrow(trace)-7,bootstraps, array_object,verbose = F)
      sampling <- lapply(sampling,.get_signal_fast,nthresh=nthresh)
      names(sampling) <- 1:length(sampling)
      sampling <- lapply(names(sampling),function(x){
        score <- sampling[[x]]
        df <- data.frame(score=score,position=trace$position,id=x)
        rownames(df) <- NULL
        return(df)
      })
      sampling <- do.call(rbind,sampling)
      trace <- rbind(trace,sampling)
      cols <- c('red',rep(alpha('grey75',.05),bootstraps))
    }
    
    p <- ggplot(trace, aes(x=position, y=score, col=id)) + 
           geom_line()
    if(estimate_backgroud==F){
      p <- p + geom_point(pch=16)
    }
    
    p <- p +
      ylim(0,ymax) +
      xlim(min.x,max.x) +
      scale_color_manual(values = cols) +
      prism() + 
      theme(aspect.ratio = .5, legend.position = 'none') +
      ylab('SWaFi score (a.u.)') + 
      xlab('Position') + 
      ggtitle(label = gsub(' \\(.*','',unlist(get_metadata(array_object, proteins=id, colnames = 'Protein.name'))),
              subtitle = paste0('(ArrayID: ',id,')'))
    
    return(p)
  })
  
  nCol <- floor(sqrt(length(trace_plots)))
  if(Plot==T){
    plot_list_of_ggplots(trace_plots, ncol = nCol)
  } else {
    return(trace_plots)
  }
  
}

plot_SWaFi_trace_w_sites <- function(array_object=NULL, 
                                     dataset=NULL, 
                                     proteins=NULL, 
                                     ymax=NULL, 
                                     confidence='low', 
                                     nCol=NULL){
  
  #traces <- get_SWaFi_traces(array_object,dataset,proteins)
  
  traces <- unlist(lapply(dataset,function(x){
    s <- quiet(get_SWaFi_traces(array_object = array_object, dataset = x, proteins = proteins))
    
    if(x=='FBBE1'){
      s <- lapply(s,function(x){x$score <- x$score*2.14912; return(x)})
    }
    
    return(s)
  }), recursive=F)
  
  traces <- split(traces,unique(names(traces)))
  
  traces <- lapply(traces,function(s){data.frame(score=Reduce('+',s)$score, position=s[[1]]$position)})
  
  ids <- rev(names(traces))
  
  ratio <- rev(unlist(lapply(traces,nrow)))
  
  if(is.null(ymax)){
    ymax <- max(unlist(lapply(traces,function(x){max(x$score)})))
  }
  
  sites <- quiet(get_binding_sites(proteins = proteins, binding = T, datasets = c('FBBE1','FBBE2'), confidence = confidence, array_object = array_object))   
  
  trace_plots <- lapply(ids,function(id){
    trace <- traces[[id]]
    
    trace$id <- '0'
    max.x <- ifelse(nrow(trace)>75,round(max(trace$position),-2),round(max(trace$position),-1))
    min.x <- ifelse(nrow(trace)>75,round(min(trace$position),-2),round(min(trace$position),-1))
    max.x <- ifelse(max.x<max(trace$position),max(trace$position),max.x)
    min.x <- ifelse(min.x>min(trace$position),min(trace$position),min.x)
    
    cols <- 'black'
    
    site <- sites[sites$UniProtKB.AC==id,]
    
    p <- ggplot()
    
    if(nrow(site)>0){
      p <- p + geom_rect(data = site, aes(xmin=start, xmax=end, ymin=rep(-Inf,nrow(site)), ymax=rep(Inf,nrow(site)), fill=1:nrow(site)), alpha=0.25)
    } else {
      p <- p + geom_rect(aes(xmin=min.x, xmax=max.x, ymin=-Inf, ymax=Inf), alpha=0)
    }
    
    p <- p + geom_line(data = trace, aes(x=position, y=score), col = '#0275d8') +
      geom_point(data = trace, aes(x=position, y=score), pch=16, col = '#0275d8') +
      ylim(0,ymax) +
      xlim(min.x,max.x) +
      scale_fill_viridis() +
      prism() + 
      theme(aspect.ratio = NULL, legend.position = 'none', panel.background = element_blank()) +
      ylab('SWaFi score (a.u.)') + 
      xlab('Position') 
    
    return(p)
  })
  
  if(length(ids)==2){
    trace_plots[[2]] <- quiet(trace_plots[[2]] + scale_fill_viridis(direction = -1) + remove_y_axis())
  } else {
    trace_plots[[2]] <- NULL
  }
  
  if(is.null(nCol)){
    nCol <- floor(sqrt(length(trace_plots)))
  }
  
  plot_grid(plotlist = trace_plots, rel_widths = ratio, align = 'hv', ncol = nCol) 
}

plot_SWaFi_regions <- function(array_object=NULL, dataset, protein){
  
  if(is.null(array_object)){
    array_object <- Find_array_object()
  }
  
  regions <- get_SWaFi_regions(array_object,dataset,protein)
  ids <- names(regions)
  
  region_plots <- lapply(ids,function(id){
    region <- regions[[id]]
    
    trace <- min(region$start):max(region$end)
    max.x <- ifelse(length(trace)>75,round(max(trace),-2),round(max(trace),-1))
    min.x <- ifelse(length(trace)>75,round(min(trace),-2),round(min(trace),-1))
    max.x <- ifelse(max.x<max(trace),max(trace),max.x)
    min.x <- ifelse(min.x>min(trace),min(trace),min.x)
    
    all  <- IRanges(region$start,region$end,names = region$pval)
    bins <- disjointBins(IRanges(start(all), end(all) + 1))
    dat  <- cbind(as.data.frame(all), bin = bins)
    
    plot_title <- gsub(' \\(.*','',unlist(get_metadata(array_object, proteins=id, colnames = 'Protein.name')))
    
    ggplot(dat) + 
      geom_rect(aes(xmin = start, xmax = end,
                    ymin = bin, ymax = bin + 0.9, fill=as.numeric(names))) +
      prism() + 
      theme(aspect.ratio = .2) + 
      remove_y_axis() + 
      xlim(min.x,max.x) +
      scale_fill_viridis_c() +
      xlab('Position') + 
      labs(fill='p-value') +
      ggtitle(label = plot_title,
              subtitle = paste0('(ArrayID: ',id,')'))
    
  })
  
  nCol <- floor(sqrt(length(region_plots)))
  do.call("grid.arrange", c(region_plots, ncol=nCol))
  
}

plot_SWaFi_score <- function(array_object = NULL, 
                             dataset  = NULL, 
                             proteins = NULL, 
                             pthresh  = NULL, 
                             cut = NA, 
                             get_min = F, 
                             lcut = NULL, 
                             pAdjust = NULL){
  res <- get_SWaFi_score(array_object,dataset,proteins,pthresh,cut,get_min,lcut,pAdjust,Plot=T)
  return(res)
}

# plot_qq_plot <- function(qq, 
#                          viridis_color='magma', 
#                          begin = 0.1, 
#                          end = 0.8, 
#                          alpha=0.3){
#   qq_plot <- ggplot(qq) + 
#     geom_point(aes(x=Empirical, y=Theoretical, col=as.factor(Terminal_length)), alpha=alpha, pch=16) + 
#     prism() +
#     theme(legend.position='none') +
#     geom_abline(slope=0,intercept = 0, lty=2) + 
#     geom_abline(slope=1,intercept = 0, col='red', lty=2, linewidth=1) + 
#     scale_color_viridis_d(option=viridis_color, begin = begin, end = end)
#   
#   qq_plot
# }
# 
# plot_SWaFi_matrix_column_summaries <- function(SWaFi_matrix,
#                                                name,
#                                                ymax=.5,
#                                                colFunction=colMeans) {
#   
#   regions <- gsub('-.*','',colnames(SWaFi_matrix))
#   
#   N  <- length(regions[regions=='N'])
#   TM <- length(regions[regions=='TM'])
#   C  <- length(regions[regions=='C'])
#   
#   if(N > 0){N_position <- seq(-(N*2)+1,-1,2)} else {N_position <- NULL}
#   if(C > 0){C_position <- seq(1,2*C,2)} else {C_position <- NULL}
#   
#   function_score <- colFunction(SWaFi_matrix,na.rm=T)
#   function_score[is.na(function_score)] <- 0
#   
#   df <- data.frame(position = c(N_position,0,rep(NA,TM-2),0,C_position),
#                    average_score = function_score, 
#                    terminal = as.factor(c(rep(1,N+1),rep(2,TM-2),rep(3,C+1))))
#   
#   # Plot
#   ggplot(na.omit(df), aes(x=position,y=average_score,col=terminal)) + 
#     geom_line(linewidth=1) + 
#     prism() + # custom theme
#     ggtitle(name) + 
#     facet_wrap(~terminal, ncol = 2, scales = 'free_x') + 
#     ylim(0,ymax) + 
#     theme(strip.text.x = element_blank(), 
#           axis.text = element_text(size=8), 
#           axis.title = element_text(size=10), 
#           title = element_text(size=14), 
#           legend.position = 'none') + 
#     scale_color_viridis_d(option = 'turbo') + 
#     ylab('Score')
# }

# plot_panther_overrepresentation <- function(overrepresentation_results, 
#                                             over_represented=T, 
#                                             use_corrected=T, 
#                                             aspect_ratio=NULL,
#                                             xlim=NULL, 
#                                             palette='viridis', 
#                                             alpha=0.8,
#                                             min_point_size=1,
#                                             max_point_size=10,
#                                             plot=T){
#   
#   if(use_corrected==F){
#     overrepresentation_results[,4] <- overrepresentation_results$pval
#   } 
#   
#   colnames(overrepresentation_results)[4] <-   'p_value'
#   
#   selected_sign <- ifelse(over_represented==T,'+','-')
#   overrepresentation_results <- overrepresentation_results[overrepresentation_results$`+/-`==selected_sign,]
#   
#   if(is.null(aspect_ratio)){
#     bins <- sort(riffle(ntile(x = 1:1000, n=200)-0.5,ntile(x = 1:1000, n=200)))
#     aspect_ratio <- bins[nrow(overrepresentation_results)]
#     if(aspect_ratio>10){
#       input <- readline("You're trying to plot a lot of terms are yo sure you want to continue? (options: 'y' or 'n')\n")
#       if(input=='n'){
#         stop()
#       }
#     }
#   }
#   
#   overrepresentation_results$term_label <- make.unique(overrepresentation_results$term_label)
#   overrepresentation_results <- na.omit(overrepresentation_results)
#   overrepresentation_results$database <- toTitleCase(overrepresentation_results$database)
#   
#   p <- ggplot(overrepresentation_results, aes(y=factor(term_label,levels = term_label),
#                                               x=-log10(p_value))) + 
#     geom_point(aes(size=fold_change, col=database), pch=16) + 
#     prism() + 
#     theme(aspect.ratio = aspect_ratio, 
#           legend.position = 'right', 
#           legend.justification = 'left', 
#           axis.title.y = element_blank()) +
#     scale_y_discrete(labels=gsub('\\..*','',overrepresentation_results$term_label)) + 
#     guides(col=guide_legend(ncol=1,byrow=F)) + 
#     scale_fill_viridis_d(option = palette, alpha = alpha) + 
#     scale_size_continuous(range = c(min_point_size,max_point_size))
#   
#   if(use_corrected==F){
#     p <- p + xlab(bquote(-log[10](p-value))) 
#   } else {
#     p <- p + xlab(bquote(-log[10](q-value))) 
#   }
#   
#   if(!is.null(xlim)){
#     p <- p + xlim(0,xlim)
#   } else {
#     p <- p + xlim(0,ceiling(max(-log10(overrepresentation_results$p_value)+5)))
#   }
#   
#   if(plot==T){
#     print(p)
#   }
#   
#   return(p)
# }
# 
# plot_panther_fold_change_vs_pvalue <- function(overrepresentation_results, 
#                                                use_corrected=T, 
#                                                show_top=10, 
#                                                thresholds=c(1,0.05), 
#                                                add_to_ylim_to_fit_labels=c(-log10(0.035),0), 
#                                                add_to_xlim_to_fit_labels=c(30,30),
#                                                change_distance_from_point_to_text=c(0,0),
#                                                palette='viridis', 
#                                                alpha=0.8, 
#                                                aspect.ratio=0.5){
#   
#   if(use_corrected==F){
#     overrepresentation_results[,4] <- overrepresentation_results$pval
#   } 
#  
#   colnames(overrepresentation_results)[4] <-   'p_value'
#   overrepresentation_results$fold_change <- log2(overrepresentation_results$fold_change)
#   overrepresentation_results$fold_change[overrepresentation_results$fold_change==-Inf] <- 0
#   
#   smallest <- floor(min(overrepresentation_results$fold_change))
#   biggest  <- ceiling(max(overrepresentation_results$fold_change))
# 
#   overrepresentation_results$label <- overrepresentation_results$term_label
#   overrepresentation_results$label[overrepresentation_results$p_value>thresholds[2] | 
#                                    abs(overrepresentation_results$fold_change) < thresholds[1]] <- ""
#   
#   if(!is.null(show_top)){
#     top_right <- sort(overrepresentation_results$p_value[overrepresentation_results$fold_change>thresholds[1]])[show_top]
#     top_left <-  sort(overrepresentation_results$p_value[overrepresentation_results$fold_change<(-thresholds[1])])[show_top]
#     overrepresentation_results$label[overrepresentation_results$p_value>top_right & overrepresentation_results$fold_change>thresholds[1] | 
#                                      overrepresentation_results$p_value>top_left  & overrepresentation_results$fold_change<(-thresholds[1])] <- ""
#   }
#   
#   OR_res_no_labels <- overrepresentation_results[overrepresentation_results$label=="",]
#   OR_res_no_labels$label_position.x <- NA
#   OR_res_no_labels$label_position.y <- NA
#   
#   OR_res_labels <- overrepresentation_results[overrepresentation_results$label!="",]
#   OR_res_labels$label_position.x <- 0
#   OR_res_labels$label_position.x[OR_res_labels$fold_change<0] <- smallest+change_distance_from_point_to_text[1]
#   OR_res_labels$label_position.x[OR_res_labels$fold_change>0] <- biggest+change_distance_from_point_to_text[2]
#   
#   smallest_p <- ceiling(-log10(min(overrepresentation_results$p_value))) + add_to_ylim_to_fit_labels[2]
#   
#   OR_res_labels <- OR_res_labels |> 
#     group_by(label_position.x) |> 
#     arrange(desc(p_value)) |> 
#     mutate(label_position.y=seq((0+add_to_ylim_to_fit_labels[1]),smallest_p,((smallest_p-add_to_ylim_to_fit_labels[1])/(length(p_value)-1)))) |> 
#     ungroup()
#   
#   overrepresentation_results <- rbind(OR_res_no_labels,OR_res_labels)
#   
#   if(smallest_p<ceiling(-log10(min(overrepresentation_results$p_value)))){
#     smallest_p <- ceiling(-log10(min(overrepresentation_results$p_value)))
#   }
#   
#   overrepresentation_results$database[overrepresentation_results$p_value>thresholds[2] | abs(overrepresentation_results$fold_change)<thresholds[1]] <- NA
#   
#   p <- ggplot(overrepresentation_results, aes(y=-log10(p_value),x=fold_change,col=database)) + 
#     geom_text(data=overrepresentation_results[overrepresentation_results$label_position.x>0,], aes(label = paste0(label,' p=',format(round(p_value, 5), nsmall = 5)), x=label_position.x, y=label_position.y), show.legend = FALSE, hjust = 0) +
#     geom_text(data=overrepresentation_results[overrepresentation_results$label_position.x<0,], aes(label = paste0('p=',format(round(p_value, 5), nsmall = 5),' ',label), x=label_position.x, y=label_position.y), show.legend = FALSE, hjust = 1) +
#     geom_point(alpha = alpha) + 
#     prism() + 
#     theme(aspect.ratio = aspect.ratio, 
#           legend.position = 'bottom', 
#           legend.justification = 'left') +
#     ylab(bquote(-log[10](p-value))) +
#     xlab(bquote(log[2](FoldChangeOverExpected))) +
#     guides(col=guide_legend(ncol=2,byrow=F)) + 
#     scale_color_manual(na.value = 'grey80', values = setNames(viridis(length(unique(overrepresentation_results$database)),option = palette),unique(overrepresentation_results$database)) ) + 
#     geom_hline(yintercept = -log10(thresholds[2]), lty=2, col='grey30') + 
#     geom_vline(xintercept = c(-thresholds[1],0,thresholds[1]), lty=2, col='grey30') + 
#     xlim(smallest-add_to_xlim_to_fit_labels[1],biggest+add_to_xlim_to_fit_labels[2]) + 
#     ylim(0,smallest_p)
# 
#     suppressWarnings(print(p))
# }

SWaFi_matrix_dendrogram <- function(SWaFi_matrix, 
                                    distance_measure='euclidean', 
                                    hclust_method='ward.D'){
  SWaFi_matrix[is.na(SWaFi_matrix)] <- 0
  d <- dist(SWaFi_matrix, method = distance_measure)
  hc <- hclust(d, method = hclust_method)
  dendrogram <- as.dendrogram(hc)
  return(dendrogram)
}

plot_SWaFi_matrix <- function(SWaFi_matrix,
                              hclust_method='ward.D',
                              distance_measure='euclidean',
                              colorRamp2_function=NULL, 
                              brewerpal_color='YlGnBu', 
                              na_color = 'white',
                              non_binding_color='lavender',
                              dendrogram=NULL,
                              dend_width=5,
                              show_axis=T,
                              add_column_labels=T,
                              column_labels_interval=50,
                              row_stats=NULL,
                              bar_width=0.01,
                              border=T,
                              show_rownames=F,
                              legend=F,
                              use_raster=T,
                              row_split=NULL,
                              ylim=NULL,
                              show_barplot_label=F,
                              heatmap_name_rot=90,
                              heatmap_name=NULL){
  
  # Define color scale ---------------------------------------------------------
  if(is.null(colorRamp2_function)){
    pos <- brewer.pal(5,brewerpal_color)
    colorRamp2_function <- colorRamp2(c(-1,-0.999999999,0,0.1,0.2,0.3,0.4,0.5,0.6), 
                                      c('white','white',non_binding_color,'grey40',pos))
    
    }
  
  # Hierarchical clustering of matrix ------------------------------------------
  heat_matrix <- SWaFi_matrix
  heat_matrix[is.na(heat_matrix)] <- 0
  
  if(!is.logical(dendrogram)){
    if(is.null(dendrogram)){
      d <- dist(heat_matrix, method = distance_measure)
      hc <- hclust(d, method = hclust_method)
      dendrogram <- as.dendrogram(hc)
    } 
  }

  
  # Add back lengths to matrix by re-introducing NAs ---------------------------
  mask_matrix <- SWaFi_matrix
  mask_matrix[!is.na(mask_matrix)] <- 1
  heat_matrix <- heat_matrix*mask_matrix
  

  # Add column labels ----------------------------------------------------------
  if(add_column_labels==T){
    pos1 <- find_SWaFi_matrix_TM_segment(SWaFi_matrix)$end
    pos_at <- seq(pos1,ncol(SWaFi_matrix),by=column_labels_interval)
    pos_labels <- pos_at-pos_at[1]
    pos_at[1] <- pos_at[1] + 1
    pos_labels <- pos_labels*2
    pos_labels[1] <- pos_labels[1] + 1
    
    neg1 <- find_SWaFi_matrix_TM_segment(SWaFi_matrix)$start
    neg_at <- seq(1,neg1,by=column_labels_interval)
    neg_at <- neg_at+neg1-max(neg_at)
    neg_labels <- neg_at[length(neg_at)]-neg_at
    neg_labels <- neg_labels*2
    neg_at[length(neg_at)] <- neg_at[length(neg_at)] - 1
    neg_labels[length(neg_at)] <- neg_labels[length(neg_at)] + 1
    neg_labels <- -neg_labels
    
    at <- c(neg_at,pos_at)
    label <- c(neg_labels,pos_labels)
    
    ha <- columnAnnotation(position = anno_mark(at = at, 
                                                labels = label, 
                                                labels_rot = 0, link_height = unit(0.25,'cm'),
                                                side = 'bottom',
                                                extend = 0.5))
    } else {ha <- NULL}
  
  
  if(!is.null(row_stats)){
    
    if(!is.function(row_stats)){
      cat("\nERROR: 'row_stats' has to be a function that can be applied to the rows of a matrix\n")
      stop()
    }
    
    hr = rowAnnotation(f=anno_barplot(x = row_stats(SWaFi_matrix,na.rm=T), show_name=F,ylim = ylim,
                       bar_width = bar_width,border = F,add_numbers = F, width = unit(dend_width,'cm'), axis = show_axis, 
                       gp = gpar(col = alpha('black',0.3),fontsize=0), 
                       axis_param = list(side = "top",labels_rot=0)),
                       show_annotation_name=show_barplot_label)
    
  } else {hr <- rowAnnotation(foo = anno_empty(border = F, show_name = F,width = unit(1,'cm')))}
  
   ht_opt$message = FALSE #turn of RStudio warning
  
  # Plot -----------------------------------------------------------------------
  ht <- Heatmap(SWaFi_matrix,
                cluster_rows = dendrogram,
                show_row_dend = T,
                cluster_columns = F, 
                show_column_names = F, 
                show_row_names = show_rownames, 
                col = colorRamp2_function, 
                na_col = na_color, 
                row_dend_width = unit(dend_width, "cm"), 
                use_raster = use_raster, 
                show_heatmap_legend = legend, 
                bottom_annotation = ha, 
                right_annotation = hr, 
                border = border,
                row_split = row_split, 
                row_title = heatmap_name,
                row_title_rot = heatmap_name_rot)
  
  return(ht)
}

plot_SWaFi_matrix_subsets <- function(SWaFi_matrix,
                                      protein_list,
                                      heatmap_names=F,...){
  cluster_heatmaps <- NULL # Empty object
  
  if(heatmap_names==T & is.null(names(protein_list))){
    cat("\nThe 'protein_list' does not have names, please add names or set 'heatmap_names'=FALSE ... \n")
    stop()
  }
  
  for(i in seq_along(protein_list)){
    
    # Subset matrix
    sub_mtrx <- subset_SWaFi_matrix(SWaFi_matrix,protein_list[[i]])
    sub_mtrx[,find_SWaFi_matrix_TM_segment(sub_mtrx)$columns] <- NA
    
    if(heatmap_names==T){
      heatmap_name <-names(protein_list)[i]
    } else {
      heatmap_name <- NULL
    }
    
    # Plot 
    ht <- plot_SWaFi_matrix(sub_mtrx,
                            add_column_labels=ifelse(i==length(protein_list),T,F),
                            show_axis=ifelse(i==1,T,F),heatmap_name=heatmap_name,...)
    
    # Create plotable ComplexHeatmap object
    cluster_heatmaps <- cluster_heatmaps %v% ht
  }
  return(cluster_heatmaps)
}

keep_regions_in_SWaFi_matrix <- function(SWaFi_matrix, SWaFi_score){
  
  TM_seg <- find_SWaFi_matrix_TM_segment(SWaFi_matrix)
  
  mask_matrix <- mask_proteins_in_SWaFi_matrix(SWaFi_matrix,rownames(SWaFi_matrix))
  
  SWaFi_score <- as.data.frame(SWaFi_score)
  
  for(i in 1:nrow(SWaFi_score)){
    if(SWaFi_score[i,'terminal']=='N'){
      
      last <- TM_seg$start-floor(SWaFi_score[i,'distance_from_TM']/2)
      first <- last-ceiling(SWaFi_score[i,'length']/2)
      
      if(first<1){first <- 1}
      
      if(last>0){
        mask_matrix[rownames(mask_matrix)==strip_terminal_specifier(SWaFi_score[i,'UniProtKB.AC']), first:last] <- mask_matrix[rownames(mask_matrix)==strip_terminal_specifier(SWaFi_score[i,'UniProtKB.AC']), first:last]+1
      }
    } else {
      
      first <- TM_seg$end+ceiling(SWaFi_score[i,'distance_from_TM']/2)
      last <- first+floor(SWaFi_score[i,'length']/2)-1
      
      if(last>ncol(mask_matrix)){last <- ncol(mask_matrix)}
      
      if(first<ncol(mask_matrix)+1){
        mask_matrix[rownames(mask_matrix)==strip_terminal_specifier(SWaFi_score[i,'UniProtKB.AC']), first:last] <- mask_matrix[rownames(mask_matrix)==strip_terminal_specifier(SWaFi_score[i,'UniProtKB.AC']), first:last]+1
      }
    }
  }
  
  mask_matrix[mask_matrix>1] <- 1
  masked <- SWaFi_matrix*mask_matrix
  return(masked)
}

# plot_AA_probability_density_fold_change <- function(positive_proteins = NULL,
#                                                     negative_proteins = NULL,
#                                                     positive_sequences = NULL,
#                                                     negative_sequences = NULL, 
#                                                     plot = T, 
#                                                     color_vector_with_names=NULL, 
#                                                     ylimits=NULL){
#   
#   pos <- get_AA_frequency(sequences = positive_sequences, proteins = positive_proteins, freq=F)
#   neg <- get_AA_frequency(sequences = negative_sequences, proteins = negative_proteins, freq=F)
#   
#   logFC <- data.frame(log2(pos/neg))
#   colnames(logFC) <- c('AA','log2FC')
#   
#   logFC <- logFC[order(logFC$log2FC),]
#   rownames(logFC) <- NULL
#   
#   if(is.null(color_vector_with_names)){
#     color_vector_with_names <- setNames(c('#ffdb58','blue','purple','red','#ffdb58','#ffdb58','purple',
#                                           'red','cyan','blue','#ffdb58','#ffdb58','blue','#ffdb58',
#                                           '#ffdb58','cyan','purple','purple','#ffdb58','#ffdb58','#ffdb58'),
#                                         c("A", "R", "N", "D", "C", "U" ,"Q", "E",
#                                           "G", "H", "I", "L", "K", "M", "F",
#                                           "P", "S", "T", "W", "Y", "V"))
#   }
#   
#   if(is.null(ylimits)){
#     ylimits <- c(-max(abs(logFC$log2FC)),max(abs(logFC$log2FC)))
#   }
#   
#   p <- ggplot(logFC,aes(x=reorder(AA,log2FC),y=log2FC, fill=reorder(AA,log2FC))) + 
#     geom_bar(stat='identity') + 
#     scale_fill_manual(values=alpha(color_vector_with_names,.75)) + 
#     prism()  + theme(legend.position='none') +
#     xlab('Amino acid 1-letter code') + 
#     ylab('log2 Fold change "log2(positive/negative)"') + ylim(ylimits[1],ylimits[2])
#   
#   if(plot==T){
#     print(p)
#   }
#   return(list(log2FC=logFC,plot=p))
# }

plot_list_of_ggplots <- function(plot_list,ncol=1){
  grid.draw(do.call("grid.arrange", c(plot_list, ncol=ncol)))
}

# used in plot_pdb_with_binding_sites function
seq2 <- Vectorize(seq.default, vectorize.args = c("from", "to"))

plot_pdb_with_binding_sites <- function(UniProtKB.AC,
                                        SWaFi_score=NULL,
                                        pthresh=0.05,
                                        protein_color="#d3d3d3",
                                        protein_opacity=0.5,
                                        site_color='red',
                                        AlphaFold_directory=NULL,
                                        plot_all_scores=T,
                                        array_object=NULL,
                                        dataset=NULL,
                                        zoom_limits=c(10,500),
                                        quality=5,
                                        colorRamp2_function=NULL,
                                        brewerpal_color='YlGnBu',
                                        non_binding_color='lavender',
                                        cartoon_style='oval',
                                        background_color='white'){
  
  #library(r3dmol)
  #library(bio3d)
  
  if(is.null(array_object)){
    array_object <- Find_array_object()
  }
  
  if(!is.null(SWaFi_score) & plot_all_scores==T){
    cat('Provided SWaFi_scores will be plotted for significant regions only because plot_all_scores=T ... \n')
    if(is.null(dataset)){
      cat('No specific dataset requested defaulting to FBBE2 ... \n')
      dataset <- 'FBBE2'
    }
    only_plot_scores_in_sites <- SWaFi_score
  }
  
  # Get default binding scores if not provided 
  if(is.null(SWaFi_score)){
    cat('No SWaFi scores provided using default settings for specified dataset ... \n')
    if(is.null(dataset)){
      cat('No specific dataset requested defaulting to FBBE1 + FBBE2 ... \n')
      dataset <- c('FBBE1','FBBE2')
    }
    only_plot_scores_in_sites <- SWaFi_score
    SWaFi_score <- get_binding_sites(dataset = dataset, confidence = pthresh, proteins = UniProtKB.AC, array_object = array_object)
  }
  
  # Identify the type of id provided and subset scores
  if(detect_identifier(UniProtKB.AC)=='ArrayID'){
    SWaFi_score <- SWaFi_score[SWaFi_score$UniProtKB.AC==UniProtKB.AC,]
  } else {
    SWaFi_score <- SWaFi_score[strip_terminal_specifier(SWaFi_score$UniProtKB.AC)==UniProtKB.AC,]
  }
  
  # # Set p-value threshold
  # SWaFi_score <- SWaFi_score[SWaFi_score$pval<=pthresh,]
  
  # Define amino acids residues (binding site) to highlight in structure
  if(nrow(SWaFi_score)<1){
    highlight <- c(0)
  } else {
    highlight <- as.vector(unlist(seq2(from = SWaFi_score$start, to = SWaFi_score$end, by = 1)))
  }
  
  # Decide whether to use provided directory of structures or the data stored in the array_object
  if(!is.null(AlphaFold_directory)){
    pdb_files <- list.files(path=AlphaFold_directory, pattern = '.pdb', full.names = T)
  
    if(length(pdb_files)<1){
      cat('To use this function you need to have downloaded the pdb files needed ... \n Maybe you made a mistake in the path? ... \n')
      stop()
    }
  
    pdb_file <- pdb_files[str_detect(pdb_files,strip_terminal_specifier(UniProtKB.AC))]
    pdb_file <- pdb_file[1]
    
    # Test if structure file exists
    if(length(pdb_file)<1){
      cat('File not found in provided directory ... \n')
      pdb <- NA
    } else {
      pdb <-read.pdb(pdb_file) #createbio3dobject
    }
    
  } else {
    pdb <- quiet(get_metadata(array_object = array_object, colnames = 'AlphaFold_pdb', proteins = UniProtKB.AC)[[1]])
  }
  
  # Test if structure has been found
  if(length(pdb)==1){
    cat('array_object or directory does not contain a structure for this protein ... \n')
    res <- NA
  } else {
  
  # loadbio3dobjectintor3dmol 
  res <- r3dmol(viewer_spec = m_viewer_spec(
                backgroundColor = background_color,
                cartoonQuality = quality,
                lowerZoomLimit = zoom_limits[1],
                upperZoomLimit = zoom_limits[2]),
                id = "demo",
                elementId = "demo") %>% 
    m_add_model(data=m_bio3d(pdb)) %>% 
    m_zoom_to() %>%  
    m_set_style(style = m_style_cartoon(color = protein_color, opacity = protein_opacity, style=cartoon_style)) 
  
  if(plot_all_scores==T){
    
    # Define color scale ---------------------------------------------------------
    if(is.null(colorRamp2_function)){
      pos <- brewer.pal(5,brewerpal_color)
      colorRamp2_function <- colorRamp2(breaks = c(0,0.1,0.2,0.3,0.4,0.5,0.6), 
                                        colors = c(non_binding_color,'grey40',pos))
      
    }
  
    scores <- unlist(lapply(dataset,function(x){
      s <- quiet(get_SWaFi_traces(array_object = array_object, dataset = x, proteins = UniProtKB.AC))
      
      if(x=='FBBE1'){
        s <- lapply(s,function(x){x$score <- x$score*2.14912; return(x)})
      }
      
      return(s)
    }), recursive=F)
    
    scores <- split(scores,unique(names(scores)))
    
    scores <- lapply(scores,function(s){data.frame(score=Reduce('+',s)$score, position=s[[1]]$position)})
      
    scores <- lapply(scores,function(x){
                     x$color <- colorRamp2_function(x$score)
                     x$color <- substring(x$color,1,7) # for some reason colorRamp adds two F's to the end of each color??? does not work with r3dmol so we remove them
                     data.frame(position=riffle(x$position,x$position+1) ,color=riffle(x$color,x$color))
                     })
    
    len <- nchar(unlist(get_metadata(array_object = array_object, colnames = "sub.sequence", proteins = names(scores))))
    
    scores <- na.omit(do.call(rbind,lapply(seq_along(len), function(i){
                    scores[[i]] <- scores[[i]][1:len[i],]
                    })))
    
    if(!is.null(only_plot_scores_in_sites)){
      scores$color[!scores$position%in%highlight] <- non_binding_color
    }
    
    for(i in 1:nrow(scores)){
      res <- res %>% m_set_style(style = m_style_cartoon(color = scores$color[i]), sel = m_sel(resi = c(scores$position[i],scores$position[i])))
    }
    
  } else {
    res <- res %>% m_set_style(style = m_style_cartoon(color = site_color), sel = m_sel(resi = highlight))
  }
  
  res %>% m_render()
  }
  return(res)
}

plot_pdb_list <- function(pdb_list,
                          ncol=2, 
                          background_color='black', 
                          control_all=F){
  m_grid(
    viewer = pdb_list,
    rows = ceiling(length(pdb_list)/ncol),
    cols = ncol,
    control_all = control_all,
    viewer_config = m_viewer_spec(
      backgroundColor = background_color
    )
  )
}

plot_multiple_pdbs <- function(UniProtKB.AC, 
                               ncol=2, 
                               background_color='black', 
                               control_all=F, ...){
  
  plot_pdb_list(lapply(UniProtKB.AC,plot_pdb_with_binding_sites, ...), 
                ncol=ncol,
                background_color=background_color,
                control_all=control_all)

}

plot_pdb_compare_datasets <- function(UniProtKB.AC, 
                                      datasets=c('FBBE1','FBBE2'), 
                                      background='black',
                                      control_all=T,
                                      ncol=1, ...){
  
  res <- quiet(lapply(datasets, function(dataset){
    p <- plot_pdb_with_binding_sites(UniProtKB.AC = UniProtKB.AC, dataset = dataset, ...)
    return(p)
  }))
  
  plot_pdb_list(res, background_color = background, control_all = control_all, ncol=ncol)
}

plot_pdb_highlight_residues <- function(r3dmol_obj, residues=1, color='red', add_label=T){
  
  if(add_label==T){
    r3dmol_obj <- r3dmol_obj %>%
      m_add_res_labels(sel = m_sel(resi = residues), 
                       style = m_style_label(showBackground = F, 
                                             fontColor = color), 
                       byframe = F) %>% 
      m_set_style(style = m_style_cartoon(color = color), 
                  sel   = m_sel(resi = residues)) %>%
      m_render()
  } else {
    r3dmol_obj <- r3dmol_obj %>%
      m_set_style(style = m_style_cartoon(color = color), 
                  sel   = m_sel(resi = residues))  %>%
      m_render()
  }
  return(r3dmol_obj)
}

save_pdb_with_binding_sites <- function(proteins, out_folder=NULL, ...){
  
  if(is.null(out_folder)){
    cat('Provide and output folder name ... \n')
    stop()
  }
  
  dir.create(here('output',out_folder))
  dir <- here('output',out_folder)
  wd_temp <- getwd()
  setwd(dir)
  
  res <- lapply(proteins, function(i){
    p <- plot_pdb_with_binding_sites(UniProtKB.AC = i, ...)
    name <- paste0('./',i,'.html')
    
    htmlwidgets::saveWidget(p,
                            name,
                            selfcontained = TRUE)
    })
  setwd(wd_temp)
  return('Done, your folder has been created in the - output - folder')
}

plot_conservation_score <- function(proteins, 
                                    conservation_score = 'WRSRscore.smooth',
                                    terminal_only = T,
                                    no_columns = NULL,
                                    higher_conservation = 'forestgreen',
                                    lower_conservation = 'darkred',
                                    line_color = 'white',
                                    Plot = T,
                                    array_object = NULL){
  
  if(is.null(array_object)){
    array_object <- Find_array_object() 
  }
  
  if(terminal_only==T) {
    metadata_to_get <-  "terminal.conservation"
  } else {
    metadata_to_get <-  "conservation"
  }
  
  fullplotlist <- unlist(lapply(proteins, function(protein){
    
    name <- unique(get_short_protein_name(protein))
    
    conservation_list <- get_metadata(array_object = array_object,
                                      colnames = metadata_to_get,
                                      proteins = protein)
    
    
    plotlist <- lapply(conservation_list, function(df){
      
      df <- df[,c('AA','residue',conservation_score)]
      colnames(df)[3] <- 'conservation_score'
      
      p <- ggplot(data=df, aes(x = residue, y = conservation_score)) +
        geom_rect(xmin=min(df$residue), 
                  xmax=max(df$residue), 
                  ymin=0, 
                  ymax=max(df$conservation_score)+1, 
                  fill=higher_conservation) +
        geom_area(fill=lower_conservation) +
        geom_line(col=line_color, linewidth=1) +
        scale_x_continuous(breaks = df$residue,
                           labels = df$AA) +
        prism() +
        theme(legend.position = 'none', 
              aspect.ratio = .25, 
              axis.text.x = element_text(size=8)) +
        xlab("Sequence Homo Sapiens") +
        scale_y_continuous(name = "Relative Substitution Score") + 
        ggtitle(paste0(name,'\n',protein))
      
      return(p)
    })
    
    return(plotlist)
  }), recursive = F, use.names = T)
  
  if(is.null(no_columns)){
    floor(sqrt(length(fullplotlist)))
  }
  
  if(Plot==T){
    plot_list_of_ggplots(fullplotlist, ncol = no_columns)
  }
  return(fullplotlist)
}

plot_trace_and_conservation <- function(protein = NULL, 
                                        dataset = 'FBBE2', 
                                        array_object = NULL){
  
  if(is.null(array_object)){
    array_object <- Find_array_object()
  }
  
  if(is.null(protein)){
    message('Please provide a protein to plot ... ')
    stop()
  }
  
  n_terminals <- quiet(length(get_metadata(colnames = 'terminal',proteins= protein, array_object = array_object)))
  
  if(n_terminals == 2){
    
    plot_list_of_ggplots(quiet(c(rev(plot_conservation_score(proteins = protein, terminal_only = T, no_columns = 1, array_object = array_object, Plot = F)),
                                plot_SWaFi_trace(proteins = paste0(protein,'_N'), dataset = 'FBBE2', ymax = 1.5, array_object = array_object, Plot=F),
                                plot_SWaFi_trace(proteins = paste0(protein,'_C'), dataset = 'FBBE2', ymax = 1.5, array_object = array_object, Plot=F))), 
                         ncol = 2)
    
  } else {
    
    plot_list_of_ggplots(quiet(c(plot_conservation_score(proteins = protein, terminal_only = T, no_columns = 1, Plot = F, array_object = array_object),
                                 plot_SWaFi_trace(proteins = protein, dataset = 'FBBE2', ymax = 1.5, array_object = array_object, Plot=F))), 
                         ncol = 1)
    
  }
  
}

#-############################################################################-#
#-############### ADD FUNCTIONS ('add' data to array_object) #################-#
#-############################################################################-#

# add_author <- function(author, array_object=NULL){
#   if(is.null(array_object)){
#     array_object <- Find_array_object()
#   }
#   array_object$About$Authors <- paste0(array_object$About$Authors,' & ',author)
#   return(array_object)
# }
# 
# add_control_metadata <- function(dataset,metadata_list, array_object=NULL){
#   
#   if(is.null(array_object)){
#     array_object <- Find_array_object()
#   }
#   
#   if(!is.list(metadata_list)){
#     print('Error, please provide a list of metadata ... ')
#     stop()
#   }
#   
#   array_object$Experiments[[dataset]]$Controls_metadata <- metadata_list
#   
#   updated_array_object <- update_date(array_object)
#   return(updated_array_object)
# }
# 
# add_experiment <- function(array_object=NULL,
#                            normalised_data,
#                            dataset,
#                            file_path=NULL,
#                            columns=NULL,
#                            column_with_peptide_sequences=NULL){
#   
#   
#   if(is.null(array_object)){
#     array_object <- Find_array_object()
#   }
#   
#   if(is.null(file_path)){
#     cat("\nPlease provide the path to the raw file used to create this experiment entry ... \n")
#     stop()
#   }
#   
#   if(dataset %in% get_colnames(array_object)){
#     cat(paste0("\nDataset or item with name '",dataset,"' already exists...!!! \n"))
#     stop()
#   }
#   
#   if(is.null(column_with_peptide_sequences)){
#     cat('\nName of column containing peptide sequences has not been specified ... \n')
#     cat("Trying with DEFAULT name: 'coresequence'\n")
#     pep_col <- 'coresequence'
#   } else {
#     pep_col <- column_with_peptide_sequences
#   }
#   
#   if(is.null(columns)){
#     columns <- c('signal','coresequence','sector','row','col','normalised')
#     cat('\nNo column names provided ...\n')
#     cat('Trying with default names: \n')
#     cat(paste0(paste(columns,collapse = ' | '),'\n'))
#   }
#   
#   if(pep_col %notin% columns){
#     cat('\nSpecified name of column containing peptide sequences does not exist in the data to be added')
#     cat(paste0('\nName provided: ',pep_col))
#     stop()
#   }
#   
#   terminals <- get_terminals(array_object)
#   
#   if(pep_col %notin% names(terminals[[1]]$PeptideProperties)){
#     cat(paste0('\nCould not find column name <',pep_col,'> in existing data ... '))
#     cat("\nReplacing column name with: 'coresequence'")
#     colnames(normalised_data)[which(colnames(normalised_data)==pep_col)] <- 'coresequence'
#     pep_col <- 'coresequence'
#   } 
#   
#   res <- pblapply(terminals,function(y){
#     
#     newdata <- normalised_data[normalised_data[,pep_col]%in%y$PeptideProperties$coresequence,]
#     newdata <- merge(x  = y$PeptideProperties[,c(pep_col,"Order")],
#                      y  = newdata[,columns], 
#                      by = pep_col, all.x=T)
#     
#     newdata <- newdata[order(newdata$Order),]
#     
#     z <- which(names(y)=='PeptideProperties')
#     y <- c(y[1:z-1],list(newdata),y[z:length(y)])
#     names(y)[z] <- dataset
#     return(y)
#   })
#   
#   updated_array_object <- array_object
#   updated_array_object$Protein_terminals <- res
#   
#   Controls <- extract_controls(updated_array_object,normalised_data,dataset,pep_col)[,columns]
#   
#   updated_array_object$Experiments[[dataset]] <- list(Added=today(), From_file=file_path, Controls=Controls)
#   
#   updated_array_object <- update_date(updated_array_object)
#   
#   return(updated_array_object)
# }
# 
# add_fits <- function(array_object=NULL, dataset, fit_list){
#   # Messages and errors
#   if(is.null(dataset)){
#     print('Please specify which dataset you want add the fits to  ...')
#     stop()
#   }
#   
#   if(is.null(array_object)){
#     array_object <- Find_array_object()
#   }
#   
#   array_object$Experiments[[dataset]]$Fits <- fit_list
#   
#   updated_array_object <- update_date(array_object)
#   
#   return(updated_array_object)
# }
# 
# add_imputation_fit <- function(array_object=NULL, dataset, impute_fit){
#   # Messages and errors --------------------------------------------------------
#   if(is.null(dataset)){
#     print('Please specify which dataset you want add the fit to  ...')
#     stop()
#   }
#   
#   if(is.null(array_object)){
#     array_object <- Find_array_object()
#   }
#   
#   array_object$Experiments[[dataset]]$Imputation_fit <- impute_fit
#   
#   updated_array_object <- update_date(array_object)
#   
#   return(updated_array_object)
# }
# 
# add_random_sampling <- function(array_object=NULL,
#                                 dataset,
#                                 random_sampling_list){
#   
#   # Messages and errors
#   if(is.null(dataset)){
#     print('Please specify which dataset you want add the fit to  ...')
#     stop()
#   }
#   
#   if(is.null(array_object)){
#     array_object <- Find_array_object()
#   }
#   
#   array_object$Experiments[[dataset]]$Random_Sampling <- random_sampling_list
#   
#   updated_array_object <- update_date(array_object)
#   
#   return(updated_array_object)
#   
# }
# 
# add_SWaFi_traces_to_array_object <- function(array_object=NULL, 
#                                              dataset=NULL,
#                                              nthresh=0.2){
#   
#   # Messages and errors --------------------------------------------------------
#   if(is.null(dataset)){
#     print('Please specify which dataset you want add SWaFi scores to  ...')
#     stop()
#   }
#   
#   if(is.null(array_object)){
#     array_object <- Find_array_object()
#   }
#   
#   # Defining input to get_signal function --------------------------------------
#   terminals <- get_terminals(array_object=array_object)
#   impute <- get_imputation_information(array_object, dataset)
#   
#   # Run get_signal function ----------------------------------------------------
#   traces <- pblapply(terminals,
#                    .get_signal,
#                    dataset=dataset,
#                    nthresh=nthresh,
#                    impute=impute$parameters,
#                    impute_distribution=impute$Function,
#                    Plot=F)
#   
#   # Add traces to array_object and return result -------------------------------
#   array_object$Experiments[[dataset]]$SWaFi_traces <- traces
#   array_object$Experiments[[dataset]]$SWaFi_traces_parameters <- list(nthresh=nthresh)
#   
#   updated_array_object <- update_date(array_object)
#   
#   return(updated_array_object)
# }
# 
# add_SWaFi_regions_to_array_object <- function(array_object=NULL,
#                                               dataset,lcut){
#   
#   if(is.null(array_object)){
#     array_object <- Find_array_object()
#   }
#   
#   traces <- get_SWaFi_traces(array_object, dataset)
#   
#   fits   <- get_fits(array_object, dataset)
#   distribution <- fits$Function
#   fits <- fits$Fits
#   
#   nthresh <- get_SWaFi_traces(array_object,dataset, parameters=T)$SWaFi_parameters$nthresh
#   all_signals <- na.omit(get_dataset(array_object,dataset)$normalised)
#   mx <- max(.get_signal_fast(y=rep(max(all_signals),length(fits)), nthresh=nthresh)) #theoretical max
#   
#   regions <- pblapply(traces,
#                       .Define_regions_and_calculate_p_vals,
#                       lcut=lcut,
#                       fits=fits,
#                       mx=mx,
#                       distribution=distribution)
#   
#   regions <- lapply(regions,function(x){
#     x$adjusted_pval <- x$pval
#     return(x)
#   })
#   
#   array_object$Experiments[[dataset]]$SWaFi_regions <- regions
#   array_object$Experiments[[dataset]]$SWaFi_regions_parameters <- list(lcut=lcut, pAdjust='none')
#   
#   updated_array_object <- update_date(array_object)
#   print(paste0('SWaFi regions successfully added to dataset: ',dataset,'!'))
#   return(updated_array_object)  
# }
# 
# add_SWaFi_binding_sites <- function(array_object=NULL,
#                                     dataset,
#                                     pthresh,
#                                     useAdjusted=F,
#                                     get_min=F,
#                                     overwrite=F){
#   
#   if(is.null(array_object)){
#     array_object <- Find_array_object()
#   }
#   
#   lcut <- get_SWaFi_regions(array_object,dataset,parameters = T)$SWaFi_parameters$lcut
#   
#   if('SWaFi_binding_sites'%in%names(get_experiment(array_object,dataset))){
#     
#     previous_calculations <- names(get_experiment(array_object,dataset)$SWaFi_binding_sites)
#     Type <- ifelse(get_min==F,'widest','lowest')
#     calculation_id <- paste0('Threshold_',pthresh,'__Adjusted_',useAdjusted,'__Type_',Type)
#     
#     if(calculation_id%in%previous_calculations & overwrite==F){
#       print(paste0("'SWaFi_binding_sites' with these parameters already exist in dataset: '",dataset,"' ..."))
#       stop()
#     } else {
#       
#       if(overwrite==T){
#         continue <- NA
#         while(continue%notin%c('yes','n','N','no','No','NO')){
#           continue <- readline('!!!! You are overwriting data !!!! Are you sure? (options: "yes" or "no") ')
#         }
#         if(continue%in%c('n','N','no','No','NO')){
#           overwrite <- F
#           stop()
#         }
#         cat('\n')
#       }
#       
#       if(overwrite==T){
#         print('Overwriting ...')
#         print(paste0("Calculating binding sites ... calculation id: '",calculation_id,"'"))
#         binding_sites <- .Calculate_SWaFi_binding_sites(array_object,dataset,proteins=NULL,pthresh,lcut,useAdjusted,get_min)
#         array_object$Experiments[[dataset]]$SWaFi_binding_sites[[calculation_id]] <- binding_sites
#       }
#     }
#   } else {
#     
#     Type <- ifelse(get_min==F,'widest','lowest')
#     calculation_id <- paste0('Threshold_',pthresh,'__Adjusted_',useAdjusted,'__Type_',Type)
#     
#     print(paste0("Calculating binding sites ... calculation id: '",calculation_id,"'"))
#     binding_sites <- .Calculate_SWaFi_binding_sites(array_object,dataset,proteins=NULL,pthresh,lcut,useAdjusted,get_min)
#     array_object$Experiments[[dataset]][['SWaFi_binding_sites']][[calculation_id]] <- binding_sites
#     
#   }
#   
#   updated_array_object <- update_date(array_object)
#   print(paste0('SWaFi binding sites successfully added to dataset: ',dataset,'!'))
#   return(updated_array_object) 
#   
# }
# 
# add_new_metadata <- function(array_object=NULL, 
#                              metadata_list=NULL, 
#                              metadata_name=NULL, 
#                              backup=T, 
#                              overwrite=F){
#   
#   if(is.null(array_object)){
#     array_object <- Find_array_object()
#   }
#   
#   arrayID <- F
#   all_ids <- get_available_proteins(array_object)
#   if(detect_identifier(names(metadata_list))=="ArrayID"){
#     arrayID <- T
#     all_ids <- sort(names(get_terminals(Array_object)))
#   }
#   
#   if(is.null(metadata_list) | is.null(metadata_name)){
#     cat('You need to provide a list of meta data as well as a name for the meta data ... \n meta data has not been updated ... \n')
#     stop()
#   }
#   
#   if(metadata_name %in% get_colnames() & overwrite==F){
#     cat('meta data with the provided name already exists ... if you want to overwrite set overwrite=T \n meta data has not been updated ...\n')
#   }
#   
#   if(!identical(sort(names(metadata_list)),all_ids)){
#     cat('meta data is required for all available proteins; add uniprot ids as names to the metadata_list \n Note: check which proteins are in the array_object using get_available_proteins() and remove duplicates \n meta data has not be added ... \n')
#     stop()
#   }
#   
#   if(overwrite==T){
#     continue <- NA
#     while(continue%notin%c('y','Y','ye','YE','Ye','yes','YES','Yes','YeS','YEs','yeS','yES','n','N','no','No','NO')){
#       continue <- readline('!!!! You are trying to overwrite meta data !!!! \n       Are you sure you want to do this?! \n Note: if you want to update all meta data from uniprot please just use the update_metadata function ... \n (options: "yes" or "no") ')
#     }
#     if(continue%in%c('n','N','no','No','NO')){
#       overwrite <- F
#       stop()
#     }
#     cat('\n')
#   }
#   
#   if(backup==T){
#     cat('Creating Backup ... saved to backup folder\n')
#     array_object <- create_array_backup(array_object)
#   }
#   
#   updated_terminals <- lapply(get_terminals(array_object),function(terminal){
#     
#     ID <- terminal$UniProtKB.AC
#     
#     if(arrayID==T){
#       ID <- paste0(terminal$UniProtKB.AC,'_',terminal$terminal)
#     }
#     
#     terminal[[metadata_name]] <- metadata_list[[ID]]
#     
#     terminal <- move_element_to_bottom_of_list(list = terminal,to_move_to_bottom = 'status')
#     return(terminal)
#     
#   })
#   
#   array_object$Protein_terminals <- updated_terminals
#   
#   updated_array_object <- update_date(array_object)
#   return(updated_array_object)
#   
# }

add_helical_confidence_to_SWaFi_binding_sites <-function(SWaFi_sites,
                                                         cassette_len = 11,
                                                         pH = 7.4,
                                                         array_object=NULL,
                                                         seq_col = 'sequence'){
  
  if(is.null(array_object)){
    array_object <- Find_array_object()
  }
  
  SWaFi_sites <- as.data.frame(SWaFi_sites)
  
  AF <- get_metadata(colnames='AlphaFold', array_object = array_object)
  
  helix_properties <- do.call(rbind,pblapply(1:nrow(SWaFi_sites), function(i){
    subAF <- AF[[SWaFi_sites[i,'UniProtKB.AC']]] %>%
      filter(residue %in% c(SWaFi_sites[i,'start']:SWaFi_sites[i,'end']))
    
    subAF$confidence[!str_detect(subAF$structure,'HELX_RH_AL_P')] <- 0
    
    if(nchar(SWaFi_sites[i,seq_col])<cassette_len){
      hm   <- hydrophobic_moment_screen(SWaFi_sites[i,seq_col], cassette_len = nchar(SWaFi_sites[i,seq_col]))
      h    <- hydrophobicity_screen(SWaFi_sites[i,seq_col],cassette_len = nchar(SWaFi_sites[i,seq_col])) 
      z    <- charge_screen(SWaFi_sites[i,seq_col],cassette_len = nchar(SWaFi_sites[i,seq_col]), pH = pH)
      mers <- SWaFi_sites[i,seq_col]
      conf <- ifelse(nrow(subAF)>0, mean(subAF$confidence),NA)
    } else {
      
      # Extract all 11mers
      mers <- combinated_letters(SWaFi_sites[i,seq_col],n = cassette_len)
      hm   <- hydrophobic_moment_screen(SWaFi_sites[i,seq_col],cassette_len = cassette_len)
      h    <- hydrophobicity_screen(SWaFi_sites[i,seq_col],cassette_len = cassette_len) 
      z    <- charge_screen(SWaFi_sites[i,seq_col],cassette_len = cassette_len, pH = pH)
      
      if(nrow(subAF)>0){
        conf <- rowMeans(rollapply(subAF$confidence,cassette_len,by=1,c))
        MIP  <- hm*(conf+1)
      } else {
        MIP  <- hm
        conf <- NA
      }
      
      # Extract maximum MIP
      mers <- mers[which(MIP==max(MIP))][1]
      hm   <- hm[which(MIP==max(MIP))][1]
      h    <- h[which(MIP==max(MIP))][1]
      z    <- z[which(MIP==max(MIP))][1]
      
      if(nrow(subAF)>0){
        conf   <- conf[which(MIP==max(MIP))][1]
      }
    }
    
    #Avoids error if no AlphaFold data is available:
    names(hm) <- NULL
    names(h)  <- NULL
    names(z)  <- NULL
    
    data.frame(helical_conf=mean(subAF$confidence), 
               helical_residues=length(subAF$structure[str_detect(subAF$structure,'HELX_RH_AL_P')]),
               max_helix=mers,
               max_helix_conf=conf,
               max_helix_H_moment=hm,
               max_helix_hydrophobicity=h,
               max_helix_charge=z)
  }))
  
  if(sum(colnames(helix_properties)%in%colnames(SWaFi_sites))>0){
    continue <- NA
    while(continue%notin%c('y','Y','ye','YE','Ye','yes','YES','Yes','YeS','YEs','yeS','yES','n','N','no','No','NO')){
      continue <- readline('Do you want to overwrite existing helicity data? \n (options: "y" or "n") ')
    }
    if(continue%notin%c('n','N','no','No','NO')){
      SWaFi_sites <- SWaFi_sites[,which(colnames(SWaFi_sites)%notin%colnames(helix_properties))]
    } 
    cat('\n')
  }
  
  cbind(SWaFi_sites,helix_properties)
  
}

add_new_binding_site_confidence <- function(datasets=NULL, 
                                            p_thresh=c(0.01,0.05,0.1), 
                                            confidence=c('high_confidence',
                                                         'medium_confidence',
                                                         'low_confidence'), 
                                            array_object=NULL, 
                                            overwrite=F, ...){
  
  if(is.null(array_object)){
    array_object <- Find_array_object()
  }
  
  item_name <- paste(datasets,collapse = '_')
  
  if(!'Binding_sites'%in%names(array_object)){
    cat('Binding_sites do not exist, creating new list element ... \n')
    array_object$Binding_sites <- list()
    array_object <- move_element_to_bottom_of_list(array_object,c('Import_functions','About'))
  }
  
  if(item_name%in%names(array_object$Binding_sites)){
    if(overwrite==F){
      cat('You are trying to overwrite existing data ... \nPlease specify that you want to do this. Use overwrite=T \n')
      stop()
    } else {
      array_object$Binding_sites <- array_object$Binding_sites[-which(names(array_object$Binding_sites)==item_name)]
    }
  }
  
  scores <- lapply(p_thresh, function(p){
    extract_overlapping_sites(datasets, proteins = species_filter(array_object = array_object), pthresh = p, array_object = array_object, ...)
  })

  names(scores) <- confidence
  
  scores_to_add <- list(scores)
  names(scores_to_add) <- item_name
  
  array_object$Binding_sites <- c(array_object$Binding_sites,scores_to_add)
  
  updated_array_object <- update_date(array_object = array_object)
  
  return(updated_array_object)
}

identify_amphipathic_helices <- function(SWaFi_sites, 
                                         p_thresh=0.05,
                                         helical_thresh=50, 
                                         non_binding_p_thresh=0.95,
                                         H_moment_thresh=NULL, 
                                         helical_conf_column=NULL, 
                                         H_moment_column=NULL, 
                                         pval_column=NULL){
  
  if(is.null(helical_conf_column)){
    helical_conf_column <- 'max_helix_conf'
  }
  
  if(is.null(H_moment_column)){
    H_moment_column <- 'max_helix_H_moment'
  }
  
  if(is.null(pval_column)){
    pval_column <- 'fisher_pval'
  }
  
  if(sum(c(helical_conf_column,H_moment_column,pval_column)%notin%colnames(SWaFi_sites))>0){
    message('Column names not found ....\n- Maybe you have not run the add_helical_confidence_to_SWaFi_binding_sites() function on the SWaFi sites dataframe?\n- If desired columns exist please specify their names ...\n')
    stop()
  }
  
  if(is.null(H_moment_thresh)){
    message('H_moment_thresh not provided, setting H_moment_thresh to the median Hydrophobic moment of non binding sites\n')
    H_moment_thresh <- SWaFi_sites[,H_moment_column]
    H_moment_thresh <- median(H_moment_thresh[SWaFi_sites[,pval_column]>non_binding_p_thresh], na.rm=T) 
  }
  
  SWaFi_sites$amphipathic_helix <- SWaFi_sites[,helical_conf_column] >= helical_thresh & 
                                   SWaFi_sites[,H_moment_column] >= H_moment_thresh 
  
  SWaFi_sites$amphipathic_helix[is.na(SWaFi_sites$amphipathic_helix)] <- T
  
  return(SWaFi_sites)
}

#-############################################################################-#
#-###### Functions used to import COMPARTMENTS database to array_object ######-#
#-############################################################################-#

# download_COMPARTMENTS_database <- function(species='human'){
# 
#   if(is.numeric(species)){
#     available_species <- setNames(c(9606,6239,7227,559292,10116,3702,10090), c('human','worm','fly','yeast','rat','arabidopsis','mouse'))
#     
#     if(species %in% available_species){
#       species <- names(available_species)[which(available_species==species)]
#     } else {
#       cat('Error: species has to be one of:\nhuman (9606), mouse (10090), rat (10116), fly (7227), worm (6239), yeast (559292) or arabidopsis (3702)\n')
#       stop()
#     }
#   }
#   
#   if(species %notin% c('human','mouse','rat','fly','worm','yeast','arabidopsis')){
#     cat('Error: species has to be one of:\nhuman (9606), mouse (10090), rat (10116), fly (7227), worm (6239), yeast (559292) or arabidopsis (3702)\n')
#     stop()
#   }
#   
#   if(here('data','COMPARTMENTS_Database') %notin% list.dirs(here('data'), recursive = F)){
#     dir.create(here('data','COMPARTMENTS_Database'))
#     dir.create(here('data','COMPARTMENTS_Database','downloads'))
#     last_download <- data.frame(date=date(),number_of_downloads=0, species='human')
#     write.csv(last_download,here('data','COMPARTMENTS_Database','last_download.csv'), row.names = F)
#     
#   }
#   
#   # read index file
#   last_download <- read.csv(here('data','COMPARTMENTS_Database','last_download.csv'), header = T)
#   index <- nrow(last_download)
#   
#   cat(paste0('Downloading ',species,' COMPARTMENTS data from https://download.jensenlab.org ... \n'))
#   COMPARTMENTS <- curl::curl_download(url = paste0('https://download.jensenlab.org/',species,'_compartment_integrated_full.tsv'),
#                                       destfile = here('data','COMPARTMENTS_Database','downloads',paste0(species,'_compartment_integrated_full_',index,'.tsv')))
# 
#   COMPARTMENTS <- read_delim(COMPARTMENTS, col_names = F, show_col_types = F)
#   colnames(COMPARTMENTS) <- c('Ensembl_id','Gene','GO_id','GO_term','Confidence')
#   attributes(COMPARTMENTS)$spec <- NULL
#   attributes(COMPARTMENTS)$problems <- NULL
#   COMPARTMENTS <- as.data.frame(COMPARTMENTS)
#   
#   # update index file
#   last_download <- rbind(last_download,data.frame(date=date(),number_of_downloads=index,species=species))
#   write.csv(last_download,here('data','COMPARTMENTS_Database','last_download.csv'), row.names = F)
#   
#   return(COMPARTMENTS)
#     
# }
# 
# update_COMPARTMENTS_database   <- function(array_object=NULL, 
#                                            overwrite=T, backup=T){
#   
#   if(is.null(array_object)){
#     array_object <- Find_array_object()
#   }
#   
#   # Checking species availability
#   species <- as.numeric(unique(unlist(get_metadata(array_object=array_object,colnames='Organism.ID'))))
#   available_species <- setNames(c(9606,6239,7227,559292,10116,3702,10090), c('human','worm','fly','yeast','rat','arabidopsis','mouse'))
#   species <- species[species %in% available_species]
#   
#   if(length(which(species%notin%available_species))>0){
#     cat(paste0('Species: ',paste(species[which(species%notin%available_species)], collapse = ', '),' are not available in the COMPARTMENTS database and will not be added ... \n'))
#   }
#   
#   # Fetching data from database
#   COMPARTMENTS_data <- lapply(species,download_COMPARTMENTS_database)
#   names(COMPARTMENTS_data) <- species
#   
#   cat('Mapping and Extracting Array proteins ... \n')
#   COMPARTMENTS_list <- pblapply(get_terminals(array_object = array_object), function(terminal){
#     
#     protein_COMPARTMENTS  <- COMPARTMENTS_data[[terminal$Organism.ID]][COMPARTMENTS_data[[terminal$Organism.ID]]$Gene %in% unlist(str_split(terminal$Genes,' ')),]
#     rownames(protein_COMPARTMENTS) <- NULL
#     return(protein_COMPARTMENTS)
#   })
#   
#   updated_array_object <- add_new_metadata(array_object = Array_object, metadata_list = COMPARTMENTS_list, metadata_name = 'COMPARTMENTS', overwrite = overwrite, backup = backup)
#   cat('Done!\n')
#   return(updated_array_object)
# }

#-############################################################################-#
#-######### Function used to remove a whole dataset from the object ##########-#
#-############################################################################-#
# remove_experiment <- function(array_object=NULL,dataset){
#   
#   if(is.null(array_object)){
#     array_object <- Find_array_object()
#   }
#   
#   if(dataset %notin% names(array_object$Experiments)){
#     print('Dataset with this identifier does not exist, please enter a valid identifier...')
#   } else {
#     array_object$Protein_terminals <- lapply(array_object$Protein_terminals,function(terminal){
#       terminal <- terminal[names(terminal)!=dataset]
#       return(terminal)
#     })
#     array_object$Experiments <- array_object$Experiments[names(array_object$Experiments)!=dataset]
#     print(paste0('Dataset with identifier: ',dataset,' has been removed!'))
#   }
#   
#   updated_array_object <- update_date(array_object)
#   
#   return(updated_array_object)
# }

#-############################################################################-#
#-############################ DOMAIN FILTERING ##############################-#
#-############################################################################-#

# Author: Junior Agenant
# Function to get the proportion of the site that falls within domains
# Required:
# 1) binding_sites: data frame with the sites and their location
# 2) merizo_domain_regions: df of domains and their location within the protein
# Check for each site if it falls within a domain
proportion_site_in_domain <- function(binding_sites, merizo_domain_regions){
  
  # For loop to get the proportion and longest NDR for each site
  for(site in 1:nrow(binding_sites)) {                                          # loop through the site
    
    if(strip_terminal_specifier(binding_sites$UniProtKB.AC[site]) %in%          # check if the protein is in the set with domain regions
       merizo_domain_regions$UniProtKB.AC) {
      
      range_site <- c(binding_sites$start[site]:binding_sites$end[site])        # get range of the site
      
      domain_regions <-                                                         # select the domain regions for this protein
        merizo_domain_regions[merizo_domain_regions$UniProtKB.AC == 
                                strip_terminal_specifier(binding_sites$UniProtKB.AC[site]),]
      
      range_regions <- unlist(sapply(1:nrow(domain_regions), function(x){       # loop through the domain regions
        domain_regions[x,'domain_start']:domain_regions[x,'domain_end']         # get range of this domain region
      }))
      
      binding_sites[site,'prop_domain'] <- 
        (length(intersect(range_site, range_regions))) / length(range_site)     # get the proportion that falls within domains
      
      
      range_ndr <- setdiff(range_site, range_regions)                           # find the regions of the site that do not fall within domains (NDR)
     
      if(length(rle(diff(range_ndr))$lengths[rle(diff(range_ndr))$values == 1])>0){
        longest_ndr <- max(rle(diff(range_ndr))$lengths[rle(diff(range_ndr))$values == 1]) # find the longest NDR
      } else {
        longest_ndr <- 0
      }
      
      if (longest_ndr > 0) {                                                    # if the longest NDR is larger than 0:
        binding_sites[site,'longest_ndr'] <- longest_ndr                        # add the longest NDR to the longest_ndr column
      } else {
        binding_sites[site,'longest_ndr'] <- 0                                  # otherwise add zero to the longest_ndr column
      }
      
    } else {                                                                    # if the protein of the site is not in the domain regions dataframe:
      binding_sites[site,'prop_domain'] <- NA                                   # add NA for the proportion within domains
      binding_sites[site,'longest_ndr'] <- NA                                   # add NA for the longest NDR
    }
  }
  
  # Add column with status for each site
  binding_sites$in_domain <- rep('untested', nrow(binding_sites))               # created a column for the whether or not within domain status
  
  binding_sites$in_domain[binding_sites$prop_domain == 1] <- 'yes'              # add 'yes' for a proportion of 1
  
  binding_sites$in_domain[binding_sites$prop_domain == 0] <- 'no'               # add 'no' for a proportion of 0
  
  binding_sites$in_domain[binding_sites$prop_domain > 0 &                       # for proportions between 0 and 1 add 'partially'
                            binding_sites$prop_domain < 1] <- 'partially'
  
  # Return the dataframe with the proportions
  return(binding_sites)
}

get_merizo_domains <- function(overlap='within', 
                               pIoU_cutoff=0.75,
                               plDDT_cutoff=70, 
                               proteins=NULL, 
                               array_object=NULL){
  
  merizo <- unique(do.call(rbind, get_metadata(array_object = array_object, 
                                               colnames = 'merizo_domains',
                                               proteins = proteins)))                 # get merizo data from array_object
  
  if(is.numeric(overlap)){
    merizo <- merizo[merizo$overlap<=overlap,]
  } else if (is.character(overlap)) {
    merizo <- merizo[merizo$within==T,]
  }
  
  merizo <- merizo[merizo$pIoU_per_domain >= pIoU_cutoff & 
                   merizo$average_plDDT >= plDDT_cutoff , ]
  
  return(merizo)
}

#-############################################################################-#
#-######## Functions used to manipulate EXPERIMENTS in array_object ##########-#
#-############################################################################-#

# correct_for_multiple_testing <- function(array_object,dataset,pAdjust){
#   regions <- get_SWaFi_regions(array_object,dataset)
#   
#   regions <- lapply(regions,function(x){
#     x$adjusted_pval <- p.adjust(x$pval,method = pAdjust) 
#     return(x)
#   })
#   
#   array_object$Experiments[[dataset]]$SWaFi_regions <- regions
#   array_object$Experiments[[dataset]]$SWaFi_regions_parameters$pAdjust <- pAdjust
#   
#   updated_array_object <- update_date(array_object)
#   return(updated_array_object)  
#   
# }

#-############################################################################-#
#-################# FUNCTIONS FOR CLUSTERING SWaFi MATIRCES ##################-#
#-############################################################################-#

# cluster_SWaFi_matrix <- function(SWaFi_matrix, 
#                                  metadata=NULL, 
#                                  neighbours=12,
#                                  min_dist=0.1, 
#                                  scale=F, 
#                                  center=F, 
#                                  npcs=10,
#                                  set_number_of_PCs=NULL, 
#                                  cluster_resolution=0.5,
#                                  array_object=NULL){
# 
#   library('Seurat')
#   
#   if(is.null(array_object)){
#     array_object <- Find_array_object()
#   }
#   
#   SWaFi_matrix <-  SWaFi_matrix[,which(colSums(is.na(SWaFi_matrix))!=0)]
#   
#   # Create Seurat Object
#   theme_set(new = theme_minimal())
#   
#   Seurat_object <- CreateSeuratObject(counts = t(SWaFi_matrix),
#                                       assay = 'Protein', 
#                                       row.names = rownames(SWaFi_matrix), 
#                                       project = 'HD_Array',
#                                       meta.data = metadata)
#   
#   # Normalise and run PCA
#   #Seurat_object <- NormalizeData(Seurat_object, normalization.method = "LogNormalize", scale.factor = 10000)
#   cat('Scaling data and running PCA ... \n')
#   Seurat_object <- ScaleData(Seurat_object, do.scale = scale, do.center = center)
#   Seurat_object <- RunPCA(Seurat_object, features = rownames(Seurat_object), approx=FALSE, npcs=npcs)
#   
#   if(!is.null(set_number_of_PCs)){
#     PCs <- set_number_of_PCs
#   } else {
#     cat('Estimating elbow point using Kneedle algorithm (Satopaa 2011) ... \n')
#     # Determine percent of variation associated with each PC
#     pct <- Seurat_object[["pca"]]@stdev / sum(Seurat_object[["pca"]]@stdev) * 100
#     
#     # calculate elbow point with Kneedle algorithm (Satopaa 2011)
#     par(mfrow=c(1,1))
#     PCs <- kneedle(pct,-1)
#     rm(pct)
#   }
#   
#   # UMAP and determine clusters 
#   cat('Finding Clusters ... \n')
#   Seurat_object <- FindNeighbors(Seurat_object, graph.name = 'Clusters', dims = 1:PCs)
#   Seurat_object <- FindClusters(Seurat_object,graph.name = 'Clusters', resolution = cluster_resolution, algorithm = 3 , n.start = 100)
#   #Seurat_object <- FindSubCluster(Seurat_object, cluster = '0',graph.name = 'Clusters', resolution = 5, algorithm = 3, subcluster.name = 'sub0')
#   #Seurat_object <- FindSubCluster(Seurat_object, cluster = '1',graph.name = 'Clusters', resolution = .075, algorithm = 3, subcluster.name = 'sub1')
#   cat('Reducing dimensions using UMAP ... \n')
#   Seurat_object <- RunUMAP(Seurat_object, dims=1:PCs, umap.method='umap-learn', n.neighbors=neighbours, min.dist=min_dist)
#   
#   print(DimPlot(Seurat_object, reduction = "umap", pt.size = 3))
#   
#   detach("package:Seurat", unload=TRUE)
#   detach("package:SeuratObject", unload=TRUE)
#   
#   return(Seurat_object)
# }

#-############################################################################-#
#-################### FUNCTIONS for Processing ORTHOLOGS #####################-#
#-############################################################################-#

# read url function using tryCatch adapted from:
# https://stackoverflow.com/questions/12193779/how-to-use-the-trycatch-function
# GET_orthologs <- function(url, .content_type='application/json') {
#   tryCatch(
#     {
#       r <- suppressWarnings(GET(url,content_type(.content_type)))
#       
#       if(is.error(stop_for_status(r))==F){
#         out <- content(r)[[1]][[1]]
#       } else {
#         out <- 'call error'
#       }
#       out
#     },
#     error = function(cond) {
#       message(paste("URL caused an error:", url))
#       message(conditionMessage(cond))
#       # Choose a return value in case of error
#       'no orthologs'
#     },
#     warning = function(cond) {
#       message(paste("URL caused a warning:", url))
#       message(conditionMessage(cond))
#       # Choose a return value in case of warning
#       NULL
#     }
#   )
# }
# 
# map_uniprot_to_ensembl     <- function(UniProtKB.AC=NULL,
#                                        array_object=NULL){
#   
#   # Find array_object in environment
#   if(is.null(array_object)){
#     array_object <- Find_array_object()
#   }
#   
#   if(is.null(UniProtKB.AC)){
#     print(paste0('No UniProt accessions provided - running for example IDs'))
#     UniProtKB.AC <- names(get_terminals(array_object = array_object))[1:10]
#   }
#   
#   ID_df <- data.frame(Array_ID  = UniProtKB.AC, 
#                       UniProtKB.AC = strip_terminal_specifier(UniProtKB.AC))
#   
#   # filter out outdated ids
#   outdated   <- ID_df$UniProtKB.AC[is.outdated(ID_df$UniProtKB.AC)]
#   up_to_date <- ID_df[!is.outdated(ID_df$UniProtKB.AC),]
#   
#   print('GETing data from uniprot servers ... ')
#   info <- retrieve_uniprot_info(unique(up_to_date$UniProtKB.AC))
#   
#   dbrefs <- do.call(rbind, pblapply(info,function(x){
#     
#     org <- x$organism$taxonomy                  # get organism
#     db  <- x$dbReferences$id                    # get mapping data
#     
#     if(org == 9606){
#       db <- db[str_detect(db,'ENSG')]           # find Ensembl identifiers
#       db <- db[!str_detect(db,'ENSGT')]          # remove GeneTree identifiers
#       db <- table(gsub('.*ENSG','ENSG',db))     # count frequency of identifiers
#     }
#     
#     if(org == 7227){
#       db <- db[str_detect(db,'FB:FBgn')]        # find Ensembl identifiers
#       db <- table(gsub('.*FBgn','FBgn',db))     # count frequency of identifiers
#     }
#     
#     if(org == 6239){
#       db <- db[str_detect(db,'WBGene')]         # find Ensembl identifiers
#       db <- table(gsub('.*WBGene','WBGene',db)) # count frequency of identifiers
#     }
#     
#     if(org == 559292){
#       db <- db[str_detect(db,'SGD:S')]          # find Ensembl identifiers
#       db <- table(gsub('SGD:S','S',db))         # count frequency of identifiers
#     }
#     
#     if(org == 10116){
#       db <- db[str_detect(db,'RGD:')]           # find Ensembl identifiers
#       db <- table(gsub('RGD:','',db))           # count frquency of identifiers
#     }
#     
#     ensembl <- names(db)[db==max(db,na.rm=T)]   # select most frequent
#     
#     return(data.frame(UniProtKB.AC=x$accession, Ensembl_ID=ensembl, TaxID=org))
#   }))
#   
#   # merge into ID_df
#   dbrefs <- merge(x=ID_df, y=dbrefs, by='UniProtKB.AC', all.x=T)
#   
#   print('Done - note that the following UniProt accessions are outdated:')
#   print(outdated)
#   
#   return(dbrefs)
# }
# 
# get_ensembl_orthologs <- function(ensembl_mapping = map_uniprot_to_ensembl(),
#                                   server = 'https://rest.ensembl.org',
#                                   .content_type = 'application/json'){
#   
#   # create urls for GETing orthologs
#   ensembl_mapping$url <- paste0(server,'/homology/id/',
#                                 ensembl_mapping$TaxID,'/',
#                                 ensembl_mapping$Ensembl_ID,
#                                 '?content-type=orthologues')
#   
#   ensembl_mapping$url[is.na(ensembl_mapping$Ensembl_ID)] <- NA
#   
#   # Convert to list
#   ensembl_mapping_list <- split(ensembl_mapping,ensembl_mapping$UniProtKB.AC)
#   
#   print(paste0('Retrieving orthologues from ',server))
#   ensembl_orthologs <- pblapply(ensembl_mapping_list, function(x){
#     
#     # GET URLs from server
#     if(!is.na(unique(x$Ensembl_ID))){
#       out <- GET_orthologs(unique(x$url), .content_type = .content_type) 
#       
#     } else {
#       out <- 'outdated'
#     }
#     
#     # Duplicate content if N and C termini are in array_object
#     if(nrow(x)==1){
#       out = list(a=out) 
#     } else {
#       out <- list(a=out,b=out)
#     }
#     names(out) <- x$Array_ID
#     
#     return(out)
#     
#   })
#   
#   ensembl_orthologs <- unlist(ensembl_orthologs, recursive = F, use.names = T)
#   names(ensembl_orthologs) <- gsub('.*\\.','',names(ensembl_orthologs))
#   
#   return(ensembl_orthologs)
#   
# }
# 
# calculate_conservation_score <- function(ensembl_orthologs, 
#                                          ensembl_tree,
#                                          substitution_matrix,
#                                          aminode_matrix,
#                                          array_object=NULL){
#   
#   # Find array_object in environment
#   if(is.null(array_object)){
#     array_object <- Find_array_object()
#   }
#   
#   conservation_scores <- pblapply(names(ensembl_orthologs), function(Array_ID){
#     
#     if(!is.outdated(Array_ID)){
#       
#       orthologs <- ensembl_orthologs[[Array_ID]]
#       
#       if(length(orthologs$homologies)>0){
#         
#         orthologs <- orthologs$homologies
#         
#         # Get names of orthologs
#         names(orthologs) <- unlist(lapply(orthologs, function(x){x$target$species}))
#         
#         # Get aa sequences of orthologs
#         targets <- gsub('-','',
#                         unlist(
#                           lapply(orthologs, function(x){x$target$align_seq})
#                         )
#         )
#         
#         # Test which orthologs are in the ensembl tree
#         targets <- targets[names(targets)%in%tolower(gsub(' ','_',ensembl_tree$tip.label))]
#         
#         # Get source sequence from array_object and add to targets
#         source <- unlist(get_metadata(array_object = Array_object, 
#                                       colnames =  c('Full.Sequence','Organism'),
#                                       proteins = Array_ID), use.names = F)
#         
#         source <- setNames(source[1],
#                            tolower(gsub('_\\(.*','',gsub(' ','_',source[2]))))
#         
#         # Ensure that homo_sapiens does not appear as an ortholog
#         targets <- targets[names(targets)!="homo_sapiens"]
#         
#         sequences <- c(source, targets)
#         sequences <- gsub('U','C',sequences)
#         sequences <- gsub('[^QWERTYIPASDFGHKLCVNM]','',sequences)
#         
#         
#         # Removing multiple orthologs: keeping most related sequence -------------------
#         cat('\nProcessing intra-species paralogs...\n')
#         dups <- unique(names(sequences)[duplicated(names(sequences))])
#         sequences_temp <- sequences[!names(sequences) %in% dups]
#         keep <- c()
#         for(organism in dups){
#           comp <- sequences[names(sequences) == organism]
#           names(comp) <- paste0(names(comp),'-',c(1:length(comp)))
#           
#           bg   <- quiet(msa(c(source,comp), 
#                             type='protein', 
#                             method = 'ClustalOmega', 
#                             substitutionMatrix = 'BLOSUM80',
#                             verbose = F))
#           
#           ref  <- as.character(AAStringSet(bg)[names(source)])
#           mis  <- mismatches(ref, bg, thresh=0.5, remove=F)
#           keep <- c(keep, comp[which(mis==min(mis,na.rm=T))[1]])
#         }
#         if(length(keep)>0){
#           names(keep) <- gsub('-.*','',names(keep))
#         }
#         
#         sequences <- c(sequences_temp,keep)
#         
#         if(length(sequences)>1){
#           
#           cat('Running multiple sequence alignment algorithm ClustalOmega with BLOSUM80...\n')
#           alignment <- quiet(msa(sequences, 
#                                  type = 'protein', 
#                                  substitutionMatrix = 'BLOSUM80',
#                                  order = 'input', method = 'ClustalOmega', verbose = F))
#           
#           # remove sequences with less than 50% alignment --------------------------------
#           
#           # Get human sequence from slignment
#           ref <- as.character(AAStringSet(alignment)[names(source)])
#           
#           cat('Removing sequences with >50% mismatches...\n')
#           mismatched <- mismatches(ref=ref, 
#                                    alignment=alignment,
#                                    thresh = 0.5, 
#                                    remove = T)
#           
#           
#           ref <- unlist(str_split(ref,""))
#           Gaps <- which(ref == '-')
#           
#           autoMasked <- alignment
#           colmask(autoMasked) <- IRanges(start = Gaps, end = Gaps)
#           rowmask(autoMasked) <- IRanges(start = mismatched)
#           
#           #cluster_msa_alignment(autoMasked)
#           
#           if(nrow(as.matrix(autoMasked))>1){
#             
#             conservation <- msaConservationScore(autoMasked, substitution_matrix)
#             conservation <- conservation[!is.na(conservation)]
#             names(conservation) <- unlist(str_split(as.character(AAStringSet(autoMasked)[1]),''))
#             
#             conservation <- data.frame(residues = names(conservation),
#                                        seq = c(1:length(conservation)),
#                                        score = conservation + min(conservation),
#                                        relative.score = (conservation+abs(min(conservation)))/mean(conservation+abs(min(conservation))))
#             
#             conservation$smooth <- smooth(smooth(smooth(conservation$relative.score,span=5),span=3),span=3) #11+7+7 cassettes
#             
#             y <- as.matrix(autoMasked)
#             
#             row.names(y) <- gsub('\\.[0-9]','',
#                                  str_replace(gsub('_',' ',row.names(y)), "^\\w{1}", toupper))
#             
#             temp_tree <- keep.tip(ensembl_tree, row.names(y))
#             edges_ref <- as.data.frame(temp_tree$edge)
#             species   <- data.frame(species=temp_tree$tip.label)
#             
#             cat('Calculating relative weighted substitution score...\n')
#             aminode <- apply(y,2,function(x){
#               
#               if(length(unique(x)) > 1){
#                 
#                 rooted_temp_tree <- multi2di(temp_tree)
#                 
#                 probs <- as.data.frame(ace(phy = rooted_temp_tree,
#                                            x = x[unique(names(x))],
#                                            model = "ER",
#                                            type = "discrete")$lik.anc)
#                 
#                 #fitER <- rerootingMethod(tree=rooted_temp_tree, x = x[unique(names(x))], model='ER')
#                 #plot(fitER)
#                 
#                 rooted_temp_tree$tip.label <- merge(x=species, 
#                                                     y=data.frame(species=names(x),aa=x), sort=F)[,2]
#                 
#                 tips <- data.frame(node=c(1:length(rooted_temp_tree$tip.label)),aa=rooted_temp_tree$tip.label)
#                 
#                 nodes <- rbind(tips,data.frame(node=rownames(probs),
#                                                aa=colnames(probs)[max.col(probs,ties.method="first")]))
#                 
#                 probs <- as.data.frame(sapply(rownames_to_column(probs),as.numeric))
#                 
#                 edges <- merge(x=merge(x=edges_ref,y=nodes,by.x='V2', by.y='node'),
#                                y=nodes,by.x='V1', by.y='node')
#                 
#                 edges$aa <- apply(edges,1,function(j){
#                   
#                   if(!as.numeric(j[2]) %in% tips$node){
#                     aa <- colnames(probs)[max.col(probs[probs$rowname == as.numeric(j[2]),][-1]*
#                                                     aminode_matrix[j[4],colnames(probs)[-1]])+1]
#                   } else { aa <- j[3]}
#                   aa
#                 })
#                 
#                 aminode <- sum(1-apply(edges,1,function(z){aminode_matrix[z[4],z[5]]}))
#                 
#               } else {
#                 aminode <- max(c(temp_tree$edge))*(1-aminode_matrix[unique(x),unique(x)])
#               }
#               aminode
#             })
#             cat('Saving...\n')
#             
#             aminode <- data.frame(aminode=aminode/mean(aminode))
#             aminode$aminode_smooth <- smooth(smooth(smooth(aminode$aminode,span=5),span=3),span=3)
#             aminode$aminode_smooth_spline <- smooth.spline(aminode$aminode,df=30)$y
#             
#             conservation <- cbind(conservation, aminode)
#             
#             colnames(conservation) <- c('AA','residue','Cscore','Cscore.relative','Cscore.smooth',
#                                         'WRSRscore','WRSRscore.smooth','WRSRscore.spline.smooth')
#             
#             
#             terminal <- get_metadata(array_object = array_object, 
#                                      colnames = c('sub.start','sub.end','sub.sequence'),
#                                      proteins = Array_ID)[[1]]
#             
#             terminal_conservation <- conservation[conservation$residue %in% c(terminal$sub.start:terminal$sub.end),]
#             
#             cat('Checking AA Sums...\n')
#             mistake <- F
#             if(paste0(terminal_conservation$AA, collapse = '') == terminal$sub.sequence){
#               cat('Sums checked!\n')
#             } else {
#               mistake <- T
#             }
#             
#             if(mistake==F){
#               out <- list(conservation = conservation,
#                           terminal.conservation = terminal_conservation,
#                           alignment = autoMasked)
#             } else {
#               #cat('Mistake in alignment identified\n')
#               out <- list(conservation = 'mistake',
#                           terminal.conservation = 'mistake',
#                           alignment = autoMasked)
#             }
#             
#           } else {
#             #cat('No orthologs found with > 50% identity!...\n')
#             out <- list(conservation = 'no orthologs with >50% identity',
#                         terminal.conservation = 'mistake',
#                         alignment = autoMasked)
#           }
#         } else {
#           #cat('No orthologs found with > 50% identity!...\n')
#           out <- list(conservation = 'no orthologs with >50% identity',
#                       terminal.conservation = 'no orthologs with >50% identity',
#                       alignment = autoMasked)
#         }
#       } else {
#         #cat('No orthologs found!...\n')
#         out <- list(conservation = 'no orthologs',
#                     terminal.conservation = 'no orthologs',
#                     alignment = 'no orthologs')
#       }
#       
#     } else {
#       #outdated
#       out <- list(conservation = 'outdated',
#                   terminal.conservation = 'outdated',
#                   alignment = 'outdated')
#     }
#     
#   })
#   
#   names(conservation_scores) <- names(ensembl_orthologs)
#   return(conservation_scores)
# }
# 
# cluster_msa_alignment <- function(proteins, method='complete', Plot=T, array_object=NULL, ...){
#   
#   if(is.null(array_object)){
#     array_object <- Find_array_object() 
#   }
#   
#   msa_alignments <- get_metadata(array_object = array_object, 
#                                  colnames = 'alignment', 
#                                  proteins = proteins)
#   
#   if(Plot==T){
#     columns <- floor(sqrt(length(msa_alignments)))
#     rows <- ceiling(length(msa_alignments)/columns)
#     par(mfrow=c(rows,columns))
#   }
#   
#   cluster_data <- lapply(names(msa_alignments),function(protein){
#     
#     name <- unique(get_short_protein_name(protein))
#     
#     msa_alignment <- msa_alignments[[protein]]
#     
#     sdist <- stringDist(as(msa_alignment,'AAStringSet'), method='hamming')
#     clust <- hclust(sdist, ...)
#   
#     if(Plot==T){
#       quiet(print(plot(clust, main = paste0(name,'\n',protein))))
#     }
#     return(list(distance_matrix = sdist, hclust_object = clust))
#   })
#   
#   return(cluster_data)
# }

#-############################################################################-#
#-############################# Hidden functions #############################-#
#-############################################################################-#

# .inflect         <- function(x, threshold = 2){
#   if(length(x)>1){
#     up   <- sapply(1:threshold, function(n) c(x[-(seq(n))], rep(NA, n)))
#     down <- sapply(-1:-threshold, function(n) c(rep(NA,abs(n)), x[-seq(length(x), length(x) - abs(n) + 1)]))
#     a    <- cbind(x,up,down)
#     return(list(minima = which(apply(a, 1, min) == a[,1]), maxima = which(apply(a, 1, max) == a[,1])))
#   } else {
#     return(list(minima=integer(),maxima=integer()))
#   }
# }
# 
# .detect_minima   <- function(SWaFi_trace,lcut,Plot=F){
#   extrema <- .inflect(SWaFi_trace$score,threshold = 2)
#   minmax  <- SWaFi_trace[unlist(extrema),]
#   minmax$id <- c(rep('min',length(extrema$minima)),
#                  rep('max',length(extrema$maxima)))
#   minmax <- minmax[order(minmax$position),]
#   minmax$diff  <- abs(c(diff(minmax$score),1)/minmax$score)
#   minmax$diff2 <- abs(c(1,c(diff(minmax$score))/minmax$score[-nrow(minmax)]))
#   minmax <- minmax[minmax$diff>=lcut|minmax$diff2>=lcut,]
#   minmax <- minmax[!rowSums(is.na(minmax[,1:5])) == 5,]
#   
#   mn <- minmax[minmax$id=='min',]
#   mn <- mn[!is.na(mn$score),]
#   mn$diff[is.infinite(mn$diff)] <- NA
#   mn$diff2[is.infinite(mn$diff2)] <- NA
#   mn <- mn[is.na(mn$diff) + is.na(mn$diff2) < 2,]
#   
#   
#   if(Plot==T){
#     plot(SWaFi_trace$position,SWaFi_trace$score,type=ifelse(nrow(mn)==0,'o','l'),
#          xlab='Position',
#          ylab='Binding score (a.u.)', main= paste0('Local minima cut: ',lcut))
#     points(mn$position,mn$score, col='red', pch = 19)
#     points(minmax$position[minmax$id=='max'],
#            minmax$score[minmax$id=='max'],col='steelblue',pch=19)
#   }
#   return(mn$position)
# }
# 
# .get_signal      <- function(df, 
#                              nthresh=0.2, 
#                              impute=NULL, 
#                              dataset=NA,
#                              impute_distribution=rBCCGo, 
#                              Plot=F){
#   
#   if(is.na(dataset)){
#     print('please provide the name of the dataset...')
#     stop()
#   }
#   
#   nthresh <- nthresh/2
#   
#   if(!is.list(df)){
#     y <- df
#     x <- 1:length(df)
#   } else {
#     y <- df[[dataset]][,'normalised']
#     x <- df[[dataset]][,'Order']
#   }
#   
#   if(length(y[is.na(y)])>0 & !is.null(impute)){
#     
#     if(length(impute)==2){
#       y[is.na(y)] <- impute_distribution(n=length(y[is.na(y)==T]),impute[1],impute[2])
#     }
#     
#     if(length(impute)==3){
#       y[is.na(y)] <- impute_distribution(n=length(y[is.na(y)==T]),impute[1],impute[2],impute[3])
#     }
#     
#     if(length(impute)==4){
#       y[is.na(y)] <- impute_distribution(n=length(y[is.na(y)==T]),mu=impute[1],impute[2],impute[3],impute[4])
#     }
#     
#     if(length(impute)>4 | length(impute)==1){
#       print('model has more or less parameters than accepted... min parameters = 2, max parameters = 4')
#     }
#     
#   }
#   
#   if(is.null(impute)){
#     for(i in which(is.na(y))){
#       
#       if(i == length(y) | i==1){
#         if(i==1){
#           position_after <- NA
#           while(is.na(position_after)){
#             index_after <- 1
#             position_after <- y[i+index_after]
#           }
#           y[i] <- position_after
#         } else {
#           position_before <- NA
#           index_before <- 0
#           while(is.na(position_before)){
#             index_before <- 1
#             position_before <- y[i-index_before]
#           }
#           y[i] <- position_before
#         }
#         
#       } else {
#         
#         position_before <- NA
#         index_before <- 0
#         while(is.na(position_before)){
#           index_before <- 1
#           position_before <- y[i-index_before]
#         }
#         
#         position_after <- NA
#         while(is.na(position_after)){
#           index_after <- 1
#           position_after <- y[i+index_after]
#         }
#         
#         y[i] <- mean(c(position_before,position_after))
#       }
#     }
#   }
#   
#   if(Plot==T){
#     par(mfrow=c(3,3))
#     plot(x,y,xlab='Position of first amino acid in 16mer peptide', 
#          ylab='Fluorescence signal a.u.', main="Normalised Data")
#     lines(x,y)
#   }
#   
#   m <- matrix(0,ncol=length(y)+7,nrow=length(y))
#   
#   for(i in 1:length(y)){
#     m[i,c(i:(i+7))] <- y[i]
#   }
#   
#   mda <- floor(median(1:ncol(m)))
#   mx <- length(x)
#   if(ncol(m)%%2==0){
#     w <- c(1:(mda),(mda):1)
#   }  else {
#     w <- c(1:(mda),(mda-1):1)
#   }
#   
#   if(mx < 8){
#     w[w>mx]<-mx
#   } else {
#     w[w>8] <- 8
#   }
#   
#   m <- colSums(as.data.frame(m))/w
#   x <- c(x,seq(from=max(x)+2,to=max(x)+15,by=2))
#   
#   if(Plot==T){
#     plot(x,m,xlab='Position', ylab='Fluorescence signal a.u.', 
#          main="Positionally Corrected Data")
#     lines(x,m)
#   }
#   
#   #noise removal by fast fourier transformation
#   noisy_signal_ft <- fft(c(rep(0,8),m,rep(0,8)))
#   x_axis <- c(1:length(noisy_signal_ft))
#   
#   if(Plot==T){
#     plot(x = x_axis, y = abs(Re(noisy_signal_ft)), type = "l", col = "blue", 
#          xlab = "Frequency", ylab = "Magnitude", lwd = 2, main="Fast Fourier Transformation")
#   }
#   
#   noisy_signal_ft[round((length(x_axis)*nthresh)):(length(x_axis)-round(length(x_axis)*nthresh))] <- 0 + 0i
#   if(Plot==T){
#     plot(x = x_axis, y = abs(Re(noisy_signal_ft)), type = "l", col = "blue", 
#          xlab = "Frequency", ylab = "Magnitude", lwd = 2, 
#          main=paste0('Keep ',(nthresh+nthresh)*100,'% of Frequencies'))
#   }
#   
#   noisy_signal_ifft <- fft(noisy_signal_ft, inverse = TRUE)/length(noisy_signal_ft)
#   cleaned <- Re(noisy_signal_ifft)[9:(length(x_axis)-8)]
#   
#   if(Plot==T){
#     plot(x=x, y=cleaned, type="l", col="blue",xlab="Position", ylab="Binding score (a.u.)", lwd=3, 
#          main = 'Retransformation')
#   }
#   
#   cleaned[cleaned<0]<-0
#   return(tibble(score=cleaned,position=x))
# }
# 
# .get_signal_fast <- function(y,
#                              nthresh=0.2,
#                              impute=NULL,
#                              impute_distribution=rBCCGo){
#   
#   nthresh <- nthresh/2
#   
#   if(length(y[is.na(y)])>0 & !is.null(impute)){
#     
#     if(length(impute)==2){
#       y[is.na(y)] <- impute_distribution(n=length(y[is.na(y)==T]),impute[1],impute[2])
#     }
#     
#     if(length(impute)==3){
#       y[is.na(y)] <- impute_distribution(n=length(y[is.na(y)==T]),impute[1],impute[2],impute[3])
#     }
#     
#     if(length(impute)==4){
#       y[is.na(y)] <- impute_distribution(n=length(y[is.na(y)==T]),mu=impute[1],impute[2],impute[3],impute[4])
#     }
#     
#     if(length(impute)>4 | length(impute)==1){
#       print('model has more or less parameters than accepted... min parameters = 2, max parameters = 4')
#     }
#     
#   }
#   
#   if(is.null(impute)){
#     for(i in which(is.na(y))){
#       
#       if(i == length(y) | i==1){
#         if(i==1){
#           position_after <- NA
#           while(is.na(position_after)){
#             index_after <- 1
#             position_after <- y[i+index_after]
#           }
#           y[i] <- position_after
#         } else {
#           position_before <- NA
#           index_before <- 0
#           while(is.na(position_before)){
#             index_before <- 1
#             position_before <- y[i-index_before]
#           }
#           y[i] <- position_before
#         }
#         
#       } else {
#       
#       position_before <- NA
#       index_before <- 0
#       while(is.na(position_before)){
#         index_before <- 1
#         position_before <- y[i-index_before]
#       }
#       
#       position_after <- NA
#       while(is.na(position_after)){
#         index_after <- 1
#         position_after <- y[i+index_after]
#       }
#       
#       y[i] <- mean(c(position_before,position_after))
#       }
#     }
#   }
#   
#   m <- matrix(0,ncol=length(y)+7,nrow=length(y))
#   
#   for(i in 1:length(y)){
#     m[i,c(i:(i+7))] <- y[i]
#   }
#   
#   mda <- floor(median(1:ncol(m)))
#   if(ncol(m)%%2==0){
#     w <- c(1:(mda),(mda):1)
#   }  else {
#     w <- c(1:(mda),(mda-1):1)
#   }
#   
#   mx <- length(y)
#   if(mx < 8){
#     w[w>mx]<-mx
#   } else {
#     w[w>8] <- 8
#   }
#   
#   #noise removal by fast fourier transformation
#   noisy_signal_ft <- fft(c(rep(0,8),colSums(as.data.frame(m))/w,rep(0,8)))
#   noisy_signal_ft[round((length(noisy_signal_ft)*nthresh)):(length(noisy_signal_ft)-round(length(noisy_signal_ft)*nthresh))] <- 0 + 0i
#   noisy_signal_ifft <- fft(noisy_signal_ft, inverse = TRUE)/length(noisy_signal_ft)
#   cleaned <- Re(noisy_signal_ifft)[9:(length(noisy_signal_ft)-8)]
#   cleaned[cleaned<0]<-0
#   return(cleaned)
# }
# 
# .get_p_val       <- function(peak, mx, fits, fl, distribution){
#   if(length(peak)>0){
#     sl <- length(peak)
#     ss <- sum((peak+1))/(mx+1)
# 
#     params <- fits[[sl]]$parameters
#     
#     if(length(params)==2){
#       p <- distribution(ss,fits[[sl]][[params[1]]],
#                         fits[[sl]][[params[2]]],
#                         lower.tail=F)
#     }
#     
#     if(length(params)==3){
#       p <- distribution(ss,fits[[sl]][[params[1]]],
#                         fits[[sl]][[params[2]]],
#                         fits[[sl]][[params[3]]], 
#                         lower.tail=F)
#     }
#     
#     if(length(params)==4){
#       p <- distribution(ss,fits[[sl]][[params[1]]],
#                         fits[[sl]][[params[2]]],
#                         fits[[sl]][[params[3]]],
#                         fits[[sl]][[params[4]]],
#                         lower.tail=F)
#     }
#     
#     if(length(params)>4 | length(params)==1){
#       print('model has more parameters than accepted... max parameters = 4')
#     }
#     
#     p <- 1-(1-p)^(1/(length(fits)/(fl/2)))
#   } else {
#     p <- NA
#   }
#   return(p)
# }
# 
# .Define_regions_and_calculate_p_vals <- function(SWaFi_trace,
#                                                  fl=NULL,
#                                                  lcut,
#                                                  fits,
#                                                  mx,
#                                                  distribution,
#                                                  Plot=F){
#   
#   if(is.null(fl)){
#     fl <- nrow(SWaFi_trace)
#   }
#   
#   if(nrow(SWaFi_trace)>1 & sum(SWaFi_trace$score)>0){
#     
#     minimas <- .detect_minima(SWaFi_trace,lcut,Plot)
#     threshold <- c(0,SWaFi_trace$score[order(SWaFi_trace$score)])
#     
#     all_possible_regions <- lapply(1:(nrow(SWaFi_trace)-1),function(x){
#       regions <- SWaFi_trace[SWaFi_trace$score>threshold[x],]
#       
#       aa_intervals <- c(2,diff(regions$position))
#       aa_intervals[regions$position%in%minimas] <- -2
#       
#       regions$ID <- cumsum(aa_intervals != 2)
#       regions <- split(regions,regions$ID)
#       
#       regions <- lapply(regions,function(site){
#         data.frame(start=min(site$position),
#                    end=max(site$position),
#                    pval=.get_p_val(site$score,mx=mx,fits=fits,
#                                   fl=fl,distribution=distribution))
#       })
#       
#       regions <- rbindlist(regions)
#       
#       return(regions)
#     })
#   
#   } else {
#     
#     all_possible_regions <- list(data.frame(start=min(SWaFi_trace$position),
#                                             end=max(SWaFi_trace$position),
#                                             pval=.get_p_val(SWaFi_trace$score,mx=mx,fits=fits,
#                                                             fl=fl,distribution=distribution)))
#     } 
#   
#   return(unique(rbindlist(all_possible_regions)))
# }
# 
# .Approximate_binding_sites <- function(regions, 
#                                        SWaFi_trace, 
#                                        pthresh, 
#                                        lcut, 
#                                        mx,
#                                        fits, 
#                                        distribution, 
#                                        get_min=F, 
#                                        pAdjust,
#                                        useAdjusted=F){
#   
#   if(useAdjusted==T){
#     regions$pval <- regions$adjusted_pval
#     print('Using adjusted p-values ... not recommended')
#   }
#   
#   regions <- regions[,c('start','end','pval')]
#   fl <- nrow(SWaFi_trace)
#   checkpoint <- ifelse(nrow(SWaFi_trace)>1 & sum(SWaFi_trace$score)>0,0,1)
#   
#   if(checkpoint==0){
#     
#     all_ranges <- IRanges(regions$start, regions$end, names = regions$pval)
#     
#     regions <- regions[regions$pval<pthresh,]
#     regions$id <- 1:nrow(regions)
#     potential_sites <- nrow(regions)
#     
#     if(potential_sites>0){
#       results <- list()
#       counter <- 0
#       while(nrow(regions)>0){
#         counter <- counter+1
#         
#         if(get_min==T){
#           
#           smallest_pval <- which(regions$pval==min(regions$pval))
#           site <- regions[smallest_pval,]
#           largest_site  <- which((site$end-site$start)==max(site$end-site$start))
#           site <- site[largest_site,]
#           sele <- site$id[1]
#           site <- site[site$id==sele,]
#           
#         } else {
#           
#           largest_site  <- which((regions$end-regions$start)==max(regions$end-regions$start))
#           site <- regions[largest_site,]
#           smallest_pval <- which(site$pval==min(site$pval))
#           site <- site[smallest_pval,]
#           sele <- site$id[1]
#           site <- site[site$id==sele,]
#           ol <- IRanges::findOverlapPairs(all_ranges,IRanges(site$start, site$end),type='within')
#           p <-  min(as.numeric(as.data.frame(ol@first)$names))
#           site$pval <- p
#           
#         }
#         
#         regions <- regions[-sele,]
#         results[[counter]] <- site #add max bindingsite to list
#         
#         peak <- IRanges(site$start, site$end)
#         regions_ranges <- IRanges(regions$start, regions$end, names=regions$ID)
#         ol <- IRanges::findOverlapPairs(peak, regions_ranges, type="any")
#         ol <- as.data.frame(ol@second)
#         regions <- regions[!regions$ID%in%ol$names,]
#         
#       }
#       
#       results <- unique(rbindlist(results))
#       results <- results[order(results$start),]
#       results <- results[,-4]
#       
#     } else {
#       
#       results <- data.frame(start=NA,end=NA,pval=NA)[-1,]
#       
#     }
#     
#   } else {
#     
#     results <- unique(regions)
#     
#   }
#   
# 
#   # Calculate scores of non-binding regions ------------------------------------
#   minimas <- .detect_minima(SWaFi_trace,lcut)
#   
#   # old site implementation 
#   #minimas <- minimas[!minimas%in%c(results$start,results$end,results$start-1,results$end+1)]
#   
#   #non_bind <- data.frame(start=sort(c(min(SWaFi_trace$position),minimas+1,results$end+1)),
#   #                       end  =sort(c(max(SWaFi_trace$position),minimas,results$start-1)))
#   
#   # new implementation ###################################################
#   peaks <- data.frame(start=sort(c(min(SWaFi_trace$position),minimas)),
#                       end  =sort(c(max(SWaFi_trace$position),minimas)))
#   
#   peaks$ID <- as.character(1:nrow(peaks))
#   
#   peak_ranges <- IRanges(peaks$start, peaks$end, names = peaks$ID)
#   results_ranges <- IRanges(results$start, results$end)
#   ol <- IRanges::findOverlapPairs(results_ranges,peak_ranges,type='within')
#   ol <- as.data.frame(ol@second)
#   non_bind <- peaks[!peaks$ID%in%ol$names,-which(colnames(peaks)=='ID')]
#   non_bind <- non_bind[!non_bind$start>non_bind$end,]
#   
#   #########################################################################
# 
#   non_bind$pval <- apply(non_bind,1,function(x){
#     
#     non_bind_region <- SWaFi_trace[SWaFi_trace$position%in%c(x[1]:x[2]),]
#     
#     if(nrow(non_bind_region)>0){
#       return(unique(
#               min(
#                 p.adjust(
#                   .Define_regions_and_calculate_p_vals(non_bind_region, 
#                                                        mx=mx, 
#                                                        lcut=lcut,fl=fl,
#                                                        Plot=F, 
#                                                        fits=fits, 
#                                                        distribution=distribution)$pval,
#                   method=pAdjust),
#               na.rm = T)
#               )
#              )
#     } else {
#       return(NA)
#     }
#   })
# 
#   if(length(results)==0){
#     results <- non_bind
#   } else {
#     if(ncol(results)==5){
#       results <- rbind(results[,c(2:5)],non_bind)
#     } else {
#       results <- rbind(results,non_bind)
#     }
#   }
#   
#   # Add area under the curve for each region -----------------------------------
#   results$score <- apply(results,1,function(x){
#     sum(SWaFi_trace$score[SWaFi_trace$position%in%c(x[1]:x[2])])
#   })
#   
#   # Add significance logical ---------------------------------------------------
#   if(!is.na(pthresh) & !length(results)==0){
#     results$signif <- results$pval<pthresh
#   } else {
#     if(!length(results)==0){
#       results$signif <- NA
#     }
#   }
#   results <- results[order(results$start),]
#   results <- unique(results)
#   results$ID <- 1:nrow(results)
# 
#   return(results)
#   
# }
# 
# .Calculate_SWaFi_binding_sites <- function(array_object=NULL,
#                                            dataset,
#                                            proteins=NULL,
#                                            pthresh,
#                                            lcut,
#                                            useAdjusted=F,
#                                            get_min=F){
#   
#   if(is.null(array_object)){
#     array_object <- Find_array_object()
#   }
#   
#   fits   <- get_fits(array_object, dataset)
#   distribution <- fits$Function
#   fits <- fits$Fits
#   mx <- get_theoretical_maximum_SWaFi_score(array_object,dataset)
#   
#   if(useAdjusted==T){
#     pAdjust <- get_SWaFi_regions(Array_object,dataset,parameters = T)$SWaFi_parameters$pAdjust
#   } else {
#     pAdjust <- 'none'
#   }
#   
#   #if(is.null(proteins)){
#     regions      <- get_SWaFi_regions(array_object,dataset,proteins)
#     SWaFi_traces <- get_SWaFi_traces(array_object,dataset,proteins)
#     proteins     <- names(SWaFi_traces)
#   # } else {
#   #   regions      <- get_SWaFi_regions(array_object,dataset,proteins)
#   #   SWaFi_traces <- get_SWaFi_traces(array_object, dataset,proteins)
#   #   proteins     <- names(SWaFi_traces)
#   # }
#   
#   bindingsites <- pblapply(proteins, function(id){
# 
#     sites <- .Approximate_binding_sites(regions[[id]],SWaFi_traces[[id]],pthresh,lcut,mx,fits,distribution,get_min,pAdjust,useAdjusted)
#     sites$UniProtKB.AC <- id
#     sites <- sites[,c('UniProtKB.AC','ID','start','end','score','pval','signif')]
#     
#     return(sites)
#   })
#   
#   names(bindingsites) <- proteins
#   
#   return(bindingsites)
#   
# }
# 
# .pivot_binding_site_longer <- function(binding_site){
#   do.call(rbind,apply(binding_site,1,function(x){
#     df <- data.frame(position=c(as.numeric(x[3]):as.numeric(x[4])))
#     df$ID <- as.numeric(x[2])
#     df$pval <- as.numeric(x[6])
#     df$signif <- x[7]
#     df$UniProtKB.AC <- x[1]
#     return(df)
#   }))
# }
# 
# .get_score <- function(protein,
#                        array_object,
#                        cut=NA, 
#                        pthresh=NULL, 
#                        Plot=F,
#                        get_min=F,
#                        lcut=NULL, 
#                        pAdjust=NULL, 
#                        dataset=NULL){
#   
#   # get input to functions -----------------------------------------------------
#   nthresh <- get_SWaFi_traces(array_object,dataset, protein, parameters=T)$SWaFi_parameters$nthresh
#   
#   fits   <- get_fits(array_object, dataset)
#   distribution <- fits$Function
#   fits <- fits$Fits
#   
#   impute <- get_imputation_information(array_object, dataset)
#   impute_distribution <- impute$Function
#   impute <- impute$parameters
#   
#   all_signals <- na.omit(get_dataset(array_object,dataset)$normalised)
#   
#   mx <- get_theoretical_maximum_SWaFi_score(array_object,dataset)
#   
#   # If plotting is TRUE SWaFi trace is re- calculated to show the process ------
#   if(Plot==T){
#     SWaFi_trace <- .get_signal(df = get_terminals(array_object,protein)[[1]], 
#                               nthresh = nthresh, 
#                               impute = impute, 
#                               dataset = dataset, 
#                               impute_distribution = impute_distribution, 
#                               Plot = T)
#   }
#   
#   # If plotting is set to FALSE fetch pre-calculated trace ---------------------
#   if(Plot==F){
#     SWaFi_trace <- get_SWaFi_traces(array_object,dataset,proteins=protein)[[1]]
#   }
#   
#   fl <- nrow(SWaFi_trace) # no of peptides sampled
#   
#   # If hard cut-off is set, everything under the cut-off is set to 0 -----------
#   if(!is.na(cut)){
#     SWaFi_trace[SWaFi_trace<=cut] <- 0
#   }
#   
#   # Calculate regions and deterimine significance ------------------------------
#   res <- .Define_regions_and_calculate_p_vals(SWaFi_trace,fl,lcut,fits,mx,distribution,Plot)
#   
#   if(Plot==T){ 
# 
#     clr <- rev(viridis::viridis(nrow(res), alpha = 0.75))
#     res$log <- -log10(res$pval)
#     plot(c(min(res$start),max(res$end)), rep(-log10(pthresh),2),
#          ylim=c(min(res$log),max(res$log)),
#          xlab='Position',ylab='-log10(p-value)', type = "l",lty=2,col='red',lwd=1.5)
#     title(main = paste0('Cut-off alpha ',pthresh))
#     for(i in 1:nrow(res)){
#       lines(c(res$start[i],res$end[i]), rep(res$log[i],2), col=clr[i])
#     }
#     res <- res[,-c('log')]
#   }
#   
#   # Correct for multiple testing -----------------------------------------------
#   res$adjusted_pval <- p.adjust(res$pval,method=pAdjust)
#   
#   # Define binding sites -------------------------------------------------------
#   res <- .Approximate_binding_sites(res,SWaFi_trace,pthresh,lcut,mx,fits,distribution,get_min,pAdjust)
#   
#   # Add identifiers and arrange final table ------------------------------------
#   res$UniProtKB.AC <- protein
#   res <- res[,c('UniProtKB.AC','ID','start','end','score','pval','signif')]
#   
#   # Plot results ---------------------------------------------------------------
#   if(Plot==T){
#     plot(SWaFi_trace$position,SWaFi_trace$score, xlab='Position', ylab='Binding score (a.u.)',
#          type = "l",col='grey80', ylim=c(-0.01,max(SWaFi_trace$score)+0.01),lty=2)
#     abline(h=cut, col='red', lwd=2,lty=2)
#     axis(3, at=unique(c(res$start[res$signif==T],
#                         res$end[res$signif==T])),
#          labels=as.character(unique(c(res$start[res$signif==T],
#                                       res$end[res$signif==T]))),las=2)
#     
#     no_of_sites <- nrow(res[res$signif==T,])
#     
#     if(no_of_sites>0){
#       clr <- viridis::viridis(no_of_sites)
#       for(j in 1:no_of_sites){
#         i <- which(res$signif==T)[j]
#         site <- SWaFi_trace[SWaFi_trace$position >= res$start[i] & SWaFi_trace$position <= res$end[i],]
#         abline(v=c(res$start[i],res$end[i]), col=clr[j], lty=2, lwd=2)
#         points(site$position,site$score, type="o", pch=19, col=clr[j])
#         text(median(site$position),max(site$score)+1,
#              labels = paste0('p-value\n~',round(res$pval[i],4)),cex=0.75, col=clr[j])
#       }
#     }
#   } 
#   
#   # Plot table of scores -------------------------------------------------------
#   if(Plot==T){
#     pres <- res
#     pres$pval  <- round(pres$pval,4)
#     pres$score <- round(pres$score,2)
#     
#     plot.new()
#     addtable2plot(0,0,pres, 
#                   xpad=.2, ypad=.5,
#                   bty='o',
#                   display.rownames = F, 
#                   hlines = TRUE,
#                   vlines = TRUE)
#   }
#   
#   # Return results -------------------------------------------------------------
#   return(res)
#   
# }

#-############################################################################-#
#-########################### GENERIC FUNCTIONS ##############################-#
#-############################################################################-#

# Kneedle Algorithm
# Find the knee/elbow point in a vector using the Kneedle algorithm.
#
# This function uses the Kneedle algorithm (Satopaa 2011)
# to find the index of the knee point in the provided vector.
# If the values are mostly increasing, use sign = 1. If they are
# mostly decreasing, use sign = -1.
#
# @param values The values to find a knee/elbow in.
# @param sign -1 if the values are mostly decreasing, 1 if the 
# values are mostly decreasing.
#
# @return The index of the knee/elbow.
#
# edited from MollahLab

dist2d  <- function(a,b,c) {
  v1 <- b - c
  v2 <- a - b
  m <- cbind(v1,v2)
  d <- det(m)/sqrt(sum(v1*v1))
  d
}
kneedle <- function(values, sign) {
  start = c(1, values[1])
  end = c(length(values), values[length(values)])
  k <- which.max(lapply(1:length(values),
                        function(idx) {
                          sign * -1 * dist2d(c(idx, values[idx]),
                                             start,
                                             end
                          )
                        })
  )
  plot(1:length(values),values,type="l", xlab='PC',ylab='Standard deviation')
  lines(c(1,length(values)),c(values[1],values[length(values)]),col="red")
  points(k, values[k],col="blue")
  abline(v=k,col="blue")
  return(k)
}

# Filter functions
SubsetList <- function(df, filter){
   df[names(df) %in% filter]
  }

# custom not in function
`%notin%` <- Negate(`%in%`)

# Get column medians function
colMedians <- function(matrix, na.rm=T){apply(matrix, 2, median,na.rm=na.rm)}

# Get row medians function
rowMedians <- function(matrix, na.rm=T){apply(matrix, 1, median,na.rm=na.rm)}

# Get column max function
colMax     <- function(matrix, na.rm=T){apply(matrix, 2, max,na.rm=na.rm)}

# Get row max function
rowMax     <- function(matrix, na.rm=T){apply(matrix, 1, max,na.rm=na.rm)}

# Split string into overlapping sequences of length n
combinated_letters <- function(string, n) {
  length_ <- str_length(string)
  str_sub(string, seq(1, length_ + 1 - n), seq(n, length_))
}

# Smooth function
smooth <- function(col, span=1, treat_start=T){
  s <- c()
  for (f in 1:length(col)){
    v <- f-span
    d <- f+span
    if(v < 0){
      v <- f-f
      d <- f+f
    }
    if(d > length(col)){
      d <- f+abs(f-length(col))
    }
    s <- c(s,mean(col[v:d]))
  }
  return(s)
}

# interlacing two vectors
riffle <- function (a,b) { 
  n  <- min(length(a),length(b)) 
  p1 <- as.vector(rbind(a[1:n],b[1:n])) 
  p2 <- c(a[-(1:n)],b[-(1:n)]) 
  c(p1,p2) 
} 

# Calculate hydrophobic moment and individual vectors
Hydrophobic_moment <- function(sequence_list, 
                               get_vectors=F, 
                               pH=7.4){
  
  AA_Properties <- AA_Properties(pH=pH)
  
  Hm_list <- lapply(sequence_list, function(sequence){
    seq <- unlist(str_split(sequence,''))
    seq <- data.frame(AA=seq,order=1:length(seq))
    seq <- merge(seq,AA_Properties,by='AA',all.x=T)
    seq <- seq[order(seq$order),]
    rownames(seq) <- NULL
    
    # Calculate location of residues assuming a perfect alpha helix
    seq$cos <- cos(seq$order*100*pi/180)
    seq$sin <- sin(seq$order*100*pi/180)
    seq$start <- 0 #start of vectors
    
    # Calculate hydrophobicity vector sums
    res <- sqrt(sum(seq$Hydrophobicity*seq$cos)^2+sum(seq$Hydrophobicity*seq$sin)^2)/nrow(seq)
    
    if(get_vectors==T){
      res <- seq
    }
    
    return(res)
  })
  
 return(Hm_list)
}

# Screen sequences for Hydrophobicity with a cassette length of n
hydrophobicity_screen <- function(sequences,
                                  cassette_len=11){
  
  H_matrix <- do.call(rbind,lapply(split(c(0.310,-1.010,-0.600,-0.770,1.540,1.540,-0.220,-0.640,0.000,0.130,1.800,1.700,-0.990,1.230,1.790,0.720,-0.040,0.260,2.250,0.960,1.220)[match(unlist(str_split(sequences,'')),
                                         c("A", "R", "N", "D", "C", "U" ,"Q", "E","G", "H", "I", "L", "K", "M", "F","P", "S", "T", "W", "Y", "V"))],rep(1:length(sequences),nchar(sequences))),
                                   function(x){rollapply(x,cassette_len,by=1,c)}))
  
  unlist(split(rowMeans(H_matrix),rep(1:length(sequences),nchar(sequences)-cassette_len+1)))
}

# Screen sequences for Charge with a cassette length of n
charge_screen <- function(sequences,
                          cassette_len=11,
                          pH=7.4){
  
  H_matrix <- do.call(rbind,lapply(split(c(0/(10^(14-pH)+1), 1/(10^(pH-12.48)+1),
                                           0,-1/(10^(3.65-pH)+1),-1/(10^(8.18-pH)+1),
                                           -1/(10^(8.18-pH)+1),
                                           0,-1/(10^(4.25-pH)+1),0,1/(10^(pH-6)+1),
                                           0,0,1/(10^(pH-10.53)+1),-1/(10^(14-pH)+1),
                                           0,0,-1/(10^(14-pH)+1),-1/(10^(14-pH)+1),
                                           0,-1/(10^(10.07-pH)+1),0)[match(unlist(str_split(sequences,'')),
                                                                           c("A", "R", "N", "D", "C", "U" ,"Q", "E","G", "H", "I", "L", "K", "M", "F","P", "S", "T", "W", "Y", "V"))],rep(1:length(sequences),nchar(sequences))),
                                   function(x){rollapply(x,cassette_len,by=1,c)}))
  
  unlist(split(rowSums(H_matrix),rep(1:length(sequences),nchar(sequences)-cassette_len+1)))
}

# Get distribution of charge from an aligned set of sequences
charge_distribution <- function(sequences,
                                cassette_len=11,
                                pH=7.4, 
                                summary_function=colMeans){
  
  if(length(unique(nchar(sequences)))>1){
    print('Sequences are not of the same length, add "-" to sequences in alignment ... ')
    stop()
  }
  
  H_matrix <- do.call(rbind,lapply(split(c(0/(10^(14-pH)+1), 1/(10^(pH-12.48)+1),
                                           0,-1/(10^(3.65-pH)+1),-1/(10^(8.18-pH)+1),
                                           -1/(10^(8.18-pH)+1),
                                           0,-1/(10^(4.25-pH)+1),0,1/(10^(pH-6)+1),
                                           0,0,1/(10^(pH-10.53)+1),-1/(10^(14-pH)+1),
                                           0,0,-1/(10^(14-pH)+1),-1/(10^(14-pH)+1),
                                           0,-1/(10^(10.07-pH)+1),0,NA)[match(unlist(str_split(sequences,'')),
                                                                           c("A", "R", "N", "D", "C", "U" ,"Q", "E","G", "H", "I", "L", "K", "M", "F","P", "S", "T", "W", "Y", "V","-"))],rep(1:length(sequences),nchar(sequences))),
                                   function(x){rollapply(x,cassette_len,by=1,c)}))
  
  if(!is.null(summary_function)){
    H_matrix <- summary_function(H_matrix, na.rm=T)
  }
  return(H_matrix)
}

# Screen sequences for Hydrophobic moment with a cassette length of n
hydrophobic_moment_screen <- function(sequences,
                                      cassette_len=11){
  
  H_matrix <- do.call(rbind,lapply(split(c(0.310,-1.010,-0.600,-0.770,1.540,1.540,-0.220,-0.640,0.000,0.130,1.800,1.700,-0.990,1.230,1.790,0.720,-0.040,0.260,2.250,0.960,1.220)[match(unlist(str_split(sequences,'')),
                                                                                                                                                                                       c("A", "R", "N", "D", "C", "U" ,"Q", "E","G", "H", "I", "L", "K", "M", "F","P", "S", "T", "W", "Y", "V"))],rep(1:length(sequences),nchar(sequences))),
                                   function(x){rollapply(x,cassette_len,by=1,c)}))
  
  unlist(split(sqrt(rowSums(H_matrix*matrix(rep(cos(c(1:cassette_len)*100*pi/180),nrow(H_matrix)),ncol=cassette_len,byrow = T))^2 + 
                      rowSums(H_matrix*matrix(rep(sin(c(1:cassette_len)*100*pi/180),nrow(H_matrix)),ncol=cassette_len,byrow = T))^2)/cassette_len,rep(1:length(sequences),nchar(sequences)-cassette_len+1)))
}

# Test if there is a region that could be an AH
# This function tests if a region with Hydrophobic moment 
# has occurred by chance or if it is a continuous region
# Idea by Tommas T. E. Nielsen
# Does not work yet!
is_amphipathic_helix_region <- function(sequence,
                                        min_len=11,
                                        thresh=0.1){
  
  Hm <- hydrophobic_moment_screen(sequence, cassette_len = min_len)
  m <- which(Hm==max(Hm))
  
  if(m==length(Hm)){
    after <- c(0,diff(Hm))[m]
  } else {
    after <- c(0,diff(Hm))[m+1]
  }
  
  if(m==1){
    before <- c(0,diff(Hm))[m+1]
  } else {
    before <- c(0,diff(Hm))[m]
  }
  
  AH <- ifelse(before<thresh | after<thresh,T,F)
  return(AH)
}

# Get cassette of length n with the highest  hydrophobic moment
max_hydrophobic_moment <- function(sequences,
                                   cassette_len=11,
                                   pH=7.4){
  sequences <- na.omit(sequences)
  rbindlist(pblapply(sequences,function(x){
    
    mers <- combinated_letters(x,cassette_len)
    Hm <- unlist(Hydrophobic_moment(mers,pH=pH))
    
    Hm <- data.frame(sequence=x,XXmers=mers,max.H.moment=Hm)
    Hm <- Hm[Hm$max.H.moment==max(Hm$max.H.moment),]
  }))
  
}

# Fast function to calculate the maximum hydrophobic moment for many proteins
Screen_for_max_hydrophobic_moment <- function(sequences,
                                              cassette_len=11){

  H_matrix <- do.call(rbind,lapply(split(c(0.310,-1.010,-0.600,-0.770,1.540,1.540,-0.220,-0.640,0.000,0.130,1.800,1.700,-0.990,1.230,1.790,0.720,-0.040,0.260,2.250,0.960,1.220)[match(unlist(str_split(sequences,'')),
                                         c("A", "R", "N", "D", "C", "U" ,"Q", "E","G", "H", "I", "L", "K", "M", "F","P", "S", "T", "W", "Y", "V"))],rep(1:length(sequences),nchar(sequences))),
                                   function(x){rollapply(x,cassette_len,by=1,c)}))
  
  unlist(lapply(split(sqrt(rowSums(H_matrix*matrix(rep(cos(c(1:cassette_len)*100*pi/180),nrow(H_matrix)),ncol=cassette_len,byrow = T))^2 + 
                             rowSums(H_matrix*matrix(rep(sin(c(1:cassette_len)*100*pi/180),nrow(H_matrix)),ncol=cassette_len,byrow = T))^2)/cassette_len,rep(1:length(sequences),nchar(sequences)-cassette_len+1)),max))
}

# Heliquest function
heliquest <- function(sequence,
                      plot=T,
                      pH=7.4){
  
  #Reading in aa hydrophobicity data
  AA_Properties <- AA_Properties(pH=pH)
  
  # Create data frame of sequence and merge with properties
  seq <- unlist(str_split(sequence,''))
  seq <- seq[seq %in% AA_Properties$AA]
  seq <- data.frame(AA=seq,order=1:length(seq))
  seq <- merge(seq,AA_Properties,by='AA',all.x=T)
  seq <- seq[order(seq$order),]
  rownames(seq) <- NULL
  
  # Calculate location of residues assuming a perfect alpha helix
  seq$cos <- cos(seq$order*100*pi/180)
  seq$sin <- -sin(seq$order*100*pi/180)
  seq$start <- 0 #start of vectors
  
  # Calculate hydrophobicity vector sums
  original_vec <- c(sum(seq$Hydrophobicity*seq$cos)/nrow(seq),sum(seq$Hydrophobicity*seq$sin)/nrow(seq))
  
  # Create rotation matrix to point hydrophobic moment vector down
  pivot <- (acos(abs(original_vec[1])/sqrt(abs(original_vec[1])^2+abs(original_vec[2])^2))*180/pi)
  if(original_vec[1]>0 & original_vec[2]>0){pivot <- 270-pivot}
  if(original_vec[1]<0 & original_vec[2]>0){pivot <- 90+pivot}
  if(original_vec[1]<0 & original_vec[2]<0){pivot <- 270-pivot-180}
  if(original_vec[1]>0 & original_vec[2]<0){pivot <- 270+pivot}
  if(original_vec[1]==0 & original_vec[2]==0){pivot <- 0}
  
  rot_matrix <- matrix(c(cos(pivot*pi/180),sin(pivot*pi/180),-sin(pivot*pi/180),cos(pivot*pi/180)),ncol = 2, byrow = F)
  
  # Rotate Hm vector and residue positions
  rotated_vec <-    rot_matrix %*% original_vec 
  orginal_positions <- seq[,c('cos','sin')]
  seq$cos <- apply(orginal_positions,1,function(x){(rot_matrix %*% c(x[1],x[2]))[1,1]})
  seq$sin <- apply(orginal_positions,1,function(x){(rot_matrix %*% c(x[1],x[2]))[2,1]})
  
  fctr <- (0.375*(sort(rep((1:ceiling(nrow(seq)/18)),18))[1:nrow(seq)]-1))+1
  
  # calculate mean hydrophobic moment
  moment <- sqrt(original_vec[1]^2 + original_vec[2]^2)
  
  # Plot helical wheel
  p <- ggplot(seq, aes(cos*fctr,sin*fctr)) + 
    geom_path(aes(col=order, linewidth=order)) + 
    scale_linewidth(range=c(0.5,3)) +
    scale_color_gradient(low='grey80',high='grey30') + 
    ggnewscale::new_scale_color() +
    geom_segment(aes(x=start, y=start, xend=cos*Hydrophobicity ,yend=sin*Hydrophobicity),
                 arrow = arrow(length = unit(0.5, "cm"),type = 'open'),linetype=2,linewidth=.5)  + 
    annotate('segment',x=0, y=0, xend=rotated_vec[1,1] ,yend=rotated_vec[2,1],
             arrow = arrow(length = unit(0.5, "cm"),type = 'closed'), linewidth=1.5) + 
    geom_point(aes(col=AA),pch=16,size=20) +
    geom_text(aes(label=AA)) +
    geom_text2(aes(label=c('N',rep('',nrow(seq)-2),'C')), nudge_x = -0.16, nudge_y = -0.1,col='red') +
    prism() + remove_x_axis() + remove_y_axis() + theme(legend.position='none') +
    xlim(-2.3,2.3) +
    ylim(-2.3,2.3) +
    scale_color_manual(values = setNames(as.character(AA_Properties$color),AA_Properties$AA)) + 
    annotate(geom="text", x=0, y=0.2, label=paste0('<µH> ',round(moment, digits = 3)), size = unit(8, "pt"))
  
  if(plot==T){
   print(p)
  }
  
  
  return(list(coordinates=seq,H.moment=moment, plot=p))

}

# Sped up version of heliquest
fast_heliquest_coords <- function(sequence, pH=7.4){
  
  # Create data frame of sequence and merge with properties
  seq <- plyr::join(data.frame(AA=unlist(str_split(sequence,'')),order=1:nchar(sequence)), AA_Properties(pH=pH),by='AA',match='first')

  # Calculate hydrophobicity vector sums
  original_vec <- c(sum(seq$Hydrophobicity*cos(seq$order*100*pi/180))/nrow(seq),sum(seq$Hydrophobicity*(-sin(seq$order*100*pi/180)))/nrow(seq))
  
  # Create rotation matrix to point hydrophobic moment vector down
  if(original_vec[1]>0 & original_vec[2]>0){pivot <- 270-(acos(abs(original_vec[1])/sqrt(abs(original_vec[1])^2+abs(original_vec[2])^2))*180/pi)}
  if(original_vec[1]<0 & original_vec[2]>0){pivot <- 90 +(acos(abs(original_vec[1])/sqrt(abs(original_vec[1])^2+abs(original_vec[2])^2))*180/pi)}
  if(original_vec[1]<0 & original_vec[2]<0){pivot <- 270-(acos(abs(original_vec[1])/sqrt(abs(original_vec[1])^2+abs(original_vec[2])^2))*180/pi)-180}
  if(original_vec[1]>0 & original_vec[2]<0){pivot <- 270+(acos(abs(original_vec[1])/sqrt(abs(original_vec[1])^2+abs(original_vec[2])^2))*180/pi)}
  if(original_vec[1]==0 & original_vec[2]==0){pivot <- 0}
  
  # Rotate Hm vector and residue positions
  seq$cos <- apply(data.frame(cos=cos(seq$order*100*pi/180), sin= -sin(seq$order*100*pi/180)),1,function(x){(matrix(c(cos(pivot*pi/180),sin(pivot*pi/180),-sin(pivot*pi/180),cos(pivot*pi/180)),ncol = 2, byrow = F) %*% c(x[1],x[2]))[1,1]})
  seq$sin <- apply(data.frame(cos=cos(seq$order*100*pi/180), sin= -sin(seq$order*100*pi/180)),1,function(x){(matrix(c(cos(pivot*pi/180),sin(pivot*pi/180),-sin(pivot*pi/180),cos(pivot*pi/180)),ncol = 2, byrow = F) %*% c(x[1],x[2]))[2,1]})
  seq$start <- 0 #start of vectors
  
  return(seq)
  
}

# Identify the optimal amphipathic helix length

# This function tries to optimize the AH based on binding site length as well as 
# size of the hydrophobic face
# The charge and polar weights tell the function how hard K, R, S, T, N and Q
# should be punished for obstructing the hydrophobic face.
# The misplaced_hydrophobic_weigth tells the function how hard to punish 
# hydrophobic residues in the hydrophilic face
# Proline_weight tells the function how har to punish prolines in the helix 
# if allow_prolines is set to TRUE.
# len_requirement determines when length of the helix should stop being 
# punished with respect to the input sequence length 
# (if it has value 1/3 (default), all potential helices short than 1/3 times 
# the input sequence length will be punished, 
# this ensures that a short helix of e.g. 11aa will not be favored to explain 
# the binding of a sequence longer than 33aa's).
# The punishable length can also be fixed if the above is not desired 
# using fixed_len_threshold. 
# This option takes a length in aa e.g. 11 as input. This would only punish 
# sequences shorter than 11aa.
# The log_base parameter is used to define how to weigh the punishment of 
# short sequences.
# They are punished on a log scale so long binding sites do not weight too much, 
# fx. log10 will favor shorter sequences while log2 will favor longer sequences.
# Choose the log_base depending on how much weight you want to give longer 
# sequences default is ln.
# If you have long sequences I suggest to limit max_len to a reasonable 
# length, fx. 33aa, else run will take longer
# The segregation cutoff is based on the the distribution of amino acids with
# hydrophobic, polar and charged nature in the helix;
# the default 2 in the perfect helix would allow for a face of minimum 
# 2 tryptophans, due to its high hydrophobicity, 
# but usually this ensures a minimum face of 3 hydrophobic residues.
# Helicity is calculated either using netsurf3.0 (needs internet connection) or
# defaults to AlphaFold if an array_object is provided (local).
# Should only be run on pre-defined regions known to bind membranes; 
# this function alone does NOT predict membrane binding

get_best_AH <- function(sequences, 
                        UniProtKB.AC = NULL, 
                        min_len = 7, 
                        speed_up = T,
                        speed_up_by = 3,
                        res = 0.2,
                        max_len = 33,
                        len_requirement = 1/3, 
                        fixed_len_threshold = NULL,
                        use_netsurf = F, 
                        plot = F, 
                        G_to_hydrophobic = F, 
                        A_to_hydrophobic = F, 
                        charge_weight = 1, 
                        polar_weight = 1, 
                        misplaced_hydrophobic_weight = 1, 
                        allow_prolines = T, 
                        proline_weight = 1, 
                        log_base = exp(1), 
                        return_all_info = T, 
                        return_top_hit = T, 
                        set_min_moment = NULL, 
                        set_segregation_cutoff = 2, 
                        array_object = NULL){
  
  # Statement to decide whether to initiate netsurf or not ---------------------
  if(!is.null(UniProtKB.AC)){
    if(use_netsurf==T){
      cat('use_netsurf set to TRUE, but will not be used because UniProtKB.AC has been provided\n')
      use_netsurf <- F
    }
    
    # if netsurf is not used an array object has to be loaded in the env. with AlphaFold or other structure predictions
    if(is.null(array_object)){
      array_object <- Find_array_object() 
    }
  }
  
  # Notify user if the Speed_up has been used, might impact precision ----------
  if(speed_up==T){
    cat(paste0('Speed_up set to TRUE only keeping sequences with the top ',10*(1/speed_up_by),'% highest hydrophobic moment for each length ... \nMight be less accurate, but can increase speed by up to 50%! ...\n'))
  }
  
  # Begin function -------------------------------------------------------------
  
  result <- lapply(seq_along(sequences), function(i){  # Running lapply to vectorize function 
  
    sequence <- sequences[i] # get sequence
    
    if(!is.na(sequence)){
    
      if(!is.null(UniProtKB.AC)){
          UniProt <- UniProtKB.AC[i] # get UniProt ID if array_object is used for helicity
      }
      
      # Setting maximum length of AHs to scan for --------------------------------
      if(is.null(max_len)){
        max_len <- nchar(sequence)
      }
      
      if(max_len>nchar(sequence)){
        max_len <- nchar(sequence)
      }
      
      # Get all possible sequences of lengths from min_len to max_len ------------
      Hm <- na.omit(rbindlist(lapply(min_len:max_len,function(l){
        
        mers <- combinated_letters(sequence,l)
        
        # remove sequences with prolines if requested
        if(allow_prolines==F){
          mers <- mers[!str_detect(mers,'P')]
        }
        
        if(length(mers)>0){
          data.frame(XXmers=mers,H.moment=unlist(Hydrophobic_moment(mers))) # return seqs and their Hm
        } else {
          data.frame(XXmers=NA,H.moment=NA) # no seqs found without prolines; NA returned
        }
      })))
      
      if(allow_prolines==F){
        allow_prolines_temp <- F
      } else {
        allow_prolines_temp <- T
      }
      
      # If allow_prolines=FALSE and no peptides without prolines have been found, prolines are included for that sequence
      if(nrow(Hm)<1){
        cat(paste0('No peptides of length ',min_len,' without prolines in sequence: ', sequence ,' ... \n prolines will be kept in sequence \n'))
        allow_prolines_temp <- T
        Hm <- rbindlist(lapply(min_len:max_len,function(l){
          data.frame(XXmers=combinated_letters(sequence,l),H.moment=unlist(Hydrophobic_moment(mers)))
        }))
      }
      
      # Speed_up will discard sequences with low hydrophobic moment --------------
      if(speed_up==T){
        Hm <- do.call(rbind,lapply(split(Hm,nchar(Hm$XXmers)), function(x){
          n <- 1+round(nrow(x)*(1/speed_up_by)) # keeps top 30% of seqs for each length based on Hm
          if(n<3){n <- nrow(x)} # Ensures that at least 5 seqs are retained
          x$keep <- x$H.moment>=(sort(x$H.moment, decreasing = T)[n])
          return(x)
        }))
      
        keep <- Hm$keep
        Hm <- Hm[Hm$keep==T,1:2]
      } else {
        keep <- rep(T,nrow(Hm)) # else all sequences are kept
      }
      
      # Calculate position of amino acids with respect to the mean hydrophobic moment for each sequence, 
      # calculate size of hydrophobic face and punish amino acids such as charges in the hydrophobic face.
      Hm$segregation <- sapply(Hm$XXmers,function(mer){
        coords <- fast_heliquest_coords(mer) # get amino acid coordinates with respect to Hm vector
        
        # Group amino acids into categories
        coords$categories <- as.character(coords$categories)
        if(G_to_hydrophobic==T){coords$categories[coords$AA=='G'] <- 'Hydrophobic'} 
        if(A_to_hydrophobic==F){coords$categories[coords$AA=='A'] <- 'Special Cases'}
        coords$categories[coords$categories %notin% c('Hydrophobic','Special Cases','Polar')] <- 'Charged'
        
        # Get largest hydrophobic face, punish based on location with respect to Hm
        face <- sapply(str_extract_all(sapply(seq(-1,1,res),function(i){
          hyd_face <- paste(substring(coords$categories, 1, 1)[coords$sin<=i],collapse='')}), "(H)\\1+"), function(x){
            ifelse(length(x)>0,max(nchar(x)),0)})*(1-seq(0,1,1/(2/res)))
  
        # Punish location of charged residues
        Charge <- sapply(seq(-1,1,res),function(i){
          sum(abs(coords$Hydrophobicity[coords$sin<=i & coords$categories== 'Charged']))
        })*(charge_weight-seq(0,charge_weight,1/((2/res)/charge_weight))) # weight of punishment
        
        # Punish location of polar residues
        Polar <- sapply(seq(-1,1,res),function(i){
          sum(abs(coords$Hydrophobicity[coords$sin<=i & coords$categories=='Polar']))
        })*(polar_weight-seq(0,polar_weight,1/((2/res)/polar_weight))) # weight of punishment
        
        # Punish location of hydrophobic residues in hydrophilic face
        hydrophobic_misplaced <- sapply(seq(-1,1,res),function(i){
          sum(coords$Hydrophobicity[coords$sin>=i & coords$categories=='Hydrophobic'])
        })*seq(0,misplaced_hydrophobic_weight,1/((2/res)/misplaced_hydrophobic_weight)) # weight of punishment
        
        # Find best score and punish for prolines in the sequence
        max(face-Charge-Polar-hydrophobic_misplaced-(proline_weight*length(unlist(str_extract_all(mer,'P')))))
      })
  
      # Punishment for short sequences -------------------------------------------
      Hm$len <- nchar(Hm$XXmers) 
      Hm$length_correction <- log(nchar(Hm$len), base=log_base)
      
      if(!is.null(fixed_len_threshold)){
        req <- fixed_len_threshold
      } else {
        req <- ifelse((nchar(sequence)*len_requirement)>11,nchar(sequence)*len_requirement,11)
      }
      
      Hm$length_correction[Hm$len>req] <- log(req, base=log_base)
    
      # Get helicity from netsurf3.0 or from predictions in array_object ---------
      heli_orig <- NULL
      
      if(use_netsurf==T){
        netsurf <- Netsurfp3(sequence,sequence) # function to access netsurf on biolib from python 
        
        # Add mean helicity for all possible peptides
        Hm$helicity <- unlist(lapply(min_len:max_len,function(l){
          m <- rollapply(netsurf$helix, l, by = 1, c)
          if(allow_prolines_temp==F){
            unlist(lapply(split(m, row(m)),mean,na.rm=T))[which(!str_detect(apply(rollapply(netsurf$seq, l, by = 1, c),1,paste,collapse=''),'P'))]
          } else {
            unlist(lapply(split(m, row(m)),mean,na.rm=T))
          }
        }))[keep]
        
        heli_orig <- 'NetSurfp3.0' # metadata
      }
      
      if(!is.null(UniProtKB.AC)){
        
        # Get position of sequence of site in protein
        start <- gregexpr(sequence,get_metadata(array_object, colnames = 'Full.Sequence', proteins = UniProt)[[1]][1])[[1]][1]
    
        # Get AlphaFold predictions for site
        AF <- get_metadata(array_object, colnames = 'AlphaFold', proteins = UniProt)[[1]]
        
        if(nrow(AF)>1){ # test if structure is available in array_object
          AF <- AF[AF$residue%in%start:(start+nchar(sequence)-1),]
          AF$confidence[AF$structure!='HELX_RH_AL_P'] <- 1 # set helical conf of non helical backbones to 1%
          
          # Add mean helicity for all possible peptides (AF conf divided by 100 to be between 0 and 1)
          Hm$helicity <- unlist(lapply(min_len:max_len,function(l){
            m <- rollapply(AF$confidence/100, l, by = 1, c)
            if(allow_prolines_temp==F){
              unlist(lapply(split(m, row(m)),mean,na.rm=T))[which(!str_detect(apply(rollapply(convert_aa(AF$aa), l, by = 1, c),1,paste,collapse=''),'P'))]
            } else {
              unlist(lapply(split(m, row(m)),mean,na.rm=T))
            }
          }))[keep]
          
          heli_orig <- 'AlphaFold2' # metadata
        } else {
          Hm$helicity <- 0.01
        }
      }
      
      # Calculation of final score -----------------------------------------------
      # Message produced if helicity has not been used for calculation!
      if(is.null(heli_orig)){
        cat('Helicity not calculated and included in optimization, specify UniProtKB.AC to use AlphaFold or set use_netsurf to T\n')
        Hm$score <- Hm$H.moment*Hm$segregation+Hm$length_correction 
      } else {
        Hm$score <- Hm$H.moment*Hm$segregation*Hm$helicity+Hm$length_correction 
      }
      
      # Order peptides by score --------------------------------------------------
      Hm <- Hm[order(Hm$score, decreasing = T),]
      Hm$is_helix <- T
      
      not_filtered <- Hm # backup prior to filtering
      not_filtered$is_helix <- F
      
      # If filters have been provided, peptides are removed here -----------------
      if(!is.null(set_min_moment)){
        Hm <- Hm[Hm$H.moment>=set_min_moment,]
      }
      
      if(!is.null(set_segregation_cutoff)){
        Hm <- Hm[Hm$segregation>=set_segregation_cutoff,]
      }
      
      # Creating object to be returned based on function input -------------------
      if(nrow(Hm)>=1){
        if(return_all_info==T){
          if(return_top_hit==F){
            res <- heliquest(Hm$XXmers[1], plot = plot)
            res <- c(res,list(XXmers=Hm, helicity_calculated_with=heli_orig))
          } else {
            res <- Hm[1,]
          }
        } else {
            res <- Hm$XXmers[1]
        }
      } else { # return empty object if no AH has been detected 
        if(return_all_info==T){
          if(return_top_hit==F){
            res <- c(list(coordinates=NA, H.moment=NA),list(XXmers=not_filtered,helicity_calculated_with=heli_orig))
          } else {
            res <- not_filtered[1,]
          }
        } else {
          res <- NA
        }
      }
    
    } else { # return empty object if no AH has been detected 
      if(return_all_info==T){
        if(return_top_hit==F){
          res <- c(list(coordinates=NA, H.moment=NA),list(XXmers=NA,helicity_calculated_with=NA))
        } else {
          res <- data.frame(XXmers=NA, H.moment=NA, segregation=NA, len=NA, length_correction=NA, score=NA, is_helix=NA)
        }
      } else {
        res <- NA
      }
    }
    
    # Returning final object ---------------------------------------------------
    return(res)
  
  })
  # Aggregating results --------------------------------------------------------  
  if(return_all_info==F){
    result <- unlist(result)
  } else {
    if(return_top_hit==T){
      result <- rbindlist(result)
    } else {
      if(length(sequences)==1){
        result <- result[[1]]
      }
    }
  }

  return(result)
  
}

# Plot Heliquest output
plot_helical_wheel <- function(heliquest_coordinates, 
                               correct_rotation=T, 
                               pH=7.4, 
                               density=F, 
                               add_noise=0, 
                               noise_sd=0.1, 
                               contour_var='count', 
                               bandwidth=0.3, 
                               adjust_bandwidth=1, 
                               n=100, 
                               plot=T, 
                               legend_position='none'){
  
  rotated_vec <- c(sum(heliquest_coordinates$Hydrophobicity*heliquest_coordinates$cos)/nrow(heliquest_coordinates),sum(heliquest_coordinates$Hydrophobicity*heliquest_coordinates$sin)/nrow(heliquest_coordinates))
  colrs <- setNames(as.character(heliquest_coordinates$color),heliquest_coordinates$AA)
  colrs <- colrs[unique(names(colrs))]
  
  if(correct_rotation==T){
    # Create rotation matrix to point hydrophobic moment vector down
    pivot <- (acos(abs(rotated_vec[1])/sqrt(abs(rotated_vec[1])^2+abs(rotated_vec[2])^2))*180/pi)
    if(rotated_vec[1]>0 & rotated_vec[2]>0){pivot <- 270-pivot}
    if(rotated_vec[1]<0 & rotated_vec[2]>0){pivot <- 90+pivot}
    if(rotated_vec[1]<0 & rotated_vec[2]<0){pivot <- 270-pivot-180}
    if(rotated_vec[1]>0 & rotated_vec[2]<0){pivot <- 270+pivot}
    
    rot_matrix <- matrix(c(cos(pivot*pi/180),sin(pivot*pi/180),-sin(pivot*pi/180),cos(pivot*pi/180)),ncol = 2, byrow = F)
    
    # Rotate Hm vector and residue positions
    rotated_vec <-    rot_matrix %*% rotated_vec 
    orginal_positions <- heliquest_coordinates[,c('cos','sin')]
    heliquest_coordinates$cos <- apply(orginal_positions,1,function(x){(rot_matrix %*% c(x[1],x[2]))[1,1]})
    heliquest_coordinates$sin <- apply(orginal_positions,1,function(x){(rot_matrix %*% c(x[1],x[2]))[2,1]})
  }
  
  # Calculate hydrophobicity vector sums ---------------------------------------
  original_vec <- c(sum(heliquest_coordinates$Hydrophobicity*heliquest_coordinates$cos)/nrow(heliquest_coordinates),sum(heliquest_coordinates$Hydrophobicity*heliquest_coordinates$sin)/nrow(heliquest_coordinates))
  moment <- sqrt(original_vec[1]^2 + original_vec[2]^2)
  
  if(density==F){
    fctr <- (0.375*(sort(rep((1:ceiling(nrow(heliquest_coordinates)/18)),18))[1:nrow(heliquest_coordinates)]-1))+1
    
    p <- ggplot(heliquest_coordinates, aes(cos*fctr,sin*fctr)) + 
      geom_path(aes(col=order,linewidth=order)) + 
      scale_linewidth(range=c(0.5,3)) +
      scale_color_gradient(low='grey80',high='grey30') + 
      ggnewscale::new_scale_color() +
      geom_segment(aes(x=start, y=start, xend=cos*Hydrophobicity ,yend=sin*Hydrophobicity),
                    arrow = arrow(length = unit(0.5, "cm"),type = 'open'),linetype=2,linewidth=.5) +
      geom_point(aes(col=AA),pch=16,size=20) +
      geom_segment(aes(x=0, y=0, xend=rotated_vec[1] ,yend=rotated_vec[2]),
                   arrow = arrow(length = unit(0.5, "cm"),type = 'closed'), linewidth=1.5) +
      geom_text(aes(label=AA)) +
      geom_text2(aes(label=c('N',rep('',nrow(heliquest_coordinates)-2),'C')), nudge_x = -0.16, nudge_y = -0.1,col='red') +
      prism() + remove_x_axis() + remove_y_axis() + theme(legend.position=legend_position) +
      xlim(-2.3,2.3) +
      ylim(-2.3,2.3) +
      scale_color_manual(values = colrs) + 
      annotate(geom="text", x=0, y=0.2, label=paste0('<µH> ',round(moment, digits = 3)), size = unit(8, "pt"))
  } else {

  if(add_noise>0){
    temp <- heliquest_coordinates
    for(i in 1:add_noise){
      heliquest_noised <- temp
      heliquest_noised$cos <- heliquest_noised$cos*(rnorm(n = 1, sd=.1)+1)
      heliquest_noised$sin <- heliquest_noised$sin*(rnorm(n = 1, sd=.1)+1)
      heliquest_coordinates <- rbind(heliquest_coordinates,heliquest_noised)
    }
  }
  
  p <- ggplot(heliquest_coordinates, aes(x=cos, y=sin))  + 
    stat_density_2d(geom = "polygon", aes(fill=color, alpha = after_stat(level), group=factor(color)), contour_var = contour_var ,adjust = adjust_bandwidth,n = n, h = bandwidth) +
    geom_segment(aes(x=0, y=0, xend=rotated_vec[1] ,yend=rotated_vec[2]),
                 arrow = arrow(length = unit(0.5, "cm"),type = 'closed'), linewidth=1.5) +
    prism() + remove_x_axis() + remove_y_axis() + 
    theme(legend.position=legend_position) +
    xlim(-2.3,2.3) +
    ylim(-2.3,2.3) +
    scale_alpha_continuous(range=c(0,1)) +
    scale_fill_manual(labels = unlist(lapply(lapply(split(heliquest_coordinates$AA,heliquest_coordinates$color),unique),paste,collapse=',')[ levels(factor(heliquest_coordinates$color))]),
                      values = levels(factor(heliquest_coordinates$color))) + 
    annotate(geom="text", x=0, y=0.2, label=paste0('<µH> ',round(moment, digits = 3)), size = unit(8, "pt"))
  }
  if(plot==T){
    print(p)
  }
  return(p)
}

# Change color of amino acids in helical wheel
change_helical_wheel_colors <- function(heliquest_coordinates,
                                        aa = NULL,
                                        colors = NULL,
                                        remove_color = F, 
                                        replacement_color = NULL){
  
  if(length(aa) > length(colors)){
    colors <- rep(colors,length(aa))
  }
  
  #levels(heliquest_coordinates$color) <- c(levels(heliquest_coordinates$color),colors)
  heliquest_coordinates$color <- as.character(heliquest_coordinates$color)
  
  for(i in 1:length(aa)){
    heliquest_coordinates$color[heliquest_coordinates$AA==aa[i]] <- colors[i]
  }
  
  if(remove_color==T & is.null(replacement_color)){
    groups <- list(c('R','K','H','E','D'),c('S','T','N','Q'),c('P','G','A'),c('I','L','V','C','M','W','F','Y'))
    groups <- lapply(groups,function(x){x[x%notin%aa]})
    grey_cols <- c('grey60','grey70','grey80','grey50')
    
    for(i in seq_along(groups)){
      heliquest_coordinates$color[heliquest_coordinates$AA%in%groups[[i]]] <- grey_cols[i]
    }
  }
  
  if(remove_color==T & !is.null(replacement_color)){
    heliquest_coordinates$color[heliquest_coordinates$AA%notin%aa] <- replacement_color
  }
  
  heliquest_coordinates$color <- as.factor(heliquest_coordinates$color)
  
  return(heliquest_coordinates)
  
}

# Move amino acids to front of helical wheel density plot
highlight_amino_acid_helical_wheel <- function(heliquest_coordinates,
                                               aa = NULL){
  target <- c(Amino_acids(aa),aa)
  heliquest_coordinates <- heliquest_coordinates[order(factor(heliquest_coordinates$AA, levels = target)),]
  return(heliquest_coordinates)
}

# Align sequences in helical wheel
helical_alignment <- function(alignment, 
                              coverage = NULL, 
                              start = NULL, 
                              end = NULL, 
                              add_noise = 0, 
                              noise_sd = 0.1, 
                              plot = F, 
                              legend_position = 'none'){
  
  terminals <- as.character(alignment)
  
  if(!is.null(coverage)){
    start <- do.call(rbind,sapply(terminals,str_split,pattern=''))
    start <- range(which(apply(start,2,function(x){sum(table(x[x!='-']))})>coverage*length(terminals)))
    
    if(!is.null(start)){
      start <- c(start,end)
    }
    
    terminals <- lapply(terminals,function(x){gsub('-','',substr(x,start[1],start[2]))})
    terminals <- terminals[unlist(lapply(terminals,nchar))>7]
  }
  
  
  
  terminals <- pblapply(terminals,heliquest,plot=F)
  terminals <- lapply(terminals,function(x){x[[1]]})
  
  merged <- do.call(rbind,terminals)
  merged <- rownames_to_column(merged, var='ID')
  merged$ID <- gsub('\\..*','',merged$ID)
  
  p <- plot_helical_wheel(merged,density = T,plot = F, add_noise = add_noise, noise_sd=noise_sd, legend_position=legend_position)
  
  if(plot==T){
    print(p)
  }
  
  Hm <- (sqrt(sum(merged$Hydrophobicity*merged$cos)^2 + sum(merged$Hydrophobicity*merged$sin)^2))/nrow(merged)
  
  return(list(coordinates=merged,H.moment=Hm,plot=p))
}

# Cluster helical alignment
cluster_helical_alignment    <- function(alignment, 
                                         centers = 3, 
                                         bins = 22, 
                                         n_components = 11, 
                                         n_components_kmeans = NULL,
                                         impute_with_constant = NULL,
                                         precision = 100,
                                         helical_plots_rows = 4,
                                         umap_config = umap.defaults,
                                         use_pacmap = T,
                                         aa = NULL,
                                         aa_color = NULL,
                                         point_size = 5,
                                         viridis_palette='H', 
                                         contour_var = 'count',
                                         ...){
      
      if(bins<max(nchar(alignment))){
        cat('\nToo few bins for all sequences, hydrophobicities that fall in the same bin will be summed! \n')
      }
      
      cat('Aligning helices ... \n')
      coords <- helical_alignment(alignment)$coordinates
      
      coords$radians <- sign(asin(coords$sin))*acos(-coords$cos)
      coords$intervals <-cut(coords$radians,breaks = bins)
      coords$diff <- c(1,diff(coords$order))
      coords$diff[coords$diff==1] <- 0
      coords$diff[!coords$diff==0] <- 1
      coords$diff <- cumsum(coords$diff)
      
      ints <- data.frame(intervals=levels(coords$intervals))
      
      cluster_coords <- split(coords,coords$diff)
      
      cat('Imputing data and creating matrix ... \n')
      cluster_coords <- do.call(rbind,pblapply(cluster_coords,function(x){
        
        x$intervals <- as.character(x$intervals)
        x <- merge(ints,x[,c(4,13)],all=T)
        
        x <- x %>% arrange(factor(intervals, levels = ints$intervals))
        
        x <- x %>%
          group_by(intervals)  %>%
          summarize(Hydrophobicity=sum(Hydrophobicity)) %>%
          ungroup()
        
        x$original <- x$Hydrophobicity
        x <- rbind(x,x,x) |> fill(Hydrophobicity, .direction = "downup")
        
        if(is.null(impute_with_constant)){
          for(i in 1:precision){
            x$Hydrophobicity <- smooth(x$Hydrophobicity,treat_start = F)
            x$Hydrophobicity[!is.na(x$original)] <- x$original[!is.na(x$original)] 
          }
        } else {
          
          x$Hydrophobicity[is.na(x$original)] <- impute_with_constant
          
        }
        
        pivot_wider(x[(bins+1):(2*bins),1:2], values_from = 'Hydrophobicity', names_from = 'intervals')
        
      }))
      
      cat('Running PCA ... \n')
      res <- prcomp(as.matrix(cluster_coords))
      
      if(is.na(n_components)){
        n_components <- bins
      }
      
      cat('Clustering PCA results using kmeans ... \n')
      if(is.null(n_components_kmeans)){
        cat('Setting no. of clusters for Kmeans to n_components \n')
        n_components_kmeans <- n_components
      }
      k_clusters <- kmeans(res$x[,1:n_components_kmeans], centers = centers, nstart = 1000, ...)
      
      if(use_pacmap==T){
        cat('Running PaCMAP ... \n')
        ump <- pacmap(res$x[,1:n_components])
        lab_x <- 'PaCMAP1'
        lab_y <- 'PaCMAP2'
      } else {
        cat('Running UMAP ... \n')
        ump <- umap::umap(res$x[,1:n_components], config = umap_config)
        ump <- as.data.frame(ump[['layout']])
        lab_x <- 'UMAP1'
        lab_y <- 'UMAP2'
      }
      
      ump$cluster <- as.factor(k_clusters$cluster)
      
      p <- ggplot(ump) + prism() + 
        geom_point(aes(x=V1,y=V2,col=cluster),size=point_size, pch=16) + 
        xlab(lab_x) + ylab(lab_y) + 
        scale_color_viridis_d(option = viridis_palette)
      
      ump$sequences <- alignment
      
      cat('Creating helical plots ... \n')
      helical_plots_only <- lapply(1:centers,function(i){
        plot_cluster <- ump$sequences[ump$cluster==i]
        plot_only <- plot_helical_alignment(plot_cluster, plot = F, aa = aa, aa_color = aa_color, contour_var = contour_var)$plot
        plot_only + annotate(geom="text", x=0, y=2, label=paste0('Cluster: ',i),
                             color=viridis(n=centers, option = viridis_palette)[i], size = unit(8, "pt"))
      })
      
      cat('Plotting ... \n')
      p2 <- do.call("arrangeGrob", c(helical_plots_only,ncol=ceiling(centers/helical_plots_rows)))
      
      p_final <- grid.arrange(arrangeGrob(p, ncol=1, nrow=1),
                              arrangeGrob(p2, ncol=1), heights=c(5,1), widths=c(3.5,1.5))
      
      return(list(input=alignment, matrix=as.matrix(cluster_coords), pca=res, 
                  kmeans=k_clusters, umap=ump, umap_plot=p, 
                  helical_plots=helical_plots_only, overview_plot=p_final))
      
}

re_cluster_helical_alignment <- function(cluster_output,
                                         centers = 3,
                                         n_components = 11,
                                         n_components_kmeans = NULL,
                                         helical_plots_rows = 4,
                                         umap_config = umap.defaults,
                                         use_pacmap = T,
                                         umap_only = F,
                                         aa = NULL,
                                         aa_color = NULL,
                                         point_size = 5,
                                         viridis_palette = 'H',
                                         contour_var='count', ...){
  
  cluster_coords <- cluster_output$matrix
  
  res <- cluster_output$pca
  
  if(is.na(n_components)){
    n_components <- ncol(cluster_coords)
  }
  
  if(n_components > ncol(cluster_coords)){
    cat('Too many n_components requested, setting to max ... \n')
    n_components <- ncol(cluster_coords)
  }
  
  cat('Clustering PCA results using kmeans ... \n')
  if(is.null(n_components_kmeans)){
    cat('Setting no. of clusters for Kmeans to n_components\n')
    n_components_kmeans <- n_components
  }
  k_clusters <- kmeans(res$x[,1:n_components_kmeans], centers = centers, nstart = 1000, ...)
  
  if(use_pacmap==T){
    cat('Running PaCMAP ... \n')
    ump <- pacmap(res$x[,1:n_components])
    lab_x <- 'PaCMAP1'
    lab_y <- 'PaCMAP2'
  } else {
    cat('Running UMAP ... \n')
    ump <- umap::umap(res$x[,1:n_components], config = umap_config)
    ump <- as.data.frame(ump[['layout']])
    lab_x <- 'UMAP1'
    lab_y <- 'UMAP2'
  }
  
  ump$cluster <- as.factor(k_clusters$cluster)
  
  p <- ggplot(ump) + prism() + 
    geom_point(aes(x=V1,y=V2,col=cluster),size=point_size, pch=16) + 
    xlab(lab_x) + ylab(lab_y) + 
    scale_color_viridis_d(option = viridis_palette)
  
  ump$sequences <- cluster_output$input
  
  if(umap_only==T){
    cat('Plotting ... \n')
    print(p)
    p_final <- NA
    helical_plots_only <- NA
  } else {
  
  cat('Creating helical plots ... \n')
  helical_plots_only <- lapply(1:centers,function(i){
    plot_cluster <- ump$sequences[ump$cluster==i]
    plot_only <- plot_helical_alignment(plot_cluster, plot = F, aa = aa, aa_color = aa_color, contour_var = contour_var)$plot
    plot_only + annotate(geom="text", x=0, y=2, label=paste0('Cluster: ',i),
                         color=viridis(n=centers, option = viridis_palette)[i], size = unit(8, "pt"))
  })
  
  cat('Plotting ... \n')
  p2 <- do.call("arrangeGrob", c(helical_plots_only,ncol=ceiling(centers/helical_plots_rows)))
  
  p_final <- grid.arrange(arrangeGrob(p, ncol=1, nrow=1),
                          arrangeGrob(p2, ncol=1), heights=c(5,1), widths=c(3.5,1.5))
  
  }
  
  return(list(input= cluster_output$input, matrix=cluster_coords, pca=res, 
              kmeans=k_clusters, umap=ump, umap_plot=p, 
              helical_plots=helical_plots_only, overview_plot=p_final))
  
  
}

# Functional analysis of helical alignment
label_cluster_helical_alignment_output_with_COMPARTMENTS <- function(clustered_alignment, 
                                                                     term, 
                                                                     confidence=4, 
                                                                     point_size = 5, 
                                                                     array_object=NULL){
  
  if('UniProtKB.AC'%notin%colnames(clustered_alignment)){
    cat('Error - Please add a column to the output with the UniProt IDs and call it "UniProtKB.AC"\n')
    stop()
  }
  
  if(is.null(array_object)){
    array_object <- Find_array_object()
  }
  
  all_comps <- get_metadata(colnames = 'COMPARTMENTS',proteins = clustered_alignment$UniProtKB.AC, array_object = array_object)
  
  all_comps <- stack(lapply(all_comps, function(x){x$GO_term [x$Confidence>=confidence]}))
  
  if(term%notin%all_comps$values){
    cat('The requested term does not exist in COMPARTMENTS database ... \n')
    stop()
  }
  
  sele <- unique(all_comps[all_comps$values==term,])
  colnames(sele) <- c('term','UniProtKB.AC')
  to_plot <- merge(x=clustered_alignment,y=sele, all.x=T, by='UniProtKB.AC')
  
  in_group <- to_plot[!is.na(to_plot$term),]
  
  p <- ggplot() + 
    geom_point(data = to_plot, aes(x=V1, y=V2, col=cluster), alpha=0.1, size = point_size,pch=16) + 
    geom_point(data=in_group, aes(x=V1, y=V2, fill=term)) +
    prism()
  print(p)
  return(list(df=to_plot,plot=p))
}

functional_enrichment_cluster_helical_alignment <- function(SWaFi_score,
                                                            cluster_helical_alignment_output,
                                                            background = species_filter(), 
                                                            COMPARTMENTS = T, 
                                                            PANTHER = F,
                                                            omit_terms_including = 'organelle|system|cytoplasm|structure',
                                                            PANTHER_database = 'GO',
                                                            PANTHER_select_sub_database = NULL,
                                                            compartments_threshold = 2, 
                                                            pthresh = 1, 
                                                            pAdjust = 'fdr',
                                                            array_object = NULL){
  
  helix_clusters <- cbind(SWaFi_score, cluster_helical_alignment_output$umap)
  
  cluster_OR <- do.call(rbind,lapply(1:max(as.numeric(helix_clusters$cluster)), function(i){
    cluster <- strip_terminal_specifier(helix_clusters$UniProtKB.AC[helix_clusters$cluster==i])
    if(PANTHER==T){
      OR <- quiet(get_panther_overrepresentation(cluster,Background = background, database = PANTHER_database, cutoff = pthresh, pAdjust = pAdjust))
      if(!is.null(PANTHER_select_sub_database)){
        OR <- OR[OR$database%in%PANTHER_select_sub_database,]
      }
    }
    
    if(COMPARTMENTS==T & PANTHER==F){
      
      if(is.null(array_object)){
        array_object <- Find_array_object()
      }
      
      OR <- COMPARTMENTS_overrepresentation(cluster,background = background, compartments_threshold = compartments_threshold, p_thresh = pthresh, pAdjust = pAdjust, array_object = array_object)
      OR <- OR[OR$adjusted<pthresh,]
    }
    
    OR <- OR[!str_detect(OR$term_label,omit_terms_including),]
    OR$database <- paste0('Cluster ',i)
    return(OR)  
  }))
  
  return(cluster_OR)
}

# Align helical wheels of aligned sequences
plot_helical_alignment <- function(alignment, 
                                   coverage=NULL, 
                                   start=NULL, 
                                   end=NULL, 
                                   add_noise=0, 
                                   noise_sd=0.1, 
                                   plot=T, 
                                   legend_position='none', 
                                   aa=NULL, 
                                   aa_color=NULL, 
                                   contour_var='count'){
  
  terminals <- as.character(alignment)
  
  if(!is.null(coverage)){
    start <- do.call(rbind,sapply(terminals,str_split,pattern=''))
    start <- range(which(apply(start,2,function(x){sum(table(x[x!='-']))})>coverage*length(terminals)))
    
    if(!is.null(start)){
      start <- c(start,end)
    }
    
    terminals <- lapply(terminals,function(x){gsub('-','',substr(x,start[1],start[2]))})
    terminals <- terminals[unlist(lapply(terminals,nchar))>7]
  }
  
  
  
  terminals <- pblapply(terminals,heliquest,plot=F)
  terminals <- lapply(terminals,function(x){x[[1]]})
  
  merged <- do.call(rbind,terminals)
  merged <- rownames_to_column(merged, var='ID')
  merged$ID <- gsub('\\..*','',merged$ID)
  
  if(!is.null(aa)){
    cat(paste0('Changing amino acid color for ',paste(aa,collapse = ', '),' to ',paste(aa_color, collapse = ', '),'\n'))
    merged <- change_helical_wheel_colors(merged,aa = aa, colors = aa_color, remove_color = F)
  }
  
  p <- plot_helical_wheel(merged,density = T,plot = F, add_noise = add_noise, noise_sd=noise_sd, legend_position=legend_position, contour_var = contour_var)
  
  if(plot==T){
    print(p)
  }
  
  Hm <- (sqrt(sum(merged$Hydrophobicity*merged$cos)^2 + sum(merged$Hydrophobicity*merged$sin)^2))/nrow(merged)
  
  return(list(coordinates=merged,H.moment=Hm,plot=p))
}

# Weighted kde2d function
# Adapted from: https://stat.ethz.ch/pipermail/r-help/2006-June/107405.html
kde2d.weighted <- function(x, y, h, n=25, lims=c(range(x),range(y)), w=1){
  nx <- length(x)
  if (length(y) != nx) { stop("data vectors must be the same length")}
  if (any(!is.finite(x)) || any(!is.finite(y))) {
    stop("missing or infinite values in the data are not allowed")
    }
  if (any(!is.finite(lims))) {
    stop("only finite values are allowed in 'lims'")
  }
  n <- rep(n, length.out = 2L)
  gx <- seq.int(lims[1L], lims[2L], length.out = n[1L])
  gy <- seq.int(lims[3L], lims[4L], length.out = n[2L])
  h <- if (missing(h)) {
    c(bandwidth.nrd(x), bandwidth.nrd(y))
  } else {rep(h, length.out = 2L)}
  if (any(h <= 0)) { stop("bandwidths must be strictly positive")}
  h <- h/4
  ax <- outer(gx, x, "-")/h[1L]
  ay <- outer(gy, y, "-")/h[2L]
  ax <- matrix(dnorm(ax), ncol = nx)*matrix(rep(w,n[1L]), nrow=n, ncol=nx, byrow=TRUE) 
  ay <- matrix(dnorm(ay), ncol = nx)
  z <- (tcrossprod(ax,ay)/(nx * h[1L] * h[2L]))
  list(x = gx, y = gy, z = z)
}

# plot 2d density plots of amino acids in aligned helices
# The input to this function is the output of 'plot_helical_alignment()'
plot_2d_helical_density <- function(alignment_1, 
                                    alignment_2 = NULL,
                                    aa = c('E','D'), 
                                    colr = 'red', 
                                    mid = '#FFFFBF', 
                                    equal_col_scale = NULL, 
                                    binwidth = 0.001,
                                    viridis_option = 'D', 
                                    viridis_direction = 1, 
                                    plot = T, 
                                    aspect_ratio = 1, 
                                    list_weights_for_alignments = NULL,
                                    weights_column = NULL, 
                                    bw.x = NULL, 
                                    bw.y = NULL, 
                                    yrng = NULL, 
                                    n = 400){
  
  b  <- alignment_1$coordinates
  b$radians <- sign(asin(b$sin))*acos(-b$cos) # minus sin is necessary to keep helix lefthanded, since we are plotting from N -> C and not C -> N 
  
  if(!is.null(alignment_2)){
    nb <- alignment_2$coordinates
    nb$radians <- sign(asin(nb$sin))*acos(-nb$cos)
    xrng <- range(c(b$radians[b$AA%in%aa], nb$radians[nb$AA%in%aa]))
    
    if(is.null(yrng)){
      cat('Calculating common range ... \n')
      yrng <- range(c(0, (1+max(c(b$order,nb$order)))))
    }
    
    if(is.null(bw.x)){
      cat('Estimating bandwidths ... \n')
      # Calculate the 2d density estimate over the common range
      bw.x <- bandwidth.nrd(c(b$radians[b$AA%in%aa],nb$radians[nb$AA%in%aa]))
      bw.y <- bandwidth.nrd(c(b$order[b$AA%in%aa],nb$order[nb$AA%in%aa]))
    }
  } else {
    
    xrng <- range(b$radians[b$AA%in%aa])
    if(is.null(yrng)){
      cat('Calculating range ... \n')
      yrng <- range(c(0,(1+max(b$order))))
    }
    
    if(is.null(bw.x)){
      cat('Estimating bandwidths ... \n')
      # Calculate the 2d density estimate over the common range
      bw.x <- bandwidth.nrd(b$radians[b$AA%in%aa])
      bw.y <- bandwidth.nrd(b$order[b$AA%in%aa])
    }
  }
  
  weights_added <- F
  if(!is.null(list_weights_for_alignments)){
    
    if(is.list(list_weights_for_alignments)==T){
      weights_added <- T
      b_weights <- list_weights_for_alignments[[1]][b$AA%in%aa]
      if(!is.null(alignment_2)){
        if(length(list_weights_for_alignments)>1){
          nb_weights <- list_weights_for_alignments[[2]][b$AA%in%aa]
        } else {
          cat('No weights provided for alignment_2 ... \n')
          stop()
        }
      }
      cat('Weights have been added!\n')
    } else {
      cat('Weights should be provided as a list containing weights for alignment_1 and if needed alignment_2 ... \n')
      stop()
    }
  } 
  
  if(!is.null(weights_column) & weights_added==F){
    
    if(weights_column%in%colnames(b)){
      b_weights <- b[b$AA%in%aa,weights_column]
      if(!is.null(alignment_2)){
        if(weights_column%in%colnames(nb)){
          nb_weights <- nb[nb$AA%in%aa,weights_column]
        } else {
          cat('Column name provided for weights does not exist in alignment_2$coordinates ... \n')
          stop() 
        }
      }
      weights_added <- T
      cat('Weights have been added!\n')
    } else {
      cat('Column name provided for weights does not exist in alignment_1$coordinates ... \n')
      stop()
    }
  }
  
  if(weights_added==F){
    b_weights <- 1
    if(!is.null(alignment_2)){
      nb_weights <- 1
    }
    cat('No weights provided, regular density will be calculated ... \n')
  }
  
  cat('Running Kernel Density Estimation ... \n')
  d1 <- kde2d.weighted(b$radians[b$AA%in%aa] , b$order[b$AA%in%aa],  lims=c(xrng, yrng), n=n, h=c(bw.x,bw.y),w=b_weights)
  if(!is.null(alignment_2)){
    d2 <- kde2d.weighted(nb$radians[nb$AA%in%aa], nb$order[nb$AA%in%aa], lims=c(xrng, yrng), n=n, h=c(bw.x,bw.y),w = nb_weights)
    xiden <- identical(d1$x, d2$x) 
    yiden <- identical(d1$y, d2$y) # test if ranges are the same
    cat(paste0('- x ranges identical: ',xiden,'\n- y ranges identical: ',yiden,'\n'))
  }
  
  cat('Converting to relative density ... \n')
  b_rel_amount <- nrow(b[b$AA%in%aa,])/nrow(b)
  d1$z <- d1$z*b_rel_amount
  
  if(!is.null(alignment_2)){
    nb_rel_amount <- nrow(nb[nb$AA%in%aa,])/nrow(nb) 
    d2$z <- d2$z*nb_rel_amount
  }
  
  # Calculate the difference between the 2d density estimates
  if(!is.null(alignment_2)){
    cat('Calculating Density Difference ... \n')
    diff <- d1
    diff$z  <- (d1$z-d2$z)/max(c(d1$z,d2$z))
  }
  
  cat('Melting Data into long format... \n')
  ## Melt data into long format
  # First, add row and column names (x and y grid values) to the z-value matrix
  rownames(d1$z) <- d1$x
  colnames(d1$z) <- d1$y
  d1 <- melt(d1$z, id.var=rownames(d1))
  names(d1) <- c("location","order","z")
  
  if(!is.null(alignment_2)){
    ## Melt data into long format
    # First, add row and column names (x and y grid values) to the z-value matrix
    rownames(d2$z) <- d2$x
    colnames(d2$z) <- d2$y
    d2 <- melt(d2$z, id.var=rownames(d2))
    names(d2) <- c("location","order","z")
    
    ## Melt data into long format
    # First, add row and column names (x and y grid values) to the z-value matrix
    rownames(diff$z) <- diff$x
    colnames(diff$z) <- diff$y
    diff <- melt(diff$z, id.var=rownames(diff))
    names(diff) <- c("location","order","z")
  }
  
  if(binwidth>(max(d1$z)/4)){
    cat('\nImportant! \nYou have set the binwidth very low, stat_contour() might not be able to display the densities ... \n\n')
  }
  
  if(!is.null(alignment_2)){
    cat('Generating plots ... \n')
    
    if(is.null(equal_col_scale)){
      lim_max <- max(c(d1$z,d2$z))
    } else {
      lim_max <- equal_col_scale
    }
    p_d1 <- ggplot() +
      geom_tile(aes(x=location, y=order, fill=z),data=d1) + 
      stat_contour(aes(x=location,y=order,z=z,colour=after_stat(level)),binwidth=binwidth,data=d1) +
      scale_fill_gradient2(  low="white",mid=mid, high=colr, midpoint=lim_max/2, limits=c(0,lim_max)) +
      scale_colour_gradient2(low="white",mid=mid, high=colr, midpoint=lim_max/2, limits=c(0,lim_max)) +
      prism() + xlab('Phase shift') + ylab(expression('Position in helix (N '  %->% 'C)')) +
      theme(legend.position = 'none', 
            aspect.ratio = aspect_ratio, 
            axis.text = element_text(size = 22), 
            axis.title = element_text(size = 22),
            title = element_text(size=24)) +
      geom_vline(xintercept = c(-acos(0),0,acos(0)),lty=2) +
      scale_x_continuous(breaks = c(-pi,-acos(0),0,acos(0),pi),
                         labels= c(expression(-pi),expression(frac(-pi,2)),expression(0),expression(frac(pi,2)),expression(pi))) +
      scale_y_continuous(breaks = round(seq(1, yrng[2]-1, by = 1),1)) + 
      guides(colour='none') 
    
    p_d2 <- ggplot() +
      geom_tile(aes(x=location, y=order, fill=z),data=d2) + 
      stat_contour(aes(x=location,y=order,z=z,colour=after_stat(level)),binwidth=binwidth,data=d2) +
      scale_fill_gradient2(low="white",mid=mid, high=colr, midpoint=lim_max/2, limits=c(0,lim_max)) +
      scale_colour_gradient2(low="white",mid=mid,high=colr,midpoint=lim_max/2, limits=c(0,lim_max)) +
      prism() + xlab('Phase shift') + ylab(expression('Position in helix (N '  %->% 'C)')) +
      theme(legend.position = 'none',
            aspect.ratio = aspect_ratio, 
            axis.text = element_text(size = 22), 
            axis.title = element_text(size = 22),
            title = element_text(size=24)) +
      geom_vline(xintercept = c(-acos(0),0,acos(0)),lty=2) +
      scale_x_continuous(breaks = c(-pi,-acos(0),0,acos(0),pi),
                         labels= c(expression(-pi),expression(frac(-pi,2)),expression(0),expression(frac(pi,2)),expression(pi))) +
      scale_y_continuous(breaks = round(seq(1, yrng[2]-1, by = 1),1)) + 
      guides(colour='none')
    
    #set colours
    diff_col <- viridis(3,option=viridis_option, direction = viridis_direction)
    #colour limits
    lims <- c(abs(min(diff$z)),max(diff$z))
    lims <- lims[which(lims==max(lims))]
    
    # Plot density difference between binding and non-binding 
    p_diff <- ggplot() +
      geom_tile(aes(x=location, y=order, fill=z),data=diff) + 
      stat_contour(aes(x=location,y=order,z=z,colour=after_stat(level)),binwidth=binwidth,data=diff) +
      scale_fill_gradient2(  low=diff_col[1], mid='white', high=diff_col[3], midpoint=0, limits=c(-lims,lims)) +
      scale_colour_gradient2(low=diff_col[1], mid='white', high=diff_col[3], midpoint=0, limits=c(-lims,lims)) +
      prism() + xlab('Phase shift') + ylab(expression('Position in helix (N '  %->% 'C)')) +
      theme(legend.position = 'none', 
            aspect.ratio = aspect_ratio, 
            axis.text = element_text(size = 22), 
            axis.title = element_text(size = 22),
            title = element_text(size=24)) +
      geom_vline(xintercept = c(-acos(0),0,acos(0)),lty=2) +
      scale_x_continuous(breaks = c(-pi,-acos(0),0,acos(0),pi),
                         labels= c(expression(-pi),expression(frac(-pi,2)),expression(0),expression(frac(pi,2)),expression(pi))) +
      scale_y_continuous(breaks = round(seq(1, yrng[2]-1, by = 1),1)) + 
      guides(colour='none') 
    
    p_diff_legend <- p_diff + theme(legend.position = 'bottom', legend.text = element_text(size=18, angle = 45, hjust = 1))
    p_d1_legend   <- p_d1   + theme(legend.position = 'bottom', legend.text = element_text(size=18, angle = 45, hjust = 1))
    p_d2_legend   <- p_d2   + theme(legend.position = 'bottom', legend.text = element_text(size=18, angle = 45, hjust = 1))
    
    p_diff_legend <- grab_legend(p_diff_legend)
    p_d1_legend   <- grab_legend(p_d1_legend)
    p_d2_legend   <- grab_legend(p_d2_legend)
    
    p <- list(plots=list(alignment1=p_d1,alignment2=p_d2,difference=p_diff),
              legends=list(alignment1=p_d1_legend,alignment2=p_d2_legend,difference=p_diff_legend),
              density_data=list(alignment1=d1,alignment2=d2,difference=diff),
              coordinates =list(b,nb),
              kde2_parameters=list(xrng=xrng,yrng=yrng,bw.x=bw.x,bw.y=bw.y,n=n))
    
    if(plot==T){
      grid.draw(do.call("grid.arrange", c(p$plots, ncol=3)))
    }
    
  } else {
    cat('Generating plot ... \n')
    
    if(is.null(equal_col_scale)){
      lim_max <- max(d1$z)
    } else {
      lim_max <- equal_col_scale
    }
    
    p_d1 <- ggplot() +
      geom_tile(aes(x=location, y=order, fill=z),data=d1) + 
      stat_contour(aes(x=location,y=order,z=z,colour=after_stat(level)),binwidth=binwidth,data=d1) +
      scale_fill_gradient2(  low="white",mid=mid, high=colr, midpoint=lim_max/2, limits=c(0,lim_max)) +
      scale_colour_gradient2(low="white",mid=mid, high=colr, midpoint=lim_max/2, limits=c(0,lim_max)) +
      prism() + xlab('Phase shift') + ylab(expression('Position in helix (N '  %->% 'C)')) +
      theme(legend.position = 'none', 
            aspect.ratio = aspect_ratio, 
            axis.text = element_text(size = 22), 
            axis.title = element_text(size = 22),
            title = element_text(size=24), 
            legend.text = element_text(size=24)) +
      geom_vline(xintercept = c(-acos(0),0,acos(0)),lty=2) +
      scale_x_continuous(breaks = c(-pi,-acos(0),0,acos(0),pi),
                         labels= c(expression(-pi),expression(frac(-pi,2)),expression(0),expression(frac(pi,2)),expression(pi))) +
      scale_y_continuous(breaks = round(seq(1, 11, by = 1),1)) + 
      guides(colour='none')
    
    p_d1_legend   <- p_d1   + theme(legend.position = 'bottom')
    p_d1_legend   <- grab_legend(p_d1_legend)
    
    p <- list(plots=list(alignment1=p_d1),
              legends=list(alignment1=p_d1_legend),
              density_data=list(alignment1=d1),
              coordinates =list(b),
              kde2_parameters=list(xrng=xrng,yrng=yrng,bw.x=bw.x,bw.y=bw.y,n=n))
    
    if(plot==T){
      grid.draw(p_d1)
    }
    
  }
  cat('Done!\n')
  return(p)
  
}

# Plot 3d helical density of helix
plot_3d_helical_density <- function(heliquest_coordinates, 
                                    aa_group_list=NULL, 
                                    group_colors=NULL,  
                                    remove_amino_acids=NULL, 
                                    add=F, 
                                    alpha_power=5, 
                                    radius=0.5, 
                                    jitter=T){
  
  if(is.null(aa_group_list)){
    aa_group_list <- list(c('E','D'), c('K','R','H'),c('W','F','Y'),c('L','I','V'),c('T','S','N','Q'),c('A','G'),c('C','M'),c('P'))
  }
  
  if(is.null(group_colors)){
    group_colors <- c('red','blue','orange','yellow','purple','grey50','turquoise','forestgreen')
  }
  
  if(length(group_colors) != length(aa_group_list)){
    print('Please provide the same number of colors and Amino acid groups ... ')
    stop()
  }
  
  o <- heliquest_coordinates
  
  for(i in 1:length(aa_group_list)){
    aa <- aa_group_list[[i]]
    aa <- aa[aa%notin%remove_amino_acids]
    
    if(length(aa)>0){
      colfunc <- colorRampPalette(c("white", group_colors[i]))
      d <- kde3d(o$cos[o$AA%in%aa], o$sin[o$AA%in%aa], o$order[o$AA%in%aa], n = c(100,100,11), lims = c(-1.5,1.5,-1.5,1.5,1,11))
      misc3d::image3d(d$d, d$x, d$y, d$z, alpha.power = alpha_power, col =  colfunc(9)  ,jitter = jitter, sprites = T, radius = radius, add=add)
      add <- T
    }
  }
  
  names(aa_group_list) <- group_colors  
  
  aa_group_list <- lapply(aa_group_list,function(x){x[x%notin%remove_amino_acids]})
  aa_group_list <- aa_group_list[unlist(lapply(aa_group_list,length))>0]
  
  print(aa_group_list)
  
}

# Function for determining electrochemical properties of peptides
PeptideProperties <- function(df, 
                              FASTAColumn, 
                              pH=7.4){
  
  df <- as.data.frame(df)
  
  #Reading in aa hydrophobicity data from file
  AA_Properties <- AA_Properties()
  AA_Properties <- AA_Properties[order(AA_Properties$AA),]
  
  z_aa <- AA_Properties$AA[AA_Properties$z != 0]
  
  Properties <- rbindlist(sapply(df[,FASTAColumn],function(x){
    
    if(is.na(x)){
      result <- list(data.frame(Hydrophobicity=NA, H.moment=NA,
                                z=NA, z.cont=NA, No.of.Charges=NA,
                                Proline.residues=NA))
      return(result)
    } else {
      AA <- unlist(str_split(x,''))
      AA <- AA[AA != 'X'] 
      charges <- table(c(AA[AA%in%z_aa],z_aa))-1
      no <- sum(charges[-3])
      z <- sum(charges*c(-1,-1,1,1))
      z.cont <- sum(charges*AA_Properties$z.cont[AA_Properties$AA%in%z_aa])
      prolines <- sum(AA == "P")
      H <- sum((table(c(AA,AA_Properties$AA))-1)*AA_Properties$Hydrophobicity)/length(AA)
      Hm <- data.frame(H=sapply(AA,function(j){
        AA_Properties$Hydrophobicity[AA_Properties$AA==j]
      }),
      position=1:length(AA))
      
      Hm$sin <- (Hm$H*sin(Hm$position*(100*(pi/180))))
      Hm$cos <- (Hm$H*cos(Hm$position*(100*(pi/180))))
      Hm <- (1/(length(AA))) * (sum(Hm$sin)^2 + sum(Hm$cos)^2)^0.5
      
      result <- list(data.frame(Hydrophobicity=H, H.moment=Hm,
                                z=z, z.cont=z.cont, No.of.Charges=no,
                                Proline.residues=prolines))
      return(result)
    }
  }))
  
  df <- cbind(df, Properties)
  return(df)
}

# Function for converting .fasta to dataframe
FASTA_to_csv <- function(FASTA){
  FASTA <- as.data.frame(FASTA)
  FASTA <- t(FASTA)
  FASTA <- rownames_to_column(as.data.frame(FASTA))
  colnames(FASTA) <- c("ID", "Sequence")
  FASTA
}

# Blast function for local protein blast
pBlast <- function(db = "swissprot", 
                   path=getwd()){
  
  # Getting packages for parallelisation
  if(!require(foreach)) install.packages("foreach")
  library(foreach)
  if(!require(doParallel)) install.packages("doParallel")
  library(doParallel)
  
  # Getting desired fasta files from folder
  infileNms <- list.files(path=path, pattern = "*.fasta")
  
  # Determining names of output files
  outnames <- paste(path,"/",unlist(sapply(infileNms, strsplit, split = "*.FASTA")), 
                    ".out", sep = "")
  
  # Creating folder for output files
  dir.create(paste0(path,"/pBlast_Outfiles"))
  
  # create commands function
  Blastp <- function(infile, outfile, database, threads){
    paste0("blastp -db ", database," -query ", infile, " -num_threads ", threads,
           " -outfmt 6 -out pBlast_Outfiles/", outfile)
  }
  
  # Determining Optimal Core usage for Parallelisation
  subcl <- makeCluster((detectCores()-1)/length(infileNms))
  
  # create actual commands
  cmds <- mapply(FUN = Blastp, infile = infileNms, outfile = outnames, 
                 database = db, threads = as.character(length(subcl)))
  
  #Defining Clusters
  cl <- makeCluster(length(cmds))
  registerDoParallel(cl)
  
  # Blasting in cmd prompt.
  foreach( i = 1:length(cmds)) %dopar% {
    sapply(cmds[i], system)
  }
  
  # Importing outfiles to R
  Hits <- foreach(i = 1:length(outnames)) %dopar% {
    temp <- read.table(paste0("pBlast_Outfiles/",outnames[i]), header = F)
    colnames(temp) <- c("query", "subject", "%identity", "alignment.length", 
                        "mismatches", "gap.opens", "q.start", "q.end", "s.start", 
                        "s.end", "evalue", "bit.score")
    temp
  }
  
  # Stopping cluster to save resources
  stopCluster(cl)
  
  # Renaming list elements
  names(Hits) <- outnames
  
  # Pulling individual data.frames from list to Env.
  list2env(Hits, envir = .GlobalEnv)
}

# Function returns the Jaccard index and Jaccard distance
jaccard <- function(df, 
                    margin) {
  if (margin == 1 | margin == 2) {
    M_00 <- apply(df, margin, sum) == 0
    M_11 <- apply(df, margin, sum) == 2
    if (margin == 1) {
      df <- df[!M_00, ]
      JSim <- sum(M_11) / nrow(df)
    } else {
      df <- df[, !M_00]
      JSim <- sum(M_11) / length(df)
    }
    JDist <- 1 - JSim
    return(c(JSim = JSim, JDist = JDist))
  } else break
}

# Making custom function to create xx-mers
AAmers <- function(string.list, 
                   resolution = 1, 
                   pep.seq){
  
  peptide.list <- foreach(i = 1:(length(string.list)-(pep.seq-1)), .combine = c,
                          .multicombine = T, .inorder = T) %do% {
                            peptide <- string.list[i:(i+(pep.seq-1))]
                            peptide <- paste(peptide, collapse = '')
                          }
}

# Gene set enrichment analysis
GSEA <- function(gene_list, 
                 GO_list, 
                 pvalue) {
  
  if ( any( duplicated(names(gene_list)) )  ) {
    warning("Duplicates in gene names")
    gene_list = gene_list[!duplicated(names(gene_list))]
  }
  if  ( !all( order(gene_list, decreasing = TRUE) == 1:length(gene_list)) ){
    warning("Gene list not sorted")
    gene_list = sort(gene_list, decreasing = TRUE)
  }
  #myGO = fgsea::gmtPathways(GO_file)
  
  myGO <- GO_list
  
  fgRes <- fgsea::fgsea(pathways = myGO,
                        stats = gene_list,
                        minSize=5, ## minimum gene set size
                        maxSize=500, ## maximum gene set size
                        nPermSimple = 10000) %>%  
    as.data.frame() %>% 
    dplyr::filter(!is.na(pval)) %>% 
    dplyr::filter(padj < pvalue) %>% 
    arrange(desc(NES))
  
  message(paste("Number of signficant gene sets =", nrow(fgRes)))
  
  message("Collapsing Pathways -----")
  concise_pathways = collapsePathways(data.table::as.data.table(fgRes),
                                      pathways = myGO,
                                      stats = gene_list[names(gene_list)%in%unlist(fgRes$leadingEdge)])
  fgRes = fgRes[fgRes$pathway %in% concise_pathways$mainPathways, ]
  message(paste("Number of gene sets after collapsing =", nrow(fgRes)))
  
  fgRes$Enrichment = ifelse(fgRes$NES > 0, "Up-regulated", "Down-regulated")
  filtRes = rbind(head(fgRes, n = 10),
                  tail(fgRes, n = 10 ))
  
  total_up = sum(fgRes$Enrichment == "Up-regulated")
  total_down = sum(fgRes$Enrichment == "Down-regulated")
  header = paste0("Top 10 (Total pathways: Up=", total_up,", Down=",    total_down, ")")
  
  colos = setNames(c("firebrick2", "dodgerblue2"),
                   c("Up-regulated", "Down-regulated"))
  
  g1= ggplot(filtRes, aes(reorder(pathway, NES), NES)) +
    geom_point( aes(fill = Enrichment, size = size), shape=21) +
    scale_fill_manual(values = colos ) +
    scale_size_continuous(range = c(2,10)) +
    geom_hline(yintercept = 0) +
    coord_flip() +
    labs(x="Pathway", y="Normalized Enrichment Score",
         title=header) + theme_minimal() 
  
  
  output = list("Results" = fgRes, "Plot" = g1)
  return(output)
}

plot_geneset_clusters = function(gs_results, 
                                 GO_list, 
                                 min.sz=4, 
                                 main="GSEA clusters"){
  
  myGO = GO_list
  df = matrix(nrow=nrow(gs_results), ncol = nrow(gs_results), data = 0)
  rownames(df) = colnames(df) = gs_results$pathway
  
  for ( i in 1:nrow(gs_results)) {
    genesI =  unlist(myGO[names(myGO) == gs_results$pathway[i] ])
    for (j in 1:nrow(gs_results)) {
      genesJ = unlist(myGO[names(myGO) == gs_results$pathway[j] ])
      ## Jaccards distance  1 - (intersection / union )
      overlap = sum(!is.na(match(genesI, genesJ )))
      jaccards = overlap / length(unique(c(genesI, genesJ) ))
      df[i,j] = 1-jaccards
    }
  }
  
  ## Cluster nodes using dynamic tree cut, for colors
  distMat = as.dist(df)
  dendro = hclust(distMat, method = "average" )
  clust = dynamicTreeCut::cutreeDynamicTree( dendro, minModuleSize = min.sz )
  ## Note: in dynamicTreeCut, cluster 0, is a garbage cluster for things that dont cluster, so we remove it
  
  gs_results$Cluster = clust
  gs_results = gs_results[gs_results$Cluster != 0, ]
  
  ## select gene sets to label for each clusters
  bests = gs_results %>%  
    group_by( Cluster ) %>% 
    top_n(wt = abs(size), n = 1) %>% 
    .$pathway
  ## determine cluster order for plotting
  clust_ords = gs_results %>% 
    group_by( Cluster ) %>% 
    reframe("Average" = NES ) %>% 
    arrange(desc(Average)) %>% 
    .$Cluster %>% 
    unique
  
  gs_results$Cluster = factor(gs_results$Cluster, levels = clust_ords)
  
  gs_results$Label = ""
  gs_results$Label[gs_results$pathway %in% bests ] = gs_results$pathway[gs_results$pathway %in% bests ]
  gs_results$Label = str_remove(gs_results$Label, "GO_")
  gs_results$Label = tolower(gs_results$Label)
  
  g1 = ggplot(gs_results, aes(x = Cluster, y = NES, label = Label )) +
    geom_jitter( aes(color = Cluster,  size = size), alpha = 0.8, height = 0, width = 0.2 ) +
    scale_size_continuous(range = c(0.5,5)) +
    geom_text_repel( force = 2, max.overlaps = Inf) +
    ggtitle(main) + prism()
  
  return(g1)
}

# overrepresentation analysis
Overrepresentation_analysis <- function(genes,
                                        term_list,
                                        p_thresh=0.05,
                                        pAdjust='fdr',
                                        hypothesis='two.sided'){
  
  term_list <- lapply(term_list,unique) # make sure all terms only have the same gene assigned once
  all_genes <- unique(unlist(term_list)) # getting all genes
  attributes(genes) <- NULL
  
  if(any(genes%notin%all_genes)){
    cat("Some genes are not present in term_list ... \n")
    stop()
  }
  
  if(any(duplicated(genes))){
    warning("Duplicate genes found, called unique ... ")
    genes <- unique(genes)
  }
  
  associated_terms <- lapply(term_list,function(y){y[y%in%genes]})
  terms_to_test <- names(associated_terms[unlist(lapply(associated_terms,length))>0])
  
  cat('Calculting p-values with fisher tests ... \n')
  res <- do.call(rbind,pblapply(terms_to_test,function(term){
    
    fg_in  <- length(associated_terms[[term]])
    fg_out <- length(genes)-fg_in
    bg_in  <- length(term_list[[term]])
    bg_out <- length(all_genes)-bg_in
    
    m <- create_fisher_matrix(fg_in,fg_out,bg_in,bg_out)
    
    fisher_res <- fisher.test(m, conf.int = T, alternative = hypothesis)
    
    expected <- (sum(m[,1])*sum(m[1,]))/sum(m)
    
    res <- data.frame(term=term, count=fg_in, expected=expected, fold_change=fg_in/expected ,pval=fisher_res$p.value, fisher_odds_ratio=fisher_res$estimate, fisher_conf_int_min=fisher_res$conf.int[1], fisher_conf_int_max=fisher_res$conf.int[2])
    rownames(res) <- NULL
    return(res)
  }))
  
  res$`+/-` <- ifelse(log2(res$fold_change)>=0,'+','-')
  res$adjusted <- p.adjust(res$pval,method = pAdjust)
  res$signif <- res$adjusted<p_thresh
  res <- res[,c(1:4,9,5,10,11,6:8)]
  
  res <- res[order(res$`+/-`, -log10(res$adjusted), decreasing = T),]
  
  return(res)
}

create_fisher_matrix <- function(fg_in, 
                                 fg_out, 
                                 bg_in, 
                                 bg_out){
  m <- matrix(c(fg_in,fg_out,bg_in,bg_out), byrow = F, ncol = 2)
  colnames(m) <- c('foreground','background')
  rownames(m) <- c('in','out')
  return(m)
}

# COMPARTMENTS database enrichment analysis
COMPARTMENTS_overrepresentation <- function(proteins, 
                                            background=NULL, 
                                            compartments_threshold=4, 
                                            p_thresh=0.05, 
                                            pAdjust = 'fdr', 
                                            array_object=NULL){
  
  if(is.null(array_object)){
    array_object <- Find_array_object()
  }
  
  background <- get_metadata(array_object = array_object, colnames = 'COMPARTMENTS', proteins = background)
  background <- lapply(background,function(x){
    
  x$GO_term[x$Confidence>=compartments_threshold]
    
  })
  names(background) <- strip_terminal_specifier(names(background))
  background <- background[!duplicated(names(background))]
  
  background <- stack(background)
  background$ind <- as.character(background$ind)
  proteins <- proteins[proteins%in%background$ind]
  background <- split(background$ind,background$values)
  
  Enrich <- Overrepresentation_analysis(genes = proteins, term_list = background, p_thresh = p_thresh, pAdjust = pAdjust)
  
  Enrich$database <- 'COMPARTMENTS'
  colnames(Enrich)[1] <- 'term_label'
  Enrich <- Enrich[,c(12,2,4,7,3,6,5,1,8,9,10,11)]
  Enrich <- Enrich[order(Enrich$pval,decreasing = F),]
  return(as_tibble(Enrich))
}

# Get 2d density
get_density <- function(x, y, ...) {
  dens <- MASS::kde2d(x, y, ...)
  ix <- findInterval(x, dens$x)
  iy <- findInterval(y, dens$y)
  ii <- cbind(ix, iy)
  return(dens$z[ii])
}

# Access netsurf via python
Netsurfp3 <- function(seqs, 
                      names=NULL, 
                      SavePath=getwd(), 
                      filename=NA){
  cat('\n------------------------------------------------------------------------------------------\nIf this is the first time you run this function on your PC it will throw an error ... \nThis is because you need to create a new conda environment; \n> The instructions have been saved to a README.txt file in the folder netsurfp3 in your directory ... \n------------------------------------------------------------------------------------------\n\n')
  
  if(!is.na(filename)){
    dato <- filename
    wd <- getwd()
    setwd(SavePath)
  } else {
    
    if(!is.null(names)==T){
      names(seqs) <- names
      if(length(unique(names))!=length(names)){
        names(seqs) <- paste0(names,'_')
      }
    } else {
      if(length(unique(names(seqs)))!=length(names(seqs))){
        names(seqs) <- paste0(names(seqs),'_')
      }
    }
    
    if(length(attributes(seqs)$names) != length(seqs)){
      cat('please give your sequences names... \n add them as an attribute to your sequences or define them with names=your_names')
      break
    }
    
    wd <- getwd()
    setwd(SavePath)
    
    if(!'./netsurfp3' %in% list.dirs(full.names = T, recursive = F)){
      cat('Creating netsurfp3 directories and files ... \n')
      dir.create('netsurfp3/')
      dir.create('netsurfp3/output/')
      dir.create('netsurfp3/input/')
      
      write_lines(c("# -*- coding: utf-8 -*-","import os","import biolib",
                    "import sys","folder = sys.argv[1]","filename = sys.argv[2]","outdir = sys.argv[3]",
                    "os.chdir(folder)","#print(biolib.__version__)","nsp3 = biolib.load('DTU/NetSurfP_3')",
                    "nsp3_results = nsp3.cli(args='-i ' + filename)","nsp3_results.save_files(outdir)"),
                  file='netsurfp3/netsurfp3.py')
      
      write_lines('conda create --name netsurfp3 \npip3 install -qU pybiolib\n\nIf you want to log in to BioLib while running the function create a BioLib account and add an API token to your environment variables named BIOLIB_TOKEN;\n see https://biolib.com/docs/using-applications/python/\n', 
                  file = 'netsurfp3/README.txt')
      
    }
    cat('Preparing sequences ... \n')
    long <- seqs[nchar(seqs)>1022]
    grps <- c(1,diff(rep(1:40, length.out = length(long))))-1
    grps[which(grps < 0)] <- 1 
    grps <- cumsum(grps)
    long <- split(long,grps)
    
    seqs <- seqs[nchar(seqs)<1022]
    grps <- c(1,diff(rep(1:5000, length.out = length(seqs))))-1
    grps[which(grps < 0)] <- 1 
    grps <- cumsum(grps)
    seqs <- split(seqs,grps)
    
    if(length(long[[1]])>0){
      seqs <- c(seqs,long)
    }
    
    cat('Creating sub-directories and writing fasta files ... \n')
    dato <- gsub(' ','_',gsub(':','',date()))
    dir.create(paste0('netsurfp3/input/',dato,'/'))
    dir.create(paste0('netsurfp3/output/',dato,'/'))
    
    for(i in seq_along(seqs)){
      seqinr::write.fasta(as.list(seqs[[i]]),
                  names = names(seqs[[i]]), 
                  file.out = paste0('netsurfp3/input/',dato,'/',i,'.fasta'),
                  as.string = T)
    }
    cat('Running netsurfp3 on BioLib server via python (this takes a while) ... \n')
    arg1 <- paste0('"',getwd(),'/netsurfp3"')
    for(i in seq_along(seqs)){
      arg2 <- paste0('input/',dato,'/',i,'.fasta')
      arg3 <- paste0('output/',dato,'/output_',i,'/')
      system(paste0('conda run -n netsurfp3 python netsurfp3/netsurfp3.py ',arg1,' ',arg2,' ',arg3), ignore.stdout = T)
    }
  }
  
  cat('Reading netsurfp3 output files ... \n')
  res <- do.call(rbind, 
                 lapply(list.dirs(paste0('netsurfp3/output/',dato,'/'),full.names=F,recursive=F),
                        function(x){
                          read.csv(paste0('netsurfp3/output/',dato,'/',x,'/results.csv'))
                        }
                 )
  )
  colnames(res)[7:9] <- c('helix','strand','coil')
  colnames(res)[11:18] <- c('3_10_helix','alpha_helix','pi_helix','beta_bridge','beta_strand','bend','H_bonded_turn','other')
  res$id <- gsub('>','',res$id)
  setwd(wd)
  cat('\nDone!\n')
  return(res)
  
}

# Access DeepTMHMM via python
DeepTMHMM <- function(seqs,
                      names=NULL, 
                      return_3line=T, 
                      return_gff=F, 
                      SavePath=getwd(), 
                      folder_name=NA, 
                      seq_limit=1000, 
                      add_sleep=120){
  cat('\n------------------------------------------------------------------------------------------\nIf this is the first time you run this function on your PC it will throw an error ... \nThis is because you need to create a new conda environment; \n> The instructions have been saved to a README.txt file in the folder DeepTMHMM in your directory ... \n------------------------------------------------------------------------------------------\n\n')
  
  if(!is.na(folder_name)){
    dato <- folder_name
    wd <- getwd()
    setwd(SavePath)
  } else {
    
    if(!is.null(names)==T){
      names(seqs) <- names
      if(length(unique(names))!=length(names)){
        names(seqs) <- make.unique(names,sep = '-')
      }
    } else {
      if(length(unique(names(seqs)))!=length(names(seqs))){
        cat('Names are not unique, creating unique IDs: Names are made made unique by adding a "-" separator followed by a number\n')
        names(seqs) <- make.unique(names(seqs),sep = '-')
      }
    }
    
    if(length(attributes(seqs)$names) != length(seqs)){
      cat('please give your sequences names... \n add them as an attribute to your sequences or define them with names=your_names\n')
      break
    }
    
    wd <- getwd()
    setwd(SavePath)
    
    if(!'./DeepTMHMM' %in% list.dirs(full.names = T, recursive = F)){
      cat('Creating DeepTMHMM directories and files ... \n')
      dir.create('DeepTMHMM/')
      dir.create('DeepTMHMM/output/')
      dir.create('DeepTMHMM/input/')
      
      write_lines(c("# -*- coding: utf-8 -*-","import os","import biolib",
                    "import sys","folder = sys.argv[1]","filename = sys.argv[2]","outdir = sys.argv[3]",
                    "os.chdir(folder)","#print(biolib.__version__)","deeptmhmm = biolib.load('DTU/DeepTMHMM')",
                    "deeptmhmm_results = deeptmhmm.cli(args='--fasta ' + filename)","deeptmhmm_results.save_files(outdir)"),
                  file='DeepTMHMM/DeepTMHMM.py')
      
      write_lines('conda create --name DeepTMHMM \nconda activate DeepTMHMM \npip3 install -qU pybiolib\n\nIf you want to log in to BioLib while running the function create a BioLib account and add an API token to your environment variables named BIOLIB_TOKEN;\n see https://biolib.com/docs/using-applications/python/\n', 
                  file = 'DeepTMHMM/README.txt')
      
    }
    # long <- seqs[nchar(seqs)>1022]
    # grps <- c(1,diff(rep(1:40, length.out = length(long))))-1
    # grps[which(grps < 0)] <- 1 
    # grps <- cumsum(grps)
    # long <- split(long,grps)
    
    # seqs <- seqs[nchar(seqs)<1022]
    # grps <- c(1,diff(rep(1:5000, length.out = length(seqs))))-1
    # grps[which(grps < 0)] <- 1 
    # grps <- cumsum(grps)
    # seqs <- split(seqs,grps)
    
    cat('Preparing sequences ... \n')
    grps <- c(1,diff(rep(1:seq_limit, length.out = length(seqs))))-1 #DeepTMHMM has a limit 
    grps[which(grps < 0)] <- 1 
    grps <- cumsum(grps)
    seqs <- split(seqs,grps)
    
    # if(length(long[[1]])>0){
    #   seqs <- c(seqs,long)
    # }
    
    cat('Creating sub-directories and writing fasta files ... \n')
    dato <- gsub(' ','_',gsub(':','',date()))
    dir.create(paste0('DeepTMHMM/input/',dato,'/'))
    dir.create(paste0('DeepTMHMM/output/',dato,'/'))
    
    for(i in seq_along(seqs)){
      seqinr::write.fasta(as.list(seqs[[i]]),
                  names = names(seqs[[i]]), 
                  file.out = paste0('DeepTMHMM/input/',dato,'/',i,'.fasta'),
                  as.string = T)
    }
    cat('Running DeepTMHMM on BioLib server via python (this takes a while) ... \n')
    arg1 <- paste0('"',getwd(),'/DeepTMHMM"')
    for(i in seq_along(seqs)){
      arg2 <- paste0('input/',dato,'/',i,'.fasta')
      arg3 <- paste0('output/',dato,'/output_',i,'/')
      system(paste0('conda run -n DeepTMHMM python DeepTMHMM/DeepTMHMM.py ',arg1,' ',arg2,' ',arg3))
      Sys.sleep(rnorm(1,add_sleep,add_sleep/2))
    }
  }
  cat('Reading DeepTMHMM output files ... \n')
  if(return_3line==T){
    res <- read_DeepTMHMM_3line(path = paste0('DeepTMHMM/output/',dato,'/'))
  }
  
  if(return_gff==T){
    res <- read_DeepTMHMM_3line(path = paste0('DeepTMHMM/output/',dato,'/'))
  }

  if(return_gff==F & return_3line==F){
    cat(paste0('Files will not be returned, but files are always stored in your current directory:\nDeepTMHMM/output/',dato,'/\n'))
    res <- NA
  }
  
  setwd(wd)
  cat('\nDone!\n')
  return(res)
  
}

# Functions to read netsurf and DeepTMHMM output
read_DeepTMHMM_3line <- function(path=NULL){
  
  if(is.null(path)){
    cat('No path provided ...\n')
    stop()
  }
  
  res <- unlist(lapply(list.dirs(path,full.names=F,recursive=F),
         function(x){
           out_files <- seqinr::read.fasta(file=paste0(path,x,'/predicted_topologies.3line'), seqtype = 'AA')
           lapply(out_files,function(y){
             data.frame(aa=y[1:(length(y)/2)], topology=y[(1+(length(y)/2)):length(y)])
           })
         }),recursive = F)
  attributes(res)$keys <- data.frame(code=c('S','I','O','M'), note=c('Signalpeptide', 'Intracellular', 'Extracellular','In Membrane'))
  return(res)
}

read_DeepTMHMM_gff <- function(path=NULL){
  
  if(is.null(path)){
    cat('No path provided ...\n')
    stop()
  }
  
  res <- do.call(rbind,lapply(list.dirs(path,full.names=F,recursive=F), function(x){
                         out_files <- na.omit(read.delim(paste0(path,x,'/TMRs.gff3'), header=F, comment.char="#", na.strings = '//')[,1:4])
                       }))
  return(res)
}

list_DeepTMHMM_folders_in_directory <- function(path=getwd()){
  
  wd <- getwd()
  setwd(path)
  
  if("./DeepTMHMM"%in%list.dirs(recursive = F)){
    
    res <- list.dirs(path = './DeepTMHMM/output', recursive = F, full.names = F)
    
  } else {
    cat("You have not run DeepTMHMM() in this directory please specify path to the desired DeepTMHMM folder ... \n")
    stop()
  }
  
  setwd(wd)
  return(res)
}

list_netsurfp3_folders_in_directory <- function(path=getwd()){
  
  wd <- getwd()
  setwd(path)
  
  if("./netsurfp3"%in%list.dirs(recursive = F)){
    
    res <- list.dirs(path = './netsurfp3/output', recursive = F, full.names = F)
    
  } else {
    cat("You have not run Netsurfp3() in this directory please specify path to the desired netsurfp3 folder ... \n")
    stop()
  }
  
  setwd(wd)
  return(res)
}

# Functions for reading alphafold .cif files
read_cif_structure <- function(cif){
  
  if(str_detect(cif,'.cif')){
    
    cif <- readLines(cif) #read whole file
    
    for(i in c('_struct_','loop_')){ #discard lines surrounding structure
      keep <- which(str_detect(cif,i))
      cif <- cif[c(min(keep):max(keep))]
    }
    
    cif <- cif[-which(str_detect(cif,'#|loop_'))] #remove remaining symbols
    header <- cif[which(str_detect(cif,'_struct_conf.'))] #extract headers
    cif <- cif[c((length(header)+1):length(cif))] #removing headers
    
    cif <- gsub('\\? ',' \\?',cif) #removing spaces
    cif <- as.data.frame(cif) #converting to dataframe 
    cif <- cif %>% separate(cif, header,"\\s+") #separating columns
    cif <- cif[,c(5:7,12:13)]
    colnames(cif) <- c('start_aa','start','structure','end_aa','end')
    cif$start <- as.numeric(cif$start)
    cif$end <- as.numeric(cif$end)
    cif
  } else {
    stop('file is not of type .cif')
  }
}
read_cif_confidence <- function(cif){
  
  if(str_detect(cif,'.cif')){
    
    cif <- readLines(cif) #read whole file
    
    #discard lines surrounding structure
    k1 <- which(str_detect(cif,'_ma_qa_metric_local'))
    k2 <- which(str_detect(cif,'_ma_qa_metric_local|#'))
    keep <- k2[c(min(which(k2 %in% k1))-1,max(which(k2 %in% k1))+1)]
    
    cif <- cif[c(keep[1]:keep[2])]
    
    cif <- cif[-which(str_detect(cif,'#|loop_'))] #remove remaining symbols
    header <- cif[which(str_detect(cif,'_ma_qa_metric_local.'))] #extract headers
    cif <- cif[c((length(header)+1):length(cif))] #removing headers
    
    cif <- str_trim(cif) #removing spaces at end
    cif <- as.data.frame(cif) #converting to dataframe 
    cif <- cif %>% separate(cif, header,"\\s+") #separating columns
    cif <- cif[,c(2:3,5)]
    colnames(cif) <- c('aa','residue','confidence')
    cif$residue <- as.numeric(cif$residue)
    cif$confidence <- as.numeric(cif$confidence)
    cif
  } else {
    stop('file is not of type .cif')
  }
}

# Returns amino acid properties
AA_Properties <- function(pH=7.4){
  
  AA_Properties <- tibble(AA= c("A", "R", "N", "D", "C", "U" ,"Q", "E",
                                "G", "H", "I", "L", "K", "M", "F",
                                "P", "S", "T", "W", "Y", "V"), 
                          Hydrophobicity = c(0.310,-1.010,-0.600,-0.770,
                                             1.540,1.540,-0.220,-0.640,0.000,
                                             0.130,1.800,1.700,-0.990,
                                             1.230,1.790,0.720,-0.040,
                                             0.260,2.250,0.960,1.220),
                          z.cont = c(0/(10^(14-pH)+1), 1/(10^(pH-12.48)+1),
                                     0,-1/(10^(3.65-pH)+1),-1/(10^(8.18-pH)+1),
                                     -1/(10^(8.18-pH)+1),
                                     0,-1/(10^(4.25-pH)+1),0,1/(10^(pH-6)+1),
                                     0,0,1/(10^(pH-10.53)+1),-1/(10^(14-pH)+1),
                                     0,0,-1/(10^(14-pH)+1),-1/(10^(14-pH)+1),
                                     0,-1/(10^(10.07-pH)+1),0), 
                          z = c(0,1,0,-1,0,0,0,-1,0,0,0,0,1,0,0,0,0,0,0,0,0))
  AA_Properties <- AA_Properties[order(AA_Properties$AA),]
  
  AA_Properties$categories <- as.factor(c(2,2,3,3,2,5,1,2,1,2,2,4,5,4,1,4,4,5,2,2,2))
  levels(AA_Properties$categories) <- c('Positively Charged','Hydrophobic','Negatively Charged','Polar','Special Cases')
  
  # Set colours
  AA_Properties$color <- as.factor(c(1,2,3,3,2,1,4,2,5,2,2,6,7,6,5,8,8,9,2,2,2))
  levels(AA_Properties$color) <- c('grey','#ffdb58','red','turquoise','blue','rosybrown1','forestgreen','purple','orange')
  
  return(AA_Properties)
}

# Get amino acids by 1 letter code or property
Amino_acids <- function(exclude=NULL, 
                        group=NULL){
  if(!is.null(group)){
    AA <- AA_Properties()
    AA <- AA$AA[AA$categories%in%group]
  } else {
    AA <- AA_Properties()$AA
  }
    AA <- AA[AA%notin%exclude]
  
  return(AA)
}

# Returns legend from a ggplot2 object
grab_legend <- function(a.gplot){ 
  tmp <- ggplot_gtable(ggplot_build(a.gplot)) 
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box") 
  legend <- tmp$grobs[[leg]] 
  legend
} 

# Functions for using ggmsa package with msa package
make_msa_object_ggmsa_compatible <- function(msa_alignment){
  seqs <- msaConvert(msa_alignment,"seqinr::alignment")$seq
  aln <- ape::as.AAbin(split(seqs, 1:length(seqs)))
  names(aln) <- rownames(msa_alignment)
  return(aln)
  
}

# Calculate dist matrix from msa
msa_dist_matrix <- function(msa_alignment){
  dist.alignment(msaConvert(msa_alignment,type="seqinr::alignment"),"identity") 
}

# Calculate mismatches in alignment
mismatches <- function(ref,alignment,thresh=0.5, remove=T){
  mismatched <- c()
  for(seq in 2:nrow(alignment)){
    d <- unlist(str_split(diag(attr(
      adist(ref,as.character(AAStringSet(alignment)[seq]),counts=T),
      "trafos")),
      ''))
    if((length(d[d!='M'])/length(d))>thresh & remove==T){
      mismatched <- c(mismatched,seq)
    }
    if(remove==F){
      mismatched <- c(mismatched,(length(d[d!='M'])/length(d)))
    }
  }
  return(mismatched)
}

# Access PaCMAP algorithm built in python
reticulate::use_condaenv("pacmap", required = TRUE)
pacmap <- function(x){
  pacmap <- reticulate::import("pacmap")

  # Initialize PaCMAP instance
  reducer <- pacmap$PaCMAP()
  
  # Perform dimensionality Reduction
  embedding <- reducer$fit_transform(x)
  
  as.data.frame(embedding)
}

# Add list names to listed dataframes as a new column
add_list_names_as_column <- function(list, 
                                     col_name){
  for(i in names(list)){
    list[[i]][,col_name] <- i
  }
  return(list)
}

# Silence a function
quiet <- function(x) { 
  sink(tempfile()) 
  on.exit(sink()) 
  invisible(force(x)) 
} 

# Move elements to end of list
move_element_to_bottom_of_list <- function(list, 
                                           to_move_to_bottom){
  
  if(class(list)=='array_object'){
    cs <- class(list)
  } else {
    cs <- class(list)
  }
  
  if(is.null(names(list))){
    
    if(!is.numeric(to_move_to_bottom)){
      cat('You have not named the list ... \nIf you want to move by index provide a numeric vector.\n')
      stop()
    } else {
      
      l <- seq_along(list)
      list <- list[c(l[!l%in%to_move_to_bottom], to_move_to_bottom)]
      
    }
  } else {
    
    if(is.numeric(to_move_to_bottom)){
      l <- seq_along(list)
      list <- list[c(l[!l%in%to_move_to_bottom], to_move_to_bottom)]
    } else {
      
      to_move_to_bottom <- sapply(to_move_to_bottom, function(x){which(names(list)==x)})
      l <- seq_along(list)
      list <- list[c(l[!l%in%to_move_to_bottom], to_move_to_bottom)]
    }
  }
  class(list) <- cs
  return(list)
}

# Convert all capital three letter AA codes to single letter
convert_aa <- function(aa){
  
  converted <- sapply(aa,function(x){
  if(nchar(x)==3){
    res <- a(toTitleCase(tolower(x)))
  }
  
  if(nchar(x)==1){
    res <- aaa(toupper(x))
  }
  
  if(nchar(x)%notin%c(1,3)){
    print('The provided string is not an Amino acid code')
    res <- NA
    print(aa)
  }
  return(res)
  })
  return(converted)
}

# Remove all variables except Array_object or other specific variable
keep_rm <- function(keep = Filter(Array_Class_Filter, ls(name = .my_env))){
  cat('\n---------------- Removing unwanted objects ----------------\n')
  message('Removing all objects')
  keep <- unique(c(keep,Filter(Array_Class_Filter, ls(name = .my_env))))
  message(paste('Exceptions:',keep, sep = ' '))
  rm(list=setdiff(ls(name = .my_env), keep), envir = .my_env)
  message('Clearing Memory')
  gc(verbose = F)
  cat('\n------------------- Importing functions -------------------\n')
  Array_object$Import_functions()
  cat('-----------------------------------------------------------\n')
}

message('Done')

# Print info about array_object ------------------------------------------------
about_array_object()