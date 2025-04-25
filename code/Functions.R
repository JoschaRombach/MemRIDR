# Functions for HD-peptide Array shiny app
# Joscha Rombach
# 24/03/2023

# Note - SWaFi stands for Smoothed Weighted-average Fluorescence intensity,
# which is referred to as "binding score" in the publication

#-############################################################################-#
#-########################### INSTALL PACKAGES ###############################-#
#-############################################################################-#

# Install packages from CRAN, Bioconductor and Github
# cran <- c('tidyverse','berryFunctions','seqinr','stringr','data.table',
#           'viridis','bio3d','ggprism','ggExtra','gridExtra','grid','zoo','RColorBrewer',
#           'shiny','DT','shiny.fluent','imola','cowplot','shiny.router','shiny.react',
#           'shinyWidgets','circlize', 'visNetwork','igraph','shinycssloaders',
#           'plyr','ggnewscale','googledrive')
# 
# bioconductor <- c('msa','ggtree')
# 
# github_packages <- c('r3dmol','dqshiny')
# github_pages    <- c('swsoyee/r3dmol','daqana/dqshiny')
# 
# # Function for installing missing packages from CRAN and Bioconductor
# install_CRAN_and_Bioconductor_packages <- function(cran,bioconductor){
#   
#   cran_missing <- cran[which(!cran %in% rownames(installed.packages()))]
#   bioconductor_missing <- bioconductor[which(!bioconductor %in% rownames(installed.packages()))]
#   n_missing <- length(cran_missing)+length(bioconductor_missing)
#   message("\n-----------------------------------------------------------")
#   message(paste0(' Installing ',n_missing,' packages from CRAN and Bioconductor'))
#   
#     if( n_missing > 0){
#     
#     continue <- NA
#     while(!continue%in%c('y','Y','ye','YE','Ye','yes','YES','Yes','YeS','YEs','yeS','yES','n','N','no','No','NO')){
#       continue <- readline(paste0(" There are ",n_missing," packages that need to be installed - Do you want to install them? options: (yes/no) - "))
#     }
#     
#     if(continue%in%c('y','Y','ye','YE','Ye','yes','YES','Yes','YeS','YEs','yeS','yES')){
#       
#       if(length(cran_missing)>0){
#         install.packages(cran_missing)
#       }
#       
#       if(length(bioconductor_missing)>0){
#         if (!require("BiocManager", quietly = TRUE))
#           install.packages("BiocManager")
#         
#         BiocManager::install(bioconductor_missing)
#       }
#     }
#     
#     cat('\n')
# 
#   } else {
#     message(' All packages already installed')
#   }
# }
# 
# # Install from CRAN and Bioconductor
# install_CRAN_and_Bioconductor_packages(cran,bioconductor)
# 
# # Install r3dmol from Github
# 
# install_github_packages <- function(github_packages, github_pages){
# 
#     missing <- which(!github_packages %in% rownames(installed.packages()))
#     
#     if(length(missing)>0){
#       
#     continue <- NA
#     while(!continue %in% c('yes','YES','Yes','YeS','YEs','yeS','no','No','NO')){
#       continue <- readline(paste0(' Install ',paste(github_packages[missing], collapse = ' & '),' from pages ',paste(github_pages[missing], collapse = ' & '),'? note: remotes is required and will also be installed. (yes/no) : '))
#     }
#     
#     if(continue%in%c('yes','YES','Yes','YeS','YEs','yeS')){
#       if(!require('remotes')){install.packages('remotes')}
#       remotes::install_github(repo = github_pages[missing], upgrade = 'never')
#     }
#     
#     rm(continue)
#   
#     } else {
#       message(' All packages already installed')
#     }
#    }
# 
# install_github_packages(github_packages, github_pages)
# 
# #-############################################################################-#
# #-########################### LOAD PACKAGES ##################################-#
# #-############################################################################-#
# 
# # Function for loading packages
# load_packages <- function(packages){
#   message(" Loading packages ...")
#   for (p in packages) {
#     if (p %in% rownames(installed.packages())) {
#       suppressPackageStartupMessages(library(p, character.only=TRUE))
#     } else {
#       message(paste0(' Error: package ',p,' is not installed ... '))
#     }
#   }
# }
# 
# # Load packages
# load_packages(c(cran, bioconductor, github_packages))
# 
# rm(cran, 
#    bioconductor, 
#    github_packages,
#    github_pages,
#    install_github_packages,
#    install_CRAN_and_Bioconductor_packages, 
#    load_packages)

#-############################################################################-#
#-############################# AESTHETICS ###################################-#
#-############################################################################-#

# Set ggplot theme for script
message(" Setting ggplot theme to prism style ...")
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

#-############################################################################-#
#-########################### SHINY APP FUNCTIONS ############################-#
#-############################################################################-#
message(" Loading functions ...")
search_for_protein <- function(search, keys){
  
  if(search==''){
    search <- 'Q01814'
  }
  
  res <- keys[str_detect(keys$Search,tolower(search)),]
  
  if(nrow(res)==0){
    res <- keys[keys$UniProtKB.AC=='Q01814',]
  }
  colnames(res) <- c('Name','Genes','UniProtKB.AC','Search','Path')
  
  return(res)
}

search_for_disease <- function(search, keys, ClinVar){

  if(search==''){
    search <- 'epilepsy'
  }

  if(tolower(search) %in% ClinVar$Phenotype){
    res <- ClinVar[ClinVar$Phenotype==tolower(search),]
  } else {
    res <- ClinVar[str_detect(ClinVar$Phenotype,tolower(search)),]
  }

  if(nrow(res)==0){
    res <- ClinVar[ClinVar$Phenotype=='epilepsy',]
  }
  
  res <- merge(x=res, y=keys, by.x='UniProt', by.y='UniProtKB.AC', all.x=T)
  res <- res[,c('Protein.name','Genes','UniProt','Phenotype','Path')]
  colnames(res) <- c('Name','Genes','UniProtKB.AC','Search','Path')
  res <- unique(res)
  
  return(res)
}

#-############################################################################-#
#-########################## GOOGLE DRIVE FUNCTIONS ##########################-#
#-############################################################################-#

Authenticate_google_drive <- function(){
  googledrive::drive_auth(path = 'data/winter-environs-456814-d2-f014d23cacc8.json')
  print(drive_user())
  cat('\nCreating temporary directory\n')
  if(!'temp_dir' %in% list.files()){dir.create('temp_dir')} # test if temporary directory exists
}

get_array_file_google_drive <- function(array_object_file, files){
  
  wd <- getwd() # get directory path
  #if(!'temp_dir' %in% list.files()){dir.create('temp_dir')} # test if temporary directory exists
  setwd('temp_dir') # navigate to temporary directory
  
  if(!array_object_file %in% list.files(recursive = F)){ # test if file has been downloaded 
    
    googledrive::drive_download(file = array_object_file, 
                                overwrite = F)
    
    # fetching file from googledrive
    # googledrive::drive_download(file = as_id(paste0("https://drive.google.com/open?id=",
    #                                                 files$id[files$name == array_object_file])), 
    #                             overwrite = F)
  } 
  
  out <- readRDS(array_object_file) # read file
  
  setwd(wd)
  return(out)
}

#-############################################################################-#
#-########################### CONVENIENCE FUNCTIONS ##########################-#
#-############################################################################-#

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

#-############################################################################-#
#-############## GET FUNCTIONS ('get' items from array_object) ###############-#
#-############################################################################-#

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

get_available_proteins <- function(array_object=NULL){
  unique(strip_terminal_specifier(names(quiet(get_terminals(array_object)))))
}

get_terminals <- function(array_object=NULL, 
                          proteins=NULL){
  
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

get_SWaFi_traces <- function(array_object=NULL,
                             dataset=NULL,
                             proteins=NULL,
                             parameters=F){

  if(!is.null(array_object)){

  experiment <- get_experiment(array_object, dataset)

  if('SWaFi_traces' %notin% names(experiment)){
    print(paste0("Error, SWaFi traces do not exist in experiment '",dataset,"' ... please add traces to continue"))
    stop()
  }

  traces <- experiment$SWaFi_traces
  if(!is.null(proteins)){
    id_type <- detect_identifier(proteins)

    if(id_type=="UniProtKB.AC"){
      #print('Returning all available terminals for the specified proteins...')
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
}

get_averaged_SWaFi_traces <- function(proteins = NULL,
                                      datasets = c('FBBE1','FBBE2'), 
                                      array_object=NULL){
  
  if(!is.null(array_object)){
  
  scores <- unlist(lapply(datasets, function(dataset){
    
    s <- quiet(get_SWaFi_traces(array_object = array_object, 
                                dataset  = dataset, 
                                proteins = proteins))
    
    if(dataset==datasets[1]){
      
      s <- lapply(s, function(x){ 
        x$score <- x$score*1.91549
        return(x)
      })
    }
    
    return(s)
  }), recursive=F)
  
  scores <- split(scores, unique(names(scores)))
  
  scores <- lapply(scores, function(s){
    
    if(!identical(s[[1]]$position,s[[2]]$position)){
      message('Positions datasets do not match! Please double check how datasets were imported!')
      stop()
    }
    
    data.frame(score = Reduce('+',s)$score, 
               position = s[[1]]$position)
  })
  
  return(scores)
}
}

# this is the final function to get binding sites
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
  if(!is.null(array_object)){

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
} 

#-############################################################################-#
#-############# PLOT FUNCTIONS ('plot' data from array_object) ###############-#
#-############################################################################-#

plot_SWaFi_trace_w_sites <- function(array_object=NULL, 
                                     traces = NULL,
                                     proteins=NULL,
                                     sites=NULL,
                                     ymax=NULL, 
                                     confidence='medium', 
                                     trace_color = '#0275d8',
                                     highlight_color = 'red'){
  
  # Find array object if missing
  if(!is.null(array_object)){
  
  if(proteins %notin% c(names(array_object$Protein_terminals), 
                        array_object$Protein_terminals[[1]]$UniProtKB.AC)){
    
    p <- ggplot() + 
      prism() + 
      theme(aspect.ratio = NULL, legend.position = 'none', panel.background = element_blank()) +
      ylab('Binding score (a.u.)') + 
      xlab('Position') 
    
    return(p)
    
  } else {
  
  #traces <- quiet(get_averaged_SWaFi_traces(proteins, dataset, array_object))
  
  ids <- rev(names(traces))
  
  ratio <- rev(unlist(lapply(traces,nrow)))
  
  if(is.null(ymax)){
    ymax <- max(do.call(rbind,traces)$score)
  }
  
  # if(is.null(sites)){
  #   sites <- quiet(get_binding_sites(proteins = proteins, 
  #                                    binding = T, 
  #                                    datasets = dataset, 
  #                                    confidence = confidence, 
  #                                    array_object = array_object))
  #   
  # }
  
  trace_plots <- lapply(ids,function(id){
    trace <- traces[[id]]
    
    #trace$id <- '0'
    max.x <- ifelse(nrow(trace)>75,round(max(trace$position),-2),round(max(trace$position),-1))
    min.x <- ifelse(nrow(trace)>75,round(min(trace$position),-2),round(min(trace$position),-1))
    max.x <- ifelse(max.x<max(trace$position),max(trace$position),max.x)
    min.x <- ifelse(min.x>min(trace$position),min(trace$position),min.x)
    
    #sites <- sites[sites$UniProtKB.AC==id,]
    
    p <- ggplot()
    
    sites <- sites[!is.na(sites$start),]
    
    if(nrow(sites)>0){
      p <- p + geom_rect(data = sites, aes(xmin=start, xmax=end, ymin=rep(-Inf,nrow(sites)), ymax=rep(Inf,nrow(sites)), fill=1:nrow(sites)), alpha=0.15)
    } 
    
    p <- p + 
      geom_line(data = trace,  aes(x=position, y=score), col = trace_color) +
      geom_point(data = trace, aes(x=position, y=score), pch=16, col = trace_color) +
      ylim(0,ceiling(ymax)) +
      xlim(min.x,max.x) +
      prism() + 
      theme(aspect.ratio = NULL, legend.position = 'none', panel.background = element_blank()) +
      ylab('Binding score (a.u.)') + 
      xlab('Position') 
    
    return(p)
  })
  
  if(length(ids)==2){
    if(nrow(sites)==1){
      trace_plots[[2]] <- trace_plots[[2]] + scale_fill_gradient(low =  highlight_color, high = highlight_color) + remove_y_axis()
      trace_plots[[1]] <- trace_plots[[1]] + scale_fill_gradient(low =  highlight_color, high = highlight_color)
    } else {
      trace_plots[[2]] <- trace_plots[[2]] + scale_fill_viridis(direction = -1) + remove_y_axis()
      trace_plots[[1]] <- trace_plots[[1]] + scale_fill_viridis()
    }
  } else {
    if(nrow(sites)==1){
      trace_plots[[1]] <- trace_plots[[1]] + scale_fill_gradient(low =  highlight_color, high = highlight_color)
    } else {
      trace_plots[[1]] <- trace_plots[[1]] + scale_fill_viridis()
    }
  }
  
  if(length(ids)==1){
    ratio <- 1
  }
  
  cowplot::plot_grid(plotlist = trace_plots, rel_widths = ratio, align = 'hv', nrow = 1) 
   }
 }
}

plot_list_of_ggplots <- function(plot_list,ncol=1){
  do.call("grid.arrange", c(plot_list, ncol=ncol))
}

plot_pdb_with_binding_sites <- function(UniProtKB.AC,
                                        SWaFi_score=NULL,
                                        pthresh=0.05,
                                        protein_color="#d3d3d3",
                                        protein_opacity=0.5,
                                        site_color='red',
                                        AlphaFold_directory=NULL,
                                        plot_all_scores=T,
                                        array_object=NULL,
                                        dataset=c('FBBE1','FBBE2'),
                                        zoom_limits=c(10,500),
                                        quality=5,
                                        colorRamp2_function=NULL,
                                        brewerpal_color='YlGnBu',
                                        non_binding_color='lavender',
                                        cartoon_style='oval',
                                        background_color='white', 
                                        fontColor='#0dc5c1'){
  if(!is.null(array_object)){

  Exists <- UniProtKB.AC %in% c(names(array_object$Protein_terminals), 
                                strip_terminal_specifier(names(array_object$Protein_terminals)))

  pdb <- quiet(get_metadata(array_object = array_object, 
                            colnames = 'AlphaFold_pdb', 
                            proteins = strip_terminal_specifier(UniProtKB.AC))[[1]])
  
  # Test if structure has been found
  if(length(pdb)==1){
    #message(paste0('array_object does not contain a structure for this protein - ',UniProtKB.AC,' ... \n'))
    res <- r3dmol(viewer_spec = m_viewer_spec(
                  backgroundColor = background_color,
                  cartoonQuality  = quality)) %>% 
           m_add_label(text = 'Structure not found', 
                       style = m_style_label(fixed = T,
                                             alignment = 'center', 
                                             fontColor = fontColor, 
                                             showBackground = F))
  } else {
  
  # load bio3d object into r3dmol 
  res <- r3dmol(viewer_spec = m_viewer_spec(
                backgroundColor = background_color,
                cartoonQuality  = quality,
                lowerZoomLimit  = zoom_limits[1],
                upperZoomLimit  = zoom_limits[2])) %>% 
         m_add_model(data=m_bio3d(pdb)) %>% 
         m_zoom_to() %>%  
         m_set_style(style = m_style_cartoon(color = protein_color, 
                                             opacity = protein_opacity, 
                                             style=cartoon_style)) 
  
  if(plot_all_scores==T & Exists){
    
    # Define color scale 
    if(is.null(colorRamp2_function)){
      pos <- brewer.pal(5,brewerpal_color)
      colorRamp2_function <- colorRamp2(breaks = c(0,0.1,0.2,0.3,0.4,0.5,0.6), 
                                        colors = c(non_binding_color,'grey40',pos))
      
    }
  
    scores <- get_averaged_SWaFi_traces(proteins = UniProtKB.AC,
                                        datasets = dataset, 
                                        array_object = array_object)
    
    scores <- lapply(scores,function(x){
                     x$color <- colorRamp2_function(x$score)
                     x$color <- substring(x$color,1,7) # for some reason colorRamp adds two F's to the end of each color??? does not work with r3dmol so we remove them
                     data.frame(position=riffle(x$position,x$position+1),
                                color   =riffle(x$color,x$color))
                     })
    
    len <- nchar(unlist(quiet(get_metadata(array_object = array_object, 
                                           colnames = "sub.sequence", 
                                           proteins = names(scores)))))
    
    scores <- na.omit(do.call(rbind,lapply(seq_along(len), function(i){
                    scores[[i]] <- scores[[i]][1:len[i],]
                    })))
    
    for(i in 1:nrow(scores)){
      res <- res %>% m_set_style(style = m_style_cartoon(color = scores$color[i]), 
                                 sel = m_sel(resi = c(scores$position[i],scores$position[i])))
    }
   } 
  }
  return(res)
  }
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

plot_vis_network <- function(UniProtKB.AC, graph=NULL, 
                             disease_col = "#800074",
                             protein_col = "#298c8c", 
                             height = NULL, width = NULL, 
                             diseases=NULL){
  
  if(is.null(graph)){
    message('No graph found, stopping ...')
    stop()
  }
  
  if(is.null(diseases)){
    id <- names(V(graph))[V(graph)$UniProtKB.AC%in%UniProtKB.AC]
  } else {
    id <- names(V(graph))[names(V(graph))%in%diseases]
  }
  
  if(length(id)<1){
    
    p <- visNetwork(nodes = data.frame(label='no pathological variants found in ClinVar'), width = width, height = height)
    
  } else {
    
    memb   <- components(graph)$membership
    find   <- unlist(lapply(split(names(V(graph)), memb), function(x){ sum(str_detect(x, pattern=paste(id, collapse = '|')))}))
    temp_g <- induced_subgraph(graph, V(graph)[names(V(graph))%in%names(memb[memb %in% names(find[find>0])])])
    
    vis_edges <- as_data_frame(temp_g)
    colnames(vis_edges) <- c("from", "to", "width")
    vis_edges$width = rescale(vis_edges$width, 1,15)
    
    vis_nodes <- data.frame(id    = names(V(temp_g)), 
                            label = names(V(temp_g)), 
                            value = V(temp_g)$Freq, 
                            group = V(temp_g)$type,
                            color = ifelse(V(temp_g)$type=='Protein', protein_col, disease_col),
                            sele  = names(V(temp_g)) %in% id)
    
    p <- visNetwork(nodes = vis_nodes, edges = vis_edges, width = width, height = height) %>%
      
      
      visOptions(nodesIdSelection = list(selected = id[1], 
                                         style = 'width: 0px; height: 0px; opacity: 0'), # hide selection toggles
                        clickToUse = F, 
                        highlightNearest = list(enabled = T, degree = 1, hover = T)) %>%
      
      visPhysics(stabilization = F) %>%
      
      visNodes(
        color = list(
          group = c(disease_col,protein_col)
        )) %>%
      
      visEdges(
        smooth = T,
        physics = T,
        shadow = F,
        color = list(color = alpha("grey",0.7), highlight = alpha('black',0.5))
      ) %>%
      # line to get id of node to etract ClinVar variants
      visEvents(select = "function(nodes) {
                Shiny.onInputChange('current_node_id', nodes.nodes);
                ;}")
  }
  p
}

# Heliquest function
heliquest <- function(sequence,
                      pH=7.4, ...){
  
  heliquest_coords <- fast_heliquest_coords(sequence,pH)

  moment <-  sqrt((sum(heliquest_coords$Hydrophobicity*heliquest_coords$cos)/nrow(heliquest_coords))^2 +
                  (sum(heliquest_coords$Hydrophobicity*heliquest_coords$sin)/nrow(heliquest_coords))^2)
  
  p <- plot_helical_wheel(heliquest_coords, ...)
  
  print(p)

  return(list(coordinates=heliquest_coords,H.moment=moment, plot=p))
  
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

# Plot Heliquest output
plot_helical_wheel <- function(heliquest_coordinates, 
                               correct_rotation=F, 
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
  moment <- sqrt(rotated_vec[1]^2 + rotated_vec[2]^2)
  
  if(density==F){
    fctr <- (0.375*(sort(rep((1:ceiling(nrow(heliquest_coordinates)/18)),18))[1:nrow(heliquest_coordinates)]-1))+1
    
    p <- ggplot(heliquest_coordinates, aes(cos*fctr,sin*fctr)) + 
      geom_path(aes(col=order,linewidth=order)) + 
      scale_linewidth(range=c(0.5,3)) +
      scale_color_gradient(low='grey80',high='grey30') + 
      ggnewscale::new_scale_color() +
      # geom_segment(aes(x=start, y=start, xend=cos*Hydrophobicity ,yend=sin*Hydrophobicity),
      #                  arrow = arrow(length = unit(0.5, "cm"),type = 'open'),linetype=2,linewidth=.5) +
      geom_point(aes(col=AA),pch=16,size=20) +
      annotate(geom = 'segment',x=0, y=0, xend=rotated_vec[1] ,yend=rotated_vec[2],
                   arrow = arrow(length = unit(0.5, "cm"),type = 'closed'), linewidth=1.5) +
      geom_text(aes(label=AA)) +
      geom_text2(aes(label=c('N',rep('',nrow(heliquest_coordinates)-2),'C')), nudge_x = -0.16, nudge_y = -0.1,col='red') +
      prism() + remove_x_axis() + remove_y_axis() + theme(legend.position=legend_position) +
      xlim(-1.5,1.5) +
      ylim(-1.5,1.5) +
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
      xlim(-1.5,1.5) +
      ylim(-1.5,1.5) +
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

#-############################################################################-#
#-############### ADD FUNCTIONS ('add' data to array_object) #################-#
#-############################################################################-#

add_helical_confidence_to_SWaFi_binding_sites <-function(SWaFi_sites,
                                                         cassette_len = 11,
                                                         pH = 7.4,
                                                         array_object=NULL,
                                                         seq_col = 'sequence'){
  
  if(is.null(array_object)){

  
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
}

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
#-########################### GENERIC FUNCTIONS ##############################-#
#-############################################################################-#

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
    # if(is.null(array_object)){
    #   array_object <- Find_array_object() 
    # }
  }
  
  # Notify user if the Speed_up has been used, might impact precision ----------
  if(speed_up==T){
    cat(paste0('Speed_up set to TRUE only keeping sequences with the top ',10*(1/speed_up_by),'% highest hydrophobic moment for each length ... \nMight be less accurate, but can increase speed by up to 50%! ...\n'))
  }
  
  # Begin function -------------------------------------------------------------
  
  result <- lapply(seq_along(sequences), function(i){  # Running lapply to vectorize function 
  
    sequence <- sequences[i] # get sequence
    
    if(!is.na(sequence) & nchar(sequence)>(min_len-1)){
    
      if(!is.null(UniProtKB.AC)){
          UniProt <- UniProtKB.AC # get UniProt ID if array_object is used for helicity
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

# Returns amino acid properties
AA_Properties <- function(pH=7.4){
  
  AA_Properties <- data.frame(AA= c("A", "R", "N", "D", "C", "U" ,"Q", "E",
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

# Silence a function
quiet <- function(x) { 
  sink(tempfile()) 
  on.exit(sink()) 
  invisible(force(x)) 
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

# Confidence converter 
conf_convert <- function(conf){
  
  if(conf == 'medium'){
    conf <- 0.05
  }
  
  if(conf == 'low'){
    conf <- 0.1
  }
  
  if(conf == 'high'){
    conf <- 0.01
  }
  return(conf)
}

message(' Done')

# Print info about array_object ------------------------------------------------
message("-----------------------------------------------------------",
    "\n Authors - Joscha Rombach & Tommas T. E. Nielsen &", 
    "\n           Junior Agenant & Kenneth L. Madsen",
    "\n-----------------------------------------------------------",
    "\n Date created:           2022-12-21", 
    "\n Last updated:           Thu Sep 19 18:49:00 2024", 
    "\n Last saved:             Thu Sep 19 21:24:06 2024",
    "\n # of Experiments:       2", 
    "\n Experiment IDs:         FBBE1 FBBE2", 
    "\n Affiliation:            University of Copenhagen, Denmark",
    "\n-----------------------------------------------------------",
    "\n If you use this database please cite:\n Rombach, J., Nielsen, T. T. E.,(2025) .....",
    "\n-----------------------------------------------------------")
