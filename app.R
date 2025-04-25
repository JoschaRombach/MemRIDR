# Author: Joscha Rombach
# Created: 12/04/2024
# Project: Web interface for array_object

# ctrl F and replace to change color
# primary color: steelblue 
# secondary color: red
# old color: #0dc5c1 

#-############################################################################-#
#-############################## DEPENDENCIES ################################-#
#-############################################################################-#
message(" Starting ...")
message(" Loading packages ...")

#CRAN
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(berryFunctions))
suppressPackageStartupMessages(library(seqinr))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(viridis))
suppressPackageStartupMessages(library(bio3d))
suppressPackageStartupMessages(library(ggprism))
suppressPackageStartupMessages(library(ggExtra))
suppressPackageStartupMessages(library(gridExtra))
suppressPackageStartupMessages(library(grid))
suppressPackageStartupMessages(library(zoo))
suppressPackageStartupMessages(library(RColorBrewer))
suppressPackageStartupMessages(library(shiny))
suppressPackageStartupMessages(library(DT))
suppressPackageStartupMessages(library(shiny.fluent))
suppressPackageStartupMessages(library(imola))
suppressPackageStartupMessages(library(cowplot))
suppressPackageStartupMessages(library(shiny.router))
suppressPackageStartupMessages(library(shiny.react))
suppressPackageStartupMessages(library(shinyWidgets))
suppressPackageStartupMessages(library(circlize))
suppressPackageStartupMessages(library(visNetwork))
suppressPackageStartupMessages(library(igraph))
suppressPackageStartupMessages(library(shinycssloaders))
suppressPackageStartupMessages(library(plyr))
suppressPackageStartupMessages(library(ggnewscale))
suppressPackageStartupMessages(library(googledrive))

#BIOCONDUCTOR
suppressPackageStartupMessages(library(msa))
suppressPackageStartupMessages(library(ggtree))

#GITHUB
suppressPackageStartupMessages(library(r3dmol))
suppressPackageStartupMessages(library(dqshiny))

#-############################################################################-#
#-############################ SOURCE FUNCTIONS ##############################-#
#-############################################################################-#
source('code/Functions.R', local = TRUE)

#-############################################################################-#
#-############## MAPPING FILES, GOOGLE AUTHENTICATE AND VECTORS ##############-#
#-############################################################################-#

message(" Loading in data ...")
#keys    <- read.csv('data/ID_mapping.csv')
keys    <- read.csv('data/ID_mapping_google_drive.csv')

Authenticate_google_drive()

ClinVar <- read.csv('data/ClinVar_variants_in_sites.csv')

opts <- structure(.Data = paste0(keys$UniProtKB.AC,'|',stringr::str_trunc(keys$Protein.name, 70)), 
                  names = keys$Genes[order(keys$UniProtKB.AC)])

dis_opts <- c(opts,structure(.Data = paste0('Disease|',sort(unique(ClinVar$Phenotype))), 
                             names= stringr::str_trunc(sort(unique(ClinVar$Phenotype)), 67)))

message(" Done")

#-############################################################################-#
#-############################# USING ENTER KEY ##############################-#
#-############################################################################-#

js_submit <- '
$(document).keyup(function(event) {
  if ($("#proteinquery").is(":focus") && (event.keyCode == 13)) {
      $("#submit").click();
  }
});
'
js_submit_data <- '
$(document).keyup(function(event) {
  if ($("#proteinquery_data").is(":focus") && (event.keyCode == 13)) {
      $("#submit_data").click();
  }
});
'
js_submit_analysis <- '
$(document).keyup(function(event) {
  if ($("#proteinquery_analysis").is(":focus") && (event.keyCode == 13)) {
      $("#submit_analysis").click();
  }
});
'
#-############################################################################-#
#-################################ LAYOUTS ###################################-#
#-############################################################################-#

my_template <- gridTemplate("myareas", "grid", areas = list(
    default = c(
        "header header header",
        "left main right",
        "footer footer footer"),
    xs = c(
        "header",
        "left",
        "main",
        "right",
        "footer")
    ),
    gap = '5px'
    
)

my_template_analysis <- gridTemplate("myareas", "grid", areas = list(
  default = c(
    "header right",
    "left right",
    "left right"),
  xs = c(
    "header",
    "left",
    "right")
),
gap = '20px'

)

#-############################################################################-#
#-################################ HEADER ####################################-#
#-############################################################################-#

header_download <- list(
    list(
        key  = 'download', 
        text = "Download", 
        iconProps = list(iconName = "Download")
    )
)

header_about <- list(
    list(
        key  = 'about', 
        text = "About", 
        iconProps = list(iconName = "Info")
    )
)

search_bar <- flexPanel(
  flex = c(0,0), 
  direction='row',
  align_items = 'center', 
  tags$script(HTML(js_submit)),
  
  dqshiny::autocomplete_input(
    id = "proteinquery",
    label = NULL,
    options = opts,  
    value = NULL ,
    placeholder = "Search for protein, e.g. ATP2B2", 
    max_options = 20, 
    create = T, 
    contains = T
  ),
  actionButton(inputId = "submit", label = '',
               icon = icon("magnifying-glass", 
                           class = 'fa-solid', 
                           lib = 'font-awesome'), 
               style = 'background-color: steelblue; 
                        border-color: steelblue; 
                        color: #555; 
                        font-size: 30px; 
                        margin-left: 110px; 
                        margin-bottom: 15px; 
                        border-radius: 40px; 
                        width: 60px; 
                        height: 60px; 
                        box-shadow: 0px 0px 15px 15px steelblue; 
                        z-index: unset !important; 
                        z-index: 10000000000 !important;',
               onclick = paste0("location.href='",route_link("data"),"';")
  )
)

search_bar_data <- flexPanel(
  flex = c(0,0), 
  direction='row',
  align_items = 'center', 
  tags$script(HTML(js_submit_data)),
  
  dqshiny::autocomplete_input(
    id = "proteinquery_data",
    label = NULL,
    options = opts,  
    value = NULL ,
    placeholder = "Search for protein, e.g. ATP2B2", 
    max_options = 20, 
    create = T, 
    contains = T
  ),
  actionButton(inputId = "submit_data", label = '',
               icon = icon("magnifying-glass", 
                           class = 'fa-solid', 
                           lib = 'font-awesome'), 
               style = 'background-color: steelblue; 
                        border-color: steelblue; 
                        color: #555; 
                        font-size: 30px; 
                        margin-left: 110px; 
                        margin-bottom: 15px; 
                        border-radius: 40px; 
                        width: 60px; 
                        height: 60px; 
                        box-shadow: 0px 0px 15px 15px steelblue; 
                        z-index: unset !important; 
                        z-index: 10000000000 !important;'
  )
)

search_bar_analysis <- flexPanel(
  flex = c(0,0), 
  direction='row',
  align_items = 'center', 
  tags$script(HTML(js_submit_analysis)),
  
  dqshiny::autocomplete_input(
    id = "proteinquery_analysis",
    label = NULL,
    options = dis_opts,  
    value = NULL ,
    placeholder = "Search for protein or disease", 
    max_options = 20, 
    create = T, 
    contains = T
  ),
  actionButton(inputId = "submit_analysis", label = '',
               icon = icon("magnifying-glass", 
                           class = 'fa-solid', 
                           lib = 'font-awesome'), 
               style = 'background-color: steelblue; 
                        border-color: steelblue; 
                        color: #555; 
                        font-size: 30px; 
                        margin-left: 110px; 
                        margin-bottom: 15px; 
                        border-radius: 40px; 
                        width: 60px; 
                        height: 60px; 
                        box-shadow: 0px 0px 15px 15px steelblue; 
                        z-index: unset !important; 
                        z-index: 10000000000 !important;'
  )
)

app_header_front <- flexPanel(
    id = "header",
    align_items = "center",
    flex = c(0,1,0),
    #div(img(src = "ATP2B2.png", style = "width: 100px; margin-left: 20px; margin-top: 20px;")), #html hyphen &#8209; , space &nbsp;
    div(Text(HTML("MemRIDR&nbsp;db"), 
             style="color: steelblue; font-size: 60px; margin-left: 10px;", 
             class = "display nowrap")),
    div(Text('')),
    div(tags$ul(shiny::a(href = route_link("about"), "About", 
                         style = "text-decoration: none; color: gray; font-size: 20px; margin-right: 80px")))
)

app_header <- flexPanel(
    id = "header",
    align_items = "center",
    flex = c(0,1,0),
    div(tags$a(href=route_link("/"), img(src = "ATP2B2.png", 
                                         style = "width: 140px; margin-left: 20px; margin-top: 20px; cursor:pointer;"))),
    div(search_bar_data, 
        style = 'margin-top: 30px; margin-left: 20px;'),
    div(downloadButton(outputId = 'downloadbinding', 
                       label = 'Download', 
                       icon = icon('angles-down')), 
        style = 'margin-right: 80px')
)

app_header_analysis <- flexPanel(
    id = "header",
    align_items = "center",
    flex = c(0,1,0),
    div(tags$a(href=route_link("/"), img(src = "ATP2B2.png", 
                                         style = "width:140px; margin-left: 20px; margin-top: 20px; cursor:pointer;"))),
    div(search_bar_analysis, style = 'margin-top: 30px; margin-left: 20px;'),
    downloadButton(outputId = 'downloadclinvar', 
                   label = 'Download', 
                   icon = icon('angles-down'))
)

app_header_about <- flexPanel(
    id = "header",
    align_items = "center",
    flex = c(0,1,0),
    div(Text(HTML("MemRIDR&nbsp;db"), 
             style="color: steelblue; font-size: 60px; margin-left: 10px;", 
             class = "title")),
    div(Text(''))
)

#-############################################################################-#
#-################################# FOOTER ###################################-#
#-############################################################################-#

app_footer <- flexPanel(
    id = "footer",
    align_items = "center",
    justify_content = 'space-evenly',
    style = 'border-radius: 20px; background-color: #fff; margin-top: 20px;',

    div(
        tags$a(img(src = "UCPH_logo.png", 
                   style = "margin-left: 70px; width: 300px; cursor:pointer;"),
               href="https://in.ku.dk/research/madsen-lab/", 
               target="_blank")#,
        #tags$a(img(src = "LF_logo.png",   style = "margin-left: 30px; width: 300px; cursor:pointer;"),href="https://lundbeckfonden.com/")
        ),
    div('If you use this database in a publication, please cite: Rombach, J., Nielsen, T. T. E., .....')
)


app_footer_data <- flexPanel(
  id = "footerdata",
  align_items = "center",
  justify_content = 'space-evenly',
  style = 'border-radius: 0px; background-color: #fff; margin-top: 20px;',
  div(
    tags$a(img(src = "UCPH_logo.png", 
               style = "margin-left: 70px; width: 300px; cursor:pointer;"),
           href="https://in.ku.dk/research/madsen-lab/", 
           target="_blank")#,
    #tags$a(img(src = "LF_logo.png",   style = "margin-left: 30px; width: 300px; cursor:pointer;"),href="https://lundbeckfonden.com/")
  ),
  div('If you use this database in a publication, please cite: Rombach, J., Nielsen, T. T. E., .....')
)

#-############################################################################-#
#-############################## FRONT PAGE ##################################-#
#-############################################################################-#

home_page <- gridPanel(
        id = 'home_page',
        template = my_template,
        rows = '1fr 3fr, 1fr',
        
        header = app_header_front,
        
        left = fluidRow(),
        
        main = flexPanel(
                id = "main",
                align_items = "center",
                gap = 20,
                justify_content = 'space-evenly', 
                flex = c(1),
                img(src = "ATP2B2.png", style = 'width: 650px;'),
                div(search_bar)),

        right = fluidRow(),
        
        footer = app_footer
)

#-############################################################################-#
#-############################### DATA PAGE ##################################-#
#-############################################################################-#

data_page <- gridPanel(
    
    id = 'data',
    template = my_template,
    columns = "1fr 1fr 2fr",
    rows = "0.5fr 4fr 0.5fr",
    
    header = app_header,
    
    left = fluidRow(id = 'datapanel',
                    
                    column(12,
                           dataTableOutput('proteintable') %>% withSpinner(color="steelblue", type = 7, hide.ui = F),
                           style='padding-left:0px; padding-right:0px; padding-top:0px; padding-bottom:10px; height:200px;' 
                           ),
                    column(6,
                           checkboxGroupInput("terminal", 
                                              "Select terminals", 
                                              c('N','C'), 
                                              inline = T, 
                                              selected = c('N','C'))
                    ),
                    column(6,
                           radioButtons("confidence", 
                                        "Set confidence", 
                                        c('low','medium','high'), 
                                        inline = T, 
                                        selected = c('medium'))
                    ),
                    column(6,
                           actionButton(inputId = "analysis",
                                        label = 'ClinVar variants network', 
                                        icon = icon("diagram-project", lib = 'font-awesome'), 
                                        style = 'background-color: steelblue; border-color: steelblue; color: #fff;',
                                        onclick = paste0("location.href='", route_link("analysis"),"';")
                                        ), 
                           style='padding-left:15px; padding-right:0px; padding-top:5px; padding-bottom:0px;' 
                           ),
                    column(6,
                           actionButton(inputId = "AH",
                                        label = 'Detect amphipathic helix', 
                                        icon = icon("magnifying-glass-chart"), 
                                        style = 'background-color: steelblue; border-color: steelblue; color: #fff;',
                                        class = "btn-success"),
                           style='padding-left:15px; padding-right:0px; padding-top:5px; padding-bottom:0px' 
                          ),
                    column(12,
                           plotOutput('helicalwheel', width = '90%') %>% withSpinner(color="steelblue", hide.ui = T, type = 7),
                           style='padding-left: 10%; padding-right:0px; padding-top:0px; padding-bottom:0px;'
                           ),
                    column(2, offset = 4,
                           uiOutput("UniProtLink_button"), 
                           style='padding-rigth:30px, padding-top:-100px; padding-bottom:10px;'
                    ),
                    ),
    
    main  =  flexPanel(r3dmolOutput('pdb', width = 'auto', height = 'auto')),
    
    right =  fluidRow(id = 'datapanel',

                      column(12,
                             div(Text(variant = "xLarge", 'Membrane binding traces'), 
                                 style='padding-left:22px; padding-right:1px; padding-top:20px; padding-bottom:0px; height: 10%;')
                             ),
                      column(12,
                             plotOutput('trace') %>% withSpinner(color="steelblue", hide.ui = F, type = 7), 
                             style='padding-left:10px; padding-right:20px; padding-top:0px; padding-bottom:0px; height: 50%;' 
                             ),
                      column(12, 
                             dataTableOutput('siteTable') %>% withSpinner(color="steelblue", hide.ui = F, type = 7), 
                             style='padding-left:0px; padding-right:0px; padding-top:0px; padding-bottom:20px; height: 40%;' 
                             )
                      ),
    
    footer = app_footer_data
)

#-############################################################################-#
#-############################# ANALYSIS PAGE ################################-#
#-############################################################################-#

analysis_page <- gridPanel(
  
  id = 'analysis',
  template = my_template_analysis,
  columns = "1.5fr 3.5fr",
  rows = "0.5fr 4.5fr",
  
  header = app_header_analysis,
  
  left = fluidRow(id = 'datapanel',
    column(12, 
           dataTableOutput('clinvartable') %>% withSpinner(color="steelblue", hide.ui = F, type = 7), 
           style='padding-left:0px; padding-right:0px; padding-top:0px; padding-bottom:0px;' 
           )
  ),

  right = fluidRow(fillPage(
    tags$style(type = "text/css", "#network {height: calc(100vh - 25px) !important;}"),
    visNetworkOutput("network") #%>% withSpinner(color="steelblue", type = 7, hide.ui = F)
    )
  )
)

#-############################################################################-#
#-############################## ABOUT PAGE ##################################-#
#-############################################################################-#

about_page <- gridPanel(
                    id = 'about',
                    template = my_template, #'holy-grail',
                    rows = '2fr 4fr, 1fr',
                    columns = '1.5fr 2fr 1.5fr',
                    
                    header = div(app_header_about,style = 'padding-bottom:50px;'),
                    
                    main = flexPanel(justify_content = 'left',
                                     direction = 'column',
                                     style = 'border-radius: 20px; box-shadow: 0 0 10px steelblue; text-align: center; width: 50hw;',
                                
                                tags$html(
                                  tags$body(
                                h1('What is this all about?', 
                                   style = 'padding-top:30px; padding-bottom:10px; padding-left:30px; padding-right:30px;'),
                                
                                p('Transmembrane proteins (TMPs) serve crucial functions including receptor-, transporter- and channel-activity. 
                                   Besides our in-depth understanding of their transmembrane domains, TMPs house extensive intrinsically disordered regions (IDRs) 
                                   whose functions remain largely unknown. IDR function is highly dependent on context, a trait particularly relevant for IDRs 
                                   anchored to membranes. Here, we used peptide arrays and fluorescent liposomes to create a proteome-wide map of 
                                   membrane-interacting regions in IDRs of all human TMPs. For more information see publication Rombach, Nielsen... (2025). 
                                   With MemRIDR (Membrane Recruited Intrinsically Disordered Regions), we provide access and a basic analysis tool, which 
                                   we hope will help the community in discovering new unconventional structural or regulatory mechanisms in TMPs.', 
                                  style = 'padding-top:30px; padding-bottom:30px; padding-left:30px; padding-right:30px; text-align: justify;')), 
                                ),
                                
                                div(tags$ul(shiny::a(href  = route_link("/"), "Search for Protein", class = "btn btn-primary", 
                                                     style = "text-decoration: none; font-size: 20px; margin-right: 40px; margin-top: 20px; text-align: center; background-color: steelblue; border-color: steelblue; border-radius: 40px; width: 200px; box-shadow: 0px 0px 15px 15px steelblue;"))),
                                
                                
                                div(downloadButton(outputId = 'downloadallbinding',    
                                                   label = 'Download all membrane binding data', 
                                                   icon  = icon('angles-down')), 
                                                   style = 'padding-top:50px; padding-bottom:10px;'),
                                
                                div(downloadButton(outputId = 'downloadallclinvar', 
                                                   label = 'Download all ClinVar variants in membrane binding sites', 
                                                   icon  = icon('angles-down')), 
                                                   style = 'padding-top:10px; padding-bottom:60px;')
                        ),
                    
                    footer = app_footer_data
)

#-############################################################################-#
#-################################### UI #####################################-#
#-############################################################################-#

ui <- fluidPage(
  
    tags$head(
        tags$title('MemRIDR db'),
        tags$link(href = "style.css", rel = "stylesheet", type = "text/css"),
        
        tags$style(type="text/css", "#proteinquery, #proteinquery_data, #proteinquery_analysis {
                                     border-radius: 30px;
                                     border: 1px solid #777;
                                     height: 50px; 
                                     width: 450px; 
                                     text-align: center; 
                                     font-size: 20px;}"),
        
        tags$style(type="text/css", "#proteinquery:active, #proteinquery_data:active, #proteinquery_analysis:active
                                     #proteinquery:hover #proteinquery_data:hover #proteinquery_analysis:hover{
                                     border-radius: 30px; 
                                     border: 1px solid steelblue !important;
                                     box-shadow: 0 0px 10px steelblue !important;
                                     height: 50px; 
                                     width: 450px; 
                                     text-align: center; 
                                     font-size: 20px;}"),
        
        tags$style(type="text/css", "#footer {opacity: 0.5;}"),
        tags$style(type="text/css", "#footer:hover{opacity: 1;}"),
        tags$style(type="text/css", "#footerdata {opacity: 0.5;}"),
        tags$style(type="text/css", "#footerdata:hover{opacity: 1}"),
        tags$style(type="text/css", "#AH .fa-magnifying-glass-chart {color:#fff}"),
        tags$style(type="text/css", "#analysis .fa-diagram-project {color:#fff}"),
        tags$style(type="text/css", "td {white-space: nowrap; max-width: 200px; overflow: hidden; text-overflow: ellipsis;}"),
        tags$style(type="text/css", "#datapanel {border-radius: 30px; box-shadow: 0 0 10px steelblue; margin: 0px; padding: 0px;}"),
        
        tags$style(type="text/css", "#downloadclinvar, #downloadallclinvar, #downloadallbinding, 
                                     #downloadclinvar:active, #downloadbinding, #downloadbinding:active {
                                      background-color:rgba(0,0,0,0);
                                      color: steelblue;
                                      font-family: Arial;
                                      border-color: transparent;
                                      -webkit-box-shadow: 0px;
                                      box-shadow: 0px;
                                    }"),
        tags$style('/* Works on Firefox */
                  * {
                    scrollbar-width: thin;
                    scrollbar-color: transparent transparent;
                    }
                  
                  /* Works on Chrome, Edge, and Safari */
                  *::-webkit-scrollbar 
                    width: 0px;
                    }
                  
                  *::-webkit-scrollbar-track {
                    background: transparent;
                    border-radius: 30px; 
                    }
                  
                  *::-webkit-scrollbar-thumb {
                    background-color: transparent;
                    border: 0px solid transparent;
                    }'),
        
        tags$style("html, body {overflow: visible !important;"),
        tags$style(HTML(
          'table.dataTable tbody tr.selected td,
           table.dataTable tbody td.selected {
              border-top-color: white !important;
              box-shadow: inset 0 0 0 9999px steelblue !important;
              }
          
            table.dataTable tbody tr:active td {
              background-color: steelblue !important;
              }
          
            :root {
              --dt-row-selected: transparent !important;
              }
          
            table.dataTable tbody tr:hover, table.dataTable tbody tr:hover td {
              background-color: steelblue !important;
              }')),
        tags$style(HTML(
          '.autocomplete {
              position: relative;
              border: 0px solid transparent;
              }
           .autocomplete-items div {
              background-color: #fff;
              border-bottom: 0px white;
              margin-left: 0px;
              border: 0px solid transparent;
              }
           .autocomplete-items div:hover {
              background-color: steelblue;
              color: #ffffff;
              border: 0px solid transparent;
              }
          .autocomplete-items {
            	position: absolute;
            	border: 0px solid transparent;
              box-shadow: 0 0px 0px transparent;
            	border-bottom: none;
            	border-top: none;
            	z-index: 999;
            	top: 100%;
            	left: 0;
            	right: 0;
            	max-height: 300px;
            	overflow: auto;
              padding-left: 22px;
              padding-right: 22px;
              border-bottom-left-radius: 0px; 
              border-top-left-radius: 0px; 
              border-bottom-right-radius: 0px; 
              border-top-right-radius: 0px;
              }
          .autocomplete-active {
              background-color: steelblue !important;
              color: #ffffff;
              padding-left: 22px;
              padding-right: 22px;
              border: 0px solid transparent;
              }'))
        ),

      router_ui(
        route("/",        home_page),
        route("data",     data_page),
        route("analysis", analysis_page),
        route("about",    about_page))
)

#-############################################################################-#
#-################################# SERVER ###################################-#
#-############################################################################-#

server <- function(input, output, session) {
  
    # Get query ----------------------------------------------------------------
    my_protein <- reactiveValues(query='')
    my_disease <- reactiveValues(query=NULL)
    highlight  <- reactiveValues(nodes=NULL)
    
    observeEvent(input$proteinquery,{
      my_protein$query <- gsub('\\|.*','',input$proteinquery)
    }, ignoreInit = T)
    
    observeEvent(input$proteinquery_data,{
      my_protein$query <- gsub('\\|.*','',input$proteinquery_data)
    }, ignoreInit = T)
    
    observeEvent(input$proteinquery_analysis,{
      
      if(gsub('\\|.*','',input$proteinquery_analysis)=='Disease'){
        my_protein$query <- input$proteinquery_analysis
      } else {
        search <- tolower(input$proteinquery_analysis)
        if(search==''){
          search <- 'epilepsy'
        }
        
        search <- str_replace(search, '\\(','\\\\(')
        search <- str_replace(search, '\\)','\\\\)')
        
        if(any(str_detect(ClinVar$Phenotype,search))){
          my_protein$query <- paste0('Disease|', search)
        } else {
          my_protein$query <- gsub('\\|.*','',input$proteinquery_analysis)
        }
      }
      
    }, ignoreInit = T)
  
    # Search for query matches -------------------------------------------------
    selectprotein <- eventReactive({c(input$submit,
                                      input$submit_data,
                                      input$submit_analysis
                                      )}, {
                                        
        if(gsub('\\|.*','',my_protein$query)=='Disease'){

          res <- search_for_disease(search = gsub('.*\\|','',my_protein$query), keys = keys, ClinVar = ClinVar)
          
          if(nrow(res)<3){
            
            empty_rows <- data.frame(Name=rep('',3-nrow(res)),
                                     Genes=rep('',3-nrow(res)),
                                     UniProtKB.AC=rep('',3-nrow(res)),
                                     Search=rep('',3-nrow(res)),
                                     Path=rep(res$Path[nrow(res)],3-nrow(res)))
            
            res <- rbind(res,empty_rows)
          }
          
          my_disease$query <- unique(res$Search)
          highlight$nodes  <- unique(res$Search)
          
        } else {
          
          res <- search_for_protein(search = isolate(my_protein$query), keys = isolate(keys))
          if(nrow(res)<3){
            
            empty_rows <- data.frame(Name=rep('',3-nrow(res)),
                                     Genes=rep('',3-nrow(res)),
                                     UniProtKB.AC=rep('',3-nrow(res)),
                                     Search=rep('',3-nrow(res)),
                                     Path=rep(res$Path[nrow(res)],3-nrow(res)))
           
            res <- rbind(res,empty_rows)
          }
          my_disease$query <- NULL
        }
        rownames(res) <- NULL
        selectprotein <- res
    }, ignoreInit = T)
    
    # Plot table with hits for search query ------------------------------------
    output$proteintable = DT::renderDT({datatable(
                                       selectprotein()[,c(3,1,2)],
                                       plugins = "ellipsis",
                                       selection = list(
                                                        mode='single',
                                                        selected=1
                                                        ),
                                       options=list(pageLength=3, 
                                                    dom = 'tp',
                                                    #autoWidth = T,
                                                    #responsive = TRUE,
                                                    #scrollX=T,
                                                    columnDefs = list(
                                                        list(
                                                        targets = c(2,3),
                                                        render = JS("$.fn.dataTable.render.ellipsis(14, false )")
                                                        ),
                                                        list(
                                                            width = '100px', 
                                                            targets = c(1)
                                                            )
                                                        )
                                                    )
                                              )
                                       })

    # Import Array_object ------------------------------------------------------
    
    # Local
    # Array_object_path <- reactive({
    #                           if(!is.null(input$proteintable_rows_selected)){
    #                             Array_object_path <- isolate(selectprotein())[input$proteintable_rows_selected,5]
    #                           }
    #                         })
    # Array_object <- reactive({
    #                           if(!is.null(input$proteintable_rows_selected)){
    #                             Array_object <- readRDS(Array_object_path())
    #                           }
    #                         })
    
    #Google drive
    Array_object_path <- reactive({
      if(!is.null(input$proteintable_rows_selected)){
        Array_object_path <- isolate(selectprotein())[input$proteintable_rows_selected,5]
      }
    })
    Array_object <- reactive({
      if(!is.null(input$proteintable_rows_selected)){
        Array_object <- get_array_file_google_drive(Array_object_path(), google_drive_files)
      }
    })

    # Select terminal to focus on ----------------------------------------------
    UniProt <- eventReactive(c(Array_object(),input$terminal),{
      if(!is.null(input$proteintable_rows_selected)){
        my_row <- input$proteintable_rows_selected
      } else {
        my_row <- 1
      }
          if(!is.null(my_row)){
            
            while(isolate(selectprotein())[my_row,1]==''){
              my_row <- my_row-1
            }
            
            if(length(input$terminal)==1){
              UniProt <- paste0(isolate(selectprotein())[my_row,3],'_',input$terminal)
            } else {
              UniProt <- isolate(selectprotein())[my_row,3]
            }
            
          } else {
                if(length(input$terminal)==1){
                   UniProt <- paste0(isolate(selectprotein())[1,3],'_',input$terminal)
                } else {
                   UniProt <- isolate(selectprotein())[1,3]
                }
          }
      }, ignoreInit = T)
    
    # Create button that links to UniProt --------------------------------------
    output$UniProtLink_button <- renderUI({
      shiny::a(
        h4(img(src='UniProtLogo.png', width='100px'),
           class = "btn btn-default action-button",
           style = "height: auto; 
                    width: auto; 
                    max-width: 120px; 
                    max-height: 60px;
                    border-color: #fff;
                    border-radius: 30px;"),
        target = "_blank",
        href = paste0("https://www.uniprot.org/uniprotkb/", selectprotein()[input$proteintable_rows_selected,3], "/entry")
      )
    })
    
    # Set confidence threshold -------------------------------------------------
    conf <- reactive({input$confidence})
    
    base_pdb <- reactive({
      input$terminal #trigger
      if(!is.null(Array_object())){
      base_pdb <- quiet(plot_pdb_with_binding_sites(UniProtKB.AC = isolate(UniProt()),
                                        array_object = Array_object(), 
                                        pthresh = isolate(conf()), 
                                        zoom_limits = c(5,1000),
                                        background_color='white'))
      }
    })
    
    # Plot protein structure with SWaFi score ----------------------------------
    pdb <- reactive({
      if(!is.null(base_pdb())){
        
        pdb <- base_pdb()
        
        if(!is.null(input$siteTable_rows_selected)){
          
          sites_to_highlight <- isolate(sites())[input$siteTable_rows_selected,1:2]
          sites_to_highlight <- sites_to_highlight[!is.na(sites_to_highlight$start),] 
          sites_to_highlight <- rbind(data.frame(start=1, end=1),sites_to_highlight)

          highlight <- unlist(c(apply(sites_to_highlight,1,function(x){
            x[1]:x[2]
          })))[-1]
          
          if(length(highlight)>0){
            pdb <- pdb %>% m_set_style(style = m_style_cartoon(color = 'red'),
                                       sel   = m_sel(resi = highlight))
          }
          
        }
        return(pdb)
      }
    })

    # Render protein structure -------------------------------------------------
    output$pdb <- renderR3dmol({pdb()})
    
    # Get SWaFi trace for download ---------------------------------------------
    trace <- reactive({
      
      trace <- get_averaged_SWaFi_traces(proteins = UniProt(),
                                         datasets = c('FBBE1','FBBE2'), 
                                         array_object = isolate(Array_object()))
    })
    

    # Plot SWaFi trace ---------------------------------------------------------
    selecttrace <- reactive({
      suppressWarnings(
              plot_SWaFi_trace_w_sites(array_object = Array_object(), 
                                       proteins     = UniProt(), 
                                       sites        = selected_sites(),
                                       traces       = isolate(trace()),
                                       confidence   = conf(), 
                                       trace_color  = 'steelblue',
                                       highlight_color = 'red')
                       )
                            })
    
    output$trace <- renderPlot(selecttrace())
    
    # Get binding all sites ----------------------------------------------------
    binding_sites <- reactive({
      binding_sites <-  quiet(get_binding_sites(proteins     = isolate(UniProt()),
                                                confidence   = conf(),
                                                array_object = Array_object())
                              ) 
    })
    
    # Select terminal to display -----------------------------------------------
    sites <- reactive({
      
      if(!is.null(binding_sites())){
          temp_sites <- binding_sites()
          
        if(detect_identifier(UniProt())=='ArrayID'){
          temp_sites <- temp_sites[temp_sites$UniProtKB.AC==isolate(UniProt()),]
        }
          
        sites <- temp_sites[,c(3:6,8,10)]  %>% 
                 dplyr::mutate(across(where(is.numeric), ~round(.x,digits = 5)))
      }
    })
    
    # Selecting rows in sites --------------------------------------------------
    rows <- reactiveValues(selected=NULL)
    observeEvent(input$siteTable_rows_selected,{
      rows$selected <- input$siteTable_rows_selected
    }, ignoreInit = T,ignoreNULL = FALSE)
    observeEvent(Array_object(),{
      rows$selected <- NULL
    }, ignoreInit = T)
    
    selected_sites <- reactive({
                            if(!is.null(rows$selected)){
                              selected_sites <-  sites()[rows$selected,]
                            } else {
                              selected_sites <- sites()[sites()$fisher_pval<conf_convert(isolate(conf())),]
                            }
                          })
    
    # Create binding site table ------------------------------------------------
    output$siteTable = DT::renderDT({datatable(sites(),
                                       plugins = "ellipsis",
                                       selection = list(
                                         mode='multiple'
                                       ),
                                       options=list(pageLength=6, 
                                                    dom = 'tp',
                                                    columnDefs = list(
                                                      list(
                                                        targets = c(6),
                                                        render = JS("$.fn.dataTable.render.ellipsis( 14, false )")
                                                      ),
                                                      list(
                                                        targets = c(1,2),
                                                        width = '30px'
                                                      )
                                                    )
                                       )
      )
    })
    
    # Select sites to focus on AH analysis -------------------------------------
    activeSite <- reactive({
                        if(!is.null(input$siteTable_rows_selected)){
                          activeSite <- isolate(sites())[input$siteTable_rows_selected,6]
                        } else {
                          activeSite <- NULL
                        }
                  })
    
    # Get sequence from selected site ------------------------------------------
    selectedsequence <- reactive({
                            if(!is.null(activeSite())){
                              if(length(activeSite()$sequence)<5){
                                selectedsequence <- quiet(get_best_AH(activeSite()$sequence,isolate(UniProt()), array_object = isolate(Array_object()), return_all_info = F))
                              } else {
                                selectedsequence <- c('too','many','sites','have','been','selected')
                              }
                            } else {
                              selectedsequence <- 'not_selected'
                            }
                        })
    
    # Triggers for plotting calculating AH -------------------------------------
    plot_AH <- reactiveValues(selected=0)
    
    observeEvent(input$AH, {
      plot_AH$selected <<- input$AH
    })
    
    observeEvent(c(input$siteTable_rows_selected,UniProt()),{
      plot_AH$selected <<- 0
    },ignoreInit=T)
    
    helicalwheel <- reactive({
      if(!is.null(plot_AH$selected) & plot_AH$selected>0){
        
        if(length(isolate(selectedsequence()))>4){
          helicalwheel <- ggplot() + 
            annotate(geom="text", x=0, y=0, label='Select a\nmaximum of\n4 sites', size = unit(12, "pt")) + 
            prism() + remove_x_axis() + remove_y_axis() + theme(legend.position='none') +
            xlim(-1.5,1.5) +
            ylim(-1.5,1.5)
        } else {
        
        if(isolate(selectedsequence())[1] %in% 'not_selected'){
          helicalwheel <- ggplot() + 
            annotate(geom="text", x=0, y=0, label='No site selected', size = unit(12, "pt")) + 
            prism() + remove_x_axis() + remove_y_axis() + theme(legend.position='none') +
            xlim(-1.5,1.5) +
            ylim(-1.5,1.5)
        } else {
        
        helicalwheel <- plot_list_of_ggplots(lapply(isolate(selectedsequence()), function(seqs){
          
          if(!is.na(seqs)){
            
            temp_wheel <- heliquest(seqs, plot=F)$plot
            
          } else {
            # plot empty graph with message
            temp_wheel <- ggplot() + 
              annotate(geom="text", x=0, y=0, label='No amphipathic\nhelix found', 
                       size = unit(12/ceiling(sqrt(length(isolate(selectedsequence())))), "pt")) + 
              prism() + remove_x_axis() + remove_y_axis() + theme(legend.position='none') +
              xlim(-1.5,1.5) +
              ylim(-1.5,1.5)
          }
          
          return(temp_wheel)
          
         }), ncol=ceiling(sqrt(length(isolate(selectedsequence())))))
        }
        }
       return(helicalwheel)
      }
    })
      
    # Render helical wheel -----------------------------------------------------
    output$helicalwheel <- renderPlot({helicalwheel()})

    # Import disease graph and data --------------------------------------------
    g <- readRDS("data/disease_graph.rds")

    # Evaluate whether to get network  -----------------------------------------
    analysis_button <- reactiveValues(selected=0)
    observeEvent(c(input$submit, input$submit_data),{
      analysis_button$selected <- 0
    }, ignoreInit = T)
    observeEvent(c(input$analysis,input$submit_analysis),{
      analysis_button$selected <- 1
    }, ignoreInit = T)
    
    # Get nodes selected in disease graph --------------------------------------
    myNode <- reactiveValues(selected=NULL)
    
    observeEvent(input$current_node_id, {
      myNode$selected <- input$current_node_id
    }, ignoreInit = T)
    
    # Get sub network ----------------------------------------------------------
    net <- reactive({
      
      Array_object()
      #input$analysis 
      #input$submit
      #input$submit_data
      #input$submit_analysis
      
      if(!is.null(Array_object())){
        uniprot <- isolate(strip_terminal_specifier(UniProt()))
        
        if(uniprot %in% ClinVar$UniProt){
          myNode$selected <- unique(ClinVar$Gene[ClinVar$UniProt == uniprot])
        }
        
        if(!is.null(analysis_button$selected) & analysis_button$selected>0){
          
          net <- plot_vis_network(UniProtKB.AC = uniprot, 
                                  diseases = isolate(my_disease$query),
                                  graph = g, 
                                  protein_col = 'steelblue', 
                                  disease_col = 'red') 
          
          myNode$selected <- isolate(my_disease$query)
          my_disease$query <- NULL
        } else {
          
          #message(' turning off interactive network')
          myNode$selected <- NULL
          net <- plot_vis_network(UniProtKB.AC = 'none', graph = g, 
                                  protein_col = 'steelblue', disease_col = 'red')
        }
        return(net)
      }
    })
    
    # Plot interactive network -------------------------------------------------
    output$network <- renderVisNetwork({net()})
    observe({
      visNetworkProxy("network") %>%
        visSelectNodes(id=c(ClinVar_df()$Phenotype[input$clinvartable_rows_selected], 
                            ClinVar_df()$Gene[input$clinvartable_rows_selected]))
    })
    
    # Create ClinVar data frame subset -----------------------------------------
    ClinVar_df <- reactive({input$analysis # triggers
                            #input$submit_analysis
                            input$current_node_id
                            
                            if(!is.null(Array_object())){
                              
                            uniprot <- isolate(strip_terminal_specifier(UniProt()))
                            
                             if(uniprot %in% ClinVar$UniProt){
                               
                               if(is.null(myNode$selected)){
                                 
                                 temp_node <- unique(ClinVar$Gene[ClinVar$UniProt == uniprot])
                                 ClinVar_df <- ClinVar[ClinVar$Gene %in% temp_node,]
                                 
                               } else {
                                 
                                 temp_node <- myNode$selected
                                 if(all(temp_node %in% ClinVar$Gene)){
                                   ClinVar_df <- ClinVar[ClinVar$Gene %in% temp_node,]
                                 } else {
                                   ClinVar_df <- ClinVar[ClinVar$Phenotype %in% temp_node,]
                                 }
                                 
                               }
                             } else {
                               ClinVar_df <- ClinVar[0,]
                             }
                            } else {
                              ClinVar_df <- ClinVar[0,]
                            }
                        })
    
    # Render ClinVar subset data frame -----------------------------------------
    output$clinvartable <- DT::renderDT({datatable(
      ClinVar_df()[,c(11,2,4,5,7,8,9,12,19)],
      plugins = "ellipsis",
      selection = list(mode='single', 
                       selected=1),
      options=list(pageLength=18, 
                   dom = 'tp',
                   autoWidth = T,
                   #responsive = TRUE,
                   scrollX=T,
                   columnDefs = list(
                     list(targets = c(3,8,9),width = '100px'),
                     list(targets = c(1,2,4),width = '40px' ),
                     list(targets = c(5,6,7),width = '25px' ),
                     list(targets = c(0),    width = '0px'  ),
                     list(targets = c(9),   render = JS("$.fn.dataTable.render.ellipsis( 12, false )")),
                     list(targets = c(3,8), render = JS("$.fn.dataTable.render.ellipsis( 12, false )")),
                     list(targets = c(4,7), render = JS("$.fn.dataTable.render.ellipsis( 6, false )" )),
                     list(targets = c(5,6), render = JS("$.fn.dataTable.render.ellipsis( 4, false )" ))
                    )
                   )
                  )
                }) 
    
    # download handlers --------------------------------------------------------
    
    # Binding site data
    output$downloadbinding <- downloadHandler(
      filename = function(){paste0(gsub('-','',Sys.Date()),'_',UniProt(),'_data.zip')},
      content = function(file){
        
        temp_directory <- file.path(tempdir(), as.integer(Sys.time()))
        dir.create(temp_directory)
        
        quiet(ggsave(filename = file.path(temp_directory, paste0(UniProt(),'_membrane_binding_plot.pdf')), 
               plot = selecttrace(), 
               device = "pdf", width = 15, height = 7))
        
        if(!is.null(plot_AH$selected) & plot_AH$selected>0){
          quiet(ggsave(filename = file.path(temp_directory, paste0(UniProt(),'_amphipathic_helices.pdf')), 
                       plot = helicalwheel(), 
                       device = "pdf", width = 10, height = 10))
        }
        
        # Work around for self contained, needs to be in current directory to work!
        wd_temp <- getwd()
        setwd(temp_directory)
        
        htmlwidgets::saveWidget(widget = pdb(),
                                file = file.path(temp_directory, paste0(UniProt(),'_pdb.html')),
                                selfcontained = TRUE)
        setwd(wd_temp)
        
        trace_download <- do.call(rbind,trace())
        readr::write_csv(trace_download[order(trace_download$position),],  file.path(temp_directory, paste0(UniProt(),'_membrane_binding_trace.csv')))
        readr::write_csv(sites(),  file.path(temp_directory, paste0(UniProt(),'_membrane_binding_sites.csv')))
        
        zip::zip(
          zipfile = file,
          files = dir(temp_directory),
          root = temp_directory
        )
        
      },
      contentType = "application/zip"
      )
    
    # ClinVar data
    output$downloadclinvar <- downloadHandler(
      filename = function(){paste0(gsub('-','',Sys.Date()),'_',UniProt(),'_ClinVar_data.zip')},
      content = function(file){
        
        temp_directory <- file.path(tempdir(), as.integer(Sys.time()))
        dir.create(temp_directory)
        
        readr::write_csv(ClinVar_df(),  file.path(temp_directory, paste0(UniProt(),'_ClinVar_data.csv')))
     
        # Work around for self contained, needs to be in current directory to work!
        wd_temp <- getwd()
        setwd(temp_directory)

        htmlwidgets::saveWidget(widget = net(),
                                file = file.path(temp_directory, paste0(UniProt(),'_ClinVar_network.html')),
                                selfcontained = TRUE)
        setwd(wd_temp)
        
        zip::zip(
          zipfile = file,
          files = dir(temp_directory),
          root = temp_directory
        )
        
      },
      contentType = "application/zip"
      )
    
    # All ClinVar data
    output$downloadallclinvar <- downloadHandler(
      filename = 'all_ClinVar_data.csv',
      content = function(file) {
        file.copy("data/ClinVar_variants_in_sites.csv", file)
      }
      )
    
    # All binding data
    output$downloadallbinding <- downloadHandler(
      filename = 'all_membrane_binding_sites.csv',
      content = function(file) {
        file.copy("data/all_sites.csv", file)
      }
      )
    
    # setup multipage app ------------------------------------------------------
    router_server()    
    
    # Stop app on window closing - remove before deploying! --------------------
    # session$onSessionEnded(function(){
    #   cat("Session Ended\n")
    #   stopApp()
    # })
    

}

#-############################################################################-#
#-################################### RUN ####################################-#
#-############################################################################-#

shinyApp(ui = ui, server = server, options = list(launch.browser = TRUE))
