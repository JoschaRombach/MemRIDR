# Author: Joscha Rombach
# Created: 12/04/2024
# Project: Web interface for array_object




################################################################################
### TO RUN THE APP PRESS "Run App" IN THE TOP RIGHT CORNER OF THIS PANEL -> ####
################################################################################




# ################################################################################
# ################################## READ DATA ###################################
# ################################################################################

suppressPackageStartupMessages(library(here))

Array_object <<- readRDS(here('data/Array_objects_single_proteins','Array_object_Q9P2U7.rds'))
#Array_object$Import_functions()
source(here('code','Functions.R'))
keys <- read.csv('data/ID_mapping.csv')

################################################################################
############################### USING ENTER KEY ################################
################################################################################

js <- '
$(document).keyup(function(event) {
  if ($("#proteinquery").is(":focus") && (event.keyCode == 13)) {
      $("#submit").click();
  }
});
'

################################################################################
################################## LAYOUTS #####################################
################################################################################

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
    gap = '20px'
)

################################################################################
################################## HEADER ######################################
################################################################################

header_download <- list(
    list(
        key  = 'download', 
        text = "Download data", 
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
                tags$script(HTML(js)),
                textInput(inputId = "proteinquery", 
                        label = NULL, 
                        value = NULL, 
                        placeholder = 'Search for protein, e.g. vGluT1'),
                actionButton(inputId = "submit", 
                             label = icon("magnifying-glass"), 
                             style = 'background-color: #0275d8; border-color: #0275d8; font-color: #fff; font-size: 30px; margin-left: 110px; margin-bottom: 15px; border-radius: 40px; width: 60px; height: 60px; box-shadow: 0px 0px 15px 15px #0275d8;',
                             onclick = paste0("location.href='",route_link("data"),"';")
                             )
                  )

app_header_front <- flexPanel(
    id = "header",
    align_items = "center",
    flex = c(0,1,0),
    div(img(src = "vglut1.png", style = "width: 100px; margin-left: 20px; margin-top: 20px;")),
    div(Text("MemBindMap", style="color: #0B51A0; font-size: 60px; margin-left: 10px;", class = "title")),
    div(tags$ul(shiny::a(href = route_link("about"), "About", style = "text-decoration: none; color: gray; font-size: 20px; margin-right: 80px")))
)

app_header <- flexPanel(
    id = "header",
    align_items = "center",
    flex = c(0,1,0),
    div(tags$a(href=route_link("/"), img(src = "vglut1.png", style = "width: 100px; margin-left: 20px; margin-top: 20px; cursor:pointer;"))),
    div(search_bar, style = 'margin-top: 30px; margin-left: 20px;'),
    div(CommandBar(items = header_download), style = 'margin-right: 80px')
)

app_header_about <- flexPanel(
    id = "header",
    align_items = "center",
    flex = c(0,1,0),
    div(tags$a(href=route_link("/"), img(src = "vglut1.png", style = "width: 100px; margin-left: 20px; margin-top: 20px; cursor:pointer;"))),
    div(Text("MemBindMap", style="color: #0B51A0; font-size: 60px; margin-left: 10px;", class = "title")),
    div(tags$ul(shiny::a(href = route_link("data"), "Search for Protein", class = "btn btn-primary", 
                         style = "text-decoration: none; font-size: 20px; margin-right: 80px; margin-top: 20px; text-align: center; border-color: #0275d8; border-radius: 40px; width: 200px; box-shadow: 0px 0px 15px 15px #0275d8;")))
)

################################################################################
################################### FOOTER #####################################
################################################################################

app_footer <- flexPanel(
    id = "footer",
    align_items = "center",
    justify_content = 'space-evenly',
    style = 'border-radius: 20px; background-color: #daedf0;',
    div(
        tags$a(img(src = "UCPH_logo.png", style = "margin-left: 30px; width: 300px; cursor:pointer;"),href="https://in.ku.dk/research/madsen-lab/")#,
        #tags$a(img(src = "LF_logo.png",   style = "margin-left: 30px; width: 300px; cursor:pointer;"),href="https://lundbeckfonden.com/")
        ),
    div('If you use this database in a publication, please cite: Rombach, J., Nielsen, T. T. E., .....')
)

app_footer_data <- flexPanel(
  id = "footerdata",
  align_items = "center",
  justify_content = 'space-evenly',
  style = 'border-radius: 20px; margin-bottom: 30px; background-color: #fff;',
  div(
    tags$a(img(src = "UCPH_logo.png", style = "margin-left: 30px; width: 300px; cursor:pointer;"),href="https://in.ku.dk/research/madsen-lab/")#,
    #tags$a(img(src = "LF_logo.png",   style = "margin-left: 30px; width: 300px; cursor:pointer;"),href="https://lundbeckfonden.com/")
  ),
  div('If you use this database in a publication, please cite: Rombach, J., Nielsen, T. T. E., .....')
)

################################################################################
################################ FRONT PAGE ####################################
################################################################################

home_page <- gridPanel(
        id = 'home_page',
        template = my_template, #'holy-grail',
        rows = '1fr 3fr, 1fr',
        
        header = app_header_front,

        main = flexPanel(
                id = "main",
                align_items = "center",
                gap = 20,
                justify_content = 'space-evenly',
                flex = c(1),
                img(src = "vglut1.png", style = 'width: 350px;'),
                div(search_bar)),

        footer = app_footer
)

################################################################################
################################# DATA PAGE ####################################
################################################################################

data_page <- gridPanel(
    
    id = 'data',
    template = my_template,
    columns = "1fr 1fr 2fr",
    rows = "0.5fr 4fr 0.5fr",
    
    header = app_header,
    
    left = fluidRow(id = 'datapanel',
                    column(12,
                           dataTableOutput('proteintable')
                           ),
                    column(12, 
                           checkboxGroupInput("terminal", 
                                              "Select terminals to display", 
                                              c('C','N'), 
                                              inline = T, 
                                              selected = c('C','N'))
                    ),
                    column(6, 
                           radioButtons("confidence", 
                                        "Set confidence", 
                                        c('low','medium','high'), 
                                        inline = T, 
                                        selected = c('medium'))
                    ),
                    column(6, 
                           uiOutput("UniProtLink_button")
                           )
                    ),
    
    main  =  flexPanel(r3dmolOutput('pdb', width = 'auto', height = 'auto')),
    
    right =  fluidRow(id = 'datapanel',
                      column(12, 
                             plotOutput('trace')
                             ),
                      column(12, 
                             dataTableOutput('siteTable')
                             ),
                      column(6, 
                             plotOutput('helicalwheel')
                             ),
                      column(6,
                             actionButton(inputId = "analysis", 
                                          label = icon("magnifying-glass-chart"), 
                                          style = 'background-color: #0275d8; border-color: #0275d8; font-color: #fff; font-size: 30px; border-radius: 40px; margin-top: 60px; width: 60px; height: 60px; box-shadow: 0px 0px 15px 15px #0275d8;',
                                          onclick = paste0("location.href='",route_link("analysis"),"';")
                             )
                             )
                      ),
    
    footer = app_footer_data
)

################################################################################
############################### ANALYSIS PAGE ##################################
################################################################################

analysis_page <- gridPanel(
  
  id = 'analysis',
  template = my_template,
  columns = "1fr 1fr 2fr",
  rows = "0.5fr 4fr 0.5fr",
  
  header = app_header,
  
  left = fluidRow(),
  
  main  =  flexPanel(),
  
  right =  fluidRow(),
  
  footer = app_footer_data
  
)

################################################################################
################################ ABOUT PAGE ####################################
################################################################################

about_page <- gridPanel(
                    id = 'about',
                    template = my_template, #'holy-grail',
                    rows = '2fr 4fr, 1fr',
                    
                    header = app_header_about,
                    
                    main = flexPanel(
                                direction = 'row',
                                div(Text('What is this all about?')),
                                style = 'border-radius: 20px; box-shadow: 0 0 10px #0275d8; text-align: center;'
                        ),
                    
                    footer = app_footer
)

################################################################################
##################################### UI #######################################
################################################################################

ui <- fluidPage(
    #tags$hr(),
    tags$head(
        tags$link(href = "style.css", rel = "stylesheet", type = "text/css"),
        tags$style(type="text/css", "#proteinquery {border-radius: 30px; height: 50px; width: 450px; text-align: center; font-size: 20px;}"),
        tags$style(type="text/css", "#footer {opacity: 0.5;}"),
        tags$style(type="text/css", "#footer:hover{opacity: 1; box-shadow: 0 0 10px #0275d8;}"),
        tags$style(type="text/css", "#footerdata {opacity: 0.5;}"),
        tags$style(type="text/css", "#footerdata:hover{opacity: 1}"),
        tags$style(type="text/css", "#submit .fa-magnifying-glass {color:#fff}"),
        tags$style(type="text/css", "#analysis .fa-magnifying-glass-chart {color:#fff}")
    ),
    
    router_ui(
        route("/",     home_page),
        route("data",  data_page),
        route("analysis", analysis_page),
        route("about", about_page))
)

################################################################################
################################### SERVER #####################################
################################################################################

server <- function(input, output, session) {
  
    
  
    # Search for query matches -------------------------------------------------
    selectprotein <- eventReactive(input$submit, {
        search_for_protein(search = input$proteinquery, keys = keys)
    })
    
    # Plot table with hits for search query ------------------------------------
    output$proteintable = DT::renderDT(selectprotein()[,c(3,1,2)],
                                       plugins = "ellipsis",
                                       selection = list(
                                                        mode='single',
                                                        selected=1
                                                        ),
                                       options=list(pageLength=3, 
                                                    dom = 'tp',
                                                    # limit cells in columns 1 and 2 to 17 characters
                                                    columnDefs = list(
                                                        list(
                                                        targets = c(2,3),
                                                        render = JS("$.fn.dataTable.render.ellipsis( 14, false )")
                                                        ),
                                                        list(
                                                            width = '100px', 
                                                            targets = c(1)
                                                            )
                                                        )
                                                    )
                                       )

    # Import Array_object ------------------------------------------------------
    Array_object_path <- reactive({selectprotein()[input$proteintable_rows_selected,5]})
    
    Array_object <- reactive({readRDS(Array_object_path())})
    
    # Select terminal to focus on ----------------------------------------------
    UniProt <- reactive({
      
          if(!is.null(input$proteintable_rows_selected)){
                if(length(input$terminal)==1){
                   UniProt <- paste0(selectprotein()[input$proteintable_rows_selected,3],'_',input$terminal)
                } else {
                   UniProt <- selectprotein()[input$proteintable_rows_selected,3]
                }
          } else {
                if(length(input$terminal)==1){
                   UniProt <- paste0(selectprotein()[1,3],'_',input$terminal)
                } else {
                   UniProt <- selectprotein()[1,3]
                }
          }
      })
    
    # Create button that links to UniProt --------------------------------------
    output$UniProtLink_button <- renderUI({
      shiny::a(
        h4(img(src='UniProtLogo.png', width='100px'),
           class = "btn btn-default action-button",
           style = "width: 120px; height 60px; border-color: #fff;"),
        target = "_blank",
        href = paste0("https://www.uniprot.org/uniprotkb/", selectprotein()[input$proteintable_rows_selected,3], "/entry")
      )
    })
    
    # Set confidence threshold -------------------------------------------------
    conf <- reactive({input$confidence})
    
    # Plot structure with SWaFi score ------------------------------------------
    output$pdb <- renderR3dmol({
      
      pdb <- quiet(plot_pdb_with_binding_sites(UniProt(),
                                         array_object=Array_object(), 
                                         pthresh = conf(), 
                                         zoom_limits = c(5,1000)))
      
      if(!is.null(input$siteTable_rows_selected)){
        # highlight selected site
        highlight <- as.vector(unlist(seq2(from = sites()[input$siteTable_rows_selected,1], to = sites()[input$siteTable_rows_selected,2], by = 1)))
        pdb <- pdb %>% m_set_style(style = m_style_cartoon(color = 'red'), sel = m_sel(resi = highlight))
      }
      
      return(pdb)
      
    })
    
    # Plot SWaFi trace ---------------------------------------------------------
    selecttrace <- reactive(plot_SWaFi_trace_w_sites(array_object = Array_object(), 
                                                     proteins=UniProt(), 
                                                     dataset = c('FBBE1','FBBE2'), 
                                                     confidence = conf(), 
                                                     nCol=2))
    
    output$trace <- renderPlot(selecttrace(), height = 300)
    
    # Get binding sites --------------------------------------------------------
    sites <- reactive({
      quiet(get_binding_sites(proteins=UniProt(),
                        #binding = T,
                        confidence = conf(),
                        array_object = Array_object())[,c(3:6,8,10)])  %>% 
        dplyr::mutate(across(where(is.numeric), round, 5))
    })
    
    # Create binding site table ------------------------------------------------
    output$siteTable = DT::renderDT(sites(),
                                       plugins = "ellipsis",
                                       selection = list(
                                         mode='single'
                                       ),
                                       options=list(pageLength=3, 
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
    
    # Select site to focus on --------------------------------------------------
    activeSite <- reactive({
                        if(!is.null(input$siteTable_rows_selected)){
                          activeSite <- sites()[input$siteTable_rows_selected,6]
                        } else {
                          activeSite <- NA
                        }
                  })
      
    
    # Get best AH --------------------------------------------------------------
    selectedsequence <- reactive({
                            if(!is.na(activeSite())){
                              selectedsequence <- quiet(get_best_AH(activeSite(),UniProt(), array_object = Array_object(), return_all_info = F))
                            } else {
                              selectedsequence <- 'not_selected'
                            }
                        })
    
    # Plot helical wheel of AH -------------------------------------------------
    output$helicalwheel <- renderPlot({
      
      
      if(!is.na(selectedsequence())){
          # plot helical wheel
          if(selectedsequence()=='not_selected'){
            helicalwheel <- ggplot() + 
              annotate(geom="text", x=0, y=0, label='No site selected', size = unit(12, "pt")) + 
              prism() + remove_x_axis() + remove_y_axis() + theme(legend.position='none') +
              xlim(-1.5,1.5) +
              ylim(-1.5,1.5)
          } else {
            helicalwheel <- heliquest(sequence = selectedsequence())
          }
            
       } else {
          # plot empty graph with message
         helicalwheel <- ggplot() + 
                    annotate(geom="text", x=0, y=0, label='No amphipathic\nhelix found', size = unit(12, "pt")) + 
                    prism() + remove_x_axis() + remove_y_axis() + theme(legend.position='none') +
                    xlim(-1.5,1.5) +
                    ylim(-1.5,1.5)
       }
      
       return(helicalwheel)
      })
    
    # setup multipage app ------------------------------------------------------
    router_server()    
    
}

################################################################################
##################################### RUN ######################################
################################################################################

shinyApp(ui = ui, server = server)