library(shiny)
#call for functions calculating the poisson score
source('./modules/poiss_calc.R')
#call for function for validation of the input
source('./modules/validateinput.R')


shinyServer( function(input, output, session){
  #defining constants
  #symbol for lambda
  lambda <- intToUtf8(955)
  pc_global <<- NULL
  
  #functionality for instructions button
  observeEvent(input$instructions, {
    updateNavbarPage(session, "nav-Page",selected = "How to use graPe")
  })
  
  #functionality for about graPe button
  observeEvent(input$about, {
    updateNavbarPage(session, "nav-Page",selected = "About")
  })
  
  #functionality for citing graPe button
  observeEvent(input$cite, {
    updateNavbarPage(session, "nav-Page",selected = "Citing graPe")
  })
  
  #functionality for submit button
  observeEvent(input$do, {
    show(selector = "#nav-Page li a[data-value=Output]")
    updateNavbarPage(session, "nav-Page",selected = "Output")
  })
  
  #functionality for the reset button
  observeEvent(input$reset, {
    reset('file1')
    disable("do")
    hide('choose_negcon')
    hide('choose_poscon')
    hide(selector = "#nav-Page li a[data-value=Output]")
  })
  
  #activate the reset button when a file is uploaded
  observeEvent(input$file1,
               if(!is.null(input$file1)){
                 enable("reset")
                 show('choose_negcon')
                 show("choose_poscon")
                 hide(selector = "#nav-Page li a[data-value=Output]")
               }, 
               ignoreNULL = TRUE
  )
  
  #get the file name
  file_name <- reactive({
    inFile <- input$file1
    if (is.null(inFile))
      return(NULL)
    return (inFile$name)
  })
  output$file1.name <- renderText({ file_name() })
  
  #get treatment names from the uploaded data
  trt_lst <- eventReactive(input$file1,{
    if (!is.null(input$file1)) {
      inFile <- input$file1
      validate(
        need(file_ext(inFile$name) %in% c(
          'csv'), "Wrong File Format try again!"))
      in.dat <- read.csv(inFile$datapath)
      cols.expected = c("run_id","plate_id","well_id","site_id", "trt","total","count" )
      validate(
        need(((setequal(colnames(in.dat),cols.expected))), 
             "Input file header row must conform to the Instructions."))
      in.dat <- read.csv(inFile$datapath)
      trt_unique_val <- in.dat %>% 
        group_by(plate_id, well_id,trt) %>% summarise(freq_site = n()) %>% 
        group_by(trt) %>% summarize(freq =  n()) %>% 
        filter(freq > 1) %>% arrange(desc(freq))
      return(unique(trt_unique_val$trt))
    }
    else
      return(NULL)
  })
  
  #render the negcon dropdown when a file is uploaded
  output$choose_negcon <- renderUI({
    req(input$file1)
    selectInput("negcon",label = "Select negative control",choices=c("NONE",as.character(unlist(trt_lst()))),selected="NONE")
  })
  
  #render the poscon dropdown and disable/enable 'Go' button when the negcon input is selected
  output$choose_poscon <- renderUI({
    req(input$file1)
    req(input$negcon)
    if(input$negcon != 'NONE'){
      negcon_val<-input$negcon
      choice_poscon<-setdiff(trt_lst(),negcon_val)
      selectInput("poscon",label = "Select positive control",choices=c("NONE",choice_poscon),selected="NONE")
    }
    else{
      #error.flag <<- 0
      #new.input <<- 1
      disable("do")
      return(NULL)
    }
  })
  
  #enable/disable the 'Go' button if the poscon input changes
  observeEvent(input$poscon,{
    req(input$file1)
    #new.input <<- 1
    if((input$negcon != 'NONE') & (input$poscon != 'NONE')){
      enable("do")
    }
    else{
      #error.flag <<- 0
      disable("do")
    }
  },priority = 1)
  
  #show/hide the output tab when the Submit button is clicked
  observe({
    hide(selector = "#nav-Page li a[data-value=Output]")
  })
  
  #functionality for the reset button in output tab
  observeEvent(input$restart, {
    reset('file1')
    disable("do")
    hide('choose_negcon')
    hide('choose_poscon')
    hide(selector = "#nav-Page li a[data-value=Output]")
    updateNavbarPage(session, "nav-Page",selected = "graPe")
  })
  
  #functionality for the reset button in output tab
  observeEvent(input$update, {
    updateNavbarPage(session, "nav-Page",selected = "graPe")
  })
  
  #error checking when the 'Go' button is pressed
  calculate.scores <- eventReactive(input$do, {
    inFile <- input$file1
    dat = read.csv(inFile$datapath)
    cols.expected = c("run_id","plate_id","well_id","site_id", "trt","total","count" )
    if((setequal(colnames(dat),cols.expected)) & (length(colnames(dat)) == length(cols.expected))) {
      nc = input$negcon
      pc = input$poscon
      pc_global <<- input$poscon
      if(!(nc %in% dat$trt)){
        return("Negative control name must occur in the input data.")
      }
      if(!(pc %in% dat$trt)){
        return("Positive control name must occur in the input data.")
      }
      if(!(is.numeric(dat$site_id) & is.numeric(dat$total) & is.numeric(dat$count))) {
        return("Numeric data columns must not contain alphanumeric data.")
      }
      #make sure each plate has at least one DMSO well and collect stats
      unique_plate_id = unique(dat$plate_id) #set of unique plates
      if(length(intersect(unique_plate_id,unique(dat[dat$trt==nc,'plate_id']))) != length(unique_plate_id)) {
        return(paste0("Negative control wells must exist on plate_id: ",setdiff(unique_plate_id,unique(dat[dat$trt==nc,'plate_id']))))
      }
      if(length(intersect(unique_plate_id,unique(dat[dat$trt==pc,'plate_id']))) != length(unique_plate_id)) {
        return(paste0("Positive control wells must exist on plate_id: ",setdiff(unique_plate_id,unique(dat[dat$trt==pc,'plate_id']))))
      }
      else {
        #calculate the scores if there are no errors
        compound.scores <- run_graPe(dat,nc)
        #compound.scores[,-c(1,2)] = round(compound.scores[,-c(1,2)],2)
        return(data.frame(compound.scores))
      }
    } else {
      return("Input file header row must conform to the Instructions.")
    }
  })
  
  #read the file and display the contents in the Uploaded Data tab
  output$contents <- DT::renderDataTable({
    req(input$file1)
    return(DT::datatable(read.csv(input$file1$datapath)))
  }, 
  #extensions = c('ColVis', 'ColReorder'),
  #options = list(dom = 'C<"clear">Rlfrtip',
  #               colReorder = list(realtime = TRUE))
  )
  
  #render data if any of the inputs change or the 'Go' button is pressed
  observeEvent({input$file1
    input$negcon
    input$poscon
    input$do
  },
  {  #recreate the error message-ui
    output$err.msg <- renderText({
      #scoredata <- calculate.scores()
      inFile <- input$file1
      if(file_ext(inFile$name) %in% c('csv')){
        return(NULL)
      }
      else{
        disable("do")
        hide(selector = "#nav-Page li a[data-value=Output]")
        return(NULL) 
      } 
    })
    
    #output the calculated Poisson z-prime(plate quality tab)
    output$qual <- DT::renderDataTable({
      scoredata <- calculate.scores()
      positive.control = pc_global
      qual_output = scoredata[(scoredata$trt == positive.control),
                              c('run_id',
                                'plate_id',
                                'trt',
                                'zprime_poiss',
                                'nc.lambdahat.per.plate',
                                'nc.upperqi.per.plate.spline',
                                'trt.lambdahat.per.plate',
                                'trt.lowerqi.per.plate.spline')]
      float_cols <- c('zprime_poiss',
                      'nc.lambdahat.per.plate',
                      'nc.upperqi.per.plate.spline',
                      'trt.lambdahat.per.plate',
                      'trt.lowerqi.per.plate.spline')
      qual_output[,float_cols] <- round(qual_output[,float_cols],2)
      colnames(qual_output) = c('run_id','plate_id','trt','poiss_zp',paste0('nc_',lambda,'hat'),'nc_upperci',paste0('pc_',lambda,'hat'),'pc_lowerci')
      rownames(qual_output) = seq(1,nrow(qual_output))
      return(DT::datatable(qual_output))
      #}
    }, 
    #extensions = c('ColVis', 'ColReorder'),
    options = list(order = list(list(1, 'asc')))
    )
    
    
    #output the calculated poisson d-score(treatment score tab)
    output$dval <- DT::renderDataTable({
      scoredata <- calculate.scores()
      dscore_output = unique(scoredata[,c('run_id','trt','dscore_poiss','fold_change','effect_size','trt.lambdahat.per.run','trt.lowerci.per.run','nc.lambdahat.per.run','nc.upperci.per.run')])
      float_cols <- c('dscore_poiss',
                      'fold_change',
                      'effect_size',
                      'trt.lambdahat.per.run',
                      'trt.lowerci.per.run',
                      'nc.lambdahat.per.run',
                      'nc.upperci.per.run')
      dscore_output[,float_cols] <- round(dscore_output[,float_cols],2)
      colnames(dscore_output) = c('run_id','trt','poiss_ds','fold_chg','eff_siz',paste0('trt_',lambda,'hat'),'trt_lowerci',paste0('nc_',lambda,'hat'),'nc_upperci')
      rownames(dscore_output) <- seq(length=nrow(dscore_output)) 
      return(DT::datatable(dscore_output))
      #}
    }, 
    options = list(order = list(list(2, 'desc')))
    )
    
    
    #render the download button only when there is no error (plate quality tab)
    output$qual_download_button <- renderUI({
      scoredata <- calculate.scores()
      #if (new.file == 1){
      #  return(NULL)
      #}
      req(input$do)
      #if (download.button.flag == 1) {
      downloadButton('qual_downloadData', 'Download', class = "btn-primary")
      #}
    })
    
    #show legend for qual panel
    observeEvent(input$showlegend_qualpanel, {
      showModal(modalDialog(
        title = "Legend",
        tags$ul(style="list-style-type:circle;padding-left: 12px;",
                tags$li("run_id = user-provided run identifier"),
                tags$li("plate_id = user-provided plate identifier"),
                tags$li("trt = treatment condition of the well (e.g., DMSO or a compound name)"),
                tags$li("poiss_zp = Poisson-based Z' factor for the plate based on user-provided definitions of control treatments"),
                tags$li(paste0('nc_',lambda,'hat = ',lambda,' hat (expectation value) of the Poisson distribution fit to user-provided negative-control wells for the plate')),
                tags$li("nc_upperci = upper 95% confidence boundary of the expectation value for negative-control wells"),
                tags$li(paste0('pc_',lambda,'hat = ',lambda,' hat (expectation value) of the Poisson distribution fit to user-provided positive-control wells for the plate')),
                tags$li("pc_lowerci = lower 95% confidence boundary of the expectation value for positive-control wells")
        ),
        easyClose = TRUE
      ))
    })
    
    #show legend for score panel
    observeEvent(input$showlegend_scorepanel, {
      showModal(modalDialog(
        title = "Legend",
        tags$ul(style="list-style-type:circle;padding-left: 12px;",
                tags$li("run_id = user-provided run identifier"),
                tags$li("trt = treatment condition of the well (e.g., DMSO or a compound name)"),
                tags$li("poiss_ds = Poisson-based D-score for the treatment the compound (cf. https://www.ncbi.nlm.nih.gov/pubmed/24464433/)"),
                tags$li("fold_chg = fold-change (ratio) of the treatment response compared to the negative control"),
                tags$li("eff_siz = effect size (difference) between the treatment response and the negative control"),
                tags$li(paste0('trt_',lambda,'hat = ',lambda,' hat (expectation value) of the Poisson distribution fit to treatment wells')),
                tags$li("trt_lowerci = lower 95% confidence boundary of the expectation value for treatment wells"),
                tags$li("trt_upperci = upper 95% confidence boundary of the expectation value for treatment wells"),
                tags$li(paste0('nc_',lambda,'hat = ',lambda,' hat (expectation value) of the Poisson distribution fit to negative-control wells in this run')),
                tags$li("nc_upperci = upper 95% confidence boundary of the expectation value for negative-control wells")
        ),
        easyClose = TRUE
      ))
    })
    
    #render the download button only when there is no error (treament score tab)
    output$dval_download_button <- renderUI({
      scoredata <- calculate.scores()
      #if (new.file == 1){
      #  return(NULL)
      #}
      req(input$do)
      #if (download.button.flag == 1) {
      downloadButton('dval_downloadData', 'Download', class = "btn-primary")
      #}
    })
  })
  
  
  #download button functionality for plate quality tab
  output$qual_downloadData <- downloadHandler(
    filename = function() { #name for the downloaded file
      inFile <- input$file1
      paste0(tools::file_path_sans_ext(inFile$name),"_plate_quality_",format(Sys.time(), "%Y%m%d_%H%M"), ".csv", sep = "")
    },
    content = function(file) { #content of downloaded file
      scoredata <- calculate.scores()
      positive.control = input$poscon
      qual_output = scoredata[(scoredata$trt == positive.control),c(2,1,7,5,6,3,4)]
      colnames(qual_output) = c('plate_id','trt','poiss_zp','nc_lambdahat','nc_upperci','pc_lambdahat','pc_lowerci')
      write.csv(qual_output[order(qual_output$plate_id),], file, row.names = FALSE)
      
    }
  )
  
  #download button functionality for treatment score tab
  output$dval_downloadData <- downloadHandler(
    filename = function() { #name for the downloaded file
      inFile <- input$file1
      paste0(tools::file_path_sans_ext(inFile$name),"_treatment_scores_",format(Sys.time(), "%Y%m%d_%H%M"), ".csv", sep = "")
    },
    content = function(file) { #content of downloaded file
      scoredata <- calculate.scores()
      dscore_output = unique(scoredata[,c('trt','dscore_poiss','fold_change','effect_size','trt.lambdahat.per.run','trt.lowerci.per.run','nc.lambdahat.per.run','nc.upperci.per.run')])
      colnames(dscore_output) = c('trt','poiss_ds','fold_chg','eff_siz','trt_lambdahat','trt_lowerci','nc_lambdahat','nc_upperci')
      write.csv(dscore_output[order(dscore_output$poiss_ds,decreasing = TRUE),], file, row.names = FALSE)
    }
  )
  
}
)
