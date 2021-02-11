#appTitle:graPe - generalized rare-event assay Poisson estimator
#this app has 2 tabs - main page and help

library(shiny)
library(shinyBS)
library(shinyWidgets)
library(shinythemes)
library(shinyjs)
library(markdown)
lambda <- intToUtf8(955)

#select the theme from the predefined list of 
#themes at https://rstudio.github.io/shinythemes/

shinyUI(
  #create the navigation bar
  navbarPage(
    " ",
    id = "nav-Page",
    theme = shinytheme("flatly"),
    collapsible = TRUE,
    windowTitle = "graPe",
    #Tab #1: Home 
    tabPanel(
      "graPe",
      fluidRow(
        column(
          width = 8,
          offset = 2,
          includeHTML("ui_html_text/home.txt"),
          column(
            width = 12,
            #offset = 2,
            #CSS style to change color of the text in validate()
            tags$head(
              tags$style(HTML("
              .shiny-output-error-validation {
              color: #ff0000;
              font-weight: bold;
              }"))
            ),
            sidebarPanel(
              style = "margin: 20px 0px 100px;",
              width = 12,
              actionButton("instructions", 
                           HTML("&nbsp;&nbsp;Instructions"), 
                           icon = icon("chalkboard-teacher"),
                           class = "btn-primary",
                           width = '150px',
                           style='margin:5px;padding:5px; font-size:90%'),
              actionButton("about", 
                           HTML("&nbsp;&nbsp;About graPe"), 
                           icon = icon("wine-glass-alt"),
                           class = "btn-primary",
                           width = '150px',
                           style='margin:5px;padding:5px; font-size:90%'),
              actionButton("cite", 
                           HTML("&nbsp;&nbsp;Citing graPe"), 
                           icon = icon("quote-left"),
                           class = "btn-primary",
                           width = '150px',
                           style='margin:5px;padding:5px; font-size:90%'),
              br(),
              br(),
              fileInput('file1', 
                        div(
                          style = "color:black", #text-align:left;
                          strong("Analyze results using graPe  "),
                          tipify(circleButton("pB2", 
                                          "", 
                                          icon = icon("info-circle"),
                                          size = "xs"), 
                                 paste0("Analyze your results by ",
                                 "uploading a CSV file. ",
                                 "Click on instructions for more details "
                                )
                          )
                          )
                        ),
              uiOutput("choose_negcon"),
              uiOutput("choose_poscon"),
              fluidRow(
                column(12,
                       useShinyjs(),
                       disabled(actionButton("do",
                                             HTML("&nbsp;&nbsp;Submit"),
                                             width = '150px',
                                             style = "color: #FFFFF; 
                                                     background-color: #994EFF; 
                                                     border-color: #994EFF;
                                                     margin:5px;
                                                     padding:5px; 
                                                     font-size:90%")),
                       disabled(actionButton("invalidate", 
                                             HTML("&nbsp;&nbsp;Validate input"),
                                             class = "btn-warning",
                                             width = '150px',
                                             style = 'margin:5px;padding:5px; 
                                                     font-size:90%')),
                       disabled(actionButton("reset", 
                                             HTML("&nbsp;&nbsp;Reset"), 
                                             class = "btn-warning",
                                             width = '150px',
                                             style = 'margin:5px;padding:5px; 
                                                     font-size:90%'))
                )
              ),
              fluidRow(
                column(10,
                       h5(textOutput('err.msg'), align = 'left', style="color:red")
                )
              )
            )
          )
        )
      ),
      absolutePanel(
        bottom = 0,
        left = 0, 
        right = 0,
        fixed=TRUE,
        div(
          style="padding: 2px; border-top: 2px solid #CCC; background: #F7F7F7;",
          HTML(
            markdownToHTML(
              fragment.only=TRUE, 
              text=c("&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;2020 
                     The Broad Institute of MIT and Harvard")
            )
          )
        )
      )
    ),
    #the output panel
    tabPanel(
      "Output",
      fluidPage(
        tags$head(
          tags$style(
          HTML(".dataTables_wrapper { overflow-x: scroll; }" )
          )
        ),
        fluidRow(
          column(
            width = 8,
            offset = 2,
            sidebarPanel(
              style = "margin: 20px 0px 0px 0px;",
              width=12,
              fluidRow(
                column(
                  4,
                  div(
                    img(src="logo.png", height="80px"),
                    style = "padding:5px; margin:5px;"
                  ),
                  align = "left"
                ),
                column(
                  8,
                  fluidRow(
                    div(h4('Uploaded filename:', style="display:inline-block;float:left;margin-left:10px;"),
                        h4(textOutput("file1.name"), style="display:inline-block;float:left;font-weight: bold;"))
                  ),
                  fluidRow(
                    useShinyjs(),
                    actionButton("update","Update input", class = "btn-primary"),
                    actionButton("restart", "Restart analysis", class = "btn-primary"),
                    tags$style(type='text/css', "button#restart {margin-left: 10px;}"),
                    tags$style(type='text/css', "button#update {margin-left: 10px;}")
                  )
                  
                )
              ),
              
              )
          )
        ),
        fluidRow(
          column(
            value = "qualpanel", 
            width = 8,
            offset = 2,
            sidebarPanel(
              style = "border: 1px solid black; padding: 10px; margin: 20px 0px 0px 0px; background-color:#FFFFFF;",
              width=12,
              h3("Plate quality", align = "center"),
              br(),
              fluidRow(
                column(12,
                       DT::DTOutput('qual'),
                       br(),
                       #buttons
                       fluidRow(
                         column(2,
                                conditionalPanel(
                                  condition = "typeof output.qual_download_button !== 'undefined'",
                                  actionButton("showlegend_qualpanel", "Show legend", class = "btn-primary")
                                )
                         ),
                         column(3, offset = 7 , align = 'right',
                                #download button
                                uiOutput("qual_download_button")
                         )
                       )
                )
              )
            )
          )
        ),
        fluidRow(
          column(
            value = "scorepanel",
            width = 8,
            offset = 2,
            sidebarPanel(
              style='border: 1px solid black; padding: 10px; margin: 20px 0px 0px 0px; background-color:#FFFFFF;',
              width = 12,
              h3("Treatment scores", align = "center"),
              br(),
              fluidRow(
                column(12,
                       DT::DTOutput('dval'),
                       br(),
                       #buttons
                       fluidRow(
                         column(2,
                                conditionalPanel(
                                  condition = "typeof output.dval_download_button !== 'undefined'",
                                  actionButton("showlegend_scorepanel", "Show legend", class = "btn-primary")
                                )
                         ),
                         column(3, offset = 7 , align = 'right',
                                #download button
                                uiOutput("dval_download_button")
                         )
                       )
                )
              )
            )
          )
        ),
        fluidRow(
          column(
            value = "inputpanel", 
            width = 8,
            offset = 2,
            sidebarPanel(
              style='border: 1px solid black; padding: 10px; margin: 20px 0px 120px 0px; background-color:#FFFFFF;',
              width = 12,
              h3("Data uploaded", align = "center"),
              br(),
              fluidRow(
                column(12,
                       DT::DTOutput('contents')
                )
              )
            )
          )
        )
      ),
      absolutePanel(
        bottom = 0,
        left = 0, 
        right = 0,
        fixed=TRUE,
        div(
          style="padding: 2px; border-top: 2px solid #CCC; background: #F7F7F7;",
          HTML(
            markdownToHTML(
              fragment.only=TRUE, 
              text=c("&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;2020 
                     The Broad Institute of MIT and Harvard")
            )
          )
        )
      )
    ),
    #the help menu/panel
    navbarMenu(
      "Help",
      tabPanel(
        "About",
        mainPanel(
          fluidRow(
            column(
              8,
              offset = 2,
              includeHTML("ui_html_text/about.txt")
              )
            ),
          width = 12
          ),
        absolutePanel(
          bottom = 0,
          left = 0, 
          right = 0,
          fixed=TRUE,
          div(
            style="padding: 2px; border-top: 2px solid #CCC; background: #F7F7F7;",
            HTML(
              markdownToHTML(
                fragment.only=TRUE, 
                text=c("&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;2020 
                     The Broad Institute of MIT and Harvard")
              )
            )
          )
        )
      ),
      tabPanel(
        "How to use graPe",
        mainPanel(
          fluidRow(
            column(
              8,
              offset = 2,
              br(),
              tags$ul(style="list-style-type:circle;padding-left: 12px;",
                      includeHTML("ui_html_text/instructions.txt"))
              )
            ),
          width = 12
          ),
        absolutePanel(
          bottom = 0,
          left = 0, 
          right = 0,
          fixed=TRUE,
          div(
            style="padding: 2px; border-top: 2px solid #CCC; background: #F7F7F7;",
            HTML(
              markdownToHTML(
                fragment.only=TRUE, 
                text=c("&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;2020 
                     The Broad Institute of MIT and Harvard")
                )
              )
            )
          )
        ),
      tabPanel(
        "Citing graPe",
        mainPanel(
          fluidRow(
            column(
              8,
              offset =2,
              br(),
              tags$ul(style="list-style-type:circle;padding-left: 12px;",
                      includeHTML("ui_html_text/citation.txt"))
            )
          ),
          width = 12
        ),
        absolutePanel(
          bottom = 0,
          left = 0, 
          right = 0,
          fixed=TRUE,
          div(
            style="padding: 2px; border-top: 2px solid #CCC; background: #F7F7F7;",
            HTML(
              markdownToHTML(
                fragment.only=TRUE, 
                text=c("&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;2020 
                     The Broad Institute of MIT and Harvard")
              )
            )
          )
        )
      )
    ),
    tags$head(tags$style(HTML('.navbar-header {margin-left:17% !important; font-size:35px; text-align:center;}')))
  )
)
