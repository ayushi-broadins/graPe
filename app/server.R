library(shiny)
#call for functions calculating the poisson score
source('./poiss.calc.R')
#call for function for validation of the input
source('./validateinput.R')


shinyServer( function(input, output, session){
  #defining constants
  #symbol for lambda
  lambda <- intToUtf8(955)
  
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
  
}
)