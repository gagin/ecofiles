library(shiny)

shinyUI(fluidPage(
        titlePanel("Uploading Files"),
        sidebarLayout(
                sidebarPanel(
                        fileInput('file1', 'Measurements TXT (TSV) File',
                                  accept=c('text/txt', 
                                           'text/tab-separated-values,text/plain', 
                                           '.tsv')),
                        tags$hr(),
                        fileInput('file2', 'Samples mapping CSV File',
                                  accept=c('text/csv', 
                                           'text/comma-separated-values,text/plain', 
                                           '.csv'))
                        
                ),
                mainPanel(
                        p("Results"),
                        tableOutput("contents")
                )
        )
))