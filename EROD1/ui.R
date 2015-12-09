library(shiny)

shinyUI(fluidPage(
        h1("Uploading Files"),
        fluidRow(
                column(6,
                        fileInput('file1', 'Measurements TXT (TSV) File',
                                  accept=c('text/txt', 
                                           'text/tab-separated-values,text/plain', 
                                           '.tsv'))
                       ),
                column(6,
                        fileInput('file2', 'Samples mapping CSV File',
                                  accept=c('text/csv', 
                                           'text/comma-separated-values,text/plain', 
                                           '.csv'))
                       )
                ),
                        h2("Results"),
                        tabsetPanel(
                                tabPanel("Slopes", tableOutput("contents")), 
                                tabPanel("Protein", plotOutput("p1")), 
                                tabPanel("Resorufin", plotOutput("p2")),
                                tabPanel("Slopes same scale", plotOutput("p3")),
                                tabPanel("Slopes diff scales", plotOutput("p4"))
                        )
                )
        )