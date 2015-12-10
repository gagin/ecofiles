library(shiny)

shinyUI(fluidPage(
        h1("Uploading Files"),
        fluidRow(
                column(4,
                       tags$a(href = "https://raw.githubusercontent.com/gagin/ecofiles/master/2015-11-26%20EROD%20pr%20.txt", "Data example"),
                       fileInput('file1', 'Measurements TXT (TSV) File',
                                 accept=c('text/txt', 
                                          'text/tab-separated-values,text/plain', 
                                          '.tsv'))
                ),
                column(4,
                       tags$a(href = "https://raw.githubusercontent.com/gagin/ecofiles/master/default.csv", "Mapping example"),
                       fileInput('file2', 'Sample mapping CSV File',
                                 accept=c('text/csv', 
                                          'text/comma-separated-values,text/plain', 
                                          '.csv'))

                ),
                column(4,
                       radioButtons('inGroup', 'Columns per sample group', 
                                    choices = 2:3, inline=TRUE)
                )
        ),
        fluidRow(
                column(5,
                       textInput('proteins', 'Protein calibration',
                                 value = '0, 0.1, 0.3, 0.75, 1.1, 1.4')
                ),
                column(5,
                       textInput('rs', 'Resorufin calibration',
                                 value = '0, 0.1, 0.2, 0.3, 0.4, 0.5')
                ),
                column(2,
                       submitButton("Recalculate")
                )
        ),
        h1("Results"),
        tabsetPanel(
                tabPanel("Slopes", 
                         fluidRow(
                                 column(6, h2("Slopes by groups"), tableOutput("slopes.narrow")),
                                 column(6, h2("Based on all points per sample"), tableOutput("slopes"))
                         )),
                tabPanel("Slopes grid", tableOutput("slopes.grid")),
                tabPanel("Protein", plotOutput("p1")), 
                tabPanel("Resorufin", plotOutput("p2")),
                tabPanel("Charts same scale", plotOutput("p3")),
                tabPanel("Charts separate scales", plotOutput("p4")),
                tabPanel("Info", 
                         h1("What is this?"),
                         p("This script automates operations on aquatic toxicity measurements as performed previously in Excel like this:"),
                         tags$a(href = "https://raw.githubusercontent.com/gagin/ecofiles/master/2015-11-26%20EROD%20calc.xlsx", "Original calculation file"),
                         h2("Assumptions"),
                         p("Names of samples are uploaded as a CSV table, and the program will try to calculate slope for every non-empty name. Topright 6x4 cells should be kept empty though, as they are used for calibration."),
                         p("Please also notice that in the current implementation protein calibration is done over averaged values for each time mark. When approximation is build over indepenedent points, result is different - at least the intercept is."),
                         h3("Warnings"),
                         p("Please note that case of three columns per group hasn't been tested."),
                         p("Groups per sample is only used for the first table, everything else uses sample name matching.")
                         )
        )
)
)