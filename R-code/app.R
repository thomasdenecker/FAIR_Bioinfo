################################################################################
# FAIR_app
# Thomas DENECKER
# 01 /2019
#
# GitHub :
# https://github.com/thomasdenecker/FAIR_Bioinfo
################################################################################

################################################################################
# Library
################################################################################

library(shiny)
library(shinydashboard)

################################################################################
# UI
################################################################################

header <- dashboardHeader(title = "FAIR_app")
sidebar <- dashboardSidebar(
  sidebarMenu(
    menuItem("Import files", tabName = "import", icon = icon("th")),
    menuItem("Description", tabName = "description", icon = icon("th")),
    menuItem("Data exploration", tabName = "dataExplo", icon = icon("th")),
    menuItem("Normalization", tabName = "normalization", icon = icon("th")),
    menuItem("Differential analysis", tabName = "diffAna", icon = icon("th")),
    menuItem("Session information", tabName = "session", icon = icon("cubes")),
    menuItem("Bibliography", tabName = "bibli", icon = icon("book"))
  )
)

body <- dashboardBody(
  tabItems(
    # First tab content
    tabItem(tabName = "import",
            fluidRow(
              box(plotOutput("plot1", height = 250)),
              
              box(
                title = "Controls",
                sliderInput("slider", "Number of observations:", 1, 100, 50)
              )
            )
    ),
    
    # Second tab content
    tabItem(tabName = "description",
            h2("Widgets tab content")
    ), 
    
    # Second tab content
    tabItem(tabName = "session",
            h2("Session"),
            p("The versions of the R software and Bioconductor packages used for this analysis are listed below. 
              It is important to save them if one wants to re-perform the analysis in the same conditions."),
            uiOutput("sessionText")
    )
  )
)


ui <- dashboardPage(skin= "red", header, sidebar, body)

server <- function(input, output, session) {
  si = sessionInfo()
  output$sessionText = renderUI({
    HTML(paste(si[[1]]$version.string,",", si[[1]]$platform, "<br>","<br>",
               "Locale : ", si[[3]], "<br>","<br>",
               "Attached base packages : ", paste(si[[5]], collapse = " , "),"<br>","<br>",
               "Other attached packages :", paste(unlist(lapply(si$otherPkgs, function(x){paste(x$Package, x$Version)})), collapse = " , "), "<br>","<br>",
               "Loaded via a namespace (and not attached) :" ,paste(unlist(lapply(si$loadedOnly, function(x){paste(x$Package, x$Version)})), collapse = " , ") 
          ))
  })
}

shinyApp(ui, server)