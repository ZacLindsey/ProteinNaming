library(shiny)

# Define UI for miles per gallon application
shinyUI(pageWithSidebar(
  
  # Application title
  headerPanel("Plasmid Backbone"),
  
  # Sidebar with controls to select the variable to plot against mpg
  # and to specify whether outliers should be included
  sidebarPanel(
    selectInput("variable", "Variable:",
                choices=names(
                  setNames(
                    unique(
                      ClusterGroups$V1),
                    unique(
                      ClusterGroups$V1)
                    )
                  )
                ),
    
    numericInput("rC2",label="Height (Pixels):",value=1000),
    textInput("ranges",label="Enter ranges to subset (1-5, 10-20, ...):"),
    
    checkboxInput("Subset", "Subset", FALSE),
    
    downloadButton('downloadData', 'Download (may take ~5 min)')
  ),
  
  # Show the caption and plot of the requested variable against mpg
  mainPanel(
    h3("Plot of this protein (may take a few min to render)"),
    
    plotOutput("proteinPlot",width="100%")
  )
))
