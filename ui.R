library(shiny)
ClusterGroups <<- read.csv(file="ClusterGroups.csv",head=FALSE,sep=",")

shinyUI(pageWithSidebar(
  
  headerPanel("Plasmid Backbone"),
  
  sidebarPanel(
    selectInput("variable", "Backbone Protein Name:",
                choices=names(
                  setNames(
                    unique(
                      ClusterGroups$V1),
                    unique(
                      ClusterGroups$V1)
                    )
                  )
                ),
    
    radioButtons("arrange","Order By:",c("Sequence Length","AA Cluster")),
    textInput("ranges",label="Select Ranges (1-5, 10-20, ...):"),
    
    checkboxInput("Subset", "Only Select These Ranges", FALSE)
    
    # To use the local download feature, uncomment the following line
    # AND add a comma after checkboxInput
    #downloadButton('downloadData', 'Download (may take ~5 min)')
  ),
  
  mainPanel(
    h3("Plot of this protein (may take a few min to render)"),
    h5("Bold labeled font = matched a reference plasmid"),
    plotOutput("proteinPlot",width="100%")
  )
))
