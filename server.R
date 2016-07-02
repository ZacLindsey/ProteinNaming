library(shiny)
library(datasets)
library(ggplot2) # load ggplot
library(seqinr)
ClusterGroups <- read.csv(file="ClusterGroups_6-28.csv",sep=",")
# Define server logic required to plot various variables against mpg
shinyServer(function(input, output) {
  
  # Compute the forumla text in a reactive function since it is 
  # shared by the output$caption and output$mpgPlot functions
  formulaText <- reactive(function() {
    paste("mpg ~", input$variable)
  })
  
  # Return the formula text for printing as a caption
  output$caption <- renderText(function() {
    formulaText()
  })
  
  # Generate a plot
  
  getHeight<-function() # NOTE: I don't know how/why this works!
  {
    print(exprToFunction(input$rC2))
    return(exprToFunction(input$rC2))
  }
  
  output$proteinPlot <- renderPlot({
    # check for the input variable
    backbone <- input$variable
    
    x=c()
    for (id in unique(ClusterGroups$V1))
    {
      x<-c(x,(nrow(ClusterGroups[which(ClusterGroups$V1==id),])))
    }
    maxx <- max(x)
    
    
      CG2 <- ClusterGroups[which(ClusterGroups$V1==backbone),]
      cols <- rainbow(length(unique(CG2$V2)))
      CG2$V6 <- c(1:nrow(CG2))
      rC<<-nrow(ClusterGroups)
      
      index = 1
      for (i in unique(CG2$V2))
      {
        CG2$V6[which(CG2$V2==i)] <- cols[index]
        index = index + 1
      }

      # Chunk to ensure correct data formatting
      # These few lines may be unnecessary
      CG2$V5[which(CG2$V5=="*")] <- "100.0"
      CG2$V5<-as.double(as.character(CG2$V5))
      CG2$V5[which(is.na(CG2$V5))] <- 0.0
      CG2$V4 <- as.double(CG2$V4)
      CG2$V5 <- CG2$V5/100.00

      CG3 <- CG2[order(CG2$V4),] # CRITICAL: Allow the reordering of colors!
      
      CG3$V7<-c(1:nrow(CG3)) # CRITICAL: Used for labeling and position!
      # Will be the same as V7 if no subset

      if (input$Subset) # For subsetting only: modify CG3
      {
        CG3<-subsetData(CG3,input$ranges)
      }
      
      #CG3$V8<-c(1:nrow(CG3)) # CRITICAL: Used for labeling and position!
      
      seqData <<- as.character(CG3$V3)
      
      print(ggplot(CG3,aes(x=c(1:nrow(CG3)),y=CG3$V4))+ # Notice X position is determined by V7
               geom_bar(stat="identity",
                        fill=CG3$V6,
                        alpha=CG3$V5
               )+
               scale_y_continuous()+
               scale_x_discrete()+ # Remove extra space
               geom_text(size=3.75,y=max(CG3$V4)/2.0,label=CG3$V3,color="black") +
               geom_text(size=3.5,y=-15,
                         aes(label=rev(CG3$V7)),
                         color="black") +
              coord_flip() +
               ggtitle(backbone)+
               theme(axis.text.y=element_blank(),
                     axis.ticks.y=element_blank(),
                     axis.title.y=element_blank()))
  },
  height=getHeight() # WEIRD!!!
  )
  
  output$downloadData <- downloadHandler(
    
    # This function returns a string which tells the client
    # browser what name to use when saving the file.
    filename = function() {
      paste(input$variable, "faa", sep = ".")
    },
    
    # This function should write data to a file given to it by
    # the argument 'file'.
    content = function(file) {

      #plasmidAA <- read.fasta("Plasmids20-200kb-6-9-2016AA.fa","AA")
      print("FINISHED FASTA READ")
      
      writeSeqs <- function(seq)
      {
        seqName <- attr(seq,"name")
        if (grepl("<unknown description>",seqName))
        {
          seqName <- strsplit(seqName," ")[[1]][1]
        }
        if (tolower(seqName) %in% seqData)
        {
          print("RETURNING...")
          return(seq)
          print("DONE WITH A RETURN")
        }
      }
      Seqs <<- lapply(plasmidAA,writeSeqs)
      Seqs<-Seqs[!sapply(Seqs, is.null)] # REMOVE NULL seqs 
                                         # (since lapply return same length)
      print("DONe")
      write.fasta(Seqs,input$variable,
                  file.out=paste(input$variable, "faa", sep = "."))
    }
  )
})