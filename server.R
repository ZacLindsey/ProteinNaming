library(shiny)
library(datasets)
library(ggplot2)
library(seqinr)
shinyServer(function(input, output) {
  
  getHeight<-function() 
  {
      return(
        exprToFunction(
          
          if(
          (50+13*nrow(
            if(input$Subset)
            {
              subsetData(
                ClusterGroups[which(ClusterGroups$V1==input$variable),],input$ranges)
            }
            else
            {
              ClusterGroups[which(ClusterGroups$V1==input$variable),]
            }
            )
           )<32000
          )
          {
            (50+13*nrow(
              if(input$Subset)
              {
                subsetData(
                  ClusterGroups[which(ClusterGroups$V1==input$variable),],input$ranges)
              }
              else
              {
                ClusterGroups[which(ClusterGroups$V1==input$variable),]
              }
            ))
          }
          else
          {
            32000
          }
          )
        )
  }
  
  output$proteinPlot <- renderPlot({
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

      if (input$arrange == "Sequence Length")
      {
        CG3 <- CG2[order(CG2$V4),] # CRITICAL: Allow the reordering of colors!
      }
      else
      {
        CG3 <- CG2
      }
      
      CG3$V7<-c(1:nrow(CG3)) # CRITICAL: Used for labeling and position!
      # Will be the same as V7 if no subset
      
      references <- tolower(c("R7K","R46","pRA3","pKJK5","RK2","RP4","pNDM-1_Dok01"))
      cids <- c()
      for (cRow in c(1:nrow(CG3)))
      {
        plasmid <- strsplit(as.character(CG3$V3[cRow]),"_")[[1]]
        if(length(plasmid)>2)
        {
          plasmid <- plasmid[2:length(plasmid)]
          plasmid <- paste(plasmid,collapse="_")
        }
        else
        {
          plasmid <- plasmid[2]
        }
        if (tolower(plasmid) %in% references)
        {
          cids <- c(cids, CG3$V2[cRow])
        }
      }
      cids2 <- unique(cids)
      CG3$V8 <- FALSE
      CG3[which(CG3$V2 %in% cids),"V8"] <- TRUE
      
      if (input$Subset) # For subsetting only: modify CG3
      {
        CG3<-subsetData(CG3,input$ranges)
      }
      
      #CG3$V8<-c(1:nrow(CG3)) # CRITICAL: Used for labeling and position!
      
      seqData <<- as.character(CG3$V3) # CRITICAL: Used to write file
      
      print(ggplot(CG3,aes(x=c(1:nrow(CG3)),y=CG3$V4))+ # Notice X position is determined by V7
               geom_bar(stat="identity",
                        fill=CG3$V6,
                        alpha=CG3$V5
               )+
               scale_y_continuous()+
               scale_x_discrete()+ # Remove extra space
               geom_text(size=3.75,y=max(CG3$V4)/2.0,
                         fontface=ifelse(CG3$V8,"bold","plain"),
                         label=CG3$V3,
                         color="black") +
               geom_label(size=3.5,y=-5,
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

      print("READING FASTA... (MAY TAKE UP TO 5 MINUTES)")
      plasmidAA <- read.fasta("Plasmids20-200kb-6-9-2016AA.fa","AA")
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
      print(input$variable)
      fourName <- input$variable
      sNames <- c(fourName) # AS of 7/15/16 NEEDS TO BE FIXED
      # THE FIRST NAME DOES NOT HAVE THE OLD NAME APPENDED
      for (i in c(1:length(Seqs)))
      {
        seqName <- attr(Seqs[i],"name")
        if (grepl("<unknown description>",seqName))
        {
          seqName <- strsplit(seqName," ")[[1]][1]
        }
        sNames <- c(sNames,paste(fourName,seqName,sep="_"))
      }
      print(sNames)
      write.fasta(Seqs,sNames,
                  file.out=paste(fourName, "faa", sep = "."))
    }
  )
})