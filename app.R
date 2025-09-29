library(shiny)
library(dplyr)
library(ggplot2)
library(stringr)
library(ggbeeswarm)
library(ggpubr)

load("Omicron_DIA.RData")

ui <- fluidPage(
  titlePanel("BAL DIA-MS Omicron Proteomics Analysis"),
  
  sidebarLayout(
    sidebarPanel(
      radioButtons(
        inputId = "biomakerID",
        label = "Please select the biomarkers based on",
        choices = c("Protein Accessions (UniProt ID)",
          "Protein Descriptions"
          
        )
      ),
      
      
      conditionalPanel(
        "input.biomakerID == 'Protein Accessions (UniProt ID)'",
        selectizeInput(
          inputId = "ProteinAccessions",
          label = "Protein Accessions (UniProt ID; Choose from the drop-down list or Type to search)",
          choices = sort(ProteinIndex$PG.ProteinAccessions),
          selected = ProteinIndex$PG.ProteinAccessions[1410],
          options = list(
            maxOptions = 2850,
            placeholder = "Type to search..."
          )
        )
      ),
      
      conditionalPanel(
        "input.biomakerID == 'Protein Descriptions'",
        selectizeInput(
          inputId = "ProteinDescriptions",
          label = "Protein Descriptions (Choose from the drop-down list or Type to search)",
          choices = sort(ProteinIndex$PG.ProteinDescriptions),
          selected = ProteinIndex$PG.ProteinDescriptions[1410],
          options = list(
            maxOptions = 2850,
            # substring search is enabled by default
            # but you can add more options here if needed
            placeholder = "Type to search..."
          )
        )
      ),
      
      
      hr(),
      br(),
      wellPanel(
        h4("Biomarker to be assessed ", style = "color:brown;"),
        tableOutput("Alltable")
      )
    ),
    
    mainPanel(
      h4(
        "1. Beeswarm plots comparing Acute (A), Recovery (R), and Convalescent(C) phases with Control (CNTL) Groups:",
        style = "color:blue;"
      ),
      plotOutput("beePlot"),
      br(),
      h4(
        "2. Beeswarm plots comparing Acute (A), Recovery (R), and Convalescent (C) phases:",
        style = "color:blue;"
      ),
      plotOutput("beePlot2")
    )
  )
)

server <- function(input, output) {

  # index table to use later
  IndexTable<- reactive({
    ID= input$biomakerID
    AccessionID =  as.character(input$ProteinAccessions)
    DescID        =  as.character(input$ProteinDescriptions)
 
    Index2 <- NULL
    
    if (ID=='Protein Accessions (UniProt ID)') {Index2= ProteinIndex  %>%filter(PG.ProteinAccessions==AccessionID) }
    else {if (ID=='Protein Descriptions') {Index2= ProteinIndex  %>%filter(PG.ProteinDescriptions==DescID) }
    }
    
    if (nrow(Index2) == 0) {
      # Return a placeholder empty data frame with correct colnames if no match
      return(data.frame(
         `Protein Descriptions` = character(0),
         `Protein Accessions` = character(0)
      ))
    }
    
    Index2$No <- NULL
    return(Index2[, c(1,2)])
    
  })
  
  
  #-----
  output$Alltable <- renderTable({
    df <- IndexTable()
    if (nrow(df) == 0) return(NULL)  # Show nothing when no match
    data.frame(Biomarker = t(df))
  })
 
  
  #-----
  
  output$beePlot <- renderPlot({

    df <- IndexTable()
    req(nrow(df) > 0)  # Stop if no data
    
    ProtID= as.character(df[,1])
    Prot_title=  as.character(df[,2])
    

    Pt_sub= allpts %>% filter(PG.ProteinAccessions==ProtID )

    Pt_sub$R.Condition =factor(Pt_sub$R.Condition, levels = c("A", "R", "C" ,"CNTL"))

    Pt_sub2A = data.frame( part="MS1",   Pt_sub[, c("ptID","R.Condition","PG.MS1Quantity")])
    Pt_sub2A<- rename(Pt_sub2A, MS_value="PG.MS1Quantity" )

    Pt_sub2B =  data.frame( part="MS2",   Pt_sub[, c("ptID","R.Condition","PG.MS2Quantity")])
    Pt_sub2B<- rename(Pt_sub2B, MS_value="PG.MS2Quantity" )

    Pt_sub2= rbind(Pt_sub2A, Pt_sub2B )

    ggplot(data= Pt_sub2, aes(x=R.Condition, y=MS_value,  col=R.Condition)) +
               geom_quasirandom(alpha=0.7, size=1.2, width=0.3, aes(group=ptID))+
               scale_y_log10() +
               labs( x="", y="", title=Prot_title) +
               scale_color_manual(values=rev(c("forestgreen", "dodgerblue" ,"orange","red"))) +   guides(color="none")+
               # add median
               stat_summary(  fun = "median", shape=15,  size = 0.7,alpha=1 ) +
               stat_summary( fun = "median", geom = "crossbar", width = 0.6, size=0.6, alpha=0.9)+   theme_bw()+
               facet_wrap(vars(part) , scales='free', ncol=2)+
               theme(axis.text=element_text( size=13)  ,
                     plot.title =    element_text( size=15, face='bold', color = 'blue'),
                     text= element_text( size=11),
                     strip.text = element_text(size = 12,face='bold', color = 'blue' ),
                     legend.text = element_text( size=12))

  },
  # height = 420, width = 600
  )

  
  output$beePlot2 <- renderPlot({
    
    df <- IndexTable()
    req(nrow(df) > 0)  # Stop if no data
    
    
    ProtID= as.character(df[,1])
    
    Prot_title=  as.character(df[,2])
      
    Pt_sub= allpts %>% filter(PG.ProteinAccessions==ProtID & R.Condition != "CNTL" )
    
    Pt_sub$R.Condition =factor(Pt_sub$R.Condition, levels = c("A", "R", "C"))
    
    Pt_sub2A = data.frame( part="MS1",   Pt_sub[, c("ptID","R.Condition","PG.MS1Quantity")])
    Pt_sub2A<- rename(Pt_sub2A, MS_value="PG.MS1Quantity" )
    
    Pt_sub2B =  data.frame( part="MS2",   Pt_sub[, c("ptID","R.Condition","PG.MS2Quantity")])
    Pt_sub2B<- rename(Pt_sub2B, MS_value="PG.MS2Quantity" )
    
    Pt_sub2= rbind(Pt_sub2A, Pt_sub2B )
    
    ggplot(data= Pt_sub2, aes(x=R.Condition, y=MS_value,  col=R.Condition)) +
               geom_quasirandom(alpha=0.7, size=1.2, width=0.3, aes(group=ptID))+
               scale_y_log10() +
               labs( x="", y="", title=Prot_title) +
               scale_color_manual(values=rev(c("dodgerblue" ,"orange","red"))) +   guides(color="none")+
               # add median
               stat_summary(  fun = "median", shape=15,  size = 0.7,alpha=1 ) +
               stat_summary( fun = "median", geom = "crossbar", width = 0.6, size=0.6, alpha=0.9)+   theme_bw()+
               facet_wrap(vars(part) , scales='free', ncol=2)+
               theme(axis.text=element_text( size=13)  ,
                     plot.title =    element_text( size=15, face='bold', color = 'blue'),
                     text= element_text( size=11),
                     strip.text = element_text(size = 12,face='bold', color = 'blue' ),
                       legend.text = element_text( size=12))
    
  },
  # height = 420, width = 600
  )
  

}




# Create Shiny app ----
shinyApp(ui = ui, server = server)

