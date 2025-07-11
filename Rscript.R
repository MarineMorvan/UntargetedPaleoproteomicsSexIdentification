install.packages("shiny")
install.packages("bslib")
install.packages("stringr")
install.packages("tidyverse")
install.packages("ggplot2")
install.packages("hrbrthemes")
install.packages("DT")
library(shiny)
library(bslib)
library(stringr)
library(tidyverse)
library(ggplot2)
library(hrbrthemes)
library(DT)

ui <- fluidPage(
  
  titlePanel("Untargeted paleoproteomics sex identification"),
  
  sidebarPanel(
    wellPanel(
      tags$h4("Upload your data:"),
      fileInput("peptides", label = h5("peptides.txt"),accept = c(".txt")),
      fileInput("msms", label = h5("msms.txt"),accept = c(".txt"))
    ),
    wellPanel(
      tags$h4("Graph"),
      actionButton("GraphButton", "Plot graph"),
      downloadButton('downloadGraph',"Download the Graph")
    ),
    wellPanel(
      tags$h4("Table"),
      actionButton("TableButton", "Plot table"),
      downloadButton('downloadTable',"Download the Table")
    ),
),
 
  mainPanel(
    tabsetPanel(type = "tabs",
                tabPanel("Graph and Table",
                         plotOutput("plot1"),  
                         DT::dataTableOutput("table1")
                ) 
    )
  )
)

options(shiny.maxRequestSize = 300*1024^2)
server <- function(input, output, session) {

observeEvent(input$peptides, {
    req(input$peptides)  # Ensure the file is uploaded
    data <- read.delim(input$peptides$datapath, header = TRUE, sep = "\t")   
    assign("peptides_data", data, envir = .GlobalEnv)    
    print("Peptides file has been saved to the global environment as 'peptides_data'")
    peptides_filter <- data[, c("Proteins", "Sequence", "Start.position", "End.position")]
    assign("peptides_filter", peptides_filter, envir = .GlobalEnv)
    AMELY <- subset(peptides_filter,Proteins=='sp|Q99218|AMELY_HUMAN' )
    AMELX <- subset(peptides_filter,Proteins=='sp|Q99217|AMELX_HUMAN' )
    AMELY <- AMELY %>% add_column(SMI = 59)
    AMELY <- AMELY%>% filter( between(AMELY$SMI, AMELY$Start.position, AMELY$End.position) )
    AMELX <- AMELX %>% add_column(SIR = 44)
    AMELX <- AMELX %>% add_column(IR = 45)
    AMELX <- AMELX%>% filter( between(AMELX$IR, AMELX$Start.position, AMELX$End.position) )
    AMELX <- AMELX%>% filter( between(AMELX$SIR, AMELX$Start.position, AMELX$End.position) )
    assign("AMELY", AMELY, envir = .GlobalEnv)
    assign("AMELX", AMELX, envir = .GlobalEnv)
})

observeEvent(input$msms, {
     req(input$msms)  # Ensure the file is uploaded 
     data <- read.delim(input$msms$datapath, header = TRUE, sep = "\t")     
     assign("msms_data", data, envir = .GlobalEnv)      
     print("Msms file has been saved to the global environment as 'msms_data'")
     msms_filter <- data[, c("Sequence", "Precursor.Intensity", "Raw.file", "Proteins")]
     assign("msms_filter", msms_filter, envir = .GlobalEnv)
     AMELY <- merge.data.frame(AMELY, msms_filter, by.x = c("Sequence", "Proteins"), by.y = c("Sequence", "Proteins"))
     AMELX <- merge.data.frame(AMELX, msms_filter, by.x = c("Sequence", "Proteins"), by.y = c("Sequence", "Proteins"))
     AMELY[nrow(AMELY) + 1,] = list(Sequence = "A", Proteins = "in case of all samples are female", Start.position=0, End.position=0, Precursor.Intersity=0, Raw.file=0)
     AMELX <- aggregate(Precursor.Intensity~Sequence+Raw.file+Proteins+Start.position+End.position,data=AMELX,FUN=sum)
     AMELY <- aggregate(Precursor.Intensity~Sequence+Raw.file+Proteins+Start.position+End.position,data=AMELY,FUN=sum)
     AMELX_intensities <- aggregate(Precursor.Intensity~Raw.file+Proteins,data=AMELX,FUN=sum)
     AMELY_intensities <- aggregate(Precursor.Intensity~Raw.file+Proteins,data=AMELY,FUN=sum)
     names(AMELX_intensities)[names(AMELX_intensities) == 'Precursor.Intensity'] <- 'Precursor.Intensity.AMELX'
     names(AMELY_intensities)[names(AMELY_intensities) == 'Precursor.Intensity'] <- 'Precursor.Intensity.AMELY'
     assign("AMELY_intensities", AMELY_intensities, envir = .GlobalEnv)
     assign("AMELX_intensities", AMELX_intensities, envir = .GlobalEnv) 
     Total_intensities <- merge(AMELX_intensities, AMELY_intensities,by = intersect(AMELX_intensities(Raw.file), AMELY_intensities(Raw.file)), by.x = c("Raw.file"), by.y = c("Raw.file"), all = TRUE)
     Total_intensities <- Total_intensities[!grepl("in case of all samples are female", Total_intensities$Proteins.y),]
     Total_intensities <- Total_intensities %>%
       add_column(logAMELX = log(Total_intensities$Precursor.Intensity.AMELX))
     Total_intensities <- Total_intensities %>%
       add_column(logAMELY = log(Total_intensities$Precursor.Intensity.AMELY))
     Total_intensities$logAMELY <- replace(Total_intensities$logAMELY, is.na(Total_intensities$logAMELY), 0)
     Total_intensities <- Total_intensities[!grepl(-Inf, Total_intensities$logAMELY),]
     Total_intensities <- Total_intensities %>%
       add_column(Biological_sex = if_else(Total_intensities$logAMELY >1, "Male", "Female"))
     Total_intensities = subset(Total_intensities, select = -c(Precursor.Intensity.AMELX, Precursor.Intensity.AMELY))
     assign("Total_intensities", Total_intensities, envir = .GlobalEnv)
}) 

observeEvent(input$GraphButton, {
  req(Total_intensities)  # Ensure the data is available
  output$plot1 <- renderPlot({
    Total_intensities %>%
      ggplot(aes(x = logAMELX, y = logAMELY, colour = Biological_sex, size = 1)) +
      geom_point() +
      scale_color_manual(values = c("#ffb627", "#60d394")) +
      theme(axis.title.x = element_text(size = 20), axis.title.y = element_text(size = 20),legend.text = element_text(size=18), legend.title = element_text(size=18),axis.text = element_text(size = 18)) +
      labs(color = "Legend", x = bquote(Log[10]("AMELX")), y = bquote(Log[10]("AMELY"))) +
      guides(size = "none")  
  })
})

output$downloadGraph <- downloadHandler(
  filename = function() {
    paste("graph-DDA-SexID-", Sys.Date(), ".tiff", sep = "")
  },
  content = function(file) {
    ggsave(file, plot = last_plot(), device = "tiff", width = 12, height = 8)
  }
)
  
observeEvent(input$TableButton, {
  req(Total_intensities)  # Ensure the data is available
  output$table1 <- DT::renderDataTable({
    DT::datatable(Total_intensities, 
                  options = list(scrollX = TRUE, buttons = c('copy', 'csv', 'excel')),
                  rownames = TRUE
    )
  })
})

output$downloadTable <- downloadHandler(
  filename = function() {
    paste("table-DDA-SexID-", Sys.Date(), ".csv", sep = "")
  },
  content = function(file) {
    write.csv(Total_intensities, file, row.names = TRUE)
  }
)
}
shinyApp(ui = ui, server = server)
