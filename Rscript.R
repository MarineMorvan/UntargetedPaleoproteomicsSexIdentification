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

options(shiny.maxRequestSize = 600*1024^2)
server <- function(input, output, session) {
  
  # Reactive storage for the processed data table
  reactive_data <- reactiveVal()
  
  observeEvent(input$peptides, {
    req(input$peptides)
    data <- read.delim(input$peptides$datapath, header = TRUE, sep = "\t")   
    peptides_filter <- data[, c("Proteins", "Sequence", "Start.position", "End.position")]
    AMELY <- subset(peptides_filter, grepl("Q99218", Proteins) & !grepl("Q99217", Proteins))
    AMELX <- subset(peptides_filter, grepl("Q99217", Proteins))
    AMELY <- AMELY %>% add_column(SMI = 59)
    AMELY <- AMELY %>% filter(between(SMI, Start.position, End.position))
    AMELX <- AMELX %>% add_column(SIR = 44)
    AMELX <- AMELX %>% add_column(IR = 45)
    AMELX <- AMELX %>% filter(between(IR, Start.position, End.position))
    AMELX <- AMELX %>% filter(between(SIR, Start.position, End.position))
    assign("AMELY", AMELY, envir = .GlobalEnv)
    assign("AMELX", AMELX, envir = .GlobalEnv)
  })
  
  observeEvent(input$msms, {
    req(input$msms)
    data <- read.delim(input$msms$datapath, header = TRUE, sep = "\t")
    msms_filter <- data[, c("Sequence", "Precursor.Intensity", "Raw.file", "Proteins")]
    
    AMELY <- get("AMELY", envir = .GlobalEnv)
    AMELX <- get("AMELX", envir = .GlobalEnv)
    
    AMELY <- merge(AMELY, msms_filter, by = c("Sequence", "Proteins"))
    AMELX <- merge(AMELX, msms_filter, by = c("Sequence", "Proteins"))
    
    AMELY[nrow(AMELY) + 1,] <- list(Sequence = "A", Proteins = "in case of all samples are female", Start.position=0, End.position=0, Precursor.Intensity=0, Raw.file=0)
    
    AMELX <- aggregate(Precursor.Intensity ~ Sequence + Raw.file + Proteins + Start.position + End.position, data = AMELX, FUN = sum)
    AMELY <- aggregate(Precursor.Intensity ~ Sequence + Raw.file + Proteins + Start.position + End.position, data = AMELY, FUN = sum)
    
    AMELX_intensities <- aggregate(Precursor.Intensity ~ Raw.file + Proteins, data = AMELX, FUN = sum)
    AMELY_intensities <- aggregate(Precursor.Intensity ~ Raw.file + Proteins, data = AMELY, FUN = sum)
    
    names(AMELX_intensities)[names(AMELX_intensities) == 'Precursor.Intensity'] <- 'Precursor.Intensity.AMELX'
    names(AMELY_intensities)[names(AMELY_intensities) == 'Precursor.Intensity'] <- 'Precursor.Intensity.AMELY'
    
    # Merge intensities on Raw.file
    Total_intensities <- merge(
      AMELX_intensities, AMELY_intensities,
      by = "Raw.file",
      all = TRUE,
      suffixes = c(".x", ".y")
    )
    
    Total_intensities <- Total_intensities[!grepl("in case of all samples are female", Total_intensities$Proteins.y), ]
    
    # Add log columns, handle zeros and NAs
    Total_intensities$logAMELX <- log(Total_intensities$Precursor.Intensity.AMELX)
    Total_intensities$logAMELY <- log(Total_intensities$Precursor.Intensity.AMELY)
    Total_intensities$logAMELY[is.na(Total_intensities$logAMELY) | is.infinite(Total_intensities$logAMELY)] <- 0
    
    # Biological sex determination
    Total_intensities$Biological_sex <- ifelse(Total_intensities$logAMELY > 1, "Male", "Female")
    
    # Add empty Label column
    Total_intensities$Label <- ""

    # Add empty Group column
    Total_intensities$Group <- ""
    
    # Reorder/select columns
    Total_intensities <- Total_intensities[, c("Raw.file", "Label", "Group", "Proteins.x", "Proteins.y", "logAMELX", "logAMELY", "Biological_sex")]
    
    # Reduce log values with 5 decimals
    Total_intensities$logAMELX<-round(Total_intensities$logAMELX,5)
    Total_intensities$logAMELY<-round(Total_intensities$logAMELY,5)

    # Save reactive data
    reactive_data(Total_intensities)
  })
  
  # Render plot reactively
  output$plot1 <- renderPlot({
    data <- reactive_data()
    req(data)
    
    ggplot(data, aes(x = logAMELX, y = logAMELY, colour = Biological_sex, shape = Group)) +
      geom_point(size=3.5) +
      geom_text(aes(label = Label), size = 3.5, check_overlap = FALSE, vjust = -1, show.legend = F) +
      scale_color_manual(values = c("#ffb627", "#60d394")) +
      theme(
        axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size = 20),
        legend.text = element_text(size = 18),
        legend.title = element_text(size = 18),
        axis.text = element_text(size = 18)
      ) +
      labs(
        color = "Legend",
        x = bquote(Log[10]("AMELX")),
        y = bquote(Log[10]("AMELY"))
      ) +
      guides(size = "none")
  })
  
  # Render editable datatable reactively
  output$table1 <- DT::renderDataTable({
    data <- reactive_data()
    req(data)
    
    DT::datatable(
      data,
      editable = list(target = 'cell', disable = list(columns = c(0, 3, 4, 5, 6, 7))),  # Only 'Label' and 'Group' editable (col 1 and 2)
      options = list(scrollX = TRUE, buttons = c('copy', 'csv', 'excel')),
      rownames = FALSE
    )
  })
  
  # Update reactive data on label and group edit
  observeEvent(input$table1_cell_edit, {
  info <- input$table1_cell_edit
  data <- reactive_data()
  req(data)
  
  row <- info$row
  col <- info$col
  value <- info$value
  
  # Column indices are 0-based
  if (col == 1) {  # Label
    data[row, "Label"] <- value
  } else if (col == 2) {  # Group
    data[row, "Group"] <- value
  }
  
  reactive_data(data)
})



  # Download plot as TIFF
  output$downloadGraph <- downloadHandler(
    filename = function() {
      paste("graph-DDA-SexID-", Sys.Date(), ".tiff", sep = "")
    },
    content = function(file) {
      data <- reactive_data()
      req(data)
      ggsave(file, plot = last_plot(), device = "tiff", width = 12, height = 8)
    }
  )
  
  # Download table as CSV
  output$downloadTable <- downloadHandler(
    filename = function() {
      paste("table-DDA-SexID-", Sys.Date(), ".csv", sep = "")
    },
    content = function(file) {
      data <- reactive_data()
      req(data)
      write.csv(data, file, row.names = FALSE)
    }
  )
}
shinyApp(ui = ui, server = server)
