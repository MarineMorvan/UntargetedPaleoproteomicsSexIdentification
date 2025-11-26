library(shiny)
library(dplyr)
library(ggplot2)
library(DT)

options(shiny.maxRequestSize = 600 * 1024^2)

ui <- fluidPage(
  titlePanel("Untargeted paleoproteomics sex identification"),
  
  sidebarPanel(
    wellPanel(
      tags$h4("Upload your data:"),
      fileInput("peptides", label = h5("peptides.txt"), accept = ".txt"),
      fileInput("msms", label = h5("msms.txt"), accept = ".txt")
    ),
    wellPanel(
      tags$h4("Graph"),
      downloadButton('downloadGraph', "Download the Graph")
    ),
    wellPanel(
      tags$h4("Table"),
      downloadButton('downloadTable', "Download the Table")
    )
  ),
  
  mainPanel(
    tabsetPanel(
      type = "tabs",
      tabPanel("Graph and Table",
               plotOutput("plot1"),
               DTOutput("table1")
      ),
      tabPanel("Peptide list",
               selectInput("raw_file_select", "Select Raw file:", choices = NULL),
               DTOutput("peptide_table")
      )
    )
  )
)

server <- function(input, output, session) {
  
  rv <- reactiveValues(
    peptide_list = NULL,
    peptides_raw = NULL,
    checkbox_states = list()
  )
  
  rv_table <- reactiveVal(NULL)
  
  # ------------------- Load peptides.txt -------------------
  observeEvent(input$peptides, {
    req(input$peptides)
    data <- read.delim(input$peptides$datapath, sep="\t", header=TRUE)
    peptides_filter <- data[, c("Proteins", "Sequence", "Start.position", "End.position")]
    
    AMELY <- peptides_filter %>%
      filter(grepl("Q99218", Proteins) & !grepl("Q99217", Proteins)) %>%
      mutate(SMI=59) %>%
      filter(between(SMI, Start.position, End.position))
    
    AMELX <- peptides_filter %>%
      filter(grepl("Q99217", Proteins)) %>%
      mutate(SIR=44, IR=45) %>%
      filter(between(IR, Start.position, End.position)) %>%
      filter(between(SIR, Start.position, End.position))
    
    rv$peptides_raw <- list(AMELX=AMELX, AMELY=AMELY)
  })
  
  # ------------------- Load msms.txt and merge -------------------
  observeEvent(input$msms, {
    req(input$msms, rv$peptides_raw)
    
    msms <- read.delim(input$msms$datapath, sep="\t", header=TRUE)
    
    clean_matches <- function(x){
      if(is.na(x) || x=="") return(NA)
      ions <- unlist(strsplit(x,";"))
      ions <- ions[grepl("^[by]\\d+$", ions)]
      paste(ions, collapse=";")
    }
    msms$Matches <- sapply(msms$Matches, clean_matches)
    
    msms_filter <- msms[, c("Sequence","Precursor.Intensity","Raw.file","Proteins")]
    matches_col <- msms[, c("Sequence","Proteins","Matches")]
    
    AMELX <- merge(rv$peptides_raw$AMELX, msms_filter, by=c("Sequence","Proteins"))
    AMELY <- merge(rv$peptides_raw$AMELY, msms_filter, by=c("Sequence","Proteins"))
    
    AMELX <- merge(AMELX, matches_col, by=c("Sequence","Proteins"))
    AMELY <- merge(AMELY, matches_col, by=c("Sequence","Proteins"), all.x=TRUE)
    
    if(nrow(AMELY)==0){
      AMELY[1,] <- list("A","in case of all samples are female",0,0,0,0,NA)
      names(AMELY) <- names(AMELX)
    }
    
    AMELX <- aggregate(Precursor.Intensity ~ Sequence+Raw.file+Proteins+Start.position+End.position+Matches, AMELX, sum)
    AMELY <- aggregate(Precursor.Intensity ~ Sequence+Raw.file+Proteins+Start.position+End.position+Matches, AMELY, sum)
    
    peptide_list <- rbind(AMELX, AMELY)
    
    count_ions <- function(m, t){
      if(is.na(m) || m=="") return(0)
      sum(grepl(paste0("^",t,"\\d+$"), unlist(strsplit(m,";"))))
    }
    
    peptide_list <- peptide_list %>%
      rowwise() %>%
      mutate(
        `b ions`=paste0(count_ions(Matches,"b"), "/", nchar(Sequence)),
        `y ions`=paste0(count_ions(Matches,"y"), "/", nchar(Sequence)),
        Use = TRUE
      ) %>% ungroup()
    
    rv$peptide_list <- peptide_list
    
    # Initialize checkbox states per Raw.file
    rv$checkbox_states <- peptide_list %>%
      split(.$Raw.file) %>%
      lapply(function(df) setNames(df$Use, seq_len(nrow(df))))
    
    updateSelectInput(session, "raw_file_select", choices=unique(peptide_list$Raw.file))
  })
  
  # ------------------- Peptide table -------------------
  peptide_proxy <- dataTableProxy("peptide_table")
  
  output$peptide_table <- renderDT(server = TRUE, {
    req(rv$peptide_list, input$raw_file_select)
    selected_file <- input$raw_file_select
    df <- rv$peptide_list %>% filter(Raw.file == selected_file)
    
    if(is.null(rv$checkbox_states[[selected_file]])){
      rv$checkbox_states[[selected_file]] <- setNames(df$Use, seq_len(nrow(df)))
    }
    
    df_display <- df
    df_display$`Raw file` <- df_display$Raw.file
    df_display$`Start position` <- df_display$Start.position
    df_display$`End position` <- df_display$End.position
    df_display$Add <- sapply(seq_len(nrow(df_display)), function(i){
      id <- paste0("chk_", selected_file, "_", i)
      checked <- ifelse(isTRUE(rv$checkbox_states[[selected_file]][i]), "checked", "")
      sprintf('<input type="checkbox" id="%s" %s onclick="Shiny.setInputValue(\'%s\', this.checked, {priority: \'event\'})">', id, checked, id)
    })
    
    datatable(df_display %>% select(Add, `Raw file`, Sequence, Proteins, `Start position`, `End position`, `b ions`, `y ions`),
              escape = FALSE,
              rownames = FALSE,
              options = list(scrollX=TRUE, pageLength=100, lengthMenu=c(10,25,50,100))
    )
  })
  
  # ------------------- Update checkbox states -------------------
  observe({
    req(rv$peptide_list, input$raw_file_select)
    selected_file <- input$raw_file_select
    df <- rv$peptide_list %>% filter(Raw.file == selected_file)
    n <- nrow(df)
    
    for(i in seq_len(n)){
      id <- paste0("chk_", selected_file, "_", i)
      if(!is.null(input[[id]])){
        rv$checkbox_states[[selected_file]][i] <- input[[id]]
      }
    }
    
    df_display <- df
    df_display$`Raw file` <- df_display$Raw.file
    df_display$`Start position` <- df_display$Start.position
    df_display$`End position` <- df_display$End.position
    df_display$Add <- sapply(seq_len(nrow(df_display)), function(i){
      id <- paste0("chk_", selected_file, "_", i)
      checked <- ifelse(isTRUE(rv$checkbox_states[[selected_file]][i]), "checked", "")
      sprintf('<input type="checkbox" id="%s" %s onclick="Shiny.setInputValue(\'%s\', this.checked, {priority: \'event\'})">', id, checked, id)
    })
    
    replaceData(peptide_proxy, df_display %>% select(Add, `Raw file`, Sequence, Proteins, `Start position`, `End position`, `b ions`, `y ions`),
                resetPaging = FALSE, rownames = FALSE, clearSelection = FALSE)
  })
  
  # ------------------- Recompute table1 / plot1 -------------------
  observe({
    req(rv$peptide_list, rv$checkbox_states)
    all_raw <- unique(rv$peptide_list$Raw.file)
    
    df_list <- lapply(all_raw, function(raw_file){
      df <- rv$peptide_list %>% filter(Raw.file == raw_file)
      chk <- rv$checkbox_states[[raw_file]]
      df$Use <- chk
      df[df$Use==TRUE, ]
    })
    
    df <- do.call(rbind, df_list)
    
    if(nrow(df)==0){
      df <- data.frame(
        `Raw file`=character(0),
        `Protein X`=character(0),
        `Protein Y`=character(0),
        Precursor.Intensity.AMELX=numeric(0),
        Precursor.Intensity.AMELY=numeric(0),
        logAMELX=numeric(0),
        logAMELY=numeric(0),
        `Biological sex`=character(0),
        Label=character(0)
      )
    } else {
      AMELX <- df[grepl("Q99217", df$Proteins), ]
      AMELY <- df[grepl("Q99218", df$Proteins), ]
      
      AMELX_int <- if(nrow(AMELX)>0) aggregate(Precursor.Intensity ~ Raw.file + Proteins, AMELX, sum) else data.frame(Raw.file=character(0), Proteins=character(0), Precursor.Intensity=numeric(0))
      AMELY_int <- if(nrow(AMELY)>0) aggregate(Precursor.Intensity ~ Raw.file + Proteins, AMELY, sum) else data.frame(Raw.file=character(0), Proteins=character(0), Precursor.Intensity=numeric(0))
      
      names(AMELX_int)[3] <- "Precursor.Intensity.AMELX"
      names(AMELY_int)[3] <- "Precursor.Intensity.AMELY"
      
      df <- merge(AMELX_int, AMELY_int, by="Raw.file", all=TRUE)
      df <- df[!grepl("in case of all samples are female", df$Proteins.y), ]
      
      df$`Raw file` <- df$Raw.file
      df$`Protein X` <- df$Proteins.x
      df$`Protein Y` <- df$Proteins.y
      df$logAMELX <- log(df$Precursor.Intensity.AMELX)
      df$logAMELY <- log(df$Precursor.Intensity.AMELY)
      df$logAMELY[!is.finite(df$logAMELY)] <- 0
      df$`Biological sex` <- ifelse(df$logAMELY>1, "Male", "Female")
      
      if(!"Label" %in% names(df)) df$Label <- ""
    }
    
    rv_table(df)
  })
  
  # ------------------- Table1 output -------------------
  output$table1 <- renderDT({
    df <- rv_table()
    req(df)
    df <- df[, c("Raw file", "Label", "Protein X", "Protein Y", "logAMELX", "logAMELY", "Biological sex")]
    datatable(df,
              rownames=FALSE,
              editable=list(target="cell", disable=list(columns=c(0,2,3,4,5,6))),
              options=list(scrollX=TRUE))
  })
  
  # ------------------- Update Label -------------------
  observeEvent(input$table1_cell_edit, {
    info <- input$table1_cell_edit
    df <- rv_table()
    req(df)
    
    row <- info$row
    col <- info$col
    value <- info$value
    
    if (col == 1) df[row, "Label"] <- value
    
    rv_table(df)
  })
  
  # ------------------- Plot output -------------------
  output$plot1 <- renderPlot({
    df <- rv_table()
    req(df)
    ggplot(df, aes(logAMELX, logAMELY, color=`Biological sex`)) +
      geom_point(size=4) +
      geom_text(aes(label=Label), vjust=-1, size=4, show.legend = FALSE) +
      scale_color_manual(values=c("#ffb627","#60d394")) +
      theme_minimal(base_size=16)
  })
  
  # ------------------- Downloads -------------------
  output$downloadGraph <- downloadHandler(
    filename=function() paste0("graph-DDA-SexID-", Sys.Date(), ".tiff"),
    content=function(file) ggsave(file, plot=last_plot(), width=12, height=8, device="tiff")
  )
  
  output$downloadTable <- downloadHandler(
    filename=function() paste0("table-DDA-SexID-", Sys.Date(), ".csv"),
    content=function(file) write.csv(rv_table(), file, row.names=FALSE)
  )
  
}

shinyApp(ui, server)
