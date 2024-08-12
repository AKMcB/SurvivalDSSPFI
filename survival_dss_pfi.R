library(shiny)
library(bslib)
library(mclust)
library(tidyverse)
library(survival)
library(survminer)
set.seed(123)

expr <- readRDS("survival_tcga.rds")
gene_data <- read.csv("expr_cutoff_genes_cancer_type_tcga.csv")


create_km <- function(data, gene_data, GENE, CANCER, analysis_type) {
  expr <- expr %>% select(c("id", GENE, `cancer type abbreviation`, all_of(analysis_type), paste(analysis_type,".time", sep = "")))
  colnames(expr)[2] <- "gene"
  expr <- subset(expr, expr$`cancer type abbreviation` == CANCER)
  info_surv <- select(expr, c("id", "gene", "cancer type abbreviation", analysis_type, paste(analysis_type,".time", sep = "")))
  
  info_surv$cancer_type <- as.factor(info_surv$`cancer type abbreviation`)
  
  # Extract cutoff based on the gene and cancer type of interest
  cutoff <- gene_data %>% filter(gene== GENE & cancer == CANCER) %>% pull(cutoff_value)
  
  if (length(cutoff) == 0) {
    message("no cutoff value found for the combination")
    return(NULL)
  } 
  
  # Define gene expression based on cutoff
  info_surv$gene_expression <- ifelse(info_surv$gene >= cutoff, 'High', "Low")
  
  if (all(is.na(info_surv$gene_expression))){ 
    message("all values are NA")
    return(NULL)
  }
  
  column <- paste(analysis_type, ".time", sep = "")
  info_surv <- subset(info_surv, info_surv[[column]] > 0)
  
  if (nrow(info_surv) == 0) {
    message("no valdig data is available")
    return(NULL)
  } 
  
  # Convert time into years
  info_surv$years <- info_surv[[column]] / 365
  info_surv <- info_surv[complete.cases(info_surv),]
  
  # Check if the gene expression groups have sufficient observations
  if (length(unique(info_surv$gene_expression)) < 2) {
    message("not enough data points in gene expression groups")
    return(NULL)
  } 
  
  column_analysis <- paste(analysis_type)
  # Define survival
  survival = Surv(time = info_surv$years, event = info_surv[[column_analysis]])
  
  survival_fit <- survfit(formula = survival ~ gene_expression, data = info_surv)
  
  p <- ggsurvplot(
    fit = survival_fit,
    pval = TRUE,
    surv.median.line = "hv", legend = c(0.1, 0.1),
    xlab = paste(analysis_type, "(Years)"),
    ylab = paste(analysis_type, "Probability"),
    title = paste(GENE, analysis_type, "for", CANCER),
    palette = c("#E69F00", "#0072B2"),
    pval.coord = c(0.1,0.2), 
    break.x.by = 5,         
    conf.int = TRUE,
    risk.table = TRUE, risk.table.title = "",
    risk.table.height = 0.15,
    ncensor.plot = TRUE,
    ncensor.plot.height = 0.15,
    legend.labs = c("High", "Low"),
    legend.title = paste(GENE), 
    tables.theme = theme(
      axis.text.x = element_blank(), 
      axis.ticks.x = element_blank(), 
      axis.title.x = element_blank(), 
      axis.title.y = element_blank(),
      axis.text.y = element_text(size = 14)
    ), 
    font.x = 14, font.y = 14, font.tickslab = 14
  )
  return(p)
}

# Shiny UI
ui <- fluidPage(
  titlePanel("Survival Plot: DSS or PFI"),
  
  sidebarLayout(
    sidebarPanel(
      selectInput("GENE", "Select Gene:", choices = sort(unique(gene_data$gene))),
      selectInput("CANCER", "Select Cancer Type:", choices = sort(unique(gene_data$cancer))),
      selectInput("analysis_type", "Survival Type", choices = c("DSS", "PFI")),
      sliderInput("y_axis", "Y Axis Limit", min = 0, max = 1, value = c(0, 1)),
      sliderInput("x_axis", "X Axis Limit (Years)", min = 0, max = 25, value = c(0, 25)),
      actionButton("update", "Update"),
      downloadButton("downloadPlot", "Download Plot")
    ),
    
    mainPanel(
      plotOutput("SurvPlot")
    )
  )
)


# Shiny Server
server <- function(input, output, session) {
  
  observeEvent(input$update, {
    result <- create_km(expr, gene_data, input$GENE, input$CANCER, input$analysis_type)
    
    output$SurvPlot <- renderPlot({
      if (!is.null(result)) {
        # Extract the main ggplot object from ggsurvplot result
        main_plot <- result$plot + 
          xlim(input$x_axis[1], input$x_axis[2]) + 
          ylim(input$y_axis[1], input$y_axis[2])
        
        # Print the main plot
        print(main_plot)
      }
    })
  })
  
  output$downloadPlot <- downloadHandler(
    filename = function() {
      paste("survival_", input$CANCER, "_", input$GENE, "_", input$analysis_type, ".pdf", sep = "")
    },
    content = function(file) {
      p <- create_km(expr, gene_data, input$GENE, input$CANCER, input$analysis_type)
      if (!is.null(p)) {
        main_plot <- p$plot + 
          xlim(input$x_axis[1], input$x_axis[2]) + 
          ylim(input$y_axis[1], input$y_axis[2])
        ggsave(file, plot = main_plot, device = "pdf")
      }
    }
  )
}

shinyApp(ui = ui, server = server)
