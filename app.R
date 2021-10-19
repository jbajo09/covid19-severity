# Before the deploy to shinyapps.io, run:
# library(BiocManager)
# options(repos = BiocManager::repositories())
# See here:
# https://community.rstudio.com/t/failing-to-deploy-shinyapp-depending-on-bioconductor-packages/6970/5


# Packages
library(shiny)
library(shinydashboard)
library(dplyr)
library(KnowSeq)
library(reshape2)
library(caret)
library(ggplot2)
library(ggalluvial)
library(DT)
library(waiter)
library(ROCR)
library(limma)
library(class)
library(M3C)
options(repos = BiocManager::repositories())

setwd(getwd())
load("index.RData")



# File to slightly modify dataPlot function
#source("www/dataPlot.R")

# Define some spinners
spinner_abrir <- tagList(
  spin_folding_cube(),
  span(br(), h4("Loading application..."), style="color:white;")
)

spinner <- tagList(
  spin_chasing_dots(),
  span(br(), h4("Loading..."), style="color:white; display: inline-block;")
)

ui <- dashboardPage(title = "COVID-19 Severity", # Title in web browser
                    ## Theme
                    skin = "black",
                    ## Header
                    dashboardHeader(title = span(
                      "COVID-19",
                      style = "font-family: Lucida Console; font-weight: bold"
                    )),
                    ## Sidebar
                    dashboardSidebar(
                      tags$head(
                        tags$link(rel = "stylesheet", type = "text/css", href = "style.css")
                      ),
                      
                      sidebarMenu(
                        menuItem("Introduction", tabName = "intro", icon = icon("file-alt")),
                        menuItem("Data loading", tabName = "datos", icon = icon("database")),
                        menuItem('Preprocessing', tabName = 'preprocessing', icon = icon('tools')),
                        menuItem("Genes selection", tabName = "genes", icon = icon("dna")),
                        menuItem("training-test", tabName = "entrenamiento", icon = icon("play")),
                        menuItem("LOOCV", tabName = "validation", icon = icon("check-circle")),
                        menuItem("T-SNE", tabName = "tsne", icon = icon("chart-line")),
                        menuItem("Code", tabName = "codigo", icon = icon("code"))
                      )
                    ),
                    ## Body
                    dashboardBody(
                      use_waiter(),
                      # Spinners to show on load, or when the application is busy
                      #waiter_show_on_load(spinner_abrir, color = "#027368"),
                      #waiter_on_busy(spinner, color = "#027368"),
                      tabItems(
                        # Tab 1
                        tabItem(tabName = "intro",
                                
                                h1("About this web application"),
                                tags$i("COVID-19 severity biomarkeRs"), "is a web application that allows users with no previous knowledge of programming to carried out COVID-19 transcriptomic RNA-Seq and developed a COVID-19 severity CDSS ",
                                br(), br(),
                                
                                "The ", tags$i("COVID-19 severity biomarkeRs"), "application is part of the", 
                                tags$a(
                                  " work entitled COVID-19 Severity Index based on Genetic Biomarker using Machine Learning",
                                  href = "https://github.com/jbajo09/covid19-severity",
                                  target="_blank"
                                ),
                                "It's developed in R-Shiny and the code is ",
                                tags$a(
                                  "open source.",
                                  href = "https://github.com/jbajo09/covid19-severity",
                                  target="_blank"
                                ),
                                
                                h2("Abstract "),
                                
                                h3(tags$b("Background")),
                                "A fundamental challenge in the fight against COVID-19 is the realization of reliable and precise tools able
                                to predict the progression of the disease in a patient. This information can be extremely useful to discern hospitalized
                                patients with higher risk of UCI requirement from patients with low severity. How SARS-CoV-2 infection will evolve remains
                                unclear. In this paper, the design of a RNA-Seq genetic biomarker classification system of COVID-19 severity is proposed. Our
                                system can be of great clinical utility for the strategy of planning, organization and administration of human/material resources, 
                                along with automatic classification of the severity of patients affected by COVID-19. A novel pipeline able to integrate RNA-Seq data from different databases, obtaining a genetic biomarker
                                COVID-19 Severity Index using artificial intelligence algorithm, was developed. CD93, RPS24, PSCA and CD300E were identified
                                as a COVID-19 severity gene signature, by carrying out a differentially gene expression analysis and a posterior feature
                                selection step. Furthermore, by using the obtained gene signature, an effective multi-class classifier able to discern between control,
                                outpatient, inpatient and ICU COVID-19 patients was optimized reaching an accuracy of 97.5%.",
                              
                                
                                
                                # Images
                                br(), br(), br(),
                                fluidRow(column(6, tags$img(src = "ugr.png", height = "100px")),
                                         column(6, tags$img(src = "knowseq.png", height = "120px")))
                                
                        ),
                        
                        # Tab 2
                        tabItem(tabName = "datos",
                                
                                # Left column
                                fluidRow(column(6, 
                                                h1("Data loading"),
                                                fileInput(inputId = "file_labels",
                                                          label = span("Select CSV file with labels (see ",
                                                                       tags$a(
                                                                         "here",
                                                                         href = "https://drive.google.com/file/d/1DTdwFhYkwTkfOBs18rtTWdneML1eNnko/view?usp=sharing",
                                                                         target="_blank"
                                                                       ),
                                                                       "an example)"),
                                                          accept = ".csv",
                                                          width = "100%"
                                                ),
                                                fileInput(inputId = "file_DEGsMatrix",
                                                          label = span("Select CSV file expression matrix (see ",
                                                                       tags$a(
                                                                         "here",
                                                                         href = "https://drive.google.com/file/d/1P5zhkstO-ezFQv7mu-yWJfaMfxfaSUXO/view?usp=sharing",
                                                                         target="_blank"
                                                                       ),
                                                                       "an example)"),
                                                          accept = ".csv",
                                                          width = "100%"
                                                ),
                  
                                                
                                                actionButton(inputId = "boton_importar",
                                                             label = "Import file",
                                                             icon = icon("fas fa-file-import", lib = "font-awesome"),
                                                             width = "100%"),
                                                br(),
                                                
                                                conditionalPanel(condition = "input.boton_importar!=0",
                                                                 
                                                                 h2("Distribution of classes"),
                                                                 
                                                                 tableOutput("tabla1")),
                
                                ),
                            
                                )
                                
                        ),
                        
                        
                        # Tab 3
                        tabItem(tabName = "preprocessing",
                                h1("Normalize between arrays"),
                                 'If your data is a collection of different series, you can carry out a normalization between arrays, scaling the log-ratios
to have the same median-absolute-deviation (MAD) across arrays',
                                br(),
                                br(),
                                actionButton(inputId = "boton_normalization",
                                             label = "Normalize between arrays",
                                             icon = icon("dna", lib = "font-awesome"),
                                             width = "50%"),
                                br(),
                                br(),
                                
                                conditionalPanel(condition = "input.boton_normalization!=0",
                                                 
                                                 'Normalization between arrays finished',
                                                 ),
                                
                                br(),
                                br(),
                                
                                h1("Batch effect treatment using SVA"),
                                'In order to treat Batch effect Surrogate Variable Analysis can be carried out',
                                br(),
                                br(),
                                actionButton(inputId = "boton_batch",
                                             label = "Batch effect SVA",
                                             icon = icon("dna", lib = "font-awesome"),
                                             width = "50%"),
                                br(),
                                br(),
                                conditionalPanel(condition = "input.boton_batch!=0",
                                                 
                                                 'Batch effect treatment with SVA completed',
                                ),
                                br(),
                                br(),
                                
                                h1('Outlier removal'),
                                'A majority voting system can be carried out in order to remove outliers. Three different tests are carried out. An interquantile range method, a Kolmogorov Smirnov (K-S) test and MA-plot. ',
                                br(),
                                br(),
                                actionButton(inputId = "boton_outlier",
                                             label = "Remove outlier",
                                             icon = icon("dna", lib = "font-awesome"),
                                             width = "50%"),
                                br(),
                                br(),
                                conditionalPanel(condition = "input.boton_outlier!=0",
                                                 
                                                 'Outliers removed',
                                                 
                                                 textOutput('outlier'),
                                ),
                                
                                br(),
                                
                                conditionalPanel(condition = "input.boton_outlier!=0",
                                                 
                                                 h2("Distribution of classes after preprocessing"),
                                                 
                                                 tableOutput("tabla2"),
                                                 
                                                 textOutput('dim')),

                        ),
                        
                      
                        # Tab 4
                        tabItem(tabName = "genes",
                                'At this point data is splitted in training and test samples with a 80%-20% ratio. In order to replicate the same results the same partition used in our paper is used. Training and test samples are indicated in the repository with the code under training-test file.',
                                h1("DEGs extraction 10-folds CV"),
                                textInput(inputId = "LFC", label = "Select the threshold LFC", value = 1, width = "50%"),
                                
                                sliderInput(inputId = "COV", label = "Select the threshold COV", value = 2, min = 1, max = 3, step = 1, width = "50%"),
                                
                                textInput(inputId = "pvalue", label = "Select the threshold p-value", value = 0.05, width = "50%"),
                                
                                actionButton(inputId = "boton_genes",
                                             label = "Extract DEGs",
                                             icon = icon("dna", lib = "font-awesome"),
                                             width = "50%"),
                                br(),
                                br(),
                                conditionalPanel(condition = "input.boton_genes!=0",
                                                 
                                                 h3("Number of DEGs extracted"),
                                                 br(),
                                                 textOutput("numberDEGs")
                                                 
                                                
                                ),
                                br(),
                                br(),
                                h1('Feature selection algorithm mRMR'),
                                actionButton(inputId = "boton_ranking",
                                             label = "Apply mRMR",
                                             icon = icon("dna", lib = "font-awesome"),
                                             width = "50%"),
                                br(),
                                
                                conditionalPanel(condition = "input.boton_ranking!=0",
                                                 br(),
                                                 'Ranking obtained with mRMR',
                                                 br(),
                                                 #fluidRow(
                                                   #column(4, h4(tags$b("  MRMR")), tableOutput("genes_mrmr")),
                                                 textOutput('genes_mrmr')
                                                 
                                ),
                        ),   
                        
                        # Tab 5
                        tabItem(tabName = "entrenamiento",
                                h1("Model training"),
                                
                                # Choose number of folds
                                selectInput("number_folds",
                                            label = "Number of folds",
                                            choices = c(3, 5, 10),
                                            selected = 5,
                                            width = "50%"),
                                
                                # Train model button
                                actionButton(inputId = "boton_model_training",
                                             label 
                                             = "Train model",
                                             icon = icon("play", lib = "font-awesome"),
                                             width = "50%"),
                                
                                br(),
                                
                                # OPtimal k
                                #textOutput("optimal_knn"),
                                
                                br(),
                                h2("Training CV plot"),
                                plotOutput("train", width = "90%", height = '600px'),
                                br(),
                                br(),
                                h2('Gene expression Boxplot'),
                                br(),
                                actionButton(inputId = "boton_boxplot",
                                             label 
                                             = "BoxPlots",
                                             icon = icon("play", lib = "font-awesome"),
                                             width = "50%"),
                                br(),
                                br(),
                                textInput(inputId = "nboxplot", label = "Select the number of genes from mRMR ranking to get boxplot", value = 8, width = "50%"),
                                br(),
                                plotOutput("boxplot", width = "90%", height = '600px'),
                                br(),
                                h2('Precict over test samples'),
                                br(),
                                sliderInput(inputId = "numero_genes_validation", label = "Select the number of genes to use:",
                                            value = 4, min = 1, max = 50, step = 1, width = "50%"),
                                br(),
                                actionButton(inputId = "test_validation",
                                             label 
                                             = "Predict over test samples",
                                             icon = icon("play", lib = "font-awesome"),
                                             width = "50%"),
                                br(),
                                br(),
                                plotOutput("results_validation", width = "90%", height = '600px'),
                                
                        ),
                        
                        # Tab 5
                        tabItem(tabName = "validation",
                                h1("Model validation"),
                              
                                sliderInput(inputId = "numero_genes_LOOCV", label = "Select the number of genes to use:",
                                            value = 4, min = 1, max = 50, step = 1, width = "50%"),
                                
                                actionButton(inputId = "boton_model_LOOCV",
                                             label = "LOOCV strategy",
                                             icon = icon("play", lib = "font-awesome"),
                                             width = "50%"),
                                
                                br(),

                                plotOutput("results_LOOCV",
                                          width = "60%"),
                                br(),
                                br(),
                                
                                h2('Area Under the Curve'),
                                actionButton(inputId = "boton_calculate_AUC",
                                             label = "Calculate AUC",
                                             icon = icon("play", lib = "font-awesome"),
                                             width = "50%"),
                                textOutput("meanAUC"),
                                plotOutput("AUC",
                                           width = "60%"),
                                
                        ),
                        
                        # Tab 7 
                        tabItem(tabName = 'tsne',
                                h1('T-Distributed Stochastic  Neighbour Embedding'),
                                br(),
                                sliderInput(inputId = "numero_genes_TSNE", label = "Select the number of genes to use:",
                                            value = 4, min = 1, max = 50, step = 1, width = "50%"),
                                br(),
                                actionButton(inputId = "boton_calculate_TSNE",
                                             label = "Carry out T-SNE without any preprocessing step",
                                             icon = icon("play", lib = "font-awesome"),
                                             width = "60%"),
                                br(),
                                conditionalPanel(condition = "input.boton_TSNE!=0",
                                                 plotOutput('tsne', width = '60%')),
                                br(),
                                actionButton(inputId = "boton_calculate_TSNE1",
                                             label = "Carry out T-SNE after normalization between series",
                                             icon = icon("play", lib = "font-awesome"),
                                             width = "60%"),
                                br(),
                                conditionalPanel(condition = "input.boton_TSNE1!=0",
                                                 plotOutput('tsne1', width = '60%')),
                                br(),
                                actionButton(inputId = "boton_calculate_TSNE2",
                                             label = "Carry out T-SNE after removing outliers and treat batch effect",
                                             icon = icon("play", lib = "font-awesome"),
                                             width = "60%"),
                                br(),
                                conditionalPanel(condition = "input.boton_TSNE!=0",
                                                 plotOutput('tsne2', width = '60%')),
                                
                                ),
                        
                        # Tab 8
                        tabItem(tabName = "codigo",
                                h1("Code"),
                                tags$h4(
                                  "In ", tags$a(href = "https://github.com/jbajo09/covid19-severity", "this repository"),
                                  "you can find the code of this web application.")
                        )
                      ) # Close tabs
                    ) # Close dashboard body
) # Close dashboard page

# Extend size of accepted files (40MB instead of the 5MB - default)
options(shiny.maxRequestSize = 600*1024^2)

server <- function(input, output){
  
  values <- reactiveValues(labels = NULL, labels1=NULL, matrix = NULL, matrix1=NULL, matrix2=NULL, labels_train =NULL, labels_test = NULL, matrix_train = NULL, matrix_test =NULL,  DEGs = NULL, ranking = NULL, DEGsMatrix_rna_batch = NULL, optimalSVM_train = NULL, optimalkNN_train = NULL)
  
  # Server of tab: Data loading ------
  
  observeEvent(input$boton_importar, {
    
    # If files are selected, they are imported
    # Read labels
    labels <- as.vector(t(read.csv2(file = input$file_labels$datapath)))
    DEGsMatrix <- read.csv(input$file_DEGsMatrix$datapath)
    rownames(DEGsMatrix) <- DEGsMatrix[,1]
    DEGsMatrix[,1] <- NULL
    DEGsMatrix <- as.matrix(DEGsMatrix)
    
    # Table
    output$tabla1 <- renderTable({
      if(is.null(input$file_labels$datapath)) return(NULL)
      
      # Message if file is correctly imported
      showModal(modalDialog(
        h3(icon("check-circle", lib = "font-awesome", class = "fa-1x"),
           " File imported"),
        easyClose = TRUE,
        footer = NULL
      ))
      
      tabla_aux <- as.data.frame(table(labels)) %>% rename(Label = labels, Samples = Freq)
      return(tabla_aux)
    })
    
    values$labels <- labels
    values$labels1 <- labels
    values$matrix <- DEGsMatrix
    
    #output$dim <- renderText(dim(values$matrix)) 
   
    
  }) # Close import button
  
  # Server of tab : Preprocessing
  
  
  observeEvent(input$boton_normalization, {
    
    w <- Waiter$new(html = tagList(spin_folding_cube(),
                                   span(br(), br(), br(), h4("Carrying out normalization between arrays..."),
                                        style="color:white;")))  
    w$show()
    
    
    
    matrix <- values$matrix
    values$matrix1 <- matrix
    matrix_scale <- normalizeBetweenArrays(matrix, method = 'scale')
    
    w$hide()
    
    values$matrix <- matrix_scale
  })
  
  observeEvent(input$boton_batch, {
    
    w <- Waiter$new(html = tagList(spin_folding_cube(),
                                   span(br(), br(), br(), h4("Carrying out SVA..."),
                                        style="color:white;")))  
    w$show()
    
    
    
    matrix <- values$matrix
    labels <- values$labels
    values$matrix2 <- matrix
    
    matrix_batch <- batchEffectRemoval(matrix,labels, method = 'sva')
    
    
    w$hide()
    
    values$matrix <- matrix_batch
  })
  
  observeEvent(input$boton_outlier, {
    
    w <- Waiter$new(html = tagList(spin_folding_cube(),
                                   span(br(), br(), br(), h4("Removing outlier..."),
                                        style="color:white;")))  
    w$show()
    
    matrix <- values$matrix
    labels <- values$labels
    outliers <- RNAseqQA(matrix,toRemoval = TRUE, toPNG = FALSE, toPDF = FALSE) 
    w$hide()
    
    output$outlier <- renderText(outliers$outliers)
    
    values$labels <- labels[-which(colnames(matrix) %in% outliers$outliers)]
    values$matrix <- outliers$matrix
    
    # Table
    output$tabla2 <- renderTable({
      
      tabla_aux1 <- as.data.frame(table(values$labels))
      return(tabla_aux1)
    })
    
    output$dim <- renderText(paste0('The final dimension matrix contain information of  ', dim(values$matrix)[1], ' genes for ', dim(values$matrix)[2], ' samples'))
  })
  
 
  # Server of tab: Genes selection ------
  
  w <- Waiter$new(html = tagList(spin_folding_cube(),
                                 span(br(), br(), br(), h4("Running DEGs 10-folds CV..."),
                                      style="color:white;")))
  
  observeEvent(input$boton_genes, {
    
    #DEGs extraction
    w$show
    
    
    
    matrix <- values$matrix
    labels <- values$labels

    #indices <- createDataPartition(labels, p = input$porcentaje_entrenamiento / 100, list = FALSE,  times = 1)
    indices = index
    particion.entrenamiento <- matrix[,indices]
    particion.test <- matrix[,-indices]
    #output$indices <- renderText(indices)
    
    # Labels
    labels_train <- labels[indices]
    labels_test  <- labels[-indices]
    
    w$hide()
    
    w <- Waiter$new(html = tagList(spin_folding_cube(),
                                   span(br(), br(), br(), h4("Running DEGs 5-folds CV..."),
                                        style="color:white;")))
    w$show()
    set.seed(200)
    folds=10
    cvIndex <- createDataPartition(labels_train, p = .80, list = FALSE, times = folds)
    output$indices <- renderText(cvIndex)
    cvResults <- list()
    cvDEGs <- list ()
    for (i in seq_len(folds)){
      cvResults[[i]] <- DEGsExtraction(particion.entrenamiento[,cvIndex[,i]], as.factor(labels_train[cvIndex[,i]]), lfc=1, cov=2, pvalue = 0.05, number = Inf)
      cvDEGs[[i]] <- rownames(cvResults[[i]]$DEG_Results$MulticlassLFC)
    }
    DEGs <- Reduce(f='intersect', cvDEGs) # lfc 1 cov 2 pvalue 0.05 #137 genes
    output$numberDEGs<- renderText(paste0( length(DEGs), " DEGs were extracted "))

    w$hide()
    values$DEGs <- DEGs
    values$matrix_train <- particion.entrenamiento
    values$matrix_test <- particion.test
    values$labels_train <- labels_train
    values$labels_test <- labels_test
    
  })
    
  
  observeEvent(input$boton_ranking, {
    # mRMR method
    w <- Waiter$new(html = tagList(spin_folding_cube(),
                                   span(br(), br(), br(), h4("Running mRMR algorithm..."),
                                        style="color:white;")))
    w$show()
    
    
    particion.entrenamiento <- values$matrix_train
    particion.test <- values$matrix_test
    labels_train <- values$labels_train
    labels_test <- values$labels_test
    DEGs <- values$DEGs

    mrmrRanking <- featureSelection(as.matrix(t(particion.entrenamiento)), as.factor(labels_train), DEGs,
                                              mode = "mrmr")
    mrmrRanking <- names(mrmrRanking)
    w$hide()
    
    values$ranking <- mrmrRanking
    
    # Ranking tables
    output$genes_mrmr <- renderText(mrmrRanking)
    
    #output$genes_mrmr <- renderTable({
      #return(mrmrRanking)
    #}, colnames = FALSE)
    
    
  }) # Close button
  
  # Server of tab: Model training ------
  
  w2 <- Waiter$new(html = tagList(spin_folding_cube(),
                                  span(br(), br(), br(), h4(""),
                                       style="color:white;")))  
  
  observeEvent(input$boton_model_training, {
    
    
    
    particion.entrenamiento <- values$matrix_train
    particion.test <- values$matrix_test
    labels_train <- values$labels_train
    labels_test <- values$labels_test
    ranking <- values$ranking
    
    w3 <- Waiter$new(html = tagList(spin_folding_cube(),
                                    span(br(), br(), br(), h4("Training kNN algorithm..."),
                                         style="color:white;")))  
    w3$show()
    set.seed(2)
    results_cv <- knn_trn(as.matrix(t(particion.entrenamiento)), labels_train, ranking,
                          numFold = as.numeric(input$number_folds))
    values$optimalkNN_train <- results_cv$bestK
    
    w3$hide()
    
    output$optimal_knn <- renderText(paste0("\nOptimal number of neighbours k = ", results_cv$bestK))
    
    output$train <- renderPlot({
      plot(results_cv$accuracyInfo$meanAccuracy[1:20], type = 'l', col= 'black', ylab='Metric Performance', xlab='Genes', lwd=2, ylim = c(0.79,1), panel.first = grid(col='gray45'), cex.axis=1.3,cex.lab=1.3)
      lines(results_cv$sensitivityInfo$meanSensitivity[1:20], col='blue', lwd=2, lty=2)
      lines(results_cv$specificityInfo$meanSpecificity[1:20], col='#FF8B00', lwd=2, lty=4)
      lines(results_cv$F1Info$meanF1[1:20], col='red', lwd=2, lty=4)
      legend(x=15.9 ,y =0.8405, c('Accuracy', 'Sensitivity','Specificity','F1-Score'), lty = c(1,2,4,5), col = c('black','blue','#FF8B00','red'), cex=1.3)
      
      
    })
    
  }) 
  
  w2 <- Waiter$new(html = tagList(spin_folding_cube(),
                                  span(br(), br(), br(), h4(""),
                                       style="color:white;")))  
  observeEvent(input$boton_boxplot, {
    w2$show()
    
    w2$hide()
    
    ranking <- values$ranking
    matrix <- values$matrix
    labels <- values$labels 
    
    
    output$boxplot <- renderPlot({
      dataPlot(matrix[which(rownames(matrix) %in% ranking[1:as.numeric(input$nboxplot)]),],labels,mode = "genesBoxplot",toPNG = FALSE, colours = c("darkred", "forestgreen", "darkorchid3", "dodgerblue3"))
      
    })
  })
  
  observeEvent(input$test_validation, {
  
  
    w3 <- Waiter$new(html = tagList(spin_folding_cube(),
                                    span(br(), br(), br(), h4("Validating kNN algorithm..."),
                                         style="color:white;")))  
    w3$show()
    
    particion.entrenamiento <- values$matrix_train
    particion.test <- values$matrix_test
    labels_train <- values$labels_train
    labels_test <- values$labels_test
    ranking <- values$ranking
    set.seed(200)
    
    
    results_validation <- knn_test(train = t(particion.entrenamiento), labels_train,
                                   test = t(particion.test), labels_test,
                                   ranking, bestK = values$optimalkNN_train)
    
    w3$hide()
    
    
    output$results_validation <- renderPlot({
      tabla <- results_validation$cfMats[[input$numero_genes_validation]]$table
      plotConfMatrix(tabla)
    })
    
  })
  
  # LOOCV
  
  w4 <- Waiter$new(html = tagList(spin_folding_cube(),
                                  span(br(), br(), br(), h4("k-NN LOOCV ..."),
                                       style="color:white;")))  
  
  observeEvent(input$boton_model_LOOCV, {
    
    ranking <- values$ranking
    matrix <- values$matrix
    labels <- values$labels 
    

    w4$show()
   
    results_LOOCV <- knn.cv(t(matrix[which(rownames(matrix) %in% ranking[1:as.numeric(input$numero_genes_LOOCV)]),]), cl = labels, k=values$optimalkNN_train)
    w4$hide()
    
    
    output$results_LOOCV <- renderPlot({
      tabla <- confusionMatrix(data = results_LOOCV, reference = as.factor(labels))$table
      plotConfMatrix(tabla)
   })
    
    
  })
  
  w5 <- Waiter$new(html = tagList(spin_folding_cube(),
                                  span(br(), br(), br(), h4("Calculating AUC..."),
                                       style="color:white;")))  
  
  
  #AUC
  observeEvent(input$boton_calculate_AUC, {
    
    
   
    ranking <- values$ranking
    matrix <- values$matrix
    labels <- values$labels 
    
    w4 <- Waiter$new(html = tagList(spin_folding_cube(),
                                    span(br(), br(), br(), h4("Calculating AUC..."),
                                         style="color:white;")))
    
    w4$show()
    
    response <- as.factor(labels)
    aucs <- rep(NA, length(levels(response))) # store AUCs
    legendLabels <- as.character()
    colours <- c('red','blue','green','black')
    a <- list()
    
    for (i in seq_along(levels(response))) {
      cur.class <- levels(response)[i]
      binaryTraining.labels <- as.factor(labels == cur.class)
      
      
      knn_OVA_LOOCV <- knn.cv(t(matrix[which(rownames(matrix) %in% ranking[1:as.numeric(input$numero_genes_LOOCV)]),]), cl = binaryTraining.labels, k=values$optimalkNN_train)


      score <- knn_OVA_LOOCV
      score <- as.vector(score)
      score[score=='FALSE'] <- 0
      score[score=='TRUE'] <- 1
      binaryTraining.labels  <- as.vector(binaryTraining.labels )
      binaryTraining.labels [binaryTraining.labels =='FALSE'] <- 0
      binaryTraining.labels [binaryTraining.labels =='TRUE'] <- 1
      pred <- prediction(as.numeric(score), as.numeric(binaryTraining.labels ))
      perf <- performance(pred, "tpr", "fpr")
      roc.x <- unlist(perf@x.values)
      roc.y <- unlist(perf@y.values)
      a[[i]] <- roc.x
      a[[i+4]] <- roc.y
      # store AUC
      auc <- performance(pred, "auc")
      auc <- unlist(slot(auc, "y.values"))
      aucs[i] <- auc
      legendLabels[i] <- paste(levels(response)[i], " AUC: ",format(round(aucs[i], 4), nsmall = 3),sep = "")
    }
    
    output$meanAUC <- renderText((paste0("Mean AUC under the precision-recall curve is: ", round(mean(aucs), 2)))) 
    
    w4$hide()
    
    output$AUC <- renderPlot({
      plot(x=NA, y=NA, xlim=c(0,1), ylim=c(0,1),ylab="Sensitivity", xlab="1 - Specificity", bty='n', cex.lab=1.3, cex.axis=1.3)
      lines(x=c(0,1), c(0,1))
      lines(a[[5]] ~ a[[1]], col = colours[1], lwd = 2)
      lines(a[[6]] ~ a[[2]], col = colours[2], lwd = 2)
      lines(a[[7]] ~ a[[3]], col = colours[3], lwd = 2)
      lines(a[[8]] ~ a[[4]], col = colours[4], lwd = 2)
      legend(x=0.38 ,y =0.305, legendLabels, lty=1, ncol= 1,inset = c(0,0),  col = colours, cex = 1.3,lwd=3)
      
      
    })
    
    
  })

  
  # TSNE
  
  observeEvent(input$boton_calculate_TSNE, {
    
    w4 <- Waiter$new(html = tagList(spin_folding_cube(),
                                    span(br(), br(), br(), h4("Carrying out T-SNE ..."),
                                         style="color:white;")))  
    
    ranking <- values$ranking
    matrix <- values$matrix1
    labels <- values$labels1 
    
    
    w4$show()
    
    output$tsne <- renderPlot({
      tsne(matrix[which(rownames(matrix)%in% ranking[1:as.numeric(input$numero_genes_TSNE)]),],labels=as.factor(labels),controlscale=TRUE, scale=3, colvec = c('red','blue','green','black'),seed = 1,axistextsize=0)
    })
    
    w4$hide()

  })
  
  observeEvent(input$boton_calculate_TSNE1, {
    
    w4 <- Waiter$new(html = tagList(spin_folding_cube(),
                                    span(br(), br(), br(), h4("Carrying out T-SNE ..."),
                                         style="color:white;")))  
    
    ranking <- values$ranking
    matrix <- values$matrix2
    labels <- values$labels1 
    
    
    w4$show()
    
    output$tsne1 <- renderPlot({
      tsne(matrix[which(rownames(matrix)%in% ranking[1:as.numeric(input$numero_genes_TSNE)]),],labels=as.factor(labels),controlscale=TRUE, scale=3, colvec = c('red','blue','green','black'),seed = 1,axistextsize=0)
    })
    
    w4$hide()
    
  })
  
  observeEvent(input$boton_calculate_TSNE2, {
    
    w4 <- Waiter$new(html = tagList(spin_folding_cube(),
                                    span(br(), br(), br(), h4("Carrying out T-SNE ..."),
                                         style="color:white;")))  
    
    ranking <- values$ranking
    matrix <- values$matrix
    labels <- values$labels 
    
    
    w4$show()
    
    output$tsne2 <- renderPlot({
      tsne(matrix[which(rownames(matrix)%in% ranking[1:as.numeric(input$numero_genes_TSNE)]),],labels=as.factor(labels),controlscale=TRUE, scale=3, colvec = c('red','blue','green','black'),seed = 1,axistextsize=0)
    })
    
    w4$hide()
    
  })
  

}

shinyApp(ui, server)