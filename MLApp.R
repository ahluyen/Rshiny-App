library(shiny)
#install.packages("shinydashboard")
library(shinydashboard)
library(imputeTS)
library(ggplot2)
library(viridis)
library(caret)
library(rattle)
library(rpart)
library(e1071)
library(class)
library(plotly)
library(DT) 
library(pscl)
library(pls)
library(randomForest)
library(shinycssloaders)

setwd("/Users/andy/Desktop/Spring2020/Shiny")

ui <- dashboardPage(
  dashboardHeader(title = "ML App"),
  dashboardSidebar(
    sidebarMenu(
      menuItem("File input", tabName = "fileinput", icon = icon("tree")),
      menuItem("KNN", tabName = "knn", icon = icon("tree")),
      menuItem("Decision Tree", tabName = "decision_tree", icon = icon("tree")),
      menuItem("Random Forest", tabName = "RandomForest", icon = icon("tree")),
      menuItem("PCR", tabName = "pcr", icon = icon("tree"))
    )
  ),
  dashboardBody(
    tabItems(
      tabItem(tabName = "fileinput",
              fluidPage(
                sidebarLayout(
                  sidebarPanel(
                    titlePanel('Uploading files'),
                    fileInput("file", HTML("Choose CSV File <br/>(Please input a file with the first column as id and last column as pheno type)"),
                              multiple = FALSE,
                              accept = c("text/csv", "text/comma-separated-values,text/plain", ".csv")),
                    selectInput('colMissingPercent', 'Drop colums with selected percentage of missing values', c(100, 90, 80, 70, 60, 50, 40, 30, 20, 10), selected = 50),
                    selectInput('rowMissingPercent', 'Drop rows with selected percentage of missing values', c(100, 90, 80, 70, 60, 50, 40, 30, 20, 10), selected = 90),
                    selectInput('replacement', 'Replace missing values with', c('mode','mean', 'median')),
                    sliderInput('trainPercent', 'Training Set Percentage',min = 1, max = 99, value = 80)),
                  mainPanel(
                    #tableOutput("table")
                    DT::dataTableOutput("table"),
                    #verbatimTextOutput("colname_consolidated"),
                    DT::dataTableOutput("test"),
                    DT::dataTableOutput("train")
                  )
                )
              )
      ),
      tabItem(tabName = "knn",
              fluidPage(
                sidebarLayout(
                  sidebarPanel(
                    titlePanel('K-nearest Neighbors Algorithm'),
                    br(),
                    downloadButton(outputId = "downloadMatrix", label = "Download Confusion Matrix")
                  ),
                  mainPanel(
                    plotOutput('plot1') %>% withSpinner(color="#0dc5c1"),
                    verbatimTextOutput(outputId = "matrix") %>% withSpinner(color="#0dc5c1")
                  )
                )
              )
      ),
      tabItem(tabName = "decision_tree",
              fluidPage(
                titlePanel(title = "RPART Decision Tree"
                ),
                sidebarLayout(
                  sidebarPanel(
                    selectizeInput('predictor1', 
                                   'Choose predictors from list: ', 
                                   choices = c(""), 
                                   options = list(maxItems = 20), 
                                   multiple = TRUE)
                  ),
                  mainPanel(
                    
                    plotOutput("plot2"),
                    verbatimTextOutput("cm_dt"),
                    verbatimTextOutput("summary_dt")
                  )
                )
              )
      ),
      tabItem(tabName = "RandomForest", 
              fluidPage(
                titlePanel('Random Forest'), 
                sidebarLayout(
                  sidebarPanel(
                    radioButtons(inputId = "rf_predictor", label = "Select predictors: ",
                                 choices = c("All predictor varialbes", "Choose predictor variables"
                                 ), selected = "Choose predictor variables"),
                    conditionalPanel(condition = "input.rf_predictor == 'Choose predictor variables'",
                                     selectizeInput(inputId = "predictor2",
                                                    label = "Choose predictors from list: ",
                                                    options = list(maxItems = 20),
                                                    multiple = TRUE,
                                                    choices = c("") )),
                    h4("Click to save confusion matrix table: "),
                    downloadButton(outputId = "download_cm", label = "Download")
                  ),
                  mainPanel(
                    verbatimTextOutput("rf_cm"),
                    verbatimTextOutput("rf_sum")
                  )
                )
              )
      ),
      tabItem(
        tabName = "pcr",
        fluidPage(
          titlePanel("Principal Component Regression"),
          sidebarLayout(
            sidebarPanel(
              selectizeInput(inputId = "pca_no",
                             label = "Choose the No. of PCs", 
                             choices = c(''),
                             options = list(create = TRUE), multiple = TRUE)
            ),
            mainPanel(
              verbatimTextOutput("pcr_cm"),
              verbatimTextOutput("pcr_fit"),
              #plotOutput("pcr_plot"),
              DT::dataTableOutput("PCR_train_label"),
              #verbatimTextOutput("PCR_train_colname"),
              DT::dataTableOutput("PCR_test")
              #verbatimTextOutput("PCR_test_colname")
              
            )
          )
        )
      )
    )
  )
)

server <- function(input, output, session){
  dataset <- reactive({
    # default file if no user input
    filePath = "two_label_with_selected_features_rn_v3.csv"
    if (!is.null(input$file)) {
      filePath <- input$file$datapath
    }
    df <- read.csv(filePath)
    # check file format
    colnames<- names(df)
    if (colnames[1] != 'names' || colnames[ncol(df)] != 'pheno') {
      return("Input file formate error: Please make sure the first column of the input file is names and the last column is pheno.")
    }
    df <- df[,2:length(df)]
    #df$pheno <- ifelse(df$pheno == 'intolerant', 1, 0)
    # removing % of col containing NAs not required
    df_remove_col <- df[ , which(colMeans(!is.na(df)) > 1 - colPercent())]
    # removing % of rows containing NAs not required
    cleaned_data<-df_remove_col[which(rowMeans(!is.na(df_remove_col)) > 1 - rowPercent()),]
    return (cleaned_data)
  })
  
  # number of column consisting NAs to be dropped saved in colPercent()
  colPercent <- reactive({as.numeric(input$colMissingPercent) / 100})
  # number of rows consisting NAs to be dropped saved in rowPercent()
  rowPercent <- reactive({as.numeric(input$rowMissingPercent) / 100})
  # Imputation method saved in replace()
  replace <- reactive({input$replacement})
  # Incorporating the training percantage
  percent <- reactive({as.numeric(input$trainPercent) / 100})
  
  #Imputing the entire dataset only for table display
  dataset_table <- reactive(na_mean(dataset(), option = replace(), maxgap = Inf))
  
  #Display the cleaned data
  output$table <- DT::renderDataTable(
    dataset_table(), options = list(scrollX = TRUE), caption = "LactoseIntolerance SNP Dataset")
  
  #row numbers to be included in train dataset saved in split()
  # percent() is the value of train test split selected by user
  split <- reactive({
    set.seed(100)
    caret::createDataPartition(dataset()$pheno, p = percent(), list = FALSE)
  })
  
  # for X_train and y_train dataset to be used without the labels column
  data_without_labs <- reactive(dataset()[, 2:length(dataset()) - 1])
  # only the label column for X_test and y_test
  data_labels <- reactive(dataset()['pheno'])
  
  # train and test split for dataset
  X_train <- reactive({data_without_labs()[split(),]})
  X_test <- reactive({data_without_labs()[-split(),]})
  
  #train and test split for label
  y_train <- reactive({data_labels()[split(),]})
  y_test <- reactive({data_labels()[-split(),]})
  
  #Imputed training dataset
  X_train_impute <- reactive(na_mean(X_train(), option = replace(), maxgap = Inf))
  #Imputed testing dataset
  X_test_impute <- reactive(na_mean(X_test(), option = replace(), maxgap = Inf))
  
  #Binding the imputed train and label together
  
  X_train_consolidated <- reactive(cbind(X_train_impute(), y_train()))
  
  #colnames(X_train_consolidated())[length(X_train_consolidated())] <- "phenotype"
  
  X_train_consolidated_phenotype <- reactive({X_hopingthisworks = X_train_consolidated()
  colnames(X_hopingthisworks)[which(names(X_hopingthisworks)=='y_train()')]='phenotype'
  return(X_hopingthisworks)
  })
  
  X_test_consolidated_phenotype <- reactive({X_hopingthisworks = X_test_consolidated()
  colnames(X_hopingthisworks)[which(names(X_hopingthisworks)=='y_test()')]='phenotype'
  return(X_hopingthisworks)
  })
  
  
  #<- reactive(colnames(X_train_consolidated())[length(X_train_consolidated())] <- "phenotype")
  
  #output$colname_consolidated <- renderPrint({colnames(X_train_consolidated())})
  
  output$test <- DT::renderDataTable(
    X_test(), options = list(scrollX = TRUE), caption = "LactoseIntolerance Test Dataset")
  
  output$train <- DT::renderDataTable(
    X_train_consolidated_phenotype(), options = list(scrollX = TRUE), caption = "LactoseIntolerance Train Dataset" )
  
  #Train your model on X_train_consolidated
  #Test your model on X_test_impute
  # prediction = model.fit(X_test)
  # predictions vs actual (y_test) will be used for cm
  # use labels columns saved in y_test for confusion matrix(cm)
  
  # Code for KNN here
  
  knnFit <- reactive({
    trctrl <- trainControl(method = "repeatedcv", number = 10, repeats = 3)
    set.seed(3333)
    train(phenotype ~ ., data = X_train_consolidated_phenotype(), method = "knn", trControl=trctrl, preProcess = c("center","scale"), tuneLength = 20)
  })
  output$plot1 <- renderPlot({
    plot(knnFit())
  })
  
  # Confusion Matrix Generating
  conMatrix<-reactive({
    knnPredict <- predict(knnFit(),newdata = X_test_impute() )
    confusionMatrix(knnPredict, reference = y_test())
  })
  
  # Show matrix in the sde bar
  output$matrix <- renderPrint({
    m<-conMatrix()
    m
  })
  
  # Download confusion matrix as csv file
  output$downloadMatrix <- downloadHandler(
    filename = function() {
      paste("confusion_matrix_", Sys.Date(), ".csv", sep="")
    },
    content = function(file) {
      #write.csv(as.data.frame.matrix(conMatrix()) , file)
      write.csv(as.data.frame.matrix(as.table(conMatrix())) , file)
    }
  )
  
  # Code for Decision Tree here
  
  observe({
    DT <- names(X_train_impute())
    updateSelectizeInput(session, "predictor1",
                         choices = DT,
                         selected = DT[1:8],
                         options = list(maxItems = 20))
  })
  
  formula_dt <- reactive({
    as.formula(paste("phenotype",
                     paste(input$predictor1, sep = ",", collapse = " + "),
                     sep = " ~ "))
  })
  
  # fit the RPART model per chosen predictors    
  modelfit <- reactive({
     fit_dt <- caret::train(data = X_train_consolidated_phenotype(), formula_dt(), method = "rpart",
          control = rpart.control(minsplit = 1, 
                                  minbucket = 1, 
                                  cp = 0.001))
    return(fit_dt)
  })
  
  # generate predictions
  #predictfit <- reactive({
  #  fit = modelfit()
  #  fit_dt_2 <- predict(fit, newdata = X_test_impute())
  #  return(fit_dt_2)
  #})
  
  output$plot2 <- renderPlot({
    fit = modelfit()
    fancyRpartPlot(fit$finalModel)
  })
  
  # Confusion matrix and Summary for Decision Tree
  
  dt_cm<- reactive({  
    fit.dt = modelfit()
    test_pred_dt <- predict(fit.dt, X_test_impute())
    cm <- confusionMatrix(test_pred_dt, reference = y_test())
    return(cm)
  })
  
  output$summary_dt <- renderPrint({
    modelfit()
  })
  
  output$cm_dt <- renderPrint({
    dt_cm()
  })
  
  # Code for Random Forest here
  
  observe({
    X <- names(X_train_impute())
    updateSelectizeInput(session, "predictor2",
                         choices = X,
                         selected = X[1:20],
                         options = list(maxItems = 20))
  })
  
  formula_rf <- reactive({
    if (input$rf_predictor == 'Choose predictor variables') {
      as.formula(paste("phenotype",
                       paste(input$predictor2, sep = ",", collapse = " + "),
                       sep = " ~ "))
    } else {
      as.formula(paste("phenotype", "~ ."))
    }
  })
  
  rf.fit <- reactive({
    rf_fit <- randomForest(formula_rf(),
                           data = X_train_consolidated_phenotype(),
                           importance = TRUE)
    return(rf_fit)
  })
  
  rf.cm <- reactive({
    fit = rf.fit()
    test_pred <- predict(fit, X_test_impute())
    cm <- confusionMatrix(test_pred, reference = y_test())
    return(cm)
    
  })
  
  output$rf_sum <- renderPrint({
    rf.fit()
  })
  
  output$rf_cm <- renderPrint({
    rf.cm()
  })
  
  output$download_cm <- downloadHandler(
    filename = function() {
      paste("confusionMatrix-", Sys.Date(), ".csv", sep = "")
    },
    content = function(file) {
      write.csv(as.data.frame.matrix(as.table(rf.cm())) , file)
    }
  )
  
  
  # Code for PCR here
  
  #centering, scaling and creating the PCs training dataset
  pr.out.train <- reactive({
    prcomp(X_train_impute(), center = T, scale = T)
  })
  
  pca.x.withoutnames.train <- reactive({
    as.data.frame(pr.out.train()$x[,1:20])
  })
  
  observe({
    Q <- names(pca.x.withoutnames.train())
    updateSelectizeInput(session, "pca_no",
                         choices = Q,
                         selected = Q[1],
                         options = list(maxItems = 20))
  })
  
  #Display the PC train data
  #output$PCR_train <- DT::renderDataTable(
  #  pca.x.withoutnames.train(), options = list(scrollX = TRUE))
  
  #centering and scaling and creating PCs for test dataset as well
  #pr.out.test <- reactive({
  #  prcomp(X_test_impute(), center = T, scale = T)
  #})
  
  #pca.x.withoutnames.test <- reactive({
  #  as.data.frame(pr.out.test()$x[,1:8])
  #})
  
  # Doing PCA on the test set
  pr.out.test <- reactive({
    predict(pr.out.train(), X_test_impute())
  })
  
  #Display the PC test data
  output$PCR_test <- DT::renderDataTable(
    pr.out.test(), options = list(scrollX = TRUE), caption = "LactoseIntolerance PCA Test Dataset")
  
  # Adding the label column to PCs dataset
  X_train_consolidated_PC <- reactive(cbind(pca.x.withoutnames.train(), y_train()))
  
  # Changing names of y-train on PCs as well
  X_train_consolidated_PC_phenotype <- reactive({X_hopingthisworks_2 = X_train_consolidated_PC()
  colnames(X_hopingthisworks_2)[which(names(X_hopingthisworks_2)=='y_train()')]='phenotype'
  return(X_hopingthisworks_2)
  })
  
  #Display the PC train data
  output$PCR_train_label <- DT::renderDataTable(
    X_train_consolidated_PC_phenotype(), options = list(scrollX = TRUE), caption = "LactoseIntolerance PCA Train Dataset")
  
  formula_pcr <- reactive({as.formula(paste('phenotype',
                                        paste(input$pca_no, sep = ",", collapse = "+"),
                                        sep = "~"))
  })
  
  #pcr_model <- reactive({
  #  pcr_model <- glm(formula(), data = X_train_consolidated_PC_phenotype(), family = binomial())
  #  #pcr_model<- glm(y_train() ~ PC1+PC2, X_train_consolidated_PC(), family = binomial())
  #  return(pcr_model)
  #})
  
  pcr_model <- reactive({
    fit <- caret::train(formula_pcr(),
                        data = X_train_consolidated_PC_phenotype(),
                        method = 'glm', family = 'binomial')
    return(fit)
  })
  
  pcr_pred <- reactive({
    pcr_test_pred <- predict(pcr_model(), pr.out.test())
    return(pcr_test_pred)
  })
  
  pcr_cm <- reactive({
  cm = confusionMatrix(pcr_pred(), reference = y_test())
  return(cm)
  })
  
  
  
  #Output for confusion matrix
  output$pcr_cm <- renderPrint({pcr_cm()})
  
  #Output for summary table
  output$pcr_fit <- renderPrint({
    summary(pcr_model())
  })
} 

shinyApp(ui = ui, server = server)