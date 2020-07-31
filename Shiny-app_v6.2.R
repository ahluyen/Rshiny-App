library(shiny)
#install.packages("shinydashboard")
library(shinydashboard)
library(imputeTS)
library(dplyr)
library(ggplot2)
library(viridis)
library(plotly)
library(DT)
library(psych)
library(caret)
library(scales)

#library(rsconnect)

setwd("/Users/andy/Desktop/Spring2020/Shiny")

ui <- dashboardPage(
  dashboardHeader(title = "EDA"),
  dashboardSidebar(
    sidebarMenu(
      menuItem("File input", tabName = "fileinput", icon = icon("tree")),
      menuItem("K-Means", tabName = "k-means", icon = icon("tree")),
      menuItem("Hierarchy Clustering", tabName = "hierarchy", icon = icon("tree")),
      menuItem("PCA", tabName = "pca", icon = icon("tree")),
      menuItem("3D-PCA", tabName = "3d-pca", icon = icon("tree"))
    )
  ),
  dashboardBody(
    tabItems(
      tabItem(tabName = "fileinput",
              tabsetPanel(
                tabPanel("Data",
                         sidebarLayout(
                           sidebarPanel(
                             h4('Uploading files'),
                             fileInput("file", HTML("Choose CSV File <br/>(Please input a file with the first column as id and last column as pheno type)"),
                                       multiple = FALSE,
                                       accept = c("text/csv", "text/comma-separated-values,text/plain", ".csv")),
                             htmlOutput(outputId = "dim"),
                             br(),
                             selectInput('colMissingPercent', 'Drop colums with selected percentage of missing values',
                                         c(100, 90, 80, 70, 60, 50, 40, 30, 20, 10),
                                         selected = 50),
                             selectInput('rowMissingPercent', 'Drop rows with selected percentage of missing values',
                                         c(100, 90, 80, 70, 60, 50, 40, 30, 20, 10),
                                         selected = 90),
                             htmlOutput(outputId = "dim_remove"),
                             br(),
                             htmlOutput(outputId = "dim_after_remove"),
                             br(),
                             downloadButton(outputId = "download", label = "Download lactose intolerance data")
                             ),
                           mainPanel(
                             #tableOutput("table")
                             DT::dataTableOutput("dt")
                           )
                         )
                         ),
                tabPanel("Imputation",
                         sidebarLayout(
                           sidebarPanel(
                             selectInput('replacement', 'Replace missing values with', c('mode','mean', 'median'))
                           ),
                           mainPanel(
                             dataTableOutput("table")
                           )
                         )
                         ),
                tabPanel("Summary",
                         sidebarLayout(
                           sidebarPanel(
                             selectInput(inputId = "data_sum",
                                         label = "Select the dataset",
                                         choices = c("original data", "filtered & imputed data"),
                                         selected = "filtered & imputed data"),
                             numericInput(inputId = "from",
                                          label = "Show data summary from column",
                                          value = 1, min = 1),
                             numericInput(inputId = "to",
                                          label = "to column",
                                          value = 50, min = 1)),
                           mainPanel(
                             tags$style(type = "text/css",
                                        ".shiny-output-error { visibility: hidden; }",
                                        ".shiny-output-error:before { visibility: hidden; }"),
                             verbatimTextOutput(outputId = "summary"))
                         )),
                tabPanel("Structure",
                         sidebarLayout(
                           sidebarPanel(
                             selectInput(inputId = "data_str",
                                         label = "Select the dataset",
                                         choices = c("original data", "filtered & imputed data"),
                                         selected = "filtered & imputed data"),
                             numericInput(inputId = "start",
                                          label = "Show data structure from column",
                                          value = 1, min = 1),
                             numericInput(inputId = "end",
                                          label = "to column",
                                          value = 50, min = 1)),
                           mainPanel(
                             verbatimTextOutput(outputId = "structure"))
                         )),
                tabPanel("Plot",
                         sidebarLayout(
                           sidebarPanel(
                             selectInput(inputId = "data_plot",
                                         label = "Select the dataset",
                                         choices = c("original data", "filtered & imputed data"),
                                         selected = "filtered & imputed data"),
                             selectizeInput(inputId = "plot_var",
                                            label = "Select the variable",
                                            choices = c(""),
                                            options = list(maxItems = 1)),
                             checkboxInput(inputId = "color",
                                           label = "color by phenotype", TRUE)
                           ),
                           mainPanel(
                             plotOutput(outputId = "bar"))
                         )
                         )
              )

      ),
      tabItem(tabName = "pca",
              fluidPage(
                sidebarLayout(
                  sidebarPanel(
                    titlePanel('Principal Component Analysis'),
                    selectizeInput('xcol', 'X Variable', choices = c('PC1','PC2','PC3','PC4','PC5','PC6','PC7','PC8'),
                                   options = list(create = TRUE), selected = c('PC1')),
                    selectizeInput('ycol', 'Y Variable', choices = c('PC1','PC2','PC3','PC4','PC5','PC6','PC7','PC8'),
                                   options = list(create = TRUE), selected = c('PC2')),
                    hr(),
                    checkboxInput("col", "Color by Phenotype", FALSE),
                    checkboxInput("stat_ellipse", "Add Ellipse", FALSE)
                  ),
                  mainPanel(
                    tabsetPanel(id = "tabs",
                                tabPanel("PCA Plot",
                                         plotlyOutput("plot1")
                                ),
                                tabPanel("Scree Plot",
                                         plotOutput("plot5"),
                                         htmlOutput(outputId = "scree_title"),
                                         verbatimTextOutput("pca_var")
                                ),
                                tabPanel("Cumulative Plot",
                                         plotOutput("plot6"),
                                         htmlOutput(outputId = "cumu_title"),
                                         verbatimTextOutput("pca_cumu_var")
                                )
                    )
                  )
                )
              )
      ),
      tabItem(tabName = "k-means",
              fluidPage(
                sidebarLayout(
                  sidebarPanel(
                    titlePanel('K-means Clustering'),
                    selectizeInput('xcol2', 'X Variable', choices = c(""), options = list(maxItems = 1)),
                    selectizeInput('ycol2', 'Y Variable', choices = c(""), options = list(maxItems = 1)),
                    hr(),
                    checkboxInput('kmeans_colcluster', 'Color by Clusters', FALSE),
                    numericInput('kmeans_clusters', 'Cluster count', 2,
                                 min = 1, max = 9),
                    checkboxInput('kmeans_pheno', 'Color by Pheno', FALSE)),
                  mainPanel(
                    verbatimTextOutput(outputId = "km_table"),
                    plotlyOutput("plot2")
                  )
                )
              )
      ),
      tabItem(tabName = "3d-pca",
              fluidPage(
                sidebarLayout(
                  sidebarPanel(
                    titlePanel('3D PCA Plot'),
                    selectizeInput('xcol3', 'X Variable', choices = c('PC1','PC2','PC3','PC4','PC5','PC6','PC7','PC8'),
                                   options = list(create = TRUE), selected = c('PC1')),
                    selectizeInput('ycol3', 'Y Variable', choices = c('PC1','PC2','PC3','PC4','PC5','PC6','PC7','PC8'),
                                   options = list(create = TRUE), selected = c('PC2')),
                    selectizeInput('zcol3', 'Z Variable', choices = c('PC1','PC2','PC3','PC4','PC5','PC6','PC7','PC8'),
                                   options = list(create = TRUE), selected = c('PC3'))
                  ),
                  mainPanel(
                    plotlyOutput("plot3")
                  )
                )
              )
      ),
      tabItem(tabName = "hierarchy",
              fluidPage(
                sidebarLayout(
                  sidebarPanel(
                    titlePanel('Hierarchical clustering'),
                    selectInput('distance', 'Distance Method', c("euclidean", "maximum", "manhattan", "canberra", "binary", "minkowski")),
                    selectInput('linkage', 'Linkage Method', c("ward.D", "ward.D2", "single", "complete",
                                                               "average", "mcquitty", "median", "centroid"),
                                selected = "ward.D2"),
                    radioButtons("hcoloring", "Color option:",
                                 c("Default color" = "na",
                                   "Colored by clusters" = "cluster",
                                   "Colored by phenotype" = "pheno")),
                    conditionalPanel(
                      condition = "input.hcoloring == 'cluster'",
                      selectInput("hclusterNum", "Cluster Number", c(1,2,3,4,5,6,7,8,9))
                    )
                  ),
                  mainPanel(
                    verbatimTextOutput(outputId = "hc_table"),
                    plotOutput('plot4')
                  )
                )
              )
      )
    )
  )
)

server <- function(input, output, session){
  options(shiny.sanitize.errors = TRUE)

  dataset <- reactive({
    # default file if no user input
    df <- read.csv("two_label_with_selected_features_rn_v3.csv")
    if (!is.null(input$file)) {
      df <- read.csv(input$file$datapath)
    }
    colnames <- names(df)
    if (colnames[1] != 'names' || colnames[ncol(df)] != 'pheno') {
      return("Input file formate error: Please make sure the first column of the input file is names and the last column is pheno.")
    }
    return(df)
  })

  output$dim <- renderText({paste("<b>This data has:</b>", "<br>", dim(dataset())[1],
                                  "rows", "<br>", dim(dataset())[2], "columns")})

    output$download <- downloadHandler(
    filename = function() {
      paste("Lactose_intolerant_data_", Sys.Date(), ".csv", sep = "")
    },
    content = function(file) {
      write.csv(dataset(), file)
    }
  )

  data_removed <- reactive({
    df2 <- dataset()
    df_filtered <- df2[ , 2:length(df2)]
    df_remove_col <- df_filtered[ , which(colMeans(!is.na(df_filtered)) > 1 - colPercent())]
    cleaned_data <- df_remove_col[which(rowMeans(!is.na(df_remove_col)) > 1 - rowPercent()),]
    return(cleaned_data)
  })

  output$dt <- renderDataTable(
    data_removed(), options = list(scrollX = TRUE))

  output$dim_remove <- renderText({
      paste(dim(dataset())[2] - dim(data_removed())[2] - 1, "<b>columns</b> and",
            dim(dataset())[1] - dim(data_removed())[1], "<b>rows</b> have been removed")
  })

  output$dim_after_remove <- renderText({
    paste("<b>The data currently has:</b>", "<br>", dim(data_removed())[1],
            "rows", "<br>", dim(data_removed())[2], "columns")
  })

  # Remove x% of rows and columns
  colPercent <- reactive({as.numeric(input$colMissingPercent) / 100})
  rowPercent <- reactive({as.numeric(input$rowMissingPercent) / 100})
  replace <- reactive({input$replacement})


  # Filtered Dataset
  dataset_filtered <- reactive({
    df2 <- dataset()
    df_filtered <- df2[ , 2:length(df2)]
    df_remove_col <- df_filtered[ , which(colMeans(!is.na(df_filtered)) > 1 - colPercent())]
    cleaned_data <- df_remove_col[which(rowMeans(!is.na(df_remove_col)) > 1 - rowPercent()),]
    dataset2 <- na_mean(cleaned_data, option = replace(), maxgap = Inf)
    return(dataset2)
  })

  output$table <- DT::renderDataTable(
    dataset_filtered(), options = list(scrollX = TRUE))

  observe({
    if (input$data_sum == "filtered & imputed data") {
      x <- length(dataset_filtered())
      updateNumericInput(session, "to", value = x, max = x)

    } else {

      x <- length(dataset())
      updateNumericInput(session, "to", value = x, max = x)
      }
    })

  # Summary
  output$summary <- renderPrint({
    if (input$data_sum == "original data") {
      summary(dataset()[, input$from:input$to])
    } else {
      summary(dataset_filtered()[, input$from:input$to])
    }
  })

  observe({
    if (input$data_str == "filtered & imputed data") {
      x <- length(dataset_filtered())
      updateNumericInput(session, "end", value = x, max = x)

    } else {

      x <- length(dataset())
      updateNumericInput(session, "end", value = x, max = x)
    }
  })

  # Structure
  output$structure <- renderPrint({
    if (input$data_str == "original data") {
      str(dataset()[, input$start:input$end])
    } else {
      str(dataset_filtered()[, input$start:input$end])
    }
  })

  # Plot
  observe({
    if (input$data_plot == "filtered & imputed data") {

    plot_var <- colnames(dataset_filtered()[, 1:(length(dataset_filtered()) - 1)])

    # Can also set the label and select items
    updateSelectizeInput(session, "plot_var",
                         choices = plot_var,
                         options = list(maxItems = 1)
    )
  } else {

      plot_var <- colnames(dataset()[, 2:(length(dataset()) - 1)])
      updateSelectizeInput(session, "plot_var",
                           choices = plot_var,
                           options = list(maxItems = 1))
  }

  })

  dataset_filtered_fc <- reactive(as.data.frame(lapply(dataset_filtered(), as.factor)))

  dataset_fc <- reactive(as.data.frame(lapply(dataset(), as.factor)))

  dataset_fc_pheno <- reactive({
    if (class(dataset()[, length(dataset())]) == "factor") {
      dataset()[, length(dataset())]
    } else {
      as.factor(dataset()[, length(dataset())])
    }
  })

  output$bar <- renderPlot({
    if (input$data_plot == "filtered & imputed data") {
      #data: remove and impute
      g <- ggplot(data = dataset_filtered_fc(), aes(y = (..count..)/sum(..count..))) +
        scale_y_continuous(labels = percent)
      p <- g + geom_bar(aes_string(x = input$plot_var), fill = "steelblue",
                        width = 0.5, alpha = 0.8) + ylab("Percentage")

      if (input$color)
        p <- g + geom_bar(aes_string(x = input$plot_var, fill = phenotype()),
                          width = 0.5, alpha = 0.8, position = "dodge") +
        labs(fill = "Class") + ylab("Percentage")
      p

    } else {
      #data: remove
      g <- ggplot(data = dataset_fc(), aes(y = (..count..)/sum(..count..))) +
        scale_y_continuous(labels = percent)
      p <- g + geom_bar(aes_string(x = input$plot_var), fill = "steelblue",
                        width = 0.5, alpha = 0.8) + ylab("Percentage")

      if (input$color)
        p <- g + geom_bar(aes_string(x = input$plot_var, fill = dataset_fc_pheno()),
                          width = 0.5, alpha = 0.8, position = "dodge") +
        labs(fill = "Class") + ylab("Percentage")
      p

    }

  })


  # Extract pheno
  phenotype <- reactive({
    if (class(dataset_filtered()[, length(dataset_filtered())]) == "factor") {
      dataset_filtered()[, length(dataset_filtered())]
    } else {
      as.factor(dataset_filtered()[, length(dataset_filtered())])
    }
  })

  data_without_labs <- reactive(dataset_filtered()[, 1:length(dataset_filtered()) - 1])

  # PCA
  pr.out <- reactive({
    prcomp(data_without_labs(), center = TRUE, scale. = TRUE)
  })

  pca.x <- reactive({
    as.data.frame(pr.out()$x[,1:8])
    })

  output$plot1 <- renderPlotly({
    g <- ggplot(pca.x(), aes_string(input$xcol, input$ycol)) +
      geom_point(shape = 20, alpha = 0.5) +
      ggtitle("PCA Plot") +
      xlab(input$xcol) + ylab(input$ycol)
    if (input$col){
      g <- g + aes(col = phenotype(), fill = phenotype())
    }

    if (input$stat_ellipse){
      g <- g + stat_ellipse(geom = "polygon", alpha = 0.3)
    }
    ggplotly(g)
  })

  # Variability of each principal component: pr.var
  pr.var <- reactive({pr.out()$sdev^2})
  # Variance explained by each principal component: pve
  pve <- reactive({pr.var()/sum(pr.var())})

  output$plot5 <- renderPlot({
      plot(pve(), xlab = "Principal Component",
           ylab = "Proportion of Variance Explained",
           main = "PCA Scree Plot",
           ylim = c(0, 1), type = "b")
  })

  output$plot6 <- renderPlot({
      plot(cumsum(pve()), xlab = "Principal Component",
           ylab = "Cumulative Proportion of Variance Explained",
           main = "PCA Cumulative Plot",
           ylim = c(0, 1), type = "l")
  })

  output$scree_title <- renderText({
      paste(h5("The proportion of variance explained by each principal component"))
    })

  pca.summary <- reactive({summary(pr.out())})

  output$pca_var <- renderPrint({
    pca.summary()$importance[2, ]
  })

  output$cumu_title <- renderText({
      paste(h5("The cumulative proportion of variance explained"))
  })

  output$pca_cumu_var <- renderPrint({
    pca.summary()$importance[3, ]
  })


  observe({
    xcol2 <- names(data_without_labs()[,1:length(data_without_labs())])
    updateSelectizeInput(session, "xcol2",
                         choices = xcol2,
                         selected = "rs4988235",
                         options = list(maxItems = 1))
  })

  observe({
    ycol2 <- names(data_without_labs()[,1:length(data_without_labs())])
    updateSelectizeInput(session, "ycol2",
                         choices = ycol2,
                         options = list(maxItems = 1),
                         selected = "rs182549")
  })

  kmeans_countclusters <- reactive({
    kmeans(data_without_labs(), input$kmeans_clusters)
  })


  output$km_table <- renderPrint({
    table(kmeans_countclusters()$cluster, phenotype())
  })

  km_cluster <- reactive({as.factor(kmeans_countclusters()$cluster)})

  km_center <- reactive({as.data.frame(kmeans_countclusters()$centers)})

  output$plot2 <- renderPlotly({
    g1 <- ggplot(data_without_labs(), aes_string(input$xcol2, input$ycol2)) +
      geom_jitter(position = position_jitter(width = 0.5, height = 0.5), shape = 20, alpha = 0.5) +
      ggtitle("K-means Clustering Plot") +
      xlab(input$xcol2) + ylab(input$ycol2)
    if (input$kmeans_colcluster){
      g1 <- g1 + aes(col = km_cluster(), fill = km_cluster()) +
        annotate("point", x = km_center()[,1], y = km_center()[,2],
                 size = 5, colour = 'red', alpha = 0.5)
    }
    if (input$kmeans_pheno){
      g1 <- g1 + aes(col = phenotype(), fill = phenotype())
    }
    ggplotly(g1)
  })

  x <- reactive({
    pca.x()[, c(input$xcol3)]
  })

  y <- reactive({
    pca.x()[, c(input$ycol3)]
  })

  z <- reactive({
    pca.x()[, c(input$zcol3)]
  })

  output$plot3 <- renderPlotly({
    plot_ly(x = x(), y = y(), z = as.matrix(z()),
            color = phenotype(),
            colors = c("#F8766D", "#7CAE00", "#C77CFF", "#00BFC4"),
            type = 'scatter3d',
            alpha = 0.9, opacity = 0.6) %>%
      add_markers(showlegend = FALSE) %>%
      layout(
        title = paste(input$xcol3, "vs", input$ycol3, "vs", input$zcol3),
        scene = list(
          xaxis = list(title = input$xcol3),
          yaxis = list(title = input$ycol3),
          zaxis = list(title = input$zcol3))
      )
  })

  distance <- reactive({
    dist(data_without_labs(), method = input$distance)
  })

  clusters <- reactive({
    hclust(distance(), method = input$linkage)
  })

  clust <- reactive({cutree(clusters(), k = input$hclusterNum)})

  output$hc_table <- renderPrint({
    table(clust(), phenotype())
  })

  output$plot4 <- renderPlot({
    hcd <- as.dendrogram(clusters())
    clusterColors <- as.vector(viridis(9))
    labelColors <- c("red", "blue")

    if (input$hcoloring == 'cluster') {
      # Color according to clusters
      clusterCut <- cutree(clusters(), input$hclusterNum)

      # function to get color labels
      colLab <- function(n) {
        if (is.leaf(n)) {
          a <- attributes(n)
          labCol <- labelColors[clusterCut[which(names(clusterCut) == a$label)]]
          attr(n, "nodePar") <- c(a$nodePar, lab.col = labCol)
        }
        n
      }
      clusDendro <- dendrapply(hcd, colLab)
      plot(clusDendro, type = "rectangle", ylab = "Height")

    } else if (input$hcoloring == 'pheno') {
      # function to get color labels
      colLab <- function(n) {
        if (is.leaf(n)) {
          a <- attributes(n)
          index <- as.numeric(as.character(a$label))
          label <- as.matrix(phenotype())[index, ncol(as.matrix(phenotype()))]
          col <- "red"
          if (label == "tolerant") {
            col <- "blue"
          }
          attr(n, "nodePar") <- c(a$nodePar, lab.col = col)
        }
        n
      }
      clusDendro <- dendrapply(hcd, colLab)
      plot(clusDendro, type = "rectangle", ylab = "Height")
      legend("topright",
             legend = c("Intolerant" , "Tolerant"),
             col = c("red", "blue"),
             pch = c(20,20,4,4,4), bty = "n",  pt.cex = 1.5, cex = 0.8 ,
             text.col = "black", horiz = FALSE, inset = c(0, 0))
    } else {
      # No coloring
      plot(hcd, type = "rectangle", ylab = "Height")
    }
  })
}

shinyApp(ui = ui, server = server)
