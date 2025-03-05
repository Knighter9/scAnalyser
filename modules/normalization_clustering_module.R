library(Seurat)
library(ggplot2)
library(shiny)
library(shinyjs)
library(future)
library(promises)

# Set up multisession plan for futures
plan(multisession)

# A helper function to create a box with a plot and download button
plotBox <- function(plotOutputId, downloadButtonId, plotTitle, ns,explanationText=NULL) {
  div(
    style = "border: 1px solid #ccc; padding: 10px; margin: 5px;",
    h4(plotTitle),
    plotOutput(ns(plotOutputId)),
    downloadButton(ns(downloadButtonId), paste0("Download ", plotTitle, " Plot")),
    if(!is.null(explanationText)){
        p(style = "font-size: 0.8em; margin-top: 10px;", explanationText)
    }
  )
}

NormalizationClusteringUI <- function(id) {
  ns <- NS(id)
  tagList(
    useShinyjs(),
    tags$head(
      tags$style(HTML("
                .custom-spinner {
                    margin: 0 auto;
                    width: 50px;
                    height: 50px;
                    border: 5px solid #f3f3f3;
                    border-top: 5px solid #007bff;
                    border-radius: 50%;
                    animation: spin 1s linear infinite;
                }
                @keyframes spin {
                    0% { transform: rotate(0deg); }
                    100% { transform: rotate(360deg); }
                }
      "))
    ),
    div(style = "text-align: center; margin-top: 30px;",
        h2("Normalization and Clustering"),
        radioButtons(ns("normalization_method"), "Normalization Method:",
                     choices = c("Standard Seurat Workflow" = "standard",
                                 "SCTransform" = "sctransform"),
                     inline = TRUE),
        br(),
        actionButton(ns("start_normalization"), "Start Normalization", class = "btn btn-success")
    ),
    hr(),
    uiOutput(ns("progress_ui")),
    hr(),
    uiOutput(ns("pca_plots_ui"))
  )
}

NormalizationClusteringServer <- function(id, filteredSeurat) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns

    normalizedSeurat <- reactiveVal(NULL)
    pcaPlots <- reactiveVal(list())
    progress_log <- reactiveVal("Waiting to start...")
    normalizationRunning <- reactiveVal(FALSE)
    progress_file <- reactiveVal(NULL)  # Keep track of temp file path

    output$progress_ui <- renderUI({
      if (normalizationRunning()) {
        div(style = "text-align: center; margin: 10px 0;",
            div(style = "font-weight: bold; margin-bottom: 5px;", progress_log()),
            div(class = "custom-spinner")
        )
      } else {
        NULL
      }
    })

    observeEvent(input$start_normalization, {
      req(filteredSeurat())
      # clear the current plots 
      pcaPlots(list())
      normalizationRunning(TRUE)
      shinyjs::disable("start_normalization")
      progress_log("Starting normalization...")

      # Rename the variable to avoid conflicts with built-in functions
      norm_method <- input$normalization_method
      obj <- filteredSeurat()

      temp_file <- tempfile()
      progress_file(temp_file)  # Save temp file path for later use

      future({
        logProgress <- function(message) {
          cat(message, "\n", file = temp_file, append = TRUE)
        }

        logProgress("Loading libraries...")
        library(Seurat)

        logProgress("Starting normalization...")
        if (norm_method == "standard") {
          logProgress("LogNormalizing data...")
          obj <- NormalizeData(obj)

          logProgress("Finding variable features...")
          obj <- FindVariableFeatures(obj, selection.method = "vst", nfeatures = 2000)

          logProgress("Scaling data...")
          obj <- ScaleData(obj)

          logProgress("Running PCA...")
          obj <- RunPCA(obj, features = VariableFeatures(obj))

          logProgress("Running JackStraw...")
          obj <- JackStraw(obj, num.replicate = 100)

          logProgress("Scoring JackStraw...")
          obj <- ScoreJackStraw(obj, dims = 1:20)
        } else if (norm_method == "sctransform") {
          logProgress("Running SCTransform...")
          obj <- SCTransform(obj, verbose = FALSE)

          logProgress("Running PCA...")
          obj <- RunPCA(obj, verbose = FALSE)
        }

        logProgress("Generating plots...")

        list(obj = obj, progress_file = temp_file, norm_method = norm_method)
      }, seed = TRUE) %...>% (function(result) {
        message("Result received:")
        str(result)

        if (!is.list(result)) {
          message("Error: result is not a list, check future output!")
          showNotification("Error: result is not a list!", type = "error")
          return(NULL)
        }

        normalizedSeurat(result$obj)

        plots <- list(
          elbow = ElbowPlot(result$obj) + ggtitle("Elbow Plot"),
          jackstraw = if (result$norm_method == "standard") {
            JackStrawPlot(result$obj, dims = 1:15) + ggtitle("JackStraw Plot")
          } else {
            NULL
          }
        )
        pcaPlots(plots)

        normalizationRunning(FALSE)
        shinyjs::enable("start_normalization")
      }) %...!% (function(err) {
        normalizationRunning(FALSE)
        shinyjs::enable("start_normalization")
        progress_log("Error during normalization.")
        message("error during normalization:", conditionMessage(err))
        showNotification(paste("Error:", conditionMessage(err)), type = "error")
      })

      # Poll the log file for updates
      observe({
        if (normalizationRunning() && !is.null(progress_file())) {
          if (file.exists(progress_file())) {
            log <- readLines(progress_file())
            if (length(log) > 0) {
              progress_log(tail(log, 1))  # Update UI with last log message
            }
          }
          invalidateLater(1000, session)
        }
      })
    })


    output$pca_plots_ui <- renderUI({
        req(pcaPlots())
      # lets just check to see if the length is greater than 0, then we will render
        if (length(pcaPlots()) > 0) {
            fluidRow(
                column(6, plotBox("elbow_plot", "download_elbow_plot", "Elbow Plot", ns,p("The elbow plot helps determine the optimal number of principal components (PCs) to retain for downstream clustering analysis by showing where the explained variance begins to level off. It is generally best to find the elbow and retain PCs up to that point. Although it is often acceptable to retain slighly more PC's then the Elbow point."))),
                column(6, if (!is.null(pcaPlots()$jackstraw)) {
                plotBox("jackstraw_plot", "download_jackstraw_plot", "JackStraw Plot", ns,"The JackStraw plot helps assess the statistical significance of each principal component (PC), allowing us to choose PCs that capture meaningful biological variation rather than noise. The plot shows the distribution of p-values for each PC, with significant PCs showing a deviation from the null hypothesis (p < 0.05).")
                } else {
                div("JackStraw Plot not available for SCTransform")
                })
            )
        }
    })

    output$elbow_plot <- renderPlot({
      req(pcaPlots())
      pcaPlots()$elbow
    })

    output$jackstraw_plot <- renderPlot({
      req(pcaPlots())
      pcaPlots()$jackstraw
    })

    output$download_elbow_plot <- downloadHandler(
      filename = function() { "ElbowPlot.png" },
      content = function(file) {
        ggsave(file, plot = pcaPlots()$elbow, width = 6, height = 4)
      }
    )

    output$download_jackstraw_plot <- downloadHandler(
      filename = function() { "JackStrawPlot.png" },
      content = function(file) {
        ggsave(file, plot = pcaPlots()$jackstraw, width = 6, height = 4)
      }
    )
  })
}
