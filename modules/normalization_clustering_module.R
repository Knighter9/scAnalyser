library(Seurat)
library(ggplot2)
library(shiny)
library(shinyjs)
library(future)
library(promises)
library(DoubletFinder)
library(dplyr)

# Set up multisession plan for futures
plan(multisession)
# increase future.globals.maxsize to 1g
options(future.globals.maxSize = 1024 * 1024 * 1024 * 10)

# mulitplet rates for doublet identification
multiplet_rates <- data.frame(
    Cells_Loaded = c(
        825, 1650, 3300, 4950, 6600, 8250, 9900, 11550,
        13200, 14850, 16500
    ),
    Cells_Recovered = c(
        500, 1000, 2000, 3000, 4000, 5000, 6000,
        7000, 8000, 9000, 10000
    ),
    Multiplet_Rate = c(
        0.004, 0.008, 0.016, 0.024, 0.032, 0.04,
        0.048, 0.056, 0.064, 0.072, 0.08
    )
)
getMultipletRate <- function(num_cells_recovered) {
    closest_row <- multiplet_rates[which.min(abs(multiplet_rates$Cells_Recovered - num_cells_recovered)), ]
    return(closest_row$Multiplet_Rate)
}
# A helper function to create a box with a plot and download button
plotBox <- function(plotOutputId, downloadButtonId, plotTitle, ns, explanationText = NULL, height = "600px", width = "600px", passedPlotDim = FALSE) {
    div(
        style = "border: 1px solid #ccc; padding: 10px; margin: 5px;",
        h4(plotTitle),
        # check to see if they passed plot dimensions
        if (passedPlotDim) {
            print("plot dimensions passed")
            plotOutput(ns(plotOutputId), height = height, width = width)
        } else {
            plotOutput(ns(plotOutputId))
        },
        downloadButton(ns(downloadButtonId), paste0("Download ", plotTitle, " Plot")),
        if (!is.null(explanationText)) {
            p(style = "font-size: 0.8em; margin-top: 10px;", explanationText)
        }
    )
}

# pc selection function
pcSelectionUI <- function(normalizedSeurat, numSelectedPCs, ns) {
    req(normalizedSeurat())
    tagList(
        h4("Select Number of Principal Components (PCs) for Clustering"),
        sliderInput(
            ns("num_pcs"),
            "Number of PCs to Use:",
            min = 5,
            max = 50,
            value = numSelectedPCs(),
            step = 1
        ),
        # add in an action button to proceed with clustering and downstream analysis with the selected pcs
        actionButton(ns("proceed_clustering"), "Proceed with Clustering", class = "btn btn-success"), # an event will handle a click of this button so we can respond accordingly
    )
}

# run doublet finder function
runDoubletFinder <- function(obj, pcs, num_cells_recovered, norm_method) {
    message("Running DoubletFinder process...")
    message("about to print obj")
    message("************")
    message(head(obj))
    message("------------")
    print("about to print num cells recovered")
    print(num_cells_recovered)
    # Get multiplet rate from 10x table
    multiplet_rate <- getMultipletRate(num_cells_recovered)
    print("about to print the mulitplet rate")
    print(multiplet_rate)
    # Base expected doublets
    num_cells <- nrow(obj@meta.data)
    print("The number of cells is")
    print(num_cells)

    nExp <- round(num_cells * multiplet_rate)
    print("nExp")
    print(nExp)
    # Use Seurat clusters as "cell types" for homotypic adjustment
    annotations <- obj@meta.data$seurat_clusters # Directly use Seurat's clustering
    homotypic.prop <- DoubletFinder::modelHomotypic(annotations)
    nExp_adj <- round(nExp * (1 - homotypic.prop))

    print("nExp_adj")
    print(nExp_adj)

    print("about to print obj uisng print")
    print(head(obj))
    message(paste("Homotypic proportion:", homotypic.prop))
    message(paste("Adjusted expected doublets:", nExp_adj))


    # Param sweep and optimal pK detection
    sweep.res.list <- NULL
    if (norm_method == "standard") {
        sweep.res.list <- DoubletFinder::paramSweep(
            obj,
            PCs = 1:pcs,
            sct = FALSE
        )
    } else if (norm_method == "sctransform") {
        sweep.res.list <- DoubletFinder::paramSweep(
            obj,
            PCs = 1:pcs,
            sct = TRUE
        )
    }
    sweep.stats <- DoubletFinder::summarizeSweep(sweep.res.list, GT = FALSE)

    bcmvn <- DoubletFinder::find.pK(sweep.stats)

    # Optimal pK
    optimal_pk <- as.numeric(as.character(bcmvn$pK[which.max(bcmvn$BCmetric)]))

    print("optimal_pk")
    print(optimal_pk)

    # Run DoubletFinder with adjusted nExp
    if (norm_method == "standard") {
        obj <- DoubletFinder::doubletFinder(
            obj,
            PCs = 1:pcs,
            pN = 0.25,
            pK = optimal_pk,
            nExp = nExp_adj, # Use adjusted value
            reuse.pANN = FALSE,
            sct = FALSE
        )
    } else if (norm_method == "sctransform") {
        print("running doublet finder with sctransform")
        obj <- DoubletFinder::doubletFinder(
            obj,
            PCs = 1:pcs,
            pN = 0.25,
            pK = optimal_pk,
            nExp = nExp_adj, # Use adjusted value
            reuse.pANN = FALSE,
            sct = TRUE
        )
    }

    # Add doublet status to metadata
    doublet_column <- grep("DF.classifications", colnames(obj@meta.data), value = TRUE)
    obj$doublet_status <- obj@meta.data[[doublet_column]]
    print("about to print the object one last time")
    print(head(obj@meta.data))

    return(obj)
}

# normalization and clustering UI
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
        div(
            style = "text-align: center; margin-top: 30px;",
            h2("Normalization and Clustering"),
            radioButtons(ns("normalization_method"), "Normalization Method:",
                choices = c(
                    "Standard Seurat Workflow" = "standard",
                    "SCTransform" = "sctransform"
                ),
                inline = TRUE
            ),
            br(),
            actionButton(ns("start_normalization"), "Start Normalization", class = "btn btn-success")
        ),
        hr(),
        uiOutput(ns("progress_ui")),
        hr(),
        uiOutput(ns("pca_plots_ui")),
        hr(),
        uiOutput(ns("pc_selection_ui")),
        hr(),
        uiOutput(ns("clustering_results_ui")),
        hr(),
        # add in ui output for the clustering results here
        hr(),
        # add in the doublet removal controls
        uiOutput(ns("doublet_controls_ui")),
        uiOutput(ns("doublet_umap_ui"))
    )
}


NormalizationClusteringServer <- function(id, filteredSeurat) {
    moduleServer(id, function(input, output, session) {
        # ************************************
        # set up the necessary reactive vals here
        # ************************************
        ns <- session$ns
        normalizedSeurat <- reactiveVal(NULL)
        pcaPlots <- reactiveVal(list())
        progress_log <- reactiveVal("Waiting to start...")
        normalizationRunning <- reactiveVal(FALSE)
        clusteringRunning <- reactiveVal(FALSE)
        doubletRemovalRunning <- reactiveVal(FALSE)
        progress_file <- reactiveVal(NULL) # Keep track of temp file path
        selectedPCs <- reactiveVal(FALSE) # set the selecting pcs reactiveVal to be a false boolean, we will use this to determine when to show the
        numSelectedPCs <- reactiveVal(10) # set the default num selected pcs to be 10
        clusteringResults <- reactiveVal(NULL) # this will be used to store the seurat object after intial clustering has been run
        doubletRemovalResults <- reactiveVal(NULL) # this will be used to store the seurat object after doublet removal function has been run

        # ************************************
        # observe events here
        # ***********************************
        observeEvent(input$start_normalization, {
            # the normalizatoin button has been pressed so we are going to start the process
            # first step is getting the seurat object, clearing the pcaPlots reactive val in case they are re running normalization
            req(filteredSeurat())
            pcaPlots(list())
            normalizationRunning(TRUE)
            shinyjs::disable("start_normalization")
            progress_log("Starting normalization...")

            norm_method <- input$normalization_method
            obj <- filteredSeurat() # get the filtered seurat object

            temp_file <- tempfile()
            progress_file(temp_file) # Save temp file path for later use

            future(
                {
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
                },
                seed = TRUE
            ) %...>% (function(result) {
                # handle the return from a successful future, we will get a normalized seurat object, the normalized object will be in result so we can updat the correct reactive val.
                message("Result received:")
                str(result)

                if (!is.list(result)) {
                    message("Error: result is not a list, check future output!")
                    showNotification("Error: result is not a list!", type = "error")
                    return(NULL)
                }
                normalizedSeurat(result$obj) # updating the reactive val with the normalized seurat object
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
                # set the selectedPCs reactiveVal to be true, this will trigger the UI to show the select the PCs for clustering UI
                selectedPCs(TRUE)
            }) %...!% (function(err) { # error handling
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
                            progress_log(tail(log, 1)) # Update UI with last log message
                        }
                    }
                    invalidateLater(1000, session)
                }
            })
        })
        # add another observerve event to handle the proceed with clustering buttoin after it has been clicked
        # after the button has been clicked we will also disable the rerun normalization button and then
        # we can just for testing purposed right now run clustering and then show a umap so we can see that it is working
        observeEvent(input$proceed_clustering, {
            message("Proceed with clustering button clicked")
            shinyjs::disable("start_normalization")
            shinyjs::disable("proceed_clustering")
            clusteringRunning(TRUE)
            # we can run clustering here and then show a umap
            # run the umap with the chosen pcs
            obj <- normalizedSeurat()
            num_pcs <- numSelectedPCs()
            temp_file <- tempfile()
            progress_file(temp_file) # Save temp file path for later use
            progress_log("Running clustering...")
            # lets run this process in a future so that we can show a spinner while the process is running
            future(
                {
                    logProgress <- function(message) {
                        cat(message, "\n", file = temp_file, append = TRUE)
                    }

                    logProgress("Running clustering...")
                    obj <- FindNeighbors(obj, dims = 1:num_pcs)
                    obj <- FindClusters(obj)
                    obj <- RunUMAP(obj, dims = 1:num_pcs)

                    # return the object make sure seed is set to true?
                    obj
                },
                seed = TRUE
            ) %...>% (function(result) {
                message("Result received:")
                # populate the clustering results reactiveVal with the result
                clusteringResults(result)
                clusteringRunning(FALSE)
                # enable the start normalization button
                shinyjs::enable("start_normalization")
                shinyjs::enable("proceed_clustering")
                # show the doublet removal controls
                shinyjs::show("doublet_controls")
            }) %...!% (function(err) {
                message("Error during clustering:", conditionMessage(err))
                showNotification(paste("Error:", conditionMessage(err)), type = "error")
            })
            # poll the log file for updates
            observe({
                if (clusteringRunning() && !is.null(progress_file())) {
                    if (file.exists(progress_file())) {
                        log <- readLines(progress_file())
                        if (length(log) > 0) {
                            progress_log(tail(log, 1)) # Update UI with last log message
                        }
                    }
                    invalidateLater(1000, session)
                }
            })
        })

        observeEvent(input$run_doublet_removal, {
            message("run doublet removal buttoin clicked")
            # disable the button and noramlization and clustering buttons
            shinyjs::disable("run_doublet_removal")
            shinyjs::disable("start_normalization")
            shinyjs::disable("proceed_clustering")
            doubletRemovalRunning(TRUE)
            progress_log("Running doublet removal...")
            # get the clustering results
            req(clusteringResults())
            obj <- clusteringResults()
            norm_method <- input$normalization_method
            message(head(obj))
            message("the object is of type")
            message(class(obj))
            print("about to check for seurat_clusters")
            print(head(obj))
            if (!("seurat_clusters" %in% colnames(obj@meta.data))) {
                stop("No seurat_clusters found in object metadata. Please make sure clustering step was successful.")
            }
            num_pcs <- input$num_pcs
            num_cells_recovered <- input$cells_recovered
            temp_file <- tempfile()
            progress_file(temp_file) # Save temp file path for later use
            doublet_rate <- getMultipletRate(num_cells_recovered)
            output$doublet_rate_display <- renderText(paste("Estimated Multiplet Rate:", round(doublet_rate * 100, 2), "%"))
            future(
                {
                    logProgress <- function(message) {
                        cat(message, "\n", file = temp_file, append = TRUE)
                    }

                    logProgress("Doublet Finder Running...")
                    runDoubletFinder(obj, num_pcs, num_cells_recovered, norm_method)
                },
                seed = TRUE
            ) %...>% (function(result) {
                message("Doublet removal result received:")
                doubletRemovalResults(result)
                doubletRemovalRunning(FALSE)
                shinyjs::enable("run_doublet_removal")
                # save the rds file
                saveRDS(result, file = "doublet_removal_results.rds") # tempory so i can bypass normalization and clustring
            }) %...!% (function(err) {
                message("Error during doublet removal:", conditionMessage(err))
                showNotification(paste("Error:", conditionMessage(err)), type = "error")
                shinyjs::enable("run_doublet_removal")
            })

            # poll the log file for updates
            observe({
                if (doubletRemovalRunning() && !is.null(progress_file())) {
                    if (file.exists(progress_file())) {
                        log <- readLines(progress_file())
                        if (length(log) > 0) {
                            progress_log(tail(log, 1)) # Update UI with last log message
                        }
                    }
                    invalidateLater(1000, session)
                }
            })
        })
        # ************************************
        # render ui and plot functions here
        # ************************************
        # Setting up the spinner for normilzation, we will have it run while the normalization is running, otherwise it will not run.
        output$progress_ui <- renderUI({
            if (normalizationRunning() || clusteringRunning() || doubletRemovalRunning()) {
                div(
                    style = "text-align: center; margin: 10px 0;",
                    div(style = "font-weight: bold; margin-bottom: 5px;", progress_log()),
                    div(class = "custom-spinner")
                )
            } else {
                NULL
            }
        })
        output$pca_plots_ui <- renderUI({
            req(pcaPlots())
            # lets just check to see if the length is greater than 0, then we will render
            if (length(pcaPlots()) > 0) {
                fluidRow(
                    column(6, plotBox("elbow_plot", "download_elbow_plot", "Elbow Plot", ns, p("The elbow plot helps determine the optimal number of principal components (PCs) to retain for downstream clustering analysis by showing where the explained variance begins to level off. It is generally best to find the elbow and retain PCs up to that point. Although it is often acceptable to retain slighly more PC's then the Elbow point."))),
                    column(6, if (!is.null(pcaPlots()$jackstraw)) {
                        plotBox("jackstraw_plot", "download_jackstraw_plot", "JackStraw Plot", ns, "The JackStraw plot helps assess the statistical significance of each principal component (PC), allowing us to choose PCs that capture meaningful biological variation rather than noise. The plot shows the distribution of p-values for each PC, with significant PCs showing a deviation from the null hypothesis (p < 0.05).")
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


        # set up the pc selection UI
        output$pc_selection_ui <- renderUI({
            req(selectedPCs())
            if (selectedPCs()) {
                pcSelectionUI(normalizedSeurat, numSelectedPCs, ns)
            }
        })
        output$umap_plot <- renderPlot({
            req(clusteringResults())
            DimPlot(clusteringResults(), reduction = "umap")
        })

        output$clustering_results_ui <- renderUI({
            req(clusteringResults())
            # we create the dim plot

            # lets just check to see if the length is greater than 0, then we will render
            if (length(clusteringResults()) > 0) {
                # use the plotBox function to create a box with the umap plot and a download button
                fluidRow(
                    column(12, plotBox("umap_plot", "download_umap_plot", "UMAP Plot", ns, explanationText = NULL, height = "600px", width = "600px", TRUE))
                )
            }
        })
        output$doublet_umap <- renderPlot({
            req(doubletRemovalResults())
            if (length(doubletRemovalResults()) > 0) {
                DimPlot(doubletRemovalResults(), reduction = "umap", group.by = "doublet_status") +
                    ggtitle("UMAP After Doublet Removal") +
                    theme(plot.title = element_text(hjust = 0.5))
            }
        })
        output$doublet_umap_ui <- renderUI({
            req(doubletRemovalResults())
            if (length(doubletRemovalResults()) > 0) {
                # use the plotBox function to create a box with the umap plot and a download button
                # we can also add in an explanation text here
                # we will use the plotBox function to create a box with the umap plot and a download button
                # we can also add in an explanation text here
                obj <- doubletRemovalResults()
                print("about to print doublet plot")
                print(head(obj))
                print("the object is of type")
                print(class(obj))
                print("About to print the meta data")
                print(head(obj@meta.data))
                fluidRow(column(12, plotBox("doublet_umap", "download_doublet_umap", "UMAP After Doublet Removal", ns, height = "600px", width = "600px", passedPlotDim = TRUE)))
            }
        })
        # set up the doublet removal controls
        output$doublet_controls_ui <- renderUI({
            # tag list of the necessary elements
            req(clusteringResults())
            if (length(clusteringResults()) > 0) {
                tagList(
                    h4("Doublet Removal"),
                    numericInput(ns("cells_recovered"), "Number of Cells Recovered:", value = 10000, min = 1000, max = 20000, step = 1000),
                    actionButton(ns("run_doublet_removal"), "Run DoubletFinder", class = "btn btn-success"),
                    textOutput("doublet_rate_display")
                )
            }
        })
        # *************************************
        # download handlers here
        # *************************************
        output$download_umap_plot <- downloadHandler(
            filename = function() {
                "UMAPPlot.png"
            },
            content = function(file) {
                ggsave(file, plot = DimPlot(clusteringResults(), reduction = "umap"), width = 6, height = 4) # nolint: indentation_linter.
            }
        )
        output$download_elbow_plot <- downloadHandler(
            filename = function() {
                "ElbowPlot.png"
            },
            content = function(file) {
                ggsave(file, plot = pcaPlots()$elbow, width = 6, height = 4)
            }
        )
        output$download_jackstraw_plot <- downloadHandler(
            filename = function() {
                "JackStrawPlot.png"
            },
            content = function(file) {
                ggsave(file, plot = pcaPlots()$jackstraw, width = 6, height = 4)
            }
        )
    })
}
