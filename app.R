library(shiny)
library(shinyjs)
library(Seurat)
library(bslib)

source("modules/upload_module.R")
source("modules/QC_Analysis_Module.R")
source("modules/normalization_clustering_module.R") # New module
source("modules/clustering2_module.R") # New module
source("modules/cell_annotation_module.R") # New module

ui <- fluidPage(
    useShinyjs(),
    uiOutput("main_ui"),
    theme = bs_theme(version = 4),
    tags$script(HTML("
        if (!sessionStorage.getItem('historyStack')) {
            sessionStorage.setItem('historyStack', JSON.stringify(['#upload']));
            console.log('Initialized historyStack:', sessionStorage.getItem('historyStack'));
        }
        if (!sessionStorage.getItem('forwardStack')) {
            sessionStorage.setItem('forwardStack', JSON.stringify([]));
        }

        history.replaceState({page: '#upload'}, '', '#upload');

        window.onpopstate = function(event) {
            console.log('Back/Forward button detected');
            let newPage = window.location.hash || '#upload';
            Shiny.setInputValue('browser_page', newPage, {priority: 'event'});
        };
    "))
)

server <- function(input, output, session) {
    current_page <- reactiveVal("upload")
    seuratData <- reactiveVal(NULL)

    # Store filtered Seurat object from QC for handoff to Normalization
    filteredSeurat <- reactiveVal(NULL)
    # store doublet removal results for handoff to Clustering2
    doubletFilteredSeurat <- reactiveVal(NULL)
    # store final Seurat object for handoff to CellAnnotation
    finalSeurat <- reactiveVal(NULL)


    observe({
        hash <- session$clientData$url_hash
        message(sprintf("URL hash on load: %s", hash))

        if (is.null(hash) || hash == "") {
            updateQueryString("#upload", mode = "replace")
            current_page("upload")
        } else {
            page <- substring(hash, 2)
            message(sprintf("Hash found - setting page to: %s", page))
            current_page(page)
        }
    })

    observeEvent(input$browser_page, {
        page <- substring(input$browser_page, 2)
        message(sprintf("Browser page detected - setting page to: %s", page))
        current_page(page)
    })

    output$main_ui <- renderUI({
        if (current_page() == "upload") {
            uploadUI("uploader")
        } else if (current_page() == "QC_analysis") {
            QC_analysisUI("QC_analysis")
        } else if (current_page() == "normalization_clustering") {
            NormalizationClusteringUI("normalization_clustering")
        } else if (current_page() == "clustering2") {
            Clustering2UI("clustering2")
        } else if (current_page() == "cell_annotation") {
            CellAnnotationUI("cell_annotation")
        } else {
            div(class = "alert alert-danger", "Page not found.")
        }
    })

    uploadedSeurat <- uploadServer("uploader")

    observe({
        req(uploadedSeurat())
        seuratData(uploadedSeurat())

        showModal(modalDialog(
            title = "Upload Successful",
            p("Your Seurat object was uploaded successfully."),
            actionButton("proceed_btn", "Proceed With Analysis", class = "btn btn-success")
        ))
    })

    observeEvent(input$proceed_btn, {
        removeModal()
        runjs("
            let historyStack = JSON.parse(sessionStorage.getItem('historyStack')) || ['#upload'];
            historyStack.push('#QC_analysis');
            sessionStorage.setItem('historyStack', JSON.stringify(historyStack));
            sessionStorage.setItem('forwardStack', JSON.stringify([]));
            history.pushState({page: '#QC_analysis'}, '', '#QC_analysis');
        ")
        current_page("QC_analysis")
        updateQueryString("#QC_analysis", mode = "replace")
    })

    # Hand off filtered Seurat object from QC to Normalization page
    QC_analysisServer("QC_analysis", seuratData, onProceed = function(filteredObj) {
        filteredSeurat(filteredObj)

        runjs("
            let historyStack = JSON.parse(sessionStorage.getItem('historyStack')) || ['#upload'];
            historyStack.push('#normalization_clustering');
            sessionStorage.setItem('historyStack', JSON.stringify(historyStack));
            sessionStorage.setItem('forwardStack', JSON.stringify([]));
            history.pushState({page: '#normalization_clustering'}, '', '#normalization_clustering');
        ")
        current_page("normalization_clustering")
        updateQueryString("#normalization_clustering", mode = "replace")
    })

    NormalizationClusteringServer("normalization_clustering", filteredSeurat)

    NormalizationClusteringServer("normalization_clustering", filteredSeurat, onProceed = function(singletObj) {
        doubletFilteredSeurat(singletObj) # Store the final Seurat object

        runjs("
            let historyStack = JSON.parse(sessionStorage.getItem('historyStack')) || ['#upload'];
            historyStack.push('#clustering2');
            sessionStorage.setItem('historyStack', JSON.stringify(historyStack));
            sessionStorage.setItem('forwardStack', JSON.stringify([]));
            history.pushState({page: '#clustering2'}, '', '#clustering2');
        ")
        current_page("clustering2")
        updateQueryString("#clustering2", mode = "replace")
    })

    #  Pass the final filtered Seurat object to Clustering2Server
    Clustering2Server("clustering2", doubletFilteredSeurat)

    # handle on proceeding object for clustering2Server to launch cell annotation page
    Clustering2Server("clustering2", doubletFilteredSeurat, onProceed = function(finalObj) {
        finalSeurat(finalObj) # Store the final Seurat object
        runjs("
            let historyStack = JSON.parse(sessionStorage.getItem('historyStack')) || ['#upload'];
            historyStack.push('#cell_annotation');
            sessionStorage.setItem('historyStack', JSON.stringify(historyStack));
            sessionStorage.setItem('forwardStack', JSON.stringify([]));
            history.pushState({page: '#cell_annotation'}, '', '#cell_annotation');
        ")
        current_page("cell_annotation")
        updateQueryString("#cell_annotation", mode = "replace")
    })
    # Pass the final Seurat object to CellAnnotationServer
    CellAnnotationServer("cell_annotation", finalSeurat)
}

shinyApp(ui, server)
