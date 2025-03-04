library(shiny)
library(Seurat)

# Source modules
source("modules/upload_module.R")
source("modules/QC_Analysis_Module.R")

ui <- fluidPage(
    uiOutput("main_ui")  # Dynamically switch pages
)

server <- function(input, output, session) {

    # This tracks which page we're on
    current_page <- reactiveVal("upload")

    # Global object to store Seurat object
    seuratData <- reactiveVal(NULL)

    # --- URL Hash Handling (default to 'upload' if no hash) ---
    observe({
        hash <- session$clientData$url_hash
        if (is.null(hash) || hash == "") {
            # No hash = default to upload page
            current_page("upload")
            updateQueryString("#upload", mode = "replace")
        } else {
            # Remove leading # and update current_page()
            current_page(substring(hash, 2))
            print("Current page set to: ")
            print(current_page())
        }
    })

    # --- Dynamically render page based on current_page() ---
    output$main_ui <- renderUI({
        if (current_page() == "upload") {
            uploadUI("uploader")
        } else if (current_page() == "QC_analysis") {
            QC_analysisUI("QC_analysis")
        }
    })

    # --- Handle Seurat upload ---
    uploadedSeurat <- uploadServer("uploader")

    # Store uploaded object into global reactiveVal
    observe({
        req(uploadedSeurat())
        seuratData(uploadedSeurat())
    })

    # --- Show "Proceed" button after successful upload ---
    observe({
        if (!is.null(seuratData())) {
            showModal(modalDialog(
                title = "Upload Successful",
                p("Your Seurat object was uploaded successfully."),
                actionButton("proceed_btn", "Proceed With Analysis", class = "btn btn-success")
            ))
        }
    })

    # --- Handle "Proceed" button click to move to QC Analysis ---
    observeEvent(input$proceed_btn, {
        removeModal()
        current_page("QC_analysis")
        updateQueryString("#QC_analysis", mode = "replace")  # Update URL
    })

    # --- Handle Back/Forward button via polling ---
    observe({
        invalidateLater(500, session)
        hash <- session$clientData$url_hash
        new_page <- ifelse(is.null(hash) || hash == "", "upload", substring(hash, 2))
        if (new_page != current_page()) {
            current_page(new_page)
        }
    })

    # --- Initialize the QC Analysis module ---
    QC_analysisServer("QC_analysis", seuratData)
}

shinyApp(ui, server)
