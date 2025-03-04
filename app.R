library(shiny)
library(shinyjs)
library(Seurat)

source("modules/upload_module.R")
source("modules/QC_Analysis_Module.R")

ui <- fluidPage(
    useShinyjs(),
    uiOutput("main_ui"),
    tags$script(HTML("
        if (!sessionStorage.getItem('historyStack')) {
            sessionStorage.setItem('historyStack', JSON.stringify(['#upload']));
            console.log('Initialized historyStack:', sessionStorage.getItem('historyStack'));
        }
        if (!sessionStorage.getItem('forwardStack')) {
            sessionStorage.setItem('forwardStack', JSON.stringify([]));
        }

        // Ensure initial state is recorded (important for browser back/forward button to work correctly)
        history.replaceState({page: '#upload'}, '', '#upload');

        window.onpopstate = function(event) {
            console.log('Back/Forward button detected');

            // Just notify Shiny that the hash changed
            let newPage = window.location.hash || '#upload';
            console.log('Navigated to page (back/forward):', newPage);

            Shiny.setInputValue('browser_page', newPage, {priority: 'event'});
        };
    "))
)

server <- function(input, output, session) {
    current_page <- reactiveVal("upload")
    seuratData <- reactiveVal(NULL)

    # On initial page load (first visit or refresh)
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

    # Listen to browser back/forward events
    observeEvent(input$browser_page, {
        page <- substring(input$browser_page, 2)  # remove leading #
        message(sprintf("Browser page detected - setting page to: %s", page))
        current_page(page)
    })

    output$main_ui <- renderUI({
        if (current_page() == "upload") {
            uploadUI("uploader")
        } else if (current_page() == "QC_analysis") {
            QC_analysisUI("QC_analysis")
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
            sessionStorage.setItem('forwardStack', JSON.stringify([]));  // Clear forward on new navigation

            history.pushState({page: '#QC_analysis'}, '', '#QC_analysis');
        ")

        current_page("QC_analysis")
        updateQueryString("#QC_analysis", mode = "replace")
    })

    QC_analysisServer("QC_analysis", seuratData)
}

shinyApp(ui, server)
