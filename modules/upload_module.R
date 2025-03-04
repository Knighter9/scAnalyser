uploadUI <- function(id) {
    ns <- NS(id)
     tagList(
        fileInput(ns("seurat_file"), "Upload Seurat Object (.rds)"),
        textOutput(ns("placeholder_text")),
        textOutput(ns("error_message"))
    )
}

uploadServer <- function(id) {
    moduleServer(id, function(input, output, session) {
        error_message <- reactiveVal(NULL)

        observeEvent(input$seurat_file, {
            file_name <- input$seurat_file$name
            if (!grepl("\\.rds$", file_name, ignore.case = TRUE)) {
                error_message("Please upload a valid Seurat object in .rds format.")
            } else {
                error_message(NULL)
            }
        })

        output$error_message <- renderText({
            error_message()
        })

        reactive({
            req(input$seurat_file)

            if (!is.null(error_message())) {
                return(NULL)
            }

            readRDS(input$seurat_file$datapath)
        })
    })
}

