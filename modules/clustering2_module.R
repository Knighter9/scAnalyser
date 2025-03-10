Clustering2UI <- function(id) {
    ns <- NS(id)
    tagList(
        useShinyjs(),
        div(
            style = "text-align: center; margin-top: 30px;",
            h2("Final Clustering After Doublet Removal"),
            p("This step refines clustering analysis after filtering out doublets."),
            plotOutput(ns("final_clustering_ui")),
            downloadButton(ns("download_final_umap"), "Download UMAP Plot")
        )
    )
}


Clustering2Server <- function(id, doubletRemovalResults) {
    moduleServer(id, function(input, output, session) {
        ns <- session$ns

        finalSeurat <- reactive({
            req(doubletRemovalResults()) # Ensure the data exists
            doubletRemovalResults() # This returns the filtered (singlet-only) Seurat object
        })

        output$final_clustering_ui <- renderPlot({
            req(finalSeurat()) # Ensure Seurat object is available
            DimPlot(finalSeurat(), reduction = "umap") +
                ggtitle("Final Clustering After Doublet Removal")
        })

        output$download_final_umap <- downloadHandler(
            filename = function() {
                "Final_UMAP.png"
            },
            content = function(file) {
                ggsave(file, plot = DimPlot(finalSeurat(), reduction = "umap"), width = 6, height = 4)
            }
        )
    })
}
