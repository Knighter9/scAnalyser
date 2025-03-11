library(clustree)
source("modules/normalization_clustering_module.R")

Clustering2UI <- function(id) {
    ns <- NS(id)
    tagList(
        useShinyjs(),

        # Standard Clustering Tree
        div(
            style = "text-align: center; margin-top: 30px;",
            h2("Clustering Tree Visualization"),
            p("Explore how clusters change across different resolutions."),
            plotOutput(ns("clustree_plot")),
            downloadButton(ns("download_clustree_plot"), "Download Clustering Tree")
        ),

        # Gene Expression Clustering Tree
        div(
            style = "text-align: center; margin-top: 30px;",
            h2("Clustering Tree Colored by Gene Expression"),
            p("Enter a gene to visualize its expression in the clustering tree."),
            textInput(ns("gene_input"), "Enter Gene of Interest:", placeholder = "E.g., CD3E"),
            actionButton(ns("update_gene_clustree"), "Update Gene Clustree"),
            plotOutput(ns("gene_clustree_plot")),
            downloadButton(ns("download_gene_clustree_plot"), "Download Gene Expression Clustering Tree")
        ),

        # Final UMAP plot after doublet removal (wrapped in plotBox)
        div(
            style = "text-align: center; margin-top: 30px;",
            h2("Final Clustering After Doublet Removal"),
            p("This step refines clustering analysis after filtering out doublets."),
            plotBox(
                plotOutputId = "final_clustering_ui",
                downloadButtonId = "download_final_umap",
                plotTitle = "Final Clustering After Doublet Removal",
                ns = ns,
                height = "600px",
                width = "600px",
                passedPlotDim = TRUE
            )
        ),

        # Clustering Resolution Selection & Proceed Button
        div(
            style = "text-align: center; margin-top: 30px;",
            h2("Proceed to Cell Annotation"),
            p("Select a clustering resolution to proceed with cell annotation."),
            selectInput(ns("resolution_select"), "Select Clustering Resolution:", choices = seq(0, 1.0, by = 0.1), selected = 0.8),
            actionButton(ns("proceed_to_annotation"), "Proceed to Cell Annotation", class = "btn btn-primary")
        )
    )
}

Clustering2Server <- function(id, doubletRemovalResults, num_pcs = 10, onProceed = NULL) {
    moduleServer(id, function(input, output, session) {
        ns <- session$ns

        finalSeurat <- reactive({
            req(doubletRemovalResults())
            doubletRemovalResults() # Get the filtered Seurat object
        })

        # Rerun UMAP and clustering at multiple resolutions without future-based parallelism
        finalSeuratClustered <- reactive({
            req(finalSeurat())
            seurat_obj <- finalSeurat()
            future::plan(sequential) # Set future plan to sequential for reproducibility
            # Rerun UMAP because cells were removed
            seurat_obj <- FindNeighbors(seurat_obj, dims = 1:num_pcs)

            # Manually loop over clustering resolutions (avoid resolution 0)
            resolutions <- seq(0.0, 1.0, by = 0.1)
            for (res in resolutions) {
                set.seed(123) # Set a fixed seed for reproducibility (optional)
                seurat_obj <- FindClusters(seurat_obj, resolution = res)
            }

            seurat_obj <- RunUMAP(seurat_obj, dims = 1:num_pcs)
            return(seurat_obj)
        })

        # Render updated UMAP plot
        output$final_clustering_ui <- renderPlot({
            req(finalSeuratClustered())
            DimPlot(finalSeuratClustered(), reduction = "umap") +
                ggtitle("Final Clustering After Doublet Removal")
        })

        # Download updated UMAP plot
        output$download_final_umap <- downloadHandler(
            filename = function() {
                "Final_UMAP.png"
            },
            content = function(file) {
                ggsave(file, plot = DimPlot(finalSeuratClustered(), reduction = "umap"), width = 6, height = 6)
            }
        )

        # Standard clustering tree visualization
        output$clustree_plot <- renderPlot({
            req(finalSeuratClustered())
            clustree(finalSeuratClustered(), prefix = "SCT_snn_res.")
        })

        # Download clustering tree plot
        output$download_clustree_plot <- downloadHandler(
            filename = function() {
                "Clustering_Tree.png"
            },
            content = function(file) {
                plot <- clustree(finalSeuratClustered(), prefix = "SCT_snn_res.")
                ggsave(file, plot = plot, width = 8, height = 6)
            }
        )

        # Reactive value for user-selected gene
        selectedGene <- reactiveVal(NULL)

        observeEvent(input$update_gene_clustree, {
            gene <- input$gene_input

            # Check if the gene exists in the Seurat object
            if (!is.null(gene) && gene %in% rownames(finalSeuratClustered())) {
                selectedGene(gene)
            } else {
                showNotification("Gene not found in dataset. Showing default clustering tree.", type = "error")
                selectedGene(NULL) # Reset to default
            }
        })

        # Render clustering tree with gene expression (if valid gene is entered)
        output$gene_clustree_plot <- renderPlot({
            req(finalSeuratClustered())

            if (!is.null(selectedGene())) {
                clustree(finalSeuratClustered(),
                    prefix = "SCT_snn_res.",
                    node_colour = selectedGene(), node_colour_aggr = "median"
                )
            } else {
                clustree(finalSeuratClustered(), prefix = "SCT_snn_res.")
            }
        })

        # Download gene expression clustering tree plot
        output$download_gene_clustree_plot <- downloadHandler(
            filename = function() {
                "Gene_Expression_Clustering_Tree.png"
            },
            content = function(file) {
                if (!is.null(selectedGene())) {
                    plot <- clustree(finalSeuratClustered(),
                        prefix = "SCT_snn_res.",
                        node_colour = selectedGene(), node_colour_aggr = "median"
                    )
                } else {
                    plot <- clustree(finalSeuratClustered(), prefix = "SCT_snn_res.")
                }
                ggsave(file, plot = plot, width = 8, height = 6)
            }
        )

        # Observe event for "Proceed to Annotation" button
        observeEvent(input$proceed_to_annotation, {
            req(finalSeuratClustered())

            # Get selected clustering resolution
            selected_resolution <- input$resolution_select
            message(paste("Proceeding with resolution:", selected_resolution))

            # Apply selected resolution to the Seurat object (without future parameters)
            clustered_obj <- finalSeuratClustered()
            clustered_obj <- FindClusters(clustered_obj, resolution = as.numeric(selected_resolution))

            # Proceed with next step (cell annotation) if function is provided
            if (!is.null(onProceed)) {
                onProceed(clustered_obj)
            }

            showNotification(paste("Proceeding with clustering resolution:", selected_resolution), type = "message")
        })
    })
}
