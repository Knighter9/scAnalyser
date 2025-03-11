library(shiny)
library(shinyjs)
library(Seurat)

# UI function for the Cell Annotation module
CellAnnotationUI <- function(id) {
    ns <- NS(id)
    fluidPage(
        useShinyjs(),
        h2("Cell Annotation Module"),
        fluidRow(
            column(
                width = 8,
                h3("UMAP Plot of Clustered Cells"),
                plotOutput(ns("umap_plot")),
                h3("UMAP Plot with Assigned Cell Identities"),
                uiOutput(ns("umap_plots_ui")), # Dynamic UI for UMAPs
                h3("Generate Gene Expression Plots"),
                textInput(ns("gene_input"), "Enter Gene:", placeholder = "e.g., CD3E"),
                actionButton(ns("plot_gene"), "Generate Plots"),
                br(), br(),
                plotOutput(ns("feature_plot")),
                br(),
                plotOutput(ns("violin_plot"))
            ),
            column(
                width = 4,
                h3("Cluster Annotation Table"),
                uiOutput(ns("cluster_table_ui")),
                br(),
                actionButton(ns("save_cluster_annotations"), "Save Cluster Annotations")
            )
        ),
        hr(),
        verbatimTextOutput(ns("annotation_info"))
    )
}

# Server function for the Cell Annotation module
CellAnnotationServer <- function(id, clustered_obj) {
    moduleServer(id, function(input, output, session) {
        ns <- session$ns
        assigned_obj <- reactiveVal(NULL) # Store updated object with cell identities

        # Render original UMAP plot
        output$umap_plot <- renderPlot({
            req(clustered_obj())
            DimPlot(clustered_obj(), reduction = "umap", label = TRUE) +
                ggtitle("UMAP Plot - Clustered Cells")
        })

        # Render the UI for UMAP plots dynamically
        output$umap_plots_ui <- renderUI({
            req(assigned_obj()) # Only render after saving annotations
            tagList(
                h3("Broad Cell Type UMAP"),
                plotOutput(ns("umap_broad")),
                h3("Subtype Cell Type UMAP"),
                plotOutput(ns("umap_subtype"))
            )
        })

        # Generate FeaturePlot and VlnPlot for the entered gene
        observeEvent(input$plot_gene, {
            req(clustered_obj(), input$gene_input)
            gene <- input$gene_input
            if (gene %in% rownames(clustered_obj())) {
                output$feature_plot <- renderPlot({
                    FeaturePlot(clustered_obj(), features = gene) +
                        ggtitle(paste("Feature Plot for", gene))
                })
                output$violin_plot <- renderPlot({
                    VlnPlot(clustered_obj(), features = gene) +
                        ggtitle(paste("Violin Plot for", gene))
                })
            } else {
                showNotification(paste("Gene", gene, "not found in the dataset."), type = "error")
            }
        })

        # Build a dynamic UI table for cluster annotation
        output$cluster_table_ui <- renderUI({
            req(clustered_obj())
            clusters <- sort(unique(clustered_obj()@meta.data$seurat_clusters))
            table_rows <- lapply(clusters, function(cl) {
                fluidRow(
                    column(4, strong(paste("Cluster", cl))),
                    column(4, textInput(ns(paste0("broad_", cl)), "Broad Cell Type", value = "")),
                    column(4, textInput(ns(paste0("subtype_", cl)), "Subtype", value = ""))
                )
            })
            tagList(table_rows)
        })

        # When the save button is pressed, update the object with new identities
        observeEvent(input$save_cluster_annotations, {
            req(clustered_obj())
            clusters <- sort(unique(clustered_obj()@meta.data$seurat_clusters))

            # Create named vectors to store cell type identities
            broad_idents <- setNames(rep(NA, length(clusters)), as.character(clusters))
            subtype_idents <- setNames(rep(NA, length(clusters)), as.character(clusters))

            # Fill in the new identities from user input
            for (cl in clusters) {
                broad_idents[as.character(cl)] <- input[[paste0("broad_", cl)]]
                subtype_idents[as.character(cl)] <- input[[paste0("subtype_", cl)]]
            }

            # Create a modified Seurat object with updated metadata
            new_obj <- clustered_obj()
            new_obj@meta.data$broad_cell_type <- factor(new_obj@meta.data$seurat_clusters,
                levels = names(broad_idents),
                labels = broad_idents
            )
            new_obj@meta.data$subtype_cell_type <- factor(new_obj@meta.data$seurat_clusters,
                levels = names(subtype_idents),
                labels = subtype_idents
            )

            print("Updated Metadata (Head):")
            print(head(new_obj@meta.data))

            # Store updated object
            assigned_obj(new_obj)

            showNotification("Cluster annotations saved! UMAPs with cell identities are now visible.", type = "message")
        })

        # Render UMAP plot with Broad Cell Types
        output$umap_broad <- renderPlot({
            req(assigned_obj())
            DimPlot(assigned_obj(), reduction = "umap", group.by = "broad_cell_type", label = TRUE) +
                ggtitle("UMAP Plot - Broad Cell Types")
        })

        # Render UMAP plot with Subtype Cell Types
        output$umap_subtype <- renderPlot({
            req(assigned_obj())
            DimPlot(assigned_obj(), reduction = "umap", group.by = "subtype_cell_type", label = TRUE) +
                ggtitle("UMAP Plot - Subtype Cell Types")
        })

        # Additional info for debugging / summary
        output$annotation_info <- renderPrint({
            req(clustered_obj())
            meta <- clustered_obj()@meta.data
            cat("Clustered Seurat Object Summary:\n")
            cat("Number of cells:", nrow(meta), "\n")
            if ("seurat_clusters" %in% colnames(meta)) {
                cat("Number of clusters:", length(unique(meta$seurat_clusters)), "\n")
            } else {
                cat("No clustering information found in metadata.\n")
            }
        })
    })
}
