library(Seurat)
library(ggplot2)

# Calculate QC Metrics (nFeature, nCount, percent.mt)
calculateQCMetrics <- function(seuratObj) {
    counts <- GetAssayData(seuratObj, layer = "counts")

    seuratObj$nFeature_RNA <- Matrix::colSums(counts > 0)
    seuratObj$nCount_RNA <- Matrix::colSums(counts)

    mito_genes <- grep("^MT-", rownames(seuratObj), value = TRUE)
    percent_mt <- Matrix::colSums(counts[mito_genes, , drop = FALSE]) / seuratObj$nCount_RNA * 100
    seuratObj$percent.mt <- percent_mt

    return(seuratObj)
}

# Reusable plot container function (for pre and post plots)
plotBox <- function(plotOutputId, downloadButtonId, plotTitle, ns) {
    div(
        style = "border: 1px solid #ccc; padding: 10px; margin: 5px; text-align: center;",
        h4(plotTitle),
        plotOutput(ns(plotOutputId)),
        downloadButton(ns(downloadButtonId), paste0("Download ", plotTitle, " Plot"))
    )
}

# QC Analysis UI
QC_analysisUI <- function(id) {
    ns <- NS(id)

    tagList(
        h3("QC Analysis"),
        p("Below are QC plots for your uploaded Seurat object. You can download each plot individually."),
        fluidRow(
            column(4, plotBox("nFeature_plot", "download_nFeature", "nFeature_RNA (Pre-filter)", ns)),
            column(4, plotBox("nCount_plot", "download_nCount", "nCount_RNA (Pre-filter)", ns)),
            column(4, plotBox("percent_mt_plot", "download_percent_mt", "percent.mt (Pre-filter)", ns))
        ),
        hr(),
        h3(style = "text-align: center;", "QC Filtering Options"),
        p(style = "text-align: center;",
          "For 10x Genomics data, typical filtering thresholds are ~250 to 2500-3000 genes per cell."),
        fluidRow(
            column(4, numericInput(ns("min_features"), "Minimum Genes Per Cell (nFeature_RNA)", value = 250)),
            column(4, numericInput(ns("max_features"), "Maximum Genes Per Cell (nFeature_RNA)", value = 2500)),
            column(4, numericInput(ns("max_percent_mt"), "Maximum Percent Mitochondrial Genes", value = 5))
        ),
        div(style = "text-align: center;",
            actionButton(ns("apply_filters"), "Apply Filters", class = "btn btn-primary")
        ),
        verbatimTextOutput(ns("filter_summary")),
        hr(),
        uiOutput(ns("filtered_plots_ui")),
        uiOutput(ns("proceed_button_ui"))
    )
}

library(Seurat)
library(ggplot2)

QC_analysisServer <- function(id, seuratData,onProceed = NULL) {
    moduleServer(id, function(input, output, session) {
        ns <- session$ns

        qc_seurat <- reactiveVal(NULL)
        pre_qc_plots <- reactiveVal(list())
        filtered_seurat <- reactiveVal(NULL)

        # Cache for post-filter plots
        post_qc_plots <- reactiveVal(list())

        observeEvent(seuratData(), {
            obj <- calculateQCMetrics(seuratData())
            qc_seurat(obj)

            pre_qc_plots(list(
                nFeature = VlnPlot(obj, features = "nFeature_RNA", pt.size = 0.1) +
                    ggtitle("nFeature_RNA (Pre-filter)"),
                nCount = VlnPlot(obj, features = "nCount_RNA", pt.size = 0.1) +
                    ggtitle("nCount_RNA (Pre-filter)"),
                percent_mt = VlnPlot(obj, features = "percent.mt", pt.size = 0.1) +
                    ggtitle("percent.mt (Pre-filter)")
            ))

            # Clear filtered data and post-filter plots
            filtered_seurat(NULL)
            post_qc_plots(NULL)
        })

        # Pre-filter plots
        output$nFeature_plot <- renderPlot({ req(pre_qc_plots()); pre_qc_plots()$nFeature })
        output$nCount_plot <- renderPlot({ req(pre_qc_plots()); pre_qc_plots()$nCount })
        output$percent_mt_plot <- renderPlot({ req(pre_qc_plots()); pre_qc_plots()$percent_mt })

        # Pre-filter download handlers
        output$download_nFeature <- downloadHandler(
            filename = function() { "PreFilter_nFeature_plot.png" },
            content = function(file) { ggsave(file, plot = pre_qc_plots()$nFeature, width = 6, height = 4) }
        )
        output$download_nCount <- downloadHandler(
            filename = function() { "PreFilter_nCount_plot.png" },
            content = function(file) { ggsave(file, plot = pre_qc_plots()$nCount, width = 6, height = 4) }
        )
        output$download_percent_mt <- downloadHandler(
            filename = function() { "PreFilter_percent_mt_plot.png" },
            content = function(file) { ggsave(file, plot = pre_qc_plots()$percent_mt, width = 6, height = 4) }
        )

        # Apply filters and cache post-filter plots
        observeEvent(input$apply_filters, {
            req(qc_seurat())

            obj <- qc_seurat()
            keep_cells <- obj$nFeature_RNA >= input$min_features &
                          obj$nFeature_RNA <= input$max_features &
                          obj$percent.mt <= input$max_percent_mt

            filtered_obj <- subset(obj, cells = colnames(obj)[keep_cells])
            calculateQCMetrics(filtered_obj)
            filtered_seurat(filtered_obj)

            # Cache the post-filter plots
            post_qc_plots(list(
                nFeature = VlnPlot(filtered_obj, features = "nFeature_RNA", pt.size = 0.1) +
                    ggtitle("nFeature_RNA (Post-filter)"),
                nCount = VlnPlot(filtered_obj, features = "nCount_RNA", pt.size = 0.1) +
                    ggtitle("nCount_RNA (Post-filter)"),
                percent_mt = VlnPlot(filtered_obj, features = "percent.mt", pt.size = 0.1) +
                    ggtitle("percent.mt (Post-filter)")
            ))

            message(sprintf("Filtered: %d cells retained out of %d", sum(keep_cells), ncol(obj)))
        })

        # Filter summary
        output$filter_summary <- renderPrint({
            req(filtered_seurat())
            cat(sprintf("Filtered dataset has %d cells (original: %d cells).",
                        ncol(filtered_seurat()), ncol(qc_seurat())))
        })

        # Render post-filter plots in UI
        output$filtered_plots_ui <- renderUI({
            req(filtered_seurat())
            tagList(
                h3("Post-Filtering QC Plots", style = "text-align: center;"),
                fluidRow(
                    column(4, plotBox("filtered_nFeature_plot", "download_filtered_nFeature", "nFeature_RNA (Post-filter)", ns)),
                    column(4, plotBox("filtered_nCount_plot", "download_filtered_nCount", "nCount_RNA (Post-filter)", ns)),
                    column(4, plotBox("filtered_percent_mt_plot", "download_filtered_percent_mt", "percent.mt (Post-filter)", ns))
                )
            )
        })

        # Post-filter plots
        output$filtered_nFeature_plot <- renderPlot({ req(post_qc_plots()); post_qc_plots()$nFeature })
        output$filtered_nCount_plot <- renderPlot({ req(post_qc_plots()); post_qc_plots()$nCount })
        output$filtered_percent_mt_plot <- renderPlot({ req(post_qc_plots()); post_qc_plots()$percent_mt })

        # Post-filter download handlers (use cached plots)
        output$download_filtered_nFeature <- downloadHandler(
            filename = function() { "PostFilter_nFeature_plot.png" },
            content = function(file) { ggsave(file, plot = post_qc_plots()$nFeature, width = 6, height = 4) }
        )
        output$download_filtered_nCount <- downloadHandler(
            filename = function() { "PostFilter_nCount_plot.png" },
            content = function(file) { ggsave(file, plot = post_qc_plots()$nCount, width = 6, height = 4) }
        )
        output$download_filtered_percent_mt <- downloadHandler(
            filename = function() { "PostFilter_percent_mt_plot.png" },
            content = function(file) { ggsave(file, plot = post_qc_plots()$percent_mt, width = 6, height = 4) }
        )

        # Proceed button (only show after filtering)
        output$proceed_button_ui <- renderUI({
            req(filtered_seurat())
            div(style = "text-align: center; margin-top: 20px;",
                actionButton(ns("proceed_normalization"), "Proceed to Normalization and Clustering", class = "btn btn-success"))
        })
        # Proceed button (only show after filtering)
        output$proceed_button_ui <- renderUI({
            req(filtered_seurat())
            div(style = "text-align: center; margin-top: 20px;",
                actionButton(ns("proceed_normalization"), "Proceed to Normalization and Clustering", class = "btn btn-success"))
        })

        # Trigger `onProceed()` callback to pass filtered object and navigate
        observeEvent(input$proceed_normalization, {
            req(filtered_seurat())
            if (!is.null(onProceed)) {
                onProceed(filtered_seurat())
            }
        })
    })
}
