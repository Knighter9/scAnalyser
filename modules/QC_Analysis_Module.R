QC_analysisUI <- function(id) {
    ns <- NS(id)
    tagList(
        h3("QC Analysis"),
        p("This page will be implemented later. Stay tuned for cool QC plots and filters!")
    )
}

# Server logic for QC analysis
# lets look at the seurat object and do some QC, 


# first lets calcualte 

QC_analysisServer <- function(id, seuratData) {
    moduleServer(id, function(input, output, session) {
        # Nothing to do yet, placeholder for future QC logic
    })
}
