#!/usr/bin/env Rscript

# Load optparse we need to check inputs

suppressPackageStartupMessages(require(optparse))

# Load common functions

suppressPackageStartupMessages(require(workflowscriptscommon))

# parse options
option_list = list(
  make_option(
    c("--object-1"),
    action = "store",
    default = NA,
    type = 'character',
    help = "A Seurat object."
  ),
  make_option(
    c("--object-2"),
    action = "store",
    default = NA,
    type = 'character',
    help = "Counts Matrix."
  ),
  make_option(
    c("--integration-type"),
    action = "store",
    default = NA,
    type = 'character',
    help = "Assay Name."
  ),
  make_option(
    c("--output-object-file"),
    action = "store",
    default = NA,
    type = 'character',
    help = "File name in which to store serialized R matrix object."
  ) 
)
opt <- wsc_parse_args(option_list)

suppressPackageStartupMessages(require(Seurat))

set.seed(1234)

if (! file.exists(opt$object_1)){
  stop((paste('File', opt$object_1, 'does not exist')))
}
if (! file.exists(opt$object_2)){
  stop((paste('File', opt$object_2, 'does not exist')))
}

obj1 <- readRDS(file = opt$object_1)
obj2 <- readRDS(file = opt$object_2)
integration_type <- opt$integration_type

obj <- list(obj1, obj2)
print(0)
obj <- lapply(X = obj, FUN = function(x) {
    x <- NormalizeData(x)
    x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})
print(0)
features <- SelectIntegrationFeatures(object.list = obj)
print(0)
anchors <- FindIntegrationAnchors(object.list = obj, anchor.features = features)
integrated <- IntegrateData(anchorset = anchors)

saveRDS(integrated, file = opt$output_object_file)


