# Cell communication analysis

## prepare the CellChat database 

Rename the LR pairs for analysis the mouse and human cell interactions.

```r

# Function to capitalize the first letter of a string
capitalize_first <- function(x) {
  paste0(toupper(substr(x, 1, 1)), tolower(substr(x, 2, nchar(x))))
}

# Function to add a prefix to all elements of a dataframe
add_prefix <- function(df, prefix, upper = TRUE) {
  if (isTRUE(upper)) {
    df <- data.frame(apply(as.matrix(df), 2, function(x) {
      ifelse(x == "", x, paste0(prefix, toupper(x)))
    }))
    rownames(df) <- paste0(prefix, toupper(rownames(df)))
  } else {
    df <- data.frame(apply(as.matrix(df), 2, function(x) {
      ifelse(x == "", x, paste0(prefix, capitalize_first(x)))
    }))
    rownames(df) <- paste0(prefix, capitalize_first(rownames(df)))
  }
  return(df)
}

# Generate mouse to human database (db1)
db1 <- CellChatDB.mouse

## Modify interaction data
db1$interaction$interaction_name <- paste0('hg_', toupper(db1$interaction$interaction_name))
db1$interaction$pathway_name <- paste0('hg_', toupper(db1$interaction$pathway_name))
rownames(db1$interaction) <- db1$interaction$interaction_name
db1$interaction$ligand <- paste0('hg_', toupper(db1$interaction$ligand))
db1$interaction$receptor <- paste0('mm_', capitalize_first(db1$interaction$receptor))

## Modify complex and cofactor data
db1$complex <- rbind(
  add_prefix(db1$complex, prefix = 'hg_', upper = TRUE),
  add_prefix(db1$complex, prefix = 'mm_', upper = FALSE)
)
db1$cofactor <- rbind(
  add_prefix(db1$cofactor, prefix = 'hg_', upper = TRUE),
  add_prefix(db1$cofactor, prefix = 'mm_', upper = FALSE)
)
db1$geneInfo$Symbol <- paste0('mm_', db1$geneInfo$Symbol)

# Generate human to mouse database (db2)
db2 <- CellChatDB.human

## Modify interaction data
db2$interaction$interaction_name <- paste0('mm_', toupper(db2$interaction$interaction_name))
db2$interaction$pathway_name <- paste0('mm_', toupper(db2$interaction$pathway_name))
rownames(db2$interaction) <- db2$interaction$interaction_name
db2$interaction$ligand <- paste0('mm_', capitalize_first(db2$interaction$ligand))
db2$interaction$receptor <- paste0('hg_', toupper(db2$interaction$receptor))

## Modify complex and cofactor data
db2$complex <- rbind(
  add_prefix(db2$complex, prefix = 'hg_', upper = TRUE),
  add_prefix(db2$complex, prefix = 'mm_', upper = FALSE)
)
db2$cofactor <- rbind(
  add_prefix(db2$cofactor, prefix = 'hg_', upper = TRUE),
  add_prefix(db2$cofactor, prefix = 'mm_', upper = FALSE)
)
db2$geneInfo$Symbol <- paste0('hg_', db2$geneInfo$Symbol)

# Merge mouse and human data
db3 <- CellChatDB.human
db3$interaction <- rbind(db1$interaction, db2$interaction)
db3$complex <- rbind(db1$complex, db2$complex)
db3$cofactor <- rbind(db1$cofactor, db2$cofactor)
db3$geneInfo <- rbind(db1$geneInfo, db2$geneInfo)

# Save the merged database
saveRDS(db3, file = 'cellchatdb_mm_hg.Rds')


```

## Implement cell communication analysis

```r
library(Seurat)
library(CellChat)

dih = db_human@assays$SCT@data
rownames(dih) <- paste0('hg_', rownames(dih))
dim = db_mm@assays$SCT@data
rownames(dim) <- paste0('mm_', rownames(dim))

meta = data.frame(labels = obj_current$group, row.names = rownames(obj_current@meta.data))

spatial.locs = Seurat::GetTissueCoordinates(obj_current, scale = NULL, cols = c("imagerow", "imagecol"))

spot.size = 65 # 10x visium
scalefactors = jsonlite::fromJSON(
  txt = file.path(
    "/mnt/raid61/TNP_spatial/analysis/alignment/2024-0601-split/data/hg38/BNI_23_HCT/BNI_23_HCT/outs/spatial/",
    'scalefactors_json.json'
  )
)


conversion.factor = spot.size/scalefactors$spot_diameter_fullres
spatial.factors = data.frame(ratio = conversion.factor, tol = spot.size/2)

d.spatial <-
  computeCellDistance(coordinates = spatial.locs,
                      ratio = spatial.factors$ratio,
                      tol = spatial.factors$tol)


cellchat <-
  createCellChat(
    object = rbind(dih, dim),
    meta = meta,
    group.by = "labels",
    datatype = "spatial",
    coordinates = spatial.locs,
    spatial.factors = spatial.factors
  )


cellchat@DB <- readRDS('cellchatdb_mm_hg.Rds')

cellchat <- subsetData(cellchat)
future::plan("multisession", workers = 4) 
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat, variable.both = F)

cellchat <-
  computeCommunProb(
    cellchat,
    type = "truncatedMean",
    trim = 0.1,
    distance.use = TRUE,
    interaction.range = 250,
    scale.distance = 0.01,
    contact.dependent = TRUE,
    contact.range = 100
  )


cellchat <- filterCommunication(cellchat, min.cells = 10)
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)
saveRDS(cellchat, 'bni_cellchat.Rds')


```


