# GenomeTrackFigureMaker

```r
GenomicFeatures::makeTxDbFromUCSC(
  genome="hg38",
  tablename="ccdsGene") %>%
  Gviz::GeneRegionTrack(.,collapseTracks=T,background.title = "black",name="Gene",fill="black") %>%
  saveRDS(.,file="hg38_txTr.rds")
```

```r
install.packages("shiny")
install.packages("shinybusy")
install("shinymeta")

if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("Gviz")
```


1.  Open R Studio.
2.  Set working directory to where your parameters file is stored.
```r
setwd("/DirectoryWhereMyParametersFileIsSaved")
```
3.  Download and run this App
```r
runGitHub( "GenomeTrackFigureMaker", "SansamLab",destdir = "../")
```
**_NOTE:_**  The defauault working directory when the app is run will be the location of the app.R file. Therefore, it is crucial that the relative paths of the app.R and parameters file are at the same level. This ensures that relationships of the paths of files listed in the parameters file are maintained.
