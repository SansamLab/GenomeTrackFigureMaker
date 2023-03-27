# GenomeTrackFigureMaker

## What does this app do?
This Shiny application creates a genome track plot based on a user-provided CSV file. The server code defines several reactive expressions and observers that update plot parameters based on user input. The output plot is generated using the makeGenomeTrackPlot() function, which is called with various input parameters based on user selections. Finally, the metaCode output provides the R code for generating the plot, which can be useful for debugging and reproducing the plot outside of the Shiny app.

## Description of columns in .csv parameters file
The .csv file should contain the following information for each track to be plotted:

- "SampleLabel": A label for the sample track.
- "SampleFilename": The filename or path of the sample data file in a format supported by Gviz (e.g., bigWig, bedGraph, BAM, etc.).
- "SampleColor": The color to use for the sample track in the plot.
- "ControlLabel": A label for the control track (if applicable).
- "ControlFilename": The filename or path of the control data file in a format supported by Gviz (e.g., bigWig, bedGraph, BAM, etc.) (if applicable).
- "ControlColor": The color to use for the control track in the plot (if applicable).

Each row in the .csv file corresponds to a single track to be plotted, and the "SampleLabel" should match the label of the track selected by the user in the Shiny app. The file should be structured such that each track has its own row, and the columns are comma-separated.



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
