makeGenomeTrackPlot <- function(GenomeAxisTrack,
                                  GenomeAxisTrack2,
                                  IdeogramTrackList,
                                  TranscriptTrack,
                                  Chromosome,
                                  Start,
                                  End,
                                  parameters,
                                  SampleLabel,
                                  SampleLabels,
                                  relativeLineWidthGenomeAxis,
                                  relativeLineWidthHist,
                                  GenomeAxisFontSize,
                                  LabelFontScale,
                                  IdeogramMarkerSize){
    options(ucscChromosomeNames=FALSE)
    rl <- GenomicRanges::GRanges(seqnames = Chromosome,
                                 ranges = IRanges::IRanges(start=Start,
                                                           end=End))
    bws <- rtracklayer::BigWigSelection(rl)
    
    
    ###
    
    parameters$SampleFilename <- as.character(parameters$SampleFilename)
    parameters$ControlFilename <- as.character(parameters$ControlFilename)
    
    bigwigParameters <- parameters[grep( ".bw$", parameters$SampleFilename, ignore.case=TRUE),]
    bedParameters <- parameters[grep( ".bed$", parameters$SampleFilename, ignore.case=TRUE),]
    
    SampleBwFiles <- lapply(
      bigwigParameters$SampleFilename,
      rtracklayer::import.bw,selection=bws)
    
    names(SampleBwFiles) <- bigwigParameters$SampleLabel
    
    ControlBwFiles <- lapply(
      bigwigParameters$ControlFilename,
      rtracklayer::import.bw,selection=bws)
    
    names(ControlBwFiles) <- bigwigParameters$SampleLabel
    
    ### add y limits
    
    parameters$yLimits <- sapply(SampleLabels,function(lbl){
      round(max(c(SampleBwFiles[[lbl]]$score,ControlBwFiles[[lbl]]$score)),0)
    })
    
    ### 
    
    SampleBedFiles <- lapply(
      bedParameters$SampleFilename,function(nme){
        gr <- rtracklayer::import.bed(nme)
        GenomicRanges::strand(gr) <- "*"
        IRanges::subsetByOverlaps(gr,rl)
      }
    )
    
    names(SampleBedFiles) <- bedParameters$SampleLabel
    
    ###
    
    
    SampleDataTracks <- lapply(
      parameters$SampleLabel,
      function(lbl){
        rw <- parameters[which(parameters$SampleLabel==lbl),]
        if(grepl( ".bw$", rw$SampleFilename, ignore.case=TRUE)){
          Gviz::OverlayTrack(
            trackList = list(
              Gviz::DataTrack(
                range = SampleBwFiles[[rw$SampleLabel]],
                genome = "hg19", 
                type = "hist", 
                col=rw$SampleColor,
                col.histogram=rw$SampleColor,
                fill.histogram=rw$SampleColor,
                ylim=c(0,rw$yLimits),
                lwd = relativeLineWidthHist,
                legend=TRUE,
                name=rw$SampleLabel,
                showAxis=FALSE,
                background.title="blue",
                cex.title=LabelFontScale),
              Gviz::DataTrack(
                range = ControlBwFiles[[rw$SampleLabel]],
                genome = "hg19", 
                type = "hist", 
                lwd = relativeLineWidthHist,
                col=rw$ControlColor,
                col.histogram=rw$ControlColor,
                fill.histogram=rw$ControlColor,
                ylim=c(0,rw$yLimits),
                legend=TRUE,
                name=rw$ControlLabel,
                showAxis=FALSE,
                background.title="blue",
                cex.title=LabelFontScale)
            ),
            name = rw$SampleLabel,
            background.title = "black",
            cex.title=LabelFontScale
          )
        }else
        {
          Gviz::AnnotationTrack(
            range = SampleBedFiles[[rw$SampleLabel]],
            name=rw$SampleLabel,
            col=rw$SampleColor,
            fill=rw$SampleColor,
            col.line=rw$SampleColor,
            background.title = "black",
            cex.title=LabelFontScale
          )
        }
      })
    
    
    ideogramTrack <- IdeogramTrackList[[Chromosome]]
    
    #ideogramTrack <- Gviz::IdeogramTrack(genome = "hg19", chromosome = Chromosome)
    
    Gviz::displayPars(GenomeAxisTrack)$lwd <- relativeLineWidthGenomeAxis 
    Gviz::displayPars(ideogramTrack)$lwd <- IdeogramMarkerSize
    Gviz::displayPars(GenomeAxisTrack)$col <- "black" 
    Gviz::displayPars(GenomeAxisTrack)$fontsize <- GenomeAxisFontSize
    Gviz::displayPars(GenomeAxisTrack)$cex.title <- LabelFontScale
    Gviz::displayPars(ideogramTrack)$cex.title <- LabelFontScale
    Gviz::displayPars(TranscriptTrack)$cex.title <- LabelFontScale
    Gviz::displayPars(GenomeAxisTrack)$fontcolor <- "black"
    LabelFontScale
    
    
    trcKList <- list(ideogramTrack,TranscriptTrack,GenomeAxisTrack,GenomeAxisTrack2)
    
    # trcKList <- lapply(trcKList,function(trck){
    #   Gviz::displayPars(trck)$fontsize <- fontSize
    #   Gviz::displayPars(trck)$lwd <- relativeLineWidth
    # })
    
    Gviz::plotTracks(c(trcKList[[1]],
                       SampleDataTracks,
                       trcKList[[2]],
                       trcKList[[3]],
                       trcKList[[4]]),
                     from=Start,
                     to=End)
  }

################################################################################

plotFunctionText <- quote(
    makeGenomeTrackPlot <- function(GenomeAxisTrack,
                                    GenomeAxisTrack2,
                                    IdeogramTrackList,
                                    TranscriptTrack,
                                    Chromosome,
                                    Start,
                                    End,
                                    parameters,
                                    SampleLabel,
                                    SampleLabels,
                                    relativeLineWidthGenomeAxis,
                                    relativeLineWidthHist,
                                    GenomeAxisFontSize,
                                    LabelFontScale,
                                    IdeogramMarkerSize){
      options(ucscChromosomeNames=FALSE)
      rl <- GenomicRanges::GRanges(seqnames = Chromosome,
                                   ranges = IRanges::IRanges(start=Start,
                                                             end=End))
      bws <- rtracklayer::BigWigSelection(rl)
      
      
      ###
      
      parameters$SampleFilename <- as.character(parameters$SampleFilename)
      parameters$ControlFilename <- as.character(parameters$ControlFilename)
      
      bigwigParameters <- parameters[grep( ".bw$", parameters$SampleFilename, ignore.case=TRUE),]
      bedParameters <- parameters[grep( ".bed$", parameters$SampleFilename, ignore.case=TRUE),]
      
      SampleBwFiles <- lapply(
        bigwigParameters$SampleFilename,
        rtracklayer::import.bw,selection=bws)
      
      names(SampleBwFiles) <- bigwigParameters$SampleLabel
      
      ControlBwFiles <- lapply(
        bigwigParameters$ControlFilename,
        rtracklayer::import.bw,selection=bws)
      
      names(ControlBwFiles) <- bigwigParameters$SampleLabel
      
      ### add y limits
      
      parameters$yLimits <- sapply(SampleLabels,function(lbl){
        round(max(c(SampleBwFiles[[lbl]]$score,ControlBwFiles[[lbl]]$score)),0)
      })
      
      ### 
      
      SampleBedFiles <- lapply(
        bedParameters$SampleFilename,function(nme){
          gr <- rtracklayer::import.bed(nme)
          GenomicRanges::strand(gr) <- "*"
          IRanges::subsetByOverlaps(gr,rl)
        }
      )
      
      names(SampleBedFiles) <- bedParameters$SampleLabel
      
      ###
      
      
      SampleDataTracks <- lapply(
        parameters$SampleLabel,
        function(lbl){
          rw <- parameters[which(parameters$SampleLabel==lbl),]
          if(grepl( ".bw$", rw$SampleFilename, ignore.case=TRUE)){
            Gviz::OverlayTrack(
              trackList = list(
                Gviz::DataTrack(
                  range = SampleBwFiles[[rw$SampleLabel]],
                  genome = "hg19", 
                  type = "hist", 
                  col=rw$SampleColor,
                  col.histogram=rw$SampleColor,
                  fill.histogram=rw$SampleColor,
                  ylim=c(0,rw$yLimits),
                  lwd = relativeLineWidthHist,
                  legend=TRUE,
                  name=rw$SampleLabel,
                  showAxis=FALSE,
                  background.title="blue",
                  cex.title=LabelFontScale),
                Gviz::DataTrack(
                  range = ControlBwFiles[[rw$SampleLabel]],
                  genome = "hg19", 
                  type = "hist", 
                  lwd = relativeLineWidthHist,
                  col=rw$ControlColor,
                  col.histogram=rw$ControlColor,
                  fill.histogram=rw$ControlColor,
                  ylim=c(0,rw$yLimits),
                  legend=TRUE,
                  name=rw$ControlLabel,
                  showAxis=FALSE,
                  background.title="blue",
                  cex.title=LabelFontScale)
              ),
              name = rw$SampleLabel,
              background.title = "black",
              cex.title=LabelFontScale
            )
          }else
          {
            Gviz::AnnotationTrack(
              range = SampleBedFiles[[rw$SampleLabel]],
              name=rw$SampleLabel,
              col=rw$SampleColor,
              fill=rw$SampleColor,
              col.line=rw$SampleColor,
              background.title = "black",
              cex.title=LabelFontScale
            )
          }
        })
      
      
      ideogramTrack <- IdeogramTrackList[[Chromosome]]
      
      Gviz::displayPars(GenomeAxisTrack)$lwd <- relativeLineWidthGenomeAxis 
      Gviz::displayPars(ideogramTrack)$lwd <- IdeogramMarkerSize
      Gviz::displayPars(GenomeAxisTrack)$col <- "black" 
      Gviz::displayPars(GenomeAxisTrack)$fontsize <- GenomeAxisFontSize
      Gviz::displayPars(GenomeAxisTrack)$cex.title <- LabelFontScale
      Gviz::displayPars(ideogramTrack)$cex.title <- LabelFontScale
      Gviz::displayPars(TranscriptTrack)$cex.title <- LabelFontScale
      Gviz::displayPars(GenomeAxisTrack)$fontcolor <- "black"
      LabelFontScale
      
      
      trcKList <- list(ideogramTrack,TranscriptTrack,GenomeAxisTrack)
      
      # trcKList <- lapply(trcKList,function(trck){
      #   Gviz::displayPars(trck)$fontsize <- fontSize
      #   Gviz::displayPars(trck)$lwd <- relativeLineWidth
      # })
      
      Gviz::plotTracks(c(trcKList[[1]],
                         SampleDataTracks,
                         trcKList[[2]],
                         trcKList[[3]]),
                       from=Start,
                       to=End)
    }
  )
