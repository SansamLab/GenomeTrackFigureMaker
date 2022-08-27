#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)

# Define UI for application that draws a histogram
ui <- fluidPage(

  fluidPage(
    shinybusy::add_busy_spinner(spin = "fading-circle",
                                height="50px",
                                width="50px",
                                timeout=1000),
    fluidRow(
      shiny::fileInput("parameterFle","Select .csv parameter file for bigwigs:",accept=".csv")
    ),
    fluidRow(
      shiny::wellPanel(
        column(2,shiny::selectInput("chrom", label = "Chromosome:",
                                    choices = (paste0("chr",1:22)), selected = "chr1")),
        column(2,shiny::numericInput("start", label = "Start:", value = 65272001)),
        column(2,shiny::numericInput("end", label = "End:",value=65370732)),
        column(2,shiny::numericInput("scaleSize",label = "Scale bar size",value=1000)),
        column(2,shiny::radioButtons("ideoChoice","Show ideogram?",choices=c("yes","no"))),
        shiny::wellPanel(
          shiny::sliderInput("start2","Refine Start and End:",
                             max=65370732,
                             min=65272001,
                             value=c(65272001,65370732),
                             width="100%")))),
    wellPanel(fluidRow(
      selectizeInput("smpIDs",
                     "Tracks to plot:",
                     choices="upload parameter file",
                     multiple=T))),
    wellPanel(
      fluidRow(
        column(2,shiny::numericInput("relLineWidthGenomeAxis",label = "Genome Axis Line Width:",value=1,min=0)),
        column(2,shiny::numericInput("relLineWidthHist",label = "Histogram Line Width:",value=1,min=0)),
        column(2,shiny::numericInput("IdeogramMarkerSize",label = "Ideogram marker size:",value=1,min=0)),
        column(2,shiny::numericInput("fontSize",label = "Font size:",value=10,min=0)),
        column(2,shiny::numericInput("fontScale",label = "Label font scale:",value=1,min=0))))
    
  )
)

# Define server logic required to draw a histogram
server <- function(input, output,session) {
  
  library(Gviz)
  library(shinymeta)
  
################################################################################
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
  

################################################################################

  gtrack <- shinymeta::metaReactive({Gviz::GenomeAxisTrack()})
  gtrack2 <- shinymeta::metaReactive({Gviz::GenomeAxisTrack(scale=..(input$scaleSize))})
  txTr <- shinymeta::metaReactive({readRDS('data/hg19_txTr.rds')})
  itrack_ls <- shinymeta::metaReactive({readRDS('data/itrack_ls.rds')})
  
  
  observeEvent(input$smpIDs,{
    cnt <- length(input$smpIDs) + 3
    updateNumericInput(session,"height",value=cnt)
  })
  
  observeEvent(input$start,{
    strt <- input$start
    nd <- input$end
    updateSliderInput(session,"start2",value=c(strt,nd),min=strt,max=nd)
  })
  
  observeEvent(input$end,{
    strt <- input$start
    nd <- input$end
    updateSliderInput(session,"start2",value=c(strt,nd),min=strt,max=nd)
  })
  
  Height_pxls <- reactive({input$height*72})
  Width_pxls <- reactive({input$width*72})
  
  #parameterFile <- reactive({read.csv(input$parameterFle$name)})
  
  parameterFile <- shinymeta::metaReactive({read.csv(..(input$parameterFle$name))})
  
  observeEvent(input$parameterFle,{
    shiny::updateSelectizeInput(inputId="smpIDs",choices=parameterFile()$SampleLabel)
  })
  
  parametersSubset <- reactive({
    req(input$parameterFle)
    parameterFile()[match(input$smpIDs,parameterFile()$SampleLabel),]
  })
  
  GenomeAxisPlot <- shinymeta::metaReactive2({
    req(parametersSubset())
    #strt <- input$start2[1]
    #nd <- input$start2[2]
    shinymeta::metaExpr({makeGenomeTrackPlot(
      ..(gtrack()),
      ..(gtrack2()),
      ..(itrack_ls()),
      ..(txTr()),
      ..(input$chrom),
      ..(input$start2[1]),
      ..(input$start2[2]),
      ..(parametersSubset()),
      ..(input$smpIDs),
      ..(input$smpIDs),
      ..(input$relLineWidthGenomeAxis),
      ..(input$relLineWidthHist),
      ..(input$fontSize),
      ..(input$fontScale),
      ..(input$IdeogramMarkerSize))})
  })
  
  
  renderPlot({
    req(parametersSubset())
    makeGenomeTrackPlot(gtrack(),
                        gtrack2(),
                        itrack_ls(),
                        txTr(),
                        input$chrom,
                        input$start2[1],
                        input$start2[2],
                        parametersSubset(),
                        input$smpIDs,
                        input$smpIDs,
                        input$relLineWidthGenomeAxis,
                        input$relLineWidthHist,
                        input$fontSize,
                        input$fontScale,
                        input$IdeogramMarkerSize)
  })
  
  output$downloadPlot <- downloadHandler(
    filename = function() { "plot.pdf" },
    content = function(file) {
      pdf(file,width=input$width,height=input$height,pointsize=16)
      makeGenomeTrackPlot(gtrack(),
                          gtrack2(),
                          itrack_ls(),
                          txTr(),
                          input$chrom,
                          input$start2[1],
                          input$start2[2],
                          parametersSubset(),
                          input$smpIDs,
                          input$smpIDs,
                          input$relLineWidthGenomeAxis,
                          input$relLineWidthHist,
                          input$fontSize,
                          input$fontScale,
                          input$IdeogramMarkerSize)
      dev.off()
    }
  )
  
  shiny::fluidPage(wellPanel(fluidRow(
    column(4,shiny::numericInput("width", label = "Width:",value=4,min=0)),
    column(4,shiny::numericInput("height",label = "Height:",value=2,min=0)),
    column(4,downloadButton("downloadPlot","Download Plot"),
    ))))
  
  renderPrint({
    shinymeta::expandChain(
      plotFunctionText,
      invisible(gtrack()),
      invisible(gtrack2()),
      invisible(txTr()),
      invisible(itrack_ls()),
      invisible(parameterFile()),
      GenomeAxisPlot()
    )
  })
    
}

# Run the application 
shinyApp(ui = ui, server = server)
