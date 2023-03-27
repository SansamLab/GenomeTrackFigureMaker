# This Shiny application creates a genome track plot based on a user-provided CSV
# file. Users can select a chromosome and genomic region and choose tracks to plot
# using the uploaded parameter file. The plot can be downloaded as a PDF. The app
# uses Gviz and shinymeta packages for visualization and reactive programming. 
# Numeric inputs control plot parameters, file input is used for the parameter 
# file, and download button and plot output are included in the UI. Reactive 
# expressions and observers update plot parameters based on user input. 
# makeGenomeTrackPlot function generates the output plot based on user selections. 
# The metaCode output provides the R code for generating the plot for debugging 
# and reproducing the plot outside the Shiny app.

# The .csv file should contain the following information for each track to be plotted:
#
# "SampleLabel": A label for the sample track.
# "SampleFilename": The filename or path of the sample data file in a format supported by Gviz (e.g., bigWig, bedGraph, BAM, etc.).
# "SampleColor": The color to use for the sample track in the plot.
# "ControlLabel": A label for the control track (if applicable).
# "ControlFilename": The filename or path of the control data file in a format supported by Gviz (e.g., bigWig, bedGraph, BAM, etc.) (if applicable).
# "ControlColor": The color to use for the control track in the plot (if applicable).
#
# Each row in the .csv file corresponds to a single track to be plotted, and the 
# "SampleLabel" should match the label of the track selected by the user in the Shiny app. 
# The file should be structured such that each track has its own row, and the columns are comma-separated.


library(shiny)

################################################################################
#                                                                              #
#                                Shiny UI Code                                 #
#                                                                              #
################################################################################

ui <- fluidPage(
  
  # Add a spinner to indicate that the app is busy
  fluidPage(
    shinybusy::add_busy_spinner(spin = "fading-circle",
                                height="50px",
                                width="50px",
                                timeout=1000),
    
    # File input to select the CSV parameter file for bigwigs
    fluidRow(
      shiny::fileInput("parameterFle","Select .csv parameter file for bigwigs:",accept=".csv")
    ),
    
    # UI for selecting genomic region
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
    
    # Selectize input for choosing tracks to plot
    wellPanel(fluidRow(
      selectizeInput("smpIDs",
                     "Tracks to plot:",
                     choices="upload parameter file",
                     multiple=T))),
    
    # UI for various plot parameters
    wellPanel(
      fluidRow(
        column(2,shiny::numericInput("relLineWidthGenomeAxis",label = "Genome Axis Line Width:",value=1,min=0)),
        column(2,shiny::numericInput("relLineWidthHist",label = "Histogram Line Width:",value=1,min=0)),
        column(2,shiny::numericInput("IdeogramMarkerSize",label = "Ideogram marker size:",value=1,min=0)),
        column(2,shiny::numericInput("fontSize",label = "Font size:",value=10,min=0)),
        column(2,shiny::numericInput("fontScale",label = "Label font scale:",value=1,min=0)))),
    
    # Plot output
    wellPanel(fluidRow(
      plotOutput("plot")
    )),
    
    # Button to download the plot
    wellPanel(fluidRow(
      column(4,shiny::numericInput("width", label = "Width:",value=4,min=0)),
      column(4,shiny::numericInput("height",label = "Height:",value=2,min=0)),
      column(4,downloadButton("downloadPlot","Download Plot"),
      ))),
    
    # R code output for generating the plot
    wellPanel(verbatimTextOutput("metaCode"))
  )
)

################################################################################
#                              Server Logic                                    #
# This is the server code for the Shiny application that creates a genome track #
# plot based on a user-provided CSV file. The server code defines several      #
# reactive expressions and observers that update plot parameters based on user #
# input. The output plot is generated using the makeGenomeTrackPlot function,   #
# which is called with various input parameters based on user selections.      #
# Finally, the metaCode output provides the R code for generating the plot,     #
# which can be useful for debugging and reproducing the plot outside of the     #
# Shiny app.                                                                   #
################################################################################

server <- function(input, output,session) {
  
  # Load necessary libraries and source helper functions
  library(Gviz) # for visualization
  library(shinymeta) # for reactive programming
  source("helpers.R") # helper functions for making genome track plot
  
  # Initialize reactive objects
  gtrack <- shinymeta::metaReactive({Gviz::GenomeAxisTrack()})
  gtrack2 <- shinymeta::metaReactive({Gviz::GenomeAxisTrack(scale=..(input$scaleSize))})
  txTr <- shinymeta::metaReactive({readRDS('data/hg38_txTr.rds')}) # transcript annotation track
  itrack_ls <- shinymeta::metaReactive({readRDS('data/itrack_ls.rds')}) # intron annotation track
  
  # Update height input based on number of selected tracks
  observeEvent(input$smpIDs,{
    cnt <- length(input$smpIDs) + 3
    updateNumericInput(session,"height",value=cnt)
  })
  
  # Update refined start/end slider based on input start value
  observeEvent(input$start,{
    strt <- input$start
    nd <- input$end
    updateSliderInput(session,"start2",value=c(strt,nd),min=strt,max=nd)
  })
  
  # Update refined start/end slider based on input end value
  observeEvent(input$end,{
    strt <- input$start
    nd <- input$end
    updateSliderInput(session,"start2",value=c(strt,nd),min=strt,max=nd)
  })
  
  # Compute height and width in pixels based on user inputs
  Height_pxls <- reactive({input$height*72})
  Width_pxls <- reactive({input$width*72})
  
  # Read in parameter file name, default to example.csv if no file selected
  parameterFileName <- shinymeta::metaReactive({
    if (is.null(input$parameterFle))
      "data/example.csv"
    else
      input$parameterFle$datapath
  })
  
  # Read in parameter file based on selected file name
  parameterFile <- shinymeta::metaReactive({read.csv(..(parameterFileName()))})
  
  # Update track choices in selectizeInput based on parameter file
  observeEvent(parameterFile(),{
    shiny::updateSelectizeInput(inputId="smpIDs",choices=parameterFile()$SampleLabel)
  })
  
  # Subset parameter file based on selected tracks
  parametersSubset <- reactive({
    req(parameterFile())
    parameterFile()[match(input$smpIDs,parameterFile()$SampleLabel),]
  })
  
  # Create genome track plot
  GenomeAxisPlot <- shinymeta::metaReactive2({
    req(parametersSubset())
    shinymeta::metaExpr({
      makeGenomeTrackPlot(
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
        ..(input$IdeogramMarkerSize))
    })
  })
  
  # Render plot output
  output$plot <- renderPlot({
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
  
  # Download plot as pdf
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
  
  # Render meta code
  output$metaCode <- renderPrint({
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
