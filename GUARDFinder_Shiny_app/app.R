
#CRISPR GUARD Finder Shiny app

#load libraries
library(shiny)
library(DT)
library(tidyverse)
library(digest)

nucleotides <- c("A", "a", "C", "c", "T", "t", "G", "g")
special_chars <- c("\\/", "\\*", "\\,", "\\\\", "\\-", " ", "\\^", "\\|", "\\>")

# Define UI for application
ui <- fluidPage(
    hr(),
    #title
    titlePanel("CRISPR GUARD Finder"),
    hr(),
    #description
    em("CRISPR GUARD is a tool to reduce off-target editing by Cas9 and base editors.",
       p("Short guide RNAs called \"GUARD RNAs\" recruit Cas9 complexes to off-target sites but do not permit nuclease activty, 
         thereby protecting them from the mismatched guide RNA by direct competition.", 
         p("The CRISPR GUARD Finder tool searches for guide RNA off-targets and designs GUARD RNAs to protect them from editing."
           )
         )
       ),
    br(),
    # Sidebar with input 
    sidebarLayout(
        sidebarPanel(
            textInput("id",
                        "guide RNA name (e.g. HBB)"),
        textInput("guide",
                  "guide RNA sequence without PAM (e.g. CTTGCCCCACAGGGCAGTAA)"),
        selectInput("genome",
                  "genome version - human", c("hg38")),
        selectInput("pam",
                  "PAM", c("NGG")),
        sliderInput(inputId = "guide_mismatches",
                  label = "mismatches for the guide RNA off-target search", 
                  value = 3, min = 0, max = 5),
        selectInput("guide_min_pvalue",
                  "probability threshold for the guide RNA off-target search (e.g. 0.1 = 10 % likely)", c("0.1", "0.05", "0.01")),
        selectInput("guard_length",
                    "GUARD RNA length (nt)", c("15", "14")),
        sliderInput(inputId = "max_guard_distance", 
        label = "distance between GUARD RNA and off-target (bp)",
        value = 10, min = 0, max = 15),
        actionButton(inputId = "submit", label = "submit job"),
        hr(),
        downloadButton("downloadData", "download results"),
        hr(),
        textInput("chr",
                  "special case: If your guide RNA maps perfectly to more than one genomic position, please specify the intended on-target locus: chr e.g. chr4"),
        textInput("start",
                  "start"),
        textInput("end",
                  "end"),
        selectInput("strand",
                  "strand", c("","+", "-"))
        ),
        
        #MAIN
        mainPanel(
            br(),
            #image
            img(height= 280, width = 450, src = "CRISPR_GUARD.png"),
            br(),
            
            #print parameters
            br(),
            textOutput("parameters"),
            hr(),
            
            #results
            DT::dataTableOutput("results"),
            
            #MIT license etc with footnotes
            br(),
            hr(),
            h5("please cite: Coelho et al., CRISPR GUARD: short guide RNAs protect off-target sites from Cas9 nuclease activity, Nature Communications, 2020"),
            a("GUARD Finder GitHub", href = "https://github.com/MatthewACoelho/GUARDfinder"),
            hr(),
            h5("License"),
            h6("Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the \"Software\"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:
            The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.
            THE SOFTWARE IS PROVIDED \"AS IS\", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE."),
            hr(),
            h6("Shiny app development by Matt Coelho, CRISPR GUARD Finder code by Mike Firth. Thanks to Donny van de Meer and Jaymin Mistry for advice.")
         )
    )
)

# Define server logic required to make results table
server <- function(input, output) {
    
    #make the params.nf file from input fields 
    parameters <- eventReactive(input$submit, {
        
        # taint-check text input. guide should be 19 or 20 nucleotides ATGC only. id should be <20 chars with no special chars
        validate(
                need(str_length(input$guide) == 19 | str_length(input$guide) == 20 & str_length(input$guide) == sum(str_count(input$guide, nucleotides)), "invalid guide RNA input. 19 or 20 nt"),
                need(str_length(input$id) < 20 & sum(str_count(input$id, special_chars)) == 0, "invalid guide RNA name. please remove special characters or reduce length")
                )
        params <- print(paste(c("params {", "\n",
          "id = ", "\"", input$id, "\"", "\n",
          "guide = ", "\"", input$guide,"\"", "\n",
          "genome = ", "\"", input$genome,"\"", "\n",
          "pam = ", "\"", input$pam,"\"", "\n",
          "guide_mismatches = ", "\"", input$guide_mismatches,"\"", "\n",
          "guide_min_pvalue = ", "\"", input$guide_min_pvalue,"\"", "\n",
          "guard_length = ", "\"", input$guard_length,"\"", "\n",
          "guard_mismatches = ", "\"", "2","\"", "\n",
          "max_guard_distance = ", "\"", input$max_guard_distance,"\"", "\n",
          "chr = ", "\"", input$chr,"\"", "\n",
          "start = ", "\"", input$start,"\"", "\n",
          "end = ", "\"", input$end,"\"", "\n",
          "strand = ", "\"", input$strand,"\"", "\n",
          "}"), collapse=""))
    folder_id <- sha1(params)
    path <- paste0(getwd(), "/", folder_id)
    dir.create(path)
    setwd(path)
    write(params, file ="params.nf")
    print(params)
    })
    
    results <- eventReactive(input$submit, {
        withProgress(message = "Finding GUARDs ... ",
                     value = 0.3, 
                     detail = "for gRNAs with hundreds of off-targets, this can take several minutes", {
        if (file.exists("params.nf"))
            system("/Users/mc32/guard_root/bin/find_guards.sh", intern=TRUE)
            results_filename <- print(paste(c(input$id, "_final.txt"), collapse = ""))
                validate(
                        need(file.exists(results_filename), "no results for this gRNA query ... please relax search criteria")
                        )
                data <- as_tibble(read.delim(file = results_filename, header = TRUE, sep = "\t"))
                data <- data %>% mutate(OffTarget=OffGuide, GuardOffTargets0Mismatch=X0, GuardOffTargets1Mismatch=X1, GuardOffTargets2Mismatch=X2)
                data <- data %>% select("OffTarget","ID", "OffGuideStrand", "Guard", "GuardGC", "GuardStrand",
                                        "ForwardGuardWithPAM", "GuardOffTargets0Mismatch" , 
                                        "GuardOffTargets1Mismatch" ,"GuardOffTargets2Mismatch", "GuardChr","GuardStart", "GuardStop")
    
                })
    })
    
    #Render functions
    #show params in the shiny app
    output$parameters <- renderText({
       parameters()
    })
    
    #show results data table in the shiny app
    output$results <- DT::renderDataTable({
        results()
    })
    
    # Downloadable txt of data
    output$downloadData <- downloadHandler(
        filename = function() {
            paste(c(input$id, "_final.txt"), collapse = "")
        },
        content = function(file) {
            write.table(results(), file, sep ="\t")
        }
    )
}

# Run the application 
shinyApp(ui = ui, server = server)

