
#GUARD Finder Shiny app

library(shiny)

# Define UI for application

ui <- fluidPage(
    hr(),
    # Application title
    titlePanel("GUARD Finder"),
    hr(),
    em("CRISPR GUARD is a tool to reduce off-target editing by Cas9 and base editors.",
       p("Short guide RNAs called \"GUARD RNAs\" recruit Cas9 complexes to off-target sites but do not permit nuclease activty, 
         thereby protecting them from the mismatched guide RNA by direct competition.", 
         p("The GUARD Finder tool searches for guide RNA off-targets and designs GUARD RNAs to protect them from editing."
           )
         )
       ),
    br(),
    # Sidebar with a input 
    sidebarLayout(
        sidebarPanel(
            textInput("id",
                        "guide RNA name"),
        textInput("guide",
                  "guide RNA sequence without PAM"),
        selectInput("genome",
                  "genome version, human or mouse", c("hg38", "mm10")),
        selectInput("pam",
                  "PAM", c("NGG", "NGRRT", "NGRRN", "TTTN")),
        sliderInput(inputId = "guide_mismatches",
                  label = "maximum # mismatches for the guide off-target search", 
                  value = 4, min = 0, max = 5),
        selectInput("guide_min_pvalue",
                  "p value filter for which off-targets to take forward to guard design", c("0.1", "0.05", "0.01")),
        selectInput("guard_length",
                    "GUARD RNA length", c("15", "14")),
        sliderInput(inputId = "guard_mismatches",
                    label = "maximum # mismatches for the GUARD RNA off-target search",
                    value = 3, min = 0, max = 4),
        sliderInput(inputId = "max_guard_distance", 
        label = "maximum distance between off-target and GUARD RNA binding",
        value = 10, min = 0, max = 15),
        actionButton(inputId = "submit", label = "submit job"),
        hr(),
        textInput("chr",
                  "Special case: If your guide RNA maps perfectly to more than one genomic position, please specify the intended on-target locus: chr e.g. chr4"),
        textInput("start",
                  "start"),
        textInput("end",
                  "end"),
        selectInput("strand",
                  "strand", c("","+", "-"))
        ),
        
        mainPanel(
        tableOutput("Results"),
        br(),
        img(height= 280, width = 450, src = "CRISPR_GUARD.png"),
        hr(),
           h5("please cite: Coelho et al., CRISPR GUARD: short guide RNAs protect off target sites from Cas9 nuclease activity, Nature Communications, 2020"),
           a("GUARD Finder GitHub", href = "https://github.com/MatthewACoelho/GUARDfinder"),
        hr(),
           h5("License"),
           h6("Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the \"Software\"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software. 
THE SOFTWARE IS PROVIDED \"AS IS\", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
The Software is not intended for clinical use.")
         )
    )
)

# Define server logic required to make results table
server <- function(input, output) {
    
    #make the params.nf file from input fields 
    parameters <- eventReactive(input$submit, {
       params<- print(paste(c("params {", "\n",
          "id = ", "\"", input$id, "\"", "\n",
          "guide = ", "\"", input$guide,"\"", "\n",
          "genome = ", "\"", input$genome,"\"", "\n",
          "pam = ", "\"", input$pam,"\"", "\n",
          "guide_mismatches = ", "\"", input$guide_mismatches,"\"", "\n",
          "guide_min_pvalue = ", "\"", input$guide_min_pvalue,"\"", "\n",
          "guard_length = ", "\"", input$guard_length,"\"", "\n",
          "guard_mismatches = ", "\"", input$guard_mismatches,"\"", "\n",
          "max_guard_distance = ", "\"", input$max_guard_distance,"\"", "\n",
          "chr = ", "\"", input$chr,"\"", "\n",
          "start = ", "\"", input$start,"\"", "\n",
          "end = ", "\"", input$end,"\"", "\n",
          "strand = ", "\"", input$strand,"\"", "\n",
          "}"), collapse=""))
       write(params, file ="params.nf")
       print(params)
    })

    #show params in the shiny app
    output$Results <- renderText({
        parameters()
    })
}

# Run the application 
shinyApp(ui = ui, server = server)

