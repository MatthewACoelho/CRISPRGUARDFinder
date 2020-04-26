
#GUARD Finder Shiny app

library(shiny)

# Define UI for application

ui <- fluidPage(
    hr(),
    # Application title
    titlePanel("GUARD Finder"),
    hr(),
    em("CRISPR GUARD is a tool to reduce off-target editing by Cas9 and base editors.",
       p("Short guides called GUARD RNAs can recruit Cas9 complexes to off-target sites but do not permit nuclease acitivty, thus protecting them from editing.", 
         p("This tool designs GUARD RNAs for a given guide RNA sequence."
           )
         )
       ),
    br(),
    # Sidebar with a input 
    sidebarLayout(
        sidebarPanel(
            textInput("id",
                        "guide name"),
        textInput("guide",
                  "guide sequence without PAM"),
        selectInput("genome",
                  "genome version", c("hg38", "mm10")),
        selectInput("pam",
                  "PAM sequence", c("NGG", "NGRRT", "NGRRN", "TTTN")),
        selectInput("guide_mismatches",
                  "maximum # mismatches for the guide off-target search", c("1", "2", "3", "4", "5")),
        selectInput("guide_min_pvalue",
                  "p value filter for which off-targets to take forward to guard design", c("0.01", "0.05", "0.1")),
        selectInput("guard_length",
                    "GUARD RNA length", c("14", "15")),
        selectInput("guard_mismatches",
                    "maximum # mismatches for the GUARD RNA off-target search", c("1", "2", "3", "4")),
        sliderInput(inputId = "max_guard_distance", 
        label = "maximum distance between off-target and GUARD RNA binding",
        value = 10, min = 0, max = 15),
        actionButton(inputId = "submit", label = "submit job"),
        hr(),
        textInput("chr",
                  "Special case: If your guide maps perfectly to more than one genomic position, please specify the intended on-target locus: chr e.g. chr4"),
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
           h5("please cite: Coelho et al. Nature Communications, 2020"),
           a("GUARD Finder GitHub", href = "https://github.com/MatthewACoelho/GUARDfinder"),
        hr(),
           h5("License"),
           h6("Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the \"Software\"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software. 
THE SOFTWARE IS PROVIDED \"AS IS\", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
The software is not intended for clinical use.")
         )
    )
)


# Define server logic required to make results table
server <- function(input, output) {
    
    parameters <- eventReactive(input$submit, {
       params<- c("params",
          "{",
          "id =", input$id,
          "guide =", input$guide,
          "genome =", input$genome,
          "pam =", input$pam,
          "guide_mismatches =", input$guide_mismatches,
          "guide_min_pvalue =", input$guide_min_pvalue,
          "guard_length =", input$guard_length,
          "guard_mismatches =", input$guard_mismatches,
          "max_guard_distance =", input$max_guard_distance,
          "chr =", input$chr,
          "start =", input$start,
          "end =", input$end,
          "strand =", input$strand,
          "}")
       write(params, file="params.nf", sep="")
       print(params)
    })

    output$Results <- renderText({
        parameters()
    })
}

# Run the application 
shinyApp(ui = ui, server = server)
