# This is the user-interface definition of a Shiny web application.
# You can find out more about building applications with Shiny here:
#
# http://shiny.rstudio.com
#
source("global.R")


ui <- dashboardPage(
  title = "Coinciding Germline Somatic Spectra",
  
  dashboardHeader(title = "Computational analysis of germline and somatic mutation spectra",
                  titleWidth = 700),
  dashboardSidebar(disable = TRUE),
  dashboardBody(fluidRow(
    shinyjs::useShinyjs(),
    
    box(
      width = 3,
      HTML("<b>You are now analysing:</b>"),
      h4(textOutput("displayTitle")),
      
      #### Somatic panel options ####
      wellPanel(
        id = "somatic_panel",
        width = 3,
        
        selectInput(
          "somatic",
          label = h4("Somatic database"),
          choices = list("ICGC" = "icgc",
                         "COSMIC" = "cosmic"),
          selected = "cosmic"
        ),
        
        selectInput(
          "cancer",
          label = "Cancer type",
          choices = list(
            "Chronic lymphocytic leukemia" = "chronic_lymphocytic_leukemia",
            "Acute myeloid leukemia" = "acute_myeloid_leukemia",
            "Acute lymphoblastic b-cell leukemia" = "acute_lymphoblastic_b-cell_leukemia",
            "Acute lymphoblastic leukemia" = "acute_lymphoblastic_leukemia",
            "Lung cancer" = "lung_cancer",
            "Breast cancer" = "breast_cancer",
            "Malignant melanoma" = "malignant_melanoma",
            "Sarcoma" = "sarcoma",
            "Stomach cancer" = "stomach_cancer",
            "Prostate cancer" = "prostate_cancer",
            "Diffuse large B-cell lymphoma" = "diffuse_large_B_cell_lymphoma",
            "Glioma" = "glioma",
            "Colorectal cancer" = "colorectal_cancer",
            "Oesophageal cancer" = "oesophageal_cancer",
            "Ovarian cancer" = "ovarian_cancer",
            "Pancreatic cancer" = "pancreatic_cancer",
            "Liver cancer" = "liver_cancer",
            "Cervical cancer" = "cervical_cancer",
            "Cholangiocarcinoma" = "cholangiocarcinoma",
            "Urothelical cancer" = "urothelical_cancer",
            "Kidney cancer" = "kidney_cancer",
            "Any cancer" = "pancancer"
          ),
          
          selected = "pancancer"
        ),
        
        radioButtons(
          "cancerFreq",
          "Variant frequency across samples",
          c(
            "Recurrent" = "recurrent",
            "Nonrecurrent" = "nonrecurrent",
            "Both recurrent and nonrecurrent" = "ALL"
          ),
          selected = "ALL"
        )
      ),
      
      
      #### Germline panel options ####
      wellPanel(
        id = "germline_panel",
        selectInput(
          "germline",
          label = h4("Germline database"),
          choices = list(
            "1000Genomes" = "oneKG",
            "ExAC" = "ExAC",
            "dbSNP" = "dbSNP"
          ),
          selected = "ExAC"
        ),
        
        selectInput(
          "oneKGpopulation",
          label = "Populations",
          choices = list(
            "African" = "AFR_AF_1KG",
            "European" = "EUR_AF_1KG",
            "South Asian" = "SAS_AF_1KG",
            "American" = "AMR_AF_1KG",
            "East Asian" = "EAS_AF_1KG"
          ),
          selected = "EUR_AF_1KG"
        ),
        
        selectInput(
          "ExACpopulation",
          label = "Populations",
          choices = list(
            "African/African American" = "AF_AFR_EXAC",
            "American" = "AF_AMR_EXAC",
            "Non-Finnish European" = "AF_NFE_EXAC",
            "Finnish" = "AF_FIN_EXAC",
            "East Asian" = "AF_EAS_EXAC",
            "South Asian" = "AF_SAS_EXAC",
            "Adjusted Global" = "AF_Adj_EXAC"
          ),
          selected = "AF_Adj_EXAC"
        ),
        
        # Select input minor allele frequency
        selectInput(
          "MAF",
          label = "Minor allele frequency",
          choices = list(
            "Common: >= 5%" = "Common",
            "Low Frequency: 5% <--> 1%" = "LowFreq",
            "Rare: 1% <--> 0.1%" = "Rare",
            "Very Rare: < 0.1%" = "VeryRare",
            "Any frequency" = "ALL"
          ),
          selected = "ALL"
        ),
        
        radioButtons(
          "dbSNPvalidation",
          label = h4("dbSNP Validation"),
          choices = list(
            "byCluster" = "byCluster",
            "byHapMap" = "byHapMap",
            "byFrequency" = "byFrequency",
            "by2Hit2Allele" = "by2Hit2Allele",
            "by1000G" = "by1000G",
            "byOtherPop" = "byOtherPop",
            "suspect" = "suspect",
            "multiple" = "multiple",
            "unknown" = "unknown",
            "Any validation method" = "ALL"
          ),
          selected = "ALL"
        ),
        
        radioButtons(
          "dbSNPsubmissions",
          label = h4("dbSNP Submission frequency:"),
          list(
            "1" = "1",
            "2-5" = "2:5",
            "6-10" = "6:10",
            ">10" = ">10",
            "Any submisson" = "ALL"
          ),
          selected = "ALL"
        )
        
      ),
      
      
      radioButtons(
        "variantEffect",
        "Variant consequence",
        list(
          "Coding" = "coding",
          "Noncoding" = "noncoding",
          "Both coding and noncoding" = "ALL"
        ),
        
        selected = "ALL"
      )
    ),
    
    #### TABS ####
    tabBox(
      width = 9,
      id = "tabs",
      
      tabPanel(h4("Overall statistics"),
               value = "statistics",
               fluidRow(
                 column(12,
                        
                        fluidRow(
                          column(
                            6,
                            htmlOutput("somaticStats"),
                            htmlOutput("germlineStats"),
                            htmlOutput("coincidingStats")
                          ),
                          
                          
                          column(6,
                                 plotOutput(
                                   "vennDiagram", height = 400, width = 400
                                 ))
                        ))
               )),
      
      tabPanel(h4("Variant types"),
               value = "variants",
               
               fluidRow(
                 column(12,
                        fluidRow(
                          
                          column(6,
                                 tags$h3("Variant type plot"),
                                 plotOutput("variantPlot")),
                          column(
                            6,
                            tags$h3("Enrichment/Depletion of coinciding variants"),
                            plotOutput("variantPlot_logarithm")
                          )
                        )
                 )
                 # ,
                 # column(
                 #   12,
                 #   column(4,
                 #          plotOutput("coinciding_transition_transversion")),
                 # 
                 #   column(4,
                 #          plotOutput("somatic_transition_transversion")),
                 # 
                 #   column(4,
                 #          plotOutput("germline_transition_transversion"))
                 # )
                 
               )
      ),
      
      tabPanel(
        h4("Sequence context"),
        value = "context",
        
        wellPanel(
          id = "tPanel",
          style = "overflow-y:scroll; max-height: 1000px;background:white",
          tags$h3("Coinciding variants"),
          plotOutput("contextPlot1"),
          
          downloadButton('downloadContext', 'Download datatable'),
          
          tags$h3("Unique somatic variants"),
          plotOutput("contextPlot_somatic"),
          
          tags$h3("Unique germline variants"),
          plotOutput("contextPlot_germline"),
          
          tags$h3("Enrichment/Depletion of coinciding variants"),
          plotOutput("contextPlot_logarithm"),
          htmlOutput("ContextPlotwarning")
          
        )
      ),
      
      tabPanel(
        h4("Variant consequences"),
        value = "consequences",
        
        wellPanel(
          id = "tPanel",
          style = "overflow-y:scroll; max-height: 1000px;background:white",
          tags$h3("Coinciding variants"),
          plotOutput("consequencePlot"),
          downloadButton('downloadConsequence', 'Download datatable'),
          
          tags$h3("Unique somatic variants"),
          plotOutput("consequencePlot_somatic"),
          
          tags$h3("Unique germline variants"),
          plotOutput("consequencePlot_germline"),
          
          tags$h3("Enrichment/Depletion of coinciding variants"),
          plotOutput("consequencePlot_logarithm")
        )
      ),
      
      tabPanel(
        h4("Coinciding variants"),
        fluidPage(
          tags$span(style = "color:red", textOutput("dataTableMessage")),
          dataTableOutput("dataTable"),
          downloadButton('downloadVariantTable', 'Download datatable')
        )
      ),
      
      
      
      tabPanel(
        h4("Mutational signatures"),
        tags$h3("Signature plot"),
        plotOutput("signatureContextPlot"),
        fluidRow(column(
          12,
          column(6,
                 htmlOutput("signature2")),
          column(6,
                 plotOutput("signaturePie")),
          tableOutput('signatureTable')
        ))
      ),
      
      tabPanel(h4("About"),
               
               HTML("Developed by Øyvind Sørby. <br/>Supervised by Sigve Nakken.<br/>"),
               HTML("Built with <a href='https://shiny.rstudio.com' target='_blank'>Shiny Framework.</a>"),
               HTML("<br/> Mutational signatures identified by <a href='https://github.com/raerose01/deconstructSigs' target='_blank'>deconstructSigs.</a>"),
               HTML("<br/>Based on VCF-files preprocessed with <a href='http://www.ensembl.org/info/docs/tools/vep/index.html' target='_blank'>The Ensembl Variant Effect Predictor</a> and <a href='https://github.com/brentp/vcfanno' target='_blank'>vcfanno.</a>")
               
      ),
      
      selectInput(
        "unique",
        label = "Compare coinciding variants with:",
        choices = list("Unique Germline" = "unique_germline" ,
                       "Unique Somatic" = "unique_somatic"),
        selected = "unique_germline"
      )
    )
  ))
)

