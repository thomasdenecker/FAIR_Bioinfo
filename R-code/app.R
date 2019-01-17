################################################################################
# FAIR_app
# Thomas DENECKER
# 01 /2019
#
# GitHub :
# https://github.com/thomasdenecker/FAIR_Bioinfo
################################################################################

################################################################################
# Library
################################################################################

library(shiny)
library(shinydashboard)
library(DESeq2)
library(DT)
library(FactoMineR)
library(plotly)
library(reshape2)
library(shinyWidgets)
library(colourpicker)
library(shinyjs)
library(shinycssloaders)

################################################################################
# UI
################################################################################

header <- dashboardHeader(title = "FAIR_app")
sidebar <- dashboardSidebar(
  uiOutput('sidebar')

)

body <- dashboardBody(
  useShinyjs(),
  tags$head(tags$style(HTML('
                            
                            /* body */
                            .content-wrapper, .right-side {
                            background-color: rgba(255,255,255,1);
                            }
                            
                            .skin-red .main-header .navbar {
                                background-color: #222d32;
                            }

                          /* logo */
                                .skin-red .main-header .logo {
                                background-color: #222d32;
                                }

                                /* logo when hovered */
                                .skin-red .main-header .logo:hover {
                                background-color: #222d32;
                                }
                          button.myBtn {
                            color: white;
                            background: #ff3300;
                            border: 2px solid #cc2900;
                            border-radius: 12px;
                            font-size: 24px;
                            text-align: center;
                            margin-top: 20px;
                            margin-bottom: 40px;
                            margin-left: auto;
                            margin-right: auto;
                            display: block;
                            width: 200px;
                          }
                          
                          p.legend{
                            text-align: center;
                            font-style: italic;
                          }

                            '))),
  
  tabItems(
    #---------------------------------------------------------------------------
    # Home section
    #---------------------------------------------------------------------------
    tabItem(tabName = "home",
            h2("Welcome in FAIR_bioinfo application"),
            p("The aim is to find features that are differentially expressed between 2 conditions. The statistical 
              analysis process includes data normalization, graphical exploration of raw and normalized data, test for differential 
              expression for each feature between the conditions, raw p-value adjustment and export of lists of features having a s
              ignificant differential expression between the conditions.")
    ),
    
    
    tabItem(tabName = "import",
            h2("Data import"),
            h3("Conditions"),
            fluidRow(
              column(3,
                     h3("Parameters"),
                     fileInput("ConditionFile",label = NULL,
                               buttonLabel = "Browse...",
                               placeholder = "No file selected"),align = "center",
                     tags$hr(),
                     
                     # Input: Checkbox if file has header
                     radioButtons("header_condition", "Header",
                                  choices = c("Yes" = TRUE,
                                              "No" = FALSE),
                                  selected = TRUE, inline=T),
                     
                     # Input: Select separator ----
                     radioButtons("sep_condition", "Separator",
                                  choices = c(Comma = ",",
                                              Semicolon = ";",
                                              Tab = "\t"),
                                  selected = "\t", inline=T),
                     
                     # Input: Select quotes ----
                     radioButtons("quote_condition", "Quote",
                                  choices = c(None = "",
                                              "Double Quote" = '"',
                                              "Single Quote" = "'"),
                                  selected = "", inline=T)
              ), 
              column(9, 
                     h3("File preview"),
                     dataTableOutput(outputId = "contents_condition"))
            ),
            
            h3("Count table"),
            fluidRow(
              column(3,
                     h3("Parameters"),
                     fileInput("countTableFile",label = NULL,
                               buttonLabel = "Browse...",
                               placeholder = "No file selected"),align = "center",
                     tags$hr(),
                     
                     # Input: Checkbox if file has header
                     radioButtons("header_CT", "Header",
                                  choices = c("Yes" = TRUE,
                                              "No" = FALSE),
                                  selected = TRUE, inline=T),
                     
                     # Input: Select separator ----
                     radioButtons("sep_CT", "Separator",
                                  choices = c(Comma = ",",
                                              Semicolon = ";",
                                              Tab = "\t"),
                                  selected = "\t", inline=T),
                     
                     # Input: Select quotes ----
                     radioButtons("quote_CT", "Quote",
                                  choices = c(None = "",
                                              "Double Quote" = '"',
                                              "Single Quote" = "'"),
                                  selected = "", inline=T)
              ), 
              column(9, 
                     h3("File preview"),
                     dataTableOutput(outputId = "contents_countTable"))
            ),
            div(actionButton("run", "Run", class="myBtn"), align = "center")
            
    ),
    
    #---------------------------------------------------------------------------
    # Description section 
    #---------------------------------------------------------------------------
    tabItem(tabName = "description",
            h2("Description of raw data"),
            p("The objective of this application is to find the differentially expressed genes after using the FAIR_Bioinfo workflow"),
            h3("Conditions"),
            p("The count data files and associated biological conditions are listed in the following table : "),
            div(withSpinner(tableOutput('table')), align = "center"),
            p(class="legend", "Table 1: Data files and associated biological conditions."),
            
            h3("Count table"),
            textOutput("rawDataText"),
            tags$br(),
            withSpinner(DTOutput("countTableDT")),
            p(class="legend", "Table 2: View of the count data table."),
            
            p("Looking at the summary of the count table provides a basic description of these raw counts (min and max values, median, etc)."),
            div(tableOutput('summaryRawData'), align = "center"),
            p(class="legend", "Table 3: Summary of the raw counts."),
            
            h3("Total read count per sample"),
            p("Next figure shows the total number of mapped reads for each sample. Reads that map on multiple locations on the transcriptome are 
              counted more than once, as far as they are mapped on less than 50 different loci. We expect total read counts to be similar within conditions, 
              they may be different across conditions. Total counts sometimes vary widely between replicates. This may happen for several reasons, including:"),
            tags$ul(
              tags$li("different rRNA contamination levels between samples (even between biological replicates);"),
              tags$li("slight differences between library concentrations, since they may be difficult to measure with high precision.;")
            ),
            div(withSpinner(plotOutput("figure1", width ="50%")), align = "center"),
            p(class="legend", "Figure 1: Number of mapped reads per sample. Colors refer to the biological condition of the sample."),
            
            h3("Proportion of null counts sample"),
            textOutput("figure2Text"),
            div(withSpinner(plotOutput("figure2", width ="50%")), align = "center"),
            p(class="legend", "Figure 2: Proportion of features with null read counts in each sample."),
            
            h3("Density of count distribution"), 
            p("Next figure shows the distribution of read counts for each sample. For sake of readability, 
              log2(counts+1) are used instead of raw counts. Again we expect replicates to have similar 
              distributions. In addition, this figure shows if read counts are preferably low, medium or high. 
              This depends on the organisms as well as the biological conditions under consideration."),
            div(withSpinner(plotOutput("figure3", width ="50%")), align = "center"),
            p(class="legend", "Figure 3: Density distribution of read counts."),
            
            h3("Propotion of reads form most expressed sequence"),
            p("It may happen that one or a few features capture a high proportion of reads (up to 20% or more). 
              This phenomenon should not influence the normalization process. The DESeq2 normalization has proved 
              to be robust to this situation [Dillies, 2012]. Anyway, we expect these high count features to be 
              the same across replicates. They are not necessarily the same across conditions. Next Figure and next 
              table  illustrate the possible presence of such high count features in the data set."),
            div(withSpinner(plotOutput("figure4", width ="50%")), align = "center"),
            p(class="legend", "Figure 4: Percentage of reads associated with the sequence having the highest count 
              (provided in each box on the graph) for each sample."),
            div(withSpinner(tableOutput('table4')), align = "center"),
            p(class="legend", "Table 4: Percentage of reads associated with the sequences having the highest counts."),
            
            h3("Pairwise scatter plot"),
            p("We may wish to assess the similarity between samples across conditions. 
              A pairwise scatter plot is produced (figure 5) to show how replicates and samples from 
              different biological conditions are similar or different (log2(counts+1) are used instead of raw count values). 
              Moreover, as the Pearson correlation has been shown not to be relevant to measure the similarity between replicates, 
              the SERE statistic has been proposed as a similarity index between RNA-Seq samples [Schulze, 2012]. It measures whether 
              the variability between samples is random Poisson variability or higher. Pairwise SERE values are printed in the 
              lower triangle of the pairwise scatter plot. The value of the SERE statistic is:"),
            tags$ul(
              tags$li("0 when samples are identical (no variability at all: this may happen in the case of a sample duplication);"),
              tags$li("1 for technical replicates (technical variability follows a Poisson distribution);"),
              tags$li("greater than 1 for biological replicates and samples from different biological conditions (biological variability 
                      is higher than technical one, data are over-dispersed with respect to Poisson). The higher the SERE value, the lower 
                      the similarity. It is expected to be lower between biological replicates than between samples of different biological 
                      conditions. Hence, the SERE statistic can be used to detect inversions between samples.")
            ),
            
            div(withSpinner(plotOutput("figure5", width ="500px", height = "500px")), align = "center"),
            p(class="legend", "Figure 5: Pairwise comparison of samples.")
            
    ), 
    
    #---------------------------------------------------------------------------
    # data exploration
    #---------------------------------------------------------------------------
    tabItem(tabName = "dataExplo",
            h2("Variability within the experiment: data exploration"),
            p("The main variability within the experiment is expected to come from biological 
              differences between the samples. This can be checked in two ways. The first one is to 
              perform a hierarchical clustering of the whole sample set.[...]"),
            div(withSpinner(plotOutput("figure6", width ="50%")), align = "center"),
            p(class="legend", "Figure 6: Sample clustering based on normalized data."),
            
            p("Another way of visualizing the experiment variability is to look at the first principal 
              components of the PCA, as shown on the figure 7. On this figure, the first principal component (PC1) 
              is expected to separate samples from the different biological conditions, meaning that the biological 
              variability is the main source of variance in the data."),
            fluidRow(
              column(6,withSpinner(plotOutput("figure7_1")), align = "center"),
              column(6,withSpinner(plotOutput("figure7_2")), align = "center")
            ),
            p(class="legend", "Figure 7: First two components of a Principal Component Analysis, with percentages of 
              variance associated with each axis.")
    ),
    
    #---------------------------------------------------------------------------
    # Normalization
    #---------------------------------------------------------------------------
    tabItem(tabName = "normalization",
            h2("Normalization"),
            p("Normalization aims at correcting systematic technical biases in the data, in order to make read counts 
              comparable across samples. The normalization proposed by DESeq2 relies on the hypothesis that most 
              features are not differentially expressed. It computes a scaling factor for each sample. 
              Normalized read counts are obtained by dividing raw read counts by the scaling factor 
              associated with the sample they belong to. Scaling factors around 1 mean (almost) no normalization is 
              performed. Scaling factors lower than 1 will produce normalized counts higher than raw ones, and the other way around. 
              "),
            div(withSpinner(tableOutput('table5')), align = "center"),
            p(class="legend", "Table 5: Normalization factors."),
            
            p("The histograms (figure 8) can help to validate the choice of the normalization parameter (“median” or “shorth”).
              Under the hypothesis that most features are not differentially expressed, each size factor represented by a red line 
              is expected to be close to the mode of the distribution of the counts divided by their geometric means across samples."),
            div(withSpinner(plotOutput("figure8",  height = "800px")), align = "center"),
            p(class="legend", "Figure 8: Diagnostic of the estimation of the size factors."),
            
            p("The figure 9 shows that the scaling factors of DESeq2 and the total count normalization factors may not perform similarly."),
            div(withSpinner(plotOutput("figure9", width ="50%")), align = "center"),
            p(class="legend", "Figure 9: Plot of the estimated size factors and the total number of reads per sample."),
            
            p("Boxplots are often used as a qualitative measure of the quality of the normalization process, as they show how distributions are 
              globally affected during this process. We expect normalization to stabilize distributions across samples. Figure 10 shows boxplots 
              of raw (left) and normalized (right) data respectively."),
            fluidRow(
              column(6,withSpinner(plotlyOutput("figure10_1")), align = "center"),
              column(6,withSpinner(plotlyOutput("figure10_2")), align = "center")
            ),
            fluidRow(
              column(6,withSpinner(plotlyOutput("figure10_3")), align = "center"),
              column(6,withSpinner(plotlyOutput("figure10_4")), align = "center")
            ),
            p(class="legend", "Figure 10: Boxplots adn barplots of raw (left) and normalized (right) read counts.")
    ),
    
    
    #---------------------------------------------------------------------------
    # Differential analysis
    #---------------------------------------------------------------------------
    tabItem(tabName = "diffAna",
            h2(" Differential analysis"),
            h3("Modelisation"),
            p("DESeq2 aims at fitting one linear model per feature. For this project, the design used is counts ~ group and the goal is to 
               estimate the models' coefficients which can be interpreted as log2(FC). These coefficients will then be tested to get p-values 
               and adjusted p-values."),
            
            h3("Outlier detection"),
            p("Model outliers are features for which at least one sample seems unrelated to the experimental or study design. For every feature 
              and for every sample, the Cook's distance [Cook, 1977] reflects how the sample matches the model. A large value of the Cook's distance 
              indicates an outlier count and p-values are not computed for the corresponding feature."),
            
            h3("Dispersions estimation"),
            p("The DESeq2 model assumes that the count data follow a negative binomial distribution which is a robust alternative to the Poisson 
              law when data are over-dispersed (the variance is higher than the mean). The first step of the statistical procedure is to estimate
              the dispersion of the data. Its purpose is to determine the shape of the mean-variance relationship. The default is to apply a GLM 
              (Generalized Linear Model) based method (fitType=“parametric”), which can handle complex designs but may not converge in some cases. 
              The alternative is to use fitType=“local” as described in the original paper [Anders, 2010]. The parameter used for this project is 
              fitType=“parametric”. Then, DESeq2 imposes a Cox Reid-adjusted profile likelihood maximization [Cox, 1987 and McCarthy, 2012] and uses 
              the maximum a posteriori (MAP) of the dispersion [Wu, 2013]."),
            fluidRow(
              column(6,withSpinner(plotOutput("figure11_1")), align = "center"),
              column(6,withSpinner(plotOutput("figure11_2")), align = "center")
            ),
            p(class="legend", "Figure 11: Dispersion estimates (left) and diagnostic of log-normality (right)."),
            
            p("The left panel on figure 11 shows the result of the dispersion estimation step. The x- and y-axes represent the mean count value and
              the estimated dispersion respectively. Black dots represent empirical dispersion estimates for each feature (from the observed counts). 
              The red dots show the mean-variance relationship function (fitted dispersion value) as estimated by the model. The blue dots are the final
              estimates from the maximum a posteriori and are used to perform the statistical test. Blue circles (if any) point out dispersion outliers. 
              These are features with a very high empirical variance (computed from observed counts). These high dispersion values fall far from the
              model estimation. For these features, the statistical test is based on the empirical variance in order to be more conservative than with 
              the MAP dispersion. These features will have low chance to be declared significant. The figure on the right panel allows to check the 
              hypothesis of log-normality of the dispersions."),
            
            h3("Statistical test for differential expression"),
            p("Once the dispersion estimation and the model fitting have been done, DESeq2 can perform the statistical testing. 
              Figure 12 shows the distributions of raw p-values computed by the statistical test for the comparison(s) done. This 
              distribution is expected to be a mixture of a uniform distribution on [0,1] and a peak around 0 corresponding to the 
              differentially expressed features."),
            div(withSpinner(plotOutput("figure12", width ="50%")), align = "center"),
            p(class="legend", "Figure 12: Distribution(s) of raw p-values."),
            
            h3("Independent filtering"),
            p("DESeq2 can perform an independent filtering to increase the detection power of differentially expressed features at the 
              same experiment-wide type I error. Since features with very low counts are not likely to see significant differences typically 
              due to high dispersion, it defines a threshold on the mean of the normalized counts irrespective of the biological condition. 
              This procedure is independent because the information about the variables in the design formula is not used [Love, 2014].

              Table 6 reports the thresholds used for each comparison and the number of features discarded by the independent filtering. 
              Adjusted p-values of discarded features are then set to NA."),
            div(tableOutput('table6'), align = "center"),
            p(class="legend", "Table 6: Number of features discarded by the independent filtering for each comparison."),
            
            h3("Final results"),
            p("A p-value adjustment is performed to take into account multiple testing and control the false positive rate to a chosen level α.
              For this analysis, a BH p-value adjustment was performed [Benjamini, 1995 and 2001] and the level of controlled false positive rate 
              was set to 0.05."),
            
            withSpinner(DTOutput("DEseqTable")),
            p(class="legend", "Table 7: Complete table of Deseq2 results."),
            
            p("Figure 13 represents the MA-plot of the data for the comparisons done, where differentially expressed features are highlighted in red.
              A MA-plot represents the log ratio of differential expression as a function of the mean intensity for each feature. Triangles correspond 
              to features having a too low/high log2(FC) to be displayed on the plot."),
            div(withSpinner(plotOutput("figure13", width ="50%")), align = "center"),
            p(class="legend", "Figure 13: MA-plot(s) of each comparison. Red dots represent significantly differentially expressed features."),
            
            p("Figure 14 shows the volcano plots for the comparisons performed and differentially expressed features are still highlighted in red. 
              A volcano plot represents the log of the adjusted P value as a function of the log ratio of differential expression."),
            div(withSpinner(plotlyOutput("figure14", width ="50%")), align = "center"),
            p(class="legend", "Figure 14: Volcano plot(s) of each comparison. Red dots represent significantly differentially expressed features."),
            
            withSpinner(DTOutput("DEseqTableSelected")),
            p(class="legend", "Table 7: Selected genes.")
            
    ),
    
    #---------------------------------------------------------------------------
    # Parameters
    #---------------------------------------------------------------------------
    tabItem(tabName = "parameters",
            h2("Parameters"),
            p("You can change settings. The display will be refreshed automatically."),
            h3("Graphics"),
            colourpicker::colourInput("colHisto", "Histograms", "deepskyblue"), 
            colourpicker::colourInput("colCondA", "Condition A", "gold"), 
            colourpicker::colourInput("colCondB", "Condition B", "royalblue"),
            h3("Differential analysis"),
            
            fluidRow(
              column(3,h4("fitType"),selectInput("fitType",label = NULL, choices = c(parametric= "parametric", local="local", mean = "mean"), selected = "parametric")),
              column(9, h4("Documentation"), 
                     p('either "parametric", "local", or "mean" for the type of fitting of dispersions to the mean intensity.'),
                     tags$ul(
                       tags$li("parametric - fit a dispersion-mean relation of the form: dispersion = asymptDisp + extraPois / mean  
                               via a robust gamma-family GLM. The coefficients asymptDisp and extraPois are given in the attribute 
                               coefficients of the dispersionFunction of the object."), 
                       tags$li("local - use the locfit package to fit a local regression of log dispersions over log base mean 
                               (normal scale means and dispersions are input and output for dispersionFunction). The points are weighted by 
                               normalized mean count in the local regression."),
                       tags$li("mean - use the mean of gene-wise dispersion estimates.")
                     )
                     
              )
            ),
            
            fluidRow(
              column(3,h4("P-value"),numericInput('pvalue', NA, 0.05,
                                                  min = 0, max = 1, step = 0.001)),
              column(9, h4("Documentation"), 
                     p('The significance cutoff used for optimizing the independent filtering (by default 0.05). 
                       If the adjusted p-value cutoff (FDR) will be a value other than 0.05, alpha should be set to that value.')
              )
            ),

            fluidRow(
              column(3,h4("logFC theshold"),numericInput('logFC', NA, 2,
                                                         min = 0, max = 100, step = 0.01)),
              column(9, h4("Documentation"), 
                     p('logFC theshold')
              )
            ),
            
            fluidRow(
              column(3,h4("logFC theshold"),selectInput("pAdjustMethod",label = NULL, choices = c(holm = "holm", hochberg = "hochberg", hommel = "hommel", bonferroni = "bonferroni", 
                                                                                                  BH = "BH", BY = "BY", fdr = "fdr", none = "none"), 
                                                        selected = "BH")),
              column(9, h4("Documentation"), 
                     p("The adjustment methods include the Bonferroni correction ('bonferroni') in which the p-values are multiplied by the number of 
                                comparisons. Less conservative corrections are also included by Holm (1979) ('holm'), Hochberg (1988) ('hochberg'), Hommel (1988) ('hommel'), 
                                Benjamini & Hochberg (1995) ('BH' or its alias 'fdr'), and Benjamini & Yekutieli (2001) ('BY'), respectively. A pass-through option ('none') 
                                is also included. The set of methods are contained in the p.adjust.methods vector for the benefit of methods that need to have the method as an 
                                option and pass it on to p.adjust.

                                The first four methods are designed to give strong control of the family-wise error rate. There seems no reason to use the unmodified Bonferroni 
                                correction because it is dominated by Holm's method, which is also valid under arbitrary assumptions.
                                
                                Hochberg's and Hommel's methods are valid when the hypothesis tests are independent or when they are non-negatively associated (Sarkar, 1998; 
                                Sarkar and Chang, 1997). Hommel's method is more powerful than Hochberg's, but the difference is usually small and the Hochberg p-values are faster to compute.
                                
                                The 'BH' (aka 'fdr') and 'BY' method of Benjamini, Hochberg, and Yekutieli control the false discovery rate, the expected proportion of false discoveries 
                                amongst the rejected hypotheses. The false discovery rate is a less stringent condition than the family-wise error rate, so these methods are more powerful than the others.
                                
                                Note that you can set n larger than length(p) which means the unobserved p-values are assumed to be greater than all the observed p for 'bonferroni' and 'holm'
                                methods and equal to 1 for the other methods.")
              )
            )
            
            
    ), 
    
    #---------------------------------------------------------------------------
    # Session section
    #---------------------------------------------------------------------------
    tabItem(tabName = "session",
            h2("R session information and parameters"),
            p("The versions of the R software and Bioconductor packages used for this analysis are listed below. 
              It is important to save them if one wants to re-perform the analysis in the same conditions."),
            uiOutput("sessionText")
    ), 
    
    #---------------------------------------------------------------------------
    # Bibliography section
    #---------------------------------------------------------------------------
    tabItem(tabName = "biblio",
            h2("Bibliography"),
            tags$ol(
              tags$li(HTML("Anders S, Huber W. <b>Differential expression analysis for sequence count data</b>. <i>Genome Biology</i>. 2010; doi:10.1186/gb-2010-11-10-r106.")), 
              tags$li(HTML("Love M, Huber W, Anders S. <b>Moderated estimation of fold change and dispersion for RNA-Seq data with DESeq2</b>. <i>Genome Biology</i>. 2014; doi:10.1186/s13059-014-0550-8.")), 
              tags$li(HTML("Robinson M, McCarthy DJ, Smyth GK. <b>edgeR: a Bioconductor package for differential expression analysis of digital gene expression data</b>. <i>Bioinformatics</i>. 2009; doi:10.1093/bioinformatics/btp616.")),
              tags$li(HTML("Anders S, Pyl TP, Huber W. <b>HTSeq - A Python framework to work with high-throughput sequencing data</b>. <i>Bioinformatics</i>. 2014; doi:10.1093/bioinformatics/btu638.")),
              tags$li(HTML("Liao Y, Smyth GK and Shi W. <b>featureCounts: an efficient general purpose program for assigning sequence reads to genomic features</b>. <i>Bioinformatics</i>, 2014; doi:10.1093/bioinformatics/btt656.")),
              tags$li(HTML("Ritchie ME, Phipson B, Wu D, et al. <b>limma powers differential expression analyses for RNA-sequencing and microarray studies</b>. <i>Nucleic Acids Research</i>. 2015; doi:10.1093/nar/gkv007.")),
              tags$li(HTML("Cook RD. <b>Detection of Influential Observation in Linear Regression</b>. <i>Technometrics</i>. 1977; DOI:10.1080/00401706.2000.10485981.")),
              tags$li(HTML("Bourgon R, Gentleman R and Huber W. <b>Independent filtering increases detection power for high-throughput experiments</b>. <i>PNAS</i>. 2010; doi:10.1073/pnas.0914005107.")),
              tags$li(HTML("Benjamini Y and Hochberg Y. <b>Controlling the false discovery rate: a practical and powerful approach to multiple testing</b>. <i>Journal of the Royal Statistical Society B</i>. 1995; doi:10.2307/2346101.")),
              tags$li(HTML("Benjamini Y and Yekutieli D. <b>The control of the false discovery rate in multiple testing under dependency</b>. <i>Annals of Statistics</i>. 2001.")),
              tags$li(HTML("Schulze SK, Kanwar R, Golzenleuchter M, et al. <b>SERE: Single-parameter quality control and sample comparison for RNA-Seq</b>. <i>BMC Genomics</i>. 2012; doi:10.1186/1471-2164-13-524."))
            )
    )
  )
)


ui <- dashboardPage(skin= "red", header, sidebar, body)



server <- function(input, output, session) {
  
  si <- sessionInfo()
  
  #-----------------------------------------------------------------------------
  # Reactive Values
  #-----------------------------------------------------------------------------
  
  deseqRV <- reactiveValues()
  color <- reactiveValues()
  data <- reactiveValues()
  
  #-----------------------------------------------------------------------------
  # Side bar
  #-----------------------------------------------------------------------------
  
  data$run = F
  
  observeEvent(data$run, {
    if(data$run){
      output$sidebar <- renderUI({
        sidebarMenu( id = "tabs",
          menuItem("Home", tabName = "home", icon = icon("home")),
          menuItem("Import", tabName = "import", icon = icon("file-import")),
          menuItem("Description", tabName = "description", icon = icon("info"), selected = T),
          menuItem("Data exploration", tabName = "dataExplo", icon = icon("search")),
          menuItem("Normalization", tabName = "normalization", icon = icon("poll")),
          menuItem("Differential analysis", tabName = "diffAna", icon = icon("sliders-h")),
          menuItem("Parameters", tabName = "parameters", icon = icon("wrench")),
          menuItem("Session information", tabName = "session", icon = icon("cubes")),
          menuItem("Bibliography", tabName = "biblio", icon = icon("book"))
        )
      })
      
      updateTabItems(session, "tabs", selected = "description")
      shinyjs::runjs("window.scrollTo(0, 0)")
      
    } else {
      output$sidebar <- renderUI({
        sidebarMenu(id = "tabs",
          menuItem("Home", tabName = "home", icon = icon("home"), selected = T),
          menuItem("Import", tabName = "import", icon = icon("file-import"))
        )
      })
      updateTabItems(session, "tabs", selected = "home")
      shinyjs::runjs("window.scrollTo(0, 0)")
    }
  })
  
  #-----------------------------------------------------------------------------
  # Parameters - Colors
  #-----------------------------------------------------------------------------
  
  observeEvent(data$conds,{
    color$colConds = rep(input$colCondA, length(data$conds))
    color$colConds[which(data$conds  == "CondB")] = input$colCondB
  })
  
  
  observeEvent(input$colHisto,{
    color$hist = input$colHisto
  })
  
  observeEvent(input$colCondA,{
    color$colCondA = input$colCondA
    color$colConds = rep(input$colCondA, length(data$conds))
    color$colConds[which(data$conds  == "CondB")] = input$colCondB
  })
  
  observeEvent(input$colCondB,{
    color$colCondB = input$colCondB
    color$colConds = rep(input$colCondA, length(data$conds)) 
    color$colConds[which(data$conds  == "CondB")] = input$colCondB
  })
  
  
  #---------------------------------------------------------------------------
  # config read
  #---------------------------------------------------------------------------
  
  output$contents_condition <-  renderDataTable({
    
    req(input$ConditionFile)
    
    df <- read.csv(input$ConditionFile$datapath,
                   header = as.logical(input$header_condition),
                   sep = input$sep_condition,
                   quote = input$quote_condition,
                   nrows=10
    )
  },  options = list(scrollX = TRUE , dom = 't'))
  
  
  output$contents_countTable <-  renderDataTable({ 
    
    req(input$countTableFile)
    
    df <- read.csv(input$countTableFile$datapath,
                   header = as.logical(input$header_CT),
                   sep = input$sep_CT,
                   quote = input$quote_CT,
                   nrows=10
    )
  },  options = list(scrollX = TRUE , dom = 't'))
  
  #---------------------------------------------------------------------------
  # Run 
  #---------------------------------------------------------------------------
  
  observe({
    if(is.null(input$ConditionFile) || is.null(input$countTableFile) ){
      shinyjs::disable("run")
    } else{
      shinyjs::enable("run")
    }
  })
  
  
  observeEvent(input$run, {
    withProgress(message = 'Making plot', value = 0, {
      n <- 6
      #-------------------------------------------------------------------------
      # Read data
      #-------------------------------------------------------------------------
      
      incProgress(1/n, detail = "Read conditions")
      data$conditions <- read.csv2(input$ConditionFile$datapath,
                                   header = as.logical(input$header_condition),
                                   sep = input$sep_condition,
                                   quote = input$quote_condition
      )
      
      incProgress(1/n, detail = "Read count table")
      data$countTable <- read.csv2(input$countTableFile$datapath,
                                   header = as.logical(input$header_CT),
                                   sep = input$sep_CT,
                                   quote = input$quote_CT
      )
      
      incProgress(1/n, detail = "Summaryze")
      data$summary <- do.call(cbind, lapply(data$countTable, summary))
      
      #-------------------------------------------------------------------------
      # DEseq2
      #-------------------------------------------------------------------------
      incProgress(1/n, detail = "DEseq preparation")
      data$conds    = factor(unlist(lapply(strsplit(as.character(colnames(data$countTable)), "_"), function(l) l[[1]])))
      data$colData  = data.frame(condition = data$conds)
      deseqRV$ddsObjet = DESeqDataSetFromMatrix(countData = data$countTable, 
                                                colData   = data$colData, formula(~ condition))
      deseqRV$ddsObjet = estimateSizeFactors(deseqRV$ddsObjet)
      
      incProgress(1/n, detail = "Normalisation")
      deseqRV$SF = sizeFactors(deseqRV$ddsObjet)
      deseqRV$SF = matrix(sizeFactors(deseqRV$ddsObjet), 1, length(deseqRV$SF),byrow = T, dimnames =list("Size factor", names(deseqRV$SF)))
      deseqRV$normCountData = counts(deseqRV$ddsObjet, normalized = TRUE)
      
      incProgress(1/n, detail = "DEseq calculation")
      deseqRV$ddsEstim = DESeq(deseqRV$ddsObjet, fitType = input$fitType)
      deseqRV$resDESeq = results(deseqRV$ddsEstim, contrast = c("condition", "CondA", "CondB"), pAdjustMethod = input$pAdjustMethod )
    })
    
    sendSweetAlert(
      session = session,
      title = "Done !",
      text = "The analysis was successful !",
      type = "success"
    )
    
    data$run = T
  })
  
  
  #-----------------------------------------------------------------------------
  # Parameters - Deseq2
  #-----------------------------------------------------------------------------
  
  observeEvent(input$fitType,{
    if(!is.null(deseqRV$ddsObjet)){
      deseqRV$ddsEstim = DESeq(deseqRV$ddsObjet, fitType = input$fitType)
      deseqRV$resDESeq = results(deseqRV$ddsEstim, contrast = c("condition", "CondA", "CondB"), pAdjustMethod = input$pAdjustMethod )
    }
  })
  
  observeEvent(input$pAdjustMethod,{
    if(!is.null(deseqRV$ddsObjet)){
      deseqRV$resDESeq = results(deseqRV$ddsEstim, contrast = c("condition", "CondA", "CondB"), pAdjustMethod = input$pAdjustMethod )
    }
  })
  
  #-----------------------------------------------------------------------------
  # Session section
  #-----------------------------------------------------------------------------
  
  output$sessionText = renderUI({
    HTML(paste(si[[1]]$version.string,",", si[[1]]$platform, "<br>","<br>",
               "Locale : ", si[[3]], "<br>","<br>",
               "Attached base packages : ", paste(si[[5]], collapse = " , "),"<br>","<br>",
               "Other attached packages :", paste(unlist(lapply(si$otherPkgs, function(x){paste(x$Package, x$Version)})), collapse = " , "), "<br>","<br>",
               "Loaded via a namespace (and not attached) :" ,paste(unlist(lapply(si$loadedOnly, function(x){paste(x$Package, x$Version)})), collapse = " , ") 
    ))
  })
  
  #---------------------------------------------------------------------------
  # Description section
  #---------------------------------------------------------------------------
  
  output$table <- renderTable({
    data$conditions
  })
  
  
  output$rawDataText = renderText({
    paste("After loading the data we first have a look at the raw data table itself. The data table contains 
            one row per annotated feature and one column per sequenced sample. Row names of this table are feature IDs (unique identifiers). 
            The table contains raw count values representing the number of reads that map onto the features. For this project, there are", nrow(data$countTable) ,"features in the count data table.")
  })
  
  output$countTableDT <- renderDT(data$countTable, server = TRUE)
  
  output$summaryRawData <- renderTable(data$summary, rownames = T) 
  
  output$figure1 <- renderPlot({
    barplot(colSums(data$countTable), ylab = "Total read count per sample",
            main = "Total read count", col = color$colConds,
            names = colnames(data$countTable))
  })
  
  output$figure2Text = renderText({
    paste("Next figure shows the proportion of features with no read count in each sample. We expect this 
            proportion to be similar within conditions. Features with null read counts in the 4 samples are 
            left in the data but are not taken into account for the analysis with DESeq2. 
            Here, ",sum(apply(data$countTable, 1, function(x) all(x==0)))," features (",sum(apply(data$countTable, 1, function(x) all(x==0)))*100/nrow(data$countTable) ,"%) are in this situation (dashed line). 
            Results for those features (fold-change and p-values) are set to NA in the results files.")
  })
  
  output$figure2 <- renderPlot({
    barplot(apply(data$countTable, 2, function(c){sum(c==0)*100/length(c)}), ylab = "Proportion of null counts",
            main = "Proportion of null counts per sample", col = color$colConds,
            names = colnames(data$countTable))
    abline(h=sum(apply(data$countTable, 1, function(x) all(x==0)))*100/nrow(data$countTable), lty= 2)
  })
  
  output$figure3 <- renderPlot({
    plot(1, type="n", xlab="log2(raw count +1)", ylab="Density", ylim=c(0, 1), xlim=c( log(min(data$countTable)+1), log(max(data$countTable)+1)))
    for( i in 1:ncol(data$countTable)){
      if(length(grep("CondA", colnames(data$countTable)[i])) == 1){
        lines(density(log(data$countTable[,i]+1)), col= color$colCondA)
      } else {
        lines(density(log(data$countTable[,i]+1)), col = color$colCondB)
      }
    }
    legend("topright", legend = c("CondA", "CondB"), col = c(color$colCondA, color$colCondB), lty = 1, inset = 0.01, box.lty = 0)
    
  })
  output$figure4 <- renderPlot({
    barplot(apply(data$countTable, 2, function(c){max(c)*100/sum(c)}), ylab = "Proportion of reads",
            main = "Proportion of reads from most expressed sequence", col = color$colConds,
            names = colnames(data$countTable))
    
  })
  
  
  output$table4 <- renderTable(apply(data$countTable, 2, function(x){x*100/sum(x)})[apply(data$countTable, 2, function(c){names(c)[which.max(c)]}), ], rownames = T) 
  
  output$figure5 <- renderPlot({
    # plot(log(data$countTable+1), pch= 20)
    NULL
  })
  
  #---------------------------------------------------------------------------
  # data exploration
  #---------------------------------------------------------------------------
  output$figure6 <- renderPlot({
    plot(hclust(dist(t(data$countTable))), hang = -1)
  })
  
  output$figure7_1 <- renderPlot({
    PCA(data$countTable, graph = T,axes = c(1, 2))
  })
  
  output$figure7_2 <- renderPlot({
    PCA(data$countTable, graph = T,axes = c(1, 3))
  })
  
  
  #-----------------------------------------------------------------------------
  # Normalization
  #-----------------------------------------------------------------------------
  output$table5 <- renderTable(deseqRV$SF, rownames = T) 
  
  output$figure8 <- renderPlot({
    par(mfrow = c(2,3))
    for(i in 1:ncol(deseqRV$normCountData)){
      hist(log2(deseqRV$normCountData[,i]/mean(deseqRV$normCountData[,i],na.rm = T)), main = paste("Size factors diagnostic -", colnames(deseqRV$normCountData)[i]),
           xlab = "log2(counts/geometic mean)", col = color$hist)
      abline(v=deseqRV$SF[1,i], col = "red")
    }
  })
  
  output$figure9 <- renderPlot({
    plot(deseqRV$SF[1,], apply(deseqRV$normCountData, 2, sum), xlab = "Size factors", ylab="Total number reads")
    abline(0,1)
  })
  
  output$figure10_1 <- renderPlotly({
    inter = melt(as.matrix(log(data$countTable+1)))
    plot_ly(inter,  y = ~value, x = ~Var2 , color = unlist(lapply(strsplit(as.character( inter$Var2 ), "_"), function(l) l[[1]]))  , colors = unique(color$colConds),
            text = ~paste("Feature: ", Var1),
            type = "box") %>%
      layout(title = 'Raw counts distribution',
             yaxis = list(zeroline = FALSE, title= "log2(raw count +1)"),
             xaxis = list(zeroline = FALSE, title= ""))
  })
  
  output$figure10_2 <- renderPlotly({
    inter = melt(as.matrix(log(deseqRV$normCountData+1)))
    plot_ly(inter,  y = ~value, x = ~Var2 , color = unlist(lapply(strsplit(as.character( inter$Var2 ), "_"), function(l) l[[1]]))  , colors = unique(color$colConds),
            text = ~paste("Feature: ", Var1),
            type = "box") %>%
      layout(title = 'Normalized counts distribution',
             yaxis = list(zeroline = FALSE, title= "log2(norm count +1)"),
             xaxis = list(zeroline = FALSE, title= ""))
    
  })
  
  output$figure10_3 <- renderPlotly({
    inter = colSums(data$countTable)
    plot_ly(
      x = names(inter),
      y = inter,
      color = unlist(lapply(strsplit(as.character( names(inter) ), "_"), function(l) l[[1]]))  , 
      colors = unique(color$colConds),
      name = unlist(lapply(strsplit(as.character( names(inter) ), "_"), function(l) l[[1]])) ,
      type = "bar"
    ) %>%
      layout(title = 'Total read count - Raw',
             yaxis = list(zeroline = FALSE, title= "Total read count raw per sample"),
             xaxis = list(zeroline = FALSE, title= ""))
    
  })
  
  output$figure10_4 <- renderPlotly({
    
    inter = colSums(deseqRV$normCountData)
    plot_ly(
      x = names(inter),
      y = inter,
      color = unlist(lapply(strsplit(as.character( names(inter) ), "_"), function(l) l[[1]]))  , 
      colors = unique(color$colConds),
      name = unlist(lapply(strsplit(as.character( names(inter) ), "_"), function(l) l[[1]])) ,
      type = "bar"
    ) %>%
      layout(title = 'Total read count - Normalized',
             yaxis = list(zeroline = FALSE, title= "Total read normalized count per sample"),
             xaxis = list(zeroline = FALSE, title= ""))
  })
  
  #-----------------------------------------------------------------------------
  # Differantial analysis
  #-----------------------------------------------------------------------------
  
  output$figure11_1 <- renderPlot({
    plotDispEsts(deseqRV$ddsEstim)
  })
  
  output$figure11_2 <- renderPlot({
    inter = log(dispersions(deseqRV$ddsEstim))
    hist(inter, 
         col=color$hist, 
         border="black",
         prob = TRUE, 
         xlab = "Feature dispersion estimate",
         main = "log-normality dispersion diagnostic",
         ylim = c(0,1))
    lines(density(inter, na.rm = T), 
          lwd = 2, 
          col = "black")
  })
  
  output$figure12 <- renderPlot({
    hist(deseqRV$resDESeq$pvalue, breaks= 100, col = color$hist)
  })
  
  output$table6 <- renderTable(metadata(deseqRV$resDESeq)$filterNumRej[which(rownames(metadata(deseqRV$resDESeq)$filterNumRej) == names(metadata(deseqRV$resDESeq)$filterThreshold)),] , rownames = T) 
  
  
  output$DEseqTable <- renderDT(as.data.frame(deseqRV$resDESeq), 
                                server = FALSE,
                                editable = F,
                                extensions = 'Buttons',
                                options = list(scrollX = TRUE, 
                                               extensions = 'Buttons', 
                                               searchHighlight = TRUE,
                                               dom = 'Bfrtip',
                                               buttons = c('csv', 'excel','print')))
  
  output$figure13 <- renderPlot({
    plotMA(deseqRV$resDESeq)
  })
  
  output$figure14 <- renderPlotly({
    
    inter = cbind(x = deseqRV$resDESeq$log2FoldChange, y = -log10(deseqRV$resDESeq$padj), feature = rownames(deseqRV$resDESeq), SE = deseqRV$resDESeq$lfcSE)
    inter = na.omit(inter)
    inter = as.data.frame(inter)
    inter[,1] = as.numeric(as.character(inter[,1]))
    inter[,2] = as.numeric(as.character(inter[,2]))
    
    color = rep("black", nrow(inter))
    pos = which(abs(inter$x) >= input$logFC & inter$y >= -log10(input$pvalue))
    color[pos] = "red"
    
    deseqRV$selected = deseqRV$resDESeq[which(abs(deseqRV$resDESeq$log2FoldChange) >= input$logFC & deseqRV$resDESeq$padj <= input$pvalue),]
    
    plot_ly(inter, x = ~x, y = ~y, type = 'scatter', mode = 'markers',
            text = ~paste("Feature: ", feature, '<br>lfcSE:', SE),
            marker = list(color = color)) %>%
      layout(title = 'Volcano plot',
             shapes=list(list(type='line', x0=min(inter$x)-1, x1= max(inter$x)+1, y0=-log10(input$pvalue), y1=-log10(input$pvalue), line=list(dash='dot', width=1)),
                         list(type='line', x0=-input$logFC, x1= -input$logFC, y0=0, y1=max(inter$y), line=list(dash='dot', width=1)),
                         list(type='line', x0=input$logFC, x1= input$logFC, y0=0, y1=max(inter$y), line=list(dash='dot', width=1))),
             yaxis = list(zeroline = FALSE, title= "-log10(adjusted pvalue)"),
             xaxis = list(zeroline = FALSE, title= "log2(fold change)"))
    
  })
  
  
  output$DEseqTableSelected <- renderDT(as.data.frame(deseqRV$selected), 
                                server = FALSE,
                                editable = F,
                                extensions = 'Buttons',
                                options = list(scrollX = TRUE, 
                                               extensions = 'Buttons', 
                                               searchHighlight = TRUE,
                                               dom = 'Bfrtip',
                                               buttons = c('csv', 'excel','print')))
  
}

shinyApp(ui, server)