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

################################################################################
# UI
################################################################################

header <- dashboardHeader(title = "FAIR_app")
sidebar <- dashboardSidebar(
  sidebarMenu(
    menuItem("Home", tabName = "home", icon = icon("home")),
    menuItem("Description", tabName = "description", icon = icon("info")),
    menuItem("Data exploration", tabName = "dataExplo", icon = icon("search")),
    menuItem("Normalization", tabName = "normalization", icon = icon("poll")),
    menuItem("Differential analysis", tabName = "diffAna", icon = icon("sliders-h")),
    menuItem("Session information", tabName = "session", icon = icon("cubes")),
    menuItem("Bibliography", tabName = "biblio", icon = icon("book"))
  )
)

body <- dashboardBody(
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

                            '))),
  
  tabItems(
    #---------------------------------------------------------------------------
    # Home section
    #---------------------------------------------------------------------------
    tabItem(tabName = "home",
            h2("Welcome in FAIR_bioinfo application"),
            
            p("The aim is to find features that are differentially expressed between 2 conditions. The statistical 
              analysis process includes data normalization, graphical exploration of raw and normalized data, test for differential 
              expression for each feature between the conditions, raw p-value adjustment and export of lists of features having a significant differential expression between the conditions.")
            
    ),
    
    #---------------------------------------------------------------------------
    # Description section 
    #---------------------------------------------------------------------------
    tabItem(tabName = "description",
            h2("Description of raw data"),
            p("The objective of this application is to find the differentially expressed genes after using the FAIR_Bioinfo workflow"),
            h3("Conditions"),
            p("The count data files and associated biological conditions are listed in the following table : "),
            div(tableOutput('table'), align = "center"),
            h3("Count table"),
            textOutput("rawDataText"),
            tags$br(),
            DTOutput("countTableDT"),
            p("Looking at the summary of the count table provides a basic description of these raw counts (min and max values, median, etc)."),
            div(tableOutput('summaryRawData'), align = "center"),
            
            h3("Total read count per sample"),
            p("Next figure shows the total number of mapped reads for each sample. Reads that map on multiple locations on the transcriptome are 
              counted more than once, as far as they are mapped on less than 50 different loci. We expect total read counts to be similar within conditions, 
              they may be different across conditions. Total counts sometimes vary widely between replicates. This may happen for several reasons, including:"),
            tags$ul(
              tags$li("different rRNA contamination levels between samples (even between biological replicates);"),
              tags$li("slight differences between library concentrations, since they may be difficult to measure with high precision.;")
            ),
            div(plotOutput("figure1", width ="50%"), align = "center"),
            
            h3("Proportion of null counts sample"),
            textOutput("figure2Text"),
            div(plotOutput("figure2", width ="50%"), align = "center"),
            
            h3("Density of count distribution"), 
            p("Next figure shows the distribution of read counts for each sample. For sake of readability, 
              log2(counts+1) are used instead of raw counts. Again we expect replicates to have similar 
              distributions. In addition, this figure shows if read counts are preferably low, medium or high. 
              This depends on the organisms as well as the biological conditions under consideration."),
            div(plotOutput("figure3", width ="50%"), align = "center"),
            
            h3("Propotion of reads form most expressed sequence"),
            p("It may happen that one or a few features capture a high proportion of reads (up to 20% or more). 
              This phenomenon should not influence the normalization process. The DESeq2 normalization has proved 
              to be robust to this situation [Dillies, 2012]. Anyway, we expect these high count features to be 
              the same across replicates. They are not necessarily the same across conditions. Next Figure and next 
              table  illustrate the possible presence of such high count features in the data set."),
            div(plotOutput("figure4", width ="50%"), align = "center"),
            div(tableOutput('table4'), align = "center"),
            
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
            
            div(plotOutput("figure5", width ="500px", height = "500px"), align = "center")
            
    ), 
    
    #---------------------------------------------------------------------------
    # data exploration
    #---------------------------------------------------------------------------
    tabItem(tabName = "dataExplo",
            h2("Variability within the experiment: data exploration"),
            p("The main variability within the experiment is expected to come from biological 
              differences between the samples. This can be checked in two ways. The first one is to 
              perform a hierarchical clustering of the whole sample set.[...]"),
            div(plotOutput("figure6", width ="50%"), align = "center"),
            
            p("Another way of visualizing the experiment variability is to look at the first principal 
              components of the PCA, as shown on the figure 7. On this figure, the first principal component (PC1) 
              is expected to separate samples from the different biological conditions, meaning that the biological 
              variability is the main source of variance in the data."),
            fluidRow(
              column(6,plotOutput("figure7_1"), align = "center"),
              column(6,plotOutput("figure7_2"), align = "center")
            )
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
            div(tableOutput('table5'), align = "center"),
            
            p("The histograms (figure 8) can help to validate the choice of the normalization parameter (“median” or “shorth”).
              Under the hypothesis that most features are not differentially expressed, each size factor represented by a red line 
              is expected to be close to the mode of the distribution of the counts divided by their geometric means across samples."),
            div(plotOutput("figure8",  height = "800px"), align = "center"),
            
            p("The figure 9 shows that the scaling factors of DESeq2 and the total count normalization factors may not perform similarly."),
            div(plotOutput("figure9", width ="50%"), align = "center"),
            
            p("Boxplots are often used as a qualitative measure of the quality of the normalization process, as they show how distributions are 
              globally affected during this process. We expect normalization to stabilize distributions across samples. Figure 10 shows boxplots 
              of raw (left) and normalized (right) data respectively."),
            fluidRow(
              column(6,plotOutput("figure10_1"), align = "center"),
              column(6,plotOutput("figure10_2"), align = "center")
            ),
            fluidRow(
              column(6,plotOutput("figure10_3"), align = "center"),
              column(6,plotOutput("figure10_4"), align = "center")
            )
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
              column(6,plotOutput("figure11_1"), align = "center"),
              column(6,plotOutput("figure11_2"), align = "center")
            ),
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
            div(plotOutput("figure12", width ="50%"), align = "center"),
            
            h3("Independent filtering"),
            p("DESeq2 can perform an independent filtering to increase the detection power of differentially expressed features at the 
              same experiment-wide type I error. Since features with very low counts are not likely to see significant differences typically 
              due to high dispersion, it defines a threshold on the mean of the normalized counts irrespective of the biological condition. 
              This procedure is independent because the information about the variables in the design formula is not used [Love, 2014].

              Table 6 reports the thresholds used for each comparison and the number of features discarded by the independent filtering. 
              Adjusted p-values of discarded features are then set to NA."),
            div(tableOutput('table6'), align = "center"),
            
            h3("Final results"),
            p("A p-value adjustment is performed to take into account multiple testing and control the false positive rate to a chosen level α.
              For this analysis, a BH p-value adjustment was performed [Benjamini, 1995 and 2001] and the level of controlled false positive rate 
              was set to 0.05."),
            
            p("Figure 13 represents the MA-plot of the data for the comparisons done, where differentially expressed features are highlighted in red.
              A MA-plot represents the log ratio of differential expression as a function of the mean intensity for each feature. Triangles correspond 
              to features having a too low/high log2(FC) to be displayed on the plot."),
            div(plotOutput("figure13", width ="50%"), align = "center"),
            
            p("Figure 14 shows the volcano plots for the comparisons performed and differentially expressed features are still highlighted in red. 
              A volcano plot represents the log of the adjusted P value as a function of the log ratio of differential expression."),
            div(plotlyOutput("figure14", width ="50%"), align = "center")
            
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
  #-----------------------------------------------------------------------------
  # Reactive Values
  #-----------------------------------------------------------------------------
  
  rv <- reactiveValues()
  
  #-----------------------------------------------------------------------------
  # Read data
  #-----------------------------------------------------------------------------
  si <- sessionInfo()
  
  conditions <- read.csv2("/home/rstudio/conditions.txt", sep ="\t", header = T)
  countTable <-  read.csv2("/home/rstudio/Project/countTable.txt", sep ="\t", header = T)
  summary <- do.call(cbind, lapply(countTable, summary))
  
  #-----------------------------------------------------------------------------
  # DEseq2
  #-----------------------------------------------------------------------------
  conds    = factor(unlist(lapply(strsplit(as.character(colnames(countTable)), "_"), function(l) l[[1]])))
  colConds = rep("gold", length(conds))
  colConds[which(conds == "CondB")] = "royalblue"
  
  colData  = data.frame(condition = conds)
  
  ddsObjet = DESeqDataSetFromMatrix(countData = countTable, 
                                    colData   = colData, formula(~ condition))
  ddsObjet = estimateSizeFactors(ddsObjet)
  SF = sizeFactors(ddsObjet)
  SF = matrix(sizeFactors(ddsObjet), 1, length(SF),byrow = T, dimnames =list("Size factor", names(SF)))
  normCountData = counts(ddsObjet, normalized = TRUE)
  
  ddsEstim = DESeq(ddsObjet)
  resDESeq = results(ddsEstim, contrast = c("condition", "CondA", "CondB"))
  
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
    conditions
  })
  
  
  output$rawDataText = renderText({
    paste("After loading the data we first have a look at the raw data table itself. The data table contains 
      one row per annotated feature and one column per sequenced sample. Row names of this table are feature IDs (unique identifiers). 
      The table contains raw count values representing the number of reads that map onto the features. For this project, there are", nrow(countTable) ,"features in the count data table.")
  })
  
  output$countTableDT <- renderDT(countTable, server = TRUE)
  
  output$summaryRawData <- renderTable(summary, rownames = T) 
  
  output$figure1 <- renderPlot({
    barplot(colSums(countTable), ylab = "Total read count per sample",
            main = "Total read count", col = c(rep("gold", 3), rep("royalblue", 3)),
            names = colnames(countTable))
  })
  
  output$figure2Text = renderText({
    paste("Next figure shows the proportion of features with no read count in each sample. We expect this 
          proportion to be similar within conditions. Features with null read counts in the 4 samples are 
          left in the data but are not taken into account for the analysis with DESeq2. 
          Here, ",sum(apply(countTable, 1, function(x) all(x==0)))," features (",sum(apply(countTable, 1, function(x) all(x==0)))*100/nrow(countTable) ,"%) are in this situation (dashed line). 
          Results for those features (fold-change and p-values) are set to NA in the results files.")
  })
  
  output$figure2 <- renderPlot({
    barplot(apply(countTable, 2, function(c){sum(c==0)*100/length(c)}), ylab = "Proportion of null counts",
            main = "Proportion of null counts per sample", col = c(rep("gold", 3), rep("royalblue", 3)),
            names = colnames(countTable))
    abline(h=sum(apply(countTable, 1, function(x) all(x==0)))*100/nrow(countTable), lty= 2)
  })
  
  output$figure3 <- renderPlot({
    plot(1, type="n", xlab="log2(raw count +1)", ylab="Density", ylim=c(0, 1), xlim=c( log(min(countTable)+1), log(max(countTable)+1)))
    for( i in 1:ncol(countTable)){
      if(length(grep("CondA", colnames(countTable)[i])) == 1){
        lines(density(log(countTable[,i]+1)), col= "gold")
      } else {
        lines(density(log(countTable[,i]+1)), col = "royalblue")
      }
    }
    legend("topright", legend = c("CondA", "CondB"), col = c("gold", "royalblue"), lty = 1, inset = 0.01, box.lty = 0)
    
  })
  output$figure4 <- renderPlot({
    barplot(apply(countTable, 2, function(c){max(c)*100/sum(c)}), ylab = "Proportion of reads",
            main = "Proportion of reads from most expressed sequence", col = c(rep("gold", 3), rep("royalblue", 3)),
            names = colnames(countTable))
    
  })
  
  
  output$table4 <- renderTable(apply(countTable, 2, function(x){x*100/sum(x)})[apply(countTable, 2, function(c){names(c)[which.max(c)]}), ], rownames = T) 
  
  output$figure5 <- renderPlot({
    # plot(log(countTable+1), pch= 20)
    NULL
  })
  
  #---------------------------------------------------------------------------
  # data exploration
  #---------------------------------------------------------------------------
  output$figure6 <- renderPlot({
    plot(hclust(dist(t(countTable))), hang = -1)
  })
  
  output$figure7_1 <- renderPlot({
    PCA(countTable, graph = T,axes = c(1, 2))
  })
  
  output$figure7_2 <- renderPlot({
    PCA(countTable, graph = T,axes = c(1, 3))
  })
  
  
  #-----------------------------------------------------------------------------
  # Normalization
  #-----------------------------------------------------------------------------
  output$table5 <- renderTable(SF, rownames = T) 
  
  output$figure8 <- renderPlot({
    par(mfrow = c(2,3))
    for(i in 1:ncol(normCountData)){
      hist(log2(normCountData[,i]/mean(normCountData[,i],na.rm = T)), main = paste("Size factors diagnostic -", colnames(normCountData)[i]),
           xlab = "log2(counts/geometic mean)", col = "cyan")
      abline(v=SF[1,i], col = "red")
    }
  })
  
  output$figure9 <- renderPlot({
    plot(SF[1,], apply(normCountData, 2, sum), xlab = "Size factors", ylab="Total number reads")
    abline(0,1)
  })
  
  output$figure10_1 <- renderPlot({
    boxplot(log(countTable+1), ylab = "log2(raw count +1)", col = colConds,
            main = "Raw counts distribution")
  })
  
  output$figure10_2 <- renderPlot({
    boxplot(log(normCountData+1), ylab = "log2(norm count +1)", col = colConds,
            main = "Normalized counts distribution")
  })
  
  output$figure10_3 <- renderPlot({
    barplot(colSums(countTable), ylab = "Total read count raw per sample",
            main = "Total read count", col = colConds,
            names = colnames(countTable))
  })
  
  output$figure10_4 <- renderPlot({
    barplot(colSums(normCountData), ylab = "Total read normalized count per sample",
            main = "Total read count", col = colConds,
            names = colnames(normCountData))
  })
  
  #-----------------------------------------------------------------------------
  # Differantial analysis
  #-----------------------------------------------------------------------------
  
  output$figure11_1 <- renderPlot({
    plotDispEsts(ddsEstim)
  })
  
  output$figure11_2 <- renderPlot({
    inter = log(dispersions(ddsEstim))
    hist(inter, 
         col="cyan", 
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
    hist(resDESeq$pvalue, breaks= 100, col = "cyan")
  })
  
  output$table6 <- renderTable(metadata(resDESeq)$filterNumRej[which(rownames(metadata(resDESeq)$filterNumRej) == names(metadata(resDESeq)$filterThreshold)),] , rownames = T) 
  
  output$figure13 <- renderPlot({
    plotMA(resDESeq)
  })
  
  output$figure14 <- renderPlotly({
    
    inter = cbind(x = resDESeq$log2FoldChange, y = -log10(resDESeq$padj), feature = rownames(resDESeq), SE = resDESeq$lfcSE)
    inter = na.omit(inter)
    inter = as.data.frame(inter)
    inter[,1] = as.numeric(as.character(inter[,1]))
    inter[,2] = as.numeric(as.character(inter[,2]))
    
    
    plot_ly(inter, x = ~x, y = ~y, type = 'scatter', mode = 'markers',
            text = ~paste("Feature: ", feature, '<br>lfcSE:', SE),
            color = "blue") %>%
        layout(title = 'Volcano plot',
               yaxis = list(zeroline = FALSE, title= "-log10(adjusted pvalue)"),
               xaxis = list(zeroline = FALSE, title= "log2(fold change)"))
    
  })

  
}

shinyApp(ui, server)