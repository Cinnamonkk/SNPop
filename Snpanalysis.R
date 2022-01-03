#Required libraries###################################################################################################################
library(shiny)
library(bslib)
library(tidyverse)
library(LDlinkR)
library(heatmaply)
library(shinycssloaders)
library(reactlog)
library(devtools)

#Funcs#################################################################################################################################

Plotmaker <- function(SNP,df) {
#Function to plot the major,minor allele frequencies in different populations##########################################################
  x <- data.frame(which(df == SNP,arr.ind =TRUE))
  rowcount <- x[['row']]

  valuesALL <-  c(df[rowcount,'Total.minor.allele.frequency'],df[rowcount,'Total.major.allele.frequency'])
  valuesAFR <-  c(df[rowcount,'African.minor.allele.frequency'],df[rowcount,'African.major.allele.frequency'])
  valuesAMR <-  c(df[rowcount,'American.minor.allele.frequency'],df[rowcount,'American.major.allele.frequency'])
  valuesEUR <-  c(df[rowcount,'European.minor.allele.frequency'],df[rowcount,'European.major.allele.frequency'])
  valuesSAS <-  c(df[rowcount,'South.Asian.minor.allele.frequency'],df[rowcount,'South.Asian.major.allele.frequency'])
  valuesEAS <-  c(df[rowcount,'East.Asian.minor.allele.frequency'],df[rowcount,'East.Asian.major.allele.frequency'])
  labels1 <- c(df[rowcount,'Minor.allele'],df[rowcount,'Major.allele'])

layout(matrix(c(1,2,3,4,5,6), 2, 3, byrow = TRUE))
   pie(valuesALL,paste(paste(round(valuesALL,2)*100,'%'),labels1),main = "Total Population Allele Frequencies",radius = 1.09)
   pie(valuesAFR,paste(paste(round(valuesAFR,2)*100,'%'),labels1),main = "African Population Allele Frequencies",radius = 1.09)
   pie(valuesAMR,paste(paste(round(valuesAMR,2)*100,'%'),labels1),main = "American Population Allele Frequencies",radius = 1.09)
   pie(valuesEUR,paste(paste(round(valuesEUR,2)*100,'%'),labels1),main = "European Population Allele Frequencies",radius = 1.09)
   pie(valuesSAS,paste(paste(round(valuesSAS,2)*100,'%'),labels1),main = "South Asian Population Allele Frequencies",radius = 1.09)
   pie(valuesEAS,paste(paste(round(valuesEAS,2)*100,'%'),labels1),main = "East Asian Population Allele Frequencies",radius = 1.09)
}

Plotmaker1 <- function(SNP,df){
  x <- data.frame(which(df == SNP,arr.ind =TRUE))
  #Function to plot the genotype frequencies in different populations##########################################################
  rowcount <- x[['row']]
  valuesALL <- c(df[rowcount,'Total.heterozygous'],df[rowcount,'Total.minor.allele.homozygous'],df[rowcount,'Total.major.allele.homozygous'])
  valuesAFR <- c(df[rowcount,'African.heterozygous'],df[rowcount,'African.minor.allele.homozygous'],df[rowcount,'African.major.allele.homozygous'])
  valuesAMR <- c(df[rowcount,'American.heterozygous'],df[rowcount,'American.minor.allele.homozygous'],df[rowcount,'American.major.allele.homozygous'])
  valuesEUR <- c(df[rowcount,'European.heterozygous'],df[rowcount,'European.minor.allele.homozygous'],df[rowcount,'European.major.allele.homozygous'])
  valuesSAS <- c(df[rowcount,'South.Asian.heterozygous'],df[rowcount,'South.Asian.minor.allele.homozygous'],df[rowcount,'South.Asian.major.allele.homozygous'])
  valuesEAS <- c(df[rowcount,'East.Asian.heterozygous'],df[rowcount,'East.Asian.minor.allele.homozygous'],df[rowcount,'East.Asian.major.allele.homozygous'])
  labels1 <- c( paste0(df[rowcount,'Minor.allele'],df[rowcount,'Major.allele']), paste0(df[rowcount,'Minor.allele'],df[rowcount,'Minor.allele']),paste0(df[rowcount,'Major.allele'],df[rowcount,'Major.allele']))

  layout(matrix(c(1,2,3,4,5,6), 2, 3, byrow = TRUE))
    pie(valuesALL,paste(paste(round(valuesALL,2)*100,'%'),labels1),main = "Total Population Genotype Frequencies",radius = 1.09)
    pie(valuesAFR,paste(paste(round(valuesAFR,2)*100,'%'),labels1),main = "African Population Genotype Frequencies",radius = 1.09)
    pie(valuesAMR,paste(paste(round(valuesAMR,2)*100,'%'),labels1),main = "American Population Genotype Frequencies",radius = 1.09)
    pie(valuesEUR,paste(paste(round(valuesEUR,2)*100,'%'),labels1),main = "European Population Genotype Frequencies",radius = 1.09)
    pie(valuesSAS,paste(paste(round(valuesSAS,2)*100,'%'),labels1),main = "South Asian Population Genotype Frequencies",radius = 1.09)
    pie(valuesEAS,paste(paste(round(valuesEAS,2)*100,'%'),labels1),main = "East Asian Population Genotype Frequencies",radius = 1.09)
}

traitablefunc <- function(haptable,row,df){

  y <- colnames(haptable)[0:(length(colnames(haptable)) -2)]

  z <- haptable[row,]

  v <- character()
  m <- character()
  
  l <- data.frame(z[0:(length(colnames(haptable)) -2)])

  for(snphap in y){
    dfrow <- which(df$SNP == snphap, arr.ind=TRUE)
    minoralleledf <- df[dfrow,][['Minor.allele']]
    majoralleledf <- df[dfrow,][['Major.allele']]
    allelehap <- z[[snphap]]

      if(allelehap == minoralleledf){
        w <- df[dfrow,][['Minor.allele.traits']]
        k <- df[dfrow,][['Clinical.Significance']]
        }
      if(allelehap == majoralleledf){
        w <- df[dfrow,][['Major.allele.traits']]
        k <- ''
        }
    v <- c(v,w)
    m <- c(m,k)
  }

  n <-v
  x <- data.frame(Traits = n,Clinical_Singificance = m,row.names = y)
  x1 <- as.data.frame(t(x))


  mylist <- list(one = x1 ,two = l)

  return(mylist)
}

hapfinderfunc <- function(snplist,population){

  snp <- sort(snplist)

  if (population == 'Total'){
    hapmap <- (LDhap(snp,pop ='ALL',token ='430a1f88f6ee'))
  }
  if (population == 'African'){
    hapmap <- (LDhap(snp,pop ='AFR',token ='430a1f88f6ee'))
  }
  if (population == 'American'){
    hapmap <- (LDhap(snp,pop ='AMR',token ='430a1f88f6ee'))
  }
  if (population == 'East Asian'){
    hapmap <- (LDhap(snp,pop ='EAS',token ='430a1f88f6ee'))
  }
  if (population == 'South Asian'){
    hapmap <- (LDhap(snp,pop ='SAS',token ='430a1f88f6ee'))
  }
  if (population == 'European'){
    hapmap <- (LDhap(snp,pop ='EUR',token ='430a1f88f6ee'))
  }

  return(hapmap)

}

r2finderfunc <- function(snp,population){
#Function to create an r2 matrix#######################################################################################################
  snp <- sort(snp)

  if (population == 'Total'){
    r2matrix <- (LDmatrix(snp,pop ='ALL',r2d = 'r2',token ='430a1f88f6ee'))
  }
  if (population == 'African'){
    r2matrix <- (LDmatrix(snp,pop ='AFR',r2d = 'r2',token ='430a1f88f6ee'))
  }
  if (population == 'American'){
    r2matrix <- (LDmatrix(snp,pop ='AMR',r2d = 'r2',token ='430a1f88f6ee'))
  }
  if (population == 'East Asian'){
    r2matrix <- (LDmatrix(snp,pop ='EAS',r2d = 'r2',token ='430a1f88f6ee'))
  }
  if (population == 'South Asian'){
    r2matrix <- (LDmatrix(snp,pop ='SAS',r2d = 'r2',token ='430a1f88f6ee'))
  }
  if (population == 'European'){
    r2matrix <- (LDmatrix(snp,pop ='EUR',r2d = 'r2',token ='430a1f88f6ee'))
  }


  return(r2matrix)
}

dfinderfunc <- function(snp,population){
  #Function to create an r2 matrix#######################################################################################################
  snp <- sort(snp)

  if (population == 'Total'){
    dmatrix <- (LDmatrix(snp,pop ='ALL',r2d = 'd',token ='430a1f88f6ee'))
  }
  if (population == 'African'){
    dmatrix <- (LDmatrix(snp,pop ='AFR',r2d = 'd',token ='430a1f88f6ee'))
  }
  if (population == 'American'){
    dmatrix <- (LDmatrix(snp,pop ='AMR',r2d = 'd',token ='430a1f88f6ee'))
  }
  if (population == 'East Asian'){
    dmatrix <- (LDmatrix(snp,pop ='EAS',r2d = 'd',token ='430a1f88f6ee'))
  }
  if (population == 'South Asian'){
    dmatrix <- (LDmatrix(snp,pop ='SAS',r2d = 'd',token ='430a1f88f6ee'))
  }
  if (population == 'European'){
    dmatrix <- (LDmatrix(snp,pop ='EUR',r2d = 'd',token ='430a1f88f6ee'))
  }

  return(dmatrix)
}

normheatmapfunc <- function(r2matrix){

  row.names(r2matrix) <- r2matrix$RS_number

  r2matrix1 <- r2matrix[,2:(length(r2matrix$RS_number)+1)]

  heatmatrix <- data.matrix(r2matrix1)

  heatmaply(heatmatrix, grid_gap = 1, distfun = dist_no_na ,na.color="black")
}


dist_no_na <- function(r2matrix1) {
#handles NA values in the heatmap
  edist <- dist(r2matrix1)
  edist[which(is.na(edist))] <- max(edist, na.rm=TRUE) * 1.1
  return(edist)
}


#Constants################################################################################################################################
thousandgenomes <- data.frame(c('ALL','AFR','AMR','EAS','EUR','SAS'),c(2504,661,347,504,503,489))


#Ui#####################################################################################################################################
ui <- fluidPage(
  theme = bs_theme(version = 4, bootswatch = "united"),
  titlePanel('SNPdata Visualization & Analysis'),
  tags$hr(),
  tabsetPanel(
    #First Panel Input#####################################################################################################################
    tabPanel('Input',
             sidebarLayout(
               wellPanel(
                 #SNP file input
                 fileInput("file1", "Input SNP data from SNPfinder.py:",
                           accept = '.csv'),
               ),
               mainPanel(
                 #SNP file output
                 DT::dataTableOutput("contents", width = 600),
               )
)),
    #Second Panel Input#####################################################################################################################
    tabPanel('Overview',
            sidebarLayout(
              wellPanel(
                #Creates the selection input
                uiOutput("SNPout"),
                tags$hr(),
                tags$b(HTML("SNP information:")),
                DT::dataTableOutput('SNPinfo',width = 300),
                tags$hr(),
                tags$b(HTML('Number of individuals per population:')),
                verbatimTextOutput('popcounts'),
                tags$hr()
              ),
              mainPanel(
                tags$i(HTML("<u>Minor allele frequencies for the SNP loaded in the total population:</u>")),
                withSpinner(plotOutput('rsidfreq',width = '1500px',height = '400px')),
                tabsetPanel(
                  tabPanel('Allele Frequencies',
                    withSpinner(plotOutput('popAlleleFreq',width = '1500px',height = '400px')),
                ),
                  tabPanel('Genotype Frequencies',
                    withSpinner(plotOutput('popGenFreq',width = '1500px',height = '400px'))
                         ),
                )
                )
)),
  #Third Panel Input#########################################################################################################################
  tabPanel('Linkage disequilibrium(LD) per Chromosome & Population',
           sidebarLayout(
             wellPanel(
              tags$div(HTML("Select a <b>chromosome</b> and a <b>population</b><br>to view pairwise Linkage Disequilibrium(LD)<br>data for the Loaded SNPs:")),
              tags$hr(),
              uiOutput("chromosomeout"),
              uiOutput("Populationout"),
              tags$hr(),
              tags$div(HTML("r2/D' Ranges between 0 and 1:</br>1 when two SNPs provide identical information,</br>0 when they are in perfect equilibrium.</br>- r2 takes into account the minor allele frequency.</br>- D' assumes that the minor allele is always present.")),
              tags$hr(),
              tags$b(HTML('Pairwaise r2 Matrix:')),
              withSpinner(DT::dataTableOutput('r2matrixtable', width = 300)),
              tags$hr(),
              tags$b(HTML("Pairwaise D' Matrix:")),
              withSpinner(DT::dataTableOutput('dmatrixtable', width = 300))
             ),
             mainPanel(
               tabsetPanel(
                tabPanel('Interactive LD heatmap:r2',
                  withSpinner(plotlyOutput('Normalheatmapr2',width = '100%',height = '1000px')),
                ),
                tabPanel("Interactive LD heatmap:D'",
                  withSpinner(plotlyOutput('NormalheatmapD',width = '100%',height = '1000px')),
                ),
                tabPanel('Full r2 pairwise matrix',
                  withSpinner(DT::dataTableOutput("r2matrixtablefull", width = 300)),
                ),
                tabPanel("Full D' pairwise matrix",
                         withSpinner(DT::dataTableOutput("dmatrixtablefull", width = 300)),
                ),
                )

)
)),
#Fourth Panel Input
  tabPanel('Haplotypes',
    sidebarLayout(
      wellPanel(
        tags$div(HTML("Select a <b>chromosome</b> , <b>2 - 30 SNPs</b> and a <b>population</b></br>to view the possible haplotypes</br>and their frequencies:")),
        tags$hr(),
        uiOutput("chromosomeout1"),
        tags$hr(),
        uiOutput('SNP1'),
        uiOutput('SNP2'),
        tags$hr(),
        uiOutput('Populations'),
        actionButton('Load','Submit'),
        actionButton('reset','Reset'),
        tags$hr(),
      ),
      mainPanel(
        tabsetPanel(
          tabPanel('Haplotype Frequencies',
            withSpinner(tableOutput('Haptable'))
          ),
          tabPanel('Haplotype to traits',
            uiOutput('Haplochoice'),
            tags$hr(),
            tags$b('Selected Haplotype:'),
            withSpinner(DT::dataTableOutput('Selectedhap')),
            tags$hr(),
            tags$b('Haplotype Traits:'),
            withSpinner(DT::dataTableOutput('Hapfreqtable')),

          ),
                )
                )

      )
)
)
)


#Server
server <- function(input, output) {
  #Panel handler#################################################################################


  #First Panel Output#########################################################################################################################
    SNPdata <- reactive({
      req(input$file1, file.exists(input$file1$datapath))
      read.csv(input$file1$datapath)
    })
    output$contents <-  DT::renderDataTable({
      req(SNPdata())
      SNPdata()
    })

  #Second Panel Output#######################################################################################################################
    output$SNPout<- renderUI({
      selectInput('SNPchoice',HTML('Select an SNP to view its<br>Allele and Genotype frequencies:'),choices = SNPdata()[['SNP']])
    })
    output$popcounts <- renderText({'Source: 1000 Genomes Project\nphase 3\nTotal: 2504\nAfrican: 661\nAmerican: 347\nEuropean: 503\nSouth Asian: 489\nEast Asian: 504'})
    choice1 <- reactive({input$SNPchoice
      })
    row <- reactive({which(SNPdata() == choice1(), arr.ind=TRUE)})
    output$SNPinfo <- DT::renderDataTable(SNPdata()[row()[1],18:22],options = list(scrollX = TRUE,  scrollY = "500px", scrollCollapse = TRUE,
                                                                    paging = FALSE ,searching = FALSE))

    output$popAlleleFreq <- renderPlot({Plotmaker(choice1(),SNPdata())
    })
    output$popGenFreq <- renderPlot({Plotmaker1(choice1(),SNPdata())
    })
    plot1 <-  reactive({ggplot(data = SNPdata())+ geom_col(mapping = aes(x = SNP, y = Total.minor.allele.frequency,alpha =Total.minor.allele.frequency,tooltip = SNP, data_id = SNP),fill = 'skyblue')})
    output$rsidfreq <-renderPlot({plot1() +theme(axis.text.x = element_text(angle = 90)) + labs(y ="In Total Population", x = 'SNP rsID')
    })

  #Third Panel Output########################################################################################################################
    chomosomes <- reactive({SNPdata()[['Chromosome']]})
    output$chromosomeout <- renderUI({
      selectInput('Chromchoice',label = NULL ,choices =chomosomes(),selected = NULL)
    })
    output$Populationout <- renderUI({
      selectInput('Popchoice',label = NULL, choices = c('Total','African','American','European','South Asian','East Asian'),selected = NULL)
    })

    choice2 <- reactive({input$Chromchoice
    })
    choice3 <- reactive({input$Popchoice
    })
    new <-reactive({which(SNPdata()$Chromosome == choice2(), arr.ind = TRUE)
    })
    newSNP <- reactive({SNPdata()[new(),][['SNP']]
    })
    newPOS <- reactive({SNPdata()[new(),][['Position']]
    })
    output$r2text <- renderText({'r2 : Ranges between 0 and 1\n1 when the two markers provide identical information.\n0 when they are in perfect equilibrium.'
    })
    r2matrix <- reactive({
      validate(need(length(newSNP()) >=2 ,'The Selected Chromosome has only one SNP in your Data,Please select another.'))
      r2finderfunc(newSNP(),choice3())
    })
    dmatrix <-reactive({
      validate(need(length(newSNP()) >=2 ,'The Selected Chromosome has only one SNP in your Data,Please select another.'))
      dfinderfunc(newSNP(),choice3())
    })
    output$r2matrixtable <- DT::renderDataTable(r2matrix(),options = list(scrollX = TRUE,  scrollY = "150px", scrollCollapse = TRUE,
                                                                       paging = FALSE ,searching = FALSE))
    output$r2matrixtablefull <- DT::renderDataTable(r2matrix(),options = list(paging = FALSE ,searching = FALSE))

    output$dmatrixtable <- DT::renderDataTable(dmatrix(),options = list(scrollX = TRUE,  scrollY = "150px", scrollCollapse = TRUE,
                                                                         paging = FALSE ,searching = FALSE))
    output$dmatrixtablefull <- DT::renderDataTable(dmatrix(),options = list(paging = FALSE ,searching = FALSE))

    output$NormalheatmapD <- renderPlotly({
    validate(need(length(newSNP()) >=2 ,'The Selected Chromosome has only one SNP in your Data,Please select another.'))
    normheatmapfunc(dmatrix())
    })
    output$Normalheatmapr2 <- renderPlotly({
      validate(need(length(newSNP()) >=2 ,'The Selected Chromosome has only one SNP in your Data,Please select another.'))
      normheatmapfunc(r2matrix())
    })
  #Fourth panel output#######################################################################################################
    output$chromosomeout1 <- renderUI({
      selectInput('Chromchoice1',label = NULL , choices =chomosomes(),selected = NULL)
    })
    output$SNP1 <- renderUI({
      selectInput('SNP1choice',label = NULL , choices = newSNP1(),selected = NULL,multiple = TRUE,selectize =TRUE,width = '400px')
    })
    output$Populations <- renderUI({
      selectInput('Popchoice1',label = NULL , choices = c('Total','African','American','European','South Asian','East Asian'),selected = NULL)
    })
    choice4 <- reactive({input$Chromchoice1})
    choice5 <- reactive({input$SNP1choice})
    choice7 <- reactive({input$Popchoice1})
    choice8 <- reactive({input$haplochoice})
    new1 <-reactive({which(SNPdata()$Chromosome == choice4(), arr.ind = TRUE)
    })
    newSNP1 <- reactive({SNPdata()[new1(),][['SNP']]
    })
    hapmap  <- reactive({hapfinderfunc(choice5(),choice7())
    })
    traits <- reactive({traitablefunc(hapmap(),traitrow(),SNPdata())
    })
    hapstring <- reactive({apply(hapmap(),1,paste,collapse= ' ')
    })
    traitrow <- reactive({grep(choice8(),hapstring())
    })
    output$Haptable <-  renderTable({
        validate(need(input$Load,''))
        validate(need(length(choice5()) < 30 , 'Input must be 1-30 SNPs'))
        hapmap()
    })
    output$Haplochoice <- renderUI({
        selectInput('haplochoice',label = 'Select a haplotype:', choices =hapstring())
    })
    output$Selectedhap <-  DT::renderDataTable(traits()$two ,options = list(paging = FALSE ,searching = FALSE))


    output$Hapfreqtable <- DT::renderDataTable(traits()$one ,options = list(paging = FALSE ,searching = FALSE))

    observeEvent(input$reset,{
      output$chromosomeout1 <- renderUI({
        selectInput('Chromchoice1',label = NULL , choices =chomosomes(),selected = NULL)
      })
      output$SNP1 <- renderUI({
        selectizeInput('SNP1choice',label = NULL , choices = newSNP1(),selected = NULL, multiple =TRUE,width = '400px')
      })
      output$Populations <- renderUI({
        selectInput('Popchoice1',label = NULL , choices = c('Total','African','American','European','South Asian','East Asian'),selected = NULL)
      })
    })
}

shinyApp(ui, server)


