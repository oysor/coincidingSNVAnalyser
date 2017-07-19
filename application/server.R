
# This is the server logic for a Shiny web application.
# You can find out more about building applications with Shiny here:
#
# http://shiny.rstudio.com
#

source("global.R")

server <- function(input, output) {
  
  output$displayTitle <- renderText({ 
    
    somatic <- switch(input$somatic,
                      "cosmic" = "The Catalogue Of Somatic Mutations In Cancer",
                      "icgc" = "The International Cancer Genome Consortium")
    
    germline <- switch(input$germline,
                       "oneKG" = "The 1000 Genomes Project",
                       "dbSNP" = "The Database of Short Genetic Variation",
                       "ExAC" = "The Exome Aggregation Consortium")
    
    paste(somatic," & ", germline, sep = " ")
  })
  
  #### shinyjs ####
  observe({
    
    if(input$tabs == "variants" | input$tabs == "context" | input$tabs == "consequences"){
      shinyjs::show("unique")
    }else{
      shinyjs::hide("unique")
    }
    
    if(input$somatic == 'cosmic'){
      shinyjs::show("cancer")
      shinyjs::show("cancerFreq")
      shinyjs::show("dataTable")
    }else{
      shinyjs::reset("cancer")
      shinyjs::reset("cancerFreq")
      shinyjs::hide("cancer")
      shinyjs::hide("cancerFreq")
      shinyjs::hide("dataTable")
    }
    
    if (input$germline == "oneKG" | input$germline == "ExAC") {
      shinyjs::show("MAF")
      shinyjs::hide("dbSNPvalidation")
      shinyjs::hide("dbSNPsubmissions")
      shinyjs::reset("dbSNPvalidation")
      shinyjs::reset("dbSNPsubmissions")
      
      if(input$germline == "oneKG"){
        shinyjs::show("oneKGpopulation")
        shinyjs::reset("ExACpopulation")
        shinyjs::hide("ExACpopulation")
        
      }else if(input$germline == "ExAC"){
        shinyjs::show("ExACpopulation")
        shinyjs::reset("oneKGpopulation")
        shinyjs::hide("oneKGpopulation")
      }
    } else {
      shinyjs::reset("ExACpopulation")
      shinyjs::reset("oneKGpopulation")
      shinyjs::reset("MAF")
      shinyjs::hide("ExACpopulation")
      shinyjs::hide("oneKGpopulation")
      shinyjs::hide("MAF")
      shinyjs::show("dbSNPvalidation")
      shinyjs::show("dbSNPsubmissions")
    }
  })
  
  #### Database INPUT ####
  # Updates every time a somatic or germline database is changed
  # Updates every time a variant consequence is changed
  databaseInput <-  reactive({
    
    #print("1 DATABASE INPUT")
    dataset <- databases[[paste(input$somatic, input$germline, sep="_")]]
    coinciding <- dataset$coinciding
    unique_germline <- dataset$germline
    unique_somatic <- dataset$somatic
    
    if(input$variantConsequence != 'ALL' ){
      coinciding <- coinciding[coinciding$coding == input$variantConsequence,]
      unique_germline <- unique_germline[unique_germline$coding == input$variantConsequence,]
      unique_somatic  <- unique_somatic[unique_somatic$coding == input$variantConsequence,]
    }
    
    dataframes <- list(
      "coinciding" = coinciding,
      "unique_germline" = unique_germline,
      "unique_somatic" = unique_somatic
    )
  })
  
  #### Coinciding INPUT ####
  # Updates every time databaseInput() is updated
  # Updates every time a germline setting or somatic setting is changed
  coincidingInput <- reactive({
    
    #  print("2 COINCIDING INPUT")
    coinciding <- databaseInput()$coinciding
    
    if(input$somatic == "cosmic"){
      variantTable <- switch(input$germline,
                             "oneKG" = variantTables$cosmic_oneKG,
                             "ExAC" = variantTables$cosmic_ExAC,
                             "dbSNP" = variantTables$cosmic_dbSNP
      )
      variantTable_DL <- switch(input$germline,
                                "oneKG" = variantTables$cosmic_oneKG_DL,
                                "ExAC" = variantTables$cosmic_ExAC_DL,
                                "dbSNP" = variantTables$cosmic_dbSNP_DL
      )
    }else {
      variantTable = NULL 
      variantTable_DL = NULL
    }
    
    if(input$variantConsequence != 'ALL' ){
      index <- with(variantTable, grepl(paste0("= ",input$variantConsequence), variantTable$info, fixed = TRUE))
      variantTable <- variantTable[index]
      index_dl <- with(variantTable_DL, grepl(paste0("=",input$variantConsequence), variantTable_DL$info, fixed = TRUE))
      variantTable_DL <- variantTable_DL[index_dl]
    }
    
    if(input$germline == "oneKG" | input$germline == "ExAC"){
      
      population <- switch (input$germline,
                            "oneKG" = input$oneKGpopulation,
                            "ExAC" = input$ExACpopulation
      )
      
      coinciding <- coinciding[coinciding$region == population,]
      
      ## ensures filtering of frequency for the right population
      if (input$MAF != "ALL"){
        coinciding <- coinciding[coinciding$alleleFreq == input$MAF,]
        population <- paste0(population,":[-+]?[0-9]*\\.?[0-9]+([eE][-+]?[0-9]+)?\\(",input$MAF,"\\)")
      }
      
      index <- with(variantTable, grepl(population , variantTable$info))
      variantTable <- variantTable[index]
      variantTable_DL <- variantTable_DL[index]
      # dbSNP  
    } else {
      
      if(input$dbSNPvalidation != "ALL"){
        
        coinciding <- coinciding[coinciding$validation == input$dbSNPvalidation,]
        index <- with(variantTable, grepl(paste0("\\(",input$dbSNPvalidation,"\\)"), variantTable$info))
        variantTable <- variantTable[index]
        variantTable_DL <- variantTable_DL[index]
        
      }
      
      if(input$dbSNPsubmissions != "ALL"){
        
        coinciding <- coinciding[coinciding$submission == input$dbSNPsubmissions,]
        
        index <- with(variantTable, grepl(paste0("\\(",input$dbSNPsubmissions,"\\)"), variantTable$info))
        variantTable <- variantTable[index]
        variantTable_DL <- variantTable_DL[index]
        
      }
    }
    
    
    
    if(input$somatic == "cosmic"){
      
      coinciding <- coinciding[coinciding$cancer == input$cancer,]
      cancerType <- input$cancer
      
      if(input$cancerFreq != "ALL"){
        coinciding <- coinciding[coinciding$cancerFreq == input$cancerFreq,]
        cancerType <-paste0(cancerType,":[0-9]+\\(",input$cancerFreq,"\\)")
      }
      index <- with(variantTable, grepl(cancerType, variantTable$cancer))
      variantTable <- variantTable[index]
      variantTable_DL <- variantTable_DL[index]
      
    }
    

    # aggregate
    coinciding <- setkeyv(coinciding, c('plot','set','type'))
    coinciding <- coinciding[, sum(n, na.rm = TRUE),by = c('plot','set','type')]
    
    if(!is.null(variantTable)){
      variantTable <- data.table(variantTable, key="gdna_pos")
      variantTable_DL <- data.table(variantTable_DL, key="gdna_pos")
    }
    
    
    # if (length(variantTable$gene) != length(variantTable_DL$gene)){
    #   print("VARIANT TABLES HAVE DIFFERENT LENGHTS!")
    # }
    
    coincidingDataList <- list(
      variantPlot = coinciding[coinciding$plot == "variantPlot",],
      consequencePlot = coinciding[coinciding$plot == "consequencePlot",],
      contextPlot = coinciding[coinciding$plot == "contextPlot",],
      variantTable = variantTable,
      variantTable_DL = variantTable_DL
    )
    
  })
  
  
  #### Somatic INPUT ####
  # Updates every time a somatic setting is changed
  somaticInput <- reactive({
    
    unique_somatic <- databaseInput()$unique_somatic
    coinciding_somatic <- databaseInput()$coinciding
    
    
    if(input$somatic == "cosmic"){
      unique_somatic <- unique_somatic[unique_somatic$cancer == input$cancer,]
      coinciding_somatic <- coinciding_somatic[coinciding_somatic$cancer == input$cancer,]
      
      if(input$cancerFreq != "ALL"){
        unique_somatic <- unique_somatic[unique_somatic$cancerFreq == input$cancerFreq,]
        coinciding_somatic <- coinciding_somatic[coinciding_somatic$cancerFreq == input$cancerFreq,]
      }
    }

    # aggregate
    unique_somatic <- setkeyv(unique_somatic, c('plot','set','type'))
    unique_somatic <- unique_somatic[, sum(n, na.rm = TRUE),by = c('plot','set','type')]
    
    coinciding_somatic <- setkeyv(coinciding_somatic, c('plot','set','type'))
    coinciding_somatic <- coinciding_somatic[, sum(n, na.rm = TRUE),by = c('plot','set','type')]
    
    
    somaticPlotsList <- list(
      
      coinciding_somatic_variantPlot = coinciding_somatic[coinciding_somatic$plot == "variantPlot",],
      coinciding_somatic_consequencePlot = coinciding_somatic[coinciding_somatic$plot == "consequencePlot",],
      coinciding_somatic_contextPlot = coinciding_somatic[coinciding_somatic$plot == "contextPlot",],
      
      variantPlot = unique_somatic[unique_somatic$plot == "variantPlot",],
      consequencePlot = unique_somatic[unique_somatic$plot == "consequencePlot",],
      contextPlot = unique_somatic[unique_somatic$plot == "contextPlot",]
      
    )
  })
  
  #### Germline INPUT ####
  # Updates every time databaseInput() is updated
  # Updates every time a germline setting is changed
  germlineInput <- reactive({
    
    # print("4 GERMLINE INPUT")
    unique_germline <- databaseInput()$unique_germline
    coinciding_germline <- databaseInput()$coinciding
    

    
    if(input$germline == "oneKG" | input$germline == "ExAC"){
      if (input$germline == "oneKG"){
        unique_germline <- unique_germline[unique_germline$region == input$oneKGpopulation,]
        coinciding_germline <- coinciding_germline[coinciding_germline$region == input$oneKGpopulation,]
      } else if (input$germline == "ExAC"){
        unique_germline <- unique_germline[unique_germline$region == input$ExACpopulation,]
        coinciding_germline <- coinciding_germline[coinciding_germline$region == input$ExACpopulation,]
      }
      if(input$MAF != "ALL"){
        unique_germline <- unique_germline[unique_germline$alleleFreq == input$MAF,]
        coinciding_germline <- coinciding_germline[coinciding_germline$alleleFreq == input$MAF,]
      }
      # dbSNP  
    } else {
      
      if(input$dbSNPvalidation != "ALL"){
        unique_germline <- unique_germline[unique_germline$validation == input$dbSNPvalidation,]
        coinciding_germline <- coinciding_germline[coinciding_germline$validation == input$dbSNPvalidation,]
      }
      if(input$dbSNPsubmissions != "ALL"){
        unique_germline <- unique_germline[unique_germline$submission == input$dbSNPsubmissions,]
        coinciding_germline <- coinciding_germline[coinciding_germline$submission == input$dbSNPsubmissions,]
      }
    }
    
    unique_germline <- setkeyv(unique_germline, c('plot','set','type'))
    unique_germline <- unique_germline[, sum(n, na.rm = TRUE),by = c('plot','set','type')]
    
    coinciding_germline <- setkeyv(coinciding_germline, c('plot','set','type'))
    coinciding_germline <- coinciding_germline[, sum(n, na.rm = TRUE),by = c('plot','set','type')]
    
    germlinePlotsList <- list(
      
      coinciding_germline_variantPlot = coinciding_germline[coinciding_germline$plot == "variantPlot",],
      coinciding_germline_consequencePlot = coinciding_germline[coinciding_germline$plot == "consequencePlot",],
      coinciding_germline_contextPlot = coinciding_germline[coinciding_germline$plot == "contextPlot",],
      
      variantPlot = unique_germline[unique_germline$plot == "variantPlot",],
      consequencePlot = unique_germline[unique_germline$plot == "consequencePlot",],
      contextPlot = unique_germline[unique_germline$plot == "contextPlot",]
      
    )
  })
  
  ### Setting up the corrext unique proportion of unique variants according to the venn diagram.
  plotInput <- reactive({
    
    germline_variantPlot <- germlineInput()$variantPlot
    germline_variantPlot$V1 <- (germlineInput()$variantPlot$V1 + germlineInput()$coinciding_germline_variantPlot$V1) - coincidingInput()$variantPlot$V1
    germline_contextPlot <- germlineInput()$context
    germline_contextPlot$V1 <- (germlineInput()$contextPlot$V1 + germlineInput()$coinciding_germline_contextPlot$V1) - coincidingInput()$contextPlot$V1
    germline_consequencePlot <- germlineInput()$consequencePlot
    germline_consequencePlot$V1 <- (germlineInput()$consequencePlot$V1 + germlineInput()$coinciding_germline_consequencePlot$V1) - coincidingInput()$consequencePlot$V1
    
    somatic_variantPlot <- somaticInput()$variantPlot
    somatic_variantPlot$V1 <- (somaticInput()$variantPlot$V1 + somaticInput()$coinciding_somatic_variantPlot$V1) - coincidingInput()$variantPlot$V1
    somatic_contextPlot <- somaticInput()$context
    somatic_contextPlot$V1 <- (somaticInput()$contextPlot$V1 + somaticInput()$coinciding_somatic_contextPlot$V1) - coincidingInput()$contextPlot$V1
    somatic_consequencePlot <- somaticInput()$consequencePlot
    somatic_consequencePlot$V1 <- (somaticInput()$consequencePlot$V1 + somaticInput()$coinciding_somatic_consequencePlot$V1) - coincidingInput()$consequencePlot$V1
    
    #### plotInput ####
    plotInputList <- list(
      
      germline_variantPlot = germline_variantPlot,
      germline_contextPlot = germline_contextPlot,
      germline_consequencePlot = germline_consequencePlot,
      somatic_variantPlot = somatic_variantPlot,
      somatic_contextPlot = somatic_contextPlot,
      somatic_consequencePlot = somatic_consequencePlot
      
    )
  })
  
  
  
  ####  Variant type PLOT  #### 
  output$variantPlot <- renderPlot({
    
    coinciding <- coincidingInput()$variantPlot
    
    validate(
      need(sum(coinciding$V1) != 0, "There are no coinciding dataset values left")
    )
    req(sum(coinciding$V1) != 0)
    
    
    
    unique <- switch(input$unique, 
                     "unique_somatic"=plotInput()$somatic_variantPlot, 
                     "unique_germline"=plotInput()$germline_variantPlot
    )
    
    
    # Ordering labels after filtering & aggregating.
    coinciding$type <- reorder.factor(coinciding$type, new.order = mutOrder$type)
    coinciding <- coinciding[order(coinciding$type), ]
    unique$type <- reorder.factor(unique$type, new.order = mutOrder$type)
    unique <- unique[order(unique$type), ]
    
    
    theTitle <- paste(input$somatic,input$germline, sep="_")
    
    types <- coinciding$type
    coinciding$V1 <- prop.table(coinciding$V1)
    unique$V1 <- prop.table(unique$V1)
    
    plotData <- rbind(unique, coinciding)
    plotData[0:6,"c"] <- c(input$unique)
    plotData[7:12,"c"] <- c('coinciding')
    
    plotData$type <- factor(plotData$type, levels = mutOrder$type)
   
    
    ggplot(plotData, aes(y=plotData$V1, x=plotData$type, fill=c)) +
      geom_bar(colour="black", position="dodge", stat="identity") +
      labs(x = "", y = "frelative requency") +
      theme_minimal() +
      theme(axis.text=element_text(size=12), axis.text.x = element_text(angle = 45, hjust = 1), axis.title=element_text(size=12)) +
      scale_fill_manual(values = c(coinciding="#0072B2", unique_germline="#E69F00",unique_somatic="#F0E442")) +
      theme(legend.text = element_text(size = 14, colour = "black"))+  
      scale_y_continuous(limits = c(0,max(plotData$V1)*1.2), expand = c(0, 0)) 
  })
  
  #### Variant enrichment/depletion PLOT ####  
  output$variantPlot_logarithm <- renderPlot({
    
    
    coinciding <- coincidingInput()$variantPlot
    
    validate(
      need(sum(coinciding$V1) != 0, "There are no coinciding dataset values left")
    )
    req(sum(coinciding$V1) != 0)
    
    
    unique <- switch(input$unique, 
                     "unique_somatic"=plotInput()$somatic_variantPlot, 
                     "unique_germline"=plotInput()$germline_variantPlot
    )

    
    # Ordering labels after filtering & aggregating.
    coinciding$type <- reorder.factor(coinciding$type, new.order = mutOrder$type)
    coinciding <- coinciding[order(coinciding$type), ]
    unique$type <- reorder.factor(unique$type, new.order = mutOrder$type)
    unique <- unique[order(unique$type), ]
    
    mutations <- coinciding$type
    mutations <- factor(mutations, levels = mutOrder$type)  
    
    test <- coinciding
    test2 <- unique
    test$V1 <- prop.table(test$V1)
    test2$V1 <- prop.table(test2$V1)
    
    values = prop.table(coinciding$V1)/prop.table(unique$V1)
    values <- log2(values)

    # 
    unique_color <- "#E69F00"
    if(input$unique == "unique_somatic"){
      unique_color <- "#F0E442"
    }  
    
    df <- data.frame(x = mutations, y = values)
    

    ggplot(df, aes(x=mutations, y=values)) +
      geom_bar(stat = "identity", position = "identity", colour="black",
               fill = ifelse(values > 0, "#0072B2", unique_color)) + 
      geom_text(aes(label = round(y, 2),
                    vjust = ifelse(y >= 0, 0, 1))) +
      theme_minimal() +
      theme(axis.text=element_text(size=12), axis.title=element_text(size=14,face="bold")) +
      scale_x_discrete(name ="") +
      scale_y_continuous(paste("log2(coinciding/",input$unique, ")", sep=""))
    
  })  
  
  #### Consequence PLOT ####  
  output$consequencePlot_coinciding <- renderPlot({
    
    tempData <- coincidingInput()$consequencePlot
    
    validate(
      need(sum(tempData$V1) != 0, "There are no coinciding dataset values left")
    )
    req(sum(tempData$V1) != 0)
    
    
    maxValue <- max(max(prop.table(coincidingInput()$consequencePlot$V1)), 
                    max(prop.table(plotInput()$somatic_consequencePlot$V1)), 
                    max(prop.table(plotInput()$germline_consequencePlot$V1)))
    
    
    # Ordering labels after filtering &  aggregating.
    inputData <- tempData
    inputData$type <- reorder.factor(inputData$type, new.order = consOrder$type)
    inputData <- inputData[order(inputData$type), ]
    
    
    order <- factor(inputData$type, levels = consequence_labels)
    values = prop.table(inputData$V1)
    df <- data.frame(x = order, y = values)
    
    # remove consequences that are zero
    df <- df[df$y > 0.0,]
    
    ggplot(df, aes(x=df$x, y=df$y, fill = df$x, width=0.70)) +
      geom_bar(stat = "identity", colour="black") +
      labs(x = "", y = "relative frequency")+
      theme_minimal() +
      theme(axis.text=element_text(size=12), axis.title=element_text(size=12)) +
      theme(legend.title=element_blank())+
      scale_y_continuous(limits = c(0,maxValue*1.05), expand = c(0, 0)) +
      scale_x_discrete(labels=df$x) +
      theme(axis.text.x = element_text(angle = 90, hjust = 1))+
      theme(legend.text = element_text(size = 10,
                                       colour = "black",
                                       angle = 0,
                                       face = "bold"))+     
      theme(legend.position = "none") 
    
  })
  
  #### Consequence PLOT germline ####  
  output$consequencePlot_germline <- renderPlot({
    
    
    tempData <- plotInput()$germline_consequencePlot
    
    validate(
      need(sum(tempData$V1) != 0, "There are no germline dataset values left")
    )
    req(sum(tempData$V1) != 0)
    
    maxValue <- max(max(prop.table(coincidingInput()$consequencePlot$V1)), 
                    max(prop.table(plotInput()$germline_consequencePlot$V1)), 
                    max(prop.table(plotInput()$somatic_consequencePlot$V1)))
    
    # Ordering labels after filtering &  aggregating.
    inputData <- tempData
    inputData$type <- reorder.factor(inputData$type, new.order = consOrder$type)
    inputData <- inputData[order(inputData$type), ]
    
    order <- factor(inputData$type, levels = consequence_labels)
    values = prop.table(inputData$V1)
    df <- data.frame(x = order, y = values)
    
    # remove consequences that are zero
    df <- df[df$y > 0.0,]
    
    ggplot(df, aes(x=df$x, y=df$y, fill = df$x, width=0.70)) +
      geom_bar(stat = "identity", colour="black") +
      labs( x = "", y = "relative frequency")+
      theme_minimal() +
      theme(axis.text=element_text(size=12), axis.title=element_text(size=12)) +
      scale_y_continuous(limits = c(0,maxValue*1.05), expand = c(0, 0)) +
      scale_x_discrete(labels=df$x) +
      theme(axis.text.x = element_text(angle = 90, hjust = 1))+
      theme(legend.text = element_text(size = 10,
                                       colour = "black",
                                       angle = 0,
                                       face = "bold"))+     
      theme(legend.position = "none") 
    
  })
  
  #### Consequence PLOT Somatic ####  
  output$consequencePlot_somatic <- renderPlot({
    
    tempData <- plotInput()$somatic_consequencePlot
    
    validate(
      need(sum(tempData$V1) != 0, "There are no somatic dataset values left")
    )
    req(sum(tempData$V1) != 0)
    
    maxValue <- max(max(prop.table(coincidingInput()$consequencePlot$V1)), 
                    max(prop.table(plotInput()$germline_consequencePlot$V1)), 
                    max(prop.table(plotInput()$somatic_consequencePlot$V1)))
    
    # Ordering labels after filtering &  aggregating.
    inputData <- tempData
    inputData$type <- reorder.factor(inputData$type, new.order = consOrder$type)
    inputData <- inputData[order(inputData$type), ]
    
    order <- factor(inputData$type, levels = consequence_labels)
    values = prop.table(inputData$V1)
    df <- data.frame(x = order, y = values)
    
    # remove consequences that are zero
    df <- df[df$y > 0.0,]
    
    ggplot(df, aes(x=df$x, y=df$y, fill = df$x, width=0.70)) +
      geom_bar(stat = "identity", colour="black") +
      labs(x = "", y = "relative frequency") +
      theme_minimal() +
      theme(axis.text=element_text(size=12), axis.title=element_text(size=12)) +
      scale_y_continuous(limits = c(0,maxValue*1.05), expand = c(0, 0)) +
      scale_x_discrete(labels=df$x) +
      theme(axis.text.x = element_text(angle = 90, hjust = 1))+
      theme(legend.text = element_text(size = 10,
                                       colour = "black",
                                       angle = 0,
                                       face = "bold"))+     
      theme(legend.position = "none") 
  })
  
  
  #### Consequence enrichment/depletion PLOT ####  
  output$consequencePlot_logarithm <- renderPlot({
    
    coinciding <- coincidingInput()$consequencePlot
    
    
    validate(
      need(sum(coinciding$V1) != 0, "There are not enought dataset values left")
    )
    req(sum(coinciding$V1) != 0)
    
    unique <- switch(input$unique, 
                     "unique_somatic"=plotInput()$somatic_consequencePlot, 
                     "unique_germline"= plotInput()$germline_consequencePlot
    )
    
    # Ordering labels after filtering & aggregating.
    coinciding$type <- reorder.factor(coinciding$type, new.order = consOrder$type)
    coinciding <- coinciding[order(coinciding$type), ]
    
    unique$type <- reorder.factor(unique$type, new.order = consOrder$type)
    unique <- unique[order(unique$type), ]
    
    ## removing consequence variants that have a zero value (form each dataset)
    remove <- coinciding[!coinciding$V1 > 0.0,]
    coinciding <- coinciding[coinciding$V1 > 0.0,]
    unique <- unique[ ! unique$type %in% remove$type, ]
    remove <- unique[!unique$V1 > 0.0,]
    unique <- unique[unique$V1 > 0.0,]
    coinciding <- coinciding[ ! coinciding$type %in% remove$type, ]
    
    order <- factor(coinciding$type, levels = consequence_labels)
    values = log2(prop.table(coinciding$V1)/prop.table(unique$V1))
    
    y1 = values
    values[!is.finite(values)] <- 0
    df <- data.frame(x = order, y = values)
    
    minValue <- min(values)*1.2
    maxValue <- max(values)*1.2
    
    ggplot(df, aes(x=df$x, y=df$y, fill = df$x, width=0.70)) +
      geom_bar(stat = "identity", colour="black") +
      labs(x = "",y = paste("log2(coinciding/",input$unique, ")", sep=""))+
      theme_minimal() +
      scale_y_continuous(limits = c(minValue,maxValue),expand = c(0, 0)) +
      scale_x_discrete(name = "",labels=df$x) +
      theme(plot.title = element_text(family = "Trebuchet MS", color="#666666", face="bold", size=16)) +
      theme(axis.text=element_text(size=12), axis.title=element_text(size=12,face="bold")) +
      theme(axis.text.x = element_text(angle = 90, hjust = 1))+
      geom_text(aes(label = ifelse(!is.na(y1), round(y1,2), 'NaN'),
                    vjust = ifelse(y >= 0, 0, 1))) +
      theme(legend.position = "none") 
  })
  
  #### Context PLOT ####  
  output$contextPlot_coinciding <- renderPlot({
    
    plotData <- coincidingInput()$contextPlot
    
    validate(
      need(sum(plotData$V1) != 0, "There are non coinciding variants left")
    )
    req(sum(plotData$V1) != 0)
    
    maxValue <- max(max(prop.table(coincidingInput()$contextPlot$V1)), 
                    max(prop.table(plotInput()$germline_contextPlot$V1)), 
                    max(prop.table(plotInput()$somatic_contextPlot$V1)))
    
    # Ordering labels after filtering & aggregating.
    plotData$type <- reorder.factor(plotData$type, new.order = contextOrder$type)
    plotData <- plotData[order(plotData$type), ]
    # Setting new labels
    plotData$b <- contextOrder$b
    plotData$type <- contextLabel
    # Relative frequency
    plotData$V1 <- prop.table(plotData$V1)  
    
    # print("contextPlot")
    # print(plotData)
    
    ggplot(data=plotData, aes(x=type, y=V1, fill = b, width=0.70))+
      geom_bar(stat="identity", 
               position = "identity")+
      labs(x = "", 
           y = "relative frequency")+
      theme_minimal() +
      theme(legend.title = element_blank(), 
            legend.position = "top", 
            legend.text = element_text(size = 15, colour = "black", angle = 0, face = "bold"), 
            legend.text.align = 0.5)+
      theme(axis.text.x  = element_text(angle=90,vjust=0.5,size=10,color="black"))+
      theme(axis.text.y = element_text(size = 12)) +
      scale_y_continuous(limits = c(0,maxValue*1.05), 
                         expand = c(0, 0)) +
      scale_fill_manual(values = contextColors) +
      scale_x_discrete(limits=plotData$type,  
                       expand = c(0, 0))+
      guides(fill=guide_legend(nrow=1,byrow=TRUE,
                               keywidth=1.0,
                               keyheight=0.2,
                               default.unit="inch",
                               label.position = "top")) 
  })
  
  #### Context PLOT Somatic ####  
  output$contextPlot_somatic <- renderPlot({
    
    plotData <- plotInput()$somatic_contextPlot
    
    validate(
      need(sum(plotData$V1) != 0, "There are no somatic dataset values left")
    )
    req(sum(plotData$V1) != 0)
    
    maxValue <- max(max(prop.table(coincidingInput()$contextPlot$V1)), max(prop.table(plotInput()$germline_contextPlot$V1)), max(prop.table(plotInput()$somatic_contextPlot$V1)))
    
    # Ordering labels after filtering & aggregating.
    plotData$type <- reorder.factor(plotData$type, new.order = contextOrder$type)
    plotData <- plotData[order(plotData$type), ]
    
    plotData$b <- contextOrder$b
    plotData$type <- contextLabel
    plotData$V1 <- prop.table(plotData$V1)  
    
    ggplot(data=plotData, aes(x=type, y=V1, fill = b, width=0.70))+
      geom_bar(stat="identity", position = "identity")+
      labs(x = "", y = "relative frequency")+
      theme_minimal() +
      theme(legend.position = "top") +
      theme(legend.text = element_text(size = 15, colour = "black", angle = 0, face = "bold"))+
      theme(legend.title.align = 0) +
      theme(legend.text.align = 0.5) +
      theme(legend.title = element_blank(), axis.text.x  = element_text(angle=90,vjust=0.5,size=10,color="black"))+
      theme(axis.text.y = element_text(size = 12)) +
      scale_y_continuous(limits = c(0,maxValue*1.05), expand = c(0, 0)) +
      scale_fill_manual(values = contextColors) +
      scale_x_discrete(limits=plotData$type,  expand = c(0, 0))+
      guides(fill=guide_legend(nrow=1,byrow=TRUE,
                               keywidth=1.0,
                               keyheight=0.2,
                               default.unit="inch",
                               label.position = "top"))
  })
  
  #### Context PLOT Germline ####
  output$contextPlot_germline <- renderPlot({
    
    plotData <-plotInput()$germline_contextPlot
    
    validate(
      need(sum(plotData$V1) != 0, "There are no germline dataset values left")
    )
    req(sum(plotData$V1) != 0)
    
    maxValue <- max(max(prop.table(coincidingInput()$contextPlot$V1)), max(prop.table(plotInput()$germline_contextPlot$V1)), max(prop.table(plotInput()$somatic_contextPlot$V1)))
    
    # Ordering labels after filtering & aggregating.
    plotData$type <- reorder.factor(plotData$type, new.order = contextOrder$type)
    plotData <- plotData[order(plotData$type), ]
    
    plotData$b <- contextOrder$b
    plotData$type <- contextLabel
    plotData$V1 <- prop.table(plotData$V1)  
    
    ggplot(data=plotData, aes(x=type, y=V1, fill = b, width=0.70))+
      geom_bar(stat="identity", position = "identity")+
      labs(x = "", y = "relative frequency")+
      theme_minimal() +
      theme(legend.position = "top") +
      theme(legend.text = element_text(size = 15, colour = "black", angle = 0, face = "bold"))+
      theme(legend.title.align = 0) +
      theme(legend.text.align = 0.5) +
      theme(legend.title = element_blank(), axis.text.x  = element_text(angle=90,vjust=0.5,size=10,color="black"))+
      theme(axis.text.y = element_text(size = 12)) +
      scale_y_continuous(limits = c(0,maxValue*1.05), expand = c(0, 0)) +
      scale_fill_manual(values = contextColors) +
      scale_x_discrete(limits=plotData$type,  expand = c(0, 0))+
      guides(fill=guide_legend(nrow=1,byrow=TRUE,
                               keywidth=1.0,
                               keyheight=0.2,
                               default.unit="inch",
                               label.position = "top"))
  })
  
  
  #### Context enrichment/depletion PLOT ####  
  output$contextPlot_logarithm <- renderPlot({
    
    coincidingData <- coincidingInput()$contextPlot
    
    validate(
      need(sum(coincidingData$V1) != 0, "There are not enought dataset values left")
    )
    req(sum(coincidingData$V1) != 0)
    
    uniqueData <- switch(input$unique, 
                         "unique_somatic"=plotInput()$somatic_contextPlot, 
                         "unique_germline"=plotInput()$germline_contextPlot
    )
    
    # Ordering labels after filtering & aggregating.
    coincidingData$type <- reorder.factor(coincidingData$type, new.order = contextOrder$type)
    coincidingData <- coincidingData[order(coincidingData$type), ]
    uniqueData$type <- reorder.factor(uniqueData$type, new.order = contextOrder$type)
    uniqueData <- uniqueData[order(uniqueData$type), ]
    
    coincidingData$V1 <- prop.table(coincidingData$V1)  
    uniqueData$V1 <- prop.table(uniqueData$V1)  
    res <- coincidingData$V1/uniqueData$V1
    
    plotData <- contextOrder
    plotData$type <- contextLabel
    
    plotData$n <- log2(res)
    n1 =  plotData$n
    plotData$n[!is.finite(plotData$n)] <- 0
    
    ggplot(data=plotData, aes(x=type, y=n, fill = b, width=0.70))+
      geom_bar(stat="identity", position = "identity")+
      labs( y = paste("log2(coinciding/", input$unique,") ", paste= "" ))+
      theme_minimal() +
      theme(plot.title = element_text(family = "Trebuchet MS", color="#666666", face="bold", size=16)) +
      theme(legend.position = "top") +
      theme(legend.text = element_text(size = 15,
                                       colour = "black",
                                       angle = 0,
                                       face = "bold"))+
      theme(axis.title = element_text(face="bold")) +
      theme(legend.title = element_blank(),axis.text.x  = element_text(angle=90,
                                                                       vjust=0.5, 
                                                                       size=10, 
                                                                       color="black"))+
      theme(axis.text.y = element_text(size = 12)) +
      scale_y_continuous(limits = c(min(plotData$n)*1.5,max(plotData$n)*1.5), expand = c(0, 0)) +
      scale_fill_manual(values = contextColors) +
      theme(legend.text.align = 0.5) +
      scale_x_discrete(name = "", limits=plotData$type,  expand = c(0, 0))+
      guides(fill=guide_legend(nrow=1,byrow=TRUE,
                               keywidth=1.0,
                               keyheight=0.2,
                               default.unit="inch",
                               label.position = "top"))+
      geom_text(aes(label = ifelse((is.finite(n1)), "", n1),
                    vjust = ifelse(n >= 0, 0, 1))) 
  })
  
  #### Data TABLES ####
  output$dataTable <- DT::renderDataTable(
    
    coincidingInput()$variantTable,
    filter = 'bottom',
    options = list(pageLength = 20),
    extensions = "Responsive",
    escape = FALSE
  )
  
  
  observe({
    
    if(sum(coincidingInput()$contextPlot$V1) != 0){
      shinyjs::show("downloadContext")
      shinyjs::show("downloadConsequence")
    } else {
      shinyjs::hide("downloadContext")
      shinyjs::hide("downloadConsequence")
    }
    
    if(sum(coincidingInput()$contextPlot$V1) == 0 | input$somatic == "icgc"){
      shinyjs::hide("downloadData")
    } else {
      shinyjs::show("downloadData")
    }
    
  })
  
  
  output$dataTableMessage <- renderText({
    
    if(input$somatic == "icgc"){
      "Only COSMIC have coinciding datatables available"
    }
    
  })
  
  
  downloadDatatable <- reactive({
    
    dl <- coincidingInput()$contextPlot
    dl <- dl[,c("type","V1")]
    
    dl$type <- reorder.factor(dl$type, new.order = contextOrder$type)
    dl <- dl[order(dl$type), ]
    dl$type <- contextLabel
    
    dl2 <- coincidingInput()$consequencePlot
    dl2 <- dl2[,c("type","V1")]
    
    list("variantDatatable" = coincidingInput()$variantTable_DL,
         "contextTable" = dl,
         "consequenceTable" = dl2
    )
  })
  
  
  downloadFileName <- reactive({
    
    
    if(input$variantConsequence != "ALL"){
      name <- paste0(input$variantConsequence,"_")
    } else {
      name <- ""
    }
    
    if(input$cancer != 'pancancer'){
      name <- paste0(name,input$cancer, "_")
    }
    if(input$cancerFreq != 'ALL'){
      name <- paste0(name,input$cancerFreq ,"_")
    }
    
    if (input$germline == "oneKG" | input$germline == "ExAC"){
      if(input$germline == "oneKG"){
        
        name <- paste0(name ,input$oneKGpopulation, "_")
        
      } else if(input$germline == "ExAC"){
        
        name <- paste0(name ,  input$ExACpopulation,"_")
        
      }
      if(input$MAF!= 'ALL'){
        name <- paste0(name , input$MAF,"_")
      }
      
    } else {
      
      if(input$dbSNPvalidation != "ALL"){
        name <- paste0(name ,input$dbSNPvalidation,"_")
      }
      if(input$dbSNPsubmissions != "ALL"){
        name <- paste0(name ,input$dbSNPsubmissions,"_")
      }
    }
    name <- paste0(name,input$somatic, "_",input$germline,"_coinciding")
    return(name)
  })
  
  output$downloadVariantTable <- downloadHandler(
    filename = function() { 
      paste0(downloadFileName(),'_variantTable.csv') 
    },
    content <- function(file) {
      write.csv(downloadDatatable()$variantDatatable, file)
    }
  )
  
  output$downloadContext <- downloadHandler(
    filename = function() { paste0(downloadFileName(), '_context.csv') },
    content <- function(file) {
      write.table(downloadDatatable()$contextTable, file, sep = ',', row.names = FALSE)
    }
  )
  
  output$downloadConsequence <- downloadHandler(
    filename = function() { paste0(downloadFileName(), '_consequence.csv') },
    content <- function(file) {
      write.table(downloadDatatable()$consequenceTable, file, sep = ',', row.names = FALSE)
    }
  )
  
  
  
  #### STATETISTICS ####
  output$somaticStats <- renderUI({
    
    somatic <- switch(input$somatic,
                      "cosmic" = paste0("<a href='http://cancer.sanger.ac.uk/cosmic/about' target='_blank'>COSMIC</a>"),
                      "icgc" = paste0("<a href='http://icgc.org/#about' target='_blank'>ICGC</a>"))
    
    if(input$variantConsequence != "ALL"){
      choices <- paste0("-Variant consequence: ",input$variantConsequence,"<br/>")
    } else {
      choices <- ""
    }
    
    if (input$somatic == "cosmic"){
      
      
      cancerType <- switch (input$cancer,
                            "chronic_lymphocytic_leukemia" = "Chronic lymphocytic leukemia",
                            "acute_myeloid_leukemia" = "Acute myeloid leukemia",
                            "acute_lymphoblastic_b-cell_leukemia" = "Acute lymphoblastic b-cell leukemia",
                            "acute_lymphoblastic_leukemia" = "Acute lymphoblastic leukemia",
                            "lung_cancer" = "Lung cancer", 
                            "breast_cancer" = "Breast cancer", 
                            "malignant_melanoma" = "Malignant melanoma",
                            "sarcoma" = "Sarcoma",
                            "stomach_cancer" = "Stomach cancer",
                            "prostate_cancer" = "Prostate cancer",
                            "diffuse_large_B_cell_lymphoma" = "Diffuse large B-cell lymphoma",
                            "glioma" = "Glioma",
                            "colorectal_cancer" = "Colorectal cancer",
                            "oesophageal_cancer" = "Oesophageal cancer",
                            "ovarian_cancer" = "Ovarian cancer",
                            "pancreatic_cancer" = "Pancreatic cancer",
                            "liver_cancer" = "Liver cancer",
                            "cervical_cancer" = "Cervical cancer",
                            "cholangiocarcinoma" = "Cholangiocarcinoma",
                            "urothelical_cancer" = "Urothelical cancer",
                            "kidney_cancer" = "Kidney cancer",
                            "pancancer" = "Any cancer")
      
      if(cancerType != 'Any cancer'){
        choices <- paste0("-Cancer type: ",cancerType, "<br/>")
      }
      if(input$cancerFreq != 'ALL'){
        choices <- paste0(choices,"-Cancer type frequency: ",input$cancerFreq ,"<br/>")
      }
    }
    outText <- paste0("<h4>Somatic variant database: ",somatic, "</h4>")
    outText <- paste0(outText,choices)
    somatic_sum = sum(plotInput()$somatic_variantPlot$V1) +  sum(coincidingInput()$variantPlot$V1)
    somatic_number <- format(somatic_sum, digits=nchar(somatic_sum), decimal.mark=",", big.mark=" ",small.mark=".", small.interval=3)
    outText  <- paste0(outText,"-Size: ", somatic_number, " variants", "<br/>")
    HTML(outText)
  })
  
  output$germlineStats <- renderUI({
    
    germline <- switch(input$germline, 
                       "oneKG"=paste0("<a href='http://www.internationalgenome.org/about' target='_blank'>1000 Genomes Project</a>"),
                       "dbSNP"=paste0("<a href='https://www.ncbi.nlm.nih.gov/snp' target='_blank'>dbSNP</a>"),
                       "ExAC"=paste0("<a href='http://exac.broadinstitute.org/about' target='_blank'>ExAC</a>"))
    
    if(input$variantConsequence != "ALL"){
      choices <- paste0("-variant consequence: ",input$variantConsequence,"<br/>")
    } else {
      choices <- ""
    }
    
    if (input$germline == "oneKG" | input$germline == "ExAC"){
      if(input$germline == "oneKG"){
        
        
        population_1kg <- switch(
          input$oneKGpopulation,
          "AFR_AF_1KG" = "African",
          "EUR_AF_1KG" = "European",
          "SAS_AF_1KG" = "South Asian",
          "AMR_AF_1KG" = "American",
          "EAS_AF_1KG" = "East Asian"
        )
        choices <- paste0(choices ,"-Populations: ", population_1kg, "<br/>")
        
        
      } else if(input$germline == "ExAC"){
        
        
        exac_population <- switch(
          input$ExACpopulation,
          "AF_AFR_EXAC" = "African/African American",
          "AF_AMR_EXAC" = "American",
          "AF_Adj_EXAC" = "Adjusted Global",
          "AF_NFE_EXAC" = "Non-Finnish European",
          "AF_FIN_EXAC" = "Finnish",
          "AF_EAS_EXAC" = "East Asian",
          "AF_SAS_EXAC" = "South Asian"
        )
        choices <- paste0(choices ,"-Populations: ", exac_population,"<br/>")
        
      }
      if(input$MAF!= 'ALL'){
        freq_choice <- switch(input$MAF, 
                              "Common"="Common", 
                              "LowFreq"="Low Frequency", 
                              "Rare"="Rare", 
                              "VeryRare"="Very Rare")
        choices <- paste0(choices , "-Allele frequency: ", freq_choice,"<br/>")
      }
    } else {
      
      if(input$dbSNPvalidation != "ALL"){
        choices <- paste0(choices , "-Validation: ",input$dbSNPvalidation,"<br/>")
      }
      if(input$dbSNPsubmissions != "ALL"){
        choices <- paste0(choices , "-Submissions: ",input$dbSNPsubmissions,"<br/>")
      }
    }
    
    outText <- paste0("<h4>Germline variant database: ",germline, "</h4>")
    outText <- paste0(outText,choices)
    germline_sum = sum(plotInput()$germline_variantPlot$V1) + sum(coincidingInput()$variantPlot$V1)
    germline_number <- format(germline_sum, digits=nchar(germline_sum), decimal.mark=",", big.mark=" ",small.mark=".", small.interval=3)
    outText  <- paste0(outText,"-Size: ", germline_number, " variants", "<br/>")
    HTML(outText)
  })
  
  #### Overall statistics ####
  output$coincidingStats <- renderText({
    
    germline = input$germline
    somatic = input$somatic
    
    # temporary filtered
    coinciding_temp_sum <- sum(coincidingInput()$variantPlot$V1)
    
    somatic_temp_sum <- sum(plotInput()$somatic_variantPlot$V1) 
    germline_temp_sum <- sum(plotInput()$germline_variantPlot$V1) 
    
    coinciding_part_of_somatic <- percent(prop.table(c(coinciding_temp_sum, somatic_temp_sum))[[1]])
    coinciding_part_of_germline <- percent(prop.table(c(coinciding_temp_sum, germline_temp_sum))[[1]])
    
    unique_part_of_somatic <- percent(prop.table(c(coinciding_temp_sum, somatic_temp_sum))[[2]])
    unique_part_of_germline <- percent(prop.table(c(coinciding_temp_sum, germline_temp_sum))[[2]])
    
    coinciding_number <- format( coinciding_temp_sum, digits=nchar(coinciding_temp_sum), decimal.mark=",", big.mark=" ",small.mark=".", small.interval=3)
    outText <- paste("<h4>Results from intersection of germline and somatic variant loci:</h4> ","1. Number of <em>coinciding</em> variants: ", coinciding_number, " variants (", coinciding_part_of_somatic," of somatic, ",coinciding_part_of_germline," of germline)", sep = "")   
    somatic_number <- format( somatic_temp_sum, digits=nchar(somatic_temp_sum), decimal.mark=",", big.mark=" ",small.mark=".", small.interval=3)
    outText <- paste0(outText,"<br/>2. Number of <em>unique somatic</em> variants: ", somatic_number," variants (",unique_part_of_somatic," of somatic)")
    germline_number <- format(germline_temp_sum, digits=nchar(germline_temp_sum), decimal.mark=",", big.mark=" ",small.mark=".", small.interval=3)
    outText <- paste0(outText,"<br/>3. Number of <em>unique germline</em> variants: ", germline_number," variants (",unique_part_of_germline," of germline)")
    
    HTML(outText)
  })
  
  
  
  
  #### deContructSigs ####
  deConstructSigs <- reactive({
    
    plotData <- coincidingInput()$contextPlot
    
    
    validate(
      need(sum(plotData$V1) != 0, "There are not enought dataset values left")
    )
    req(sum(plotData$V1) != 0)
    
    # Ordering labels after filtering & aggregating.
    plotData$type <- reorder.factor(plotData$type, new.order = contextOrder$type)
    plotData <- plotData[order(plotData$type), ]
    
    plotData$V1 <- prop.table(plotData$V1)
    plotData$type <- contextLabel
    
    df <- plotData[,c("type","V1")]
    
    contextPlotWeights <- matrix(df$V1,nrow=1)
    contextPlotWeights <- data.frame(contextPlotWeights)
    colnames(contextPlotWeights) <- df$type 
    contextPlotWeights <- contextPlotWeights[1,]
    
    # print("mutationPlot")
    # print(contextPlotWeights)
    
    res <- whichSignatures(tumor.ref = contextPlotWeights, 
                           signatures.ref = signatures.cosmic,
                           sample.id = 1,
                           contexts.needed = TRUE,
                           tri.counts.method = 'default'
    )
    
    list("results" = res)
    
  })
  
  
  output$signatures <- renderUI({
    
    res <- deConstructSigs()$results
    req(!is.null(res))
    

    
    signatures <- paste0(colnames(res$weights[1]), " : ", paste0(round(res$weights[1]*100, 3), "%"))
    if(length(res$weights) > 1){
      for(sig in 2:length(res$weights)){
        if(res$weights[sig] > 0){
          signatures <- paste(signatures, paste0(colnames(res$weights[sig]), " : ", paste0(round(res$weights[sig]*100, 3),"%" )), sep = " & ")
        }
      }
    }
    HTML(paste0("<strong>",signatures,"</strong>"))
  })
  
  #### Signature Pie ####
  output$signaturePie <- renderPlot({
    
    res <- deConstructSigs()$results
    req(!is.null(res))
    names(res$unknown) <- "coinciding context pattern"
    rownames(res$weights) <- "Resembling signatures"
    makePie(res)
  })
  
  output$signatures<- renderUI({
    res <- deConstructSigs()$results
    req(!is.null(res))
    signatures <- paste0(colnames(res$weights[1]), " : ", paste0(round(res$weights[1]*100, 2), "%"))
    if(length(res$weights) > 1){
      for(sig in 2:length(res$weights)){
        if(res$weights[sig] > 0){
          signatures <- paste(signatures, paste0(colnames(res$weights[sig]), " : ", paste0(round(res$weights[sig]*100, 2),"%" )), sep = "</br>")
        }
      }
    }
    signatures <- paste(signatures,paste0("</br>Unknown", " : ", paste0(round((1-sum(res$weights))*100, 2), "%")))
    
    signatureCosmicSite <- "You can read more information about signatures here: <a href='http://cancer.sanger.ac.uk/cosmic/signatures' target='_blank'>http://cancer.sanger.ac.uk/cosmic/signatures</a>"
    HTML(paste0("<strong>",signatures,"</strong></br>",signatureCosmicSite))
  })
  
  #### signature PLOT ####
  output$signatureContextPlot <- renderPlot({
    
    res <- deConstructSigs()$results
    req(!is.null(res))
    
    tumor <- res$product
    signature <- data.frame("type" = colnames(tumor), "V1" = tumor[0:96])
    
    
    plotData <- signature["type"]
    plotData["signature"] <- signature$V1
    plotData$b <- contextOrder$b
    
  
    
    ggplot(data=plotData, aes(x=type, y=plotData["signature"], fill = b, width=0.70))+
      geom_bar(stat="identity", position = "identity")+
      labs(x = "", y = "relative frequency")+
      theme_minimal() +
      theme(legend.position = "top") +
      theme(legend.text = element_text(size = 15, colour = "black", angle = 0, face = "bold"))+
      theme(legend.title.align = 0) +
      theme(legend.text.align = 0.5) +
      theme(legend.title = element_blank(), axis.text.x  = element_text(angle=90,vjust=0.5,size=10,color="black"))+
      theme(axis.text.y = element_text(size = 12)) +
      scale_y_continuous(limits = c(0,max(plotData["signature"])*1.05), expand = c(0, 0)) +
      scale_fill_manual(values = contextColors) +
      scale_x_discrete(limits=plotData$type,  expand = c(0, 0))+
      guides(fill=guide_legend(nrow=1,byrow=TRUE,
                               keywidth=1.0,
                               keyheight=0.2,
                               default.unit="inch",
                               label.position = "top")) 
  })
  
  
  
  output$ContextPlotwarning <- renderUI({
    
    outText <- c("")
    coincidingValues <-coincidingInput()$contextPlot$V1
    req(sum(coincidingValues) != 0)
    if (0 %in% coincidingValues){
      
      outText <- paste0(outText, "Warning: coinciding dataset contains missing values</br>")
    }  
    if (0 %in% plotInput()$somatic_contextPlot$V1 == TRUE){
      
      outText <- paste0(outText, "Warning: somatic dataset contains missing values</br>")
    }
    if (0 %in% plotInput()$germline_contextPlot$V1 == TRUE){
      
      outText <- paste0(outText, "Warning: germline dataset contains missing values</br>")
    } 
    # This text is <span style="color:red">red</span></div>
    HTML(paste0("<div><span style='color:red'>",outText,"</span></div>"))
    
  })
  
  #### venn diagram ####
  output$vennDiagram <- renderPlot({
    
    germline_db <- switch (input$germline,
                           "oneKG" = "1000Genomes",
                           "ExAC" = "ExAC",
                           "dbSNP" = "dbSNP"
    )
    
    somatic_db <- switch (input$somatic,
                          "cosmic" = "COSMIC",
                          "icgc" = "ICGC"
    )
    
    coinciding_temp_sum <- sum(coincidingInput()$variantPlot$V1)
    somatic_temp_sum <- sum(plotInput()$somatic_variantPlot$V1) + sum(coincidingInput()$variantPlot$V1)
    germline_temp_sum <- sum(plotInput()$germline_variantPlot$V1) + sum(coincidingInput()$variantPlot$V1)
    
    req(coinciding_temp_sum > 0)
    
    vennData <- prop.table(c(somatic_temp_sum, germline_temp_sum, coinciding_temp_sum))
    
    vennData <- vennData * 100
    grid.newpage()
    draw.pairwise.venn(vennData[1], vennData[2], vennData[3], 
                       category = c(somatic_db, germline_db), 
                       lty = rep("blank", 2), 
                       fill = c("light blue", "pink"), 
                       alpha = rep(0.5, 2),
                       cat.pos = c(0, 0), 
                       cat.dist = rep(0.025, 2),
                       print.mode = FALSE,
                       ext.text = FALSE,
                       sigdigs=4
    )
  })
  
  
  output$signatureTable <- renderTable({
    signatures <- deConstructSigs()$results
    index <- which(signatures$weights != 0)
    signatureTable[index]
    
  })
  
}
