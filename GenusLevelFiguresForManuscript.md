Skin Microbiome Dynamics Following *B. subtilis* formulation Challange
----------------------------------------------------------------------

##### Load necassary libraries.

    library(RUVSeq)
    library(DESeq2)
    library(ggplot2)
    library(RColorBrewer)
    library(pheatmap)
    library(grid)
    library(reshape2)
    library(gtools) 
    library(dplyr)
    library(EDASeq)
    library(vsn)
    library(tidyverse)

The workspace includes the abundance data from two seperate sequencing
experiments as well as the design of the experiment as is described in
the manuscript.

The abundance matrix used here excludes Day 6.

### Load input data

    #The experimental design 
    load(file = paste0("/Bigdata/Dropbox (Technion Dropbox)/Rina_Benel/Home/BKU/Mizrahi_veronica/data/GenusMergedDesignFirstSecond_NoSixthDayNoNew.RData"))

    #The abundance data 
    load(file = paste0("/Bigdata/Dropbox (Technion Dropbox)/Rina_Benel/Home/BKU/Mizrahi_veronica/data/GenusMergedDataFirstSecond_NoSixthDayNoNew.RData"))

#### Filter the data

    #counts the number of rows in the final matrix 
    nrow(FinalMergedData)

    ## [1] 710

    ##we decided to filter for 50 reads in at least 2 samples
    filter <- apply(FinalMergedData, 1, function(x) length(x[x>50])>= 2)

    ##subset the raw abundance data 
    filteredGenusLevel <- FinalMergedData[filter, ]

    #check now the number of rows
    nrow(filteredGenusLevel)

    ## [1] 286

    #for downstream processes convert the input data into an expression set 
    #create an expression set using the assay and pheno data
    ExpressionSet <- newSeqExpressionSet(as.matrix(filteredGenusLevel), 
                                         phenoData = mergedDesign, 
                                         row.names = colnames(filteredGenusLevel))

### RUVseq: Estimating the factors of unwanted variation using replicate samples

    #construct a matrix that specifies the replicate samples 
    #we will use the ruvseq column for this, each row indicates a set of replicate samples 
    differences <- RUVSeq::makeGroups(mergedDesign$RUVseq) #use the ruvseq column for the design 

    #Run ruvseq with the differences and the number of desired k's 
    ExpressionSet2 <- RUVSeq::RUVs(ExpressionSet, rownames(filteredGenusLevel), k = 40, differences)

    #examine all data now included in the phenotypic data 
    #pData(ExpressionSet2)

### PCA

    #PCA divided by treatments 
    colorsDay <- c("#000000", rep(c("#1805f0", "#fc0f03", "#1fa66b"), 5))

    #Convert the Colors into a numeric vector which matches the color assignment 
    colorsDayNumeric <- colorsDay[as.numeric(mergedDesign$RUVseq)]

    #Add the shapes to define the days
    shapes = c(1, rep(18, 3), rep(17, 3), rep(15,3), rep(8, 3), rep(16, 3))

    #Convert the shape vector into a numeric vector which matches the shape assignment 
    shapesNumeric <- shapes[as.numeric(mergedDesign$RUVseq)]

    #The margins of the plot 
    par(mar=c(11, 4.1, 2, 5.1), xpd=TRUE)

    #the letter k represents how many PC's you want displayed. 
    pca_data <- EDASeq::plotPCA(ExpressionSet2, col=colorsDayNumeric, k = 2, 
                                cex=2.5, labels = F, pch=shapesNumeric)
    #this is for manaul input into legend... 
    legend("bottom",inset=c(0, -0.8), legend=levels(mergedDesign$RUVseq), pch=shapes, 
           title="Treatment Groups",cex=0.8, ncol = 3, col=colorsDay)

![](GenusLevelFiguresForManuscript_files/figure-markdown_strict/unnamed-chunk-5-1.png)

### Differential expression with DESeq2

    #take input and make a deseq2 object
    #the design columns are incoporates into the ExpressionSet2 when using RUVseq 
    dds <- DESeqDataSetFromMatrix(countData = counts(ExpressionSet2),
                                  colData = pData(ExpressionSet2),
                                  design = ~ W_1+W_2+W_3+W_4+W_5+W_6+W_7+W_8+W_9+W_10+W_11+W_12+W_13+W_14+
                               W_15+W_16+W_17+W_18+W_19+W_20+W_21+W_22+W_23+W_24+W_25+W_26+W_27+W_28+W_29+
                                    W_30+W_31+W_32+W_33+W_34+W_35+W_36+W_37+W_38+W_39+W_40+RUVseq)

    #run deseq function
    dds <- DESeq(dds)
    #clean the rows that are no converging in beta
    ddsClean <- dds[which(mcols(dds)$betaConv),]

    #find which result names can be used for contrasts
    res <- resultsNames(ddsClean)

    #find which result names can be used for contrasts
    res <- resultsNames(ddsClean)

    #open a df for all of the DE results
    DEresultsDF <- data.frame()

    #get all the contrasts that start with RUVseq
    contrasts <- grep("RUVseq", res)

    for(i in contrasts){
      #split the contrasts names to a vector
      splitContrast <- unlist(strsplit(res[i], "_"))
      splitContrast <- splitContrast[splitContrast != "vs"] #remove the word vs so we have a three element vector
      res_df <- results(ddsClean, contrast = splitContrast) #call the results function 
      
      #covert to DF from deseq2 object 
      res_df <- as.data.frame(res_df)
      
      #this will be used for downstream analysis
      #if the padj is more than 0.05 input 0 in the FC column
      res_df$log2FoldChange <- ifelse(res_df$padj > 0.05, 0, res_df$log2FoldChange) 
      #If the FC, is less than abs (1) we will input 0 in FC column
      res_df$log2FoldChange <- ifelse(abs(res_df$log2FoldChange) < 1, 0, res_df$log2FoldChange) 
      
      #change the colnames so they reflect each contrast 
      colnames(res_df) <- paste0(colnames(res_df), "_", res[i])
      #add to final data frame   
      if(nrow(DEresultsDF) == 0) {
        DEresultsDF <- res_df
      } else {
        DEresultsDF <- cbind(DEresultsDF, res_df)
      }
    }

### Manipulate final results DF

    #retrieve the row that reflect the contrasts of b. subtilis
    Bacillus_ResData <-  DEresultsDF[rownames(DEresultsDF) == "Bacillus", ]

    #retrieve only the padj values for all contrasts for Bacillus subtilis
    BacillusPadj.df <- as.data.frame(t(dplyr::select(Bacillus_ResData, matches("padj"))))

    #change the name to padj
    colnames(BacillusPadj.df) <- "padj"

    #add an additional column that defines if significant 
    BacillusPadj.df$significant <- ifelse(BacillusPadj.df$padj < 0.05, "sig", "non-sig")

    #move the rownames to a column 
    BacillusPadj.df <- tibble::rownames_to_column(BacillusPadj.df, "treatment")
    #clean up the names for the graph 
    BacillusPadj.df$treatment <- gsub("_vs_notreatment", "", BacillusPadj.df$treatment)
    BacillusPadj.df$treatment <- gsub("padj_RUVseq_", "", BacillusPadj.df$treatment)

### Get the normalizes counts

    #get normalized counts through EDAseq
    NormCounts <- as.data.frame(EDASeq::normCounts(ExpressionSet2))

    #let's make the column names simpler and easier to manipulate 
    samples_name <- as.list(paste(mergedDesign$RUVseq, mergedDesign$mouse, sep = "_"))

    #now change colnames
    colnames(NormCounts) <- samples_name

    #remove the no treatment samples from the DF
    #will add them to the graph after calculating the mean and sd 
    NormCounts_Treatment <- NormCounts[, -grep("notreatment*", colnames(NormCounts))]

### Find the mean for the no treatment samples

    #get the "control" DF
    NormCounts_NoTreatment <-  NormCounts[, grep("notreatment*", colnames(NormCounts))]

    #Extract just the Bacillus normalized counts 
    bacillusControl <- NormCounts_NoTreatment[rownames(NormCounts_NoTreatment) == "Bacillus subtilis", ]

    #get Mean
    bacillusMean <- rowSums(bacillusControl)/ncol(bacillusControl)

### Get count data

    #manipulate the count data 
    #add the rownames as a column 
    counts <- tibble::rownames_to_column(NormCounts_Treatment, "id")

    #change to long format 
    #gather by treamtent type
    counts_long <- counts %>% tidyr::gather(sample, count, 2:ncol(counts))

    #remove the mouse/replicate so that we can find mean, sd 
    counts_long$sample = gsub("_A|_B|_C|_D|_E|_F", "", counts_long$sample)

    #group now by ID as well as treatment type, and get mean
    new.df <- counts_long %>% group_by(id, sample) %>% summarize(avg = mean(count))

    ## `summarise()` regrouping output by 'id' (override with `.groups` argument)

    # extract the data for a specific genus
    ind = grep("Bacillus", new.df$id) 

    df = new.df[ind, ]

    #merge the count data with the padj and significance 
    df.pvalue <- merge(df, BacillusPadj.df, by.x = "sample", by.y = "treatment")

    #add the following columns in order to graph the seperate treatments 
    #add column with the day and relvel 
    df.pvalue$day <- sapply(strsplit(as.character(df.pvalue$sample), ".", fixed = TRUE), '[', 1)
    df.pvalue$day <-factor(df.pvalue$day, levels=c("day2", "day4", "day8", "day11", "day14"))

    #add column with the treament info 
    df.pvalue$treatment <- sapply(strsplit(as.character(df.pvalue$sample), ".", fixed = TRUE), '[', 2)

### Graph of *Bacillus* counts

    #because there are many points near each other
    pd <- position_dodge(0.05) # move them .05 to the left and right

    #plot data
    df.pvalue.plot <- ggplot(df.pvalue, aes(x=day, y=avg)) + 
      geom_line(position = pd, size =1, aes(group = treatment, color = treatment)) +
      geom_point(position = pd, size = 5, aes(shape = significant)) + 
      ylab("Bacillus Normalized Counts") + xlab("") + 
      geom_hline(yintercept=bacillusMean, linetype="dashed", color = "black", size = 1.5) + 
      theme(legend.text = element_text(colour="black", size = 13, face = "plain")) +
      theme( axis.title.x = element_text(family="sans",size = 16, face="bold", hjust=0.5, vjust=-0.5),
             axis.title.y = element_text(family="sans",size = 16, angle=90, face="bold", hjust=0.5, vjust=1)) +
      theme( axis.text.x = element_text(family = "sans",size = 14, angle=45, face='plain', colour="#353535",   
                                        hjust=1, vjust=1) ) +
      theme( axis.text.y = element_text(family = "sans",size = 14, face='plain', colour="#353535",  vjust=0.5)) +
      theme(axis.line.x = element_line(color="black", size = 0.5),
            axis.line.y = element_line(color="black", size = 0.5)) +
      theme(legend.background = element_rect()) + 
      theme(legend.position="bottom") +
      theme(legend.direction = "vertical", legend.box = "horizontal") +
      theme(panel.border = element_blank(), panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),  panel.background = element_blank()) +
      scale_shape_manual(values = c(1, 16))


    print(df.pvalue.plot)

![](GenusLevelFiguresForManuscript_files/figure-markdown_strict/unnamed-chunk-11-1.png)

### Barplot of Most Abundant Species

#### Prepare the data for the graph

    #get the sum of the row 
    NormCounts_Treatment$sum <- apply(NormCounts_Treatment, 1, sum)

    #seperate the df by highly expressed bacteria
    #define what "high" expression threshold will be  
    threshold <- 10000
    threshNormCounts <- NormCounts_Treatment[(NormCounts_Treatment$sum) > threshold, ]

    #Now remove the "Sum" coulmn, as no longer needed 
    threshCountsNoSum <- threshNormCounts[, !names(threshNormCounts) %in% c("sum")]

    #add the rownames as a column 
    threshCountsNoSum <- tibble::rownames_to_column(threshCountsNoSum, "id")

    #change to long format 
    #gather by treatment type
    threshCounts_long <- threshCountsNoSum %>% tidyr::gather(sample, count, 2:ncol(threshCountsNoSum))

    #get rid of the mouse/replicate so that we can find mean, sd 
    threshCounts_long$sample = gsub("_A|_B|_C|_D|_E|_F", "", threshCounts_long$sample)

    #group now by ID as well as treatment type, and get mean
    DFplot <- threshCounts_long %>% group_by(sample, id) %>% summarise(mean = round(mean(count), 2))

    #group by sample and calculate relative abundance
    DFplot2 <- DFplot  %>% group_by(sample) %>% mutate(relAbundByPath = (100 * (mean / sum(mean))))

#### Stacked barplot

    #list of the days of the experiment 
    days <- c("day2.", "day4.", "day8.", "day11.", "day14.")

    #order the x axis by each treatment
    orderXaxis <- c(paste0(days, "pluronicBacteria"), paste0(days, "bacteria"), paste0(days, "pluronic"))

    #plot 
    bar <- ggplot(DFplot2, aes(x = sample, y = relAbundByPath, fill = id)) + 
      geom_bar(stat = "identity") + 
      xlab("") + ylab("Relative Abundance [%]") +
      theme(legend.title = element_blank()) +
      theme(legend.text = element_text(colour="black", size = 14, face = "plain")) +
      theme(axis.title.x = element_text(family="sans",size = 14, face="bold", hjust=0.5, vjust=-0.5),
            axis.title.y = element_text(family="sans",size = 14, angle=90, face="bold", hjust=0.5, vjust=1)) +
      theme(axis.text.x = element_text(family = "sans",size = 12, angle=60, face='plain', colour="#353535",   hjust=1, vjust=1) ) +
      theme(axis.text.y = element_text(family = "sans",size = 14, face='plain', colour="#353535",  vjust=0.5) ) +
      theme(axis.line.x = element_line(color="black", size = 0.5),
            axis.line.y = element_line(color="black", size = 0.5)) +
      theme(legend.background = element_rect()) + 
      theme(panel.border = element_blank(), panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),  panel.background = element_blank()) 

    #Specify the ordered x axis and manual fill colors 
    bar + scale_x_discrete(limits = orderXaxis) + 
      scale_fill_manual(values=c("#3a7751", "#d8732f", "#e00652", "#876ba0", "#0b00a5", 
                                 "#edc021", "#593f07", "#6d6d6c"))

![](GenusLevelFiguresForManuscript_files/figure-markdown_strict/unnamed-chunk-13-1.png)

### Heatmap of FC for Pluronic Bacteria

#### Clean and Organize the data

    #for the heatmap we want to access only the columns with FC data
    FC_cols <- c("log2FoldChange")

    #select all relevant columns
    onlyFC.df <- dplyr::select(DEresultsDF, matches(FC_cols))

    #clean up column names
    colnames(onlyFC.df) <- gsub('.*RUVseq_', '', colnames(onlyFC.df))

    #we are interested in displaying only the pluronic bacteria treatment
    pluronicBacteria <- onlyFC.df[ , grepl("pluronicBacteria", names(onlyFC.df))]

    #clean up the column names some more
    colnames(pluronicBacteria) <- gsub('.pluronicBacteria_vs_notreatment*', '', colnames(pluronicBacteria))

    #remove all the NAs, and turn to zero. 
    #Note:There are NAs bec this DF is initially from the total analysis
    pluronicBacteria[is.na(pluronicBacteria)] <- 0

    #check how many rows we have 
    nrow(pluronicBacteria)

    ## [1] 184

    #filter for all rows, with a sum of 0, this way there must be at least one change per row
    pluronicBacteria <- pluronicBacteria[which(rowSums(pluronicBacteria) != 0), ] 

    #lets remove any major FC outlier 
    logicalOutlier <- abs(rowSums(pluronicBacteria)) <= 20

    pluronicBacteria <- pluronicBacteria[logicalOutlier, ]

    #check how many rows there are again 
    nrow(pluronicBacteria)

    ## [1] 58

#### Heatmap

##### Load external function

    #load a function that was found here:https://github.com/raivokolde/pheatmap/issues/48
    #Used to make names bold 
    make_bold_names <- function(mat, rc_fun, rc_names) {
      bold_names <- rc_fun(mat)
      ids <- rc_names %>% match(rc_fun(mat))
      ids %>%
        walk(
          function(i)
            bold_names[i] <<-
            bquote(bold(.(rc_fun(mat)[i]))) %>%
            as.expression()
        )
      bold_names
    }

    #chose breakpoints after exmaining the range of the FC data 
    my_breaks <- c(-6.05, -5.05, -4.05, -3.05, -2.05, -1.05,  0.90, 1.35, 1.95,  2.95, 3.9, 4.9, 6.05)

    #chose custom colors for the heatmap 
    my_palette <- c(colorRampPalette(colors = c("darkblue", "lightblue"))(n = 5),
                    "white",
                    c(colorRampPalette(colors = c("tomato1", "darkred"))(n = 6)))

    #use grid to change the width and height of the plot so there is room for the label
    setHook("grid.newpage", function() pushViewport(viewport(x=1,y=1,width=0.95, height=0.95, name="vp", just=c("right","top"))), action="prepend")

    #draw the heatmap 
    heatmap <- pheatmap(pluronicBacteria,
                        clustering_method = "complete", 
                        clustering_distance_rows= "correlation",
                        border_color = "grey80",
                        cluster_cols = F,
                        cluster_rows = T,
                        show_rownames = T, 
                        show_colnames = T,
                        cellwidth = NA,
                        cellheight = NA,
                        angle_col = 45,
                        breaks = my_breaks,
                        treeheight_row = 0,
                        annotation_legend = F, 
                        color = my_palette,
                        labels_row = make_bold_names(pluronicBacteria, rownames, c("Bacillus"))) 

    #add a label to the x axis 
    setHook("grid.newpage", NULL, "replace")
    grid.text("Pluronic Bacteria", y=0, gp=gpar(fontsize=16))

![](GenusLevelFiguresForManuscript_files/figure-markdown_strict/unnamed-chunk-15-1.png)
