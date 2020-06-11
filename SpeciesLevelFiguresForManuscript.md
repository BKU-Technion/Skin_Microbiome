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
    library(tidyr)

The workspace includes the abundance data from two seperate sequencing
experiments as well as the design of the experiment as is described in
the manuscript.

The abundance matrix used here excludes Day 6.

### Load input data

    #The abundance data 
    load(file = paste0(getwd(), "/SpeciesMergedDesignFirstSecond_NoSixthDayNoNew.RData"))

    #The experimental design 
    load(file = paste0(getwd(), "/SpeciesMergedDataFirstSecond_NoSixthDayNoNew.RData"))

#### Filter the data:

    #count the number of rows are in the final matrix
    nrow(FinalMergedData)

    ## [1] 2014

    #we decided to filter for 50 reads in at least 2 samples
    filter <- apply(FinalMergedData, 1, function(x) length(x[x > 50]) >= 2)

    #subset the raw abundance data 
    filteredSpeciesLevel <- FinalMergedData[filter, ]

    #how many rows remain in the matrix?
    nrow(filteredSpeciesLevel)

    ## [1] 378

    #for downstream processes convert the input data into an expression set 
    #create an expression set using the assay and pheno data
    ExpressionSet <-EDASeq::newSeqExpressionSet(as.matrix(filteredSpeciesLevel), 
                                         phenoData = mergedDesign, 
                                         row.names = colnames(filteredSpeciesLevel))

### RUVseq: Estimating the factors of unwanted variation using replicate samples

    #construct a matrix that specifies the replicate samples 
    #we will use the ruvseq column for this, each row indicates a set of replicate samples 
    differences <- RUVSeq::makeGroups(mergedDesign$RUVseq) 

    #Run ruvseq with the differences and the number of desired k's 
    ExpressionSet2 <- RUVSeq::RUVs(ExpressionSet, rownames(filteredSpeciesLevel), k = 15, differences)

    #examine all data now included in the phenotypic data 
    #pData(ExpressionSet2)

### PCA

    #PCA divided by treatment groups 
    TreatmentColors <- c("#000000", rep(c("#250885", "#850816","#035940"), 5))

    #Convert the Colors into a numeric vector which matches the color assignment 
    TreatmentColorsNumeric <- TreatmentColors[as.numeric(mergedDesign$RUVseq)]

    #Define the shape 
    TreatmentShapes <- c(rep(16, 16))

    #Convert the Treatment into a numeric vector which matches the shape assignment 
    TreatmentShapesNumeric <- TreatmentShapes[as.numeric(mergedDesign$RUVseq)]

    #The margins of the plot 
    par(mar=c(11, 4.1, 2, 5.1), xpd=TRUE)

    #the letter k here represents which PC's to display 
    pca_data <- EDASeq::plotPCA(ExpressionSet2, col=TreatmentColorsNumeric, k = 2, 
                                cex=2.5, labels = F, pch=TreatmentShapesNumeric)
    #this is the input into legend... 
    legend("bottom",inset=c(0, -0.8), legend=levels(mergedDesign$RUVseq), pch=TreatmentShapes, 
           title="Treatment Groups",cex=0.8, ncol = 3, col=TreatmentColors)

![](SpeciesLevelFiguresForManuscript_files/figure-markdown_strict/unnamed-chunk-5-1.png)

#### Differential expression with DESeq2

    ## estimating size factors

    ## estimating dispersions

    ## gene-wise dispersion estimates

    ## mean-dispersion relationship

    ## final dispersion estimates

    ## fitting model and testing

    ## 133 rows did not converge in beta, labelled in mcols(object)$betaConv. Use larger maxit argument with nbinomWaldTest

### Manipulate final results DF

    #retrieve the row that reflect the contrasts of b. subtilis
    Bacillus_ResData <- DEresultsDF[rownames(DEresultsDF) == "Bacillus subtilis", ]

    #retrieve only the padj values for all contrasts for Bacillus subtilis
    BacillusPadj.df <- as.data.frame(t(dplyr::select(Bacillus_ResData, matches("padj"))))

    #change the name to padj
    colnames(BacillusPadj.df) <- "padj"

    #add an additional column that defines if significant 
    BacillusPadj.df$significant <- ifelse(BacillusPadj.df$padj < 0.01, "sig", "non-sig") 

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

    #remove the mouse/replicate so that we can find mean
    counts_long$sample = gsub("_A|_B|_C|_D|_E|_F", "", counts_long$sample)

    #group now by ID as well as treatment type, and get mean
    new.df <- counts_long %>% group_by(id, sample) %>% summarize(avg = mean(count))

    ## `summarise()` regrouping output by 'id' (override with `.groups` argument)

    # extract the data for a specific genus
    ind = grep("Bacillus subtilis", new.df$id) 

    df = new.df[ind, ]

    #merge the count data with the padj and significance 
    df.pvalue <- merge(df, BacillusPadj.df, by.x = "sample", by.y = "treatment")

    #add the following columns in order to graph the seperate treatments 
    #add column with the day and relvel 
    df.pvalue$day <- sapply(strsplit(as.character(df.pvalue$sample), ".", fixed = TRUE), '[', 1)
    df.pvalue$day <-factor(df.pvalue$day, levels=c("day2", "day4", "day8", "day11", "day14"))

    #add column with the treament info 
    df.pvalue$treatment <- sapply(strsplit(as.character(df.pvalue$sample), ".", fixed = TRUE), '[', 2)

### Graph *B. subtillus* counts

    #because there are many points near each other
    pd <- position_dodge(0.05) # move them .05 to the left and right

    #plot data
    df.pvalue.plot <- ggplot(df.pvalue, aes(x=day, y=avg)) + 
      geom_line(position = pd, size =1, aes(group = treatment, color = treatment)) +
      geom_point(position = pd, size = 5, aes(shape = significant)) + 
      ylab("Bacillus subtilis Normalized Counts") + xlab("") + 
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

![](SpeciesLevelFiguresForManuscript_files/figure-markdown_strict/unnamed-chunk-11-1.png)
