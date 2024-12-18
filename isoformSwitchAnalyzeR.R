library(IsoformSwitchAnalyzeR)

salmonQuant <- importIsoformExpression(
  parentDir = "/lustre/alice3/scratch/vasccell/cs806/brbSeq/testSeq/salmonQuant"
)

head(salmonQuant$abundance, 2)

ids <- colnames(salmonQuant$counts[-1])
myDesign <- data.frame(sampleID = ids,
                       condition = ifelse(str_detect(ids, "untr"), "UNTR", "TNFa"))

### Create switchAnalyzeRlist
huvecSwitchList <- importRdata(
  isoformCountMatrix   = salmonQuant$counts,
  isoformRepExpression = salmonQuant$abundance,
  designMatrix         = myDesign,
  isoformExonAnnoation = "/scratch/vasccell/cs806/exprPhenoData/gencode.v44.primary_assembly.annotation.gtf.gz",
  isoformNtFasta       = "/scratch/vasccell/cs806/exprPhenoData/gencode.v44.transcripts.fa.gz",
  fixStringTieAnnotationProblem = FALSE,
  showProgress = FALSE,
  removeNonConvensionalChr = TRUE
)

summary(huvecSwitchList)

huvecSwitchListFiltered <- preFilter(
  switchAnalyzeRlist = huvecSwitchList,
  geneExpressionCutoff = nrow(huvecSwitchList$designMatrix)/2,
  isoformExpressionCutoff = 3,
  removeSingleIsoformGenes = TRUE
)


huvecSwitchListAnalRes <- isoformSwitchAnalysisPart1(
  switchAnalyzeRlist   = huvecSwitchListFiltered,
  pathToOutput = "/scratch/vasccell/cs806/brbSeq/testSeq/isaFastaOutput",
  outputSequences      = TRUE, # change to TRUE whan analyzing your own data 
  prepareForWebServers = TRUE  # change to TRUE if you will use webservers for external sequence analysis
)

extractSwitchSummary(huvecSwitchList)


# Testing for Isoform Switches via DEXSeq

# Perform test
huvecSwitchListDexSeq <- isoformSwitchTestDEXSeq(
  switchAnalyzeRlist = huvecSwitchListFiltered,
  reduceToSwitchingGenes=TRUE
)

extractSwitchSummary(huvecSwitchListDexSeq)

dexSeqRes <- huvecSwitchListDexSeq$isoformSwitchAnalysis
annot <- huvecSwitchList$isoformFeatures

dexSeqResAnnot <- merge(annot[, c(3,4,7)], dexSeqRes[, -c(1,2)], by = "isoform_id")

write.csv(dexSeqResAnnot, "testSeq/tables/Differential_Exon_usage.csv", row.names = F)

head(annot)
head(dexSeqRes)
#######################################################

# Perform test
huvecSwitchListSatuRn <- isoformSwitchTestSatuRn(
  switchAnalyzeRlist = huvecSwitchListFiltered,
  reduceToSwitchingGenes=TRUE
)

extractSwitchSummary(huvecSwitchListSatuRn)

SatuRnRes <- huvecSwitchListSatuRn$isoformSwitchAnalysis
annot <- huvecSwitchList$isoformFeatures

SatuRnResAnnot <- merge(annot[, c(3,4,7)], SatuRnRes[, -c(1,2)], by = "isoform_id")

write.csv(SatuRnResAnnot, "testSeq/tables/Differential_Exon_usage.csv", row.names = F)

head(annot)
head(SatuRnRes)




# External sites
# Best to upload files instead of copy and paste
# https://cpc2.gao-lab.org/
# https://www.ebi.ac.uk/Tools/hmmer/search/hmmscan
# https://iupred2a.elte.hu/                             
# http://www.cbs.dtu.dk/services/SignalP/  https://services.healthtech.dtu.dk/services/SignalP-5.0/



huvecSwitchList2 <- isoformSwitchAnalysisPart2(
  switchAnalyzeRlist        = huvecSwitchList, 
  n                         = 5,    # if plotting was enabled, it would only output the top 10 switches
  removeNoncodinORFs        = TRUE,
  pathToCPC2resultFile      = "testSeq/isaFastaOutput/cpc2_result.txt",
  pathToPFAMresultFile      = "testSeq/isaFastaOutput/pfam_results.txt",
  pathToIUPred2AresultFile  = "testSeq/isaFastaOutput/iupred2a_result.txt.gz",
  pathToSignalPresultFile   = "testSeq/isaFastaOutput/signalP_results.txt",
  outputPlots               = FALSE  # keeps the function from outputting the plots from this example
)




#################################################################################

data("exampleSwitchList")




exampleSwitchList <- isoformSwitchAnalysisPart2(
  switchAnalyzeRlist        = exampleSwitchList, 
  n                         = 10,    # if plotting was enabled, it would only output the top 10 switches
  removeNoncodinORFs        = TRUE,
  pathToCPC2resultFile      = system.file("extdata/cpc2_result.txt"         , package = "IsoformSwitchAnalyzeR"),
  pathToPFAMresultFile      = system.file("extdata/pfam_results.txt"        , package = "IsoformSwitchAnalyzeR"),
  pathToIUPred2AresultFile  = system.file("extdata/iupred2a_result.txt.gz"  , package = "IsoformSwitchAnalyzeR"),
  pathToSignalPresultFile   = system.file("extdata/signalP_results.txt"     , package = "IsoformSwitchAnalyzeR"),
  outputPlots               = FALSE  # keeps the function from outputting the plots from this example
)
