suppressPackageStartupMessages({
  library(dplyr)
  library(bedtoolsr)
  library(argparse)
  library(stringr)
})
options(scipen=999)
parsingArgs <- function(){
  parser <- ArgumentParser(prog = 'Rscript RegFeatAct.R',
                           description='Activity of regulatory features')
  parser$add_argument('currentAnnotation', metavar='currentAnnotation', 
                      type="character", default=NULL,
                      nargs=1, help='Full path to the gff3 file with
                      regulatory annotation.')
  parser$add_argument('species', metavar='S', type="character", nargs=1,
                      help='Species name in lowercase format and using
                      underscore (e.g xxx_yyy).')
  parser$add_argument('assembly', metavar='A', type="character", nargs=1,
                      help='Assembly name.')
  parser$add_argument('outPath', metavar='output', type="character", 
                      nargs=1, help = 'Full path to the output folder.')
  parser$add_argument('annotatedTSSs', metavar='TSSsFile', type="character", 
                      nargs=1, help='Full path to the bed file containing 
                      TSSs (chr, start, end, strand).')
  parser$add_argument('path2EpigenomesOCFiles', metavar='OCPFiles', type="character", 
                      nargs=1, help='Full path to the folder where 
                      atac-seq/dnase-seq peaks files are.')
  parser$add_argument('--epiPattern', dest='epiPattern', action='store',
                      default ='', help='Pattern characterizing the subset of 
                      atac-seq/dnase-seq bed files to be used.')
  parser$add_argument('--TSSwindowUp', dest='TSSwindowUp', action='store', 
                      default=90, help='Window size (bp) upstream TSSs for 
                      promoters identification (default 90).')
  parser$add_argument('--TSSwindowDown', dest='TSSwindowDown', 
                      action='store', 
                      default=10, help='Window size (bp) downstream TSSs for 
                      promoters identification (default 10).')
  parser$add_argument('--coreWUp', dest='coreWUp', action='store', 
                      default=490, help='Upstream core size (bp) for potential
                      promoters (default 490).')
  parser$add_argument('--coreWDown', dest='coreWDown', action='store', 
                      default=10, help='Downstream core size (bp) for potential
                      promoters (default 10).')
  parser$add_argument('--featureOverlap', dest='featureOverlap', 
                      action='store', default = 0.2, help='Fraction of overlap
                      between RFs and peaks (default 0.2).')
  return(parser)
}

loadEpiPeaks <- function(path2EpigenomesOCFiles, epiBedFiles, chrs, assembly){
  epigenomes <- do.call(rbind, strsplit(epiBedFiles, split = assembly))[,2]
  epigenomes <- do.call(rbind, strsplit(epigenomes, split = '-seq-'))[,2]
  epigenomes <- unique(do.call(rbind, strsplit(epigenomes, split = '_'))[,1])
  epiPeaks <- list()
  for(epi in epigenomes){
    bedFiles <- epiBedFiles[grepl(paste('seq-', epi, '_', sep =''), 
                                  epiBedFiles)]
    peaks <- NULL
    for(i in 1:length(bedFiles)){
      peaks %>%
        bind_rows(read.delim(paste(path2EpigenomesOCFiles, bedFiles[i],
                                   sep = "/"), header = F)) %>%
        bt.sort() %>%
        bt.merge() %>%
        filter(V1 %in% chrs) -> peaks
    }
    epiPeaks[[epi]] <- peaks
  }
  return(epiPeaks)
}
addingTSSsInfo <- function(annotatedTSSsAndInfo, TSSwindowUp, TSSwindowDown,
                           coreWUp, coreWDown){
  annotatedTSSsAndInfo %>% 
    select(V1, V2, V3, V4) %>%
    unique() -> annotatedTSSs
  names(annotatedTSSs) <- c('TSSChr', 'TSSStart', 'TSSEnd', 'TSSStrand')
  # for TSSs in '+' strand
  annotatedTSSs %>%
    mutate(wStart = if_else(TSSStrand == '+',
                            TSSStart - TSSwindowUp, TSSStart - TSSwindowDown ),
           wEnd= if_else(TSSStrand == '+',
                         TSSEnd + TSSwindowDown, TSSEnd + TSSwindowUp )) -> annotatedTSSs
  
  annotatedTSSs %>%
    mutate(cStart = if_else(TSSStrand == '+',
                            TSSStart - coreWUp, TSSStart - coreWDown ),
           cEnd = if_else(TSSStrand == '+',
                          TSSEnd + coreWDown, TSSEnd + coreWUp )) -> annotatedTSSs
  
  
  annotatedTSSs %>%
    mutate(across(c(wStart, cStart), 
                  function(x){if_else( x < 0, 0, x )})) -> annotatedTSSs
  
  annotatedTSSs %>%
    select(TSSChr, wStart, wEnd, TSSStrand, TSSStart, TSSEnd, cStart, cEnd) ->
    annotatedTSSs
  return(annotatedTSSs)
}
rfsInEpigenomes <- function(rfs, epiPeaks, annotatedTSSs, featureOverlap, eVal){
  # the count will be computed differently for promoters than for enhancers and
  #open chromatin regions
  # For the latter, we request minimum overlapping fraction (featureOverlap)
  
  # exclude EMARs
  rfs %>%
    filter(type != 'EMAR') -> rfs
  
  notPromRFs <- NULL
  for(epi in names(epiPeaks)){
    notPromRFs <- rbind(notPromRFs, 
                        cbind(bt.intersect(
                          a = rfs[rfs[,'type'] != 'promoter',
                                  c('chr', 'start', 'end', 'type')],
                          b = epiPeaks[[epi]], u = TRUE, wa = TRUE,
                          f = featureOverlap, F = featureOverlap,
                          e = eVal), V5=epi))
  }
  notPromRFs %>%
    mutate(V6 = paste(V4, V1, paste(V2, V3, sep = "-"), 
                      sep = ':')) -> notPromRFs
  
  # For promoters we request minimum overlapping. Because we extended UPs that
  # overlap with windows defined around TSSs,
  # it could happen we identified promoters but then they do not overlap 
  # with union peaks if they are overlapping two TSS ('-' and '+') that are
  # too close
  rfs %>%
    filter(type == 'promoter') %>%
    select(chr, start, end, type ) -> promsToOverlap
  
  annotatedTSSs %>%
    bt.intersect(promsToOverlap,
                 u = T, wa = T) -> annotatedTSSsInP
  
  promsInepis <- NULL
  for(epi in names(epiPeaks)){
    promsInEpi <- NULL
    ocp <- epiPeaks[[epi]]
    OCPInTSS <- bt.intersect(ocp, annotatedTSSsInP, wa = T, wb = T)
    if(nrow(OCPInTSS) > 0){
      OCPInTSS %>% 
        group_by(V1,V2,V3) %>% 
        summarise(chr=unique(V1), 
                  start=min(min(V5), V2), 
                  end=max(max(V6), V3), .groups = 'drop') %>% 
        as.data.frame() -> OCPInTSSExtended
      ov <- bt.intersect(a = promsToOverlap, 
                         b = OCPInTSSExtended[,c('chr', 'start', 'end')],
                         u = TRUE, wa = TRUE)
      if(nrow(ov) > 0){
        promsInEpi <- rbind(promsInEpi, cbind(ov, V5=epi)) 
      }
    }
    OCPNotInTSS <- bt.intersect(a = ocp, b = annotatedTSSsInP, wa = T, v=T)
    OCPInP <- bt.intersect(a = promsToOverlap, 
                           b = OCPNotInTSS,
                           u = TRUE, wa = TRUE)
    if(nrow(OCPInP) > 0){
      promsInEpi <- rbind(promsInEpi, 
                          cbind(OCPInP, V5=epi)) 
    }
    promsInepis <- rbind(promsInepis, unique(promsInEpi))
  } 
  promsInepis %>%
    mutate(V6 = paste(V4, V1, paste(V2, V3, sep = "-"),  sep = ':')) %>%
    bind_rows(notPromRFs) %>%
    select(-V4) %>%
    rename(chr=V1, start=V2, end=V3, epi=V5, loc=V6)-> rfsInEpiLong
  return(rfsInEpiLong)
}
exportEpigenomeActivity <- function(outPath, species, assembly, rfsGFF, 
                                    epigenomes, rfsInEpiLong){
  date <- paste0(strsplit(as.character(Sys.Date()), 
                          split = "-")[[1]], collapse = "")
  rfsGFF %>%
    pull(ID) -> stable_ID
  for(epi_i in epigenomes){
    rfsInEpiLong %>%
      filter(epi == epi_i) %>%
      pull(loc) -> rfInEpi
    
    actLev <- data.frame(stable_ID = stable_ID, status = 'NO EVIDENCE',
                         row.names = rownames(rfsGFF))
    actLev[rfInEpi, 'status'] <- 'EVIDENCE'
    write.csv(actLev, file = paste(outPath,
                                   paste(species, assembly, epi_i,
                                         "Regulatory_activity",
                                         date, "csv", sep = "."),
                                   sep = "/"),
              row.names = F, quote=F)
  }
}


main <- function(parser){
  args <- parser$parse_args()
  currentAnnotation <- args$currentAnnotation
  rfs <- tryCatch(read.delim(currentAnnotation, header =F, skip = 1),
                  warning = function(w) NULL)
  stopifnot("File with previous regulatory annotation doesn't exist" = 
              !is.null(currentAnnotation))
  names(rfs) <- c('chr', 'source', 'type', 'start', 'end', 'score', 
                     'strand', 'phase', 'attributes')
  rfs %>%
    select(chr, start, end, type, attributes) %>%
    mutate(start = start - 1,
           ID = str_split(str_split(attributes, 'ID=', simplify = T)[,2],
                          ';', simplify = T)[,1]) -> rfs 
  rfs %>%
    mutate(loc = paste(type, chr, paste(start, end, sep = '-'), 
                       sep = ':')) %>%
    pull(loc) -> rownames(rfs)
  # we only report activity for promoter, enhancer and open_chromatin_region
  #features
  rfs %>%
    filter(type %in% c('promoter', 'enhancer', 'open_chromatin_region')) -> rfs

  species <- args$species
  assembly <- args$assembly
  outPath <- args$outPath
  if(!file.exists(outPath)){
    system(paste('mkdir ', outPath, sep=''))
  }
  annotatedTSSsAndInfo <- tryCatch(read.delim(args$annotatedTSSs, header=F, 
                                              sep= "\t"), 
                                   warning = function(w) NULL) 
  stopifnot("File with annotated TSSs doesn't exist" = 
              !is.null(annotatedTSSsAndInfo))
  annotatedTSSsAndInfo %>%
    pull(V1) %>% unique() -> chrs
  
  path2EpigenomesOCFiles <- args$path2EpigenomesOCFiles
  
  epiPattern <- args$epiPattern
  TSSwindowUp <- as.numeric(args$TSSwindowUp)
  TSSwindowDown <- as.numeric(args$TSSwindowDown)
  coreWUp <- as.numeric(args$coreWUp)
  coreWDown <- as.numeric(args$coreWDown)
  featureOverlap <- args$featureOverlap
  eVal <- FALSE
  if(!is.null(featureOverlap)){
    featureOverlap <- as.numeric(featureOverlap)
    eVal <- TRUE 
  }
  # The activity of each regulatory features across all the epigenomes will be 
  # assigned following the next steps  
  
  # 1. Load all ATAC-seq/DNase-seq peaks.   
  epiBedFiles <- list.files(path = path2EpigenomesOCFiles, 
                            pattern = epiPattern)
  
  epiPeaks <- loadEpiPeaks(path2EpigenomesOCFiles, epiBedFiles, chrs, assembly)
  print(paste("Number of epigenomes: ", length(epiPeaks), sep = ""))
  
  # 2. Adding information of windows around TSSs in order to be able to compute
  # different overlapping for promoters
  annotatedTSSs <- addingTSSsInfo(annotatedTSSsAndInfo, TSSwindowUp,
                                  TSSwindowDown, coreWUp, coreWDown)
  # Computing RFs activity
  rfsInEpiLong <- rfsInEpigenomes(rfs, epiPeaks, annotatedTSSs, 
                                  featureOverlap, eVal)
  exportEpigenomeActivity(outPath, species, assembly, rfs, names(epiPeaks),
                          rfsInEpiLong)  
  
}

parser <- parsingArgs()
main(parser)
