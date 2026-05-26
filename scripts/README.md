# regulationRegBuild

This repository contains a set of scripts used in the Regulatory Annotation Pipeline of the Ensembl Regulation Team, EMBL-EBI.

The main script is [RegulatoryBuild.R](./RegulatoryBuild.R), aimed to identify the set of Regulatory Features across the whole genome of a determined species (`S`, based on assembly `A`). For this, a set of open chromatin peaks (obtained from ATAC-seq or DNase-seq experiments), H3K4me1, H3K27ac, and H3K4me3 peaks are used. In addition, TSSs and CDSs annotations are required to the identification of promoters and enhancers, respectively.

The strategy for identifying regulatory elements is detailed [here](https://regulation.ensembl.org/help/regulatory_build).

Required and optional parameters can be checked by using the `-h` option:

```         
Rscript RegulatoryBuild.R -h
usage: Rscript RegulatoryBuild.R [-h] [--path2K4me1Files PATH2K4ME1FILES]
                                 [--path2K27acFiles PATH2K27ACFILES]
                                 [--path2K4me3Files PATH2K4ME3FILES]
                                 [--oldAnnotation OLDANNOTATION]
                                 [--ocPattern OCPATTERN]
                                 [--k4me1Pattern K4ME1PATTERN]
                                 [--k27acPattern K27ACPATTERN]
                                 [--k4me3Pattern K4ME3PATTERN]
                                 [--TSSwindowUp TSSWINDOWUP]
                                 [--TSSwindowDown TSSWINDOWDOWN]
                                 [--coreWUp COREWUP] [--coreWDown COREWDOWN]
                                 [--featureOverlap FEATUREOVERLAP]
                                 [--histoneOverlap HISTONEOVERLAP]
                                 [--OCRsHistoneOverlap OCRSHISTONEOVERLAP]
                                 [--cOCRsCDSOverlap COCRSCDSOVERLAP]
                                 [--maxDistOCR MAXDISTOCR]
                                 [--maxDistH MAXDISTH]
                                 [--maxPromExtLen MAXPROMEXTLEN]
                                 [--minOCRLen MINOCRLEN]
                                 [--maxEnhancerLen MAXENHANCERLEN]
                                 [--maxOCRLen MAXOCRLEN]
                                 [--exportActivity EXPORTACTIVITY]
                                 [--path2EpigenomesOCFiles EPIOCFILES]
                                 [--epiPattern EPIPATTERN]
                                 S A EPrefix startID TSSsFile CDSFile output
                                 OCPFiles

Regulatory features annotation

positional arguments:
  S                     Species name in lower case and using underscore as a
                        separator (e.g Xxx_yyy).
  A                     Assembly name.
  EPrefix               Ensembl prefix for the species.
  startID               Integer specifying the starting value for feature IDs.
  TSSsFile              Full path to the bed file containing TSSs annotation
                        (chr, start, end, strand).
  CDSFile               Full path to the bed file containing CDSs annotation
                        (chr, start, end).
  output                Full path to the output folder.
  OCPFiles              Full path to the folder where atac-seq/dnase-seq peak
                        files are.

options:
  -h, --help            show this help message and exit
  --path2K4me1Files PATH2K4ME1FILES
                        Full path to the folder where H3K4me1 ChiP-seq peaks
                        files are.
  --path2K27acFiles PATH2K27ACFILES
                        Full path to the folder where H3K27ac ChiP-seq peaks
                        files are.
  --path2K4me3Files PATH2K4ME3FILES
                        Full path to the folder where H3K4me3 ChIP-seq peaks
                        files are.
  --oldAnnotation OLDANNOTATION
                        Full path to the gff3 file with previous regulatory
                        annotation.
  --ocPattern OCPATTERN
                        Pattern characterizing the subset of atac-seq/dnase-
                        seq bed files to be used.
  --k4me1Pattern K4ME1PATTERN
                        Pattern characterizing the subset of H3K4me1 bed files
                        to be used.
  --k27acPattern K27ACPATTERN
                        Pattern characterizing the subset of H3K27ac bed files
                        to be used.
  --k4me3Pattern K4ME3PATTERN
                        Pattern characterizing the subset of H3K4me3 bed files
                        to be used.
  --TSSwindowUp TSSWINDOWUP
                        Window size (bp) upstream TSSs for promoters
                        identification (default 90).
  --TSSwindowDown TSSWINDOWDOWN
                        Window size (bp) downstream TSSs for promoters
                        identification (default 10).
  --coreWUp COREWUP     Upstream core size (bp) for potential promoters
                        (default 490).
  --coreWDown COREWDOWN
                        Downstream core size (bp) for potential promoters
                        (default 10).
  --featureOverlap FEATUREOVERLAP
                        Fraction of overlapping between RFs and open chromatin
                        peaks for activity annotation (default 0.2).
  --histoneOverlap HISTONEOVERLAP
                        Fraction of overlap between histone and open chromatin
                        peaks, for features identification (default 0.5).
  --OCRsHistoneOverlap OCRSHISTONEOVERLAP
                        Fraction of overlap between cOCRs and histone peaks
                        for enhancers identification (default 0.5).
  --cOCRsCDSOverlap COCRSCDSOVERLAP
                        Fraction of overlap between cOCRs and CDSs for
                        enhancers identification (default 0.1).
  --maxDistOCR MAXDISTOCR
                        Maximum distance between peaks to be merged when
                        defining cOCRs (default 100).
  --maxDistH MAXDISTH   Maximum distance between histone peaks to be merged
                        (default 100).
  --maxPromExtLen MAXPROMEXTLEN
                        Maximum size for promoters extended regions (default
                        1000).
  --minOCRLen MINOCRLEN
                        Minimum size for cOCRs (default 100).
  --maxEnhancerLen MAXENHANCERLEN
                        Maximum length for enhancers (default 2500).
  --maxOCRLen MAXOCRLEN
                        Maximum length for open chromatin regions (default
                        2500).
  --exportActivity EXPORTACTIVITY
                        Should activity files be generated?(default FALSE)
  --path2EpigenomesOCFiles EPIOCFILES
                        Full path to the folder containing the open chromatin
                        peak files to be used for activity annotation.
  --epiPattern EPIPATTERN
                        Pattern characterizing the subset of bed files to be
                        used for activity annotation.
```

There are helper scripts that complement the regulatory features definition.

The [RegFeatAct.R.R](./RegFeatAct.R) is aimed to predict RFs activity in a set of epigenomes using a pre-existent set of RFs.

Required and optional parameters can be checked by using the `-h` option:

```
Rscript RegFeatAct.R -h
usage: Rscript RegFeatAct.R [-h] [--epiPattern EPIPATTERN]
                            [--TSSwindowUp TSSWINDOWUP]
                            [--TSSwindowDown TSSWINDOWDOWN]
                            [--coreWUp COREWUP] [--coreWDown COREWDOWN]
                            [--featureOverlap FEATUREOVERLAP]
                            currentAnnotation S A output TSSsFile OCPFiles

Activity of regulatory features

positional arguments:
  currentAnnotation     Full path to the gff3 file with regulatory annotation.
  S                     Species name in lowercase format and using underscore
                        (e.g xxx_yyy).
  A                     Assembly name.
  output                Full path to the output folder.
  TSSsFile              Full path to the bed file containing TSSs (chr, start,
                        end, strand).
  OCPFiles              Full path to the folder where atac-seq/dnase-seq peaks
                        files are.

options:
  -h, --help            show this help message and exit
  --epiPattern EPIPATTERN
                        Pattern characterizing the subset of atac-seq/dnase-
                        seq bed files to be used.
  --TSSwindowUp TSSWINDOWUP
                        Window size (bp) upstream TSSs for promoters
                        identification (default 90).
  --TSSwindowDown TSSWINDOWDOWN
                        Window size (bp) downstream TSSs for promoters
                        identification (default 10).
  --coreWUp COREWUP     Upstream core size (bp) for potential promoters
                        (default 490).
  --coreWDown COREWDOWN
                        Downstream core size (bp) for potential promoters
                        (default 10).
  --featureOverlap FEATUREOVERLAP
                        Fraction of overlap between RFs and peaks (default
                        0.2).

```

The [AnnotationExtractor.R](./AnnotationExtractor.R) is useful to extract TSSs and CDSs annotation from Ensembl GFF3 files.

Usage:

```
Rscript AnnotationExtractor.R -h
usage: Rscript AnnotationExtractor.R [-h]
                                     [--transcriptsBiotypes TRANSCRIPTSBIOTYPES]
                                     [--flagsToKeep FLAGSTOKEEP]
                                     [--flagsToRemove FLAGSTOREMOVE]
                                     ChrFile annotationGFF3 isGENCODE

Extract TSSs and CDSs annotation from Ensembl GFF3 file

positional arguments:
  ChrFile               Full path to the file listing the canonical
                        chromosomes / chromosomes of interest (one per line).
  annotationGFF3        Full path to the Ensembl gff3 file containing gene
                        annotation.
  isGENCODE             Logical indicating if the genome annotation file has
                        been produced by GENCODE.

options:
  -h, --help            show this help message and exit
  --transcriptsBiotypes TRANSCRIPTSBIOTYPES
                        Comma-separated list of transcript biotypes to keep.
                        Default value is `NULL` indicating that all biotypes
                        should be kept.
  --flagsToKeep FLAGSTOKEEP
                        Comma-separated list of transcript flags to be keep.
                        Default value is `NULL` indicating that all
                        transcripts with any flag are kept.
  --flagsToRemove FLAGSTOREMOVE
                        Comma-separated list of flags indicating transcripts
                        having them should be removed. Default value is NULL,
                        all transcripts are kept.
```