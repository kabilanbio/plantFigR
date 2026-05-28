# ==========================================================================
# plantFigR_functions.R
# Gene‑peak correlations (DORC calling) and TF‑DORC associations for plants
# ==========================================================================

#' Gene-Peak correlations used for calling domains of regulatory chromatin (DORCs)
#'
#' Function to compute correlation between RNA expression and peak accessibility
#' for peaks falling within a window around each gene, across single cells,
#' using background peak correlations for significance testing. Works for any
#' plant species – provide the appropriate BSgenome and TSS GRanges.
#'
#' @param ATAC.se SummarizedExperiment object of the scATAC-seq reads in peak counts.
#' @param RNAmat Matrix object of the normalized scRNA-seq counts.
#' @param TSSg GRanges of TSS positions (1 bp each), named by gene ID.
#' @param bsgenome A BSgenome object for the target species (used for GC bias).
#' @param geneList Optional vector of gene symbols (subset of genes).
#' @param windowPadSize numeric, base pairs padded on either side of TSS (default 50000).
#' @param normalizeATACmat logical, whether to normalize ATAC counts (default TRUE).
#' @param nCores numeric, number of cores for parallelisation.
#' @param keepPosCorOnly logical, keep only positive correlations (default TRUE).
#' @param keepMultiMappingPeaks logical, keep multi‑mapping peaks (default FALSE).
#' @param n_bg numeric, number of background iterations (default 100).
#' @param p.cut numeric, p‑value cutoff for filtering results (default NULL).
#' @return data.frame of gene-peak correlations.
#' @export
runGenePeakcorr_plant <- function(
    ATAC.se,
    RNAmat,
    TSSg,
    bsgenome,
    geneList = NULL,
    windowPadSize = 50000,
    normalizeATACmat = TRUE,
    nCores = 4,
    keepPosCorOnly = TRUE,
    keepMultiMappingPeaks = FALSE,
    n_bg = 100,
    p.cut = NULL
) {

  stopifnot(inherits(ATAC.se, "RangedSummarizedExperiment"))
  stopifnot(inherits(RNAmat, c("Matrix", "matrix")))
  stopifnot(is(bsgenome, "BSgenome"))

  if (!all.equal(ncol(ATAC.se), ncol(RNAmat)))
    stop("Input ATAC and RNA objects must have same number of cells")

  message("Assuming paired scATAC/scRNA-seq data ..")

  peakRanges.OG <- granges(ATAC.se)

  rownames(ATAC.se) <- paste0("Peak", 1:nrow(ATAC.se))
  ATACmat <- assay(ATAC.se)

  if (normalizeATACmat)
    ATACmat <- centerCounts(ATACmat)

  if (is.null(rownames(RNAmat)))
    stop("RNA matrix must have gene names as rownames")

  # Remove peaks/genes with zero accessibility/expression
  if (any(Matrix::rowSums(assay(ATAC.se)) == 0)) {
    message("Peaks with 0 accessibility across cells exist ..")
    message("Removing these peaks prior to running correlations ..")
    peaksToKeep <- Matrix::rowSums(assay(ATAC.se)) != 0
    ATAC.se <- ATAC.se[peaksToKeep, ]
    ATACmat <- ATACmat[peaksToKeep, ]
    message("Important: peak indices in returned gene-peak maps are relative to original input SE")
  }

  peakRanges <- granges(ATAC.se)

  if (any(Matrix::rowSums(RNAmat) == 0)) {
    message("Genes with 0 expression across cells exist ..")
    message("Removing these genes prior to running correlations ..")
    genesToKeep <- Matrix::rowSums(RNAmat) != 0
    RNAmat <- RNAmat[genesToKeep, ]
  }

  cat("Number of peaks in ATAC data:", nrow(ATACmat), "\n")
  cat("Number of genes in RNA data:", nrow(RNAmat), "\n")

  stopifnot(is(TSSg, "GRanges"))
  names(TSSg) <- as.character(TSSg$gene_id)

  if (!is.null(geneList)) {
    if (length(geneList) == 1)
      stop("Please specify more than 1 valid gene symbol")
    if (any(!geneList %in% names(TSSg))) {
      cat("One or more of the gene names supplied is not present in the TSS annotation specified: \n")
      cat(geneList[!geneList %in% names(TSSg)], sep = ", ")
      cat("\n")
      stop()
    }
    TSSg <- TSSg[geneList]
  }

  genesToKeep <- intersect(names(TSSg), rownames(RNAmat))
  cat("\nNum genes overlapping TSS annotation and RNA matrix being considered: ",
      length(genesToKeep), "\n")

  RNAmat <- RNAmat[genesToKeep, ]
  TSSg <- TSSg[genesToKeep]

  TSSflank <- GenomicRanges::flank(TSSg, width = windowPadSize, both = TRUE)

  cat("\nTaking peak summits from peak windows ..\n")
  peakSummits <- resize(peakRanges, width = 1, fix = "center")

  cat("Finding overlapping peak-gene pairs ..\n")
  genePeakOv <- findOverlaps(query = TSSflank, subject = peakSummits)
  numPairs <- length(genePeakOv)
  cat("Found ", numPairs, "total gene-peak pairs for given TSS window ..\n")

  # GC bias correction
  if (is.null(rowData(ATAC.se)$bias)) {
    ATAC.se <- chromVAR::addGCBias(ATAC.se, genome = bsgenome)
  }

  # Background peaks
  set.seed(123)
  cat("Determining background peaks ..\n")
  if (any(Matrix::rowSums(assay(ATAC.se)) == 0)) {
    ATAC.mat <- assay(ATAC.se)
    ATAC.mat <- cbind(ATAC.mat, 1)
    ATAC.se.new <- SummarizedExperiment::SummarizedExperiment(
      assays   = list(counts = ATAC.mat),
      rowRanges = granges(ATAC.se)
    )
    bg <- chromVAR::getBackgroundPeaks(ATAC.se.new, niterations = n_bg)
  } else {
    bg <- chromVAR::getBackgroundPeaks(ATAC.se, niterations = n_bg)
  }

  cat("Computing gene-peak correlations ..\n")
  pairsPerChunk <- 500
  largeChunkSize <- 5000
  startingPoint <- 1
  chunkStarts <- seq(startingPoint, numPairs, largeChunkSize)
  chunkEnds <- chunkStarts + largeChunkSize - 1
  chunkEnds[length(chunkEnds)] <- numPairs

  dorcList <- list()
  for (i in seq_along(chunkStarts)) {
    cat("Running pairs: ", chunkStarts[i], " to ", chunkEnds[i], "\n")
    ObsCor <- FigR::PeakGeneCor(
      ATAC = ATACmat,
      RNA = RNAmat,
      OV = genePeakOv[chunkStarts[i]:chunkEnds[i]],
      chunkSize = pairsPerChunk,
      ncores = nCores,
      bg = bg
    )
    dorcList[[i]] <- ObsCor
    gc()
  }

  cat("\nMerging results ..\n")
  dorcTab <- bind_rows(dorcList)

  cat("Performing Z-test for correlation significance ..\n")
  permCols <- 4:(ncol(bg) + 3)

  if (keepPosCorOnly) {
    cat("Only considering positive correlations ..\n")
    dorcTab <- dorcTab %>% dplyr::filter(rObs > 0)
  }

  if (!keepMultiMappingPeaks) {
    cat("Keeping max correlation for multi-mapping peaks ..\n")
    dorcTab <- dorcTab %>% dplyr::group_by(Peak) %>% dplyr::filter(rObs == max(rObs))
  }

  dorcTab$Gene <- as.character(TSSg$gene_name)[dorcTab$Gene]
  dorcTab$Peak <- as.numeric(splitAndFetch(rownames(ATACmat)[dorcTab$Peak], "Peak", 2))

  dorcTab$rBgSD <- matrixStats::rowSds(as.matrix(dorcTab[, permCols]))
  dorcTab$rBgMean <- rowMeans(dorcTab[, permCols])
  dorcTab$pvalZ <- 1 - stats::pnorm(
    q = dorcTab$rObs,
    mean = dorcTab$rBgMean,
    sd = dorcTab$rBgSD
  )

  cat("\nFinished!\n")

  if (!is.null(p.cut)) {
    cat("Using significance cut-off of ", p.cut, " to subset to resulting associations\n")
    dorcTab <- dorcTab[dorcTab$pvalZ <= p.cut, ]
  }

  dorcTab$PeakRanges <- paste(
    as.character(seqnames(peakRanges.OG[dorcTab$Peak])),
    paste(
      start(peakRanges.OG[dorcTab$Peak]),
      end(peakRanges.OG[dorcTab$Peak]),
      sep = "-"
    ),
    sep = ":"
  )

  return(as.data.frame(
    dorcTab[, c("Peak", "PeakRanges", "Gene", "rObs", "pvalZ")],
    stringsAsFactors = FALSE
  ))
}


#' Infer FigR TF-DORC associations for any plant species
#'
#' Function to run TF motif-to-gene associations using reference DORC peak-gene
#' mappings and TF RNA expression levels. Works for any plant species – provide
#' the appropriate BSgenome and PWM library.
#'
#' @param ATAC.se SummarizedExperiment object of scATAC peak counts.
#' @param dorcK numeric, number of DORC kNNs to pool (default 30).
#' @param dorcTab data.frame of peak-gene associations (from runGenePeakcorr_plant).
#' @param n_bg numeric, number of background peaks for motif enrichment (default 50).
#' @param bsgenome A BSgenome object for the target species.
#' @param pwm A PWMatrixList of TF motifs.
#' @param dorcMat Smoothed DORC matrix (genes x cells).
#' @param rnaMat Smoothed RNA matrix (genes x cells).
#' @param dorcGenes Optional character vector of DORC genes to test.
#' @param nCores numeric, number of cores.
#' @return data.frame with TF-DORC motif enrichment and correlation associations.
#' @export
runFigRGRN_plant <- function(
    ATAC.se,
    dorcK = 30,
    dorcTab,
    n_bg = 50,
    bsgenome,
    pwm,
    dorcMat,
    rnaMat,
    dorcGenes = NULL,
    nCores = 4
) {
  stopifnot(all.equal(ncol(dorcMat), ncol(rnaMat)))

  if (!all(c("Peak", "Gene") %in% colnames(dorcTab)))
    stop("Expecting fields Peak and Gene in dorcTab data.frame")

  if (all(grepl("chr", dorcTab$Peak, ignore.case = TRUE))) {
    usePeakNames <- TRUE
    message("Detected peak region names in Peak field")
    if (!all(dorcTab$Peak %in% rownames(ATAC.se)))
      stop("Peak regions in dorcTab not found in input SE")
  } else {
    usePeakNames <- FALSE
    message("Assuming peak indices in Peak field")
    if (max(dorcTab$Peak) > nrow(ATAC.se))
      stop("Found DORC peak index outside range of input SE")
  }

  if (is.null(dorcGenes)) {
    dorcGenes <- rownames(dorcMat)
  } else {
    cat("Using specified list of dorc genes ..\n")
    if (!all(dorcGenes %in% rownames(dorcMat))) {
      cat("One or more of the gene names supplied is not present in the DORC matrix provided: \n")
      cat(dorcGenes[!dorcGenes %in% rownames(dorcMat)], sep = ", ")
      cat("\n")
      stop()
    }
  }

  DORC.knn <- FNN::get.knn(data = t(scale(Matrix::t(dorcMat))), k = dorcK)$nn.index
  rownames(DORC.knn) <- rownames(dorcMat)

  if (is.null(rowData(ATAC.se)$bias)) {
    ATAC.se <- chromVAR::addGCBias(ATAC.se, genome = bsgenome)
  }

  if (missing(pwm) || !is(pwm, "PWMatrixList")) {
    stop("Please provide a PWMatrixList of TF motifs via the 'pwm' argument.")
  }

  if (all(grepl("_", names(pwm), fixed = TRUE)))
    names(pwm) <- FigR::extractTFNames(names(pwm))

  message("Removing genes with 0 expression across cells ..\n")
  rnaMat <- rnaMat[Matrix::rowSums(rnaMat) != 0, ]
  myGeneNames <- gsub("-", "", rownames(rnaMat))
  rownames(rnaMat) <- myGeneNames

  motifsToKeep <- intersect(names(pwm), myGeneNames)

  cat("Getting peak x motif matches ..\n")
  motif_ix <- motifmatchr::matchMotifs(
    subject = ATAC.se,
    pwms = pwm[motifsToKeep],
    genome = bsgenome
  )
  motif_ix <- motif_ix[, Matrix::colSums(assay(motif_ix)) != 0]

  cat("Determining background peaks ..\n")
  cat("Using ", n_bg, " iterations ..\n\n")
  if (any(Matrix::rowSums(assay(ATAC.se)) == 0)) {
    ATAC.mat <- assay(ATAC.se)
    ATAC.mat <- cbind(ATAC.mat, 1)
    ATAC.se.new <- SummarizedExperiment::SummarizedExperiment(
      assays = list(counts = ATAC.mat),
      rowRanges = granges(ATAC.se)
    )
    set.seed(123)
    bg <- chromVAR::getBackgroundPeaks(ATAC.se.new, niterations = n_bg)
  } else {
    set.seed(123)
    bg <- chromVAR::getBackgroundPeaks(ATAC.se, niterations = n_bg)
  }

  cat("Testing ", length(motifsToKeep), " TFs\n")
  cat("Testing ", nrow(dorcMat), " DORCs\n")
  if (nCores > 1)
    message("Running FigR using ", nCores, " cores ..\n")

  opts <- list()
  pb <- txtProgressBar(min = 0, max = length(dorcGenes), style = 3)
  progress <- function(n) setTxtProgressBar(pb, n)
  opts <- list(progress = progress)

  cl <- parallel::makeCluster(nCores)
  clusterEvalQ(cl, .libPaths())
  doSNOW::registerDoSNOW(cl)

  mZtest.list <- foreach(
    g = dorcGenes,
    .options.snow = opts,
    .packages = c("FigR", "dplyr", "Matrix", "Rmpfr")
  ) %dopar% {
    DORCNNpeaks <- unique(dorcTab$Peak[dorcTab$Gene %in%
                                         c(g, rownames(dorcMat)[DORC.knn[g, ]])])

    if (usePeakNames)
      DORCNNpeaks <- which(rownames(ATAC.se) %in% DORCNNpeaks)

    mZ <- FigR::motifPeakZtest(
      peakSet = DORCNNpeaks,
      bgPeaks = bg,
      tfMat = assay(motif_ix)
    )
    mZ <- mZ[, c("gene", "z_test")]
    colnames(mZ)[1] <- "Motif"
    colnames(mZ)[2] <- "Enrichment.Z"
    mZ$Enrichment.P <- 2 * pnorm(abs(mZ$Enrichment.Z), lower.tail = FALSE)
    mZ$Enrichment.log10P <- sign(mZ$Enrichment.Z) * -log10(mZ$Enrichment.P)
    mZ <- cbind("DORC" = g, mZ)

    corr.r <- cor(dorcMat[g, ], t(as.matrix(rnaMat[mZ$Motif, ])),
                  method = "spearman")
    stopifnot(all.equal(colnames(corr.r), mZ$Motif))

    mZ$Corr <- corr.r[1, ]
    mZ$Corr.Z <- scale(mZ$Corr, center = TRUE, scale = TRUE)[, 1]
    mZ$Corr.P <- 2 * pnorm(abs(mZ$Corr.Z), lower.tail = FALSE)
    mZ$Corr.log10P <- sign(mZ$Corr.Z) * -log10(mZ$Corr.P)
    return(mZ)
  }

  cat("Finished!\n")
  cat("Merging results ..\n")

  TFenrich.d <- do.call('rbind', mZtest.list)
  rownames(TFenrich.d) <- NULL

  TFenrich.d <- TFenrich.d %>%
    dplyr::mutate(
      "Score" = sign(Corr) *
        as.numeric(-log10(1 - (1 - Rmpfr::mpfr(Enrichment.P, 100)) *
                            (1 - Rmpfr::mpfr(Corr.P, 100))))
    )
  TFenrich.d$Score[TFenrich.d$Enrichment.Z < 0] <- 0
  TFenrich.d
}
