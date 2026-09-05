## =============================================================================
##' @title Preprocessing JAR and Liking Data for CLUSCATA-liking
##'
##' @description
##' Preprocesses Just-About-Right (JAR) and liking data for use with
##' CLUSCATA-liking. JAR responses are converted into binary data, where
##' 1 indicates a JAR response and 0 indicates a non-JAR response.
##' The function also reshapes the JAR and liking data into the formats
##' required by CLUSCATA-liking.
##'
##' @usage
##' preprocess_JAR_liking(Data, nprod, nsub, liking_col,
##'                       levelsJAR = 3, scale = FALSE)
##'
##' @param Data A data frame where the first column contains the subjects,
##' the second column contains the products, and the remaining columns contain
##' the JAR attributes and the liking variable.
##'
##' @param nprod Integer. Number of products.
##'
##' @param nsub Integer. Number of subjects.
##'
##' @param liking_col Character. Name of the column containing the liking
##' scores.
##'
##' @param levelsJAR Integer. Number of levels of the JAR scale. Must be
##' either 3 or 5. For a 3-level scale, level 2 is considered JAR.
##' For a 5-level scale, level 3 is considered JAR. All other levels are
##' considered non-JAR.
##'
##' @param scale Logical. Should the liking data be scaled when combining
##' JAR and liking data? Default is FALSE.
##'
##' @return A list containing:
##' \itemize{
##'   \item \code{Datafinal}: combined JAR and liking data ready for
##'   CLUSCATA-liking.
##'   \item \code{CATA}: binary JAR data in CLUSCATA format, with
##'   1 = JAR and 0 = non-JAR.
##'   \item \code{liking}: liking matrix with products in rows and subjects
##'   in columns.
##'   \item \code{JAR_binary}: binary JAR data before conversion to
##'   CLUSCATA format.
##'   \item \code{NameSub}: subject names in the order used in the analysis.
##'   \item \code{NameProd}: product names in the order used in the analysis.
##' }
##'
##' @keywords JAR liking clustering
##'
##' @references
##' Llobell, F. & Guksch, T. (2026). Beyond penalty analysis:
##' Joint clustering of consumers using JAR and liking data.
##' EuroSense, Oslo, Norway.
##'
##' @examples
##' data(croissant)
##'
##' # Use a subset of 40 subjects for a faster example
##' subjects <- unique(croissant[[1]])[1:40]
##' croissant40 <- croissant[croissant[[1]] %in% subjects, ]
##'
##' prepro <- preprocess_JAR_liking(
##'   Data = croissant40,
##'   nprod = 6,
##'   nsub = 40,
##'   liking_col = "OVL_Overall Liking",
##'   levelsJAR = 5,
##'   scale = FALSE
##' )
##'
##' CATA <- prepro$CATA
##' liking <- prepro$liking
##' Data <- prepro$Datafinal
##'
##' res <- cluscata_liking(
##'   Data,
##'   nblo = 40,
##'   NameBlocks = prepro$NameSub,
##'   printlevel = FALSE
##' )
##'
##' summary(res)
##'
##' plot(
##'   res,
##'   prepro$CATA,
##'   prepro$liking,
##'   scale = FALSE
##' )
##'
##' @seealso \code{\link{cluscata_liking}}
##'
##' @export

preprocess_JAR_liking <- function(Data,
                                  nprod,
                                  nsub,
                                  liking_col,
                                  levelsJAR = 3,
                                  scale = FALSE) {

  Data <- as.data.frame(Data)

  # ------------------------------------------------------------
  # Checks
  # ------------------------------------------------------------

  if (!levelsJAR %in% c(3, 5)) {
    stop("levelsJAR must be 3 or 5.")
  }

  if (!liking_col %in% colnames(Data)) {
    stop("The liking variable was not found in Data.")
  }

  NameSub  <- unique(Data[[1]])
  NameProd <- unique(Data[[2]])

  if (length(NameSub) != nsub) {
    stop(
      "nsub does not correspond to the number of subjects. Found ",
      length(NameSub), "."
    )
  }

  if (length(NameProd) != nprod) {
    stop(
      "nprod does not correspond to the number of products. Found ",
      length(NameProd), "."
    )
  }

  if (anyDuplicated(Data[, 1:2])) {
    stop("Duplicated subject-product combinations were found.")
  }


  # ------------------------------------------------------------
  # Reorder data consistently
  # Subject first, then product
  # ------------------------------------------------------------

  ord <- order(
    match(Data[[1]], NameSub),
    match(Data[[2]], NameProd)
  )

  Data <- Data[ord, , drop = FALSE]


  # ------------------------------------------------------------
  # Liking matrix: products x subjects
  # ------------------------------------------------------------

  liking <- matrix(
    NA_real_,
    nrow = nprod,
    ncol = nsub,
    dimnames = list(
      as.character(NameProd),
      as.character(NameSub)
    )
  )

  product_position <- match(Data[[2]], NameProd)
  subject_position <- match(Data[[1]], NameSub)

  liking[
    cbind(product_position, subject_position)
  ] <- Data[[liking_col]]


  # ------------------------------------------------------------
  # Extract JAR attributes
  # ------------------------------------------------------------

  JAR_cols <- setdiff(
    seq_len(ncol(Data)),
    c(
      1,
      2,
      match(liking_col, colnames(Data))
    )
  )

  if (length(JAR_cols) == 0) {
    stop("No JAR variables were found.")
  }

  JAR <- Data[, JAR_cols, drop = FALSE]

  NameAttr <- colnames(JAR)


  # ------------------------------------------------------------
  # Check JAR values
  # ------------------------------------------------------------

  valid_levels <- if (levelsJAR == 3) {
    c(1, 2, 3)
  } else {
    c(1, 2, 3, 4, 5)
  }

  observed <- unique(unlist(JAR))
  observed <- observed[!is.na(observed)]

  if (!all(observed %in% valid_levels)) {
    stop("Unexpected values were found in the JAR variables.")
  }


  # ------------------------------------------------------------
  # Binary coding:
  #
  # 1 = JAR
  # 0 = non-JAR
  # ------------------------------------------------------------

  if (levelsJAR == 3) {

    JAR_binary <- as.data.frame(
      lapply(
        JAR,
        function(x) {
          ifelse(
            is.na(x),
            NA,
            as.integer(x == 2)
          )
        }
      )
    )

  } else {

    JAR_binary <- as.data.frame(
      lapply(
        JAR,
        function(x) {
          ifelse(
            is.na(x),
            NA,
            as.integer(x == 3)
          )
        }
      )
    )
  }

  colnames(JAR_binary) <- NameAttr


  # ------------------------------------------------------------
  # Convert into CLUSCATA format
  # ------------------------------------------------------------

  CATA <- change_cata_format(
    JAR_binary,
    nprod,
    ncol(JAR_binary),
    nsub = nsub,
    1,
    NameProds = NameProd,
    NameAttr = NameAttr
  )


  # ------------------------------------------------------------
  # Combine JAR and liking
  # ------------------------------------------------------------

  Datafinal <- combinCATALiking(
    CATA,
    liking,
    scale = scale
  )


  # ------------------------------------------------------------
  # Outputs
  # ------------------------------------------------------------

  return(
    list(
      Datafinal = Datafinal,
      CATA = CATA,
      liking = liking,
      JAR_binary = JAR_binary,
      NameSub = NameSub,
      NameProd = NameProd
    )
  )
}

