qpcr <-
function(data, response, gene, control_gene, fixed1, fixed2=NULL, rep_id, block=FALSE, splitplot=FALSE, contrasts=TRUE, interaction=FALSE, adjusted=FALSE, df=NULL){
  require(lme4)
  require(multcomp)
  if (!is.null(gene)) data[,gene] <- as.factor(data[,gene])
  if (!is.null(gene)) if (!(control_gene %in% data[,gene])) stop(paste("control_gene", control_gene, "is not a level within the factor", gene))
  data[,fixed1] <- as.factor(data[,fixed1])
  if (!is.null(fixed2)) data[,fixed2] <- as.factor(data[,fixed2])
  data[,rep_id] <- as.factor(data[,rep_id])

  if (is.null(gene)){
    ###############################################################################
    ################# one gene ####################################################
   # fixed effects formula
    if (is.null(fixed2)){
      fform <- paste(response, " ~ ", fixed1, " -1", sep="")
      fform2 <- paste(response, " ~ 1", sep="")
    } else {
      fform <- paste(response, " ~ ", fixed1, ":", fixed2, " -1", sep="")
      fform2 <- paste(response, " ~ 1", sep="")
    }
  # random effects formula
    if (splitplot){
      if (is.null(fixed2)) stop("fixed2 has to be defined for a splitplot design")
      data$FR <- data[,rep_id]:data[,fixed1]
      if (block){
        rfstr <- paste(rep_id, "/", fixed1, "/", fixed2)
      } else {
        rfstr <- paste("FR/", fixed2, sep="")
      }
      rform <- paste("+ (1|",rfstr, ")", sep="")
    } else {
      if (is.null(fixed2)){
        data$ff <- data[,fixed1]
      } else {
        data$ff <- data[,fixed1]:data[,fixed2]
      }
      if (block){
        rfstr <- paste(rep_id, "/ff")
      } else {
        rfstr <- paste(rep_id, ":ff", sep="")
      }
      rform <- paste("+ (1|",rfstr, ")", sep="")
    }
    ##########################################
    ### contrasts
    if (contrasts){
      if (is.null(fixed2)){
        fl1 <- levels(data[,fixed1])
        K <- contrMat(table(data[,fixed1]), "Tukey")
      } else {
        fl1 <- levels(data[,fixed1])
        fl2 <- levels(data[,fixed2])

        k1 <- contrMat(table(data[,fixed1]), "Tukey")
        k2 <- rbind(rep(1/length(fl2), times=length(fl2)))
        K1 <- kronecker(k2, k1)
        rownames(K1) <- rownames(k1)

        k1 <- rbind(rep(1/length(fl1), times=length(fl1)))
        k2 <- contrMat(table(data[,fixed2]), "Tukey")
        K2 <- kronecker(k2, k1)
        rownames(K2) <- rownames(k2)

        k1 <- contrMat(table(data[,fixed1]), "Tukey")
        k2 <- diag(length(fl2))
        K12 <- kronecker(k2, k1)
        rownames(K12) <- paste(rownames(k1), " | ", rep(fl2, each=nrow(k1)), sep="")

        k1 <- diag(length(fl1))
        k2 <- contrMat(table(data[,fixed2]), "Tukey")
        K21 <- kronecker(k2, k1)
        rownames(K21) <- paste(rownames(k2), " | ", rep(fl1, each=nrow(k2)), sep="")

        if (interaction){
          k1 <- contrMat(table(data[,fixed1]), "Tukey")
          k2 <- contrMat(table(data[,fixed2]), "Tukey")
          KI <- kronecker(k2, k1)
          rownames(KI) <- paste(rownames(k1), " | ", rownames(k2), sep="")
          K <- rbind(K1, K2, K12, K21, KI)
        } else {
          K <- rbind(K1, K2, K12, K21)
        }
      }
    }
    
  } else {
    ########################################################################################
    ######################## multiple genes ###################################################
    # dummy variables for pdDiag definition
    X <- data.frame(model.matrix(as.formula(paste("~", gene, "-1")), data=data))
    data <- cbind(data, X)

    # fixed effects formula
    if (is.null(fixed2)){
      fform <- paste(response, " ~ ", fixed1, ":", gene, " -1", sep="")
      fform2 <- paste(response, " ~ ", fixed1, "+", gene, sep="")
    } else {
      fform <- paste(response, " ~ ", fixed1, ":", fixed2, ":", gene, " -1", sep="")
      fform2 <- paste(response, " ~ ", fixed1, "+", fixed2, "+", gene, sep="")
    }

  # random effects formula
    if (splitplot){
      if (is.null(fixed2)) stop("fixed2 has to be defined for a splitplot design")
      data$FR <- data[,rep_id]:data[,fixed1]
      if (block){
        rfstr <- paste(rep_id, "/", fixed1, "/", fixed2)
      } else {
        rfstr <- paste("FR/", fixed2, sep="")
      }
      rform <- paste("+ (1|",rfstr, ")", paste(sapply(1:ncol(X), function(i) paste(" + (0+", names(X)[i], "|", rfstr, ")", sep="")), collapse=""), sep="")
    } else {
      if (is.null(fixed2)){
        data$ff <- data[,fixed1]
      } else {
        data$ff <- data[,fixed1]:data[,fixed2]
      }
      if (block){
        rfstr <- paste(rep_id, "/ff")
      } else {
        rfstr <- paste(rep_id, ":ff", sep="")
      }
      rform <- paste("+ (1|",rfstr, ")", paste(sapply(1:ncol(X), function(i) paste(" + (0+", names(X)[i], "|", rfstr, ")", sep="")), collapse=""), sep="")
    }

    ##########################################
    ### contrasts
    if (contrasts){
      if (is.null(fixed2)){
        fl1 <- levels(data[,fixed1])
        glev <- levels(data[,gene])
        wcg <- which(glev == control_gene)
        gk <- -1*contrMat(table(data[,gene]), "Dunnett", base=wcg)

        k1 <- contrMat(table(data[,fixed1]), "Tukey")
        K <- kronecker(gk, k1)
        rownames(K) <- paste(rownames(k1), " | ", rep(rownames(gk), each=nrow(k1)), sep="")
      } else {
        fl1 <- levels(data[,fixed1])
        fl2 <- levels(data[,fixed2])
        glev <- levels(data[,gene])
        wcg <- which(glev == control_gene)
        gk <- -1*contrMat(table(data[,gene]), "Dunnett", base=wcg)

        k1 <- contrMat(table(data[,fixed1]), "Tukey")
        k2 <- rbind(rep(1/length(fl2), times=length(fl2)))
        K1 <- kronecker(gk, kronecker(k2, k1))
        rownames(K1) <- paste(rownames(k1), " | ", rep(rownames(gk), each=nrow(k1)), sep="")

        k1 <- rbind(rep(1/length(fl1), times=length(fl1)))
        k2 <- contrMat(table(data[,fixed2]), "Tukey")
        K2 <- kronecker(gk, kronecker(k2, k1))
        rownames(K2) <- paste(rownames(k2), " | ", rep(rownames(gk), each=nrow(k2)), sep="")

        k1 <- contrMat(table(data[,fixed1]), "Tukey")
        k2 <- diag(length(fl2))
        K12 <- kronecker(gk, kronecker(k2, k1))
        rownames(K12) <- paste(rownames(k1), " | ", rep(fl2, each=nrow(k1)), " | ", rep(rownames(gk), each=nrow(k1)*nrow(k2)), sep="")

        k1 <- diag(length(fl1))
        k2 <- contrMat(table(data[,fixed2]), "Tukey")
        K21 <- kronecker(gk, kronecker(k2, k1))
        rownames(K21) <- paste(rep(rownames(k2), each=nrow(k1)), " | ", rep(fl1, times=nrow(k2)), " | ", rep(rownames(gk), each=nrow(k1)*nrow(k2)), sep="")

        if (interaction){
          k1 <- contrMat(table(data[,fixed1]), "Tukey")
          k2 <- contrMat(table(data[,fixed2]), "Tukey")
          KI <- kronecker(gk, kronecker(k2, k1))
          rownames(KI) <- paste(rownames(k1), " | ", rownames(k2), " | ", rep(rownames(gk), each=nrow(k1)*nrow(k2)), sep="")
          K <- -1*rbind(K1, K2, K12, K21, KI)
        } else {
          K <- -1*rbind(K1, K2, K12, K21)
        }
      }
    }
  }

  form <- as.formula(paste(fform, rform, sep=""))
  fm <- lmer(form, data=data)
  if (contrasts){    
    gg <- glht(fm, K)
    if (!is.null(df)){
      if (!is.integer(df)) warning(paste("df is not an integer! A df of",floor(as.integer(df)) ,"is used instead", sep=" "))
      gg$df <- as.integer(floor(df))   
    }
    pvalues <- if (adjusted) summary(gg) else summary(gg, test=univariate())
    degfree <- gg$df
  } else {
    pvalues <- NULL
    K <- NULL
    degfree <- NULL
  }

  out <- list(model=fm, pvalues=pvalues, contrasts=K, df=degfree)
  class(out) <- "qpcrmacro"
  return(out)
}
