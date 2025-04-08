#' Function to calculate phylogenetic gamma, alpha, beta diversity and dissimilarity measures
#'
#' \code{iNEXTseq}: This function calculates interpolated and extrapolated phylogenetic gamma, alpha, beta diversity, and dissimilarity measures based on abundance or incidence data, using the framework developed by Chiu et al.
#'
#' @param data OTU data can be input as a \code{matrix/data.frame} (species by assemblages), or a \code{list} of \code{matrices/data.frames}, each matrix represents a species-by-assemblages abundance or incidence matrix.\cr
#' @param q a numerical vector specifying the diversity orders. Default is \code{c(0, 1, 2)}.
#' @param base sample-size-based rarefaction and extrapolation for gamma and alpha diversity (\code{base = "size"}) or coverage-based rarefaction and extrapolation for gamma, alpha, and beta diversity (\code{base = "coverage"}). Default is \code{"coverage"}.
#' @param datatype type of input data: \code{"abundance"} (default) or \code{"incidence_raw"}. If \code{"incidence_raw"}, the first column of each matrix must be the total number of sampling units in each assemblage.
#' @param level A numerical vector specifying the particular value of sample coverage (between 0 and 1 when \code{base = "coverage"}) or sample size (\code{base = "size"}). See Details.
#' @param nboot a positive integer specifying the number of bootstrap replications when assessing sampling uncertainty and constructing confidence intervals. Use \code{0} to skip the bootstrap procedure. Default is \code{10}.
#' @param conf a number < 1 specifying the level of confidence interval. Default is \code{0.95}.
#' @param PDtree a phylogenetic tree in Newick format for all observed species in the pooled assemblage.
#' @param PDtype Type of PD to report: \code{"meanPD"} (default)  or \code{"PD"}
#' @param PDreftime a numerical value specifying the reference time for PD. Default is \code{NULL} (i.e., the age of the root of PDtree).
#'
#' @import tidyverse
#' @import magrittr
#' @import ggplot2
#' @import abind
#' @import phytools
#' @import phyclust
#' @import future.apply
#' @import ade4
#' @import tibble
#' @import stringr
#' @import forcats
#' @import dplyr
#' @import RColorBrewer
#' @import iNEXT.3D
#' @import iNEXT.beta3D
#'
#' @importFrom tidyr gather
#' @importFrom phyclust get.rooted.tree.height
#' @importFrom stats rmultinom rbinom qnorm sd optimize
#' @importFrom grDevices hcl
#' @importFrom utils combn
#' @importFrom stats quantile
#' @importFrom utils tail
#' @return If \code{base = "coverage"}, returns a list of seven data frames with three coverage-based diversity (gamma, alpha, and beta) and four dissimilarity measures. If \code{base = "size"}, returns a list of two data frames with gamma and alpha diversity.
#'
#' @examples
#' # Abundance data example
#' data("fungi")
#' data("fungi_tree")
#' output <- iNEXTseq(fungi[1], q = c(0,1,2), level = seq(0.5, 1, 0.05), nboot = 10,
#'                    conf = 0.95, PDtree = fungi_tree, PDreftime = NULL, PDtype = 'meanPD')
#'
#' # Incidence data example (assuming incidence_data is formatted correctly)
#' # output <- iNEXTseq(incidence_data, q = c(0,1,2), datatype = "incidence_raw",
#' #                    level = seq(0.5, 1, 0.05), PDtree = incidence_tree, PDtype = 'meanPD')
#'
#' @export
iNEXTseq <- function(data, q=c(0,1,2), base = "coverage", datatype = "abundance",level = NULL, nboot = 10, conf = 0.95, PDtree = NULL, PDreftime = NULL, PDtype = 'meanPD'){


  out = iNEXTbeta3D(data, diversity = "PD", q = q, datatype = datatype, base = base,
                    level = level, nboot = nboot, conf = conf, PDtype = PDtype,
                    PDtree = PDtree, PDreftime = PDreftime)
  out
  # UniFrac_out = list()
  # if(class(data)=="list"){
  #   for(i in 1:length(data)){
  #     UniFrac_out[[i]] = list(C = out[[i]]$C, U = out[[i]]$U)
  #   }
  #   if(is.NULL(names(data))){
  #     names(UniFrac_out) = paste0("Region_", 1:length(data))
  #   }else{
  #     names(UniFrac_out) = names(data)
  #   }
  #
  # }else{
  #   UniFrac_out[[1]] = list(C = out[[1]]$C, U = out[[1]]$U)
  #   names(UniFrac_out) = "Region_1"
  # }
  #UniFrac_out
}

#' ggplot2 extension for an iNEXT.seq object
#'
#' \code{ggiNEXTseq}: the \code{ggplot} extension for \code{\link{iNEXTseq}} object to plot coverage- or sample-sized-based rarefaction/extrapolation curves for phylogenetic diversity decomposition and dissimilarity measure.
#'
#' @param output the output from iNEXTseq.
#' @param type (required only when \code{base = "coverage"}), selection of plot type : \cr
#' \code{type = "B"} for plotting the gamma, alpha, and beta diversity;  \cr
#' \code{type = "D"} for plotting 4 turnover dissimilarities.
#'
#' @return a figure for phylogenetic diversity decomposition or dissimilarity measure.
#'
#' @examples
#'
#' data("fungi")
#' data("fungi_tree")
#' output <- iNEXTseq(fungi[1], q = c(0,1,2), nboot = 10, PDtree = fungi_tree)
#' ggiNEXTseq(output, type = "B")
#'
#' @export
ggiNEXTseq <- function(output, type = "B"){

  ggiNEXTbeta3D(output, type = type)

  # cbPalette <- rev(c("#999999", "#E69F00", "#56B4E9", "#009E73", "#330066", "#CC79A7", "#0072B2", "#D55E00"))
  # ylab = "UniFrac distance"
  #
  #   C = lapply(output, function(y) y[["C"]]) %>% do.call(rbind,.) %>% mutate(div_type = "1-CqN") %>% as_tibble()
  #   U = lapply(output, function(y) y[["U"]]) %>% do.call(rbind,.) %>% mutate(div_type = "1-UqN") %>% as_tibble()
  #   C = C %>% filter(Method != 'Observed')
  #   U = U %>% filter(Method != 'Observed')
  #   C[C == 'Observed_alpha'] = U[U == 'Observed_alpha'] = 'Observed'
  #
  #   df = rbind(C, U)
  #   for (i in unique(C$Order.q)) df$Order.q[df$Order.q == i] = paste0('q = ', i)
  #   df$div_type <- factor(df$div_type, levels = c("1-CqN", "1-UqN"))
  #
  #   id_obs = which(df$Method == 'Observed')
  #
  #   for (i in 1:length(id_obs)) {
  #     new = df[id_obs[i],]
  #     new$SC = new$SC - 0.0001
  #     new$Method = 'Rarefaction'
  #
  #     newe = df[id_obs[i],]
  #     newe$SC = newe$SC + 0.0001
  #     newe$Method = 'Extrapolation'
  #
  #     df = rbind(df, new, newe)
  #   }
  #
  #   lty = c(Rarefaction = "solid", Extrapolation = "dashed")
  #   df$Method = factor(df$Method, levels = c('Rarefaction', 'Extrapolation', 'Observed'))
  #
  #   double_size = unique(df[df$Method == "Observed",]$Size)*2
  #   double_extrapolation = df %>% filter(Method == "Extrapolation" & round(Size) %in% double_size)
  #
  #   ggplot(data = df, aes(x = SC, y = Estimate, col = Region)) +
  #     geom_ribbon(aes(ymin = LCL, ymax = UCL, fill = Region, col = NULL), alpha = transp) +
  #     geom_line(data = subset(df, Method != 'Observed'), aes(linetype = Method), size=1.1) + scale_linetype_manual(values = lty) +
  #     # geom_line(lty=2) +
  #     geom_point(data = subset(df, Method == 'Observed' & div_type == "Gamma"), shape = 19, size = 2) +
  #     geom_point(data = subset(df, Method == 'Observed' & div_type != "Gamma"), shape = 1, size = 2, stroke = 1.2)+
  #     geom_point(data = subset(double_extrapolation, div_type == "Gamma"), shape = 17, size = 2) +
  #     geom_point(data = subset(double_extrapolation, div_type != "Gamma"), shape = 2, size = 2, stroke = 1.2) +
  #     scale_colour_manual(values = cbPalette) +
  #     scale_fill_manual(values = cbPalette) +
  #     facet_grid(div_type ~ Order.q, scales = scale) +
  #     theme_bw() +
  #     theme(legend.position = "bottom", legend.title = element_blank()) +
  #     labs(x = 'Sample coverage', y = ylab)
}


#' function to calculate observed or asymptotic phylogenetic gamma, alpha, beta diversity and dissimilarity of order q
#'
#' \code{ObsAsyPD} computes observed and asymptotic diversity of order q between 0 and 2 (in increments of 0.2) for phylogenetic gamma, alpha, beta diversity and dissimilarity; these values with different order q can be used to depict a q-profile in the \code{ggObsAsyPD} function.
#'
#' @param data OTU data can be input as a \code{matrix/data.frame} (species by assemblages), or a \code{list} of \code{matrices/data.frames}, each matrix represents species-by-assemblages abundance matrix.\cr
#' @param q a numerical vector specifying the diversity orders. Default is \code{seq(0, 2, 0.2)}.
#' @param weight (required only when \code{type = "mle"} and \code{decomposition = "relative"}) weight for relative decomposition empirical estimate. Select size-weighted \code{("size")}, equal-weighted \code{("equal")} or a numerical vector for weight. Default is \code{"size"}.
#' @param nboot a positive integer specifying the number of bootstrap replications when assessing sampling uncertainty and constructing confidence intervals. Bootstrap replications are generally time consuming. Enter \code{0} to skip the bootstrap procedures. Default is \code{10}.
#' @param conf a positive number < 1 specifying the level of confidence interval. Default is \code{0.95}.
#' @param PDtree a phylogenetic tree in Newick format for all observed species in the pooled assemblage.
#' @param type estimate type: empirical \code{(type = "mle")} or asymptotic estimate \code{(type = "est")}. Default is \code{"mle"}.
#' @param decomposition decomposition type: relative \code{(decomposition = "relative")} or absolute decomposition \code{(decomposition = "absolute")}. Default is \code{"relative"}.
#'
#' @return a data frames with observed or asymptotic phylogenetic diversity (gamma, alpha, and beta) and four types dissimilarity measure for each dataset.
#'
#' @examples
#'
#' data("fungi")
#' data("fungi_tree")
#' ObsAsyPD_output <- ObsAsyPD(fungi[1], q = seq(0, 2, 0.2), PDtree = fungi_tree)
#'
#' @export
ObsAsyPD <- function(data, q = seq(0, 2, 0.2), weight = "size", nboot = 10, conf = 0.95,
                     PDtree, type = "mle", decomposition = "relative") {

  if (inherits(data, "data.frame") | inherits(data, "matrix"))
    data = list(Dataset_1 = data)
  if (inherits(data, "list")) {
    if (is.null(names(data)))
      dataset_names = paste0("Dataset_", 1:length(data))
    else dataset_names = names(data)
    Ns = sapply(data, ncol)
    data_list = data
  }
  if (type == "mle" & inherits(weight, "numeric")){
    if (length(data_list) != 1 & length(unique(Ns)) != 1)
      stop("Please select datasets with same number of assemblages N.")
  }


  for_each_dataset = function(data, dataset_name){

    dat <- data[rowSums(data)>0, ]
    if (decomposition == "relative"){
      method = phy.rel
    }else if (decomposition == "absolute"){
      method = phy.abs
    }

    rtip <- PDtree$tip.label[!PDtree$tip.label %in% rownames(dat)]
    rtree <- drop.tip(PDtree, rtip)
    tmp <- TranMul(dat, rtree)
    rtreephy <- newick2phylog(convToNewick(rtree))

    if (type == "mle"){
      if (inherits(weight, "numeric")){
        wk <- weight/sum(weight)
      } else if (weight == "size"){
        wk <- colSums(dat)/sum(dat) #size weight
      } else if (weight == "equal"){
        wk <- rep(1/ncol(dat), ncol(dat)) #equal weight
      }
    }else if (type == "est"){
      wk <- colSums(dat)/sum(dat) #size weight
    }else {
      wk <- colSums(dat)/sum(dat) #size weight
    }

    est <- sapply(q, function(i) method(dat, tmp, i, rtreephy, wk, type))

    if(nboot != 0){
      boot.est <- bootstrap.q.Beta_OBS(data = dat, rtree = rtree, tmp = tmp, q = q, nboot = nboot, wk = wk, type, method)

      ##
      # test <- boot.est[seq_len(2), , ]
      # #is.infinite(sum(boot.est))
      # #test <- boot.est[head(seq_len(dim(boot.est)[1]),H),1,]
      # id <- apply(test, 1:2, function(x) {
      #   bb <- x
      #   q1 <- quantile(bb, 0.25)
      #   q3 <- quantile(bb, 0.75)
      #   q1 <- q1-1.5*(q3-q1)
      #   q3 <- q1+1.5*(q3-q1)
      #   which(bb >= q1 & bb <= q3)
      # })
      # index <- Reduce(function(x,y) {intersect(x,y)}, id)
      # boot.est <- boot.est[,,index]
      ##

      #boot.est[tail(seq_len(dim(boot.est)[1]),-2*H), ][boot.est[tail(seq_len(dim(boot.est)[1]),-2*H), ]<0] <- 0
      #boot.est[tail(seq_len(dim(boot.est)[1]),-2*H), ][boot.est[tail(seq_len(dim(boot.est)[1]),-2*H), ]>1] <- 1
      #dim(boot.est)
      #diff.boot.est <- as.data.frame(apply(boot.est,3, tail, n = H))
      out = lapply(seq_along(q), function(i){
        x = transconf(Bresult = boot.est[,i,], est = est[,i], conf)
        colnames(x) = c("Estimator", "Bootstrap S.E.", "LCL", "UCL")
        rownames(x) = rownames(est)
        return(x)}) %>% do.call(rbind,.)

      if (type == "spader"){
        Order.q = rep(q, each = 2)
        Method = rep(rownames(out)[1:2], length(q))
        rownames(out) = NULL
      }
      else {
        Order.q = rep(q, each = 7)
        Method = rep(rownames(out)[1:7], length(q))
        rownames(out) = NULL
      }
      out = cbind(Dataset = dataset_name, Method, Order.q, as.data.frame(out), Decomposition = decomposition)

      out
      #
      # CL = t(sapply(seq_len(nrow(boot.est)), function(j) transconf(matrix(boot.est[j,,], nrow=length(q)), est[j,], conf)))
      # rownames(CL) <- rownames(est)
      # colnames(CL) <- c(sapply(paste0("q=",q), function(k) paste(k,c("est", "bt.sd", "LCL", "UCL"), sep = "_")))
      #
    }else{
      out <- lapply(seq_along(q), function(i){
        x = data.frame("Estimator" = est[,i], Bootstraps.e = NA, LB = NA, UB = NA)
        colnames(x) = c("Estimator", "Bootstrap S.E.", "LCL", "UCL")
        return(x)}) %>% do.call(rbind,.)

      if (type == "spader"){
        Order.q = rep(q, each = 2)
        Method = rep(rownames(out)[1:2], length(q))
        rownames(out) = NULL
      }
      else {
        Order.q = rep(q, each = 7)
        Method = rep(rownames(out)[1:7], length(q))
        rownames(out) = NULL
      }

      out = cbind(Dataset = dataset_name, Method, Order.q, as.data.frame(out), Decomposition = decomposition)
      out$'Bootstrap S.E.' = as.numeric(out$'Bootstrap S.E.')
      out$LCL = as.numeric(out$LCL)
      out$UCL = as.numeric(out$UCL)

      out
    }
  }

  output = lapply(1:length(data_list), function(i)
    for_each_dataset(data = data_list[[i]], dataset_name = dataset_names[i]))

  names(output) = dataset_names
  class(output) <- c("ObsAsyPD")
  return(output)
}


#' ggplot2 extension for an ObsAsyPD object
#'
#' \code{ggObsAsyPD}: the \code{ggplot} extension for \code{\link{ObsAsyPD}} object to plot order q against phylogenetic diversity decomposition or dissimilarity measure.
#'
#' @param output The output from \code{\link{ObsAsyPD}}.
#' @param type Selection of plot type: \cr
#'   \code{"B"} for plotting gamma, alpha, and beta diversity; \cr
#'   \code{"D"} for plotting four turnover dissimilarities.
#'
#' @return A ggplot object showing phylogenetic diversity or dissimilarity.
#'
#' @examples
#' data("fungi")
#' data("fungi_tree")
#' ObsAsyPD_output <- ObsAsyPD(fungi[1], q = seq(0, 2, 0.2), PDtree = fungi_tree)
#' ggObsAsyPD(ObsAsyPD_output, type = "B")
#'
#' @export
ggObsAsyPD <- function(output, type = "B"){

  cbPalette <- rev(c("#999999", "#E69F00", "#56B4E9", "#009E73", "#330066",
                     "#CC79A7", "#0072B2", "#D55E00"))

  plot_output <- do.call(rbind, output)
  if (type == "B") {
    plot_output = plot_output[c(grep(c("Gamma"), plot_output$Method),
                                grep(c("Alpha"), plot_output$Method),
                                grep(c("Beta"), plot_output$Method)) %>% sort(), ]
  }
  else if (type == "D") {
    plot_output = plot_output[grep("1-", plot_output$Method), ]
  }

  out = ggplot(plot_output, aes(x = Order.q, y = Estimator, colour = fct_inorder(Dataset), fill = fct_inorder(Dataset))) +
    geom_line(size = 1.5) +
    geom_ribbon(aes(ymin = LCL, ymax = UCL), linetype = 0, alpha = 0.2) +
    facet_grid(fct_inorder(Method) ~ ., scales = "free") +
    scale_colour_manual(values = cbPalette) +
    scale_fill_manual(values = cbPalette) +
    theme_bw() +
    theme(legend.position = "bottom", legend.box = "vertical",
          legend.key.width = unit(1.2, "cm"), legend.title = element_blank(),
          legend.margin = margin(0, 0, 0, 0), legend.box.margin = margin(-10, -10, -5, -10),
          text = element_text(size = 16), plot.margin = unit(c(5.5, 5.5, 5.5, 5.5), "pt"),
          strip.text = element_text(size = 15, face = "bold"),
          axis.title = element_text(hjust = 0.5, size = 15, face = "bold"),
          axis.text.x = element_text(size = 12),
          axis.text.y = element_text(size = 12),
          legend.text = element_text(size = 13)) +
    guides(linetype = guide_legend(keywidth = 2.5))

  if (type == "B") {
    out = out + ylab("Phylogenetic diversity") + xlab("Order q")
  }
  else if (type == "D") {
    out = out + ylab("Phylogenetic dissimilarity") + xlab("Order q")
  }

  return(out)
}



#' Compute phylogenetic diversity and dissimilarity using Routledge's method
#'
#' This function calculates phylogenetic gamma, alpha, and beta diversity and associated dissimilarity indices
#' using Routledge's weighted-mean approach under either coverage- or size-based rarefaction/extrapolation.
#'
#' @param data OTU abundance data: a matrix/data.frame (species by assemblages), or a list of such objects.
#' @param q A numeric vector specifying diversity orders. Default is \code{c(0, 1, 2)}.
#' @param base The type of rarefaction/extrapolation: \code{"coverage"} (default) or \code{"size"}.
#' @param level A vector specifying target sample coverage (when \code{base = "coverage"}) or sample size (when \code{base = "size"}). If \code{NULL}, the function determines appropriate levels automatically.
#' @param nboot Number of bootstrap replications. Default is \code{10}; set to \code{0} to skip uncertainty estimation.
#' @param conf Confidence level for intervals. Default is \code{0.95}.
#' @param PDtree A phylogenetic tree of class \code{phylo} with tips matching row names of \code{data}.
#' @param PDreftime Reference time for PD. If \code{NULL}, tree height is used.
#' @param PDtype Type of PD to report: \code{"meanPD"} (default).
#'
#' @return A list with elements \code{gamma}, \code{alpha}, \code{beta}, and dissimilarity metrics \code{1-C^*}, \code{1-U^*}, \code{1-V^*}, \code{1-S^*} (for coverage-based); or \code{gamma}, \code{alpha}, and \code{joint} (for size-based).
#'
#' @import tidyverse
#' @import magrittr
#' @import ggplot2
#' @import abind
#' @import phytools
#' @import phyclust
#' @import future.apply
#' @import ade4
#' @import tibble
#' @import stringr
#' @import forcats
#' @import RColorBrewer
#' @import iNEXT.3D
#' @import iNEXT.beta3D
#' @importFrom ape node.depth.edgelength
#' @importFrom tidytree drop.tip keep.tip
#' @importFrom dplyr where
#' @importFrom magrittr extract
#' @importFrom tidyr gather
#' @importFrom phyclust get.rooted.tree.height
#' @importFrom stats rmultinom rbinom qnorm sd optimize
#' @importFrom grDevices hcl
#' @importFrom utils combn
#' @importFrom utils getFromNamespace
#' @examples
#' # Abundance data example
#' data("fungi")
#' data("fungi_tree")
#' output <- iNEXT_seq_Relative(fungi[1], q = c(0,1,2), level = seq(0.5, 1, 0.05), nboot = 10,
#'                    conf = 0.95, PDtree = fungi_tree, PDreftime = NULL, PDtype = 'meanPD')
#'
#' @export
iNEXT_seq_Relative<-function(data, q = c(0, 1, 2),  base = 'coverage', level = NULL, nboot = 10, conf = 0.95,
                              PDtree = NULL, PDreftime = NULL, PDtype = 'meanPD'){
  
  PhD.m.est <- .internal_fetch("PhD.m.est", "iNEXT.3D")
  phyBranchAL_Abu <- .internal_fetch("phyBranchAL_Abu", "iNEXT.3D")
  bootstrap_population_multiple_assemblage <- .internal_fetch("bootstrap_population_multiple_assemblage", "iNEXT.beta3D")
  
  if (inherits(data, "data.frame") | inherits(data, "matrix"))
    data = list(Dataset_1 = data)
  if (inherits(data, "list")) {
    if (is.null(names(data)))
      dataset_names = paste0("Dataset_", 1:length(data))
    else dataset_names = names(data)
    Ns = sapply(data, ncol)
    data_list = data
  }

  pool.name <- lapply(data_list, function(x) rownames(x)) %>% unlist %>% unique

  if (sum(c(duplicated(PDtree$tip.label), duplicated(PDtree$node.label[PDtree$node.label!=""])))>0)
    stop("The phylo tree should not contains duplicated tip or node labels, please remove them.", call. = FALSE)

  if ( is.null(pool.name) )
    stop("Row names of data must be the species names that match tip names in tree and thus can not be empty.", call. = FALSE)

  if (sum(pool.name %in% PDtree$tip.label) != length(pool.name))
    stop("Data and tree tip label contain unmatched species", call. = FALSE)

  if (is.null(conf)) conf = 0.95
  tmp = qnorm(1 - (1 - conf)/2)
  trunc = ifelse(is.null(level), T, F)

  if (base == 'coverage' ) {

    if ( is.null(level) ) {

      level <- lapply(1:length(data_list), function(i) {


        level = seq(0.5, 1, 0.025)


        n = sum(data_list[[i]])
        data_gamma = rowSums(data_list[[i]])
        data_gamma = data_gamma[data_gamma>0]
        data_alpha = as.matrix(data_list[[i]]) %>% as.vector

        ref_gamma = Coverage(data_gamma, 'abundance', n)
        ref_alpha = Coverage(data_alpha, 'abundance', n)
        ref_alpha_max = Coverage(data_alpha, 'abundance', n*2)
        ref_gamma_max = Coverage(data_gamma, 'abundance', n*2)


        c(level, ref_gamma, ref_alpha, ref_alpha_max, ref_gamma_max) %>% sort %>% unique



      })

    } else {

      if (inherits(level, "numeric") | inherits(level, "integer") | inherits(level, "double")) {
        level <- list(level = level)
      }

      if ( (length(level) != length(data_list))) level <- lapply(1:length(data_list), function(x) level[[1]])


    }

  }
  else if ( base == 'size' ) {

    if ( is.null(level) ) {

      endpoint <- sapply(data_list, function(x) 2*sum(x))


      level <- lapply(1:length(data_list), function(i) {



        ni <- sum(data_list[[i]])


        mi <- floor(c(seq(1, ni-1, length.out = 20), ni, seq(ni+1, endpoint[i], length.out = 20)))
      })



    } else {

      if (inherits(level, "numeric") | inherits(level, "integer") | inherits(level, "double")) {
        level <- list(level = level)
      }

      if ( (length(level) != length(data_list))) level <- lapply(1:length(data_list), function(x) level[[1]])



    }
  }


  if (length(data_list) > 1) {

    pool.data = data_list[[1]] %>% data.frame %>% rownames_to_column()
    for (i in 2:length(data_list))
      pool.data = full_join(pool.data, data_list[[i]] %>% data.frame %>% rownames_to_column(), 'rowname')
    pool.data[is.na(pool.data)] = 0
    pool.data = pool.data %>% column_to_rownames() %>% rowSums

  } else pool.data = do.call(cbind, data_list) %>% rowSums



  pool.name = names(pool.data[pool.data>0])
  tip = PDtree$tip.label[-match(pool.name, PDtree$tip.label)]
  mytree = ape::drop.tip(PDtree, tip)

  # H_max = get.rooted.tree.height(mytree)
  H_max = max(ape::node.depth.edgelength(mytree))
  if(is.null(PDreftime)) { reft = H_max
  } else if (PDreftime <= 0) { stop("Reference time must be greater than 0. Use NULL to set it to pooled tree height.", call. = FALSE)
  } else { reft = PDreftime }


  for_each_dataset = function(data, dataset_name, N, level) {

    #data

    n = sum(data)
    Routledge_x <- data[(rowSums(data) != 0), ]

    data_gamma = rowSums(data)
    data_gamma = data_gamma[data_gamma>0]
    data_alpha = as.matrix(data) %>% as.vector
    wk <- colSums(Routledge_x) / sum(Routledge_x)

    ref_gamma = Coverage(data_gamma, 'abundance', n)
    ref_alpha = Coverage(data_alpha, 'abundance', n)
    ref_alpha_max = Coverage(data_alpha, 'abundance', n*2)
    ref_gamma_max = Coverage(data_gamma, 'abundance', n*2)

    # level = c(level, ref_gamma, ref_alpha, ref_alpha_max, ref_gamma_max) %>% sort %>% unique
    # level = level[level<1]

    m_gamma = sapply(level, function(i) coverage_to_size(data_gamma, i, datatype='abundance'))
    m_alpha = sapply(level, function(i) coverage_to_size(data_alpha, i, datatype='abundance'))

    m_Rou_alpha <- lapply(1:ncol(Routledge_x), function(j) {
      sapply(level, function(i) coverage_to_size(Routledge_x[, j], i, datatype = 'abundance'))
    })


    aL = phyBranchAL_Abu(phylo = PDtree, data = data_gamma, rootExtend = T, refT = reft)
    aL$treeNabu$branch.length = aL$BLbyT[,1]
    aL_table_gamma = aL$treeNabu %>% select(branch.abun, branch.length, tgroup)
    #
    aL_table_k = lapply(1:ncol(Routledge_x), function(k){

      x_k = Routledge_x[,k]
      names(x_k) = rownames(Routledge_x)

      aL = phyBranchAL_Abu(phylo = PDtree, data = x_k, rootExtend = T, refT = reft)
      aL$treeNabu$branch.length = aL$BLbyT[,1]
      aL_table = aL$treeNabu %>% select(branch.abun, branch.length, tgroup)
    })
    
    
    gamma = PhD.m.est(ai = aL_table_gamma$branch.abun, Lis = as.matrix(aL_table_gamma$branch.length), m = m_gamma, q = q, nt = n, reft = reft, cal = "PD") %>% t %>% as.data.frame %>%
      set_colnames(q) %>% gather(Order.q, Estimate) %>%
      mutate(SC = rep(level, length(q)), Size = rep(m_gamma, length(q)))


    aL_table_alpha = c()

    for (i in 1:N){

      x = data[data[,i]>0,i]
      names(x) = rownames(data)[data[,i]>0]

      aL = phyBranchAL_Abu(phylo = PDtree, data = x, rootExtend = T, refT = reft)
      aL$treeNabu$branch.length = aL$BLbyT[,1]
      aL_table = aL$treeNabu %>% select(branch.abun, branch.length, tgroup)

      aL_table_alpha = rbind(aL_table_alpha, aL_table)

    }

    tmp_alphaR_list <- lapply(1:ncol(Routledge_x), function(k) {
      ai_k <- aL_table_k[[k]]$branch.abun
      Lis_k <- as.matrix(aL_table_k[[k]]$branch.length)
      nt_k <- sum(Routledge_x[,k])  # or colSums(Routledge_x)[k]


      sapply(1:length(level), function(i) {
        m_i <- m_Rou_alpha[[k]][i]

        PhD.m.est(ai = ai_k,
                             Lis = Lis_k,
                             m = m_i,
                             q = q,
                             nt = nt_k,
                             reft = reft,
                             cal = "PD")
      }, simplify = FALSE)
    })
    ###
    alphaR_df <- do.call(rbind, lapply(1:length(tmp_alphaR_list), function(k) {
      do.call(rbind, lapply(1:length(tmp_alphaR_list[[k]]), function(i){
        tmp_1 <- tmp_alphaR_list[[k]][[i]] %>% t %>% as.data.frame()
        colnames(tmp_1) <- q
        tmp_1 %>%
          gather(Order.q, Estimate) %>%
          mutate(
            Sample = k,
            SC = level[i],                     # coverage level
            Size = m_Rou_alpha[[k]][i]        # sample size
          )
      }))
    }))
    qPDm = PhD.m.est(ai = aL_table_alpha$branch.abun, Lis = as.matrix(aL_table_alpha$branch.length), m = m_alpha, q = q, nt = n, reft = reft, cal = "PD")
    # qPDm = qPDm/N #Routledge joint

    # alpha = qPDm %>% t %>% as.data.frame %>%
    #   set_colnames(q) %>% gather(Order.q, Estimate) %>%
    #   mutate(SC = rep(level, length(q)), Size = rep(m_alpha, length(q)))



    alpha_results <- alphaR_df %>% rename(qPD = Estimate)


    alpha_summary <- alpha_results %>%
      group_by(Order.q, SC) %>%
      summarise(
        Estimate = if (unique(Order.q) == "1") {
          exp(sum(log(qPD) * wk))
        } else {
          (sum((qPD^(1 - as.numeric(unique(Order.q)))) * wk))^(1 / (1 - as.numeric(unique(Order.q))))
        },
        .groups = "drop"
      )



    level_len <- length(level)


    m_matrix <- do.call(cbind, m_Rou_alpha)  # 每一欄是一個 Assemblage 的 m 對應

    m_matrix_repeated <- m_matrix[rep(1:nrow(m_matrix), length(q)), ]

    colnames(m_matrix_repeated) <- paste0("Assemblage", 1:ncol(m_matrix_repeated))

    alpha <- alpha_summary %>%
      arrange(Order.q, SC) %>%
      mutate(Size = rep(m_alpha, length(q))) %>%
      bind_cols(as.data.frame(m_matrix_repeated)) %>%
      as.data.frame()


    PD_joint = qPDm %>% t %>% as.data.frame %>%
      set_colnames(q) %>% gather(Order.q, Estimate) %>%
      mutate(SC = rep(level, length(q)), Size = rep(m_alpha, length(q)))


    assem_colnames <- grep("^Assemblage", colnames(alpha), value = TRUE)


    gamma = (gamma %>%
               mutate(Method = ifelse(SC >= ref_gamma, ifelse(SC == ref_gamma, 'Observed', 'Extrapolation'), 'Rarefaction')))[,c(2,1,5,3,4)] %>%
      bind_cols(as.data.frame(m_matrix_repeated))%>%
      select(Estimate, Order.q, Method, SC, Size, all_of(assem_colnames))
    colnames(gamma) <- c('Estimate', 'Order.q', 'Method', 'SC', 'Size', assem_colnames)

    gamma$Order.q = as.numeric(gamma$Order.q)



    assem_colnames <- grep("^Assemblage", colnames(alpha), value = TRUE)


    alpha <- alpha %>%
      mutate(Method = ifelse(SC >= ref_alpha,
                             ifelse(SC == ref_alpha, 'Observed', 'Extrapolation'),
                             'Rarefaction')) %>%
      select(Estimate, Order.q, Method, SC, Size, all_of(assem_colnames))


    colnames(alpha) <- c('Estimate', 'Order.q', 'Method', 'SC', 'Size', assem_colnames)

    alpha$Order.q = as.numeric(alpha$Order.q)



    PD_joint = (PD_joint %>%
                  mutate(Method = ifelse(SC >= ref_alpha, ifelse(SC == ref_alpha, 'Observed', 'Extrapolation'), 'Rarefaction')))[,c(2,1,5,3,4)] %>%
      set_colnames(c('Estimate', 'Order.q', 'Method', 'SC', 'Size'))

    PD_joint$Order.q = as.numeric(PD_joint$Order.q)

    if (PDtype == 'meanPD') {
      gamma$Estimate = gamma$Estimate/reft
      alpha$Estimate = alpha$Estimate/reft
      PD_joint$Estimate = PD_joint$Estimate/reft
    }

    beta = alpha
    beta$Estimate = gamma$Estimate/alpha$Estimate
    PD_beta_Rmax = PD_joint
    PD_beta_Rmax$Estimate <- PD_joint$Estimate / alpha$Estimate

    C = beta %>% mutate(Estimate = ifelse(Order.q == 1, log(Estimate) / log(PD_beta_Rmax$Estimate), (Estimate^(1 - Order.q) - 1) / (PD_beta_Rmax$Estimate^(1 - Order.q) - 1)))
    U = beta %>% mutate(Estimate = ifelse(Order.q == 1, log(Estimate) / log(PD_beta_Rmax$Estimate), (Estimate^(Order.q - 1) - 1) / (PD_beta_Rmax$Estimate^(Order.q - 1) - 1)))
    V = beta %>% mutate(Estimate = (Estimate - 1)/(PD_beta_Rmax$Estimate - 1))
    S = beta %>% mutate(Estimate = (1/Estimate - 1)/(1/PD_beta_Rmax$Estimate - 1))


    if(nboot>1){

      #待改
      # cl = makeCluster(cluster_numbers)
      # clusterExport(cl, c("bootstrap_population_multiple_assemblage","data","data_gamma", 'data_gamma_freq',"level","N",'under_max_alpha',
      #                     'datatype', 'data_2D'))
      # clusterEvalQ(cl, library(tidyverse, magrittr))

      # plan(sequential)
      # plan(multiprocess)

      # se = parSapply(cl, 1:nboot, function(i){

      # start = Sys.time()
      se = future_lapply(1:nboot, function(i){


        tree_bt = PDtree

        bootstrap_population = bootstrap_population_multiple_assemblage(data, data_gamma, 'abundance')
        p_bt = bootstrap_population
        unseen_p = p_bt[-(1:nrow(data)),] %>% matrix(ncol = ncol(data))

        if ( nrow(p_bt) > nrow(data) & sum(unseen_p) > 0 ){

          unseen = unseen_p[which(rowSums(unseen_p) > 0),]
          unseen = matrix(unseen, ncol = ncol(unseen_p))
          p_bt = rbind(p_bt[(1:nrow(data)),], unseen)
          unseen_name = sapply(1:nrow(unseen), function(i) paste0('unseen_', i))
          rownames(p_bt) = c(rownames(data), unseen_name)

          bootstrap_sample = sapply(1:ncol(data), function(k) rmultinom(n = 1, size = sum(data[,k]), prob = p_bt[,k]))
          x_bt = bootstrap_sample

          rownames(x_bt) = rownames(p_bt)

          if ( sum(x_bt[-(1:nrow(data)),])>0 ){

            g0_hat = apply(data, 2, function(x){

              n = sum(x)
              f1 = sum(x == 1)
              f2 = sum(x == 2)

              aL = phyBranchAL_Abu(phylo = PDtree, data = x, rootExtend = T, refT = reft)

              aL$treeNabu$branch.length = aL$BLbyT[,1]
              aL = aL$treeNabu %>% select(branch.abun,branch.length)
              g1 = aL$branch.length[aL$branch.abun == 1] %>% sum
              g2 = aL$branch.length[aL$branch.abun == 2] %>% sum
              g0_hat = ifelse( g2 > ((g1*f2)/(2*f1)) , ((n-1)/n)*(g1^2/(2*g2)) , ((n-1)/n)*(g1*(f1-1)/(2*(f2+1))) )
              if(is.na(g0_hat)) {g0_hat <- 0 }
              g0_hat

            })

            te = (x_bt[1:nrow(data),]*(data == 0))>0
            used_length = sapply(1:ncol(data), function(i) {

              if (sum(te[,i]) == 0) return(0) else {

                phyBranchAL_Abu(phylo = PDtree, data = x_bt[1:nrow(data),i], rootExtend = T, refT = reft)$treeNabu %>%
                  subset(label %in% names(which(te[,i] == TRUE))) %>% select(branch.length) %>% sum

              }

            })

            g0_hat = g0_hat - used_length
            g0_hat[g0_hat < 0] = 0

            unseen_sample = x_bt[-(1:nrow(data)),]
            if (is.vector(unseen_sample)) unseen_sample = matrix(unseen_sample, ncol = ncol(x_bt))

            L0_hat = sapply(1:length(g0_hat), function(i) if(sum(unseen_sample[,i] > 0) > 0) (g0_hat[i] / nrow(unseen)) else 0 )

            L0_hat = rowSums((matrix(L0_hat, nrow(unseen_sample), ncol(unseen_sample), byrow = T) * unseen_sample)) / rowSums(unseen_sample)
            L0_hat[which(rowSums(unseen_sample) == 0)] = 0

            for (i in 1:length(L0_hat)){

              tip = list(edge = matrix(c(2,1),1,2),
                         tip.label = unseen_name[i],
                         edge.length = L0_hat[i],
                         Nnode = 1)
              class(tip) = "phylo"

              tree_bt = tree_bt + tip

            }

          } else {

            x_bt = x_bt[1:nrow(data),]
            p_bt = p_bt[1:nrow(data),]

          }

        } else {

          p_bt = p_bt[1:nrow(data),]
          x_bt = sapply(1:ncol(data), function(k) rmultinom(n = 1, size = sum(data[,k]), prob = p_bt[,k]))
          rownames(x_bt) = rownames(data)

        }

        bootstrap_data_gamma = rowSums(x_bt)
        bootstrap_data_gamma = bootstrap_data_gamma[bootstrap_data_gamma>0]
        bootstrap_data_alpha = as.matrix(x_bt) %>% as.vector
        bootstrap_data_alpha = bootstrap_data_alpha[bootstrap_data_alpha>0]

        Routledge_x <- x_bt[(rowSums(x_bt) != 0), ]

        m_gamma = sapply(level, function(i)coverage_to_size(bootstrap_data_gamma, i, datatype='abundance'))
        m_alpha = sapply(level, function(i) coverage_to_size(bootstrap_data_alpha, i, datatype='abundance'))

        m_Rou_alpha <- lapply(1:ncol(Routledge_x), function(j) {
          sapply(level, function(i) coverage_to_size(Routledge_x[, j], i, datatype = 'abundance'))
        })

        aL = phyBranchAL_Abu(phylo = tree_bt, data = bootstrap_data_gamma, rootExtend = T, refT = reft)
        aL$treeNabu$branch.length = aL$BLbyT[,1]
        aL_table_gamma = aL$treeNabu %>% select(branch.abun, branch.length, tgroup)

        gamma = as.vector(PhD.m.est(ai = aL_table_gamma$branch.abun, Lis = as.matrix(aL_table_gamma$branch.length), m = m_gamma, q = q, nt = n, reft = reft, cal = "PD") %>% t)

        #alpha_R
        aL_table_k = lapply(1:ncol(Routledge_x), function(k){

          x_k = Routledge_x[,k]
          names(x_k) = rownames(Routledge_x)

          aL = phyBranchAL_Abu(phylo = tree_bt, data = x_k, rootExtend = T, refT = reft)
          aL$treeNabu$branch.length = aL$BLbyT[,1]
          aL_table = aL$treeNabu %>% select(branch.abun, branch.length, tgroup)
        })


        aL_table_alpha = c()

        for (i in 1:N){

          # x = x_bt[x_bt[,i]>0,i]
          # names(x) = rownames(p_bt)[x_bt[,i]>0]

          x = x_bt[,i]
          names(x) = rownames(p_bt)
          x = x[x_bt[,i]>0]

          aL = phyBranchAL_Abu(phylo = tree_bt, data = x, rootExtend = T, refT = reft)
          aL$treeNabu$branch.length = aL$BLbyT[,1]
          aL_table = aL$treeNabu %>% select(branch.abun, branch.length, tgroup)

          aL_table_alpha = rbind(aL_table_alpha, aL_table)

        }

        # alpha = as.vector((iNEXT.3D:::PhD.m.est(ai = aL_table_alpha$branch.abun, Lis = as.matrix(aL_table_alpha$branch.length), m = m_alpha, q = q, nt = n, reft = reft, cal = "PD")/N) %>% t)
        joint = as.vector((PhD.m.est(ai = aL_table_alpha$branch.abun, Lis = as.matrix(aL_table_alpha$branch.length), m = m_alpha, q = q, nt = n, reft = reft, cal = "PD")) %>% t)

        # beta_obs = (iNEXT.3D:::PD.Tprofile(ai = aL_table_gamma$branch.abun, Lis = as.matrix(aL_table_gamma$branch.length), q = q, nt = n, cal = "PD") /
        #               (iNEXT.3D:::PD.Tprofile(ai = aL_table_alpha$branch.abun, Lis = as.matrix(aL_table_alpha$branch.length), q = q, nt = n, cal = "PD") / N)) %>% unlist()

        tmp_alphaR_list <- lapply(1:ncol(Routledge_x), function(k) {
          ai_k <- aL_table_k[[k]]$branch.abun
          Lis_k <- as.matrix(aL_table_k[[k]]$branch.length)
          nt_k <- sum(Routledge_x[,k])  # or colSums(Routledge_x)[k]

          sapply(1:length(level), function(i) {
            m_i <- m_Rou_alpha[[k]][i]

           PhD.m.est(ai = ai_k,
                                 Lis = Lis_k,
                                 m = m_i,
                                 q = q,
                                 nt = nt_k,
                                 reft = reft,
                                 cal = "PD")
          }, simplify = FALSE)
        })
        ###
        alphaR_df <- do.call(rbind, lapply(1:length(tmp_alphaR_list), function(k) {
          do.call(rbind, lapply(1:length(tmp_alphaR_list[[k]]), function(i){
            tmp_1 <- tmp_alphaR_list[[k]][[i]] %>% t %>% as.data.frame()
            colnames(tmp_1) <- q
            tmp_1 %>%
              gather(Order.q, Estimate) %>%
              mutate(
                Sample = k,
                SC = level[i],                     # coverage level
                Size = m_Rou_alpha[[k]][i]        # sample size
              )
          }))
        }))


        alpha_results <- alphaR_df %>% rename(qPD = Estimate)


        alpha_summary <- alpha_results %>%
          group_by(Order.q, SC) %>%
          summarise(
            Estimate = if (unique(Order.q) == "1") {
              exp(sum(log(qPD) * wk))
            } else {
              (sum((qPD^(1 - as.numeric(unique(Order.q)))) * wk))^(1 / (1 - as.numeric(unique(Order.q))))
            },
            .groups = "drop"
          ) |> as.data.frame()

        alpha = alpha_summary$Estimate

        if (PDtype == 'meanPD') {
          gamma = gamma/reft
          alpha = alpha/reft  #第二種寫法的alpha是還沒除的
          joint = joint/reft
        }

        beta = gamma/alpha
        PD_beta_Rmax = joint/alpha


        Order.q = rep(q, each = length(level))

        beta = data.frame(Estimate = beta, Order.q)

        C = (beta %>% mutate(Estimate = ifelse(Order.q == 1, log(Estimate) / log(PD_beta_Rmax), (Estimate^(1 - Order.q) - 1) / (PD_beta_Rmax^(1 - Order.q) - 1))))$Estimate
        U = (beta %>% mutate(Estimate = ifelse(Order.q == 1, log(Estimate) / log(PD_beta_Rmax), (Estimate^(Order.q - 1) - 1) / (PD_beta_Rmax^(Order.q - 1) - 1))))$Estimate
        V = (beta %>% mutate(Estimate = (Estimate - 1)/(PD_beta_Rmax - 1)))$Estimate
        S = (beta %>% mutate(Estimate = (1/Estimate - 1)/(1/PD_beta_Rmax - 1)))$Estimate



        beta = beta$Estimate

        cbind(gamma, alpha, beta, C, U, V, S) %>% as.matrix

        # }, simplify = "array") %>% apply(., 1:2, sd) %>% data.frame
      }, future.seed = TRUE) %>% abind(along = 3) %>% apply(1:2, sd)
      # end = Sys.time()
      # end - start

      # stopCluster(cl)
      # plan(sequential)

    } else {

      se = matrix(NA, ncol = 7, nrow = nrow(beta))
      colnames(se) = c("gamma", "alpha", "beta", "C", "U", 'V', 'S')
      se = as.data.frame(se)

    }




    se = as.data.frame(se)


    # if (PDtype == "PD") index = "PD"
    if (PDtype == "meanPD") index = "meanPD"


    assem_colnames <- grep("^Assemblage", colnames(alpha), value = TRUE)

    final_col_order <- c(
      "Dataset", "Order.q", "SC", "Size",
      assem_colnames,
      "Estimate", "Method", "s.e.", "LCL", "UCL", "Diversity"
    )

    gamma <- gamma %>%
      mutate(
        s.e. = se$gamma,
        LCL = Estimate - tmp * se$gamma,
        UCL = Estimate + tmp * se$gamma,
        Dataset = dataset_name,
        Diversity = index,
        Estimate = Estimate  # 重新命名 Estimate 為 Alpha
      ) %>%
      select(all_of(final_col_order)) %>%
      arrange(Order.q, SC)

    assem_colnames <- grep("^Assemblage", colnames(alpha), value = TRUE)


    final_col_order <- c(
      "Dataset", "Order.q", "SC", "Size",
      assem_colnames,
      "Estimate", "Method", "s.e.", "LCL", "UCL", "Diversity"
    )



    alpha <- alpha %>%
      mutate(
        s.e. = se$alpha,
        LCL = Estimate - tmp * se$alpha,
        UCL = Estimate + tmp * se$alpha,
        Dataset = dataset_name,
        Diversity = index,
        Estimate = Estimate  # 重新命名 Estimate 為 Alpha
      ) %>%
      select(all_of(final_col_order)) %>%
      arrange(Order.q, SC)


    final_col_order <- c(
      "Dataset", "Order.q", "SC", "Size",
      assem_colnames,
      "Estimate", "Method", "s.e.", "LCL", "UCL", "Diversity"
    )

    beta <- beta %>%
      mutate(
        s.e. = se$beta,
        LCL = Estimate - tmp * se$beta,
        UCL = Estimate + tmp * se$beta,
        Dataset = dataset_name,
        Diversity = index
      ) %>%
      select(all_of(final_col_order)) %>%
      arrange(Order.q, SC)
    C <- C %>%
      mutate(
        s.e. = se$C,
        LCL = Estimate - tmp * se$C,
        UCL = Estimate + tmp * se$C,
        Dataset = dataset_name,
        Diversity = index
      ) %>%
      select(all_of(final_col_order)) %>%
      arrange(Order.q, SC)
    U <- U %>%
      mutate(
        s.e. = se$U,
        LCL = Estimate - tmp * se$U,
        UCL = Estimate + tmp * se$U,
        Dataset = dataset_name,
        Diversity = index
      ) %>%
      select(all_of(final_col_order)) %>%
      arrange(Order.q, SC)
    V <- V %>%
      mutate(
        s.e. = se$V,
        LCL = Estimate - tmp * se$V,
        UCL = Estimate + tmp * se$V,
        Dataset = dataset_name,
        Diversity = index
      ) %>%
      select(all_of(final_col_order)) %>%
      arrange(Order.q, SC)
    S <- S %>%
      mutate(
        s.e. = se$S,
        LCL = Estimate - tmp * se$S,
        UCL = Estimate + tmp * se$S,
        Dataset = dataset_name,
        Diversity = index
      ) %>%
      select(all_of(final_col_order)) %>%
      arrange(Order.q, SC)


    if (trunc) {

      gamma = gamma %>% filter(!(Order.q==0 & round(Size)>2*n))

      alpha = alpha %>% filter(!(Order.q==0 & round(Size)>2*n))

      beta  = beta  %>% filter(!(Order.q==0 & round(Size)>2*n))

      C    =  C    %>% filter(!(Order.q==0 & round(Size)>2*n))

      U    =  U    %>% filter(!(Order.q==0 & round(Size)>2*n))

      V    =  V    %>% filter(!(Order.q==0 & round(Size)>2*n))

      S    =  S    %>% filter(!(Order.q==0 & round(Size)>2*n))

    }


    gamma = gamma %>% mutate(Reftime = reft)

    alpha = alpha %>% mutate(Reftime = reft)

    beta  = beta  %>% mutate(Reftime = reft)

    C     =  C    %>% mutate(Reftime = reft)

    U     =  U    %>% mutate(Reftime = reft)

    V     =  V    %>% mutate(Reftime = reft)

    S     =  S    %>% mutate(Reftime = reft)






    alpha$Method[alpha$SC == ref_gamma_max] =
      beta$Method[beta$SC == ref_gamma_max] =
      C$Method[C$SC == ref_gamma_max] =
      U$Method[U$SC == ref_gamma_max] =
      V$Method[V$SC == ref_gamma_max] =
      S$Method[S$SC == ref_gamma_max] = "Extrap_SC(2n, gamma)"

    gamma$Method[gamma$SC == ref_alpha_max] =
      alpha$Method[alpha$SC == ref_alpha_max] =
      beta$Method[beta$SC == ref_alpha_max] =
      C$Method[C$SC == ref_alpha_max] =
      U$Method[U$SC == ref_alpha_max] =
      V$Method[V$SC == ref_alpha_max] =
      S$Method[S$SC == ref_alpha_max] = "Extrap_SC(2n, alpha)"

    gamma$Method[gamma$SC == ref_gamma_max] = "Extrap_SC(2n, gamma)"


    alpha$Method[alpha$SC == ref_gamma] =
      beta$Method[beta$SC == ref_gamma] =
      C$Method[C$SC == ref_gamma] =
      U$Method[U$SC == ref_gamma] =
      V$Method[V$SC == ref_gamma] =
      S$Method[S$SC == ref_gamma] = "Observed_SC(n, gamma)"

    gamma$Method[gamma$SC == ref_alpha] =
      alpha$Method[alpha$SC == ref_alpha] =
      beta$Method[beta$SC == ref_alpha] =
      C$Method[C$SC == ref_alpha] =
      U$Method[U$SC == ref_alpha] =
      V$Method[V$SC == ref_alpha] =
      S$Method[S$SC == ref_alpha] = "Observed_SC(n, alpha)"

    gamma$Method[gamma$SC == ref_gamma] = "Observed_SC(n, gamma)"




    list(gamma = gamma, alpha = alpha, beta = beta, `1-C` = C, `1-U` = U, `1-V` = V, `1-S` = S)

  }

  for_each_dataset.size = function(data, dataset_name, N, level) {

    #data


    n = sum(data)
    data_gamma = rowSums(data)
    data_gamma = data_gamma[data_gamma>0]
    data_alpha = as.matrix(data) %>% as.vector
    Routledge_x <- data[(rowSums(data) != 0), ]
    wk <- colSums(Routledge_x) / sum(Routledge_x)
    # m_alpha_R = apply(Routledge_x,2,function(x)iNEXT.3D:::Coverage(x, "abundance", level))

    ref_gamma = n
    ref_alpha = n

    aL = phyBranchAL_Abu(phylo = PDtree, data = data_gamma, rootExtend = T, refT = reft)
    aL$treeNabu$branch.length = aL$BLbyT[,1]
    aL_table_gamma = aL$treeNabu %>% select(branch.abun, branch.length, tgroup)

    gamma = PhD.m.est(ai = aL_table_gamma$branch.abun, Lis = as.matrix(aL_table_gamma$branch.length), m = level, q = q, nt = n, reft = reft, cal = "PD") %>% t %>% as.data.frame %>%
      set_colnames(q) %>% gather(Order.q, Estimate) %>%
      mutate(Coverage_real = rep(Coverage(data_gamma, "abundance", level), length(q)), Size = rep(level, length(q)), Size = rep(level, length(q)))


    aL_table_alpha = c()

    for (i in 1:N){

      x = data[data[,i]>0,i]
      names(x) = rownames(data)[data[,i]>0]

      aL = phyBranchAL_Abu(phylo = PDtree, data = x, rootExtend = T, refT = reft)
      aL$treeNabu$branch.length = aL$BLbyT[,1]
      aL_table = aL$treeNabu %>% select(branch.abun, branch.length, tgroup)

      aL_table_alpha = rbind(aL_table_alpha, aL_table)

    }


    qPDm = PhD.m.est(ai = aL_table_alpha$branch.abun, Lis = as.matrix(aL_table_alpha$branch.length), m = level, q = q, nt = n, reft = reft, cal = "PD")
    # qPDm = qPDm/N
    # alpha = qPDm %>% t %>% as.data.frame %>%
    #   set_colnames(q) %>% gather(Order.q, Estimate) %>%
    #   mutate(Coverage_real = rep(iNEXT.3D:::Coverage(data_alpha, "abundance", level), length(q)), Size = rep(level, length(q)))
    #


    joint = qPDm %>% t %>% as.data.frame %>%
      set_colnames(q) %>% gather(Order.q, Estimate) %>%
      mutate(Coverage_real = rep(Coverage(data_alpha, "abundance", level), length(q)), Size = rep(level, length(q)))


    # gamma = (gamma %>%
    #            mutate(Method = ifelse(Size >= ref_gamma, ifelse(Size == ref_gamma, 'Observed', 'Extrapolation'), 'Rarefaction')))[,c(2,1,5,3,4)] %>%
    #   set_colnames(c('Estimate', 'Order.q', 'Method', 'SC', 'Size'))
    #
    # gamma$Order.q = as.numeric(gamma$Order.q)


    # alpha = (alpha %>%
    #            mutate(Method = ifelse(Size >= ref_alpha, ifelse(Size == ref_alpha, 'Observed', 'Extrapolation'), 'Rarefaction')))[,c(2,1,5,3,4)] %>%
    #   set_colnames(c('Estimate', 'Order.q', 'Method', 'SC', 'Size'))
    #
    # alpha$Order.q = as.numeric(alpha$Order.q)


    #alpha

    aL_table_k = lapply(1:ncol(Routledge_x), function(k){

      x_k = Routledge_x[,k]
      names(x_k) = rownames(Routledge_x)

      aL = phyBranchAL_Abu(phylo = PDtree, data = x_k, rootExtend = T, refT = reft)
      aL$treeNabu$branch.length = aL$BLbyT[,1]
      aL_table = aL$treeNabu %>% select(branch.abun, branch.length, tgroup)
    })



    tmp_alphaR_list <- lapply(1:ncol(Routledge_x), function(k) {
      ai_k <- aL_table_k[[k]]$branch.abun
      Lis_k <- as.matrix(aL_table_k[[k]]$branch.length)
      nt_k <- sum(Routledge_x[, k])

      sapply(1:length(level), function(i) {
        PhD.m.est(
          ai = ai_k,
          Lis = Lis_k,
          m = level[i],
          q = q,
          nt = nt_k,
          reft = reft,
          cal = "PD"
        )
      }, simplify = FALSE)
    })


    alphaR_df <- do.call(rbind, lapply(1:length(tmp_alphaR_list), function(k) {
      do.call(rbind, lapply(1:length(tmp_alphaR_list[[k]]), function(i){
        tmp <- tmp_alphaR_list[[k]][[i]] %>% t %>% as.data.frame()
        colnames(tmp) <- q
        tmp %>%
          gather(Order.q, Estimate) %>%
          mutate(
            Sample = k,
            Size = level[i]
          )
      }))
    }))

    # Rename Estimate to qPD for weighted calculation
    alpha_results <- alphaR_df %>% rename(qPD = Estimate)

    # Weighted mean calculation
    alpha_summary <- alpha_results %>%
      group_by(Order.q, Size) %>%
      summarise(
        Estimate = if (unique(Order.q) == "1") {
          exp(sum(log(qPD) * wk))
        } else {
          (sum((qPD^(1 - as.numeric(unique(Order.q)))) * wk))^(1 / (1 - as.numeric(unique(Order.q))))  # 加權 power mean
        },
        .groups = "drop"
      )

    coverage_list <- lapply(1:ncol(Routledge_x), function(j) {
      Coverage(Routledge_x[, j], datatype = "abundance", m = level)
    })
    names(coverage_list) <- paste0("Coverage_Assemblage", 1:ncol(Routledge_x))


    coverage_df <- as.data.frame(coverage_list)
    coverage_df_expanded <- coverage_df[rep(1:nrow(coverage_df), length(q)), ]


    alpha <- bind_cols(alpha_summary, coverage_df_expanded) |> as.data.frame()

    alpha <- alpha %>%
      arrange(Order.q, Size) %>%
      mutate(SC = rep(Coverage(data_alpha, "abundance", level), length(q)))


    coverage_cols <- grep("^Coverage_Assemblage", colnames(alpha), value = TRUE)


    final_col_order <- c("Estimate", "Order.q", "Method", "SC", coverage_cols, "Size")


    alpha <- (alpha %>%mutate(Method = ifelse(Size >= ref_alpha, ifelse(Size == ref_alpha, 'Observed', 'Extrapolation'), 'Rarefaction')))    %>%
      select(all_of(final_col_order))

    # alpha = (alpha %>%
    #            mutate(Method = ifelse(Size >= ref_alpha, ifelse(Size == ref_alpha, 'Observed', 'Extrapolation'), 'Rarefaction')))[,c(2,1,5,3,4)] %>%
    #   set_colnames(c('Estimate', 'Order.q', 'Method', 'SC', 'Size'))
    #
    alpha$Order.q = as.numeric(alpha$Order.q)

    final_col_order <- c("Estimate", "Order.q", "Method", "Coverage_real", coverage_cols, "Size")


    #gamma
    gamma = (gamma %>%
               mutate(Method = ifelse(Size >= ref_gamma, ifelse(Size == ref_gamma, 'Observed', 'Extrapolation'), 'Rarefaction')))[,c(2,1,5,3,4)] %>%
      bind_cols(coverage_df_expanded) %>%
      select(all_of(final_col_order))

    gamma$Order.q = as.numeric(gamma$Order.q)
    colnames(gamma) = c('Estimate', 'Order.q', 'Method', 'SC', coverage_cols, 'Size')

    #joint
    final_col_order <- c("Estimate", "Order.q", "Method", "Coverage_real", coverage_cols, "Size")

    joint = (joint %>%
               mutate(Method = ifelse(Size >= ref_alpha, ifelse(Size == ref_alpha, 'Observed', 'Extrapolation'), 'Rarefaction')))[,c(2,1,5,3,4)] %>%
      bind_cols(coverage_df_expanded) %>%
      select(all_of(final_col_order))

    joint$Order.q = as.numeric(joint$Order.q)
    colnames(joint) = c('Estimate', 'Order.q', 'Method', 'SC', coverage_cols, 'Size')


    if (PDtype == 'meanPD') {
      gamma$Estimate = gamma$Estimate/reft
      alpha$Estimate = alpha$Estimate/reft
      joint$Estimate = joint$Estimate/reft
    }

    # beta = alpha
    # beta$Estimate = gamma$Estimate/alpha$Estimate
    #
    # C = beta %>% mutate(Estimate = ifelse(Order.q == 1, log(Estimate)/log(N), (Estimate^(1 - Order.q) - 1)/(N^(1 - Order.q) - 1)))
    # U = beta %>% mutate(Estimate = ifelse(Order.q == 1, log(Estimate)/log(N), (Estimate^(Order.q - 1) - 1)/(N^(Order.q - 1) - 1)))
    # V = beta %>% mutate(Estimate = (Estimate - 1)/(N - 1))
    # S = beta %>% mutate(Estimate = (1/Estimate - 1)/(1/N - 1))

    if(nboot>1){

      # cl = makeCluster(cluster_numbers)
      # clusterExport(cl, c("bootstrap_population_multiple_assemblage","data","data_gamma", 'data_gamma_freq',"level","N",'under_max_alpha',
      #                     'datatype', 'data_2D'))
      # clusterEvalQ(cl, library(tidyverse, magrittr))

      # plan(sequential)
      # plan(multiprocess)

      # se = parSapply(cl, 1:nboot, function(i){

      # start = Sys.time()
      se = future_lapply(1:nboot, function(i){



        tree_bt = PDtree

        bootstrap_population = bootstrap_population_multiple_assemblage(data, data_gamma, 'abundance')
        p_bt = bootstrap_population
        unseen_p = p_bt[-(1:nrow(data)),] %>% matrix(ncol = ncol(data))

        if ( nrow(p_bt) > nrow(data) & sum(unseen_p) > 0 ){

          unseen = unseen_p[which(rowSums(unseen_p) > 0),]
          unseen = matrix(unseen, ncol = ncol(unseen_p))
          p_bt = rbind(p_bt[(1:nrow(data)),], unseen)
          unseen_name = sapply(1:nrow(unseen), function(i) paste0('unseen_', i))
          rownames(p_bt) = c(rownames(data), unseen_name)

          bootstrap_sample = sapply(1:ncol(data), function(k) rmultinom(n = 1, size = sum(data[,k]), prob = p_bt[,k]))
          x_bt = bootstrap_sample

          rownames(x_bt) = rownames(p_bt)

          if ( sum(x_bt[-(1:nrow(data)),])>0 ){

            g0_hat = apply(data, 2, function(x){

              n = sum(x)
              f1 = sum(x == 1)
              f2 = sum(x == 2)

              aL = phyBranchAL_Abu(phylo = PDtree, data = x, rootExtend = T, refT = reft)

              aL$treeNabu$branch.length = aL$BLbyT[,1]
              aL = aL$treeNabu %>% select(branch.abun,branch.length)
              g1 = aL$branch.length[aL$branch.abun == 1] %>% sum
              g2 = aL$branch.length[aL$branch.abun == 2] %>% sum
              g0_hat = ifelse( g2 > ((g1*f2)/(2*f1)) , ((n-1)/n)*(g1^2/(2*g2)) , ((n-1)/n)*(g1*(f1-1)/(2*(f2+1))) )
              g0_hat

            })

            te = (x_bt[1:nrow(data),]*(data == 0))>0
            used_length = sapply(1:ncol(data), function(i) {

              if (sum(te[,i]) == 0) return(0) else {

                phyBranchAL_Abu(phylo = PDtree, data = x_bt[1:nrow(data),i], rootExtend = T, refT = reft)$treeNabu %>%
                  subset(label %in% names(which(te[,i] == TRUE))) %>% select(branch.length) %>% sum

              }

            })

            g0_hat = g0_hat - used_length
            g0_hat[g0_hat < 0] = 0

            unseen_sample = x_bt[-(1:nrow(data)),]
            if (is.vector(unseen_sample)) unseen_sample = matrix(unseen_sample, ncol = ncol(x_bt))

            L0_hat = sapply(1:length(g0_hat), function(i) if(sum(unseen_sample[,i] > 0) > 0) (g0_hat[i] / nrow(unseen)) else 0 )

            L0_hat = rowSums((matrix(L0_hat, nrow(unseen_sample), ncol(unseen_sample), byrow = T) * unseen_sample)) / rowSums(unseen_sample)
            L0_hat[which(rowSums(unseen_sample) == 0)] = 0

            for (i in 1:length(L0_hat)){

              tip = list(edge = matrix(c(2,1),1,2),
                         tip.label = unseen_name[i],
                         edge.length = L0_hat[i],
                         Nnode = 1)
              class(tip) = "phylo"

              tree_bt = tree_bt + tip

            }

          } else {

            x_bt = x_bt[1:nrow(data),]
            p_bt = p_bt[1:nrow(data),]

          }

        } else {

          p_bt = p_bt[1:nrow(data),]
          x_bt = sapply(1:ncol(data), function(k) rmultinom(n = 1, size = sum(data[,k]), prob = p_bt[,k]))
          rownames(x_bt) = rownames(data)

        }

        bootstrap_data_gamma = rowSums(x_bt)
        bootstrap_data_gamma = bootstrap_data_gamma[bootstrap_data_gamma>0]
        bootstrap_data_alpha = as.matrix(x_bt) %>% as.vector
        bootstrap_data_alpha = bootstrap_data_alpha[bootstrap_data_alpha>0]
        Routledge_x <- x_bt[(rowSums(x_bt) != 0), ]

        aL = phyBranchAL_Abu(phylo = tree_bt, data = bootstrap_data_gamma, rootExtend = T, refT = reft)
        aL$treeNabu$branch.length = aL$BLbyT[,1]
        aL_table_gamma = aL$treeNabu %>% select(branch.abun, branch.length, tgroup)

        gamma = as.vector(PhD.m.est(ai = aL_table_gamma$branch.abun, Lis = as.matrix(aL_table_gamma$branch.length), m = level, q = q, nt = n, reft = reft, cal = "PD") %>% t)


        aL_table_alpha = c()

        for (i in 1:N){

          # x = x_bt[x_bt[,i]>0,i]
          # names(x) = rownames(p_bt)[x_bt[,i]>0]

          x = x_bt[,i]
          names(x) = rownames(p_bt)
          x = x[x_bt[,i]>0]

          aL = phyBranchAL_Abu(phylo = tree_bt, data = x, rootExtend = T, refT = reft)
          aL$treeNabu$branch.length = aL$BLbyT[,1]
          aL_table = aL$treeNabu %>% select(branch.abun, branch.length, tgroup)

          aL_table_alpha = rbind(aL_table_alpha, aL_table)

        }

        # alpha = as.vector((iNEXT.3D:::PhD.m.est(ai = aL_table_alpha$branch.abun, Lis = as.matrix(aL_table_alpha$branch.length), m = level, q = q, nt = n, reft = reft, cal = "PD")/N) %>% t)

        joint = as.vector((PhD.m.est(ai = aL_table_alpha$branch.abun, Lis = as.matrix(aL_table_alpha$branch.length), m = level, q = q, nt = n, reft = reft, cal = "PD")) %>% t)

        #alpha
        aL_table_k = lapply(1:ncol(Routledge_x), function(k){

          x_k = Routledge_x[,k]
          names(x_k) = rownames(Routledge_x)

          aL = phyBranchAL_Abu(phylo = tree_bt, data = x_k, rootExtend = T, refT = reft)
          aL$treeNabu$branch.length = aL$BLbyT[,1]
          aL_table = aL$treeNabu %>% select(branch.abun, branch.length, tgroup)
        })



        tmp_alphaR_list <- lapply(1:ncol(Routledge_x), function(k) {
          ai_k <- aL_table_k[[k]]$branch.abun
          Lis_k <- as.matrix(aL_table_k[[k]]$branch.length)
          nt_k <- sum(Routledge_x[, k])

          sapply(1:length(level), function(i) {
            PhD.m.est(
              ai = ai_k,
              Lis = Lis_k,
              m = level[i],
              q = q,
              nt = nt_k,
              reft = reft,
              cal = "PD"
            )
          }, simplify = FALSE)
        })


        alphaR_df <- do.call(rbind, lapply(1:length(tmp_alphaR_list), function(k) {
          do.call(rbind, lapply(1:length(tmp_alphaR_list[[k]]), function(i){
            tmp <- tmp_alphaR_list[[k]][[i]] %>% t %>% as.data.frame()
            colnames(tmp) <- q
            tmp %>%
              gather(Order.q, Estimate) %>%
              mutate(
                Sample = k,
                Size = level[i]
              )
          }))
        }))

        # Rename Estimate to qPD for weighted calculation
        alpha_results <- alphaR_df %>% rename(qPD = Estimate)

        # Weighted mean calculation
        alpha_summary <- alpha_results %>%
          group_by(Order.q, Size) %>%
          summarise(
            Estimate = if (unique(Order.q) == "1") {
              exp(sum(log(qPD) * wk))  # 加權幾何平均 for q=1
            } else {
              (sum((qPD^(1 - as.numeric(unique(Order.q)))) * wk))^(1 / (1 - as.numeric(unique(Order.q))))  # 加權 power mean
            },
            .groups = "drop"
          )

        coverage_list <- lapply(1:ncol(Routledge_x), function(j) {
          Coverage(Routledge_x[, j], datatype = "abundance", m = level)
        })
        names(coverage_list) <- paste0("Coverage_Assemblage", 1:ncol(Routledge_x))


        coverage_df <- as.data.frame(coverage_list)
        coverage_df_expanded <- coverage_df[rep(1:nrow(coverage_df), length(q)), ]


        alpha <- bind_cols(alpha_summary, coverage_df_expanded) |> as.data.frame()

        alpha <- alpha %>%
          arrange(Order.q, Size) %>%
          mutate(SC = rep(Coverage(data_alpha, "abundance", level), length(q)))


        coverage_cols <- grep("^Coverage_Assemblage", colnames(alpha), value = TRUE)


        final_col_order <- c("Estimate", "Order.q", "Method", "SC", coverage_cols, "Size")



        alpha <- (alpha %>%mutate(Method = ifelse(Size >= ref_alpha, ifelse(Size == ref_alpha, 'Observed', 'Extrapolation'), 'Rarefaction')))    %>%
          select(all_of(final_col_order))

        alpha = alpha$Estimate

        if (PDtype == 'meanPD') {
          gamma = gamma/reft
          alpha = alpha/reft
          joint = joint/reft
        }

        # beta = gamma/alpha
        #
        # Order.q = rep(q, each = length(level))
        #
        # beta = data.frame(Estimate = beta, Order.q)
        #
        # C = (beta %>% mutate(Estimate = ifelse(Order.q == 1, log(Estimate)/log(N), (Estimate^(1 - Order.q) - 1)/(N^(1 - Order.q) - 1))))$Estimate
        # U = (beta %>% mutate(Estimate = ifelse(Order.q == 1, log(Estimate)/log(N), (Estimate^(Order.q - 1) - 1)/(N^(Order.q - 1) - 1))))$Estimate
        # V = (beta %>% mutate(Estimate = (Estimate - 1)/(N - 1)))$Estimate
        # S = (beta %>% mutate(Estimate = (1/Estimate - 1)/(1/N - 1)))$Estimate
        #
        # beta = beta$Estimate
        #
        # cbind(gamma, alpha, beta, C, U, V, S) %>% as.matrix
        cbind(gamma, alpha,joint) %>% as.matrix

        # }, simplify = "array") %>% apply(., 1:2, sd) %>% data.frame
      }, future.seed = TRUE) %>% abind(along = 3) %>% apply(1:2, sd)
      # end = Sys.time()
      # end - start

      # stopCluster(cl)
      # plan(sequential)

    } else {

      # se = matrix(0, ncol = 7, nrow = nrow(beta))
      # colnames(se) = c("gamma", "alpha", "beta", "C", "U", 'V', 'S')
      # se = as.data.frame(se)

      se = matrix(NA, ncol = 3, nrow = nrow(gamma))
      colnames(se) = c("gamma", "alpha","joint")
      se = as.data.frame(se)
      #
    }




    se = as.data.frame(se)

    # if (diversity == "TD") index = "TD"
    # if (diversity == "PD" & PDtype == "PD") index = "PD"
    if (PDtype == "meanPD") index = "meanPD"

    assem_colnames <- grep("^Coverage_Assemblage", colnames(alpha), value = TRUE)


    final_col_order <- c(
      "Dataset", "Order.q", "Size", "SC",
      assem_colnames,
      "Gamma", "Method", "s.e.", "LCL", "UCL", "Diversity"
    )


    gamma <- gamma %>%
      mutate(
        s.e. = se$gamma,
        LCL = Estimate - tmp * se$gamma,
        UCL = Estimate + tmp * se$gamma,
        Dataset = dataset_name,
        Diversity = index,
        Gamma = Estimate  # 重新命名 Estimate 為 Alpha
      ) %>%
      select(all_of(final_col_order)) %>%
      arrange(Order.q, SC)

    # gamma = gamma %>% mutate(s.e. = se$gamma,
    #                          LCL = Estimate - tmp * se$gamma,
    #                          UCL = Estimate + tmp * se$gamma,
    #                          Dataset = dataset_name,
    #                          Diversity = index) %>%
    #   arrange(Order.q, Size) %>% .[,c(9, 2, 5, 4, 1, 3, 6, 7, 8, 10)] %>% rename("Gamma" = "Estimate")
    #
    # alpha = alpha %>% mutate(s.e. = se$alpha,
    #                          LCL = Estimate - tmp * se$alpha,
    #                          UCL = Estimate + tmp * se$alpha,
    #                          Dataset = dataset_name,
    #                          Diversity = index) %>%
    #   arrange(Order.q, Size) %>% .[,c(9, 2, 5, 4, 1, 3, 6, 7, 8, 10)] %>% rename("Alpha" = "Estimate")
    #
    assem_colnames <- grep("^Coverage_Assemblage", colnames(alpha), value = TRUE)


    final_col_order <- c(
      "Dataset", "Order.q", "Size", "SC",
      assem_colnames,
      "Alpha", "Method", "s.e.", "LCL", "UCL", "Diversity"
    )


    alpha <- alpha %>%
      mutate(
        s.e. = se$alpha,
        LCL = Estimate - tmp * se$alpha,
        UCL = Estimate + tmp * se$alpha,
        Dataset = dataset_name,
        Diversity = index,
        Alpha = Estimate
      ) %>%
      select(all_of(final_col_order)) %>%
      arrange(Order.q, SC)


    final_col_order <- c(
      "Dataset", "Order.q", "Size", "SC",
      assem_colnames,
      "Estimate", "Method", "s.e.", "LCL", "UCL", "Diversity"
    )

    joint <- joint %>%
      mutate(
        s.e. = se$joint,
        LCL = Estimate - tmp * se$joint,
        UCL = Estimate + tmp * se$joint,
        Dataset = dataset_name,
        Diversity = index
      ) %>%
      select(all_of(final_col_order)) %>%
      arrange(Order.q, SC)


    gamma = gamma %>% mutate(Reftime = reft)

    alpha = alpha %>% mutate(Reftime = reft)

    joint = joint %>% mutate(Reftime = reft)






    list(gamma = gamma, alpha = alpha, joint = joint)

  }


  if (base == 'coverage') output = lapply(1:length(data_list), function(i) for_each_dataset(data = data_list[[i]], dataset_name = dataset_names[i], N = Ns[i], level = level[[i]]))

  if (base == 'size') output = lapply(1:length(data_list), function(i) for_each_dataset.size(data = data_list[[i]], dataset_name = dataset_names[i], N = Ns[i], level = level[[i]]))

  names(output) = dataset_names

  class(output) <- c("iNEXT_seq_Relative")
  return(output)
}



#' Plot Relative Gamma, Alpha, Beta or Dissimilarity from iNEXT_seq_Relative output
#'
#' This function creates ggplot-based plots for the output of \code{iNEXT_seq_Relative}.
#'
#' @param output Output list from \code{iNEXT_seq_Relative}
#' @param type Type of plot. "B" for diversity (gamma, alpha, beta), "D" for dissimilarity (1-C, 1-U, 1-V, 1-S).
#'
#' @return A ggplot object showing diversity or dissimilarity profiles.
#'
#' @examples
#' data("fungi")
#' data("fungi_tree")
#' output <- iNEXT_seq_Relative(fungi[1], q = c(0,1,2), nboot = 10, PDtree = fungi_tree)
#' ggiNEXT_seq_Relative(output, type = "B")
#' @export
ggiNEXT_seq_Relative = function(output, type = 'B'){

  transp = 0.4

  # Check if the number of unique 'Assemblage' is 8 or less
  if (length(output) <= 8){
    cbPalette <- rev(c("#999999", "#E69F00", "#56B4E9", "#009E73",
                       "#330066", "#CC79A7", "#0072B2", "#D55E00"))
  }else{
    # If there are more than 8 assemblages, start with the same predefined color palette
    # Then extend the palette by generating additional colors using the 'ggplotColors' function
    cbPalette <- rev(c("#999999", "#E69F00", "#56B4E9", "#009E73",
                       "#330066", "#CC79A7", "#0072B2", "#D55E00"))
    cbPalette <- c(cbPalette, ggplotColors(length(output)-8))
  }

  if ((length(output[[1]]) == 7 & type == "B") | (length(output[[1]]) == 3)) {

    if (unique(output[[1]]$gamma$Diversity) == 'TD') { ylab = "Taxonomic diversity" }
    if (unique(output[[1]]$gamma$Diversity) == 'PD') { ylab = "Phylogenetic diversity" }
    if (unique(output[[1]]$gamma$Diversity) == 'meanPD') { ylab = "Mean phylogenetic diversity" }
    if (unique(output[[1]]$gamma$Diversity) == 'FD_tau') { ylab = "Functional diversity (given tau)" }
    if (unique(output[[1]]$gamma$Diversity) == 'FD_AUC') { ylab = "Functional diversity (AUC)" }
  }

  if (length(output[[1]]) == 7 & type == "D") {

    if (unique(output[[1]]$gamma$Diversity) == 'TD') { ylab = "Taxonomic dissimilarity" }
    if (unique(output[[1]]$gamma$Diversity) == 'PD') { ylab = "Phylogenetic dissimilarity" }
    if (unique(output[[1]]$gamma$Diversity) == 'meanPD') { ylab = "Mean phylogenetic dissimilarity" }
    if (unique(output[[1]]$gamma$Diversity) == 'FD_tau') { ylab = "Functional dissimilarity (given tau)" }
    if (unique(output[[1]]$gamma$Diversity) == 'FD_AUC') { ylab = "Functional dissimilarity (AUC)" }
  }

  if (length(output[[1]]) == 7) {
    if (type == 'B'){

      # gamma = lapply(output, function(y) {
      #
      #   tmp = y[["gamma"]]
      #   if ('mT' %in% colnames(tmp)) tmp = tmp %>% rename("Size" = "mT")
      #
      #   tmp[tmp == 'Observed_SC(n, gamma)'] = tmp[tmp == 'Observed_SC(T, gamma)'] = 'Observed'
      #   tmp[tmp == 'Extrap_SC(2n, gamma)'] = tmp[tmp == 'Extrap_SC(2T, gamma)'] = 'Extrapolation'
      #
      #   tmp[tmp == 'Observed_SC(n, alpha)'] = ifelse( unique(tmp[tmp$Method == 'Observed_SC(n, alpha)',"Size"] < unique(tmp[tmp$Method == "Observed", "Size"]) ),
      #                                                 "Rarefaction", "Extrapolation")
      #   tmp[tmp == 'Extrap_SC(2n, alpha)'] = ifelse( unique(tmp[tmp$Method == 'Extrap_SC(2n, alpha)',"Size"] < unique(tmp[tmp$Method == "Observed", "Size"]) ),
      #                                                "Rarefaction", "Extrapolation")
      #
      #   tmp[tmp == 'Observed_SC(T, alpha)'] = ifelse( unique(tmp[tmp$Method == 'Observed_SC(T, alpha)',"Size"] < unique(tmp[tmp$Method == "Observed", "Size"]) ),
      #                                                 "Rarefaction", "Extrapolation")
      #   tmp[tmp == 'Extrap_SC(2T, alpha)'] = ifelse( unique(tmp[tmp$Method == 'Extrap_SC(2T, alpha)',"Size"] < unique(tmp[tmp$Method == "Observed", "Size"]) ),
      #                                                "Rarefaction", "Extrapolation")
      #
      #   tmp
      #
      # }) %>% do.call(rbind,.) %>% rename("Estimate" = "Gamma") %>% mutate(div_type = "Gamma") %>% as_tibble()
      #
      gamma = lapply(output, function(y) {

        tmp = y[["gamma"]]
        if ('mT' %in% colnames(tmp)) tmp = tmp %>% rename("Size" = "mT")

        tmp[tmp == 'Observed_SC(n, gamma)'] = tmp[tmp == 'Observed_SC(T, gamma)'] = 'Observed'
        tmp[tmp == 'Extrap_SC(2n, gamma)'] = tmp[tmp == 'Extrap_SC(2T, gamma)'] = 'Extrapolation'

        tmp[tmp == 'Observed_SC(n, alpha)'] = ifelse( unique(tmp[tmp$Method == 'Observed_SC(n, alpha)',"Size"] < unique(tmp[tmp$Method == "Observed", "Size"]) ),
                                                      "Rarefaction", "Extrapolation")
        tmp[tmp == 'Extrap_SC(2n, alpha)'] = ifelse( unique(tmp[tmp$Method == 'Extrap_SC(2n, alpha)',"Size"] < unique(tmp[tmp$Method == "Observed", "Size"]) ),
                                                     "Rarefaction", "Extrapolation")

        tmp[tmp == 'Observed_SC(T, alpha)'] = ifelse( unique(tmp[tmp$Method == 'Observed_SC(T, alpha)',"Size"] < unique(tmp[tmp$Method == "Observed", "Size"]) ),
                                                      "Rarefaction", "Extrapolation")
        tmp[tmp == 'Extrap_SC(2T, alpha)'] = ifelse( unique(tmp[tmp$Method == 'Extrap_SC(2T, alpha)',"Size"] < unique(tmp[tmp$Method == "Observed", "Size"]) ),
                                                     "Rarefaction", "Extrapolation")

        tmp

      }) %>% do.call(rbind,.) %>% mutate(div_type = "Gamma") %>% as_tibble()

      # alpha = lapply(output, function(y) {
      #
      #   tmp = y[["alpha"]]
      #   if ('mT' %in% colnames(tmp)) tmp = tmp %>% rename("Size" = "mT")
      #
      #   tmp[tmp == 'Observed_SC(n, alpha)'] = tmp[tmp == 'Observed_SC(T, alpha)'] = 'Observed'
      #   tmp[tmp == 'Extrap_SC(2n, alpha)'] = tmp[tmp == 'Extrap_SC(2T, alpha)'] = 'Extrapolation'
      #
      #   tmp[tmp == 'Observed_SC(n, gamma)'] = ifelse( unique(tmp[tmp$Method == 'Observed_SC(n, gamma)',"Size"] < unique(tmp[tmp$Method == "Observed", "Size"]) ),
      #                                                 "Rarefaction", "Extrapolation")
      #   tmp[tmp == 'Extrap_SC(2n, gamma)'] = ifelse( unique(tmp[tmp$Method == 'Extrap_SC(2n, gamma)',"Size"] < unique(tmp[tmp$Method == "Observed", "Size"]) ),
      #                                                "Rarefaction", "Extrapolation")
      #
      #   tmp[tmp == 'Observed_SC(T, gamma)'] = ifelse( unique(tmp[tmp$Method == 'Observed_SC(T, gamma)',"Size"] < unique(tmp[tmp$Method == "Observed", "Size"]) ),
      #                                                 "Rarefaction", "Extrapolation")
      #   tmp[tmp == 'Extrap_SC(2T, gamma)'] = ifelse( unique(tmp[tmp$Method == 'Extrap_SC(2T, gamma)',"Size"] < unique(tmp[tmp$Method == "Observed", "Size"]) ),
      #                                                "Rarefaction", "Extrapolation")
      #
      #   tmp
      #
      # }) %>% do.call(rbind,.) %>% rename("Estimate" = "Alpha") %>% mutate(div_type = "Alpha") %>% as_tibble()
      #
      alpha = lapply(output, function(y) {

        tmp = y[["alpha"]]
        if ('mT' %in% colnames(tmp)) tmp = tmp %>% rename("Size" = "mT")

        tmp[tmp == 'Observed_SC(n, alpha)'] = tmp[tmp == 'Observed_SC(T, alpha)'] = 'Observed'
        tmp[tmp == 'Extrap_SC(2n, alpha)'] = tmp[tmp == 'Extrap_SC(2T, alpha)'] = 'Extrapolation'

        tmp[tmp == 'Observed_SC(n, gamma)'] = ifelse( unique(tmp[tmp$Method == 'Observed_SC(n, gamma)',"Size"] < unique(tmp[tmp$Method == "Observed", "Size"]) ),
                                                      "Rarefaction", "Extrapolation")
        tmp[tmp == 'Extrap_SC(2n, gamma)'] = ifelse( unique(tmp[tmp$Method == 'Extrap_SC(2n, gamma)',"Size"] < unique(tmp[tmp$Method == "Observed", "Size"]) ),
                                                     "Rarefaction", "Extrapolation")

        tmp[tmp == 'Observed_SC(T, gamma)'] = ifelse( unique(tmp[tmp$Method == 'Observed_SC(T, gamma)',"Size"] < unique(tmp[tmp$Method == "Observed", "Size"]) ),
                                                      "Rarefaction", "Extrapolation")
        tmp[tmp == 'Extrap_SC(2T, gamma)'] = ifelse( unique(tmp[tmp$Method == 'Extrap_SC(2T, gamma)',"Size"] < unique(tmp[tmp$Method == "Observed", "Size"]) ),
                                                     "Rarefaction", "Extrapolation")

        tmp

      }) %>% do.call(rbind,.)  %>% mutate(div_type = "Alpha") %>% as_tibble()


      # beta =  lapply(output, function(y) {
      #
      #   tmp = y[["beta"]]
      #   if ('mT' %in% colnames(tmp)) tmp = tmp %>% rename("Size" = "mT")
      #
      #   tmp[tmp == 'Observed_SC(n, alpha)'] = tmp[tmp == 'Observed_SC(T, alpha)'] = 'Observed'
      #   tmp[tmp == 'Extrap_SC(2n, alpha)'] = tmp[tmp == 'Extrap_SC(2T, alpha)'] = 'Extrapolation'
      #
      #   tmp[tmp == 'Observed_SC(n, gamma)'] = ifelse( unique(tmp[tmp$Method == 'Observed_SC(n, gamma)',"Size"] < unique(tmp[tmp$Method == "Observed", "Size"]) ),
      #                                                 "Rarefaction", "Extrapolation")
      #   tmp[tmp == 'Extrap_SC(2n, gamma)'] = ifelse( unique(tmp[tmp$Method == 'Extrap_SC(2n, gamma)',"Size"] < unique(tmp[tmp$Method == "Observed", "Size"]) ),
      #                                                "Rarefaction", "Extrapolation")
      #
      #   tmp[tmp == 'Observed_SC(T, gamma)'] = ifelse( unique(tmp[tmp$Method == 'Observed_SC(T, gamma)',"Size"] < unique(tmp[tmp$Method == "Observed", "Size"]) ),
      #                                                 "Rarefaction", "Extrapolation")
      #   tmp[tmp == 'Extrap_SC(2T, gamma)'] = ifelse( unique(tmp[tmp$Method == 'Extrap_SC(2T, gamma)',"Size"] < unique(tmp[tmp$Method == "Observed", "Size"]) ),
      #                                                "Rarefaction", "Extrapolation")
      #
      #   tmp
      #
      # })  %>% do.call(rbind,.) %>% rename("Estimate" = "Beta") %>% mutate(div_type = "Beta")  %>% as_tibble()
      #
      beta =  lapply(output, function(y) {

        tmp = y[["beta"]]
        if ('mT' %in% colnames(tmp)) tmp = tmp %>% rename("Size" = "mT")

        tmp[tmp == 'Observed_SC(n, alpha)'] = tmp[tmp == 'Observed_SC(T, alpha)'] = 'Observed'
        tmp[tmp == 'Extrap_SC(2n, alpha)'] = tmp[tmp == 'Extrap_SC(2T, alpha)'] = 'Extrapolation'

        tmp[tmp == 'Observed_SC(n, gamma)'] = ifelse( unique(tmp[tmp$Method == 'Observed_SC(n, gamma)',"Size"] < unique(tmp[tmp$Method == "Observed", "Size"]) ),
                                                      "Rarefaction", "Extrapolation")
        tmp[tmp == 'Extrap_SC(2n, gamma)'] = ifelse( unique(tmp[tmp$Method == 'Extrap_SC(2n, gamma)',"Size"] < unique(tmp[tmp$Method == "Observed", "Size"]) ),
                                                     "Rarefaction", "Extrapolation")

        tmp[tmp == 'Observed_SC(T, gamma)'] = ifelse( unique(tmp[tmp$Method == 'Observed_SC(T, gamma)',"Size"] < unique(tmp[tmp$Method == "Observed", "Size"]) ),
                                                      "Rarefaction", "Extrapolation")
        tmp[tmp == 'Extrap_SC(2T, gamma)'] = ifelse( unique(tmp[tmp$Method == 'Extrap_SC(2T, gamma)',"Size"] < unique(tmp[tmp$Method == "Observed", "Size"]) ),
                                                     "Rarefaction", "Extrapolation")

        tmp

      })  %>% do.call(rbind,.) %>% mutate(div_type = "Beta")  %>% as_tibble()


      # # Dropping out the points extrapolated over double reference size
      # gamma1 = data.frame() ; alpha1 = data.frame() ; beta1 = data.frame()
      #
      # for(i in 1:length(unique(gamma$Dataset))){
      #
      #   Gamma <- gamma %>% filter(Dataset==unique(gamma$Dataset)[i]) ; ref_size = unique(Gamma[Gamma$Method=="Observed",]$Size)
      #   Gamma = Gamma %>% filter(!(Order.q==0 & round(Size)>2*ref_size))
      #
      #   Alpha <- alpha %>% filter(Dataset==unique(gamma$Dataset)[i]) ; Alpha = Alpha %>% filter(!(Order.q==0 & round(Size)>2*ref_size))
      #   Beta <- beta %>% filter(Dataset==unique(gamma$Dataset)[i]) ; Beta = Beta %>% filter(!(Order.q==0 & round(Size)>2*ref_size))
      #
      #   gamma1 = rbind(gamma1,Gamma) ; alpha1 = rbind(alpha1,Alpha) ; beta1 = rbind(beta1,Beta)
      #
      # }
      #
      # gamma = gamma1 ; alpha = alpha1 ; beta= beta1

      df = rbind(gamma, alpha, beta)
      for (i in unique(gamma$Order.q)) df$Order.q[df$Order.q == i] = paste0('q = ', i)
      df$div_type <- factor(df$div_type, levels = c("Gamma","Alpha","Beta"))

      id_obs = which(df$Method == 'Observed')

      if (length(id_obs) > 0) {

        for (i in 1:length(id_obs)) {

          new = df[id_obs[i],]
          new$SC = new$SC - 0.0001
          new$Method = 'Rarefaction'

          newe = df[id_obs[i],]
          newe$SC = newe$SC + 0.0001
          newe$Method = 'Extrapolation'

          df = rbind(df, new, newe)

        }
      }


    }

    if (type == 'D'){

      C = lapply(output, function(y) {

        tmp = y[["1-C"]]
        if ('mT' %in% colnames(tmp)) tmp = tmp %>% rename("Size" = "mT")

        tmp[tmp == 'Observed_SC(n, alpha)'] = tmp[tmp == 'Observed_SC(T, alpha)'] = 'Observed'
        tmp[tmp == 'Extrap_SC(2n, alpha)'] = tmp[tmp == 'Extrap_SC(2T, alpha)'] = 'Extrapolation'

        tmp[tmp == 'Observed_SC(n, gamma)'] = ifelse( unique(tmp[tmp$Method == 'Observed_SC(n, gamma)',"Size"] < unique(tmp[tmp$Method == "Observed", "Size"]) ),
                                                      "Rarefaction", "Extrapolation")
        tmp[tmp == 'Extrap_SC(2n, gamma)'] = ifelse( unique(tmp[tmp$Method == 'Extrap_SC(2n, gamma)',"Size"] < unique(tmp[tmp$Method == "Observed", "Size"]) ),
                                                     "Rarefaction", "Extrapolation")

        tmp[tmp == 'Observed_SC(T, gamma)'] = ifelse( unique(tmp[tmp$Method == 'Observed_SC(T, gamma)',"Size"] < unique(tmp[tmp$Method == "Observed", "Size"]) ),
                                                      "Rarefaction", "Extrapolation")
        tmp[tmp == 'Extrap_SC(2T, gamma)'] = ifelse( unique(tmp[tmp$Method == 'Extrap_SC(2T, gamma)',"Size"] < unique(tmp[tmp$Method == "Observed", "Size"]) ),
                                                     "Rarefaction", "Extrapolation")

        tmp

      }) %>% do.call(rbind,.) %>% mutate(div_type = "1-CqN") %>% as_tibble()

      U = lapply(output, function(y) {

        tmp = y[["1-U"]]
        if ('mT' %in% colnames(tmp)) tmp = tmp %>% rename("Size" = "mT")

        tmp[tmp == 'Observed_SC(n, alpha)'] = tmp[tmp == 'Observed_SC(T, alpha)'] = 'Observed'
        tmp[tmp == 'Extrap_SC(2n, alpha)'] = tmp[tmp == 'Extrap_SC(2T, alpha)'] = 'Extrapolation'

        tmp[tmp == 'Observed_SC(n, gamma)'] = ifelse( unique(tmp[tmp$Method == 'Observed_SC(n, gamma)',"Size"] < unique(tmp[tmp$Method == "Observed", "Size"]) ),
                                                      "Rarefaction", "Extrapolation")
        tmp[tmp == 'Extrap_SC(2n, gamma)'] = ifelse( unique(tmp[tmp$Method == 'Extrap_SC(2n, gamma)',"Size"] < unique(tmp[tmp$Method == "Observed", "Size"]) ),
                                                     "Rarefaction", "Extrapolation")

        tmp[tmp == 'Observed_SC(T, gamma)'] = ifelse( unique(tmp[tmp$Method == 'Observed_SC(T, gamma)',"Size"] < unique(tmp[tmp$Method == "Observed", "Size"]) ),
                                                      "Rarefaction", "Extrapolation")
        tmp[tmp == 'Extrap_SC(2T, gamma)'] = ifelse( unique(tmp[tmp$Method == 'Extrap_SC(2T, gamma)',"Size"] < unique(tmp[tmp$Method == "Observed", "Size"]) ),
                                                     "Rarefaction", "Extrapolation")

        tmp

      }) %>% do.call(rbind,.)  %>% mutate(div_type = "1-UqN") %>% as_tibble()

      V = lapply(output, function(y) {

        tmp = y[["1-V"]]
        if ('mT' %in% colnames(tmp)) tmp = tmp %>% rename("Size" = "mT")

        tmp[tmp == 'Observed_SC(n, alpha)'] = tmp[tmp == 'Observed_SC(T, alpha)'] = 'Observed'
        tmp[tmp == 'Extrap_SC(2n, alpha)'] = tmp[tmp == 'Extrap_SC(2T, alpha)'] = 'Extrapolation'

        tmp[tmp == 'Observed_SC(n, gamma)'] = ifelse( unique(tmp[tmp$Method == 'Observed_SC(n, gamma)',"Size"] < unique(tmp[tmp$Method == "Observed", "Size"]) ),
                                                      "Rarefaction", "Extrapolation")
        tmp[tmp == 'Extrap_SC(2n, gamma)'] = ifelse( unique(tmp[tmp$Method == 'Extrap_SC(2n, gamma)',"Size"] < unique(tmp[tmp$Method == "Observed", "Size"]) ),
                                                     "Rarefaction", "Extrapolation")

        tmp[tmp == 'Observed_SC(T, gamma)'] = ifelse( unique(tmp[tmp$Method == 'Observed_SC(T, gamma)',"Size"] < unique(tmp[tmp$Method == "Observed", "Size"]) ),
                                                      "Rarefaction", "Extrapolation")
        tmp[tmp == 'Extrap_SC(2T, gamma)'] = ifelse( unique(tmp[tmp$Method == 'Extrap_SC(2T, gamma)',"Size"] < unique(tmp[tmp$Method == "Observed", "Size"]) ),
                                                     "Rarefaction", "Extrapolation")

        tmp

      }) %>% do.call(rbind,.)  %>% mutate(div_type = "1-VqN") %>% as_tibble()

      S = lapply(output, function(y) {

        tmp = y[["1-S"]]
        if ('mT' %in% colnames(tmp)) tmp = tmp %>% rename("Size" = "mT")

        tmp[tmp == 'Observed_SC(n, alpha)'] = tmp[tmp == 'Observed_SC(T, alpha)'] = 'Observed'
        tmp[tmp == 'Extrap_SC(2n, alpha)'] = tmp[tmp == 'Extrap_SC(2T, alpha)'] = 'Extrapolation'

        tmp[tmp == 'Observed_SC(n, gamma)'] = ifelse( unique(tmp[tmp$Method == 'Observed_SC(n, gamma)',"Size"] < unique(tmp[tmp$Method == "Observed", "Size"]) ),
                                                      "Rarefaction", "Extrapolation")
        tmp[tmp == 'Extrap_SC(2n, gamma)'] = ifelse( unique(tmp[tmp$Method == 'Extrap_SC(2n, gamma)',"Size"] < unique(tmp[tmp$Method == "Observed", "Size"]) ),
                                                     "Rarefaction", "Extrapolation")

        tmp[tmp == 'Observed_SC(T, gamma)'] = ifelse( unique(tmp[tmp$Method == 'Observed_SC(T, gamma)',"Size"] < unique(tmp[tmp$Method == "Observed", "Size"]) ),
                                                      "Rarefaction", "Extrapolation")
        tmp[tmp == 'Extrap_SC(2T, gamma)'] = ifelse( unique(tmp[tmp$Method == 'Extrap_SC(2T, gamma)',"Size"] < unique(tmp[tmp$Method == "Observed", "Size"]) ),
                                                     "Rarefaction", "Extrapolation")

        tmp

      }) %>% do.call(rbind,.) %>% mutate(div_type = "1-SqN") %>% as_tibble()

      # C = C %>% filter(Method != 'Observed')
      # U = U %>% filter(Method != 'Observed')
      # V = V %>% filter(Method != 'Observed')
      # S = S %>% filter(Method != 'Observed')


      # # Dropping out the points extrapolated over double reference size
      # c1 = data.frame() ; u1 = data.frame() ; v1 = data.frame() ; s1 = data.frame()
      #
      # for(i in 1:length(unique(C$Dataset))){
      #
      #   CC <- C %>% filter(Dataset==unique(C$Dataset)[i]) ; ref_size = unique(CC[CC$Method=="Observed",]$Size)
      #   CC = CC %>% filter(!(Order.q==0 & round(Size)>2*ref_size))
      #
      #   UU <- U %>% filter(Dataset==unique(C$Dataset)[i]) ; UU = UU %>% filter(!(Order.q==0 & round(Size)>2*ref_size))
      #   VV <- V %>% filter(Dataset==unique(C$Dataset)[i]) ; VV = VV %>% filter(!(Order.q==0 & round(Size)>2*ref_size))
      #   SS <- S %>% filter(Dataset==unique(C$Dataset)[i]) ; SS = SS %>% filter(!(Order.q==0 & round(Size)>2*ref_size))
      #
      #   c1 = rbind(c1,CC) ; u1 = rbind(u1,UU) ; v1 = rbind(v1,VV) ; s1 = rbind(s1,SS)
      #
      # }
      #
      # C = c1 ; U = u1 ; V = v1 ; S = s1


      df = rbind(C, U, V, S)
      for (i in unique(C$Order.q)) df$Order.q[df$Order.q == i] = paste0('q = ', i)
      df$div_type <- factor(df$div_type, levels = c("1-CqN", "1-UqN", "1-VqN", "1-SqN"))

      id_obs = which(df$Method == 'Observed')

      if (length(id_obs) > 0) {

        for (i in 1:length(id_obs)) {

          new = df[id_obs[i],]
          new$SC = new$SC - 0.0001
          new$Method = 'Rarefaction'

          newe = df[id_obs[i],]
          newe$SC = newe$SC + 0.0001
          newe$Method = 'Extrapolation'

          df = rbind(df, new, newe)

        }
      }


    }

    lty = c(Rarefaction = "solid", Extrapolation = "dashed")
    df$Method = factor(df$Method, levels = c('Rarefaction', 'Extrapolation', 'Observed'))

    double_size = unique(df[df$Method == "Observed",]$Size)*2
    double_extrapolation = df %>% filter(Method == "Extrapolation" & round(Size) %in% double_size)

    if ("Pair" %in% colnames(df)) stop("Due to too much data points under 'by_pair = TRUE', 'ggiNEXTbeta3D' don't plot the figure.", call. = FALSE)

    fig = ggplot(data = df, aes(x = SC, y = Estimate, col = Dataset)) +
      geom_line(data = subset(df, Method != 'Observed'), aes(linetype = Method), size=1.1) + scale_linetype_manual(values = lty) +
      # geom_line(lty=2) +
      geom_point(data = subset(df, Method == 'Observed' & div_type == "Gamma"), shape = 19, size = 3) +
      geom_point(data = subset(df, Method == 'Observed' & div_type != "Gamma"), shape = 1, size = 3, stroke = 1.5)+
      geom_point(data = subset(double_extrapolation, div_type == "Gamma"), shape = 17, size = 3) +
      geom_point(data = subset(double_extrapolation, div_type != "Gamma"), shape = 2, size = 3, stroke = 1.5) +
      scale_colour_manual(values = cbPalette) +
      scale_fill_manual(values = cbPalette) +
      facet_grid(div_type ~ Order.q, scales = 'free') +
      theme_bw() +
      theme(legend.position = "bottom",
            legend.title = element_blank(),
            strip.text = element_text(size = 15, face = 'bold'),
            axis.title = element_text(hjust = 0.5, size = 15, face = 'bold'),
            axis.text.x = element_text(size = 12),
            axis.text.y = element_text(size = 12),
            legend.box = "vertical",
            legend.margin = margin(0, 0, 0, 0),
            legend.box.margin = margin(-10, -10, -5, -10),
            legend.text = element_text(size = 13),
            plot.margin = unit(c(5.5, 5.5, 5.5, 5.5), "pt")) +
      labs(x = 'Sample coverage', y = ylab) +
      guides(linetype = guide_legend(keywidth = 2.5))

    if (sum(is.na(df$LCL) == 0)) fig = fig + geom_ribbon(aes(ymin = LCL, ymax = UCL, fill = Dataset, col = NULL), alpha = transp)

  } else if (length(output[[1]]) == 3) {

    gamma = lapply(output, function(y) y[["gamma"]]) %>% do.call(rbind,.) %>% rename("Estimate" = "Gamma") %>% mutate(div_type = "Gamma") %>% as_tibble()
    alpha = lapply(output, function(y) y[["alpha"]]) %>% do.call(rbind,.) %>% rename("Estimate" = "Alpha") %>% mutate(div_type = "Alpha") %>% as_tibble()

    if ('mT' %in% colnames(gamma)) {

      xlab = 'Number of sampling units'
      colnames(gamma)[colnames(gamma) == 'mT'] = 'Size'
      colnames(alpha)[colnames(alpha) == 'mT'] = 'Size'

    } else xlab = 'Number of individuals'

    df = rbind(gamma, alpha)
    for (i in unique(gamma$Order.q)) df$Order.q[df$Order.q == i] = paste0('q = ', i)
    df$div_type <- factor(df$div_type, levels = c("Gamma","Alpha"))

    id_obs = which(df$Method == 'Observed')

    if (sum(id_obs) > 0) {
      for (i in 1:length(id_obs)) {

        new = df[id_obs[i],]
        new$Size = new$Size - 0.0001
        new$Method = 'Rarefaction'

        newe = df[id_obs[i],]
        newe$Size = newe$Size + 0.0001
        newe$Method = 'Extrapolation'

        df = rbind(df, new, newe)

      }
    }

    lty = c(Rarefaction = "solid", Extrapolation = "dashed")
    df$Method = factor(df$Method, levels = c('Rarefaction', 'Extrapolation', 'Observed'))

    double_size = unique(df[df$Method == "Observed",]$Size)*2
    double_extrapolation = df %>% filter(Method == "Extrapolation" & round(Size) %in% double_size)

    if ("Pair" %in% colnames(df)) stop("Due to too much data points under 'by_pair = TRUE', 'ggiNEXTbeta3D' don't plot the figure.", call. = FALSE)

    fig = ggplot(data = df, aes(x = Size, y = Estimate, col = Dataset)) +
      geom_line(data = subset(df, Method != 'Observed'), aes(linetype = Method), size=1.1) + scale_linetype_manual(values = lty) +
      geom_point(data = subset(df, Method == 'Observed' & div_type == "Gamma"), shape = 19, size = 3) +
      geom_point(data = subset(df, Method == 'Observed' & div_type != "Gamma"), shape = 1, size = 3, stroke = 1.5) +
      geom_point(data = subset(double_extrapolation, div_type == "Gamma"), shape = 17, size = 3) +
      geom_point(data = subset(double_extrapolation, div_type != "Gamma"), shape = 2, size = 3, stroke = 1.5) +
      scale_colour_manual(values = cbPalette) +
      scale_fill_manual(values = cbPalette) +
      facet_grid(div_type ~ Order.q, scales = 'free') +
      theme_bw() +
      theme(legend.position = "bottom",
            legend.title = element_blank(),
            strip.text = element_text(size = 15, face = 'bold'),
            axis.title = element_text(hjust = 0.5, size = 15, face = 'bold'),
            axis.text.x = element_text(size = 12),
            axis.text.y = element_text(size = 12),
            legend.box = "vertical",
            legend.margin = margin(0, 0, 0, 0),
            legend.box.margin = margin(-10, -10, -5, -10),
            legend.text = element_text(size = 13),
            plot.margin = unit(c(5.5, 5.5, 5.5, 5.5), "pt")) +
      labs(x = xlab, y = ylab) +
      guides(linetype = guide_legend(keywidth = 2.5))

    if (sum(is.na(df$LCL) == 0)) fig = fig + geom_ribbon(aes(ymin = LCL, ymax = UCL, fill = Dataset, col = NULL), alpha = transp)
  }

  return(fig)
}





#' function to calculate hierarchical phylogenetic gamma, alpha, beta diversity and dissimilarity measure
#'
#' \code{hierPD}: function to calculate empirical estimates for hierarchical phylogenetic gamma, alpha, beta diversity and dissimilarity measure.
#'
#' @param data data should be input as a \code{matrix/data.frame} (species by assemblages).
#' @param mat hierarchical structure of data should be input as a \code{matrix}.
#' @param q a numerical vector specifying the diversity orders. Default is \code{seq(0, 2, 0.2)}.
#' @param weight (required only when \code{type = "mle"} and \code{decomposition = "relative"}) weight for relative decomposition empirical estimate. Select size-weighted \code{("size")}, equal-weighted \code{("equal")} or a numerical vector for weight. Default is \code{"size"}.
#' @param nboot a positive integer specifying the number of bootstrap replications when assessing sampling uncertainty and constructing confidence intervals. Bootstrap replications are generally time consuming. Enter \code{0} to skip the bootstrap procedures. Default is \code{10}.
#' @param conf a positive number < 1 specifying the level of confidence interval. Default is \code{0.95}.
#' @param PDtree a phylogenetic tree in Newick format for all observed species in the pooled assemblage.
#' @param type estimate type: empirical \code{(type = "mle")} or asymptotic estimate \code{(type = "est")}. Default is \code{"mle"}.
#' @param decomposition decomposition type: relative \code{(decomposition = "relative")} or absolute decomposition \code{(decomposition = "absolute")}. Default is \code{"relative"}.
#'
#' @return a data frames with hierarchical phylogenetic diversity (gamma, alpha, and beta) and four types dissimilarity measure.
#' @importFrom ape reorder.phylo
#' @importFrom ade4 newick2phylog
#' @examples
#'
#' data("global")
#' data("global_tree")
#' data("global_mat")
#' hier_output <- hierPD(global, mat = global_mat, q = seq(0, 2, 0.2), PDtree = global_tree)
#'
#' @export
hierPD <- function(data, mat, q = seq(0, 2, 0.2), weight = "size", nboot = 10, conf = 0.95,
                   PDtree, type = "mle", decomposition = "relative"){

  hier_method = c("qPD", "1-C", "1-U", "1-V", "1-S")
  out = hier.phylogeny(data, mat, tree = PDtree, q = q, weight = weight, nboot = nboot,
                       conf = conf, type = type, decomposition = decomposition)
  out = out[str_sub(out$Method, 1, 3) %in% hier_method, ]

  if (decomposition == "relative"){
    diss_method = out$Method[str_sub(out$Method, 1, 3) %in% hier_method[-1]]
    out$Method[str_sub(out$Method, 1, 3) %in% hier_method[-1]] = paste0(str_sub(diss_method, 1, 5), "*", str_sub(diss_method, 6))
  }

  out
}


#' ggplot2 extension for an hierPD object
#'
#' \code{gghierPD}: the \code{ggplot} extension for \code{\link{hierPD}} object to plot order q against to hierarchical phylogenetic diversity decomposition and dissimilarity measure.
#'
#' @param output the output from hierPD.
#' @param type selection of plot type : \cr
#' \code{(type = "A")} for alpha and gamma diversity; \cr
#' \code{(type = "B")} for beta diversity; \cr
#' \code{(type = "D")} for dissimilarity measure based on multiplicative decomposition.
#'
#' @return a figure for hierarchical phylogenetic diversity decomposition or dissimilarity measure.
#'
#' @examples
#'
#' data("global")
#' data("global_mat")
#' data("global_tree")
#' hier_output <- hierPD(global, mat = global_mat, q = seq(0, 2, 0.2), PDtree = global_tree)
#' gghierPD(hier_output, type = "A")
#'
#' @export
gghierPD <- function(output, type = "A"){

  m = ifelse(type=="A", 4,
             ifelse(type=="B", 5,
                    ifelse(type=="D", 6, NA)))

  if (m == 6) {
    output = output[grep("1-", output$Method), ]
    if (unique(output$Decomposition) == "relative") {
      output$group = "1-CqN*"
      output$group[grep("1-UqN*", output$Method)] = "1-UqN*"
      output$group[grep("1-SqN*", output$Method)] = "1-SqN*"
      output$group[grep("1-VqN*", output$Method)] = "1-VqN*"
    }
    else if (unique(output$Decomposition) == "absolute") {
      output$group = "1-CqN"
      output$group[grep("1-UqN", output$Method)] = "1-UqN"
      output$group[grep("1-SqN", output$Method)] = "1-SqN"
      output$group[grep("1-VqN", output$Method)] = "1-VqN"
    }
    output$Method = fct_inorder(output$Method)
    plot_out = ggplot(output, aes(x = Order.q, y = Estimator,
                                  colour = Method, fill = Method)) + facet_grid(fct_inorder(group) ~ .) +
      geom_line(size = 1.5) +
      geom_ribbon(aes(ymin = LCL, ymax = UCL, fill = Method), linetype = 0, alpha = 0.2) +
      theme_bw() +
      theme(legend.position = "bottom", legend.box = "vertical", legend.key.width = unit(1.2, "cm"),
            legend.title = element_blank(), legend.margin = margin(0, 0, 0, 0),
            legend.box.margin = margin(-10, -10, -5, -10), text = element_text(size = 16),
            plot.margin = unit(c(5.5, 5.5, 5.5, 5.5), "pt")) +
      guides(linetype = guide_legend(keywidth = 2.5))
  }
  else {
    plot_out = gghier_phylogeny(output, method = m)
  }
  plot_out + xlab("Order q") + ylab("Estimate") +
    theme(strip.text = element_text(size = 15, face = "bold"),
          axis.title = element_text(hjust = 0.5, size = 15, face = "bold"),
          axis.text.x = element_text(size = 12),
          axis.text.y = element_text(size = 12),
          legend.text = element_text(size = 13))
}

## ========== no visible global function definition for R CMD check ========== ##
utils::globalVariables(c(
  ".", "Dataset", "Estimate", "Estimator", "LCL", "Method", "Order.q", "SC", "Size", "UCL",
  "branch.abun", "branch.length", "div_type", "ggplotColors", "label",
  "qPD", "tgroup", "sd"
))

