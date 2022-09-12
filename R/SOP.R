##' Nonparametric estimators for Current-status data from multistate models
##'
##' Utilizes theory on nonparametric regression to estimate the counting process and at-risk sets for current status data from
##' general multistate models. The counting process and at-risk sets are necessary functionals for estimating the state occupation
##' probabilities.
##' @title Nonparametric estimators for Current-status data from multistate models
##' @param data.type a character string specifying the structure of the observed data. This can either be 'uncorrelated' or 'cluster-correlated'.
##' @param dat a data frame with columns 'id', 'time', 'state' indicating unique subject ID, inspection time and state occupation at inspection time. If data.type='cluster-correlated',
##' the data frame should also include 'cID' and 'csize' which are variables that specifies the unique cluster ID and cluster size.
##' @param tree a graphNEL graph with the nodes and edges of the multistate model.
##' @param start.probs a named numeric vector specifying the initial occupation probabilities of each state. The names of the elements must correspond to the nodes of tree.
##' @param weight a character vector to specify the weight for each subject if data.type='cluster-correlated'. Weights can be either 'ICS' or 'none'.
##' @param bw a numeric value specifying a bandwidth sequence for the inspection times. Default is NULL, where the bandwidth sequence is chosen by the criteria of \insertCite{wand1994;textual}{mspack2}.
##' @param ngrid an integer specifying specifying the number of grid points from 0 to max of the inspection times. Default is 1000.
##' @param pavY a logical argument specifying whether the pooled adjascent violator algorithms should be used to estimate the at-risk set. Default is TRUE.
##' @param cond a character argument specifying a single continuous covariate for conditional estimators. This is based on the method by \insertCite{lan2017;textual}{mspack2}. Default is 'none' which indicates marginal estimators are of interest.
##' @param Zval a numeric value indicating the value of the covariate at which the conditional estimator is evaluated at. Default is NULL where the conditional estimator is evaluated at the first quartile of the single covariate.
##' @param cutoffs a numeric vector specifying the time points at which the estimates are returned. Default is NULL, which returns the SOPs at all inspection time points.
##' @return a data frame of state occupation probabilities of each state in the multistate model.
##' @author Samuel Anyaso-Samuel, Dipankar Bandyopadhyay, Somnath.Datta
##' @import graph
##' @import Rgraphviz
##' @import isotone
##' @import KernSmooth
##' @import stats
##' @export
##' @references{
##'   \insertRef{datta2006}{mspack2}
##' }
##' @importFrom Rdpack reprompt
##' @examples
##' ## Three-state tracking model
##' Nodes <- c('1','2','3')
##' Edges <- list('1'=list(edges=c('2')), '2'=list(edges=c('3')),'3'=list(edges=NULL))
##'
##' ## Constructs the tree
##' tree <- new('graphNEL', nodes=Nodes, edgeL=Edges, edgemode='directed')
##'
##' ## simple data simulation
##' m <- 5
##' cdat <- lapply(1:m, function(x){
##'   ni <- rpois(n=1, lambda=10)+2 # cluster size
##'   cID <- rep(x=x, times=ni) # cluster ids
##'   time <- rweibull(n=ni, shape=3, scale=5) # inspection time
##'   state <- sample(x=c(as.numeric(Nodes)), size=ni, replace=TRUE) # state at inspection time
##'   csize <- rep(x=ni, times=ni) # cluster sizes
##'   data.frame(cID=cID, time=time, state=state, csize=csize)
##' })
##' cdat <- do.call(rbind, cdat)
##' cdat$id <- 1:nrow(cdat)
##'
##' ## initial probabilities - under the assumption that everyone states in state 1
##' start.probs <- c(1, 0, 0)  # assumes all subjects start from state 1 and time 0
##' names(start.probs) <- Nodes
##'
##' ## Estimate the state occupation probabilities
##' SOP(data.type='cluster-correlated',dat=cdat, tree=tree, start.probs=start.probs,
##'     ngrid=1000, weight='none', pavY=TRUE, cond='none', Zval=NULL,
##'     cutoffs = seq(from=min(cdat$time), to=max(cdat$time), length.out=10))
SOP <- function(data.type, dat, tree, start.probs, weight = "none", bw = NULL, ngrid = 1000, pavY = TRUE, cond = "none", Zval = NULL, cutoffs = NULL) {

    # figure out the data type
    arg0 <- c("uncorrelated", "cluster-correlated")
    arg0chk <- charmatch(data.type, arg0)
    if (is.na(arg0chk)) {
        stop("data.type should be either 'uncorrelated' or 'cluster-correlated'")
    }

    # figure out the weight
    arg1 <- c("ICS", "none")
    arg1chk <- charmatch(weight, arg1)
    if (is.na(arg1chk)) {
        stop("weight should be either 'ICS' or 'none'")
    }

    if (any(unique(dat$state) == 0)) {
        stop("Do not name any state as '0', the initial node should be '1' ")
    }

    if (!all(as.character(unique(dat$state)) %in% nodes(tree))) {
        stop("states in 'dat' object do not match with states in the 'tree' object ")
    }

    if (weight == "none")
        weight <- 1
    if (cond == "none")
        cond <- NULL

    if (data.type == "uncorrelated") {
        res <- curr_stat_unclust(dat = dat, tree = tree, start.probs = start.probs, bw = bw, ngrid = ngrid, pavY = pavY, cond = cond, Zval = Zval)
    } else if (data.type == "cluster-correlated") {
        res <- curr_stat(dat = dat, tree = tree, start.probs = start.probs, weight = weight, bw = bw, ngrid = ngrid, pavY = pavY, cond = cond, Zval = Zval)
    }

    # get the estimates at specific time points
    res <- data.frame(time = as.numeric(row.names(res)), res)
    row.names(res) <- NULL

    if (is.null(cutoffs)) {
        return(res)
    } else {
        res <- getProbs(fit = res, cutoffs = cutoffs)
        row.names(res) <- NULL
        return(res)
    }
}

