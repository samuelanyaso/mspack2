##################################################################### estimates of occupation probabilities at specific
##################################################################### time points
getProbs <- function(fit, cutoffs) {

    # time grid
    et <- fit[, c("time")]

    # indices for the specified cutoffs
    indx <- sapply(cutoffs, function(x) max(which(et <= x)), simplify = TRUE)
    fit <- fit[indx, ]

    return(fit)
}


##################################################################### Estimate state occupation probability for
##################################################################### unclustered data
curr_stat_unclust <- function(dat, tree, start.probs, bw = NULL, ngrid = NULL, pavY = TRUE, cond = NULL, Zval = NULL) {

    colnames(dat)[which(colnames(dat) == "stage")] <- "state"

    vars <- c("id", "time", "state")
    varchk <- match(vars, colnames(dat), nomatch = -1)
    if (any(varchk < 1)) {
        stop("colnames of dat should include 'id', 'time', 'state'")
    }

    # functions for kernel smoothing
    Kh <- function(C, t, h) {
        h^(-1) * dnorm((C - t)/h)
    }

    # data with sorted inspection time
    df <- with(dat, dat[order(dat$time), ])
    cstimes <- df$time
    id <- df$id
    nind <- length(unique(id))  # total no. of subjects
    n <- length(unique(id))  # total no. of clusters

    # create time grid
    maxt <- max(c(cstimes))
    if (is.null(ngrid)) {
        ngrid <- 1000
    }
    timegrid <- c(0, seq(from = (0 + 1e-04), to = maxt, length.out = ngrid))
    timegrid <- sort(timegrid)
    timegrid <- round(timegrid, 6)

    # calculate the bandwidth for the covariate
    if (is.null(bw)) {
        bw <- dpik(cstimes, kernel = "normal", gridsize = length(timegrid), range.x = range(cstimes))
    }

    # calculate the weights
    wICS <- rep(1, n)
    if (!is.null(cond)) {
        Zcov <- df[, as.character(cond)]

        if (is.null(Zval)) {
            Zval <- as.numeric(quantile(Zcov, 0.25))
        }

        # calculate the bandwidth for the covariate
        h_Z <- dpik(Zcov, gridsize = length(cstimes))
        # calculate the kernel density based on the covariates
        wZ <- as.numeric(apply(as.matrix(Zcov), 2, Kh, t = Zval, h = h_Z))

    } else {
        wZ <- rep(1, n)
    }
    # if the computed weights are 0, assign the min weights to such weights.
    if (sum(wZ == 0) > 0)
        wZ[wZ == 0] <- min(wZ[which(wZ != 0)])


    # possible transitions
    nt.states <- which(sapply(edgeL(tree), function(x) length(x$edges) > 0))  # gets the indices for the non-absorbing state
    lng <- sapply(edges(tree)[nodes(tree) %in% names(nt.states)], length)  # gets the no. of outgoing edges for each non-absorbing state
    trans <- paste(rep(nodes(tree)[nodes(tree) %in% names(nt.states)], lng), unlist(edges(tree)[nodes(tree) %in% names(nt.states)]),
        sep = "")  # names the possible transitions

    ## Indicators I(Ujj' <= C) -- this is the indicator that the j to j' transition has taken place by time C.
    Is <- matrix(0, nrow = length(cstimes), ncol = length(trans))
    rownames(Is) <- cstimes
    colnames(Is) <- paste("I", trans, sep = "")

    for (i in nodes(tree)) {
        if (length(inEdges(tree)[[i]]) == 0)
            next  # skip this iteration, if it is in the initial state
        ld <- inEdges(tree)[[i]]  #nodes from
        ex <- edges(tree)[[i]]  #nodes to
        later.states <- names(acc(tree, i)[[1]])  # gets the first node accessible from the current node
        states <- c(i, later.states)  # vector of current state and later state
        b <- paste("I", inEdges(tree)[[i]], i, sep = "")  # create name for incoming transition to the current node
        row.idx <- which(df$state %in% states)  # get indices for subjects who when inspected were at the current state or the later states, thus, they had alread made the bth transition
        col.idx <- which(colnames(Is) %in% b)
        Is[row.idx, col.idx] <- 1
    }

    # storage for initial estimates, row names are the inspection times
    Nstar <- matrix(0, nrow = length(cstimes), ncol = length(trans))
    colnames(Nstar) <- colnames(Is)
    # rownames(Nstar) <- rownames(Is)

    # performs isotonic regression
    Nstar[, colnames(Is)] <- apply(Is, 2, function(y) gpava(z = cstimes, y = y, weights = (wICS * wZ))$x)
    # Nstar <- apply(Nstar, 2, function(x) round(x*nind)) # Datta and Sundaram rounded these figures

    # (kernel estimate) of the functions for estimating N and Y
    tmp <- outer(timegrid, cstimes, Kh, h = bw)

    # perform kernel smoothing -- counting process
    Ns <- apply(Nstar, 2, function(pav) {
        gc <- (tmp %*% wICS)/nind
        (tmp %*% (pav * wICS))/gc
    })
    rownames(Ns) <- timegrid

    ## Computing the at-risk set
    Ys <- matrix(0, nrow = length(timegrid), ncol = length(nt.states))
    colnames(Ys) <- paste("Y", nodes(tree)[nt.states], sep = "")

    # compute state indicator, I(Si(Ci) = j)
    nstate <- length(unique(nodes(tree)))  # number of states in the multistate model
    Ij <- t(sapply(df$state, function(x) {
        a <- rep(0, nstate)
        a[x] = 1
        return(a)
    }, simplify = TRUE))
    colnames(Ij) <- paste("Y", nodes(tree), sep = "")
    rownames(Ij) <- rownames(Is)
    Ij <- Ij[, colnames(Ys)]

    # computes the at-risk set via kernel smoothing include a PAV step for the transitions out of the initial state.
    if (pavY) {
        Ij <- apply(Ij, 2, function(y) (-gpava(z = cstimes, y = -y, weights = (wICS * wZ))$x))
        Ys <- apply(Ij, 2, function(y) {
            gc <- (tmp %*% (wICS * wZ))/nind
            (tmp %*% (y * (wICS * wZ)))/gc
        })
    } else {
        Ys <- apply(Ij, 2, function(y) {
            gc <- (tmp %*% (wICS * wZ))/nind
            (tmp %*% (y * (wICS * wZ)))/gc
        })
    }

    colnames(Ys) <- paste("Y", nodes(tree)[nt.states], sep = "")
    rownames(Ys) <- timegrid
    rm(tmp)

    # counts the number of actual transitions as defined in the column
    dNs <- apply(Ns, 2, function(x) diff(x))

    # adjust the at risk set
    Ys <- Ys[-1, ]

    # assign column names to dNs replicates the non-absorbing states by the number of edges out of each state gets the
    # edges proceding out of each non-absorbing state
    ds <- paste("dN", rep(nodes(tree)[nodes(tree) %in% names(nt.states)], lng), unlist(edges(tree)[nodes(tree) %in% names(nt.states)]))
    colnames(dNs) <- ds

    # remove rows with no transitions
    rm.idx <- which(rowSums(dNs) <= 0)
    if (length(rm.idx) != 0) {
        dNs <- dNs[-rm.idx, ]  # remove rows with no transitions
        Ys <- Ys[-rm.idx, ]
    }

    ## compute state occupation probability
    cum.tm <- diag(nstate)
    colnames(cum.tm) <- rownames(cum.tm) <- nodes(tree)

    ps <- matrix(NA, nrow = nrow(dNs), ncol = length(nodes(tree)))  # to store the state occupation prob. for each state at each event time
    rownames(ps) <- rownames(dNs)
    colnames(ps) <- paste("p", nodes(tree), sep = "")

    # creates an array to store teh matrices of the estimates at each event time
    all.dA <- all.I_dA <- all.ajs <- array(dim = c(nstate, nstate, nrow(dNs)), dimnames = list(rows = nodes(tree), cols = nodes(tree),
        dim = rownames(dNs)))

    for (i in 1:nrow(dNs)) {
        I_dA <- diag(nstate)
        dA <- matrix(0, nrow = nstate, ncol = nstate)
        colnames(I_dA) <- rownames(I_dA) <- colnames(dA) <- rownames(dA) <- nodes(tree)

        idx <- which(dNs[i, ] > 0)
        t.nam <- colnames(dNs)[idx]
        tmp <- strsplit(t.nam, " ")
        start <- sapply(tmp, function(x) x[2])
        end <- sapply(tmp, function(x) x[3])
        idxs <- matrix(as.numeric(c(start, end)), ncol = 2)
        idxs2 <- matrix(as.numeric(c(start, start)), ncol = 2)

        # computes the integrated transition hazard matrix
        tmp <- dNs[i, idx]/Ys[i, paste("Y", start, sep = "")]
        tmp[which(tmp > 1)] <- 1
        dA[idxs] <- tmp

        if (length(idx) == 1) {
            tmp <- -dNs[i, idx]/Ys[i, paste("Y", start, sep = "")]
            tmp[which(abs(tmp) > 1)] <- -1
            dA[start, start] <- tmp
        } else {
            dA[idxs2] <- -rowSums(dA[start, ])
        }
        rm(tmp)

        I_dA <- I_dA + dA  # I + dA matrix
        all.dA[, , i] <- dA  # stores the dA matrix at time i
        all.I_dA[, , i] <- I_dA  # stores the I_dA matrix at ime

        cum.tm <- cum.tm %*% I_dA  # the product integral; the transition probability matrix at time i
        all.ajs[, , i] <- cum.tm  # stores the transition probability matrix at time i

        ps[i, ] <- start.probs %*% all.ajs[, , i]  # the state occupation prob at time i
    }
    rm.idx <- apply((is.na(ps) | is.infinite(ps)), 1, any)
    ps <- ps[!rm.idx, ]
    rm.idx <- apply((ps < 0 | ps > 1), 1, any)
    ps <- ps[!rm.idx, ]
    return(ps)
}


##################################################################### Estimate state occupation probability for
##################################################################### cluster-correlated data
curr_stat <- function(dat, tree, start.probs, weight = "ICS", bw = NULL, ngrid = NULL, pavY = TRUE, cond = NULL, Zval = NULL) {

    colnames(dat)[which(colnames(dat) == "stage")] <- "state"

    vars <- c("cID", "id", "time", "state", "csize")
    varchk <- match(vars, colnames(dat), nomatch = -1)
    if (any(varchk < 1)) {
        stop("colnames of dat should include 'cID', 'id', 'time', 'state', 'csize'")
    }

    if (any(unique(dat$state) == 0)) {
        stop("Do not name any state as '0', the initial node should be '1' ")
    }

    if (!all(as.character(unique(dat$state)) %in% nodes(tree))) {
        stop("states in 'dat' object do not match with states in the 'tree' object ")
    }

    # functions for kernel smoothing
    Kh <- function(C, t, h) {
        h^(-1) * dnorm((C - t)/h)
    }

    # data with sorted inspection times
    df <- with(dat, dat[order(dat$time), ])
    cstimes <- df$time
    # cstimes <- cstimes - min(cstimes)
    nind <- length(unique(df$id))  # total no. of subjects
    m <- length(unique(df$cID))  # total no. of clusters
    csize <- df$csize
    cid <- df$cID

    # create time grid
    maxt <- max(c(cstimes))
    if (is.null(ngrid)) {
        ngrid <- 1000
    }
    timegrid <- c(0, seq(from = (0 + 1e-04), to = maxt, length.out = ngrid))
    timegrid <- sort(timegrid)
    timegrid <- round(timegrid, 6)

    # calculate the bandwidth for the covariate
    if (is.null(bw)) {
        bw <- dpik(cstimes, kernel = "normal", gridsize = length(timegrid), range.x = range(cstimes))
    }
    # bw <- 0.5 # for this example, this bandwidth was chosen in the paper.

    # calculate the weights
    if (!is.null(cond)) {
        Zcov <- df[, as.character(cond)]

        if (is.null(Zval)) {
            Zval <- as.numeric(quantile(Zcov, 0.25))
        }

        # calculate the bandwidth for the covariate
        tmpZ <- rep(0, length(unique(cid)))
        for (i in 1:length(unique(cid))) {
            tmpZ[i] <- sample(Zcov[cid == unique(cid)[i]], 1)
        }
        h_Z <- dpik(tmpZ, gridsize = length(cstimes))

        if (weight == "ICS") {
            wICS <- 1/csize
            wZ <- as.numeric(apply(as.matrix(Zcov), 2, Kh, t = Zval, h = h_Z))
        } else if (weight == 1) {
            wICS <- rep(1, length(csize))
            wZ <- as.numeric(apply(as.matrix(Zcov), 2, Kh, t = Zval, h = h_Z))
        }
    } else {
        if (weight == "ICS") {
            wICS <- 1/csize
            wZ <- rep(1, length(csize))
        } else if (weight == 1) {
            wICS <- rep(1, length(csize))
            wZ <- rep(1, length(csize))
        }
    }
    # if the computed weights are 0, assign the min weights to such weights.
    if (sum(wZ == 0) > 0)
        wZ[wZ == 0] <- min(wZ[which(wZ != 0)])


    # possible transitions
    nt.states <- which(sapply(edgeL(tree), function(x) length(x$edges) > 0))  # gets the indices for the non-absorbing state
    lng <- sapply(edges(tree)[nodes(tree) %in% names(nt.states)], length)  # gets the no. of outgoing edges for each non-absorbing state
    trans <- paste(rep(nodes(tree)[nodes(tree) %in% names(nt.states)], lng), unlist(edges(tree)[nodes(tree) %in% names(nt.states)]),
        sep = "")  # names the possible transitions

    ## Indicators I(Ujj' <= C) -- this is the indicator that the j to j' transition has taken place by time C.
    Is <- matrix(0, nrow = length(cstimes), ncol = length(trans))
    rownames(Is) <- cstimes
    colnames(Is) <- paste("I", trans, sep = "")

    for (i in nodes(tree)) {
        if (length(inEdges(tree)[[i]]) == 0)
            next  # skip this iteration, if it is in the initial state
        ld <- inEdges(tree)[[i]]  #nodes from
        ex <- edges(tree)[[i]]  #nodes to
        later.states <- names(acc(tree, i)[[1]])  # gets the first node accessible from the current node
        states <- c(i, later.states)  # vector of current state and later state
        b <- paste("I", inEdges(tree)[[i]], i, sep = "")  # create name for incoming transition to the current node
        row.idx <- which(df$state %in% states)  # get indices for subjects who when inspected were at the current state or the later states, thus, they had alread made the bth transition
        col.idx <- which(colnames(Is) %in% b)
        Is[row.idx, col.idx] <- 1
    }

    # storage for initial estimates, row names are the inspection times
    Nstar <- matrix(0, nrow = length(cstimes), ncol = length(trans))
    colnames(Nstar) <- colnames(Is)
    # rownames(Nstar) <- rownames(Is)

    # performs isotonic regression
    Nstar[, colnames(Is)] <- apply(Is, 2, function(y) gpava(z = cstimes, y = y, weights = (wICS * wZ))$x)
    # Nstar <- apply(Nstar, 2, function(x) round(x*nind)) # Datta and Sundaram rounded these figures

    # (kernel estimate) of the functions for estimating N and Y
    tmp <- outer(timegrid, cstimes, Kh, h = bw)

    # perform kernel smoothing -- counting process
    Ns <- apply(Nstar, 2, function(pav) {
        gc <- (tmp %*% wICS)/nind
        (tmp %*% (pav * wICS))/gc
    })
    rownames(Ns) <- timegrid

    ## Computing the at-risk set
    Ys <- matrix(0, nrow = length(timegrid), ncol = length(nt.states))
    colnames(Ys) <- paste("Y", nodes(tree)[nt.states], sep = "")

    # compute state indicator, I(Si(Ci) = j)
    nstate <- length(unique(nodes(tree)))  # number of states in the multistate model
    Ij <- t(sapply(df$state, function(x) {
        a <- rep(0, nstate)
        a[x] = 1
        return(a)
    }, simplify = TRUE))
    colnames(Ij) <- paste("Y", nodes(tree), sep = "")
    rownames(Ij) <- rownames(Is)
    Ij <- Ij[, colnames(Ys)]

    # computes the at-risk set via kernel smoothing include a PAV step for the transitions out of the initial state.
    if (pavY) {
        Ij <- apply(Ij, 2, function(y) (-gpava(z = cstimes, y = -y, weights = (wICS * wZ))$x))
        Ys <- apply(Ij, 2, function(y) {
            gc <- (tmp %*% (wICS * wZ))/nind
            (tmp %*% (y * (wICS * wZ)))/gc
        })
    } else {
        Ys <- apply(Ij, 2, function(y) {
            gc <- (tmp %*% (wICS * wZ))/nind
            (tmp %*% (y * (wICS * wZ)))/gc
        })
    }

    colnames(Ys) <- paste("Y", nodes(tree)[nt.states], sep = "")
    rownames(Ys) <- timegrid
    rm(tmp)

    # counts the number of actual transitions as defined in the column
    dNs <- apply(Ns, 2, function(x) diff(x))

    # adjust the at risk set
    Ys <- Ys[-1, ]

    # assign column names to dNs replicates the non-absorbing states by the number of edges out of each state gets the
    # edges proceding out of each non-absorbing state
    ds <- paste("dN", rep(nodes(tree)[nodes(tree) %in% names(nt.states)], lng), unlist(edges(tree)[nodes(tree) %in% names(nt.states)]))
    colnames(dNs) <- ds

    # remove rows with no transitions
    rm.idx <- which(rowSums(dNs) <= 0)
    if (length(rm.idx) != 0) {
        dNs <- dNs[-rm.idx, ]  # remove rows with no transitions
        Ys <- Ys[-rm.idx, ]
    }

    ## compute state occupation probability
    cum.tm <- diag(nstate)
    colnames(cum.tm) <- rownames(cum.tm) <- nodes(tree)

    ps <- matrix(NA, nrow = nrow(dNs), ncol = length(nodes(tree)))  # to store the state occupation prob. for each state at each event time
    rownames(ps) <- rownames(dNs)
    colnames(ps) <- paste("p", nodes(tree), sep = "")

    # creates an array to store teh matrices of the estimates at each event time
    all.dA <- all.I_dA <- all.ajs <- array(dim = c(nstate, nstate, nrow(dNs)), dimnames = list(rows = nodes(tree), cols = nodes(tree),
        dim = rownames(dNs)))

    for (i in 1:nrow(dNs)) {
        I_dA <- diag(nstate)
        dA <- matrix(0, nrow = nstate, ncol = nstate)
        colnames(I_dA) <- rownames(I_dA) <- colnames(dA) <- rownames(dA) <- nodes(tree)

        idx <- which(dNs[i, ] > 0)
        t.nam <- colnames(dNs)[idx]
        tmp <- strsplit(t.nam, " ")
        start <- sapply(tmp, function(x) x[2])
        end <- sapply(tmp, function(x) x[3])
        idxs <- matrix(as.numeric(c(start, end)), ncol = 2)
        idxs2 <- matrix(as.numeric(c(start, start)), ncol = 2)

        # computes the integrated transition hazard matrix
        tmp <- dNs[i, idx]/Ys[i, paste("Y", start, sep = "")]
        tmp[which(tmp > 1)] <- 1
        dA[idxs] <- tmp

        if (length(idx) == 1) {
            tmp <- -dNs[i, idx]/Ys[i, paste("Y", start, sep = "")]
            tmp[which(abs(tmp) > 1)] <- -1
            dA[start, start] <- tmp
        } else {
            dA[idxs2] <- -rowSums(dA[start, ])
        }
        rm(tmp)

        I_dA <- I_dA + dA  # I + dA matrix
        all.dA[, , i] <- dA  # stores the dA matrix at time i
        all.I_dA[, , i] <- I_dA  # stores the I_dA matrix at ime

        cum.tm <- cum.tm %*% I_dA  # the product integral; the transition probability matrix at time i
        all.ajs[, , i] <- cum.tm  # stores the transition probability matrix at time i

        ps[i, ] <- start.probs %*% all.ajs[, , i]  # the state occupation prob at time i
    }
    rm.idx <- apply((is.na(ps) | is.infinite(ps)), 1, any)
    ps <- ps[!rm.idx, ]
    rm.idx <- apply((ps < 0 | ps > 1), 1, any)
    ps <- ps[!rm.idx, ]
    return(ps)
}

