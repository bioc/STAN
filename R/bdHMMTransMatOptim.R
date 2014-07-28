
#' Internal function that evaluates constraints of bidirectional transition matrix.
#' @param transMat Transition matrix.
#' @param xsi_sum Sum over all Xi(t,i,j).
#' @param constraints Constraints of a bidirectional transition matrix as a list.
#' @param Stationary distribution of the transition matrix.
#' 
#' @keywords internal
#' @noRd
eval_constr <- function(transMat, xsi_sum, constraints, statDistr) {
    # , statDistr
    z = c()
    
    
    ## symetric probs
    for (k in 1:length(constraints)) {
        first_prob = transMat[constraints[[k]][[1]][1], constraints[[k]][[1]][2]]
        
        sec_prob = transMat[constraints[[k]][[2]][1], constraints[[k]][[2]][2]]
        
        mult = statDistr[constraints[[k]][[2]][1]]/statDistr[constraints[[k]][[2]][2]]
        z[length(z) + 1] = first_prob - sec_prob * mult
    }
    
    return(z)
}




#'  
#' The function is called from C++ to optimize transitions.
#' 
#' @title Optimize transitions
#' 
#' @param pars Parameters for optimization.
#' 
#' @return optimized transitions
#' @usage c2optimize(pars)
#' 
#' @keywords internal
#' @noRd
c2optimize = function(pars) {
    out = NULL
    
    if (pars$method == "rsolnp") {
        out = c2solnp2(pars)
    } else {
        stop("Must specify method for numerical optimization of transition matrix!\n")
    }
    out
}


#' Internal function used to calculate the full transition matrix from the parameter vector.
#' @param params Transitions and initial state probabilities (parameters for opimization ).
#' @param xsi_sum Sum over all Xi(t,i,j).
#' @param constraints Constraints of a bidirectional transition matrix as a list.
#' @param nStates Number of states of the HMM.
#' @param statD Stationary distribution of transitions.
#' 
#' @keywords internal
#' @noRd
bidirGetTransMat_fixed2 = function(params, xsi_sum, constraints, nStates, statD) {
    
    newMat = matrix(0, nrow = nStates, ncol = nStates)
    for (i in 1:length(constraints)) {
        newMat[constraints[[i]][[1]][1], constraints[[i]][[1]][2]] = params[i]
        newMat[constraints[[i]][[2]][1], constraints[[i]][[2]][2]] = newMat[constraints[[i]][[1]][1], 
            constraints[[i]][[1]][2]] * statD[constraints[[i]][[1]][1]]/statD[constraints[[i]][[1]][2]]
    }
    
    newMat
}




#' Internal function that checks Rsolnp constraints during optimization.
#' @param params Transitions and initial state probabilities (parameters for opimization ).
#' @param xsi_sum Sum over all Xi(t,i,j).
#' @param constraints Constraints of a bidirectional transition matrix as a list.
#' @param nStates Number of states of the HMM.
#' @param stateLabel  Indicates directinality of states. States can be forward (F1, F2, ..., Fn), reverse (R1, R2, ..., Rn) or undirectional (U1, U2, ..., Um). Number of F and R states must be equal and twin states are indicated by integers in id (e.g. F1 and R1 and twins).
#' @param initGamma Gamma(0,i).
#' 
#' @keywords internal
#' @noRd
bidirEqB2 = function(params, xsi_sum, constraints, nStates, stateLabel, initGamma) {
    stateLabel = stateLabel[-grep("R", stateLabel)]
    myF = grep("F", stateLabel)
    U = grep("U", stateLabel)
    statD = c(params[length(constraints) + myF], params[length(constraints) + myF], 
        params[length(constraints) + U])
    
    newMat = bidirGetTransMat_fixed2(params, xsi_sum, constraints, nStates, statD)
    # newMat[which(newMat < 1e-100)] = 1e-100
    return(1 - c(apply(newMat, 1, sum), sum(statD)))
}


#'  
#' The function calls R densitiy function of inverse Wishart (Prior of Multivariate Gaussian) from C++.
#' 
#' @title Fit a Hidden Markov Model
#' 
#' @param hyperparams Hyperparameters for Iverse Wishart as list(cov=Covariance-matrix, S=Scale-Matrix, v=Degrees_of_freedom).
#' 
#' @return Log density of the Inverse-Wishart distribution.
#' @usage calldiwish(hyperparams)
#' 
#' @keywords internal
#' @noRd
calldiwish = function(hyperparams) {
    hyperparams$cov = matrix(hyperparams$cov, nrow = dim(hyperparams$S)[1], ncol = dim(hyperparams$S)[2])
    val = logdiwish(hyperparams$cov, hyperparams$v, hyperparams$S)
    out = val
    return(out)
}


#' Evaluates Q(theta, theta^old) during optimization of transitions.
#' @param params Transitions and initial state probabilities (parameters for opimization ).
#' @param xsi_sum Sum over all Xi(t,i,j).
#' @param constraints Constraints of a bidirectional transition matrix as a list.
#' @param nStates Number of states of the HMM.
#' @param stateLabel  Indicates directinality of states. States can be forward (F1, F2, ..., Fn), reverse (R1, R2, ..., Rn) or undirectional (U1, U2, ..., Um). Number of F and R states must be equal and twin states are indicated by integers in id (e.g. F1 and R1 and twins).
#' @param initGamma Gamma(0,i).
#' 
#' @keywords internal
#' @noRd
bidirObjective2 = function(params, xsi_sum, constraints, nStates, stateLabel, initGamma) {
    stateLabel = stateLabel[-grep("R", stateLabel)]
    myF = grep("F", stateLabel)
    U = grep("U", stateLabel)
    statD = c(params[length(constraints) + myF], params[length(constraints) + myF], 
        params[length(constraints) + U])
    out = 0
    newMat = bidirGetTransMat_fixed2(params, xsi_sum, constraints, nStates, statD)
    newMat[which(newMat < 1e-300)] = 1e-300
    for (l in 1:length(constraints)) {
        i = constraints[[l]][[1]][1]
        j = constraints[[l]][[1]][2]
        i2 = constraints[[l]][[2]][1]
        j2 = constraints[[l]][[2]][2]
        if (i != i2 & j != j2) {
            l1 = log(newMat[i, j]) * xsi_sum[i, j]
            l2 = log(newMat[i2, j2]) * xsi_sum[i2, j2]
            out = out + l1 + l2
        } else {
            out = out + log(newMat[i, j]) * xsi_sum[i, j]
        }
        
        
    }
    
    ## steady-state distribution
    for (i in 1:length(statD)) {
        out = out + log(statD[i]) * initGamma[i]
    }
    
    
    return(-out)
    
}

#' Internal function that is called from C++ and optimizes transitions using the Rsolnp package.
#' @param pars Parameters for numerical optimization using Rsolnp.
#' 
#' @keywords internal
#' @noRd
c2solnp2 = function(pars) {
    
    nStates = pars$nStates
    out = list(transMat = pars$transMat, x0 = pars$x0, statD = pars$statD, doit = as.integer(0))
    oldval = bidirObjective2(pars$x0, pars$xsi_sum, pars$constraints, pars$nStates, 
        pars$stateLabel, pars$initGamma)
    # pars$xsi_sum = pars$xsi_sum + pars$pcount
    pars$couples = pars$couples + 1
    couples = pars$couples
    g = apply(pars$xsi_sum, 1, sum)
    statD = g
    for (i in 1:length(couples)) {
        if (couples[i] != i) {
            statD[i] = (g[i] + g[couples[i]])/2
        }
    }
    statD = statD/sum(statD)
    
    stateLabel = pars$stateLabel
    myF = grep("F", stateLabel)
    U = grep("U", stateLabel)
    statD = statD[c(myF, U)]
    newx0 = c(rep(1/nStates, length(pars$constraints)), statD)
    xsi_sum = pars$xsi_sum
    
    ### 
    constraints = pars$constraints
    xsi_sum_p = pars$xsi_sum
    gamma_p = apply(xsi_sum_p, 1, sum)
    gamma_m = apply(xsi_sum_p, 2, sum)[couples]
    # print(couple) print(gamma_p+gamma_m - (gamma_m[couples]+gamma_p[couples]))
    
    xsi_sum_m = matrix(0, nrow = nStates, ncol = nStates)
    for (i in 1:length(constraints)) {
        xsi_sum_m[constraints[[i]][[1]][1], constraints[[i]][[1]][2]] = xsi_sum_p[constraints[[i]][[2]][1], 
            constraints[[i]][[2]][2]]
        xsi_sum_m[constraints[[i]][[2]][1], constraints[[i]][[2]][2]] = xsi_sum_p[constraints[[i]][[1]][1], 
            constraints[[i]][[1]][2]]
    }
    
    newMat = matrix(NA, nrow = nStates, ncol = nStates)
    newx0 = c()
    for (i in 1:length(constraints)) {
        newMat[constraints[[i]][[1]][1], constraints[[i]][[1]][2]] = (xsi_sum_p[constraints[[i]][[1]][1], 
            constraints[[i]][[1]][2]] + xsi_sum_m[constraints[[i]][[1]][1], constraints[[i]][[1]][2]])/(gamma_p[constraints[[i]][[1]][1]] + 
            gamma_m[constraints[[i]][[1]][1]])
        newMat[constraints[[i]][[2]][1], constraints[[i]][[2]][2]] = (xsi_sum_p[constraints[[i]][[2]][1], 
            constraints[[i]][[2]][2]] + xsi_sum_m[constraints[[i]][[2]][1], constraints[[i]][[2]][2]])/(gamma_p[constraints[[i]][[2]][1]] + 
            gamma_m[constraints[[i]][[2]][1]])
        newx0[i] = newMat[constraints[[i]][[1]][1], constraints[[i]][[1]][2]]
        
    }
    pars$x0 = newx0
    x = newMat
    y = x
    numiter = 0
    # print('hiaa')
    while (numiter < 10000) {
        numiter = numiter + 1
        y <- y %*% x
    }
    statD = y[1, ]
    # statD = (statD[couple]+statD)/2
    statD = statD/sum(statD)
    # print(max(abs(eval_constr( newMat, xsi_sum, constraints, statD )))) print(statD
    # %*% newMat - statD)
    statD = statD[c(myF, U)]
    newx0 = c(pars$x0, statD)
    ### 
    est = t(apply(pars$xsi_sum, 1, function(x) x/sum(x)))
    # #newx0 = c(estimate, statD) #newx0 = pars$x0 for(i in
    # 1:length(pars$constraints)) { newx0[i] = est[pars$constraints[[i]][[1]][1],
    # pars$constraints[[i]][[1]][2]] } print(newx0)
    rem = c()
    # n = sum(pars$xsi_sum[pars$constraints[[c]][[1]][1],
    # pars$constraints[[c]][[1]][2]]) transProbEM = t(apply(pars$xsi_sum, 1,
    # function(x) x/sum(x)))
    nobs = sum(xsi_sum)
    nobs = 1
    for (c in 1:length(pars$constraints)) {
        if (pars$xsi_sum[pars$constraints[[c]][[1]][1], pars$constraints[[c]][[1]][2]] < 
            nobs * 1e-06 & pars$xsi_sum[pars$constraints[[c]][[2]][1], pars$constraints[[c]][[2]][2]] < 
            nobs * 1e-06) {
            # if(est[pars$constraints[[c]][[1]][1], pars$constraints[[c]][[1]][2]] < 1e-6 &
            # est[pars$constraints[[c]][[2]][1], pars$constraints[[c]][[2]][2]] < 1e-6) {
            rem = c(rem, c)
        }
    }
    removed = list()
    if (length(rem) > 0) {
        # print(rem)
        removed = pars$constraints[rem]
        pars$constraints = pars$constraints[setdiff(1:length(pars$constraints), rem)]
        pars$x0 = pars$x0[setdiff(1:length(pars$x0), rem)]
        newx0 = c(pars$x0[1:length(pars$constraints)], statD)
    }
    constraints = pars$constraints
    # oldval = bidirObjective2(pars$x0, pars$xsi_sum, pars$constraints, pars$nStates,
    # pars$stateLabel, pars$initGamma)
    
    UB = c(pars$UB[1:length(pars$constraints)], 1.25 * statD)  #
    LB = c(pars$LB[1:length(pars$constraints)], 0.75 * statD)  # rep(0, length(UB))#
    # print(newx0)
    res = solnp(newx0, fun = bidirObjective2, eqfun = bidirEqB2, eqB = pars$eqB, 
        xsi_sum = pars$xsi_sum, initGamma = pars$initGamma, stateLabel = pars$stateLabel, 
        constraints = constraints, nStates = pars$nStates, LB = LB, UB = UB, control = pars$control)
    
    
    
    out[["nrm"]] = as.integer(length(rem))
    # cat('objective diff:', res$values[length(res$values)], '-',
    # pars$objective[length(pars$objective)], '=',
    # res$values[length(res$values)]-pars$objective[length(pars$objective)], '\n')
    if (FALSE) {
        cat("Removing ", length(rem), " constraints.\n", sep = "")
        
        cat("\nCalling Rsolnp for optimization of transition matrix.\n")
        cat("Min. value of sum(xsi)=", min(pars$xsi_sum), "\n")
        cat("Max. value of sum(xsi)=", max(pars$xsi_sum), "\n")
        cat("Removing ", length(rem), " constraints.\n", sep = "")
        
        if (res$convergence == 0) {
            cat("SUCCESS: Rsolnp converged with desired accuracy!\n")
        } else {
            cat("ERROR: Rsolnp did not converge!\n")
        }
    }
    
    rel_increase = (pars$objective[length(pars$objective)] - res$values[length(res$values)])/abs(pars$objective[length(pars$objective)])
    if (pars$objective[length(pars$objective)] == Inf) {
        rel_increase = Inf
    }
    # print(rel_increase)
    
    # if(rel_increase < 0 & res$convergence == 0 & pars$nrm[length(pars$nrm)] ==
    # length(rem)) { ## BAD CONVERGENCE out[['objective']] =
    # as.numeric(pars$objective[length(pars$objective)]) if(FALSE) { cat('Objective
    # function (optimization) already converged to optimal solution in previous
    # iteration. (=> using old estimate) \n') #cat('objective diff:',
    # res$values[length(res$values)], '-', pars$objective[length(pars$objective)],
    # '=', d, '\n') } }
    
    
    # print(oldval-res$values[length(res$values)]) GOOD CONVERGENCE
    if (res$convergence == 0 & rel_increase > 0) {
        out$doit = as.integer(1)
        params = res$pars
        # cat('(+) ')
        stateLabel = pars$stateLabel
        stateLabel = stateLabel[-grep("R", stateLabel)]
        myF = grep("F", stateLabel)
        U = grep("U", stateLabel)
        statD = c(params[length(constraints) + myF], params[length(constraints) + 
            myF], params[length(constraints) + U])
        newMat = bidirGetTransMat_fixed2(params, pars$xsi_sum, pars$constraints, 
            pars$nStates, statD)
        d = res$values[length(res$values)] - pars$objective[length(pars$objective)]
        out$transMat = newMat
        out$x0 = res$pars
        out[["objective"]] = as.numeric(res$values[length(res$values)])
        out[["statD"]] = statD
        # print(newMat)
        
        if (FALSE) {
            cat("Rsolnp converged after ", length(res$values), " iterations.\n", 
                sep = "")
            cat("objective diff: ", res$values[length(res$values)], "-", pars$objective[length(pars$objective)], 
                "=", d, "\n")
            cat("Max. violation of row sum constraints: ", max(abs(1 - apply(newMat, 
                1, sum))), "\n", sep = "")
            cat("Violation of sum(stat. distribution): ", abs(1 - sum(statD)), "\n", 
                sep = "")
            cat("Max. violation of steady-state property: ", max(abs(statD %*% newMat - 
                statD)), "\n", sep = "")
            cat("Max. violation of symmetry constraints: ", max(abs(eval_constr(newMat, 
                xsi_sum, constraints, statD))), "\n", sep = "")
            cat("Stationary distribution: ", paste(round(statD, 3), collapse = " "), 
                "\n")
        }
        
        
        
        # cat('Stationary distribution - fixed_statD (computed from gamma): ',
        # paste(round(statD, 5)-round(fixed_statD, 5), collapse=' '), '\n') NO/BAD
        # CONVERGENCE
    } else {
        if (FALSE) {
            cat("Rsolnp did not converge => Using transition matrix from previous step.\n")
        }
        out[["objective"]] = as.numeric(pars$objective[length(pars$objective)])
    }
    # print(out$doit)
    out
} 
