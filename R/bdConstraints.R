
#' Internal function that generates constraints for a bidirection transition matrix.
#' @param state.type  Indicates directinality of states. States can be forward (F1, F2, ..., Fn), reverse (R1, R2, ..., Rn) or undirectional (U1, U2, ..., Um). Number of F and R states must be equal and twin states are indicated by integers in id (e.g. F1 and R1 and twins).
#' 
#' @return A list containing the constraints
#' @keywords internal
#' @noRd

getConstraints = function(state.type) {
    fow = grep("F", state.type)
    nfow = sapply(strsplit(state.type[fow], "F"), function(x) as.integer(x[2]))
    rev = grep("R", state.type)
    nrev = sapply(strsplit(state.type[rev], "R"), function(x) as.integer(x[2]))
    u = grep("U", state.type)
    nu = sapply(strsplit(state.type[u], "U"), function(x) as.integer(x[2]))
    fow2rev = rep(0, length(nfow))
    for (i in 1:length(nfow)) {
        fow2rev[i] = rev[which(nfow[i] == nrev)]
    }
    rev2fow = rep(0, length(nrev))
    for (i in 1:length(nrev)) {
        rev2fow[i] = fow[which(nrev[i] == nfow)]
    }
    if (length(fow) != length(rev)) {
        stop("Number of foward and reverse states need to be equal!\n")
    }
    constraints = list()
    for (i in 1:length(state.type)) {
        for (j in 1:length(state.type)) {
            if (i %in% fow & j %in% fow) {
                constraints[[length(constraints) + 1]] = list(FF = c(fow[fow == i], 
                  fow[fow == j]), RR = c(fow2rev[fow[fow == j]], fow2rev[fow[fow == 
                  i]]))
            } else if (i %in% fow & j %in% rev) {
                if (fow[fow == i] <= rev2fow[rev == j]) 
                  constraints[[length(constraints) + 1]] = list(FR = c(fow[fow == 
                    i], rev[rev == j]), FR = c(rev2fow[rev == j], fow2rev[fow == 
                    i]))
            } else if (i %in% rev & j %in% fow) {
                if (rev[rev == i] <= fow2rev[fow == j]) 
                  constraints[[length(constraints) + 1]] = list(RF = c(rev[rev == 
                    i], fow[fow == j]), RF = c(fow2rev[fow == j], rev2fow[rev == 
                    i]))
            } else if (i %in% u & j %in% fow) {
                constraints[[length(constraints) + 1]] = list(UF = c(u[which(u == 
                  i)], fow[which(fow == j)]), RU = c(fow2rev[fow[fow == j]], u[which(u == 
                  i)]))
            } else if (i %in% fow & j %in% u) {
                constraints[[length(constraints) + 1]] = list(FU = c(fow[which(fow == 
                  i)], u[which(u == j)]), UR = c(u[which(u == j)], fow2rev[which(fow == 
                  i)]))
            } else if (i %in% u & j %in% u) {
                if (j >= i) {
                  constraints[[length(constraints) + 1]] = list(UU = c(u[which(u == 
                    i)], u[which(u == j)]), UU = c(u[which(u == j)], u[which(u == 
                    i)]))
                }
            }
        }
    }
    constraints
    constraints
} 
