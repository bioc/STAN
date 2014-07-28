
#' Internal function which calculates log density of the Inverse-Wishart distribution.
#' 
#' @keywords internal
#' @noRd

logdiwish = function(W, v, S) {
    nRows = nrow(S)
    myGamma = 1
    for (i in 1:nRows) {
        myGamma <- myGamma * gamma((v + 1 - i)/2)
    }
    denom = log(myGamma * 2^(v * nRows/2) * pi^((nRows * (nRows - 1)/4)))
    detS = det(S)
    detW = det(W)
    matProd = S %*% solve(W)
    myTrace = sum(matProd[row(matProd) == col(matProd)])
    numer = log(detS) * (v/2) + log(detW) * (-(v + nRows + 1)/2) + (-1/2 * myTrace)
    return(numer - denom)
} 
