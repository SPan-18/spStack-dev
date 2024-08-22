# internal function: finds spatial correlation matrix
spCor <- function(D, cor.fn, spParams){

  if(cor.fn == "exponential"){
    if (length(spParams) > 1)
      warning("Only first element of spParams used as phi.")
      phi <- spParams[1]
      R <- exp(-phi * D)
    }else if(cor.fn == "matern"){
      if (!length(spParams) == 2)
        stop("spParams must contain phi and nu.")
        phi <- spParams[1]
        nu <- spParams[2]
        R <- (D * phi)^nu/(2^(nu - 1) * gamma(nu)) *
             besselK(x = D * phi, nu = nu)
        diag(R) <- 1
    }else{
      stop("error: in spCor, specified cor.fn '", cor.fn, "' is not a valid
           option; choose, from exponential, matern.")
    }

    return(R)

}
