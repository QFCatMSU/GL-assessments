# selectivity
selectivity <- function(sel_ctl, input, sel_pars, n_year, n_age) 
{
  sel <- matrix(0, nrow = n_year, ncol = n_age)
  
  # logistic
  if(sel_ctl == 0) 
  {
    sel <- 1 / (1 + exp(-sel_pars[1] * (input - sel_pars[2])))
  } 
  else if(sel_ctl == 1)  # gamma
  {
    p <- 0.5 * (sqrt(sel_pars[2]^2 + 4 * sel_pars[1]^2) - sel_pars[2])
    sel <- (input / sel_pars[2])^(sel_pars[2] / p) * exp((sel_pars[2] - input) / p)
  }

  return(sel)
}
