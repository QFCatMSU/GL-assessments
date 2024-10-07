# selectivity
selectivity = function(sel_ctl, input, sel_pars) 
{
  # logistic
  if(sel_ctl == "logistic") 
  {
    sel = 1 / (1 + exp(-sel_pars[1] * (input - sel_pars[2])))
  } 
  else if(sel_ctl == "gamma")  # gamma
  {
    p = 0.5 * (sqrt(sel_pars[2]^2 + 4 * sel_pars[1]^2) - sel_pars[2])
    sel = (input / sel_pars[2])^(sel_pars[2] / p) * exp((sel_pars[2] - input) / p)
  }

  return(sel)
}
