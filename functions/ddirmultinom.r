# Dirichelt multinomial likelihood
ddirmultinom = function(obs, pred, input_n, theta) 
{
  dir_param = theta * input_n
  nll = lgamma(input_n + 1) + lgamma(dir_param) - lgamma(input_n + dir_param)
  nll2 = 0
  for (a in 1:length(pred)) 
  {
    nll2 = nll2 - lgamma(obs[a] * input_n + 1) -
      lgamma(obs[a] * input_n + dir_param * pred[a]) +
      lgamma(dir_param * pred[a])
  }
  nll = nll - nll2
  return(nll)
}
