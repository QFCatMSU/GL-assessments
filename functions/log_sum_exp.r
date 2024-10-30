log_sum_exp = function(x) 
{
  sum_exp = sum(exp(x))   
  out = log(sum_exp)      

  return(out)       
}