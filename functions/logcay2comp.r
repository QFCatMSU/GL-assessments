logcay2comp = function(logcay) 
{
  cay = exp(logcay)
  out = cay / sum(cay)
  
  return(out)
}
