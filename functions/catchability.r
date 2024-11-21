catchability = function(qt_ctl,
                        log_q = NULL,
                        log_qt = NULL,
                        qt_sd = NULL,
                        n_year)
{
  if(qt_ctl == "time-varying")
  {
    nll = 0;
    for(t in 2:n_year)
    {
      nll = nll - dnorm(log_qt[t], log_qt[t - 1], qt_sd, log = TRUE);
    }
    out = nll;
    
  }else if(qt_ctl == "single")
  {
    log_qt = rep(0, n_year);
    for(t in 1:n_year)
    {
      log_qt[t] = log_q;
    }
    out = log_qt;
  }
    
  return(out);
}
