recruitment = function(rec_ctl, log_r, r_sd, phi = NULL) 
{
  if (rec_ctl == "WN") # WN
  {
    r_nll = -sum(dnorm(log_r, 0, r_sd, log = TRUE))
  }else if (rec_ctl == "RW") # RW
  {
    for (t in 2:(length(log_r) - 1))
    {
      r_nll = -dnorm(log_r[t], log_r[t - 1], r_sd, log = TRUE)
    }
    }else if (rec_ctl == "AR1" | !is.null(phi)) # AR-1
    {
      stationary_sd = sqrt(r_sd * r_sd / (1 - phi * phi))
      r_nll = -dautoreg(log_r, phi = phi, scale = stationary_sd, log = TRUE)
  }
  
  return(r_nll)
}
