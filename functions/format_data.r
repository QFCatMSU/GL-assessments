format_data = function(d, type_fleet) {
  obs = list()
  n_data = length(d)
  # index for data type
  ct_grep = grep("ct|cpe", type_fleet$type)
  aa_grep = grep("iaa|pa", type_fleet$type)

  # loop through datasets
  for (i in 1:n_data)
  {
    if (i %in% ct_grep) 
    {
      obs[[i]] = d[[i]]
      obs[[i]] = data.frame(obs[[i]])
      names(obs[[i]]) = "obs"
      obs[[i]]$year = rownames(d[[i]])
      rownames(obs[[i]]) = NULL
      obs[[i]]$age = NA
      obs[[i]]$type = type_fleet$type[i]
      obs[[i]]$fleet = type_fleet$fleet[i]
      obs[[i]]$obs = log(obs[[i]]$obs)
    } else if (i %in% aa_grep) 
    {
      obs[[i]] = reshape(
        data = d[[i]],
        direction = "long",
        varying = c(1:ncol(d[[i]])),
        v.name = "obs",
        times = colnames(d[[i]]),
        timevar = "age",
        ids = rownames(d[[i]]),
        idvar = "year"
      )
      rownames(obs[[i]]) = NULL
      obs[[i]]$type = type_fleet$type[i]
      obs[[i]]$fleet = type_fleet$fleet[i]
      obs[[i]] = obs[[i]][, c("obs", "year", "age", "type", "fleet")]
      if (type_fleet$type[i] == "iaa") 
      {
        obs[[i]]$obs = log(obs[[i]]$obs)
      }
    }
  }

  # combine
  obs_all = do.call(rbind, obs)
  num_idx = which(colnames(obs_all) %in% c("obs", "year", "age", "fleet"))
  for (i in num_idx) 
  {
    obs_all[, i] = as.numeric(obs_all[, i])
  }

  return(obs_all)
}
