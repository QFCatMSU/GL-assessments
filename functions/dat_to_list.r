dat_to_list = function(path,
                          model_years,
                          model_ages,
                          model_fleets)
{
  ## get full path of all files in data directory
  dat_files = list.files(path = path, 
                        full.names = TRUE)

  ## create a list to hold data
  dat = list()

  ## model information
  dat$aux$n_year = length(model_years)
  dat$aux$n_age = length(model_ages)
  dat$aux$n_fleet = length(model_fleets)
  dat$aux$years = model_years
  dat$aux$ages = model_ages
  dat$aux$fleet = model_fleets

  ## for each data file
  # name the list items using filenames (without the extension)
  item_names = tools::file_path_sans_ext(basename(dat_files))
  for(i in 1:length(dat_files))
  {
    ## read data, there is no header and data is seperated by tabs
    tab = read.csv(dat_files[i], sep="\t", header=FALSE)
    ## add data to next position on the list
    idx = length(dat)+1
    dat[[idx]] = as.matrix(tab)
    rownames(dat[[idx]]) = model_years
    if(ncol(dat[[idx]]) == dat$aux$n_age)
    {
      colnames(dat[[idx]]) = model_ages
    }
    if(ncol(dat[[idx]]) == dat$aux$n_fleet)
    {
      colnames(dat[[idx]]) = model_fleets
    }
    names(dat)[idx] = item_names[i]
  }

  return(dat)
}