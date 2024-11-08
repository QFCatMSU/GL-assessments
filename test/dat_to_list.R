rm(list = ls());

## get full path of all files in data directory
datFiles = list.files(path="1836_Treaty_Waters/dat_files", 
                      full.names = TRUE);

## create a list to hold data
all_data = list();

## for each data file
for(i in 1:length(datFiles))
{
  ## read data, there is no header and data is seperated by tabs
  dat = read.csv(datFiles[i], sep="\t", header=FALSE);
  ## add data to next position on the list
  all_data[[length(all_data)+1]] = as.matrix(dat);
}

## name the list items using filenames (without the extension)
names(all_data) = tools::file_path_sans_ext(basename(datFiles));