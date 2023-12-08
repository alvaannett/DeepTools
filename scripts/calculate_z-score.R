
library(tidyverse)
#library(dplyr)
#library(magrittr)
#library(tibble)
#library(stringr)

calc_z_score = function(cont, cond, mark, df){

  z = df %>%
      transmute(chr, start, end,
                mean_cont = base::rowMeans(dplyr::select(., cont), na.rm = T),
                mean_cond = base::rowMeans(dplyr::select(., cond), na.rm = T)) %>%
      mutate(z_score = (mean_cond - mean_cont)/(sqrt(mean_cond + mean_cont))) %>%
     dplyr::select(-c(mean_cont, mean_cond))

  names(z)[4] = paste("Z", mark, sep="_")
  z = z %>% mutate(key = paste(chr, start, end, sep='_'))

  return(z)
}

args = commandArgs(trailingOnly=TRUE)

print(args)

path_data = args[1]
out = args[2]
control = stringr::str_split(args[3], pattern = " ")[[1]]
condition = stringr::str_split(args[4], pattern = " ")[[1]]
mark = stringr::str_split(args[5], pattern = " ")[[1]]

print(path_data)
print(out)
print(control)
print(condition)
print(mark)

data = data.table::fread(path_data,
             sep="\t",
             header = T)

names = gsub("['|#]", "", names(data))
names(data) = names

#head(data)

z_score = data %>% dplyr::select(chr, start, end)

for(m in mark){
  samples_cont = names(data)[grepl(paste0(paste0(control, ".", m), collapse = '|'), names(data))]
  samples_cond = names(data)[grepl(paste0(paste0(condition, ".", m), collapse = '|'), names(data))]

  if(samples_cond %>% length() > 0){

      message(paste('@ calculating Z-score:', m))

      message("@ files control:")
      message(paste0(samples_cont, collapse = '\n'), "\n")

      message("@ files condition:")
      message(paste0(samples_cond, collapse = '\n'), "\n")

      tmp = calc_z_score(samples_cont, samples_cond, m, data)
      z_score = cbind(z_score, tmp %>% select(starts_with('Z_'), key))

      message("@ everything ok: ")
      message(setequal(z_score$key, paste(z_score$chr, z_score$start, z_score$end, sep = "_")))

      z_score = z_score %>% select(!key)
  }
}


head(z_score)

z_score

# save tsv and granges object
write.table(z_score,
            out,
            quote=F,
            row.names=F)

sessionInfo()



data %>% select(samples_cont, samples_cond)


data %>%
  transmute(chr, start, end,
            mean_cont = base::rowMeans(dplyr::select(., samples_cont), na.rm = T),
            mean_cond = base::rowMeans(dplyr::select(., samples_cond), na.rm = T)) %>%
  mutate(z_score = (mean_cond - mean_cont)/(sqrt(mean_cond + mean_cont))) %>%
  dplyr::select(-c(mean_cont, mean_cond))
