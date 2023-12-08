
library(tidyverse)
#library(dplyr)
#library(magrittr)
#library(tibble)
#library(stringr)

# mutate_at(vars(-c('chr', 'start', 'end')), ~(.+1)) %>%

calc_log2fc = function(cont, cond, mark, df){

  log2fc = df %>%
       transmute(chr, start, end,
                mean_cont = base::rowMeans(dplyr::select(., cont), na.rm = T),
                mean_cond = base::rowMeans(dplyr::select(., cond), na.rm = T)) %>%
      mutate(log2fc = log2(mean_cond) - log2(mean_cont)) %>%
      dplyr::select(-c(mean_cont, mean_cond))

  names(log2fc)[4] = paste("log2fc", mark, sep="_")
  log2fc = log2fc %>% mutate(key = paste(chr, start, end, sep='_'))

  return(log2fc)
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

log2fc = data %>% dplyr::select(chr, start, end)

for(m in mark){
  samples_cont = names(data)[grepl(paste0(paste0(control, ".", m), collapse = '|'), names(data))]
  samples_cond = names(data)[grepl(paste0(paste0(condition, ".", m), collapse = '|'), names(data))]

  if(samples_cond %>% length() > 0){

      message(paste('@ calculating log2fc', m))

      message("@ files control:")
      message(paste0(samples_cont, collapse = '\n'), "\n")

      message("@ files condition:")
      message(paste0(samples_cond, collapse = '\n'), "\n")

      tmp = calc_log2fc(samples_cont, samples_cond, m, data)
      log2fc = cbind(log2fc, tmp %>% select(starts_with('log2fc_'), key))

      message("@ everything ok: ")
      message(setequal(log2fc$key, paste(log2fc$chr, log2fc$start, log2fc$end, sep = "_")))

      log2fc = log2fc %>% select(!key)
  }
}


head(log2fc)

log2fc

# save tsv and granges object
write.table(log2fc,
            out,
            quote=F,
            row.names=F)

sessionInfo()

