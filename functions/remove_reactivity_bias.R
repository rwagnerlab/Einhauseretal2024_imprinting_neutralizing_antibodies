#this script uses functions and code from "BA.2 and BA.5 omicron differ immunologically from both BA.1 omicron and pre-omicron variants" by Roessler, A., Netzl, A., et al. 2022
#available on repository: https://github.com/acorg/roessler_netzl_et_al2022/
#DOI: 10.5281/zenodo.7341691..
# remove reactivity bias from titertable
remove_reactivity_bias_logtiter <- function(table) {
  row_names <- rownames(table)
  table <- sapply(as.data.frame(table), function(x) {
    x - (mean(x, na.rm = T) - mean(table, na.rm = T))
  })
  
  rownames(table) <- row_names
  
  return(table)
}