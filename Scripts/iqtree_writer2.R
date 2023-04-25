# To run: Rscript path_to/iqtree_writer2.R PIPE_PATH
# where PIPE_PATH is directory containing "storage", and "in" folders for executing GMPipeline 
# Not designed to be run on its own
# Purpose is to read results from optimization trees using just master_seq and outgroup sequences 
# to determine optimal settings and ensure trees look correct. Identifies best substitution model.

library(stringr)

# Setting filepaths
PIPE_PATH <- commandArgs(trailingOnly=TRUE)[1]

# Calculating best model using model finder results for unconstrained trees
n=1
best_models <- vector()
while (n<=10){
 model_file <- paste0(PIPE_PATH, "/storage/opt_ML/main_unconstr", n, ".log")
 best_line <- grep("Best-fit", readLines(model_file), value=TRUE)
 best_models[length(best_models)+1] <- str_split(str_split(best_line, ": ")[[1]][2], " ")[[1]][1]
 n=n+1
}
best_model <- names(sort(table(best_models),decreasing=TRUE))[1]

# Parsing comparison file
comparison_lines <- readLines(paste0(PIPE_PATH, "/storage/opt_ML/concat/main_comparison.iqtree"))
AU_index <- grep("p-AU", comparison_lines)[1]
index <- AU_index+2
unconstr_scores <- ""
while (index <= AU_index+11) {
  unconstr_line <- str_remove_all(comparison_lines[index], " ")
  unconstr_scores <- paste0(unconstr_scores, str_sub(unconstr_line, start = -1))
  index <- index+1
}
constr_scores <- ""
while (index <= AU_index+21) {
  constr_line <- str_remove_all(comparison_lines[index], " ")
  constr_scores <- paste0(constr_scores, str_sub(constr_line, start = -1))
  index <- index+1
}

# Verifying that constrained and unconstrained trees pass statistical tests 
unconstr_pass <- str_count(unconstr_scores, "[+]")
constr_pass <- str_count(constr_scores, "[+]")
if(isTRUE(constr_pass<9) | isTRUE(unconstr_pass<9)) {stop('reference tree failed statistical tests: ingroups and outgroups may not significantly separated in tree, please check storage/opt_ML/main_unconstr/*.iqtree files and adjust reference and/or outgroup sequences')}

# Reporting best model 
write(best_model, file=paste0(PIPE_PATH, "/storage/opt_ML/best_model.txt"))
