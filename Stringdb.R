# Load the libraries
 library(STRINGdb)

# load string database
#9940 for sheep
string_db <- STRINGdb$new(version = "11.5", species = 9940, score_threshold = 200, input_directory="")
class(string_db)

# map the gene names to the STRING database identifiers using the map method. The map function adds an additional column with STRING identifiers to the dataframe that is passed as first parameter.
example1_mapped <- string_db$map(resSig, "gene", removeUnmappedRows = TRUE )
dim(example1_mapped)
hits <- example1_mapped$STRING_id[1:200]
string_db$plot_network(hits)
