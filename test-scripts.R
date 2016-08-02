# Test scripts to get PubMed data from NCBI into a usable format.

############################################################
#                                                          #
#                      Load packages                       #
#                                                          #
############################################################
library(rentrez)
library(purrr)

############################################################
#                                                          #
#                 Get PubMed search fields                 #
#                                                          #
############################################################
entrez_db_searchable(db = 'pubmed', config = NULL)

############################################################
#                                                          #
#                  Set search parameters                   #
#                                                          #
############################################################
# NLM journal [JOUR] abbreviations: Pain (est: 1975), J Pain (est: 2000),
# Eur J Pain (est: 1997), Clin J Pain (est: 1985)
# Year of publication [PDAT]
# Publication type [PTYP]
srch_terms <- c('((Pain[Journal] NOT Pain Suppl[Journal])
                OR J Pain[Journal] OR Eur J Pain[Journal]
                OR Clin J Pain[Journal])
                AND 1975[PDAT] : 2015[PDAT]
                AND Journal Article[PTYP]')

############################################################
#                                                          #
#                    Explore search results                #
#                                                          #
############################################################
# Find out how many articles are returned by search
explore_srch <- entrez_search(db = 'pubmed',
                             term = srch_terms,
                             retmode = 'xml',
                             retmax = 0)
# Use 'count' data from 'explore_srch' to set the 'retmax' in a new search
id_srch <- entrez_search(db = 'pubmed',
                          term = srch_terms,
                          retmode = 'xml',
                          retmax = explore_srch$count)

############################################################
#                                                          #
#               Get and sort publication IDs               #
#                                                          #
############################################################
# Extract all publication ids into numeric vector 'ids'
ids <- id_srch$ids
# Define a sequence to split the vector of ids into n = 100 chunks
splitter <- seq(from = 1, to = id_srch$count, by = 100)
# Define an empty list of length 'splitter'
id_list <- vector(mode = 'list', length = length(splitter))
# For loop to split 'ids'
for(i in seq_along(splitter)) {
    id_list[[i]] <- ids[splitter[[i]]:(splitter[[i]] + 99)]
    }

############################################################
#                                                          #
#                  Download full records                   #
#                                                          #
############################################################
# For testing use a subsample of id_list
id_test <- id_list[1:2]
id_test[[1]] <- id_test[[1]][1:20]
id_test[[2]] <- id_test[[2]][1:20]
test_rec <- id_test %>%
    map(entrez_fetch, db = 'pubmed', rettype = 'xml') %>%
    map(parse_pubmed_xml) %>%
    map(1:2)
