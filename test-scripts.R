library(rentrez)
# Get PubMed search fields
entrez_db_searchable(db = 'pubmed', config = NULL)

# NLM journal [JOUR] abbreviations: Pain (est: 1975), J Pain (est: 2000),
# Eur J Pain (est: 1997), Clin J Pain (est: 1985)
# Year of publication [PDAT]
# Publication type [PTYP]
srch_terms <- c('((Pain[Journal] NOT Pain Suppl[Journal])
                OR J Pain[Journal] OR Eur J Pain[Journal]
                OR Clin J Pain[Journal])
                AND 2015[PDAT] : 2015[PDAT]
                AND Journal Article[PTYP]')
# Search
pubmed_srch <- entrez_search(db = 'pubmed',
                             term = srch_terms,
                             retmode = 'xml',
                             retmax = 0,
                             use_history = TRUE)

# fetch
pubmed_fetch <- list()
for(seq_start in seq(1, 1000 , 100)) {
    for(i in seq_along(seq(1, 1000 , 100))) {
    pubmed_fetch[[i]] <- entrez_fetch(db = 'pubmed',
                                web_history = pubmed_srch$web_history,
                                rettype = 'xml',
                                parsed = TRUE,
                                retmax = 100,
                                restart = seq_start)
    }
}


# Parse
pubmed_parsed <- purrr::map(.x = x, .f = parse_pubmed_xml)
