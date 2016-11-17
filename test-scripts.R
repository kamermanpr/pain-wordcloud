# Test scripts to get PubMed data from NCBI into a usable format.

############################################################
#                                                          #
#                      Load packages                       #
#                                                          #
############################################################
library(rentrez)
library(stringr)
library(tidyverse)

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
search_terms <- c('((Pain[Journal] NOT Pain Suppl[Journal]) OR J Pain[Journal] OR Eur J Pain[Journal] OR Clin J Pain[Journal]) AND 1975[PDAT] : 2015[PDAT] AND Journal Article[PTYP]')

############################################################
#                                                          #
#                    Explore search results                #
#                                                          #
############################################################
# Find out how many articles are returned by search
explore_search <- entrez_search(db = 'pubmed',
                             term = search_terms,
                             retmode = 'xml',
                             retmax = 0)$count

# Use 'count' data from 'explore_srch' to set the 'retmax' in a new search
id_search <- entrez_search(db = 'pubmed',
                          term = search_terms,
                          retmode = 'xml',
                          retmax = explore_search)

############################################################
#                                                          #
#               Get and sort publication IDs               #
#                                                          #
############################################################
# Extract all publication ids into numeric vector 'ids'
ids <- id_search$ids

# Define a sequence to split the vector of ids into n = 100 chunks
splitter <- seq(from = 1, to = id_search$count, by = 100)

# Define an empty list of length 'splitter'
id_list <- vector(mode = 'list', length = length(splitter))

# For loop to split 'ids'
for(i in 1:length(splitter)) {
    id_list[[i]] <- ids[splitter[[i]]:(splitter[[i]] + 99)]
}

############################################################
#                                                          #
#                  Download full records                   #
#                                                          #
############################################################
# For testing use a subsample of id_list
id <- id_list[1:2]
id[[1]] <- id[[1]][1:20]
id[[2]] <- id[[2]][1:20]
record <- id %>%
    map(entrez_fetch, db = 'pubmed', rettype = 'xml') %>%
    map(parse_pubmed_xml) %>%
    at_depth(2, unlist) %>%
    at_depth(2, t) %>%
    at_depth(2, as_data_frame)

# Create dummy abstract column to add to all dfs
# (some records don't have an abstract column, which
# creates problems in the next step)
df_dummy <- data_frame(abstract = NA)

# Add dummy abstract column, collapse abstract columns,
# bind rows within each sub-list, bind rows between
# level 1 list elements and convert into df
record <- record %>%
    at_depth(2, bind_cols, df_dummy) %>%
    at_depth(2, unite, abstract, starts_with('abstract'), sep = ' ') %>%
    at_depth(1, bind_rows) %>%
    map_df(as_data_frame)

# Remove keyword columns and rows with empty abstract ('NA')
record <- record %>%
    select(-starts_with('key_words')) %>%
    filter(abstract != 'NA') %>%
    # Add study number as a reference
    mutate(study_id = 1:nrow(.))

# Separate out authors names and clean-up initials
authors <- record %>%
    select(starts_with('authors'), study_id) %>%
    gather(key = author_number, value = author, -study_id) %>%
    # Separate surname and first names and trim white space
    separate(author, into = c('surname', 'first'), sep = ',') %>%
    mutate(first = str_trim(first, side = 'both')) %>%
    # Separate initials ('first') into up to five individual initials
    separate(first, into = c('first', 'second', 'third', 'forth', 'fifth'),
             sep = ' ', remove = TRUE) %>%
    # Trim text to starting uppercase letter, unless the string
    # contains hyphenated uppercase letters (e.g., H-Y)
    mutate(first = ifelse(str_detect(first,
                                     pattern = '[A-Z]-[A-Z]'),
                          str_sub(first, 1, 3),
                          str_sub(first, 1, 1))) %>%
    mutate(second = ifelse(str_detect(second,
                                     pattern = '[A-Z]-[A-Z]'),
                          str_sub(second, 1, 3),
                          str_sub(second, 1, 1))) %>%
    mutate(third = ifelse(str_detect(third,
                                     pattern = '[A-Z]-[A-Z]'),
                          str_sub(third, 1, 3),
                          str_sub(third, 1, 1))) %>%
    mutate(forth = ifelse(str_detect(forth,
                                     pattern = '[A-Z]-[A-Z]'),
                          str_sub(forth, 1, 3),
                          str_sub(forth, 1, 1))) %>%
    mutate(fifth = ifelse(str_detect(fifth,
                                     pattern = '[A-Z]-[A-Z]'),
                          str_sub(fifth, 1, 3),
                          str_sub(fifth, 1, 1))) %>%
    # Unite initials columns
    unite(col = 'initials', first, second, third, forth, fifth, sep = ' ') %>%
    # Remove 'NA' from initial strings and trim any white space
    mutate(initials = str_replace_all(initials,
                                      pattern = 'NA',
                                      replacement = '')) %>%
    mutate(initials = str_trim(initials, side = 'both')) %>%
    # Arrange by study id
    arrange(study_id) %>%
    # Remove NA author name rows
    filter(!is.na(surname))

# Organise 'record' into same df structure as 'authors'
record <- record %>%
    gather(key = author_number, value = author, starts_with('author')) %>%
    # Get ride of empty author rows
    filter(!is.na(author)) %>%
    # Arrange by study id
    arrange(study_id)

# Join 'record' and 'authors'
record <- record %>%
    left_join(authors) %>%
    # Convert 'authors#' string to numeric
    mutate(author_number = as.numeric(str_sub(author_number, 8, 12))) %>%
    # Add new columns (organised by study_id) with author
    # counts and positions
    group_by(study_id) %>%
    mutate(author_count = n()) %>%
    mutate(author_position = ifelse(author_number - author_count == 0,
                                    yes = 'last author',
                                    no = ifelse(author_number == 1,
                                           yes = 'first author',
                                           no = 'middle author'))) %>%
    # Select required columns
    select(study_id, surname, initials,
           author_number, author_position, author_count,
           title, year, journal, volume, pages,
           pmid, doi, abstract)

# Final text and number clean-up
record <- record %>%
    mutate(surname = str_to_lower(surname),
           initials = str_to_lower(initials),
           title = str_to_lower(title),
           year = as.numeric(year),
           volume = as.numeric(volume),
           pmid = as.numeric(pmid),
           abstract = str_to_lower(abstract)) %>%
    mutate(journal = str_replace(journal,
                                 pattern = 'European journal of pain \\(London, England\\)',
                                 replacement = 'european journal of pain')) %>%
    mutate(journal = str_replace(journal,
                                 pattern = 'The journal of pain \\: official journal of the American Pain Society',
                                 replacement = 'journal of pain')) %>%
    mutate(journal = str_replace(journal,
                                 pattern = 'The Clinical journal of pain',
                                 replacement = 'clinical journal of pain')) %>%
    mutate(journal = str_replace(journal,
                                 pattern = 'Pain',
                                 replacement = 'pain'))

write_rds(record, './data/record.rds')

library(tm)
abstract_df <- record %>%
    select(abstract) %>%
    group_by(abstract) %>%
    summarise()

# Process
other_words <- c('may',
                 'thus',
                 'although',
                 'nevertheless',
                 'however',
                 'still',
                 'might',
                 'also',
                 'well',
                 'use',
                 'used',
                 'using',
                 'just',
                 'can',
                 'due',
                 'although',
                 'difference',
                 'differences',
                 'compare',
                 'compared',
                 'regarding',
                 'regards',
                 'will',
                 'include',
                 'included',
                 'associated',
                 'effect',
                 'effects',
                 'study',
                 'studies',
                 'studied',
                 'investigate',
                 'investigated',
                 'investigation',
                 'introduction',
                 'experiment',
                 'experiments',
                 'group',
                 'groups',
                 'show',
                 'showed',
                 'shown',
                 'showing',
                 'outcome',
                 'outcomes',
                 'find',
                 'finding',
                 'findings',
                 'found',
                 'compare',
                 'compared',
                 'aim',
                 'aims',
                 'aiming',
                 'aimed',
                 'goal',
                 'goals',
                 'measure',
                 'measured',
                 'methods',
                 'method',
                 'material',
                 'materials',
                 'result',
                 'results',
                 'conclusion',
                 'conclusions',
                 'concluded',
                 'suggest',
                 'suggested',
                 'suggesting',
                 'consequence',
                 'consequences')
abstract_tm <- DataframeSource(abstract_df)
abstract_corpus <- Corpus(abstract_tm)
abstract_corpus <- tm_map(abstract_corpus, content_transformer(tolower))
abstract_corpus <- tm_map(abstract_corpus, removePunctuation)
abstract_corpus <- tm_map(abstract_corpus, removeNumbers)
abstract_corpus <- tm_map(abstract_corpus, stripWhitespace)
abstract_corpus <- tm_map(abstract_corpus, removeWords,
                          stopwords('english'))
abstract_corpus <- tm_map(abstract_corpus, removeWords,
                          other_words)
#abstract_corpus <- tm_map(abstract_corpus, stemDocument)
abstract_corpus <- tm_map(abstract_corpus, PlainTextDocument)

# analyse
abstract_mtrx <- DocumentTermMatrix(abstract_corpus)
#abstract_mtrx2 <- removeSparseTerms(abstract_mtrx, 0.1)
freq = slam::col_sums(abstract_mtrx)
df <- data.frame(word = names(freq), freq = freq)
df



df2 <- df %>%
    filter(freq >= 10) %>%
    arrange(desc(freq))

wordcloud2::wordcloud2(df2, fontFamily = 'helvetica')

# cluster

dt <- abstract_mtrx
dt2 <- removeSparseTerms(abstract_mtrx, 0.8)
d <- dist(t(dt2), method="euclidian")
fit <- hclust(d=d, method="ward.D2")
rad_fit <- networkD3::as.radialNetwork(fit)
networkD3::radialNetwork(rad_fit)
networkD3::dendroNetwork(fit)
networkD3::diagonalNetwork(rad_fit)

ig <- igraph::graph.adjacency(d)
ig_nam <- attr(d, 'Labels')
V(ig)$name <- ig_nam
edgebundleR::edgebundle(ig)
visNetwork(vis$nodes, vis$edges, height = "500px") %>%
    visIgraphLayout(layout = "layout_in_circle") %>%
    visEdges(width = 0.5, selectionWidth = 2, color = list(color = "#97C2FC", highlight = "red")) %>%
    visNodes(size = 20) %>%
    visOptions(highlightNearest = list(enabled = T, hover = T))

