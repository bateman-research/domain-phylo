# Analyze the evolution of a protein domain across species
# Aleix Lafita - March 2021

library(seqinr)
library(dplyr)
library(tidyr)

# Columns names of the HMMER table
hmmer.cols = unlist(strsplit("qid,qlen,tid,tlen,evalue,bitscore,qstart,qend,start,end", ","))

############################ Parse data #################################

# Parse protein to proteome mapping
proteomes = read.csv(
  "proteomes.tsv",
  sep = "\t",
  col.names = c("uniprotid", "proteome")
) %>% mutate(
  uniprotid = gsub("\\|.*", "", gsub(">..\\|", "", uniprotid)),
  proteome = gsub("_.*", "", proteome)
)

# Parse proteome to species name mapping
species.tol = read.csv(
  "tol_species.tsv",
  sep = "\t",
  header = F,
  comment.char = "#",
  col.names = c("species", "proteome")
)

# Parse the domain hits
domain.hmmer = read.csv(
  "BRF1_proteomes_table.tsv",
  header = F,
  sep = "\t",
  stringsAsFactors = F,
  col.names = hmmer.cols
) %>% mutate(
  uniprotid = gsub("\\|.*", "", gsub("^..\\|", "", qid)),
  alnlen = end - start
) %>% group_by(uniprotid) %>%
  top_n(1, bitscore)

############################ Domain hits #################################

# Merge the species information with the hits
domain.species = domain.hmmer %>%
  merge(proteomes) %>% 
  merge(species.tol) %>%
  group_by(proteome) %>%
  top_n(1, bitscore) %>%
  filter(bitscore > 10) %>% # filter low scoring hits
  ungroup() %>% mutate(domain = tid)

write.table(
  domain.species %>% select(species, domain, bitscore),
  "brf1_species.tsv",
  sep = "\t",
  quote = F,
  row.names = F
)

