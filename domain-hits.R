# Analyze the evolution of a protein domain across species
# Aleix Lafita - March 2021

suppressPackageStartupMessages(library(argparse))
suppressPackageStartupMessages(library(seqinr))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(tidyr))

# Columns names of the HMMER table
hmmer.cols = unlist(strsplit("qid,qlen,tid,tlen,evalue,bitscore,qstart,qend,start,end", ","))

###################### Argparse #############################

proteomes = "proteomes.tsv"
output = "domain_species.tsv"
species = "tol_species.tsv"
hmmer = "BRF1_proteomes_table.tsv"
number = 1
bit_thr = 20

# create parser object
parser = ArgumentParser(description = 'Analyze the evolution of a protein domain across species.')

# specify our desired options 
parser$add_argument("-p", "--proteomes", default=proteomes,
                    help="Input file of proteome IDs mapped to sequence IDs [default \"%(default)s\"]")
parser$add_argument("-s", "--species", default=species,
                    help="Table of species and proteome IDs [default \"%(default)s\"]")
parser$add_argument("-t", "--hmmer", default=hmmer,
                    help="Table of HMMER hits for the domain of interest [default \"%(default)s\"]")
parser$add_argument("-o", "--output", default=output,
                    help="Name of the output file of domain hits for species [default \"%(default)s\"]")
parser$add_argument("-b", "--bit_thr", default=bit_thr,
                    help="Score threshold in bits for domain hits \"%(default)s\"]")
parser$add_argument("-n", "--number", default=number,
                    help="Number of expected hits per proteome, in case of paralogs \"%(default)s\"]")

# get command line options, if help option encountered print help and exit,
# otherwise if options not found on command line then set defaults, 
args = parser$parse_args()

proteomes = args$proteomes
species = args$species
output = args$output
hmmer = args$hmmer
bit_thr = as.integer(args$bit_thr)
number = as.integer(args$number)

############################ Parse data #################################

# Parse protein to proteome mapping
proteomes.info = read.csv(
  proteomes,
  sep = "\t",
  col.names = c("uniprotid", "proteome")
) %>% mutate(
  uniprotid = gsub("\\|.*", "", gsub(">..\\|", "", uniprotid)),
  proteome = gsub("_.*", "", proteome),
  proteome = gsub("\\.fasta", "", proteome),
)

# Parse proteome to species name mapping
species.tol = read.csv(
  species,
  sep = "\t",
  header = F,
  comment.char = "#",
  col.names = c("species", "proteome")
)

# Parse the domain hits
domain.hmmer = read.csv(
  hmmer,
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
  merge(proteomes.info) %>% 
  merge(species.tol) %>%
  group_by(proteome) %>%
  top_n(number, bitscore) %>% # select top scoring match
  top_n(number, -evalue) %>% # break any ties using the e-value
  top_n(number, uniprotid) %>% # in case of duplicated proteins
  filter(bitscore > bit_thr) %>% # filter low scoring hits
  ungroup() %>% mutate(domain = tid)

write.table(
  domain.species %>% select(species, domain, uniprotid, bitscore),
  output,
  sep = "\t",
  quote = F,
  row.names = F
)

