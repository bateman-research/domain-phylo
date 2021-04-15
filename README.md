# Pyhlogenetic analysis of protein domains

Code to create phylogenetic analyses of protein domains across a selected subset of species in the tree of life.

### Software requirements

- HMMER v3 : http://hmmer.org/
- `R` v3.6 with packages `seqinr`, `argparse`, `dplyr` and `tidyr`
- iTOL (website): https://itol.embl.de/


### Input files

1. Download (e.g. from Pfam) or create HMM models (with HMMER) for the protein domains of interest, such as the example file `BRF1.hmm`.
It is also possible to use a single sequence in FASTA format by converting it to an HMM using HMMER with:

```
hmmbuild sequence.hmm sequence.fasta
```

2. Create a list of species and associated proteome IDs from UniProt (https://www.uniprot.org/proteomes/), similar to the one provided in this repository for a subset of representative species in `tol_species.tsv`.
You can generate the list of species from the TSV file using:
```
grep -v "#" tol_species.tsv | cut -f1 > tol_species.id
```

3. Download the taxonomic tree of species from the NCBI Common Taxonomy Browser (https://www.ncbi.nlm.nih.gov/Taxonomy/CommonTree/wwwcmt.cgi) by uploading the list of species (only species names: `tol_species.id`) and saving it as a phylip format:  `tol_species_tree.phy`.

### Steps

4. Download proteome sequences of the selected subset of species from UniProt using `wget` as following. 
If the proteome is not a 'reference proteome' download it manually from the UniProt website (https://www.uniprot.org/proteomes/).
```
grep -v "#" tol_species.tsv | cut -f2 | while read id; do echo $id; wget "ftp://ftp.ebi.ac.uk/pub/databases/uniprot/current_release/knowledgebase/reference_proteomes/Eukaryota/${id}/${id}*.fasta.gz"; gunzip ${id}*.fasta.gz; done

rm UP0*DNA.fasta UP0*additional.fasta # remove other files often associated with proteomes (DNA and additional sequences)
```

5. Concatenate all proteome sequences into a single FASTA file and extract a list of protein IDs from all proteomes, to be used later:
```
cat UP0*.fasta > proteomes.fasta
for f in UP0*.fasta; do grep ">" $f | awk -v var="$f" '{print $1"\t"var}'; done > proteomes.tsv
```

6. Search for domain hits in all proteomes with `hmmsearch`. 
This command generates a `BRF1_proteomes_table.tsv` file for each proteome with the domain hits, and other HMMER related files such as the domain hits alignment `BRF1_proteomes.sto` and program raw output `BRF1_proteomes.log`.
Repeat this step for every domain of interest.
```
hmmsearch -o BRF1_proteomes.log -A BRF1_proteomes.sto --domtblout BRF1_proteomes.tsv BRF1.hmm proteomes.fasta; grep -v "#" BRF1_proteomes.tsv | awk '{print $1"\t"$3"\t"$4"\t"$6"\t"$13"\t"$14"\t"$16"\t"$17"\t"$20"\t"$21}' > BRF1_proteomes_table.tsv
```

7. Use the domain analysis script to create a final table of the list of species where the protein domain is present.
Use the `-h` option to see all the options.
```
Rscript domain-hits.R -p proteomes.tsv -s tol_species.tsv -t BRF1_proteomes_table.tsv -o BRF1_species.tsv -b 20 -n 1
```

8. Create a new tree on the iTOL website (https://itol.embl.de) by uploading the `tol_species_tree.phy` file and annotate species with domain hits using the Datasets tab.

### References

This repository contains code used for the phylogenetic analysis of domains in the RNA pol III enzyme, described in the following publication:

> Girbig, M., Misiaszek, A.D., Vorländer, M.K., Lafita A., Grötsch H., Baudin F., Bateman A & Müller C.W. Cryo-EM structures of human RNA polymerase III in its unbound and transcribing states. Nature Struct Mol Biol 28, 210–219 (2021). https://doi.org/10.1038/s41594-020-00555-5


