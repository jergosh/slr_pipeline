CLADE=Eutheria
# CLADE=Primates
# CLADE=yeast
SPECIES_LIST=data/specieslist.txt

# General setup and files to be donwloaded from Ensembl
ENS_VERSION=78
EMF_URL=ftp://ftp.ensembl.org/pub/release-$(ENS_VERSION)/emf/ensembl-compara/homologies/Compara.$(ENS_VERSION).protein.nhx.emf.gz
PR_ROOT=$(shell pwd)
ENS_ROOT=$(PR_ROOT)/data/ens/$(ENS_VERSION)
EMF_GZIP=$(ENS_ROOT)/Compara.$(ENS_VERSION).protein.nhx.gz
# EMF_FILE=$(ENS_ROOT)/Compara.$(ENS_VERSION).protein.nhx
EMF_FILE=$(ENS_ROOT)/Compara.78.protein.nhx
SPECIES_CACHE=$(ENS_ROOT)/clades.pk
TREE_ROOT=$(ENS_ROOT)/trees_greg_63
IMG_ROOT=$(ENS_ROOT)/img
SEQSETS_ROOT=$(ENS_ROOT)/seqsets_fixed
SEQSETS_CDS_ROOT=$(ENS_ROOT)/seqsets_cds
SEQSETS_PEP_ROOT=$(ENS_ROOT)/seqsets_pep
DIVERGE_ROOT=$(ENS_ROOT)/diverge
DIVERGE_SEQSETS_ROOT=$(ENS_ROOT)/diverge_seqsets_cds
DIVERGE_ALN_ROOT=$(ENS_ROOT)/aln_diverge
HMMDIVERGE_ROOT=$(ENS_ROOT)/hmm_diverge
CDS_DIR=$(ENS_ROOT)/cds
PEP_DIR=$(ENS_ROOT)/pep
PEP_HUMAN=$(PEP_DIR)/Homo_sapiens.GRCh38.pep.all.fa

ALN_ROOT=$(ENS_ROOT)/aln
PRANK_LOG_DIR=$(PR_ROOT)/log/prank
DIVERGE_PRANK_LOG_DIR=$(PR_ROOT)/log/prank_diverge
HMMDIVERGE_LOG_ROOT=$(PR_ROOT)/log/hmmdiverge

SLR_ROOT=$(ENS_ROOT)/slr
SLR_LOG_DIR=$(PR_ROOT)/log/slr

SLR_ALL=$(ENS_ROOT)/$(CLADE)_all.tab

# TODO Should this be versioned as well?
SIFTS_URL=ftp://ftp.ebi.ac.uk/pub/databases/msd/sifts/flatfiles/tsv/pdb_chain_uniprot.tsv.gz
SIFTS_GZIP=$(PR_ROOT)/data/pdb_chain_uniprot.tsv.gz
SIFTS_FILE=$(PR_ROOT)/data/pdb_chain_uniprot.tsv
UNIPROT_SEQ_URL=ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz
UNIPROT_SEQ_GZIP=$(PR_ROOT)/data/uniprot_sprot.fasta.gz
UNIPROT_SEQ_FILE=$(PR_ROOT)/data/uniprot_sprot.fasta
PDB_MAP_FILE=$(ENS_ROOT)/pdb_map_$(CLADE).tab
PDB_STATS=$(PR_ROOT)/results/pdb_stats_$(CLADE).pdf
PDB_LOG=$(PR_ROOT)/log/pdb_download.log
PDB_MISMATCH_LOG=$(PR_ROOT)/log/pdb_mismatch.log
PDB_MISSING_LOG=$(PR_ROOT)/log/pdb_missing.log

PDB_DIR=$(PR_ROOT)/data/pdb

PDB_MASTER_TABLE=$(PR_ROOT)/data/pdb_master_table_$(CLADE).tab

# YEAST_ROOT=$(PR_ROOT)/data/yeast
# TREE_ROOT=$(YEAST_ROOT)/trees
# SEQSETS_ROOT=$(YEAST_ROOT)/seqsets/
# ALN_ROOT=$(YEAST_ROOT)/aln
# SLR_ROOT=$(YEAST_ROOT)/slr
# SLR_LOG_DIR=$(PR_ROOT)/log/yeast/slr
# SLR_ALL=$(YEAST_ROOT)/$(CLADE)_all.tab

# Auxiliary output directories
SEQSETS_ROOT_UNFILTERED=$(ENS_ROOT)/seqsets_all


BSUB=bsub -I

$(EMF_GZIP):
	curl $(EMF_URL) > $(EMF_GZIP)

$(EMF_FILE): $(EMF_GZIP)
	gunzip -c $(EMF_GZIP) > $(EMF_FILE)

# TODO Write a target for any genome file? %
pep_files: # TODO depend on a wildcarded genome files
	mkdir -p $(PEP_DIR)
	python bin/get_ens_ftp.py --clade $(CLADE) --release $(ENS_VERSION) --outroot $(PEP_DIR) \
	--species_cache $(SPECIES_CACHE) --pep
	gunzip $(PEP_DIR)/*.gz

cds_files: # TODO depend on a wildcarded genome files
	mkdir -p $(CDS_DIR)
	python bin/get_ens_ftp.py --clade $(CLADE) --release $(ENS_VERSION) --outroot $(CDS_DIR) \
	--species_cache $(SPECIES_CACHE) --cds
	gunzip $(CDS_DIR)/*.gz

split_trees: # $(EMF_FILE)
	mkdir -p $(SEQSETS_ROOT)
	mkdir -p $(IMG_ROOT)
	mkdir -p $(TREE_ROOT)
	# $(BSUB) 
	python bin/split_trees.py --emf $(EMF_FILE) --clade $(CLADE) --treeroot $(TREE_ROOT) \
	--imgroot $(IMG_ROOT) --outroot $(SEQSETS_ROOT) --species_cache $(SPECIES_CACHE) \
	--species_list $(SPECIES_LIST)

split_trees_old: # $(EMF_FILE)
	mkdir -p $(SEQSETS_ROOT)
	mkdir -p $(IMG_ROOT)
	mkdir -p $(TREE_ROOT)
	# $(BSUB) 
	python bin/split_trees_old.py --emf $(EMF_FILE) --clade $(CLADE) --treeroot $(TREE_ROOT) \
	--imgroot $(IMG_ROOT) --outroot $(SEQSETS_ROOT) --species_cache $(SPECIES_CACHE)

# NOTE This used to be 'align_all'
make_seqsets: # split_trees $(wildcard $(SEQSETS_ROOT)/*/*)
	python bin/make_seqsets.py --clade $(CLADE) --cds $(CDS_DIR) --pep $(PEP_DIR) \
	--inroot $(SEQSETS_ROOT) --outroot $(SEQSETS_CDS_ROOT) --species_cache $(SPECIES_CACHE)

make_seqsets_pep: # split_trees $(wildcard $(SEQSETS_ROOT)/*/*)
	python bin/make_seqsets_pep.py --clade $(CLADE) --pep $(PEP_DIR) \
	--inroot $(SEQSETS_ROOT) --outroot $(SEQSETS_PEP_ROOT) --species_cache $(SPECIES_CACHE) \
	--specieslist data/specieslist.txt

# Nonessential
split_trees_nofilter: $(EMF_FILE)
	$(BSUB) python bin/split_trees.py --emf $(EMF_FILE) --clade $(CLADE) --treeroot $(TREE_ROOT) \
	--imgroot $(IMG_ROOT) --outroot $(SEQSETS_ROOT_UNFILTERED) --thr 0.0

# Nonessential
split_trees_diverge:
	python bin/split_trees_diverge.py --emf $(EMF_FILE) --clade $(CLADE) --treeroot $(TREE_ROOT) \
	--imgroot $(IMG_ROOT) --outroot $(DIVERGE_ROOT)

# Nonessential
align_all_diverge: # split_trees $(wildcard $(SEQSETS_ROOT)/*/*)
	python bin/align_all_pep.py --clade $(CLADE) --pep $(PEP_DIR) \
	--specieslist $(SPECIES_LIST) --inroot $(DIVERGE_ROOT) \
	--outroot $(DIVERGE_SEQSETS_ROOT)

# Nonessential
prank_diverge: 
	python bin/prank_bsub.py --clade $(CLADE) --inroot $(DIVERGE_SEQSETS_ROOT) \
	--treeroot $(DIVERGE_ROOT) --outroot $(DIVERGE_ALN_ROOT) \
	--logdir $(DIVERGE_PRANK_LOG_DIR) --mode protein

# Nonessential
rewrite_labels:
	python bin/rewrite_labels.py --clade $(CLADE) --treeroot $(DIVERGE_ROOT) \
	--alnroot $(DIVERGE_ALN_ROOT) --outroot $(HMMDIVERGE_ROOT)

# Nonessential
hmmdiverge:
	python bin/hmmdiverge_bsub.py --clade $(CLADE) --inroot $(HMMDIVERGE_ROOT) \
	--outroot $(HMMDIVERGE_ROOT) --logdir $(HMMDIVERGE_LOG_ROOT)

# Yeast analysis
# Make Fasta files using the genomes and YGOB 'pillars' as guide
# (Make alignments using PRANK)
# (Run RAxML)
parse_yeast:
	$(BSUB) python bin/yeast_make_seqsets.py 

align_yeast:
	python bin/prank_bsub.py --clade yeast --inroot $(PR_ROOT)/data/yeast/seqsets \
	--outroot $(PR_ROOT)/data/yeast/aln1 --logdir $(PR_ROOT)/log/yeast/prank1

reformat_yeast:
	$(BSUB) python bin/fasta2phylip.py --clade yeast \
	--inroot $(PR_ROOT)/data/yeast/aln1 --outroot $(PR_ROOT)/data/yeast/aln1

phylo_yeast:
	python bin/raxml_bsub.py --clade yeast --inroot $(PR_ROOT)/data/yeast/aln1 \
	--outroot $(PR_ROOT)/data/yeast/trees --logdir $(PR_ROOT)/log/yeast/raxml

realign_yeast:
	python bin/prank_bsub.py --clade yeast --inroot $(PR_ROOT)/data/yeast/seqsets \
	--treeroot $(PR_ROOT)/data/yeast/trees \
	--outroot $(PR_ROOT)/data/yeast/aln --logdir $(PR_ROOT)/log/yeast/prank

# Here the Ensembl and yeast analyses merge
# TODO Add dependencies
prank: 
	mkdir -p $(PRANK_LOG_DIR)
	python bin/prank_bsub.py --clade $(CLADE) --inroot $(SEQSETS_CDS_ROOT) --treeroot $(SEQSETS_ROOT) \
	--outroot $(ALN_ROOT) --logdir $(PRANK_LOG_DIR)

prepare_slr:
	$(BSUB) python bin/prepare_slr.py --clade $(CLADE) --treeroot $(SEQSETS_ROOT) \
	--alnroot $(ALN_ROOT) --outroot $(SLR_ROOT)
# $(SEQSETS_ROOT) \

slr:
	python bin/slr_bsub.py --clade $(CLADE) --slrroot $(SLR_ROOT) --logdir $(SLR_LOG_DIR)

process_slr:
	$(BSUB) -M 15000 -R "rusage[mem=15000]" python bin/process_slr.py --clade $(CLADE) --slrroot $(SLR_ROOT) --alnroot $(ALN_ROOT) \
	--outfile $(SLR_ALL)_tmp

all: align_all

$(SIFTS_GZIP):
	curl $(SIFTS_URL) > $(SIFTS_GZIP)

$(SIFTS_FILE): $(SIFTS_GZIP)
	gunzip -c $(SIFTS_GZIP) > $(SIFTS_FILE)

match_pdb: $(SIFTS_FILE)
	python bin/match_pdb.py --indir $(SLR_ROOT) --siftsfile $(SIFTS_FILE) --clade $(CLADE) --outfile $(PDB_MAP_FILE)

pdb_stats:
	python bin/pdb_stats.py --pep $(PEP_HUMAN) --mapfile $(PDB_MAP_FILE) --outfile $(PDB_STATS)

get_pdb:
	$(BSUB) python bin/get_pdb.py --pdbdir $(PDB_DIR) --pdbmap $(PDB_MAP_FILE) --log $(PDB_LOG) \
	--mismatch_log $(PDB_MISMATCH_LOG) --missing_log $(PDB_MISSING_LOG)

$(UNIPROT_SEQ_GZIP):
	curl $(UNIPROT_SEQ_URL) > $(UNIPROT_SEQ_GZIP)

$(UNIPROT_SEQ_FILE): $(UNIPROT_SEQ_GZIP)
	gunzip -c $(UNIPROT_SEQ_GZIP) > $(UNIPROT_SEQ_FILE)

pdb_master_table: $(UNIPROT_SEQ_FILE)
	$(BSUB) python bin/pdb_master_table.py --clade $(CLADE) --pdbmap $(PDB_MAP_FILE) --pdbdir $(PDB_DIR) \
	--slrroot $(SLR_ROOT) --outfile $(PDB_MASTER_TABLE)

# python bin/parse_sifts.py --clade Eutheria --pdbmap /nfs/research2/goldman/gregs/slr_pipeline/data/ens/78/pdb_map_Eutheria.tab --pdbdir /nfs/research2/goldman/gregs/slr_pipeline/data/pdb --slrroot /nfs/research2/goldman/gregs/slr_pipeline/data/ens/78/slr --outfile /nfs/research2/goldman/gregs/slr_pipeline/data/pdb_master_table_Eutheria.tab --siftsdir data/sifts

# bsub -I -M 20000 -R "rusage[mem=30000]" python bin/run_dssp.py --clade Eutheria --siftsmap /nfs/research2/goldman/gregs/slr_pipeline/data/pdb_master_table_Eutheria.tab --pdbdir /nfs/research2/goldman/gregs/slr_pipeline/data/pdb --slrroot /nfs/research2/goldman/gregs/slr_pipeline/data/ens/78/slr --outfile /nfs/research2/goldman/gregs/slr_pipeline/data/pdb_master_table_DSSP.tab --refgenome data/ens/78/pep/Homo_sapiens.GRCh38.pep.all.fa
