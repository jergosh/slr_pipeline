CLADE=Eutheria
# CLADE=Primates
# CLADE=yeast
SPECIES_LIST=data/specieslist.txt

PR_ROOT=`pwd`
ENS_ROOT=$(PR_ROOT)/data/ens/73
EMF_FILE=$(PR_ROOT)data/Compara.73.protein.nhx.emf
TREE_ROOT=$(ENS_ROOT)/trees
IMG_ROOT=$(ENS_ROOT)/img
SEQSETS_ROOT=$(ENS_ROOT)/seqsets
SEQSETS_CDS_ROOT=$(ENS_ROOT)/seqsets_cds
CDS_DIR=$(ENS_ROOT)/cds/canonical
PEP_DIR=$(ENS_ROOT)/pep
PEP_HUMAN=$(PEP_DIR)/Homo_sapiens.GRCh37.73.pep.all.fa

ALN_ROOT=$(ENS_ROOT)/aln
PRANK_LOG_DIR=$(PR_ROOT)/log/prank

SLR_ROOT=$(ENS_ROOT)/slr
SLR_LOG_DIR=$(PR_ROOT)/log/slr

SLR_ALL=$(ENS_ROOT)/$(CLADE)_all.tab

SIFTS_FILE=$(PR_ROOT)/data/sifts/pdb_chain_uniprot.tsv
PDB_MAP_FILE=$(PR_ROOT)/data/pdb_map_$(CLADE).tab
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


BSUB=bsub -I

split_trees: $(EMF_FILE)
	$(BSUB) python bin/split_trees.py --emf $(EMF_FILE) --clade $(CLADE) --treeroot $(TREE_ROOT) \
	--imgroot $(IMG_ROOT) --outroot $(SEQSETS_ROOT)


align_all: # split_trees $(wildcard $(SEQSETS_ROOT)/*/*)
	python bin/align_all.py --clade $(CLADE) --cds $(CDS_DIR) --pep $(PEP_DIR) \
	--specieslist $(SPECIES_LIST) --inroot $(SEQSETS_ROOT) \
	--outroot $(SEQSETS_CDS_ROOT)

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
prank: 
	python bin/prank_bsub.py --clade $(CLADE) --inroot $(SEQSETS_CDS_ROOT) --treeroot $(SEQSETS_ROOT) \
	--outroot $(ALN_ROOT) --logdir $(PRANK_LOG_DIR)

prepare_slr:
	$(BSUB) python bin/prepare_slr.py --clade $(CLADE) --treeroot $(TREE_ROOT) \
	--alnroot $(ALN_ROOT) --outroot $(SLR_ROOT)
# $(SEQSETS_ROOT) \

slr:
	python bin/slr_bsub.py --clade $(CLADE) --slrroot $(SLR_ROOT) --logdir $(SLR_LOG_DIR)

process_slr:
	$(BSUB) -M 15000 -R "rusage[mem=15000]" python bin/process_slr.py --clade $(CLADE) --slrroot $(SLR_ROOT) --alnroot $(ALN_ROOT) \
	--outfile $(SLR_ALL)

all: align_all

match_pdb:
	python bin/match_pdb.py --indir $(SLR_ROOT) --siftsfile $(SIFTS_FILE) --clade $(CLADE) --outfile $(PDB_MAP_FILE)

pdb_stats:
	python bin/pdb_stats.py --pep $(PEP_HUMAN) --mapfile $(PDB_MAP_FILE) --outfile $(PDB_STATS)

get_pdb:
	$(BSUB) python bin/get_pdb.py --pdbdir $(PDB_DIR) --pdbmap $(PDB_MAP_FILE) --log $(PDB_LOG) \
	--mismatch_log $(PDB_MISMATCH_LOG) --missing_log $(PDB_MISSING_LOG)

pdb_master_table:
	$(BSUB) python bin/pdb_master_table.py --clade $(CLADE) --pdbmap $(PDB_MAP_FILE) --pdbdir $(PDB_DIR) \
	--slrroot $(SLR_ROOT) --outfile $(PDB_MASTER_TABLE)
