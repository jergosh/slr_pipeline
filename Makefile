PR_ROOT=.
CLADE=Primates
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

ALN_ROOT=$(ENS_ROOT)/aln
PRANK_LOG_DIR=$(PR_ROOT)/log/prank

SLR_ROOT=$(ENS_ROOT)/slr
SLR_LOG_DIR=$(PR_ROOT)/log/slr

SLR_ALL = $(ENS_ROOT)/$(CLADE)_all.tab

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
# Make alignments using trees from Kevin


# Here the Ensembl and yeast analyses merge
prank: 
	python bin/prank_bsub.py --clade $(CLADE) --inroot $(SEQSETS_CDS_ROOT) --treeroot $(SEQSETS_ROOT) \
	--outroot $(ALN_ROOT) --logdir $(PRANK_LOG_DIR)

prepare_slr:
	$(BSUB) python bin/prepare_slr.py --clade $(CLADE) --treeroot $(SEQSETS_ROOT) --alnroot $(ALN_ROOT) \
	--outroot $(SLR_ROOT)

slr:
	python bin/slr_bsub.py --clade $(CLADE) --slrroot $(SLR_ROOT) --logdir $(SLR_LOG_DIR)

process_slr:
	$(BSUB) python bin/process_slr.py --clade $(CLADE) --slrroot $(SLR_ROOT) --alnroot $(ALN_ROOT) \
	--outfile $(SLR_ALL)

all: align_all
