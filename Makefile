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

BSUB=bsub -I

split_trees: $(EMF_FILE)
	$(BSUB) python bin/split_trees.py --emf $(EMF_FILE) --clade $(CLADE) --treeroot $(TREE_ROOT) \
	--imgroot $(IMG_ROOT) --outroot $(SEQSETS_ROOT)


align_all: # split_trees $(wildcard $(SEQSETS_ROOT)/*/*)
	python bin/align_all.py --clade $(CLADE) --cds $(CDS_DIR) --pep $(PEP_DIR) \
	--specieslist $(SPECIES_LIST) --inroot $(SEQSETS_ROOT) \
	--outroot $(SEQSETS_CDS_ROOT)

prank: 
	python bin/prank_bsub.py --clade $(CLADE) --inroot $(SEQSETS_CDS_ROOT) --treeroot $(TREE_ROOT) \
	--outroot $(ALN_ROOT) --logdir $(PRANK_LOG_DIR)

prepare_slr:

slr:

process_slr:

all: align_all
