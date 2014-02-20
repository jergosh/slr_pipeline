PR_ROOT=.
CLADE=Primates
SPECIES_LIST=data/specieslist.txt

ENS_ROOT=data/ens/73
EMF_FILE=data/Compara.73.protein.nhx.emf
TREE_ROOT=$(ENS_ROOT)/trees
IMG_ROOT=$(ENS_ROOT)/img
SEQSETS_ROOT=$(ENS_ROOT)/seqsets
SEQSETS_CDS_ROOT=$(ENS_ROOT)/seqsets_cds
CDS_DIR=$(ENS_ROOT)/cds
PEP_DIR=$(ENS_ROOT)/pep

BSUB=bsub -I

split_trees: $(EMF_FILE)
	$(BSUB) python bin/split_trees.py --emf $(EMF_FILE) --clade $(CLADE) --treeroot $(TREE_ROOT) \
	--imgroot $(IMG_ROOT) --outroot $(SEQSETS_ROOT)


align_all: # split_trees $(wildcard $(SEQSETS_ROOT)/*/*)
	python bin/align_all.py --cds $(CDS_DIR) --pep $(PEP_DIR) \
	--specieslist $(SPECIES_LIST) --inroot $(SEQSETS_ROOT) \
	--outroot $(SEQSETS_CDS_ROOT)

prank: 

prepare_slr:

slr:

process_slr:

all: align_all
