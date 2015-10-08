import glob
import argparse
import pickle
from os import path
from Bio.Phylo.PAML import codeml

# Model A (alternative): model = 2, NSsites = 2,  fix_omega = 0 
# Model A1 (null): model = 2, NSsites = 2,  fix_omega = 1, omega = 1 


argparser = argparse.ArgumentParser()

argparser.add_argument('--indir', metavar="dir", type=str, required=True)
argparser.add_argument('--dataset', metavar="str", type=str, required=True)
argparser.add_argument('--sample', metavar="str", type=str, required=True)

def main():
    args = argparser.parse_args()
    
    aln = path.join(args.indir, "..", args.dataset+'.paml')
    tree = path.join(args.indir, args.dataset+'.nh')

    work_dir_1 = path.join(args.indir, '1')
    out_file_1 = path.join(work_dir_1, 'A1.out')
    cml_1 = codeml.Codeml(alignment=aln, tree=tree, out_file=out_file_1, working_dir=work_dir_1)
    cml_1.set_options(method=1, model=2, NSsites=[2], fix_omega=1, omega=1, verbose=1, seqtype=1, CodonFreq=2)
    results_1 = cml_1.run()
    pickle.dump(results_1, open(path.join(work_dir_1, "results_1.pk"), 'w'))

    work_dir_2 = path.join(args.indir, '2')
    out_file_2 = path.join(work_dir_2, 'A.out')
    cml_2 = codeml.Codeml(alignment=aln, tree=tree, out_file=out_file_2, working_dir=work_dir_2)
    cml_2.set_options(method=1, model=2, NSsites=[2], fix_omega=0, omega=1, verbose=1, seqtype=1, CodonFreq=2)
    results_2 = cml_2.run()
    pickle.dump(results_2, open(path.join(work_dir_2, "results_2.pk"), 'w'))

    LRT = 2*(results_2["NSsites"][2]["lnL"] - results_1["NSsites"][2]["lnL"])
    print LRT
    # find the sites if applicable & write them out
    # write out the lnL

if __name__ == "__main__":
    main()
