import graphtools as gt
import meld
import numpy as np
import pandas as pd
import argparse

parser = argparse.ArgumentParser()
# parser.add_argument("X_pca", help="csv file storing PCA data")
# parser.add_argument("coldata", help="csv file storing coldata")
parser.add_argument("seed", help = "Seed 4 batch effect")
parser.add_argument("pop", help = "Which cell type is DA?")
parser.add_argument("--pop_enrichment", dest='pop_enr', default=0.7,
                    help = "Max condition probability in DA population")
parser.add_argument("--batchEffect_sd", dest="be_sd", default='0',
                    help = "SD of batch effect")
parser.add_argument("--k", default=50,
                    help = "K parameter")
parser.add_argument("--data_id", default="embryo",
                    help = "ID for dataset")
parser.add_argument("--MNN_correct", default="no",
                    help = "use corrected PCA?")
args = parser.parse_args()

def run_meld(X_red_dim, sample_labels, conditions, k=15):
    '''
    Run MELD
    - X_red_dim: c x d matrix of dimensionality reduction to use for graph construction
    - sample_labels: assignment of cells to samples
    - conditions: vector of condition names
    '''
    ## Make graph
    graph = gt.Graph(X_red_dim, knn=int(k))
    ## Make MELD object
    meld_op = meld.MELD()
    meld_op.graph = graph
    ## Compute density
    meld_fit = meld_op.transform(sample_labels=np.array(sample_labels))
    
    ## Mean density per replicates
    mean_density = pd.DataFrame(
            np.zeros(shape=(meld_fit.shape[0], len(conditions))),
            index=meld_fit.index,
            columns=conditions,
        )
    
    for c in conditions:
      c_mean = meld_fit.loc[:,[c in x for x in meld_fit.columns]].mean(1)
      mean_density[c] = c_mean
    
    ## From density to likelihood per condition
    likelihoods = meld.utils.normalize_densities(mean_density)
    likelihoods.columns = [col.split("_")[0] for col in likelihoods.columns]
    return(likelihoods)
    # return(mean_density)

## Read input
outdir = '/nfs/team205/ed6/data/milo_benchmark/synthetic_data/'
outprefix = "benchmark_" + args.data_id + "_pop_" + args.pop + '_enr' + str(args.pop_enr) + "_seed" + str(args.seed)
if args.MNN_correct == "no":
  X_pca = pd.read_csv(outdir + outprefix + "_batchEffect" + str(args.be_sd) + ".pca.csv", index_col=0)
else:
  X_pca = pd.read_csv(outdir + outprefix + "_batchEffect" + str(args.be_sd) + ".MNNcorrected.pca.csv", index_col=0)
coldata = pd.read_csv(outdir + outprefix + ".coldata.csv", index_col=0)

## Run MELD
print("Running MELD...")
out = run_meld(X_pca, coldata["synth_samples"], ["Condition1", "Condition2"], k=args.k)
out_df = pd.DataFrame(out["Condition2"])
out_df.index = coldata.index
out_df["method"] = "meld"
out_df["true"] = coldata["true_labels"]
out_df["true_prob"] = coldata["Condition2_prob"]

## Save output
print("Saving MELD output...")
bm_outdir = '/nfs/team205/ed6/data/milo_benchmark/'
if args.MNN_correct == "no":
  out_df.to_csv(bm_outdir+outprefix+"_batchEffect"+str(args.be_sd)+".DAresults.meld.csv")
else:
  out_df.to_csv(bm_outdir+outprefix+"_batchEffect"+str(args.be_sd)+".MNNcorrected.DAresults.meld.csv")
