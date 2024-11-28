import pandas as pd
import dask.dataframe as dd

# Path to where the GRN results are stored
results_path="./pySCENIC_auxiliary/grn_results"

# List of all the files from the 100 runs
grn_files = [f"{results_path}/grn_results_AG_run_{i}.csv" for i in range(1, 101)]

threshold_runs = int(0.9 * 100)

ddf = dd.read_csv(grn_files)
aggregated_ddf = ddf.groupby(['TF', 'target']).agg(
  count=('TF', 'size'),
  importance_sum=('importance', 'sum')
).reset_index()

aggregated = aggregated_ddf.compute()
aggregated['avg_importance'] = aggregated['importance_sum'] / aggregated['count']
filtered = aggregated[aggregated['count'] >= threshold_runs]

filtered.to_csv("patchseq_filtered_grn_results.csv", index=False)
df=filtered[['TF','target','avg_importance']]
df['importance']=df['avg_importance']
df=df.drop(columns=['avg_importance'])
df.to_csv("patchseq_consensus_grn_results.csv", index=False)
