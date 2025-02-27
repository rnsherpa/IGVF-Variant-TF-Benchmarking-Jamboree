import pandas as pd
import numpy as np
import os

predictions = np.load('example_data/adastra_variant_scores.npy')

variants = pd.read_csv('example_data/adastra_variant_list.vcf',sep="\t")
tfs = pd.read_csv('example_data/adastra_tf_name_ensembl_ids.tsv',sep="\t")

targets_file = 'example_data/targets.txt'
df = pd.read_csv(targets_file,sep="\t")

desc = df['description'].to_numpy()

#and get chromosome sizes
chrom_sizes = {}
with open("example_data/hg38.chrom.sizes", "r") as f:
    for line in f:
        line = line.strip()
        if line:  # skip any empty lines
            chrom, size = line.split()
            chrom_sizes[chrom] = int(size)
# print(chrom_sizes)

context_window = 196608 // 2

folder_path = 'predictions/'
#create the folder if it does not exist
if not os.path.exists(folder_path):
    os.makedirs(folder_path)

for TF in tfs['hgnc_name']:
    TF_indices = [i for i,c in enumerate(desc) if TF in c.split(":")[1] and "CHIP:.:" not in c]
    print(f'TF: {TF}, number of cell types: {len(TF_indices)}')
    if len(TF_indices) == 0:
        print('No cell types found, skipping')
        continue
    
    #subset to the cell types that are relevant of shape (90172, cell_types, 3, 2)
    predictions_tf = predictions[:,TF_indices]
    
    #now get the 2 separate ones to calculate logfold change across cell types
    predictions_ref = predictions_tf[:,:,:,0]
    predictions_alt = predictions_tf[:,:,:,1]
    #now calculate the log2foldchange, alt over reference, and average over cell types, now shape: (90172,3)
    log2foldchange = np.log2(predictions_alt / predictions_ref).mean(axis=1)
    
    # Average predictions across the selected cell types; predictions_tf shape: (90172, 3, 2)
    predictions_tf = predictions_tf.mean(axis=1)
    
    # Bin 0: Primary scores
    ref_score = predictions_tf[:, 2, 0]
    alt_score = predictions_tf[:, 2, 1]
    variant_effect_score = log2foldchange[:,2]
    # Placeholder for pvalue (replace with your own computation if available)
    pvalue = np.full(ref_score.shape, 'N/A')
    
    # Bin 1: 2bins scores
    ref_score_2bins = predictions_tf[:, 0, 0]
    alt_score_2bins = predictions_tf[:, 0, 1]
    variant_effect_score_2bins = log2foldchange[:,0]
    
    # Bin 2: 16bins scores
    ref_score_16bins = predictions_tf[:, 1, 0]
    alt_score_16bins = predictions_tf[:, 1, 1]
    variant_effect_score_16bins = log2foldchange[:,1]
    
    
    
    # Create a new DataFrame starting from your variants dataframe.
    # Assuming 'variants' has columns: chr, pos, spdi, ref, alt, ref_seq_context, alt_seq_context
    df_out = variants.copy()
    
    #now we get the contexxt that was used by the model
    # Calculate the start and end positions while enforcing the limits:
    start = (df_out['pos'] - context_window).clip(lower=0)
    end = (df_out['pos'] + context_window).clip(upper=df_out['chr'].map(chrom_sizes))

    # Build the ref_seq_context string using these clipped values:
    df_out['ref_seq_context'] = (
        df_out['chr'].astype(str) + ':' +
        start.astype(str) + '-' +
        end.astype(str)
    )
    
    df_out['alt_seq_context'] = df_out['ref_seq_context'] #because we are using the same context for both ref and alt
    
    # Add prediction columns
    df_out['ref_score'] = ref_score
    df_out['alt_score'] = alt_score
    df_out['variant_effect_score'] = variant_effect_score
    df_out['pvalue'] = pvalue
    df_out['Enformer.ref_score_2bins'] = ref_score_2bins
    df_out['Enformer.alt_score_2bins'] = alt_score_2bins
    df_out['Enformer.variant_effect_score_2bins'] = variant_effect_score_2bins
    df_out['Enformer.ref_score_16bins'] = ref_score_16bins
    df_out['Enformer.alt_score_16bins'] = alt_score_16bins
    df_out['Enformer.variant_effect_score_16bins'] = variant_effect_score_16bins
    
    # Save the DataFrame to a TSV file with the naming scheme: Enformer_{TF}_ADASTRA_predictions.tsv
    filename = f"Enformer_{TF}_ADASTRA_predictions.tsv"
    # break
    df_out.to_csv(folder_path + filename, sep="\t", index=False)