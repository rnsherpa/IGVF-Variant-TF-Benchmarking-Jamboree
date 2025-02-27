#this script will score variants (reference and alternate) for many TFs, it will actually score all tracks, but the outputs will be saved
#the values wil be averaged across bins, so we will do one that is 2 bins, one that is 16 bins, and one that is all bins, and see which is most informative
#this is only for adastra set!

import numpy as np
import pandas as pd
import tensorflow as tf
from tqdm import tqdm
import tensorflow_hub as hub

genome_np = 'example_data/hg38_tokenized.npz'
with np.load(genome_np) as data:
    genome = {key: np.array(data[key]) for key in data}

variants_path = 'example_data/adastra_variant_list.vcf'
variants = pd.read_csv(variants_path,sep="\t")

mapping = {'A':(1,0,0,0),'C':(0,1,0,0),'G':(0,0,1,0),'T':(0,0,0,1),'N':(0,0,0,0)}

length = 196608
output_array = np.zeros((variants.shape[0],5313,3,2))

enformer = hub.load("https://kaggle.com/models/deepmind/enformer/frameworks/TensorFlow2/variations/enformer/versions/1").model

for i in tqdm(range(variants.shape[0])):
    #initialize values and load in data
    padleft = 0
    padright = 0
    chrom = variants['chr'][i]
    pos = variants['pos'][i]-1
    start = pos - length//2
    end = pos + length//2

    #if the start or end is out of bounds, we need to pad
    if start < 0:
        padleft = -start
        start = 0
    elif end > len(genome[chrom]):
        padright = end - len(genome[chrom])
        end = len(genome[chrom])
    
    #get values and pad
    seq = genome[chrom][start:end]
    if padleft:
        #now we pad the left and eight by the amounts
        seq_pad = np.ones(length, dtype=int)*11
        seq_pad[padleft:padleft+seq.size] = seq
        seq = seq_pad
    elif padright:
        seq_pad = np.ones(length, dtype=int)*11
        seq_pad[:seq.size] = seq
        seq = seq_pad

    #now we convert to one hot
    seq_convert = (seq-7)%4
    seq_one_hot = np.zeros((seq.size,4))
    seq_one_hot[np.arange(seq.size),seq_convert] = 1
    if np.any(seq == 11):
        seq_one_hot[seq==11,:] = 0
    
    #now we need to create the full input for the model and get the alternate
    full_input = np.zeros((393216,4))
    start_index = (393216 - 196608) // 2  # Equals 98304
    end_index = start_index + 196608
    full_input[start_index:end_index,:] = seq_one_hot
    full_input_altered = full_input.copy()
    full_input_altered[393216//2] = mapping[variants['alt'][i]]

    #now combine together and predict
    full_input_combined = np.stack([full_input,full_input_altered],axis=0)
    predictions = enformer.predict_on_batch(full_input_combined)

    #output is shape (2,896,5313)
    outputs = predictions['human'].numpy()
    val_2bins = outputs[:,447:449,:].mean(axis=1).T
    val_16bins = outputs[:,440:456,:].mean(axis=1).T
    val_full = outputs.mean(axis=1).T
    #now we save out the results, each of these is now 5313x2
    output_array[i,:,0] = val_2bins
    output_array[i,:,1] = val_16bins
    output_array[i,:,2] = val_full

#now save out our output array
np.save('example_data/adastra_variant_scores.npy',output_array)