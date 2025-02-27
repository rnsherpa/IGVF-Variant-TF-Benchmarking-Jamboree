# Enformer benchmarking Variant-TF IGVF Jamboree

This folder contains instructions on how to benchmark Enformer variants and format them in a way that works for the IGVF Jamboree, it includes all necessary files and examples

Currently, the model and the results have only been run for the adastra variant list and TFs, and instructions to replicate these results are what are outlined with the tutorial below

## Environment
Setting up the environment can be done with either Conda or some other python virtual environment. Python 3.8 is recommended

```
cd Enformer #make sure you're in the Enformer directory
### Activate your environment somehow ###
pip install -r requirements.txt
```

## Data
Example data is shown in the data folder, required files are
1. a variant list vcf file (example provided)
2. a list of TFs and their names (example provided)
3. HG38 fasta file (must be downloaded)
4. targets file specifying cell type for enformer track (provided)
5. list of chromosome sizes (provided)

Examples of these files can be seen in the data folder, only HG38 must be downloaded for the tutorial

## preparing the HG38 genome
The HG38 genome can be downloaded from here:

```
wget https://www.encodeproject.org/files/GRCh38_no_alt_analysis_set_GCA_000001405.15/@@download/GRCh38_no_alt_analysis_set_GCA_000001405.15.fasta.gz -O example_data/hg38.fa.gz

yes n | gunzip example_data/hg38.fa.gz
```

For speed, I precomputed the HG38 genome into a numpy array which is tokenized. I used a tokenized version of HG38 (where A->7, C->8, G->9, T->10, N->11). This is faster to load and to one hot encode for the large amount of sequences. Given a fasta file, the hg38 genome can be computed using

```python tokenize_HG38.py```

This takes a few minutes, all you need to do is give it the fasta file input and then the output will be the tokenized sequence saved by chromosome.

## running model evaluations
The first step is to run the model evaluations, this runs it for all the reference sequences and the alternate sequence with the snp. This will then save the outputs for reference and alternate sequence, and it does it by averaging across various bin levels, the surrounding 2, surrounding 16, and all bins. These bin levels were chosen to provide different resolutions to view the changes at, and make it consistent with its own and chrombpnet benchmarking. The results are saved to a np array of shape:

```len(vcf_file), 5313, 3, 2```

The length of the vcf file is the results per variant. The 5313 is the number of tracks enformer outputs, not all are needed, as there are many DNase or CAGE tracks that are not useful for this benchmarking task. the 3 are the 3 different resolutions, going from 2 bins to 16 bins and finally all 896 bins. The 2 is reference and alternate scores.

To run the model on all the sequences, the varaints list and the genome are needed. Simply modify the script to point to the right VCF file and the tokenized genome, and the output file path, then run:

```python score_TF_variants.py```

It will save your results in the defined output file path. Note that it may take several hours on even a large GPU. All tests were performed on an A100 GPU with 50GB of RAM. Running on ~90000 variants took ~14 hours. The output file size for this set was 20GB, but this can get larger as you get more variants to test.

One solution would be to not save it for all the cell types, or simply chunk it into smaller sets of variants.

## Formatting model outputs
The formatting shapes it in the specific format defined for the Jamboree. Note that this is cell type agnostic in that it will average across cell types. If there is only one cell type with the TF, it will still use that. The results are saved for each cell type individually, so the output could easily be made cell type specific.

It loops over the list of TFs, finds which cell tyypes have that TF, and then gets and averages the scores for all those cell types. Note that the log2fold change is put in the variant_effect_score, and is calculated per cell type before being averaged. It saves out the results per TF, you must choose an output dir, and also define where the predictions are, where the list of tfs is, the variant list, a targets file specifying what the enformer outputs relate to which cell types, and a file that specifies the chromosome sizes. To run the script, simply wait for the previous model evaluations to run and then run

```python format_TF_variants.py```

This will print out the number of cell types for each TF then saves each TF file separately (assuming that there is at least 1 cell type with that TF). This script is relatively quick and may take a few minutes.

### Column Descriptions

The columns are as follows:

- **chr**: Chromosome from the VCF file.
- **pos**: Position from the VCF file.
- **spdi**: Information from the VCF file.
- **ref**: Reference nucleotide from the VCF file (can be used for validation).
- **alt**: Alternate nucleotide from the VCF file.
- **ref_seq_context**: The reference sequence context, tells you the chromosome and which locations are read by Enformer. For Enformer, it is calculated as `nucleotide position - (196608/2)`, so that we get the full 196608 context length.  
  *Note: If the context extends beyond the chromosome boundaries, Enformer appends `N` nucleotides to maintain context length.*
- **alt_seq_context**: Identical to `ref_seq_context`.
- **ref_score**: The score averaged over all 896 bins for all cell types.
- **alt_score**: The score averaged over all 896 bins for all cell types after substituting the alternative nucleotide.
- **variant_effect_score**: The log₂ fold change calculated as  
  $$ \log_2 \left( \frac{\text{score}_{alt}}{\text{score}_{ref}} \right) $$  
  Positive values imply the alternate score is higher. The score is calculated per cell type and then averaged over all cell types.
- **pvalue**: Not applicable for this model.
- **Enformer.ref_score_2bins**: The score averaged over surrounding **2 bins** for all cell types.
- **Enformer.alt_score_2bins**: The score averaged over surrounding **2 bins** for all cell types after substituting the alternative nucleotide.
- **Enformer.variant_effect_score_2bins**: The log₂ fold change calculated as  
  $$ \log_2 \left( \frac{\text{score}_{alt}}{\text{score}_{ref}} \right) $$  
  Positive values imply the alternate score is higher. The score is calculated per cell type and then averaged over all cell types.
- **Enformer.ref_score_16bins**: The score averaged over surrounding **16 bins** for all cell types.
- **Enformer.alt_score_16bins**: The score averaged over surrounding **16 bins** for all cell types after substituting the alternative nucleotide.
- **Enformer.variant_effect_score_16bins**: The log₂ fold change calculated as  
  $$ \log_2 \left( \frac{\text{score}_{alt}}{\text{score}_{ref}} \right) $$  
  Positive values imply the alternate score is higher. The score is calculated per cell type and then averaged over all cell types.
