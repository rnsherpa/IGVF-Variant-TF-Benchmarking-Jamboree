# IGVF Variant-TF Benchmarking Jamboree

Resources and code for the IGVF Variant-TF Benchmarking jamboree. Hosted by the ML/AI Focus Group.

## Links
- [Planning doc](https://docs.google.com/document/d/1YXKNBrn4QEk1z6EpZZQyhAFzY1BjN3HZYWddFkJpYlc/edit?usp=sharing)
- [Data/model spreadsheet](https://docs.google.com/spreadsheets/d/12GV6qoquMtgg02OgVpT-_Su9ERSNltz3vPBX1IU7EoE/edit?usp=sharing)
- File formats:
	- [Variant-TF binding predictions](https://docs.google.com/document/d/1DLVO6eMUvyei3Qh3CI6JZ6CisxEi0WVLsZB9NQAYgKM/edit?usp=sharing)
	- [Variant-biochemical effect predictions](https://docs.google.com/document/d/1BtxbFQReBtysS4kO_CPpmdhvOvjVjYtjA5E_jLBexcA/edit?usp=sharing)
- [Jamboree synapse](https://www.synapse.org/Synapse:syn64024472)
	- [Variant lists](https://www.synapse.org/Synapse:syn64888527)
	- [TF lists](https://www.synapse.org/Synapse:syn64891521)
	- [Prediction outputs](https://www.synapse.org/Synapse:syn64862601)
- Synapse tutorials/docs
	- [Installing API clients](https://help.synapse.org/docs/Installing-Synapse-API-Clients.1985249668.html)
	- [Python docs](https://python-docs.synapse.org/)
- [Ontology Lookup Service](https://www.ebi.ac.uk/ols4/)

## Contribution guidlines
If you are working on model predictions, please branch this repository and commit any code/documentation into the corresponding model folder. This includes any pre-processing, running models, formatting outputs, etc.

Example:

```
git clone https://github.com/rnsherpa/IGVF-Variant-TF-Benchmarking-Jamboree.git

git checkout -b SEMpl

git commit -a -m "Code for running SEMpl on jamboree variants"

git push -u origin SEMpl 
```

When satisfactory, submit pull request to merge branch into main.
