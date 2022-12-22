# Multi-omics-characterization-of-RNA-binding-proteins
## Multi-omics characterization of RNA binding proteins reveals disease comorbidities and potential drugs in COVID-19
Jiwei Pan1, Yueying Gao1, Huirui Han1, Tao Pan1, Jing Guo1, Si Li1, Juan Xu2, Yongsheng Li1

### Highlights:
	Transcriptomic and proteomic profiles of 178 COVID-19 patients were profiled and analyzed.
	RNA binding proteins were likely to be perturbed in infection and interacting with viral proteins.
	Interactome analysis revealed that RBPs were likely to locate in central regions of network.
	Network analysis revealed disease comorbidities and potential drugs in COVID-19.
![流程图](https://user-images.githubusercontent.com/91582097/209040232-31fa4674-7ad9-4db2-b32f-c9051a893e9c.png)

Scripts</br>
We shared the scripts we used for the analyses of the data.

Scripts: this folder contains several scripts used for the analyses.</br>
test: this folder contains adjusted scripts and test files for disease comorbidities and drug screening</br>

The script for disease comorbidities was run by the command line as follow:</br>
```python disease_distance.py --g1 RBP_SARS_id.txt```</br>
--g1 file containing EntreZ ID for the key genes interacting with SARS-Cov-2</br>
--g2 file containing EntreZ ID for diseases-related genes, this parameter defaults to the file in the disease folder</br>
-n file containing the network edgelist, this parameter defaults to 'interactome.tsv'</br>

The script for drug screening was run by the command line as follow:</br>
```Rscript drug_optimization.R drug_humanP.txt drugID.txt virus_mcode_select```</br>
drug_humanP.txt: file containing interaction between drugs and targets</br>
drugID.txt: file containing all drugs included in drug_humanP.txt</br>
virus_mcode_select: file containing SARS-Cov-2 related proteins in this folder</br>
