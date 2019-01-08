# Run RRHO on RNA seq data using brainSpan RNAseq

A script to compare _in-vitro_ RNAseq data from cell culture to _in-vivo_ RNAseq data from brainSpan modified from Stein, de la Torre-Ubieta et al. (2014).

The biggest difference between the published method in Stein, de la Torre-Ubieta et al. (2014) and this method is the choice of statistic to use for comparison. While the published mehtod above ranked genes by p-values this method uses logFC as it is more stable across different methods (e.g. DEseq2, edgeR, limma-voom)

The brainSpan data here only includes the cortical areas and was analyzed using limma-voom.

### Input
A folder containing CSVs of all the comparisons you are interested in in your data.
Each file must contain a column with ensembl gene IDs titled "geneID" and a column with the calculated log fold changes titled "logFC". The genes within each file don't have to be ordered by logFC.

### Output
A map of the overlap between any two systems across development.
An example of the map generated:
![alt text][results_example]

### Script

To run:

```bash
Rscript runRRHO path/to/csvs path/to/output
```

### Citation
If you use this analysis in a publication please cite the following [paper][Paper link]:

>Stein JL, de la Torre-Ubieta L, Tian Y, Parikshak NN, Hernandez IA, Marchetto MC, Baker DK, Lu D, Hinman CR, Lowe JK, Wexler EM, Muotri AR, Gage FH, Kosik KS, Geschwind DH (2014) A quantitative framework to evaluate modeling of cortical development by neural stem cells. Neuron 83:69-86.






[Paper link]: https://www.sciencedirect.com/science/article/pii/S0896627314004504?via%3Dihub
[results_example]: ../blob/master/picture1.png
