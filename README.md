## Repository for scripts and other data analysis protocols for *Aegilops tauschii* diversity project.

This paper has been published and can be found here https://www.frontiersin.org/articles/10.3389/fpls.2019.00009/full

1. Run `Tassel5 GBSv2` pipeline using the provided keyfile in `SNP_calling_scripts` folder to call SNPs.
2. `required_files` contain files with some required metadata.
3. Use the `R` scripts to perform a specific analysis.

	* `1.File.formatting.R` - format hapfile for further analyses
	* `2.SNP.Filtering.R` - filter SNPs based on criteria given in the paper
	* `3.PCA_Cluster.R` - plot PCA plots, biplots, and violin plots
	* `4.Genetic.diversity.R` - compute genetic diversity statistics
	* `5_Hybrid_graphical.R` - plot hybrids graphically
	
