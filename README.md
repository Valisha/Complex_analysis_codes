# Complex_analysis_codes
### Overview

These scripts represent some of the work I’m most proud of. They provide a streamlined framework for analyzing single-cell datasets, offering an immediate overview of differentially expressed genes (DEGs) across clusters when comparing tissues or conditions. The same logic enables flexible and intuitive clonotype-tracking visualizations, with options to highlight specific TCR sequences on demand.

### Why These Scripts Matter

Anyone who has worked in bioinformatics knows the familiar rhythm: biologists constantly need new figure aesthetics, updated comparisons, or custom subsetting—often at the last minute. These scripts were built as a “one-stop shop” to quickly generate high-quality summary plots so collaborators can immediately see the biological story. Once these initial visualizations guide the narrative, more targeted and publication-ready figures can be created based on the insights revealed.

### Within each script
1. dotplot_deg - this helps you create a dotplot for the top 25 degs for each comparison for each defined cluster. So it can help you compare the DEGs between tissue1 vs tissue2 within cluster Cd4 for example
2. top_clonotypes_alluvial - this helps you fined the shared clonotypes across different clusters/conditions and find the ones shared, along with that also "bold" the sequences that you are interested. Like I was interested in GP33-related sequences
3. editable html degs - this script will convert your DEG tabes into editable html outputs, so just add this block to you .Rmd, pass the DEG results here and the biologists will be able to use the rendered html outputs to edit the DEGs of their interest

### What You Can Do With This Toolkit
	•	Rapid DEG exploration for each cluster between tissue/condition groups
	•	Seamless clonotype tracking across clusters, samples, or phenotypic states
	•	Select and visualize specific TCR sequences as needed
	•	Produce clean, customizable figures that can be adapted for manuscripts 

