10-02-2023

GWAS variants for IBD accessed via the GWAS catalog:
https://www.ebi.ac.uk/gwas/efotraits/EFO_0003767

vars_overlapping_testable_gene_windows.txt

Filtered to only keep variants which:
1.) are genome-wide significant (pval <= 1e-08)
2.) trait is UC, CD or IBD (no pairwise comparisons)
3.) map to a gene we have ENCODE coordinates for
4.) map to our ATAC peaks in cis-regulatory region of gene (n = 117 variants)


Sample data:

For each gene:
selected_peaks.csv: number of peaks selected for that model under "nPeaks"; the R2 is the discovery R2. not cross validation R2

Example:

Result,nPeaks,R2
ABI1_RFR_rf_ranker_results.txt,32,0.8434657279896888
ABI1_RFR_perm_ranker_results.txt,23,0.8394269554402887
ABI1_RFR_dropcol_ranker_results.txt,22,0.8443856713692368



ABI1_RFR_rf_ranker_results.txt --> go to the rf_ranker directory, get the list of ATAC peaks for 32 peaks: rf_ranker/ABI1_32peaks_rfranker_importance.csv
This rf_ranker/ABI1_32peaks_rfranker_importance.csv will have the list of peaks for the selected model:

Importance,Peak
0.17451643244429596,chr10:26860658-26861158
0.07077262512628227,chr10:26697432-26697932
0.05939505355671294,chr10:26859974-26860474
0.05930999080275433,chr10:26696907-26697407

Do any variants in the vars_overlapping_testable_gene_windows.txt file exist in the rf_ranker/ABI1_32peaks_rfranker_importance.csv file?

We also want to know how many overlap peaks we did not select for the model. For example, ABI1 32 peaks total; for ABI1_RFR_perm_ranker_results.txt, 23 peaks are selected. Are any vairants in the 9 unselected peaks?

May be worth checking by only choosing variants on the same chromosome --> some variants map be "mapped" to 1 gene in the GWAS catalog, but could be in cis-regulatory windows for other genes in our dataset.

Generate a output which records which variant(s) (if any) in vars_overlapping_testable_gene_windows.txt overlaps a peak, and whether that peak is 
If it occurs more than once because it maps to more than one gene, record both instances:

Suggested columns (one output per "Result/Test/Pipeline version"):

Gene:
SNP Coord:
ATAC Coord:
Selected ATAC peak (Y/N)?: 




