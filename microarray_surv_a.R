library(affy) # to load CEL format
library(annotate) # to obtain gene names
library(hgu133a.db) # to obtain gene names
library(survival)

datadir = 'clinical_data'
file_clinical = 'GSE2034clinical.txt'
outputdir = 'output'

########################################################################
# read clinical datasets (GSE2034clinical)
# read CELL (expression) datasets (one CEL file for each sample)
#
print('reading datafiles')
# clinical data (286 samples in rows)
df_clin = read.table(file.path(datadir, file_clinical), header=FALSE)
colnames(df_clin) = c( 'pID ', 'filename ', 'nodes ', 'PFStime ', 'Progression ', 'ERstatus ', 'brainmets ')
# df_expr: expression data with shape: (22283, 286)
expr_data = justRMA(filenames=paste(df_clin$filename, ".CEL.gz", sep=""))
## attach gene symbols
ID = featureNames(expr_data)
gene_symbols = getSYMBOL(ID,"hgu133a.db")
fData(expr_data) = data.frame(ID=ID, Symbol=gene_symbols)
df_expr = exprs(expr_data)
print(dim(df_expr))

########################################################################
# survival analysis
# Cox proportional hazards method
#
survival_obj = Surv(df_clin$PFStime, df_clin$Progression)
#
## create a list of fitted models for each probe
cox_models = list()
for(i in 1:nrow(df_expr))
{
 cox_models[[i]] = coxph(survival_obj ~ df_expr[i,])
}
## get p-values and hazard ratios for each probe
## hazard ratio: positive value means that higher expression corresponds with an increased risk of disease progression
p_values = numeric(length(cox_models))
hazard_ratios = numeric(length(cox_models))
for(i in 1:length(cox_models))
{
	p_values[i] = summary(cox_models[[i]])$coefficients[5]
	hazard_ratios[i] = summary(cox_models[[i]])$coefficients[1]
}
p_adj_values = p.adjust(p_values,method="fdr")
p_ordered_row_idxs = order(p_adj_values, decreasing=FALSE)
p_ordered_IDs = featureNames(expr_data)[p_ordered_row_idxs]
p_ordered_genesymbols = fData(expr_data)$Symbol[p_ordered_row_idxs]
# collect results into a dataframe
df_result = cbind(fData(expr_data), hazard_ratios, p_values, p_adj_values)[p_ordered_row_idxs, ]
colnames(df_result) = c("ID", "Symbol", "HR", "p_value", "adj_p_value")
write.table(df_result, file=file.path(outputdir, "df_result.txt"), sep="\t", quote=FALSE, row.names=FALSE)
