# GsVec

GsVec (Gene signature Vector) is an analysis method that supports the biological interpretation of Gene signature obtained by gene expression analysis of Bioinformatics. The association between the gene signature to be interpreted and the gene signature of the Pathway / Gene Ontology data base is performed by natural language processing.

Conventional so-called Pathway enrichment analysis (Fisher exact test, GSEA, etc.) is a numerical comparison such as overlap with Gene signature and bias of distribution, and it was not possible to compare the specificity and importance of the genes involved as humans do.

GsVec is a unique method different from those methods. This method is a semantic method by taking a distributed representation method of sentences in natural language processing, extracting the characteristics of a gene and its biological gene signature, and comparing it with the gene signature to be interpreted to clarify the relevance.

In addition, this method can execute all processes only in the main language of Bioinformatics, R.

Please refer to the following paper for details.
>*Yuumi Okuzono, Takashi Hoshino. "Comprehensive Biological Interpretation of Gene Signatures using Semantic Distributed Representation" 2019 (Manuscript in preparation)*

----

## Download
- GsVec.tools.R
  - Download the file “GsVec.tools_v1.R” in this directory.
  > [GsVec.tools_v1.R](GsVec.tools_v1.R)
- Training data from MSigDB
  - **Please follow the MSigDB rules for licensing.** http://software.broadinstitute.org/gsea/license_terms_list.jsp
  - Please download the gmt files indicated by Gene Symbol of “CP: Canonical pathways” of “C2: curated gene sets” and “BP: GO biological process” of “C5: GO gene sets”. However, some CPs require a license, so be careful when downloading.
  > http://software.broadinstitute.org/gsea/downloads.jsp

----

## Install packages
GsVec uses the following CRAN packages. Please install it in advance.
~~~
install.packages("data.table")
install.packages("ggplot2")
install.packages("reshape2")
install.packages("RColorBrewer")
install.packages("Rtsne")
install.packages("fastTextR")
install.packages("mclust")
install.packages("maptpx")
~~~

----
## Usage
All subsequent steps are done in R language.


### 1. load GsVec tools
First, load the GsVec codes.
~~~
source("./GsVec.tools_v05.R")
~~~

    
### 2. Preparate training and validation data
Store the gmt file(s) of Training data from MSigDB and the gmt file(s) of your own gene signature that you want to perform biological interpretation in **different folders** and execute the following functions.
- Multiple gmt files can be included in a folder.
- Please refer to the gmt file of MSigDB for the format of the gmt file.
- The gene signature to be analyzed can be adjusted to the number of genes more or less than the number of genes by the options of "gene.number_min", "gene.number_max". (The default is recommended for training data.)
~~~
train.data <- make_train.data(
  train.data_dir.path ="./train.set",
  gene.number_min = 10,
  gene.number_max = 500,
  export_name = "gsvec"
)
~~~
After execution, "'export_name'_train.data.txt", "'export_name'_train.data_gene.txt", and "'export_name'_train.data_name.txt" will be created in the current directory.
~~~
val.data <- make_validation.data(
  val.data_dir.path ="./val.set/",
  gene.number_min = 10,
  gene.number_max = 500,
  export_name = "deg"
)
~~~
After execution, "'export_name'_val.data.txt", "'export_name'_val.data_gene.txt", and "'export_name'_val.data_name.txt" will be created in the current directory.

> ### Tip:  
> If the data you want to compare (validation) is not in the gmt format, but is in the following list where the gene name, P-value, Fold change, etc. are described, the function "**fromat_val.data_from.deg**" Can be converted to validation data format.
> 
> |Gene|P-value|log2FC|
> |:------|------:|------:| 
> | Gene1 | 0.0015 |-1.53|
> | Gene2 | 0.0063 |4.32|
> | Gene3 | 0.0072 |1.34|
> | ...  | ... |...|
> 
> ~~~ 
> fromat_val.data_from.deg(
>   deg_txt = c("deg.txt"), #DEG file > 1, 
>   gene.symbol_col.no = 1, 
>   fold.change_col.no = 3,  # if No FC data = NA
>   fold.change_type = "log" # or "linear"
> )
> ~~~
    
### 3. Create Gene-topic vector of training data
First, create a Gene vector (in short, Word2Vec) with the "**gs.train_genevec**" function.
- Vector size and number of epochs can be changed. (Default is recommended)
- Please change cpu.threads according to your environment.
~~~
train.fm <- gs.train_genevec(
  train_gene.txt="./gsvec_train.data_gene.txt", 
  train_vectorsize = 100, 
  train_epochs = 100, 
  cpu.threads =10
)
~~~
Next, use the "**estimate_cluster_size**" function to give a rough estimate of the number of topics (clusters) included in the training data.
- The clustering method can be selected from GMM and LDA.
- Otherwise, the default is recommended.
- This step may take several hours for LDA.
~~~
estimate_cluster.n(
	genevec = train.fm,
	cluster_method = c("gmm","lda"), #"g","m","GMM","LDA"
	check_cluster.numbers = c(5,10,15,20,30,40,50,60,70,80,90,100),
	plot.score = T,
	plot.export_name = feature.name
)
~~~
> ### CAUTION!
> - **The BIC (Bayesian information criterion) of GMM by mclust is 2 * log likelihood unlike normal BIC, so the largest value is the optimal cluster. When using GMM, select the number of clusters with the largest BIC.**
> - Select the lowest number of clusters as usual for LDA BICs.

Use this value to create a Gene-topic vector with the function "**gs.train_topicvec**".
- Enter the number of clusters selected in "cluster_n".
- Select the same “cluster_method” as before (GMM or LDA).
~~~
train.tv <- gs.train_topicvec(
  gene.vec = train.fm, 
  train.data = train.data, 
  cluster_n = cluster_n, 
  cluster_method = c("gmm","lda"), #"g","m","GMM","LDA"
  save.cluster.data = T, 
  save.name = "feature.name")
out <- data.table(id = rownames(train.fm), train.fm, stringsAsFactors = F)
fwrite(out,paste0("train.fm_",feature.name,".txt"),sep="\t",quote=F,row.names=F)
~~~
> ### Note:
> Once you create a Gene-topic vector, you only need to execute 1. make_val.data and the next 4. run GsVec for the second and subsequent times.
    
### 4. Run GsVec
Execute "**GSVEC**" function on "train.data", "val.data", "train.tv" (Gene-topic vector) created so far, and find the similarity between Gene signature of Training data and Validation data.
- To visualize this result with tSNE, set the "export_predict.gs.vector" option to T (TRUE).
- If the "calc.Fisher" option is set to T (TRUE), the Fisher excat test can be executed simultaneously. It is useful to interpret Fisher results side by side.
- This step may take several hours depending on the number of val.data.
~~~
gsvec <- GSVEC(
  train.data = train.data,
  val.data = val.data,
  train.gsvec = train.tv,
  feature.name = "feature.name",
  export_predict.gs.vector = T,
  calc.Fisher = T #F
)
~~~


### (Optional method:) Vidualization by tSNE
Training data and Validating data can be visualized by tSNE.
- As gsvec.mat, use "pred.val_freature.name.txt" created with the "export_predict.gs.vector" option in 4.
- By default, PCA is performed before tSNE, and tSNE is performed with more than 95% of the principal components (high speed). If you want to perform tSNE directly, set the "pca.thres" option to "NA".
- Defaults are recommended for other options.
~~~
tmp <- as.data.frame(fread("./pred.val_feature.name.txt",stringsAsFactors = F))
pred.val <- as.matrix(tmp[,-1])
rownames(pred.val) <- tmp[,1]

pca.tsne_GsVec <- function(
  gsvec.mat = pred.val,	
  gsvec.train.data = train.data,
  gsvec.val.data = val.data,
  gs.group_id.group.mat = group.mat, 	#id,Group matrix
  pca.thres = 0.95, # or NA
  out_name = "pca.tsne",
  centering = F,
  tsne_max.iter = 500,
  tsne_sta = 200
)
~~~
> ### Tip:
> - By default, Train.data and val.data are color-coded into two colors.
> - If you want to separate colors for each signature, create a Data.frame with the signature name in the ID column and the group name to be color coded in the Group column, and specify "gs.group_id.group.mat" option.
> 
> | ID | Group |
> |:---|:---| 
> | sig1 | group1 |
> | sig2 | group1 |
> | sig3 | group2 |
> | ...  | ... |

