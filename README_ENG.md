# GsVec

GsVec (Gene signature Vector) is an analysis method that supports the biological interpretation of Gene signature obtained by gene expression analysis of Bioinformatics. The association between the gene signature to be interpreted and the gene signature of the Pathway / Gene Ontology data base is performed by natural language processing.

Conventional so-called Pathway enrichment analysis (Fisher exact test, GSEA, etc.) is a numerical comparison such as overlap with Gene signature and bias of distribution, and it was not possible to compare the specificity and importance of the genes involved as humans do.

GsVec is a unique method different from those methods. This method is a semantic method by taking a distributed representation method of sentences in natural language processing, extracting the characteristics of a gene and its biological gene signature, and comparing it with the gene signature to be interpreted to clarify the relevance.

In addition, this method can execute all processes only in the main language of Bioinformatics, R.

Please refer to the following paper for details.
> [Yuumi Okuzono, Takashi Hoshino. "" XXXXXXX p.- 2019](https://www.google.co.jp/)

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
-Vector size and number of epochs can be changed. (Default is recommended)
- Please change cpu.threads according to your environment.
~~~
train.fm <- gs.train_genevec(
  train_gene.txt="./gsvec_train.data_gene.txt", 
  train_vectorsize = 100, 
  train_epochs = 100, 
  cpu.threads =10
)
~~~
次に、"**estimate_cluster_size**"関数で、trainingデータに含まれるトピック（クラスター）数の目安をつけます。
- クラスタリングの方法はGMMとLDAから選択することができます。
- それ以外はデフォルトをおすすめします。
- このステップにはLDAでは数時間かかることがあります。
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
> - **mclustによるGMMのBIC（Bayesian information criterion）は、通常のBICと異なり（2*対数尤度 になっている）、最も大きい値が最適なクラスターです。GMMを用いる場合はBICが最も大きいクラスター数を選択してください。**
> - LDAのBICは通常通り、最も低いクラスター数を選択してください。

この値を用いて、Gene-topic vectorを"**gs.train_topicvec**"の関数で作成します。
- "cluster_n"に選択したクラスター数を入力してください。
- "cluster_method"は先程と同じものを選択してください（GMM or LDA）。
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
> 一度Gene-topic vectorを作れば、次回以降は1.のmake_val.dataと次のGsVecの実行だけです。
    
### 4. Run GsVec
"**GSVEC**"関数を用いて、ここまでに作成した"train.data","val.data","train.tv"(Gene-topic vector)を使って、Training dataとValidation dataのGene signature間の類似度を求めます。
- この結果を用いてtSNEで可視化するためには、"export_predict.gs.vector"のオプションをT(TRUE)にしておいてください。
- "calc.Fisher"のオプションT(TRUE)で、Fisher excat test(以降、Fisher)を同時に実行できます。Fisherの結果も並べて解釈することは有用です。
- このステップには、val.dataの数により、数時間かかる場合があります。
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


### Optional method: Vidualization by tSNE
Training dataとValidating dataをtSNEにより可視化をすることができます。
- gsvec.matには、4.の"export_predict.gs.vector"のオプションで作成された、"pred.val_freature.name.txt"を使用します。
- デフォルトではtSNEの前にPCAを行い、その95%以上の主成分でtSNEを行います（高速）。直接tSNEを行いたい場合は、"pca.thres"のオプションを"NA"にしてください。
- 他はデフォルトを推奨します。
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
> - デフォルトではTrain.dataとval.dataを2色に色分けします。
> - Signatureごとに色を分けたい場合は、IDの列にSignature nameを、Groupの列に色分けの対象となるグループ名を記載した、Data.frameを作成し、"gs.group_id.group.mat"のオプションで指定してください。
> 
> | ID | Group |
> |:---|:---| 
> | sig1 | group1 |
> | sig2 | group1 |
> | sig3 | group2 |
> | ...  | ... |

