# GsVec

GsVec（Gene signature Vector）は、Bioinformaticsの遺伝子発現解析等で得られたGene signature(複数の遺伝子名からなるリスト)に自然言語処理を用いてPathway・Gene Ontology等のGene signatureとの関連付けを行う事で、Gene signatureの生物学的解釈を支援する解析手法です。

従来のいわゆるPathway enrichment analysis（Fisher exact test, GSEA等）では、Gene signatureとのOverlapや分布の偏りと言った数値的な比較であり、ヒトが行うような含まれる遺伝子の特異性や重要度の比較を行う事は出来ませんでした。

GsVecはそれらの方法と異なり、自然言語処理における文章の分散表現の手法を取り入れ、遺伝子およびその集合である生物学的なGene signatureの特徴を抽出し、それを手元のsignatureと比較することで意味的な関連性を明らかにするという独自の手法となっています。

また、この手法はBioinformaticsの主流言語であるR言語のみで全工程を実行可能です。

詳細は下記の論文を参照してください。
> Yuumi Okuzono, Takashi Hoshino. "Comprehensive Biological Interpretation of Gene Signatures using Semantic Distributed Representation" 2019 (Manuscript in preparation)

----

## Download
- GsVec.tools.R
  - このフォルダの「GsVec.tools_v1.R」のファイルをダウンロードしてください。
  > [GsVec.tools_v1.R](https://github.com/yuumio/GsVec/blob/master/GsVec.tools_v1.R)
- Training data from MSigDB
  - **ライセンスはMSigDBの規約に従ってください。** http://software.broadinstitute.org/gsea/license_terms_list.jsp
  - 「C2: curated gene sets」の「CP: Canonical pathways」と、「C5: GO gene sets」の「BP: GO biological process」のGene Symbolで表記されたgmtファイルをダウンロードしてください。ただし、CPには一部ライセンスが必要なものがありますので、ダウンロード時は注意してください。
  > http://software.broadinstitute.org/gsea/downloads.jsp

----

## Install packages
GsVecはCRANの以下のパッケージを使用します。事前にインストールしておいてください。
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
以降のステップは全てR言語で行います。


### 1. load GsVec tools
まず、GsVecのコードをロードしてください。
~~~
source("./GsVec.tools_v05.R")
~~~

    
### 2. Preparate training and validation data
Training data from MSigDBのgmtファイルと、生物学的解釈を行いたい独自のGene signatureのgmtファイルを別のフォルダに格納し、以下の関数を実行します。
- フォルダ内には複数のgmtファイルを含めることが可能です。
- gmtファイルの形式は、MSigDBのgmtファイルを参考にしてください。
- "gene.number_min", "gene.number_max"により、何遺伝子以上・以下のSignatureを解析対象とするか、調節できます。（trainデータはデフォルトがおすすめです。）
~~~
train.data <- make_train.data(
  train.data_dir.path ="./train.set",
  gene.number_min = 10,
  gene.number_max = 500,
  export_name = "gsvec"
)
~~~
実行が終わると、current dir内に「"export_name"_train.data.txt」、「"export_name"_train.data_gene.txt」、「"export_name"_train.data_name.txt」が作成されます。
~~~
val.data <- make_validation.data(
  val.data_dir.path ="./val.set/",
  gene.number_min = 10,
  gene.number_max = 500,
  export_name = "deg"
)
~~~
実行が終わると、current dir内に「"export_name"_val.data.txt」、「"export_name"_val.data_gene.txt」、「"export_name"_val.data_name.txt」が作成されます。

> ### Tip:  
> 比較したいデータ（validation）データがgmt形式ではなく、遺伝子名とP-value、Fold change等が記載されている下記のようなリストであれば、"**fromat_val.data_from.deg**"の関数で、validation dataの形式に変換することができます。
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
まず、"**gs.train_genevec**"関数で、Gene vector(要はWord2Vec)を作ります。
- vector sizeやepoch数も変更可能です。（デフォルトがおすすめ）
- cpu.threadsはご自身の環境に合わせて変更してください。
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

