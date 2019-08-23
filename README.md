# GsVec

GsVec（Gene signature Vector）は、Bioinformaticsの遺伝子発現解析等で得られたGene signature(複数の遺伝子名からなるリスト)に自然言語処理を用いてPathway・Gene OntologyのData baseの情報と手元のGene signatureの関連付けを行う事で、Gene signatureの生物学的解釈を支援する解析手法です。

従来のいわゆるPathway enrichment analysis（Fisher exact test, GSEA等）では、Gene signatureとのOverlapや分布の偏りと言った数値的な比較であり、ヒトが行うような含まれる遺伝子の特異性や重要度の比較を行う事は出来ませんでした。

GsVecはそれらの方法と異なり、自然言語処理における文章の分散表現の手法を取り入れ、遺伝子およびその集合である生物学的なGene signatureの特徴を抽出し、それを手元のsignatureと比較することで意味的な関連性を明らかにするという独自のGene signatureの解釈のための手法となっています。

詳細は下記のペーパーを参照してください。
> [Yuumi Okuzono, Takashi Hoshino. "" XXXXXXX p.- 2019](https://www.google.co.jp/)

-----

## Download
- GsVec.tools.R
  - このフォルダの「GsVec.tools_v5.R」のファイルをダウンロードしてください。
  > [GsVec.tools_v5.R](https://github.com/yuumio/GsVec/GsVec.tool_v5.R)
- Training data from MSigDB
  - **ライセンスはMSigDBの規約に従ってください。** http://software.broadinstitute.org/gsea/license_terms_list.jsp
  - 「C2: curated gene sets」の「CP: Canonical pathways」と、「C5: GO gene sets」の「BP: GO biological process」のGene Symbolで表記されたgmtファイルをダウンロードしてください。ただし、CPには一部ライセンスが必要なものがありますので、ダウンロード時は注意してください。
  > http://software.broadinstitute.org/gsea/downloads.jsp

----
## Workflow
以降のステップは全てR言語で行います。

### 1. load GsVec tools
まず、GsVecのコードをロードしてください。
~~~
source("//XXX/XXX/GsVec.tools_v05.R")
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
> ~~~ 
> fromat_val.data_from.deg(
>   deg_txt = c("deg.txt"), #DEG file > 1, 
>   gene.symbol_col.no = 2, 
>   fold.change_col.no = 4,  # if No FC data = NA
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
- このステップにはLDAでは数時間以上かかることがあります。
~~~
estimate_cluster.n(
	genevec = train.fm,
	cluster_method = c("gmm","lda"), #"g","m","GMM","LDA"
	check_cluster.numbers = c(1,5,10,15,20,30,40,50,60,70,80,90,100),
	plot.score = T,
	plot.export_name = feature.name
)
~~~
BIC（Bayesian information criterion）または、BF（Baysian Factor）が最も高いトピック数を選択してください。

この値を用いて、Gene-topic vectorを"**gs.train_topicvec**"の関数で作成します。
- "cluster_n"に選択したクラスター数を入力してください。
- "cluster_method"は先の"estimate_cluster_size"で指定したものを入力してください。
~~~
train.fm <- gs.train_topicvec(
  gene.vec = train.fm, 
  train.data = train.data, 
  cluster_n = cluster_n, 
  cluster_method = c("gmm","lda"), #"g","m","GMM","LDA"
  save.cluster.data = T, 
  save.name = "feature.name")
out <- data.table(id = rownames(train.fm), train.fm, stringsAsFactors = F)
fwrite(out,paste0("train.fm_",feature.name,".txt"),sep="\t",quote=F,row.names=F)
~~~

