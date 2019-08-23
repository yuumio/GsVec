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
make_train.data(
  train.data_dir.path ="./train.set",
  gene.number_min = 10,
  gene.number_max = 500,
  export_name = "gsvec"
)

make_validation.data(
  val.data_dir.path ="./val.set/",
  gene.number_min = 10,
  gene.number_max = 500,
  export_name = "deg"
)
~~~

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

