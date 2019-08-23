# GsVec

GsVec（Gene signature Vector）は、Bioinformaticsの遺伝子発現解析等で得られたGene signature(複数の遺伝子名からなるリスト)に自然言語処理を用いてPathway・Gene OntologyのData baseの情報と手元のGene signatureの関連付けを行う事で、Gene signatureの生物学的解釈を支援する解析手法です。

従来のいわゆるPathway enrichment analysis（Fisher exact test, GSEA等）では、Gene signatureとのOverlapや分布の偏りと言った数値的な比較であり、ヒトが行うような含まれる遺伝子の特異性や重要度の比較を行う事は出来ませんでした。

GsVecはそれらの方法と異なり、自然言語処理における文章の分散表現の手法を取り入れ、遺伝子およびその集合である生物学的なGene signatureの特徴を抽出し、それを手元のsignatureと比較することで意味的な関連性を明らかにするという独自のGene signatureの解釈のための手法となっています。

詳細は下記のペーパーを参照してください。
> [Yuumi Okuzono, Takashi Hoshino. "" XXXXXXX p.- 2019](https://www.google.co.jp/)

-----

## Download
- GsVec.tools.R
このフォルダの「GsVec.tools_v5.R」のファイルをダウンロードしてください。
> [GsVec.tools_v5.R](https://github.com/yuumio/GsVec/GsVec.tool_v5.R)
-Training data from MSigDB
ライセンスはMSigDBの規約に従ってください。
CP

----
## Workflow
以降のステップは全てR言語で行います。

### 1. load GsVec tools
まず、GsVecのコードをロードしてください。
~~~
source("//XXX/XXX/GsVec.tools_v05.R")
~~~

### 2. Preparation training and test data
