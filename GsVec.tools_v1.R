#------------------------------------------------------------------------------------------------------ Require packages

	library("data.table")
	library("ggplot2")
	library("reshape2")
	library("RColorBrewer")
	library("Rtsne")
	library("fastTextR")
	library("mclust")
	library("maptpx")


## ---------------------------------------------------------------------------------------------------------------------------------
## -------------------------------------------------------------------------------------------------- format train/val.data from gmt

format_gs.data_from.gmt <- function(x, population = F, gene.population = NA){
		tmp <- strsplit(x,"\t")[[1]]
		tmp2 <- unique(tmp[3:length(tmp)])
		if(population == T){
			tmp2 <- tmp2[tmp2 %in% gene.population]
		}
		out <- c(tmp[1],tmp[2],length(tmp2),paste(tmp2[order(tmp2)],collapse=" "))
		return(out)
	}

## -------------------------------------------------------------------------------------------------- make_train.data from gmt

make_train.data <- function(
 	train.data_dir.path ="./",
 	gene.number_min = 10,
 	gene.number_max = 500,
	export_name = "gsvec"
){
	train.files <- paste0(gsub("\\/$","",train.data_dir.path),"/",dir(train.data_dir.path,pattern=".*\\.gmt$"))
	for(i in 1:length(train.files)){
		tmp.td <- as.character(as.matrix(fread(train.files[i],stringsAsFactors = F,sep=",",header=F)[,1]))
		tmp.td <- t(sapply(tmp.td,format_gs.data_from.gmt))
		if(i == 1){
			train.data <- tmp.td
		}else{
			train.data <- rbind(train.data,tmp.td)
		}
	}
	colnames(train.data) <- c("name","detail","number","gene")
	train.data <- as.data.frame(train.data,stringsAsFactors = F)
	
	train.data$number <- as.numeric(train.data$number)
	train.data <- train.data[train.data$number >= gene.number_min & train.data$number <= gene.number_max,]
	train.data <- data.frame(id = paste0("t",seq(1,nrow(train.data))),train.data,stringsAsFactors = F)

	fwrite(as.data.table(train.data$gene),paste0(export_name,"_train.data_gene.txt"),sep="\t",quote=F,row.names=F,col.names=F)
	fwrite(as.data.table(train.data$name),paste0(export_name,"_train.data_name.txt"),sep="\t",quote=F,row.names=F,col.names=F)
	fwrite(as.data.table(train.data),paste0(export_name,"_train.data.txt"),sep="\t",quote=F,row.names=F,col.names=T)

	return(train.data)
}

## -------------------------------------------------------------------------------------------------- make_validation.data from.gmt

make_validation.data <- function(
 	val.data_dir.path ="./",
 	gene.number_min = NA,
 	gene.number_max = NA,
	export_name = "gsvec"
){
	train.files <- paste0(gsub("\\/$","",val.data_dir.path),"/",dir(val.data_dir.path,pattern=".*\\.gmt$"))
	for(i in 1:length(train.files)){
		tmp.td <- as.character(as.matrix(fread(train.files[i],stringsAsFactors = F,sep=",",header=F)[,1]))
		tmp.td <- t(sapply(tmp.td,format_gs.data_from.gmt))
		if(i == 1){
			val.data <- tmp.td
		}else{
			val.data <- rbind(val.data,tmp.td)
		}
	}
	colnames(val.data) <- c("name","detail","number","gene")
	val.data <- as.data.frame(val.data,stringsAsFactors = F)
	
	val.data$number <- as.numeric(val.data$number)
	if(is.na(gene.number_min)!=T){
		val.data <- val.data[val.data$number >= gene.number_min,]
	}
	if(is.na(gene.number_max)!=T){
		val.data <- val.data[val.data$number <= gene.number_max,]
	}
	val.data <- data.frame(id = paste0("v",seq(1,nrow(val.data))),val.data,stringsAsFactors = F)

	fwrite(as.data.table(val.data$gene),paste0(export_name,"_val.data_gene.txt"),sep="\t",quote=F,row.names=F,col.names=F)
	fwrite(as.data.table(val.data$name),paste0(export_name,"_val.data_name.txt"),sep="\t",quote=F,row.names=F,col.names=F)
	fwrite(as.data.table(val.data),paste0(export_name,"_val.data.txt"),sep="\t",quote=F,row.names=F,col.names=T)

	return(val.data)
}


#------------------------------------------------------------------------------------------------------ IDF

gs.train_idf <- function(
	train.data = train.data
){

	gtm <- itoken(as.character(as.matrix(train.data$gene)), preprocess_function = toupper, tokenizer = space_tokenizer, ids = train.data$id)
	genes <- create_vocabulary(gtm)
	genes <- genes[order(genes$term),]
	vectorizer <- vocab_vectorizer(genes)
	gtm <- create_dtm(gtm, vectorizer)

	train_gtm <- t(as.matrix(gtm))
	train_n <- ncol(gtm)
	train_idf <- data.frame(genes[,c("term","doc_count")], idf =log((train_n + 1) / (genes$doc_count + 1)) + 1)
	colnames(train_idf) <- c("gene","train_count","idf")

	return(train_idf)
}

## -------------------------------------------------------------------------------------------------- GV by fastText

gs.train_genevec <- function(
	train_gene.txt = "train.data_gene.txt",
	train_vectorsize = 100,
	train_epochs = 50,
	train_learning.rate =0.05,
	cpu.threads = 20
){
	cntrl <- ft.control(
		loss = "ns",
		word_vec_size = train_vectorsize, 
		learning_rate = 0.05, 	#0.025
		max_len_ngram = 10000,  
		min_count = 1, 
		min_ngram = 10000,
		max_ngram = 10000,
		nbuckets = 2000000L, 
		epoch = train_epochs, 			#100
		nthreads = cpu.threads,
		verbose=2)

	model <- fasttext(input = train_gene.txt, method = "skipgram", control = cntrl)

	tmp <- as.data.frame(fread(".vec",sep=" ",header=F,skip=1))
	train.fm <- tmp[-1,-1]
	rownames(train.fm)<- tmp[-1,1]
	train.fm <- train.fm[order(rownames(train.fm)),]

	return(train.fm)
}

## -------------------------------------------------------------------------------------------------- estimate_cluster_size

estimate_cluster.n <- function(
	genevec = train.fm,
	cluster_method = c("gmm","lda"), #"g","m","GMM","LDA"
	check_cluster.numbers = c(5,10,15,20,30,40,50,60,70,80,90,100),
	gmm.modelNames = "EEI",
	plot.score = T,
	plot.export_name = feature.name
){
	if(tolower(substr(cluster_method,1,1)) == "g"){
		out <- mclustBIC(train.fm, G=check_cluster.numbers, modelNames ="EEI")
	}
	if(tolower(substr(cluster_method,1,1)) == "l"){
		data = as.simple_triplet_matrix(t(train.fm+abs(min(train.fm))))
		out <- list()
		cat("k = ")
		for(i in 1:length(check_cluster.numbers)){
			cat(paste0(i," "))
			out[[i]] <- topics(data, K = check_cluster.numbers[i], bf=T)
		}
		save(out,file=paste0("estimate_cluster_LDA.output_",feature.name,".RData"))
	}
	if(plot.score == TRUE){
		if(tolower(substr(cluster_method,1,1)) == "g"){
			png(paste0("estimate_cluster.n__",feature.name,".png"), width = 800, height = 600)
			plot(out)
			dev.off()
		}
		if(tolower(substr(cluster_method,1,1)) == "l"){
			BIC <- rep(0.0,length(check_cluster.numbers))
			names(BIC) <- check_cluster.numbers
			T <- dim(out[[1]]$theta)[1][1]
			for(i in 1:length(out)){
				BIC[i] <- T*log(out[[i]]$D[[1]]^2)+check_cluster.numbers[i]*log(T)
			}
			png(paste0("estimate_cluster.n__",feature.name,".png"), width = 800, height = 500)
				matplot(names(BIC),BIC,type="b",pch=19,xlab="Number of components",xlim = range(check_cluster.numbers),xaxt = "n")
				axis(side = 1, at = check_cluster.numbers)
			dev.off()
		}
	}
	return(out)
}

#------------------------------------------------------------------------------------------------------ GeneVec2GsVec

predict_GsVec_from.GeneVec <- function(	
	train.gene.feature.mat = train_gene.mat,	# r: gene, c:feature
	train.data = train.data,
	val.data = val.data
){
	
	val.data <- rbind(train.data,val.data)
	feature_sum <- function(val.sig, train.gene.feature.mat){
		return(colMeans(train.gene.feature.mat[rownames(train.gene.feature.mat) %in% unique(strsplit(val.sig," ")[[1]]),]))
	}
	pred.vector <- t(sapply(val.data$gene,function(x) feature_sum(x, train.gene.feature.mat)))
	rownames(pred.vector) <- val.data$id

	return(pred.vector)
}

## -------------------------------------------------------------------------------------------------- Topic vector by GMM or LDA

gs.train_topicvec <- function(
	gene.vec = train.fm,
	train.data = train.data,
	cluster_n = 60,
	cluster_method = c("gmm","lda"), #"g","m","GMM","LDA"
	gmm.modelNames = "EEI",
	lda.model = train.cluster, #or NA or ""
	save.cluster.data = F,
	save.name = feature.name,
	lda.model_RData.path = F

){
	# IDF
	train.idf <- gs.train_idf(train.data)

	# topic
	if(tolower(substr(cluster_method,1,1)) == "g"){
		train.cluster <- Mclust(gene.vec, G = cluster_n, modelNames = gmm.modelNames)
		if(save.cluster.data == T){
			save(train.cluster,file=paste0("train.cluster_",feature.name,".RData"))
		}
		gene.vec_cluster.p <- train.cluster$z
	}
	if(tolower(substr(cluster_method,1,1)) == "l"){
		if(is.na(lda.model_RData.path)== T | lda.model_RData.path!=""){
			data = as.simple_triplet_matrix(t(gene.vec+abs(min(gene.vec))))
			train.cluster <- topics(data, K = cluster_n, bf=T)
			if(save.cluster.data == T){
				save(train.cluster,file=paste0("train.cluster_",feature.name,".RData"))
			}
		}
		gene.vec_cluster.p <- train.cluster$theta
	}
	gene.vec_n <- nrow(gene.vec)
	gene.vec_dim <- ncol(gene.vec)

	# Gene-topic vector
	train.topicvec <- matrix(numeric(), nrow = gene.vec_n, ncol= gene.vec_dim * cluster_n)
	rownames(train.topicvec) <- rownames(gene.vec)

#	for (i in 0:(cluster_n-1)) {
#	    train.topicvec [, (i * gene.vec_dim + 1):((i + 1) * gene.vec_dim)] <- gene.vec * gene.vec_cluster.p[, i + 1]
#	}
	for (i in 1:(cluster_n)) {
	    if(i == 1){
	        train.topicvec <- gene.vec * gene.vec_cluster.p[,i]
	    }else{
	        train.topicvec <- cbind(train.topicvec,gene.vec * gene.vec_cluster.p[,i])
	    }
		cat(i)
	}
	train.topicvec <- train.topicvec * train.idf$idf

	#
	return(train.topicvec)

}

## -------------------------------------------------------------------------------------------------- predict_GsVec_from.TopicVec

predict_GsVec_from.TopicVec <- function(
	train.topicvec = train.topicvec,
	train.data = train.data,
	val.data = val.data,
	sparse_thres.p = 0	# or 0.04
){
	val.data <- rbind(train.data,val.data)

	# GS vector
	gsvec <- matrix(numeric(), nrow = nrow(val.data), ncol= ncol(train.topicvec))
	rownames(gsvec) <- val.data$id
	for(i in 1:nrow(val.data)){
		tmp <- colSums(train.topicvec[rownames(train.topicvec) %in% unique(strsplit(val.data$gene[i]," ")[[1]]),])
		norm <- sqrt(sum(tmp**2))
		if(norm != 0){
			tmp <- tmp/norm
		}
		gsvec[i,] <- tmp
	}

	# Sparse Composition
	gsvec_n <- nrow(gsvec)
	a_min <- sum(apply(gsvec,1,min))*1/gsvec_n
	a_max <- sum(apply(gsvec,1,max))*1/gsvec_n
	th <- (abs(a_min)+abs(a_max))/2 * sparse_thres.p

	gsvec[abs(gsvec) < th] <- 0

	#
	return(gsvec)
}
## -------------------------------------------------------------------------------------------------- fisher exact test

gs.enrich_fisher <- function(
	train.data = train.data,
	val.data = val.data,
	use_gene.population.txt = NA,	### if use_gene.population.txt = NA, make population.data from train.data and validation data
	export_fisher.result = F,
	export_detail = F,
	export_name = export_name
){

	# prep train.data
	colnames(train.data) <- c("tid","train.data_name","train.data_detail","train.data_number","train.data_gene")
	rownames(train.data) <- train.data$tid

	# prep validation.data
	rownames(val.data) <- val.data$vid
	colnames(val.data) <- c("vid","val.data_name","val.data_detail","val.data_number","val.data_gene")

	# set gene population for Fisher exact test
	if(!is.na(use_gene.population.txt) == T){
		gene.pop <- as.character(as.matrix(fread(use_gene.population.txt, stringsAsFactors = F,sep="\t",header=F)))
	}else{
		gene.pop <- unique(strsplit(paste(train.data$train.data_gene,val.data$val.data_gene,collapse=" ")," ")[[1]])
	}
	gene.pop <- gene.pop[order(gene.pop)]

	# make matrix
	mat_ts <- sapply(train.data$train.data_gene,function(x) gene.pop %in% strsplit(x," ")[[1]])
	colnames(mat_ts) <- train.data$tid
	mat_vs <- sapply(val.data$val.data_gene,function(x) gene.pop %in% strsplit(x," ")[[1]])
	colnames(mat_vs) <- val.data$vid

	# calc fisher
	format_fisher <- function(gset,tset,gene.pop){
		indx <- gset & tset
		ov <- sum(indx)
		return(c(
			round(-log10(fisher.test(matrix(c(ov, sum(gset)-ov, sum(tset)-ov, length(gene.pop)-(sum(gset)+sum(tset)-ov)),ncol=2), alternative = "greater")$p.value),3),
			ov,
			paste(gene.pop[indx],collapse=",")
		))
	}
	cat(paste0("val.data (", ncol(mat_vs),") : "))
	for(i in 1:ncol(mat_vs)){
		tmp <- t(apply(mat_ts,2, function(y) format_fisher(mat_vs[,i],y,gene.pop)))
		tmp <- cbind(paste0(colnames(mat_vs)[i],rownames(tmp)),colnames(mat_vs)[i],rownames(tmp),tmp)
		if(i==1){
			out <- tmp
		}else{
			out <- rbind(out,tmp)
		}
		cat(paste0(i," "))
	}
	colnames(out) <- c("vtid","vid","tid","Fisher_P.-log10","Overlap_number","Overlap_gene")

	# export_fisher.result
	if(export_fisher.result == T){
		out <- merge(merge(out,val.data[,c("vid","val.data_name","val.data_detail","val.data_number")],by="vid"),
				train.data[,c("tid","train.data_name","train.data_detail","train.data_number")],by="tid")
		if(export_detail == T){
			indx <- c("vtid","vid","tid","val.data_name","val.data_detail","val.data_number","train.data_name","train.data_detail","train.data_number",
					"Fisher_P.-log10", "Overlap_number", "Overlap_gene")
		}else{
			indx <- c("vtid","vid","tid","val.data_name","val.data_number","train.data_name","train.data_number",
					"Fisher_P.-log10", "Overlap_number", "Overlap_gene")
		}
		fwrite(out[,indx], paste0("Fisher.result__",export_name,".txt"),sep="\t",quote=F,row.names=F)
	}
	return(out)
}

## -------------------------------------------------------------------------------------------------- similarity_vectors by cosine dist.

similarity_vectors <- function(
	train_feature.matrix = train.mat,
	val_feature.matrix = val.mat,
	feature.name = "gsvec",
	convert.by.melt=T
){

	numer <- as.matrix(train_feature.matrix) %*% t(as.matrix(val_feature.matrix))
	denom1 <- apply(train_feature.matrix,1,function(x) sqrt(sum(x^2)))
	denom2 <- apply(val_feature.matrix,1,function(x) sqrt(sum(x^2)))
	denom <- t(sapply(denom1,function(x) sapply(denom2,function(y) x*y)))
	dist.mat <- t(numer/denom)

	if(convert.by.melt == T){
		out <- melt(dist.mat,factorsAsStrings=F, value.name = feature.name)
		out$vtid <- paste0(out[,1],out[,2])
		out$vid <- as.character(out[,1])
		out$tid <- as.character(out[,2])
		return(out[,c(4,5,6,3)])
	}else{
		return(dist.mat)
	}
}

## -------------------------------------------------------------------------------------------------- fromat_val.data_from.deg

fromat_val.data_from.deg <- function(
  deg_txt = c("deg.txt"), #DEG file > 1,
  gene.symbol_col.no = 2,
  fold.change_col.no = 4,  # if No FC data = NA
  fold.change_type = "log" # or "linear"
){
  out <- data.frame(name=gsub("","",deg_txt),detail=deg_txt,number=numeric(),gene=character())
  if(is.na(fold.change_type)==F){
    out_up <- out
    out_up$name <- paste0(out_up$name,"__up")
    out_down <- out
    out_down$name <- paste0(out_up$name,"__down")
  }
  for(i in 1:length(deg_txt)){
    data <- fread(deg_txt[i],stringsAsFactors = F,sep="\t",header=T)
    genes <- unique(data[,gene.sybol_col.no])
    genes <- genes[genes[!is.na(genes)] & genes !=""]
    out[i,3] <- length(genes)
    out[i,4] <- genes
    if(is.na(fold.change_col.no)==T){
      next
    }
    if(fold.change_type == "linear"){
      data[,2] <- log2(data[,fold.change_col.no])
    }
    indx <- data[,fold.change_col.no] >0
    if(sum(indx)!=0){
      genes <- unique(data[indx,gene.sybol_col.no])
      genes <- genes[genes[!is.na(genes)] & genes !=""]
      out_up[i,3] <- length(genes)
      out_up[i,4] <- genes
    }
    if(sum(!indx)!=0){
      genes <- unique(data[indx,gene.sybol_col.no])
      genes <- genes[genes[!is.na(genes)] & genes !=""]
      out_down[i,3] <- length(genes)
      out_down[i,4] <- genes
    }
    if(i == length(deg_txt)){
      out <- rbind(out,out_up,out_down)
    }
    cat(i)
  }
  val.data <- data.frame(paste0("v",seq(1,nrow(out))),out,stringsAsFactors = F)
  colnames(val.data) <- c("id","name","detail","number","gene") 
  return(val.data)
}

## -------------------------------------------------------------------------------------------------- GSVEC (main function; Predict + Fisher)

GSVEC <- function(
  train.data = train.data,
  val.data = val.data,
  train.gsvec = train.fm,
  vec_method = c("gene","topic"), #"g","t","genevec","topicvec"
  feature.name = "gsvec",
  export_predict.gs.vector = T,
  sparse_thres.p = 0, #if sparse_thres.p=0, non-sparse 
  calc.Fisher = T, #F
  cosine.thres = 0.5
  ){
  
  # make gsvec
  cat(paste0("\n - predict gene signature vectors"))
  if(tolower(substr(vec_method,1,1) == "g")){
    out.mat <- predict_GsVec_from.GeneVec(	
      train.gene.feature.mat = train.gsvec,
      train.data = train.data,
      val.data = val.data
    )
  }
  if(tolower(substr(vec_method,1,1) == "t")){
    out.mat <- predict_GsVec_from.TopicVec(
      train.topicvec = train.gsvec,
      train.data = train.data,
      val.data = val.data,
      sparse_thres.p = sparse_thres.p
      )    
  }
  if(export_predict.gs.vector == T){
    cat(paste0("\n - export signature vector to pred.val_",feature.name,".txt"))
    out <- data.table(id = rownames(out.mat), out.mat, stringsAsFactors = F)
    fwrite(out,paste0("pred.val_",feature.name,".txt"),sep="\t",quote=F,row.names=F)

  }

  train.mat <- out.mat[1:nrow(train.data),]
  val.mat <- out.mat[(nrow(train.data)+1):nrow(out.mat),]
  
  ## calc similarity between signature vector
  cat(paste0("\n - calc similarity between signature vector"))
  out_sim.vec <- similarity_vectors(
    train_feature.matrix = train.mat,
    val_feature.matrix = val.mat,
    feature.name = feature.name
  )
  
  colnames(train.data) <- paste0("t",colnames(train.data))
  colnames(val.data) <- paste0("v",colnames(val.data))

  ## calc Fisher and marged results
  if(calc.Fisher == T){
    cat(paste0("\n - calc Fisher exact test - "))
    out_fisher <- gs.enrich_fisher(
      train.data = train.data,
      val.data = val.data,
      use_gene.population.txt = NA,	
      export_fisher.result = F
    )
    out <- merge(val.data[,c(1,2,4)],merge(train.data[,c(1,2,4)],merge(out_fisher,out_sim.vec[,c(1,4)],by="vtid"),by="tid"),by="vid")
    out <- out[,c(7,1,4,2,3,5,6,9,8,11,10)]
  }else{
    out <- merge(val.data[,c(1,2,4)],merge(train.data[,c(1,2,4)],out_sim.vec,by="tid"),by="vid")
    out <- out[,c(7,1,4,2,3,5,6,8)]
  }
    colnames(out)[4] <- "val.data_name"
    colnames(out)[5] <- "val.data_number"
    colnames(out)[6] <- "train.data_name"
    colnames(out)[7] <- "train.data_number"

    cat(paste0("\n - Write results to ",paste0("gsvec_",feature.name,".txt")))
    fwrite(out,paste0("gsvec_",feature.name,".txt"),sep="\t",quote=F,row.names=F)

    cat(paste0("\n - completed."))

  return(out)
}

#------------------------------------------------------------------------------------------------------ PCA.tSNE

pca.tsne_GsVec <- function(
	gsvec.mat = pred.val_tfidf,		# col = feature, row=sigs
	gsvec.train.data = train.data,
	gsvec.val.data = val.data,
	gs.group_id.group.mat = group.mat, 	#id,Group matrix
	colors = c(brewer.pal(9, "Set1"),brewer.pal(8, "Set2"),brewer.pal(8, "Set3"),brewer.pal(8, "Dark2")),
	pca.thres = 0.95,#NA
	out_name = "pca.tsne",
	centering = F,
	tsne_max.iter = 500,
	tsne_sta = 200
){

# 	colors <- c(brewer.pal(9, "Set1"),brewer.pal(8, "Set2"),brewer.pal(8, "Set3"),brewer.pal(8, "Dark2")),

	gsvec.val.data <- rbind(train.data,gsvec.val.data)
	group.mat <- data.frame(ID = gsvec.val.data$name, Group = gsvec.val.data$id)
	group.mat$Group <- gsub("t.*","training",group.mat$Group)
	group.mat$Group <- gsub("v.*","validation",group.mat$Group)

	if(is.na(pca.thres)){
		pcs <- t(gsvec.mat)
	}else{
		pcadata <- prcomp(gsvec.mat,center = centering)
		indx <- !(summary(pcadata)$importance[3,] > pca.thres)
		if(sum(indx)==0){indx[1:2] <- T}
		pcs <- pcadata$x[,indx]
	}

	tsnedata <- Rtsne(pcs, dims = 2, initial_dims = ncol(pcs), perplexity = 30,
		theta = 0.5, check_duplicates = F, pca = F, max_iter = tsne_max.iter, verbose = T, is_distance = FALSE, Y_init = NULL, eta = tsne_sta, exaggeration_factor = 12)

	colnames(gs.group_id.group.mat) <- c("ID","Group")
	data <- data.frame(id=rownames(pcs), TSNE_1 = tsnedata$Y[,1], TSNE_2 = tsnedata$Y[,2], stringsAsFactors=FALSE)
	data <- merge(merge(data,gsvec.val.data,by="id"),gs.group_id.group.mat,by.x="name",by.y="ID")
	 
	g <- ggplot(data, aes(x = TSNE_1, y = TSNE_2, color = Group, size=number))
	g <- g + geom_point(alpha=0.3)
	g <- g + scale_color_manual(values = colors)
	g <- g + ggtitle(paste0("tSNE__",out_name))
	ggsave(file = paste0("tSNE__",out_name,".png"), device="png", units="in", dpi = 150, width = 7, height = 5)

	return(data)

}

## -------------------------------------------------------------------------------------------------- EOD
