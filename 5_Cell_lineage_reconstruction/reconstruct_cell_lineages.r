library(qgraph)
library(igraph)

args = commandArgs(trailingOnly=TRUE) ### input

data_dir <- paste("$HOME/target_PCR_data/mut/",args[1],sep="")
theta <- 5
epsilon2 <- 0.01

#read data
freq <- as.matrix(read.csv(paste(data_dir,"/lineage/data_VAF.txt",sep=""), row.names=1, sep="\t"))
freq[is.na(freq)] <- 0

#add root data (freq=0.5)
freq <- rbind(freq, 0.5)
rownames(freq)[dim(freq)[1]] <- "root"

#remove low allele frequency data
for(i in dim(freq)[1]:1){
	if(max(freq[i,]) < 0.05){
		freq <- freq[-i,]
	}
}

mnames <- rownames(freq)
mutation_num <- dim(freq)[1]
sample_num <- dim(freq)[2]

#distance between mutations
Distance_matrix <- matrix(0, nrow=mutation_num, ncol=mutation_num)
for(i in 1:mutation_num){
	for(j in 1:mutation_num){
		Distance_matrix[i,j] <- sum((freq[i,]-freq[j,])**2) / (sample_num * 0.5 * (mean(freq[i,])+mean(freq[j,])))
	}
}
diag(Distance_matrix) <- 0

minus_log_Distance_matrix <- matrix(nrow=mutation_num, ncol=mutation_num)
for(i in 1:mutation_num){
	for(j in 1:mutation_num){
		if(Distance_matrix[i,j] == 0){
			minus_log_Distance_matrix[i,j] <- Inf
		}
		else{
			minus_log_Distance_matrix[i,j] <- -log(Distance_matrix[i,j])
		}

		#remove edge that satisfy -log(d(i,j)) < theta
		if(minus_log_Distance_matrix[i,j] < theta){
			minus_log_Distance_matrix[i,j] <- 0
		}
	}
}

#clustering mutation into mutation nodes based on connected components of the above graph
g <- graph.adjacency(minus_log_Distance_matrix, weighted=TRUE, mode="undirected")
cls <- clusters(g, "weak")
mutation_node_num <- cls$no
mutation_node_list <- list()
for(i in 1:mutation_node_num){
	mutation_node_list[[i]] <- which(cls$membership==i)
}

mutation_node_mnames <- rep("", mutation_node_num)
for(i in 1:mutation_node_num){
	mutation_node_mnames[i] <- rownames(freq)[mutation_node_list[[i]]][[1]]
	if(length(mutation_node_list[[i]]) >= 2){
		for(j in 2:length(mutation_node_list[[i]])){
			mutation_node_mnames[i] <- paste(mutation_node_mnames[i], rownames(freq)[mutation_node_list[[i]]][[j]], sep="\n")
		}
	}
}

#mean allele frequency of mutation node
mean_freq <- matrix(0, nrow=mutation_node_num, ncol=dim(freq)[2])
for(i in 1:mutation_node_num){
	if(length(mutation_node_list[[i]]) > 1){
		mean_freq[i,] <- colSums(freq[mutation_node_list[[i]],]) / length(mutation_node_list[[i]])
	}
	else{
		mean_freq[i,] <- freq[mutation_node_list[[i]][1],]
	}
}

#inference of one mother two daughter nodes relationships
epsilon1 <- 0.005
Issum_matrix <- matrix(rep(0,mutation_node_num**2), nrow=mutation_node_num, ncol=mutation_node_num)
while(epsilon1 > 0.0){
	print(paste("current epsilon1 =",epsilon1,sep=" "))
	Issum_matrix <- diag(rep(0,mutation_node_num))

	for(i in 1:mutation_node_num){
		for(j in 1:(mutation_node_num-1)){
			for(k in (j+1):mutation_node_num){
				if(i == j || i == k){
					next
				}

				Y <- mean_freq[i,]
				X <- mean_freq[j,] + mean_freq[k,]

				tmp <- sum((Y-X)**2/sample_num)
				tmp <- tmp / mean(X)

				if(tmp < epsilon1){
					Issum_matrix[j,i] <- 1
					Issum_matrix[k,i] <- 1
				}
			}
		}
	}

	if(sum(apply((Issum_matrix!=0), MARGIN=2, sum) > 2) > 0 || (sum(apply((Issum_matrix!=0), MARGIN=1, sum) > 1)) > 0){
		epsilon1 <- epsilon1 - 0.00025
	}
	else{
		break
	}
}

jpeg(paste(data_dir,"/lineage/lineage_step1.jpeg",sep=""),height=960, width=960, res=144)
q1 <- qgraph(t(Issum_matrix!=0), layout="spring", labels=mutation_node_mnames, edge.color="red")
dev.off()


#additional lineage inference
Total_lineage_tree <- (Issum_matrix!=0)

leaf_nodes <- c()
orphaned_nodes <- c()

old_Total_lineage_tree <- Total_lineage_tree
mutation_node_mnames2 <- mutation_node_mnames
mutation_node_num2 <- mutation_node_num
mean_freq2 <- mean_freq

#find leaf nodes
root_id <- which(mutation_node_mnames=="root")
check_leaf_node <- function(i){
	if(sum(Total_lineage_tree[,i]) > 0){
		check_leaf_node(which(Total_lineage_tree[,i]==1)[1])
		check_leaf_node(which(Total_lineage_tree[,i]==1)[2])
	}
	else{
		leaf_nodes <<- c(leaf_nodes, i)
	}
}
check_leaf_node(root_id)

#find orphaned nodes
for(i in 1:mutation_node_num){
	if(i == root_id){
		next
	}

	if(sum(Total_lineage_tree[i,]) == 0){
		orphaned_nodes <<- c(orphaned_nodes, i)
	}
}

#investigate orphaned nodes in allele frequency ascending order
tmp <- order(rowSums(mean_freq2[orphaned_nodes,]), decreasing=TRUE)
orphaned_nodes <- orphaned_nodes[tmp]

tmp_count <- 1
for(i in 1:length(orphaned_nodes)){
	child_id <- orphaned_nodes[i]

	tmp_parent <- rep(0, mutation_node_num2)
	for(j in 1:length(leaf_nodes)){
		parent_id <- leaf_nodes[j]
		
		if(sum((mean_freq2[parent_id,]-mean_freq2[child_id,]+epsilon2)<0) == 0){
			tmp_parent[parent_id] <- 1
		}
	}

	#only if the orphaned node has one possible mother node
	if(sum(tmp_parent) == 1){
		parent_id <- which(tmp_parent==1)
		Total_lineage_tree[child_id,parent_id] <- 1

		#add pseudo node (push pseudo node for mean_freq and tree)
		tmp_tree <- matrix(0, nrow=dim(Total_lineage_tree)[1]+1, ncol=dim(Total_lineage_tree)[2]+1)
		tmp_tree[1:(dim(Total_lineage_tree)[1]),1:(dim(Total_lineage_tree)[2])] <- Total_lineage_tree
		tmp_tree[dim(tmp_tree)[1],parent_id] <- 1

		mutation_node_mnames2 <- c(mutation_node_mnames2, paste("Pseudo_",tmp_count,sep=""))
		tmp_count <- tmp_count + 1

		tmp_mean_freq <- matrix(0, nrow=dim(mean_freq2)[1]+1, ncol=dim(mean_freq2)[2])
		tmp_mean_freq[1:dim(mean_freq2)[1],] <- mean_freq2
		tmp_mean_freq[dim(mean_freq2)[1]+1,] <- mean_freq2[parent_id,] - mean_freq2[child_id,]
		
		Total_lineage_tree <- tmp_tree
		mean_freq2 <- tmp_mean_freq
		mutation_node_num2 <- mutation_node_num2+1

		#remake assigned leaf node
		leaf_nodes <<- c()
		check_leaf_node(root_id)
	}
	else if(sum(tmp_parent) > 1){
        print(paste("Node :", mutation_node_mnames2[child_id], "cannot be assigned uniquely."))
		for(j in 1:sum(tmp_parent)){
			print(paste("  ", mutation_node_mnames2[child_id], " < ", mutation_node_mnames2[which(tmp_parent==1)[j]]))
		}
	}
	else{
		print(paste("Node :", mutation_node_mnames2[child_id], "cannot be assigned any leaf node."))
	}
}


setwd(paste(data_dir, "/lineage", sep=""))
step1 <- gsub("\n", ", ", mutation_node_mnames2)
step1 <- gsub("Mos", "", step1)
step2 <- gsub("(([^,]+, ){3}[^,]+), ", "\\1\n", step1)
qgraph(t(Total_lineage_tree), layout = "spring", labels = step2, edge.color = "red", filetype = "pdf", label.cex=1.2,vsize=6,vsize2=6)

jpeg(paste(data_dir,"/lineage/lineage_step2.jpeg",sep=""),height=960, width=960, res=600)
q2 <- qgraph(t(Total_lineage_tree), layout = "spring", labels = step2, edge.color = "red", label.cex=1.2,vsize=8,vsize2=8)
dev.off()

node <- mutation_node_mnames2
node <- gsub("\n", " ", node)
write.table(node, paste(data_dir,"/lineage/node.txt",sep=""))



library(ggtree)
library(ggforce)

# get relationship of nodes
adj_to_newick <- function(adj, labels) {
	build_newick <- function(node) {
		children <- which(adj[node, ] == 1)
		if (length(children) == 0) {
			return(labels[node])
		} else {
			subtrees <- sapply(children, build_newick)
			return(paste0("(", paste(subtrees, collapse=","), ")", labels[node]))
		}
	}
	root <- which(mutation_node_mnames2 == "root")
	if (length(root) != 1) stop("Error: Cannot identify unique root node")
	newick <- paste0(build_newick(root), ";")
	return(newick)
}

newick_string <- adj_to_newick(t(Total_lineage_tree), mutation_node_mnames2)
library(ape)
tree <- read.tree(text = newick_string)
tree$edge.length <- rep(1, nrow(tree$edge))


library(dplyr)

node_df <- data.frame(
	node = mutation_node_mnames2,
	vaf = mean_freq2[ ,1],
	vaf2 = 2 * mean_freq2[ ,1],
	mut_count = c(
	sapply(mutation_node_list, length),
	rep(0, length(mutation_node_mnames2) - length(mutation_node_list)))
)

node_df$vaf2[node_df$vaf2 > 1] <- 1


# plot by ggtree
p <- ggtree(tree, layout = "dendrogram") 

p_data <- p$data
df_plot <- left_join(node_df, p_data, by = c("node" = "label"))

df_plot <- df_plot %>%
	mutate(
		vaf2 = ifelse(grepl("^Pseudo_|^root$", node), 0, vaf2),
		fill_color = ifelse(vaf2 == 0, "white", "cornflowerblue")
)

df_pie <- df_plot %>%
	rowwise() %>%
	mutate(
		start1 = pi/2,                       
		end1 = pi/2 - vaf2 * 2 * pi,       
		start2 = end1,
		end2 = pi/2 - 2 * pi ,     
		is_white = vaf2 == 0
)

df_label <- df_pie[df_pie$vaf2 > 0,]

df_pie$circle_color <- ifelse(df_pie$mut_count > 4, "#dc3023", "black")
df_pie$circle_size <- ifelse(df_pie$mut_count > 4, 1.5, 0.6)

pic<- p +
	geom_arc_bar(
		data = df_pie,
		aes(x0 = x, y0 = y, r0 = 0, r = 0.4, start = start1, end = end1, fill = "#C4D0E4"),
		color = "#C4D0E4", size = 0.4, inherit.aes = FALSE
	) +
	geom_arc_bar(
		data = df_pie,
		aes(x0 = x, y0 = y, r0 = 0, r = 0.4, start = start2, end = end2, fill = "white"),
		color = "white", size = 0, inherit.aes = FALSE
	) +
	geom_circle(
		data = df_pie,
		aes(x0 = x, y0 = y, r = 0.4, color = circle_color, size=circle_size), fill = NA
	)+
	scale_color_identity()+
	scale_size_identity()+
	geom_text(
		data = df_label,
		aes(x = x - 0.12, y = y, label = mut_count),
		size = 8.5,
		hjust = 0.5,
		color = "gray40"
	) +
	geom_text(
		data = df_label,
		aes(x = x + 0.12, y = y, label = ifelse(mut_count > 1, "muts", "mut")),
		size = 8.5,
		hjust = 0.5,
		color = "gray30"
	) +
	# geom_text(
		#   data = df_label,
		#   aes(x = x + 0.5, y = y, label = sprintf("%.1f%%", vaf2 * 100)),
		#   size = 8
	# ) +
	scale_fill_identity() +
	theme(legend.position = "none")

gb <- ggplot_build(p)
xrange <- max(gb$data[[1]]$y)
yrange <- max(gb$data[[1]]$x/1.5+1)
width <- xrange * 1.6
height <- yrange * 2.1

ggsave(paste(data_dir, "/lineage/tree_plot.pdf", sep=""),pic, width = width, height = height, dpi = 600)
