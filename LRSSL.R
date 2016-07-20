###########################read training data#######################
chem.mat <-read.table("./drug_pubchem_mat.txt", sep = "\t", header = T, row.names = 1)
dom.mat <- read.table("./drug_target_domain_mat.txt", sep = "\t", header = T, row.names = 1)
go.mat <- read.table("./drug_target_go_mat.txt", sep = "\t", header = T, row.names = 1)
phar.mat <- read.table("./drug_dis_mat.txt", sep = "\t", header = T, row.names = 1)
dis.sim <- read.table("./disease_similarity.txt", sep = "\t", header = T, row.names = 1)

X1 <- as.matrix(chem.mat)
X2 <- as.matrix(dom.mat)
X3 <- as.matrix(go.mat)
Y <- as.matrix(phar.mat)

n <- nrow(Y)
c <- ncol(Y)

Xs <- list(X1, X2, X3)
###########################caculate knn graph########################
k <- 10

S1 <- X1%*%t(X1)/sqrt(rowSums(X1)%*%t(rowSums(X1)))
S2 <- X2%*%t(X2)/sqrt(rowSums(X1)%*%t(rowSums(X1)))
S3 <- X3%*%t(X3)/sqrt(rowSums(X3)%*%t(rowSums(X3)))

