##########   Clustering Methods
##########   Nathaniel E. Helwig (helwig@umn.edu)
##########   Updated: 27-Mar-2017


### load 'maps' package
if(!require(maps)){
  install.packages("maps")
  library(maps)
}

### load 'RColorBrewer' package
if(!require(RColorBrewer)){
  install.packages("RColorBrewer")
  library(RColorBrewer)
}


#####   HIERARCHICAL CLUSTERING   #####

# look at states data
?state.x77
vars <- c("Income","Illiteracy","Life Exp","HS Grad")
head(state.x77[,vars])
apply(state.x77[,vars], 2, mean)
apply(state.x77[,vars], 2, sd)

# create distance (raw and standardized)
distraw <- dist(state.x77[,vars])
diststd <- dist(scale(state.x77[,vars]))

# look at distances
dim(as.matrix(distraw))
dim(as.matrix(diststd))
as.matrix(distraw)[1:4,1:4]
as.matrix(diststd)[1:4,1:4]

# hierarchical clustering (raw data)
hcrawSL <- hclust(distraw, method="single")
hcrawCL <- hclust(distraw, method="complete")
hcrawAL <- hclust(distraw, method="average")

# plot results (raw data)
dev.new(width=12, height=4, noRStudioGD=TRUE)
par(mfrow=c(1,3))
plot(hcrawSL)
plot(hcrawCL)
plot(hcrawAL)

# hierarchical clustering (standardized data)
hcstdSL <- hclust(diststd, method="single")
hcstdCL <- hclust(diststd, method="complete")
hcstdAL <- hclust(diststd, method="average")

# plot results (standardized data)
dev.new(width=12, height=4, noRStudioGD=TRUE)
par(mfrow=c(1,3))
plot(hcstdSL)
plot(hcstdCL)
plot(hcstdAL)

# plot results (standardized data close-up)
dev.new(width=10, height=6, noRStudioGD=TRUE)
plot(hcstdCL)


#####   NON-HIERARCHICAL (K MEANS) CLUSTERING   #####

# look at states data
?state.x77
vars <- c("Income","Illiteracy","Life Exp","HS Grad")
head(state.x77[,vars])
apply(state.x77[,vars], 2, mean)
apply(state.x77[,vars], 2, sd)

# fit k means for k = 2, ..., 10 (raw data)
kmlist <- vector("list", 9)
for(k in 2:10){
  set.seed(1)
  kmlist[[k-1]] <- kmeans(state.x77[,vars], k, nstart=5000)
}

# scree plot (raw data)
tot.withinss <- sapply(kmlist, function(x) x$tot.withinss)
dev.new(width=6,height=4,noRStudioGD=TRUE)
plot(2:10, tot.withinss / kmlist[[1]]$totss, type="b", xlab="# Clusters", 
     ylab="SSW / SST", main="Scree Plot: Raw Data")

# plot results (raw data)
dev.new(width=9, height=6, noRStudioGD=TRUE)
par(mfrow=c(2,2))
for(k in 3:6){
  map(database = "state")
  title(paste0("K=",k," Clusters: Raw Data"))
  cols <- brewer.pal(k, "Paired")
  for(j in 1:k){
    ix <- names(which(kmlist[[k-1]]$cluster==j))
    if(length(ix) > 1) map(database = "state", regions = ix, col = cols[j], fill=T, add=TRUE)
  }
}

# fit k means for k = 2, ..., 10 (standardized data)
Xs <- scale(state.x77[,vars])
kmlist.std <- vector("list", 9)
for(k in 2:10){
  set.seed(1)
  kmlist.std[[k-1]] <- kmeans(Xs, k, nstart=5000)
}

# scree plot (standardized data)
tot.withinss.std <- sapply(kmlist.std, function(x) x$tot.withinss)
dev.new(width=6,height=4,noRStudioGD=TRUE)
plot(2:10, tot.withinss.std / kmlist.std[[1]]$totss, type="b", 
     xlab="# Clusters", ylab="SSW / SST", main="Scree Plot: Std. Data")

# plot results (standardized data)
dev.new(width=9, height=6, noRStudioGD=TRUE)
par(mfrow=c(2,2))
for(k in 3:6){
  map(database = "state")
  title(paste0("K=",k," Clusters: Std. Data"))
  cols <- brewer.pal(k, "Paired")
  for(j in 1:k){
    ix <- names(which(kmlist.std[[k-1]]$cluster==j))
    if(length(ix) > 1) map(database = "state", regions = ix, col = cols[j], fill=T, add=TRUE)
  }
}

