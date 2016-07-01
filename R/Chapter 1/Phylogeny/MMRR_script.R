# MMRR

# The read.matrix function requires {tseries} package to be installed and loaded.
#install.packages("tseries")
install.packages("tseries")
install.packages("car")
library(tseries)
library(car)

# Read the matrices from files.
# If the files have a row as a header (e.g. column names), then specificy 'header=TRUE', default is 'header=FALSE'.
alk = as.matrix(read.csv(row.names = "pixel_ID", file = "alk.csv", header = TRUE))
antalk = as.matrix(read.csv(row.names = "pixel_ID", file = "antalk.csv", header = TRUE))
antgut = as.matrix(read.csv(row.names = "pixel_ID", file = "antgut.csv", header = TRUE))
geo = as.matrix(read.csv(row.names = "pixel_ID", file = "geo.csv", header=TRUE))

diag(alk) = NA
diag(antalk) = NA
diag(antgut) = NA
diag(geo) = NA

# Plot pairwise distances
plot = matplot(antalk, alk, type="p",lty=1,lwd=1,lend=par("lend"),pch=20,col=1:1,cex=NULL,bg=NA, xlab="Alkaloid-bearing ant dissimilarity",ylab="Alkaloid dissimilarity")
# to plot line:
alk = as.dist(alk)
antalk = as.dist(antalk)
regLine(lm(alk~antalk, data=plot), lwd=3, col=palette()[1])

# Same for other predictors
plot = matplot(antgut, alk, type="p",lty=1,lwd=1,lend=par("lend"),pch=20,col=1:1,cex=NULL,bg=NA, xlab="Eaten ant dissimilarity",ylab="Alkaloid dissimilarity")
alk = as.dist(alk)
antgut = as.dist(antgut)
regLine(lm(alk~antgut, data=plot), lwd=3, col=palette()[1])

plot = matplot(geo, alk, type="p",lty=1,lwd=1,lend=par("lend"),pch=20,col=1:1,cex=NULL,bg=NA, xlab="Geographic distance",ylab="Alkaloid dissimilarity")
alk = as.dist(alk)
geo = as.dist(geo)
regLine(lm(alk~geo, data=plot), lwd=3, col=palette()[1])

plot = matplot(geo, antalk, type="p",lty=1,lwd=1,lend=par("lend"),pch=20,col=1:1,cex=NULL,bg=NA, xlab="Geographic distance",ylab="Alkaloid-bearing ant dissimilarity")
antalk = as.dist(antalk)
geo = as.dist(geo)
regLine(lm(antalk~geo, data=plot), lwd=3, col=palette()[1])

plot = matplot(geo, antgut, type="p",lty=1,lwd=1,lend=par("lend"),pch=20,col=1:1,cex=NULL,bg=NA, xlab="Geographic distance",ylab="Eaten ant dissimilarity")
antgut = as.dist(antgut)
geo = as.dist(geo)
regLine(lm(antgut~geo, data=plot), lwd=3, col=palette()[1])


# Make a list of the explanatory (X) matrices.
# Names are optional.  Order doesn't matter.
# Can include more than two matrices, if desired.
Xmats = list(antalk = antalk)
Xmats = list(antgut = antgut)
Xmats = list(geo = geo)

# Run MMRR function using genMat as the response variable and Xmats as the explanatory variables.
# nperm does not need to be specified, default is nperm=999)
MMRR(alk,Xmats,nperm=10000)
MMRR(antalk,Xmats,nperm=10000)
MMRR(antgut,Xmats,nperm=10000)

