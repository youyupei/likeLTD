# k means at a single locus
kmeanslocus = function(means,data)
	{
	# no result if no data
	if(length(data)==0) return(NULL)
	diffs = sapply(means,FUN=function(x) abs(x-as.numeric(data)))
	# best group if one data point
	if(length(data)==1) return(which.min(diffs))
	cols = col(diffs)
	rows = row(diffs)
	sorted = order(unlist(diffs))
	# sort everything
	sortedCols = cols[sorted]
	sortedRows = rows[sorted]
	sortedVals = c(diffs)[sorted]
	# assign best fitting groups
	out = vector(length=nrow(diffs))
	for(i in 1:5)
		{
		# best fit
		colIndex = sortedCols[1]
		rowIndex = sortedRows[1]
		out[rowIndex] = colIndex
		# remove all possible assignment for that peak
		toRemove = which(sortedRows==rowIndex)
		# remove all possible assignments for that individual (max=2 per person)
		if(length(which(out==colIndex))==2) toRemove = unique(c(toRemove,which(sortedCols==colIndex)))
		# resort
		sortedCols = sortedCols[-toRemove]
		sortedRows = sortedRows[-toRemove]
		sortedVals = sortedVals[-toRemove]
		}
	return(out)
	}

# single iteration of k means
kmeansIt = function(means,heights)
	{
	groups = sapply(heights,FUN=function(x) sapply(names(heights[[1]]),FUN=function(y) kmeanslocus(means,x[[y]])),simplify=FALSE)
	newMeans = sapply(1:length(means),FUN=function(x) mean(as.numeric(unlist(heights)[which(unlist(groups)==x)])))
	return(newMeans)
	}

# single start kMeans
startKmeans = function(heights,nU,nIt=10)
	{
	means = runif(nU,0,5000)
	for(i in 1:nIt)
		{
		means = kmeansIt(means,heights)
		}
	return(means)
	}

# multiple kmeans
multiKmeans = function(heights,nU,nStart=10,nIt=10)
	{
	multi = replicate(nIt,startKmeans(heights,nU,nIt))
	multi = apply(multi,MARGIN=2,FUN=function(x) sort(x,na.last=TRUE))
	return(rowMeans(multi,na.rm=TRUE))
	}


