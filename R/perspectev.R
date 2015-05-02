mcpRange <- 
function(t1,t2,vlist){ #calculate log minimum convex polygon range
	t1 = as.numeric(t1)
	t2 = as.numeric(t2)
	l = unique(data.frame(cbind(t1,t2))) 
	if(nrow(l) < 3){return(NA)} #mcp needs at least three unique values
	range = attributes(hr.mcp(cbind(l[,1],l[,2]),plot=FALSE,n.min=3))$area
	if(range == 0){return(NA)} #close points can be rounded down to zero
	return(log(range))
}

t1Range <- 
function(t1,t2,vlist){ #calculate lattitude breadth in decimal
	if(length(unique(t1)) < 2){return(NA)} #need at least two unique values
	t1 = as.numeric(t1)
	return(max(t1) - min(t1)) #T1 = latitude, T2 = longitude
}

t2Range <- 
function(t1,t2,vlist){ #calculate longitude breadth
	if(length(unique(t1)) < 2){return(NA)} #need at least two unique values
	t2 = as.numeric(t2)
	return(max(t2) - min(t2)) #T1 = latitude, T2 = longitude
}

t1Var <- 
function(t1,t2,vlist){ #calculate variance of trait 1
	if(length(unique(t1)) < 2){return(NA)} #need at least two unique values
	t1 = as.numeric(t1)
	return(var(t1))
}

t2Var <- 
function(t1,t2,vlist){ #calculate variance of trait 2
	if(length(unique(t2)) < 2){return(NA)} #need at least two unique values
	t2 = as.numeric(t2)
	return(var(t2))
}

t1t2Covar <- 
function(t1,t2,vlist){ #calculate covariance between trait 1 and trait 2
	if(nrow(unique(cbind(t1,t2))) < 2){return(NA)} #need at least two unique values
	t1 = as.numeric(t1)
	t2 = as.numeric(t2)
	return(cov(t1,t2))
}

t1Unique <- 
function(t1,t2,vlist){ #calculate number of unique trait 1 values
	if(length(t1) < 1){return(NA)} #need at least 1 unique value
	return(unique(t1))
}

t2Unique <- 
function(t1,t2,vlist){ #calculate number of unique trait 2 values
	if(length(t2) < 1){return(NA)} #need at least 1 unique value
	return(unique(t2))
}

t1t2Unique <- 
function(t1,t2,vlist){ #calculate unique combinations of trait 1 and trait 2
	if(nrow(unique(cbind(t1,t2))) < 2){return(NA)} #need at least two unique values
	l = unique(data.frame(cbind(t1,t2))) 
	return(nrow(l))
}

dnaTheta <-
function(t1,t2,vlist){ #calculate mean pairwise difference for set of sequences
	if(length(t1) < 2){return(NA)} #need at least two sequences
	if(is.na(vlist)){return(NA)}#must provide global distance matrix
	if(class(vlist) != "list"){return(NA)} #vlist must be a list
	compMat = vlist[[1]][t1,t1]
	return(mean(compMat[upper.tri(compMat)]))
}

hr.mcp <-
function(x,y=NULL,n.min=50,plot=TRUE,add=FALSE,ID=NULL,...){
  if(is.null(y))xy <- x
  if(!is.null(y))xy <- cbind(x,y)
  if(inherits(x,what="data.frame"))xy <- as(xy,"matrix")
  #Missing values allowed but dropped and trigger a warning 
  xy. <- xy
  xy <- na.omit(xy)
  if(nrow(xy.)!=nrow(xy))warning("Missing values dropped")
  #IF XY DOES NOT SUPPORT ESTIMATION, STOP HERE
  if(nrow(xy)<n.min){
    warning("Fewer than 'n.min' qualifying observations")
    return(NA)
  }
  mcp <- xy[chull(xy),]
  mcp <- rbind(mcp,mcp[1,])
  mcp <- Polygon(mcp,hole=FALSE)
  if(is.null(ID))ID <- as.logical(NA)
  mcp <- Polygons(list(mcp),ID=ID)
  mcp
}

perspectev.read <-
function(data,extinctionAge,occurrenceAge,upper,lower,t1,t2,traitfun=mcpRange,vlist=NULL,trim=TRUE,projection=FALSE){
	p = data
	
	cn = colnames(p)
	if(!occurrenceAge %in% cn){stop(paste("Occurrence age column",occurrenceAge,"not in dataframe!"))}
	if(!upper %in% cn){stop(paste("Upper level column",upper,"not in dataframe!"))}
	if(!lower %in% cn){stop(paste("Lower level column",lower,"not in dataframe!"))}
	if(!t1 %in% cn){stop(paste("t1 column",t1,"not in dataframe!"))}
	if(!t2 %in% cn){stop(paste("t2 column",t2,"not in dataframe!"))}

	if(projection == TRUE){ #convert coordinates to Gall equal area projection
		proj = mapproject(p[,t2],p[,t1],projection="mollweide")
		p[,t2] = proj$x * 100
		p[,t1] = proj$y * 100
	} 
	before = p[p[,occurrenceAge] >= extinctionAge,]
	after = p[p[,occurrenceAge] < extinctionAge,]	
	beforesp = unique(paste(before[,upper],before[,lower]))
	aftersp = paste(after[,upper],after[,lower])
	survsp = unique(aftersp[aftersp %in% beforesp]) #species found across extinction age
	final = data.frame(matrix(NA,nrow=length(beforesp),ncol=5))
	colnames(final) = c("Upper","Lower","Survivorship","T1","T2")
	rownames(final) = beforesp
	for(i in 1:dim(before)[1]){ #fill in data matrix - t1 = latitude, t2 = longitude
		spec1 = paste(before[i,upper],before[i,lower])
		
		if(!is.na(final[spec1,]$T1)){
			final[spec1,]$T1 = paste(final[spec1,]$T1,before[i,t1],sep=";")
		}else{final[spec1,]$T1 = as.character(before[i,t1])}

		if(!is.na(final[spec1,]$T2)){
			final[spec1,]$T2 = paste(final[spec1,]$T2,before[i,t2],sep=";")
		}else{final[spec1,]$T2 = as.character(before[i,t2])}
		if(is.na(final[spec1,]$Lower)){final[spec1,]$Lower = spec1}
		if(is.na(final[spec1,]$Upper)){final[spec1,]$Upper = as.character(before[i,upper])}
		if(is.na(final[spec1,]$Survivorship)){final[spec1,]$Survivorship = sum(spec1 %in% survsp)}
	}
	if(trim == TRUE){ #remove lower units that give NA for trait calculation
		nfinal = data.frame()
		for(i in 1:nrow(final)){
			occurrences = cbind((unlist(strsplit(final[i,]$T1,split=";"))),(unlist(strsplit(final[i,]$T2,split=";"))))
			if(!is.na(traitfun(occurrences[,1],occurrences[,2],vlist))){
				nfinal = rbind(nfinal,final[i,])
			}
		}
		final = nfinal
		colnames(final) = c("Upper","Lower","Survivorship","T1","T2")
	}
	return(final)
}


perspectev.calc <-
function(data,traitfun,vlist=list()){ #Calculate correlation between range and survivorship for a given set of genera
	low = data.frame()
	for(i in 1:nrow(data)){
		t1 = unlist(strsplit(data[i,]$T1,split=";"))
		t2 = unlist(strsplit(data[i,]$T2,split=";"))
		value = traitfun(t1,t2,vlist)
		low = rbind(low,c(data[i,]$Survivorship,value))
	}
	up = data.frame()
	uppers = unique(data$Upper)
	for(g in uppers){
		temp = data[data$Upper == g,]
		t1 = unlist(strsplit(temp$T1,split=";"))
		t2 = unlist(strsplit(temp$T2,split=";"))
		value = traitfun(t1,t2,vlist)
		survivorship = 1-prod(1-temp$Survivorship)#probability that any species survives
		up = rbind(up,c(survivorship,value))
	}
	q = quantile(up[,2])[2:4] 		 #get quantiles 
	up[,2] = scale(up[,2]) 			 #scale value
	low[,2] = scale(low[,2])
	if(var(up[,1]) == 0 || var(up[,2]) == 0){
		correlation = NA
	}else{
		correlation = cor(up[,1],up[,2]) #calculate correlation
	}
	rownames(low) = rownames(data)
	rownames(up) = uppers
	colnames(up) = c("Survivorship","Trait")
	colnames(low) = c("Survivorship","Trait")
	total = list()
	total$lower = low
	total$upper = up
	total$stats = c(correlation,q)
	return(total)
}


perspectev.permute <-
function(data,permutations,traitfun=mcpRange,vlist){
	simulations = c(); count = 0; fcount = 0
	while(count < permutations){
		temp = data 		
		allg = temp$Upper
		temp$Upper = allg[sample(1:length(allg),replace=FALSE)] 
		p = perspectev.calc(temp,traitfun,vlist)
		if(!is.na(p$stats[1])){
			simulations = rbind(simulations,p$stats)
			count = count + 1
		}
		fcount = fcount + 1
		if(fcount >  100 & (fcount-count)/fcount > 0.9){
			stop('Error: More than 90% of permutations yeild either no victims or no survivors')
		}
	}
	rownames(simulations)=NULL
	return(simulations)
}


perspectev.test <-
function(data,iterations=1000,cores=1,traitfun=mcpRange,vlist=NULL){
	sim = c(); 	i = NULL
	real = perspectev.calc(data,traitfun,vlist)
	if(cores > 1){ 
		div = floor(iterations/cores)
		plist = rep(div,length=cores)
		plist[length(plist)] = plist[length(plist)] + iterations %% cores
		cl = makeCluster(cores)
		registerDoParallel(cl)
		sim = foreach(i = 1:cores,.packages=c('perspectev'),.combine='rbind') %dopar% perspectev.permute(data,plist[i],traitfun,vlist)
		stopCluster(cl)
	}else{
		sim = perspectev.permute(data,iterations,traitfun,vlist)
	}
	total = list() #calculate genus and species p value
	total$correlation_permuted = sim[,1]
	total$correlation_observed= real$stats[1]
	total$pvalue = sum(sim[,1] >= real$stats[1])/nrow(sim)
	permutedqs = cbind(sim[,2],sim[,3],sim[,4])
	colnames(permutedqs) = c("25%","50%","75%")
	total$permuted_quantiles = permutedqs
	return(total)
}


perspectev.simulation.permute <-
function(data,simulations,intercept,slope,level,traitfun,vlist,binary,noise){
	sim = c(); count = 0; fcount = 0
	pc = perspectev.calc(data,traitfun,vlist)
	if(level=="lower"){
		while(count < simulations){
			temp = data
			survProb = inv.logit(slope * pc$lower[,2] + rnorm(nrow(temp), sd=noise) + intercept)
			if(binary==TRUE){
				temp$Survivorship = as.numeric(survProb >= runif(nrow(temp)))
			 	while(sum(temp$Survivorship) == 0 || sum(temp$Survivorship) == nrow(temp)){
			 		temp$Survivorship = as.numeric(survProb >= runif(nrow(temp)))}
			}else{temp$Survivorship = survProb}
			p = perspectev.test(temp,1,1,traitfun,vlist) 
			if(!is.na(p$correlation_permuted) & !is.na(p$correlation_observed) ){
				sim = rbind(sim,c(p$correlation_permuted,p$correlation_observed))
				count = count + 1}
			fcount = fcount + 1
			if(fcount >  100 & (fcount-count)/fcount > 0.9){
				stop('>90% of simulations yeild either no victims or no survivors')}
		}
	}else{
		while(count < simulations){
			temp = data; addOrder = c()
			survProb = inv.logit(slope * pc$upper[,2] + rnorm(nrow(pc$upper), sd=noise) + intercept)
			names(survProb) = rownames(pc$upper)
			Espec = sum(data$Survivorship)
			temp$Survivorship = 0
			uppernames = unique(temp$Upper)
			for(u in uppernames){
				st = data[data$Upper == u,]
				if(sum(st$Survivorship > 0)){
					surv = sample(st$Lower,size=1,prob=st$Survivorship)
				}else{
					surv = sample(st$Lower,size=1)
				}
				temp[surv,]$Survivorship = survProb[u]
			}
			Ospec = sum(temp$Survivorship)
			leftnames = temp[temp$Survivorship == 0,]$Lower
			leftprob = data[temp$Survivorship == 0,]$Survivorship
			if(length(leftnames[leftprob])){
				addOrder = sample(leftnames[leftprob > 0],size=length(leftnames[leftprob > 0]),prob=leftprob[leftprob > 0])
				addOrder = c(addOrder,sample(leftnames[leftprob == 0]))
			}else{addOrder=sample(leftnames)}
			counter = 1
			while(Ospec < Espec && counter <= length(addOrder)){
				spsurv = addOrder[counter]
				supper = temp[spsurv,]$Upper
				utemp = temp[temp$Upper == supper & temp$Survivorship != 0,]#get existing non-zero species in upper level
				nprob = 1 - (1-survProb[supper])^(1/(nrow(utemp)+1)) #recalculate survivorship probability 
				temp[utemp$Lower,]$Survivorship = nprob
				temp[spsurv,]$Survivorship = nprob
				counter = counter + 1
				Ospec = sum(temp$Survivorship)
			}
			if(binary == TRUE){
				sProb = temp$Survivorship
				temp$Survivorship = as.numeric(sProb >= runif(nrow(temp)))
			 	while(sum(temp$Survivorship) == 0 || sum(temp$Survivorship) == nrow(temp)){
			 		temp$Survivorship = as.numeric(sProb >= runif(nrow(temp)))}
			}
			p = perspectev.test(temp,1,1,traitfun,vlist) #occasionally get NA values from all or no genera surviving
			if(!is.na(p$correlation_observed) & !is.na(p$correlation_permuted) ){
				sim = rbind(sim,c(p$correlation_permuted,p$correlation_observed))
				count = count + 1}
			fcount = fcount + 1
			if(fcount >  100 & (fcount-count)/fcount > 0.9){
				stop('>90% of simulations yeild either no victims or no survivors')}
		}
	}
	return(sim)
}


perspectev.simulate <-
function(data,simulations,cores,traitfun=mcpRange,vlist=NULL,binary=NA,intercept=NA,slope=NA,level=NA,noise=0,fit=FALSE){
	if(is.na(intercept) || is.na(slope)){
		print('Intercept and slope parameters not specified. Assigning model fit through logistic regression.')
		fit=TRUE
	}
	if(is.na(level)){
		print("Level parameter not specified (can be 'lower' or 'upper'). Performing upper level simulation.")
		level = 'upper'
	}
	if(level != 'lower' && level != 'upper'){
		stop("Level parameter can only be 'lower' or 'upper'")
	}
	if(is.na(binary)){
		if(sum(data$Survivorship == 1) + sum(data$Survivorship == 0) == nrow(data)){
			binary=TRUE}else{binary=FALSE}
			print(paste("Binary/continuous surivorship not specified, setting binary =",binary))
	}
	g = perspectev.calc(data,traitfun,vlist)
	lower = g$lower
	upper = g$upper

	if(level == 'lower'){
		if(fit == TRUE){ #base parameters on fitting logistic regression to data
			#Note: in continuous case this involves fitting a binomial regression to 
			#non-binary data. Should give same slope and intercept.
			gl = suppressWarnings(glm(lower[,1] ~ lower[,2],family='binomial'))
			intercept = gl$coefficients[1]
			slope = gl$coefficients[2]
		}
		print(paste("Lower Level Selection Intercept",intercept))
		print(paste("Lower Level Selection Slope",slope))
	}else{
		if(fit == TRUE){
			gl = suppressWarnings(glm(upper[,1] ~ upper[,2],family='binomial'))
			intercept = gl$coefficients[1]
			slope = gl$coefficients[2]
		}
		print(paste("Upper Level Selection Intercept",intercept))
		print(paste("Upper Level Selection Slope",slope))
	}
	names(intercept) = NULL
	names(slope) = NULL
	if(cores > 1){
		div = floor(simulations/cores)
		plist = rep(div,length=cores)
		plist[length(plist)] = plist[length(plist)] + simulations %% cores
		cl = makeCluster(cores)
		registerDoParallel(cl)
		i=NULL
		sim = foreach(i = 1:cores,.packages=c('perspectev'),.combine="rbind") %dopar% perspectev.simulation.permute(data,plist[i],intercept,slope,level,traitfun,vlist,binary,noise)
		stopCluster(cl)
	}else{
		sim = perspectev.simulation.permute(data,simulations,intercept,slope,level,traitfun,vlist,binary,noise)
	}
	total = list()
	total$level = level
	total$intercept = intercept
	total$slope = slope
	total$noise = noise
	total$fitted_model = fit
	total$binary = binary
	total$pvalue = mean(sim[,1] >= sim[,2])
	total$correlation_permuted = sim[,1]
	total$correlation_observed = sim[,2]
	return(total)
}


perspectev.plot <-
function(observed,simulated,names,title){
	if(missing(observed)){
		stop("Must provide results from perspectev.test as 'observed'!")
	}
	if(missing(simulated)){
		simulated=list()
		names = c()
	}
	if(class(simulated) != 'list'){
		stop("'simulated' must be of type 'list'")
	}
	if(missing(names)){
		stop("Must provide vector of names to simulations using 'names' parameter!")
	}
	if(missing(title)){
		title = 'Difference between observed and permutated genera'
	}
	if(length(names) < 4){
		cbbPalette <- c("#2c7bb6","#abd9e9","#d7191c","#fdae61")
	}else if(length(names) < 9 && length(names) < 6){
		cbbPalette = c("#4575b4","#91bfdb","#e0f3f8","#fee090","#fc8d59","#d73027")
	}else if(length(names) < 9){
		cbbPalette <- c("#4575b4","#74add1","#abd9e9","#e0f3f8","#ffffbf","#fee090","#fdae61","#f46d43","#d73027")
	}else{
		stop("Cannot graph more than eight simulations at a time")
	}
	obs = observed$correlation_observed - observed$correlation_permuted
	tm = cbind('Observed',obs)
	for(l in 1:length(names)){
		temp1 = simulated[[l]]$correlation_observed -  simulated[[l]]$correlation_permuted
		temp2 = cbind(names[l],temp1)
		tm = rbind(tm,temp2)
	}
	tm = data.frame(tm)
	names(tm) = c('Simulation','RT')
	Simulation=RT=NULL
	tm$Simulation = factor(tm$Simulation,levels=c(names,'Observed'))
	tm$RT = as.numeric(as.character(tm$RT))

	ggplot(tm,aes(x=RT,fill=Simulation)) + 
		geom_density(alpha=0.8) + 
		xlab("R - S") +  ylab("Density") + 
		ggtitle(title) +
		guides(fill=guide_legend(title=NULL)) +
		scale_fill_manual(values = cbbPalette) +
		geom_vline(xintercept=0)
}
