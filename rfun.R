##' This function performs a combinatorial analysis, creating all
##' possible PHREEQC models with a number "len" of equilibrium phases
##' across a pool of primary and secondary minerals.
##'
##' Requires packages doParallel and foreach (for parallelization),
##' but they should already be installed as dependency of RedModRphree
##' @title Perform combinatorial analysis of a given pool of
##'     equilibrium minerals
##' @param initsol the base PHREEQC script
##' @param primary character vector containing minerals which are
##'     inserted in the models as primary minerals
##' @param secondary character vector containing minerals which are
##'     inserted in the models as secondary minerals (i.e., initial 0
##'     amount)
##' @param len number of phases in each model
##' @param procs integer, how many CPUs can we use? Defaults to 4
##' @return a list containing the parsed results (list of many blocks)
##' @author MDL

# The DoCombSim function is used here for performing combinatorial simulations.
#DoCombSim <- function(initsol, db, primary, secondary, len, procs = 4L): This line 
#defines a function named DoCombSim that takes several parameters - initsol, db, primary, 
#secondary, len, and procs, with a default value of 4 for procs.
#require(parallel), require(doParallel), require(foreach): These lines ensure that the 
#required packages (parallel, doParallel, foreach) are loaded. If not installed, these 
#packages need to be installed before using the function.

DoCombSim <- function(initsol, db, primary, secondary, len, procs=4L) {
  require(parallel)
  require(doParallel)
  require(foreach)
  
  ## Create all combinations of primary and secondary phases of length len
#  In this code, the variable "phases" is created by combining the "primary" 
#  and "secondary" vectors. The function combn is then used to generate all 
#  possible combinations of elements from "phases" with a specified length "len."
#  The result is a list of all combinations. The cat function is employed to print 
#  the number of simulations that will be conducted, which is determined by the 
#  length of the "combs" list
  phases <- c(primary, secondary)
  combs <- combn(phases, len, FUN = NULL, simplify=FALSE)
  
  cat(":: Going to do ", length(combs), " simulations\n")
  
  ## Create the phreeqc scripts.
## The function .addPhaseComb takes a vector x as an input and modifies the 
#  initial solution 'initsol' by adding properties related to the specified phases. 
#  For each phase in the input vector x, it uses the AddProp function to append 
#  properties to the solution. The properties added depend on whether the phase 
#  is in the "primary" or "secondary" category. If it's in the "primary" category,
#  the values are set to "0.0 2," and if it's in the "secondary" category, the 
#  values are set to "0.0 0." The function then returns the modified solution.
  .addPhaseComb <- function(x) {
    inp <- initsol
    for (phase in x) {
      inp <- AddProp(inp, name=phase, values=ifelse(phase %in% primary, "0.0 2", "0.0 0"), cat="pphases")
    }
    return(inp)
  }
  
  ## The lapply function is used to apply the .addPhaseComb function to each 
 # element (combination of phases) in the list combs. This results in a list 
 # biginp where each element is the modified solution obtained by adding properties 
 # based on the specific combination of phases.
  biginp <- lapply(combs, .addPhaseComb)
  
 # workhorse function to form a vector of all total elements + pH
 # The 'ExtractComponents' function in R is used here to extract and organize components 
 # from a linear model object (lin).
 # concs <- as.numeric(lin$tot[, 1]): Extract numeric concentrations from the first 
 # column ([, 1]) of the "tot" component of the linear model (lin).
 # names(concs) <- rownames(lin$tot): Assign names to the concentrations using row 
 # names from the "tot" component. This assumes that the row names are meaningful 
 # identifiers for the concentrations.
 # pH <- as.numeric(lin$desc["pH",1]): Extract the numeric pH value from the 
 # "desc" component for the "pH" row.
 # final <- c(concs, pH=pH): Combine concentrations and pH into a named
 # numeric vector (final).
  
  ExtractComponents <- function(lin) {
    concs <- as.numeric(lin$tot[,1])
    names(concs) <- rownames(lin$tot)
    pH <- as.numeric(lin$desc["pH",1])
    final <- c(concs, pH=pH)
    return(final)
  }
  
  ## Workhorse function to run simulations
 # The .runPQC function is a wrapper function that uses the phreeqc and RedModRphree packages 
 # in R to run a PHREEQC simulation.
 # phreeqc::phrSetOutputStringsOn(TRUE): This command turns on the output strings 
 # in the PHREEQC simulation, indicating that the simulation output will be captured.
 # phreeqc::phrRunString(input): This command runs the PHREEQC simulation using the 
 # provided input string. The input contains the PHREEQC script.
 # tmpout <- phreeqc::phrGetOutputStrings(): This line retrieves the output strings 
 # generated during the PHREEQC simulation.
 # res <- RedModRphree::ReadOut(tmpout)[[1]]: The output strings are processed using 
 # the RedModRphree::ReadOut function, and the result is extracted. The [[1]] indexing 
 # shows the first element of the output.
 # return(res): The function returns the result obtained from the PHREEQC simulation.
  
  .runPQC <- function(input) {
    phreeqc::phrSetOutputStringsOn(TRUE)
    phreeqc::phrRunString(input)
    tmpout <- phreeqc::phrGetOutputStrings()
    res <- RedModRphree::ReadOut(tmpout)[[1]]
    return(res)
  }
  
 # The function below checks the number of processors (procs) available. With more than 
 # one processor, it initializes a parallel computing cluster using the doParallel and parallel packages.
 # Depending on the operating system, it chooses between a PSOCK (Process Sockets) cluster 
 # or a Fork cluster for parallel processing.
 # It registers the parallel cluster using doParallel::registerDoParallel.
 # It loads a database using the phreeqc::phrLoadDatabase function on each worker in the 
 # cluster.
 # It uses the foreach::foreach loop to parallelize the computation for each element in the
 # biginp list, where .runPQC is a function that processes the input.
 # After the parallel computation, it stops the parallel cluster using 
 # parallel::stopCluster.
  #If there's only one processor available, it reverts to sequential computation 
#  using lapply.
  
  if (procs > 1) {
    if (Sys.info()[["sysname"]]=="Windows") {
      ThisRunCluster <- parallel::makePSOCKcluster(procs)
    } else {
      ThisRunCluster <- parallel::makeForkCluster(procs)
    }
    
    doParallel::registerDoParallel(ThisRunCluster)
    cat(":: Registered default doParallel cluster with ", procs, "nodes")
    parallel::clusterCall(cl=ThisRunCluster, phreeqc::phrLoadDatabase, db)
    msg(":: Database loaded on each worker")
    
    res <- foreach::foreach(i = seq_along(biginp)) %dopar% 
      .runPQC(biginp[[i]])
    cat("[ DONE ]\n")
    parallel::stopCluster(ThisRunCluster)
    
  } else {
    ## revert to sequential computation
    cat(":: Firing up PHREEQC onsingle CPU...")
    res <- lapply(biginp, .runPQC)
    cat("[ DONE ]\n")
    
  }
  
  return(res)
}

##' Computes a specific metric allowing for selection of the
##' components to be included
##'
##' @title Compute metric selecting components
##' @param data the matrix or data.frame containing all the results
##'     from the PHREEQC simulations. Its columns need to be named!
##' @param target the named vector with the target concentrations
##' @param FUN the name of the metric function
##' @param comp optional, a char vector with the names of the
##'     components. If unspecified, all components are selected
##' @param ... further parameter passed to FUN, such as "na.rm"
##' @return numeric vector with the computed metric
##' @author Marco

#The ComputeMetric function in R compute a specified metric 
#(e.g., root mean squared error - "rmse") for each row of a dataset (data) based on a 
#target vector (target).
#.Fun <- match.fun(FUN): The match.fun function is used to find the actual metric 
#function (FUN) to be applied. This provides flexibility to use different metrics.
#tmp <- subset(data, select = comp): The function selects only the columns specified 
#in the comp argument from the dataset (data).
#cvec <- target[comp]: It extracts the subset of the target vector (target) corresponding 
#to the specified columns.
#res <- apply(tmp, 1, function(x) .Fun(cvec, x, ...)): The specified metric function is 
#applied row-wise using the apply function. The result is a vector of computed metrics for 
#each row.
#return(res): The function returns the vector of computed metrics.

ComputeMetric <- function(data, target, FUN="rmse", comp=colnames(data), ...){
  ## find the metric function
  .Fun <- match.fun(FUN)
  
  ## retain only the columns given as argument
  tmp <- subset(data, select = comp)
  
  cvec <- target[comp]
  ## compute stuff using apply
  res <- apply(tmp, 1, function(x) .Fun(cvec, x, ...))
  return(res)
}

##### Filtering
#The Filter function in R is used here to filter a list of objects (lin) based on a 
#specified threshold (delta) for the absolute values of the pphases$delta component in 
#each object.
#Filter <- function(lin, delta = 0.5): This line defines a function named Filter that takes
#a list of objects (lin) and an optional parameter delta with a default value of 0.5.
#excluded <- sapply(lin, function(x) any(abs(x$pphases$delta) > delta)): This line uses 
#sapply to iterate over each object (x) in the list lin. For each object, it checks if 
#there exists any element in the absolute values of x$pphases$delta that is greater than 
#the specified delta. The result is a logical v'ector (excluded) indicating whether each 
#object should be excluded (TRUE) or not (FALSE).
#out <- lin[which(!excluded)]: This line uses 'which' to find the indices of objects that 
#are not excluded (where excluded is FALSE). It then subsets the original list lin to 
#include only those objects that are not excluded, and assigns the result to the variable out.
#return(out): Finally, the function returns the filtered list of objects (out).

Filter <- function(lin, delta=0.5) {
  excluded <- sapply(lin, function(x) any(abs(x$pphases$delta)>delta))
  out <- lin[which(!excluded)]
  return(out)
}

Filter2 <- function(lin, delta=0.5) {
  excluded <- sapply(lin, function(x) any(abs(x$pphases$delta)>delta))
  return(which(excluded))
}

FilterAll <- function(lin, delta=0.5) {
  retain <- sapply(lin, function(x) all(abs(x$pphases$delta)<delta))
  out <- which(!retain)
  return(out)
}

## Some metrics
rrmse <- function(y_true, y_pred, na.rm=TRUE)
  sqrt(mean(((y_true - y_pred)/y_true)^2, na.rm = na.rm))

## mean absolute percent error
rmape <- function(y_true, y_pred, na.rm=TRUE)
  mean(abs((y_true- y_pred)/y_true), na.rm = na.rm) * 100

##Additional errors relative Mean absolute error (mae)
rmae <- function(y_true, y_pred, na.rm=TRUE)
  mean(abs((y_true- y_pred)/y_true), na.rm = na.rm)


##The PlotComb function in R is used for visualizing and comparing compositional data 
#between two datasets (res and samples).
#comp: This argument allows the user to specify a subset of columns (comp) to be used for 
#plotting. If not provided, it defaults to the intersection of column names between res and samples.
#tmp <- subset(res, select=comp): The function subsets the relevant columns from both res and
#samples based on the specified or default set of columns.
#mins <- apply(sam, 2, min, na.rm=TRUE): The minimum (mins), maximum (maxs), median (meds),
#and mean (meas) values are calculated for each column in the samples dataset.
#colors <- heat.colors(nrow(tmp)): Colors for the barplot are generated using the heat.colors 
#function.
#out <- barplot(tmp, beside=TRUE, ylab="", log="y", col=colors, las=1, ...): The barplot 
#function is used to create a barplot for the subset of res. The bars are 
#grouped (beside=TRUE), and the y-axis is set to logarithmic scale (log="y").
#Overlaying Elements: For each column in the specified or default set, 
#rectangles (rect), dashed lines (segments) representing the median, and dotted 
#lines representing the mean are overlaid on the barplot.

PlotComb <- function(res, samples, comp, ...) {
  if (missing(comp)) {
    comp <- intersect(colnames(res), colnames(samples))
  }
  tmp <- subset(res, select=comp)
  sam <- subset(samples, select=comp)
  
  mins <- apply(sam, 2, min, na.rm=TRUE)
  maxs <- apply(sam, 2, max, na.rm=TRUE)
  meds <- apply(sam, 2, median, na.rm=TRUE)
  meas <- apply(sam, 2, mean, na.rm=TRUE)
  
  colors <- heat.colors(nrow(tmp))
  
  out <- barplot(tmp, beside=TRUE, ylab="", log="y",
                 col=colors, las=1, ...)
  for (i in seq_along(comp)) {
    rect(out[1,i]-0.6, mins[i],out[nrow(out),i]+0.6, maxs[i],col=rgb(0,0,1.0,alpha=0.5))
    segments(out[1,i]-0.6, meds[i], out[nrow(out),i]+0.6, meds[i],col="red", lwd=2, lty="dashed")
    segments(out[1,i]-0.6, meas[i], out[nrow(out),i]+0.6, meas[i],col="grey", lwd=2, lty="dotted")
  }
}



