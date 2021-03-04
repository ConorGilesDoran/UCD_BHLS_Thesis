# Final Year Project - Signature extrapolation functions
# Bioinformatic Analysis of Copy Number Alteration Signature Mechanisms in Ovarian Carcinoma
# Function code from - "Copy number signatures and mutational processes in ovarian carcinoma", Macintyre et al 2018.
# Code available - https://bitbucket.org/britroc/cnsignatures/src/master/
# Conor Giles Doran - UCD - Sept. 2020 - Feb. 2021

# Functions required for signature extrapolation script - 'sig_extrapolation.R'

#### HELPER FUNCTIONS ####

# getSegtable
# Used within each get-*Feature* function below
# Returns segment tables (one for each sample) with the following column headers: "chromosome", "start", "end", "segVal"
# Argument - QDNAseq class object
# Uses Biobase package to extract assay and feature Data of QDNAseq class

getSegTable <- function(x)
{
  dat <- x
  sn <- Biobase::assayDataElement(dat,"segmented")
  fd  <-  Biobase::fData(dat)
  # accessing the 'use' data, assigned to 'use' object. TRUE/FALSE....is this sample being used for analysis or not
  fd$use -> use
  # fd - feature data, accessing just row of data actually being used for analysis
  fdfiltfull <- fd[use,]
  # "segmented" matrix iterated over to now become just the data being used
  sn <- sn[use,]
  segTable <- c()
  for(c in unique(fdfiltfull$chromosome))
  {
    snfilt <- sn[fdfiltfull$chromosome==c]
    fdfilt <- fdfiltfull[fdfiltfull$chromosome==c,]
    sn.rle <- rle(snfilt)
    starts  <-  cumsum(c(1, sn.rle$lengths[-length(sn.rle$lengths)]))  # identify start and end points for each chromosome segment
    ends  <-  cumsum(sn.rle$lengths)
    lapply(1:length(sn.rle$lengths), function(s) {
      from  <-  fdfilt$start[starts[s]]
      to  <-  fdfilt$end[ends[s]]
      segValue  <-  sn.rle$value[s]
      c(fdfilt$chromosome[starts[s]], from, to, segValue)            # segment copy number included
    }) -> segtmp
    segTableRaw  <-  data.frame(matrix(unlist(segtmp), ncol=4, byrow=T),stringsAsFactors=F)
    segTable <- rbind(segTable,segTableRaw)
  }
  colnames(segTable)  <-  c("chromosome", "start", "end", "segVal")
  segTable
}


# getSampNames
# Returns the name of each sample from the QDNAseq class/segment tables data (abs_profiles)

# abs_profiles = CN_data/list of segment tables
getSampNames <- function(abs_profiles){
  # if CN profiles in QDNASeqCN class
  if(class(abs_profiles)=="QDNAseqCopyNumbers")
  {
    # assign column names (sample names) to 'samps'
    samps <- colnames(abs_profiles)
  }
  # otherwise
  # assign the names of sample CN profiles to 'samps' (if in segment tables perhaps?)
  else
  {
    samps <- names(abs_profiles)
  }
  return(samps)
}

# getPloidy
# used for preprocessing to identify ploidy of each chromosome segment
# uses start and endpoint data as well as copy number of each segment
# abs_profiles - QDNAseq object class or list of segment tables

getPloidy  <-  function(abs_profiles){
  out <- c()
  samps <- getSampNames(abs_profiles)
  for(i in samps)
  {
    if(class(abs_profiles)=="QDNAseqCopyNumbers")
    {
      segTab <- getSegTable(abs_profiles[,which(colnames(abs_profiles)==i)])
    }
    else
    {
      segTab <- abs_profiles[[i]]
      colnames(segTab)[4] <- "segVal"
      segTab <- segTab[!is.na(segTab$segVal),]
    }
    segLen <- (as.numeric(segTab$end)-as.numeric(segTab$start))
    ploidy <- sum((segLen/sum(segLen))*as.numeric(segTab$segVal))
    out <- c(out,ploidy)
  }
  data.frame(out,stringsAsFactors = F)
}


# normaliseMatrix
# used by the function - 'quantifySignatures' for various forms of normalisation

normaliseMatrix <- function(signature_by_sample,sig_thresh=0.01)
{
  norm_const <- colSums(signature_by_sample)                                     # sum of columns
  sample_by_signature <- apply(signature_by_sample,1,function(x){x/norm_const})  # values in rows divided by sum of column, 1 indicates rows
  sample_by_signature <- apply(sample_by_signature,1,lower_norm,sig_thresh)      # assigning 0 values using below function
  signature_by_sample <- t(sample_by_signature)                                  # transpose
  norm_const <- apply(signature_by_sample,1,sum)                                 # sum of rows
  sample_by_signature <- apply(signature_by_sample,2,function(x){x/norm_const})  # values in columns divided by sum of rows
  signature_by_sample <- t(sample_by_signature)                                  # trasnspose again
  signature_by_sample
}

# lower-norm - used in above function - 'normaliseMatrix'
# assigning 0 to values less than 0.01

lower_norm <- function(x,sig_thresh=0.01)
{
  new_x <- x
  for(i in 1:length(x))
  {
    if(x[i]<sig_thresh)
    {
      new_x[i] <- 0
    }
  }
  new_x
}


#### COPY NUMBER FEATURES ####

# getSegsize
# Argument(abs_profiles) input as QDNAseq object or segment tables
# getSampNames and getSegTable used within
# Uses start and end points of chromosome segments to determine length of each segment
# Returns table with sample ID and corresponding length of each chromosome segment within
# Used within the extractCopynumberFeatures function

getSegsize <- function(abs_profiles){
  out <- c()
  samps <- getSampNames(abs_profiles)
  for(i in samps)
  {
    if(class(abs_profiles)=="QDNAseqCopyNumbers")
    {
      segTab <- getSegTable(abs_profiles[,which(colnames(abs_profiles)==i)])
    }
    else
    {
      segTab <- abs_profiles[[i]]
      colnames(segTab)[4] <- "segVal"
      segTab <- segTab[!is.na(segTab$segVal),]
    }
    segTab$segVal[as.numeric(segTab$segVal)<0] <- 0
    seglen <- (as.numeric(segTab$end)-as.numeric(segTab$start))
    seglen <- seglen[seglen>0]
    out <- rbind(out,cbind(ID=rep(i,length(seglen)),value=seglen))
    # giving a segment ID to each value of segment length data for same sample ID (hence - rep())
  }
  rownames(out) <- NULL
  data.frame(out,stringsAsFactors = F)
}


# getBPnum - Breakpoint count per 10 megabases
# Argument(abs_profiles) input as QDNAseq object or segment tables
# chrlen - chromosome lengths data
# getSampNames and getSegTable used within
# Divides chromosome length data into 10MB intervals and counts number of chromsome segments in each
# Returns table with sample ID and corresponding breakpoint count per 10MB.
# Used within the extractCopynumberFeatures function

getBPnum <- function(abs_profiles,chrlen)
{
  out <- c()
  samps <- getSampNames(abs_profiles)
  for(i in samps)
  {
    if(class(abs_profiles)=="QDNAseqCopyNumbers")
    {
      segTab <- getSegTable(abs_profiles[,which(colnames(abs_profiles)==i)])
    }else
    {
      segTab <- abs_profiles[[i]]
      colnames(segTab)[4] <- "segVal"
      segTab <- segTab[!is.na(segTab$segVal),]
    }
    chrs <- unique(segTab$chromosome)               # chromosomes 1-22
    allBPnum <- c()                                 # empty vector called 'allBPnum'
    for(c in chrs)                                  # for each chromosome
    {
      print(paste0("chr ", c))
      currseg <- segTab[segTab$chromosome==c,]
      # all data for each chromosome , first time it will be all chr 1 data...and so on

      intervals <- seq(1,chrlen[chrlen[,1]==paste0("chr",c),2]+10000000,10000000)
      # creating the 10MB intervals, sequence from 1-length of the chromosome +10MB, interval then set to 10MB
      # paste0() just sets the sep = "" so its chr1, chr2....
      currsegmax <- max(currseg$end[-nrow(currseg)])
      if(currsegmax > rev(intervals)[1]){
        intervals <- c(intervals, max(intervals)+10000000)
      }
      res  <-  hist(as.numeric(currseg$end[-nrow(currseg)]),breaks=intervals,plot=FALSE)$counts
      # does not plot histogram, uses data from current chromosome segment of for loop,
      # removes the end point, breaks are the 10MB intervals, counts the breaks in each 10MB interval

      allBPnum <- c(allBPnum,res)
      # stores this data in allBPnum object
    }
    out <- rbind(out,cbind(ID=rep(i,length(allBPnum)),value=allBPnum))
      # creates column in 'out' for sample ID and corresponding BPnum
  }
  rownames(out) <- NULL                   # no row names
  data.frame(out,stringsAsFactors = F)  # convert to data frame
}


# getOscillation - lengths of chains of oscillating copy number
# Argument(abs_profiles) input as QDNAseq object or segment tables
# getSampNames and getSegTable used within
# Uses rounded copy number data to identify oscillating copy number states
# Returns table with sample ID and corresponding oscillating copy number chain length
# Used within the extractCopynumberFeatures function

getOscilation <- function(abs_profiles,chrlen){
  out <- c()
  samps <- getSampNames(abs_profiles)
  for(i in samps)
  {
    if(class(abs_profiles)=="QDNAseqCopyNumbers")
    {
      segTab <- getSegTable(abs_profiles[,which(colnames(abs_profiles)==i)])
    }else
    {
      segTab <- abs_profiles[[i]]
      colnames(segTab)[4] <- "segVal"
      segTab <- segTab[!is.na(segTab$segVal),]
    }
    chrs <- unique(segTab$chromosome)                    # all the same as previous^
    oscCounts <- c()
    for(c in chrs)
    {
      currseg <- segTab[segTab$chromosome==c,"segVal"]   # seg values for each chromosome
      currseg <- round(as.numeric(currseg))              # seg values rounded
      if(length(currseg)>3)                              # if seg value vector has > 3 values
      {
        prevval <- currseg[1]                            # prevval is the first element of seg values
        count=0                                          # start the count at 0
        for(j in 3:length(currseg))
        {
          if(currseg[j]==prevval&currseg[j]!=currseg[j-1])
            # First Iter: if 3rd element is == to 1st element and not equal previous element (2nd), count it
            # otherwise the oscillating counts == 0
          {
            count <- count+1
          }else{
            oscCounts <- c(oscCounts,count)
            count=0
          }
          prevval <- currseg[j-1]
           # prevval reset to 2 elements before current seg value
        }
      }
    }
    # add it to out vector, sample ID for each osc count
    out <- rbind(out,cbind(ID=rep(i,length(oscCounts)),value=oscCounts))  
    if(length(oscCounts)==0)                 # if it has no counts
    {
      out <- rbind(out,cbind(ID=i,value=0))  # assign ID to value of 0
    }
  }
  rownames(out) <- NULL
  data.frame(out,stringsAsFactors = F)
}


# getCentromereDistCounts - Breakpoint count per chromosome arm
# Argument(abs_profiles) input as QDNAseq object or segment tables
# getSampNames and getSegTable used within
# Uses centromere and chromosome length data to identify number of segments per chromosome arm - p and q
# Returns table with sample ID and corresponding breakpoint count per chromosome arm
# Used within the extractCopynumberFeatures function

getCentromereDistCounts <- function(abs_profiles,centromeres,chrlen)
{
  out <- c()
  samps <- getSampNames(abs_profiles)
  for(i in samps)
  {
    if(class(abs_profiles)=="QDNAseqCopyNumbers")
    {
      segTab <- getSegTable(abs_profiles[,which(colnames(abs_profiles)==i)])
    }else
    {
      segTab <- abs_profiles[[i]]
      colnames(segTab)[4] <- "segVal"      # same as above^
      segTab <- segTab[!is.na(segTab$segVal),]
    }
    chrs <- unique(segTab$chromosome)      # chromosomes 1:22
    all_dists <- c()                       # empty vector
    for(c in chrs)                         # for each chromosome
    {
      if(nrow(segTab)>1)                   # if the no. of rows in the segment table > 1
      {
        starts <- as.numeric(segTab[segTab$chromosome==c,2])[-1]   # all start data for that number chromosome [minus the first value]
        segstart <- as.numeric(segTab[segTab$chromosome==c,2])[1]  # segment start point is the first seg value for that no. chromosome
        ends <- as.numeric(segTab[segTab$chromosome==c,3])         # ends - all end data(col 3) for that number chromosome
        segend <- ends[length(ends)]                               # seg end is the last value of ends
        ends <- ends[-length(ends)]                                # get rid of last value in ends, as it has been assigned to segend object

        centstart <- as.numeric(centromeres[substr(centromeres[,2],4,5)==c,3])
        # centromere start - column 3 (chromstart) data according to col 2 (chrom) - chromosome number (c)

        centend <- as.numeric(centromeres[substr(centromeres[,2],4,5)==c,4])
        # centromere end - column 4 (chromeENd) data according to col 2 (chrom) - chromosome number (c)

        chrend <- chrlen[substr(chrlen[,1],4,5)==c,2]
        # chromosome end data from column 2 of chrlen data according to col 1 - chromosome number (c)

        ndist <- cbind(rep(NA,length(starts)),rep(NA,length(starts)))
        # creating row and column structure for all distance data

        ndist[starts<=centstart,1]  <-  (centstart-starts[starts<=centstart])/(centstart-segstart)*-1
        # start values < = centromere start is (centromere start - segment start points that are less than centstart) divided by
        # (centromere start minus the segment start value) multiplied by -1

        ndist[starts>=centend,1]  <-  (starts[starts>=centend]-centend)/(segend-centend)
        # start values > = centromere end is ((starts > = centromere end) - centromere end point) divided by
        # (segment end value - centromere end)

        ndist[ends<=centstart,2]  <-  (centstart-ends[ends<=centstart])/(centstart-segstart)*-1
        # same but for column 2 and using 'ends' values
        ndist[ends>=centend,2]  <-  (ends[ends>=centend]-centend)/(segend-centend)

        ndist <- apply(ndist,1,min)                  # getting minimum values

        all_dists <- rbind(all_dists,sum(ndist>0))   # sum of all distances > 0
        all_dists <- rbind(all_dists,sum(ndist<=0))  # sum of all distances <=0
      }
    }
    if(nrow(all_dists)>0)  # if number of rows >0, bind to the out dataframe, with IDs and ct1
    {
      out <- rbind(out,cbind(ID=i,ct1=all_dists[,1]))
    }
  }
  rownames(out) <- NULL
  data.frame(out,stringsAsFactors = F)
}


# getChangepointCN -  the absolute difference in copy number between adjacent segments across the genome
# Argument(abs_profiles) input as QDNAseq object or segment tables
# getSampNames and getSegTable used within
# Uses segment copy number data to identify copy number changes
# Returns table with sample ID and corresponding changepoint
# Used within the extractCopynumberFeatures function

getChangepointCN <- function(abs_profiles)
{
  out <- c()
  samps <- getSampNames(abs_profiles)
  for(i in samps)
  {
    if(class(abs_profiles)=="QDNAseqCopyNumbers")
    {
      segTab <- getSegTable(abs_profiles[,which(colnames(abs_profiles)==i)])
    }
    else
    {
      segTab <- abs_profiles[[i]]
      colnames(segTab)[4] <- "segVal"                                 # ^all same as above^
      segTab <- segTab[!is.na(segTab$segVal),]
    }
    segTab$segVal[as.numeric(segTab$segVal)<0] <- 0                   # any seg values < 0 are given value of 0
    chrs <- unique(segTab$chromosome)                                 # 1:22
    allcp <- c()                                                      # empty vector
    for(c in chrs)                                                    # for chr number in 1:22
    {
      currseg <- as.numeric(segTab[segTab$chromosome==c,"segVal"])    # all segvalue data for each chromosome , first time it will be all chr 1 data...and so on
      allcp <- c(allcp,abs(currseg[-1]-currseg[-length(currseg)]))    # subtracts first from second, second from third, third from fourth....and so on
    }
    if(length(allcp)==0)
    {
      allcp <- 0 # if there are no changepoints
    }
    out <- rbind(out,cbind(ID=rep(i,length(allcp)),value=allcp))
  }
  rownames(out) <- NULL
  data.frame(out,stringsAsFactors = F)
}


# getCN -  the absolute copy number of chromosome segments
# Argument(abs_profiles) input as QDNAseq object or segment tables
# getSampNames and getSegTable used within
# Extracts copy number values from segment copy number data
# Returns table with sample ID and corresponding copy number
# Used within the extractCopynumberFeatures function

getCN <- function(abs_profiles)
{
  out <- c()
  samps <- getSampNames(abs_profiles)
  for(i in samps)
  {
    if(class(abs_profiles)=="QDNAseqCopyNumbers")
    {
      segTab <- getSegTable(abs_profiles[,which(colnames(abs_profiles)==i)])
    }
    else
    {
      segTab <- abs_profiles[[i]]
      colnames(segTab)[4] <- "segVal"                       # ^all same as above^
      segTab <- segTab[!is.na(segTab$segVal),]
    }
    segTab$segVal[as.numeric(segTab$segVal)<0] <- 0         # any seg values < 0 are given value of 0
    cn <- as.numeric(segTab$segVal)                         # copy number = segvalues for each segment
    out <- rbind(out,cbind(ID=rep(i,length(cn)),value=cn))  # added to out with sample ID for each value
  }
  rownames(out) <- NULL
  data.frame(out,stringsAsFactors = F)
}

##below is example of how an R package function is annotated (in 'roxygen2' format):
##also addded some print statements to see how things are going in there!

#' extractCopyNumberfeatures - used to derive genome wide CN Feature Distributions within the data
#' @param CN_data - QDNAseq object class or list of segment tables
#' @param chrom_size_file file showing chromosome sizes
#' @param chrom_gap_file file showing chromosome gaps
#' @return a list of copy-number features extracted from samples

extractCopynumberFeatures <- function(CN_data, chrom_size_file, chrom_gap_file, cores = 1){
  #get chromosome lengths
  chrlen <- read.table(chrom_size_file, sep="\t",stringsAsFactors = F)[1:24,]
  #get centromere locations
  gaps <- read.table(chrom_gap_file, sep="\t", header = F, stringsAsFactors = F)
  centromeres <- gaps[gaps[,8]=="centromere",]

  if(cores > 1) {
    require(foreach)  # loads package - foreach
    doMC::registerDoMC(cores)

    temp_list = foreach::foreach(i=1:6) %dopar% {
      if(i == 1){
        print("Get Segment Size")
        list(segsize = getSegsize(CN_data) )
      } else if (i == 2) {
        print("Get Breakpoint Number")
        list(bp10MB = getBPnum(CN_data,chrlen) )
      } else if (i == 3) {
        print("Get Oscillations")
        list(osCN = getOscilation(CN_data,chrlen) )
      } else if (i == 4) {
        print("Get Centromere Distance Counts")
        list(bpchrarm = getCentromereDistCounts(CN_data,centromeres,chrlen) )
      } else if (i == 5) {
        print("Get Changepoint Copynumbers")
        list(changepoint = getChangepointCN(CN_data) )
      } else {
        print("Get Copynumbers")
        list(copynumber = getCN(CN_data) )
      }

    }
    unlist( temp_list, recursive = FALSE )
  } else {
    print("Get Segment Size")
    segsize <- getSegsize(CN_data)
    print("Get Breakpoint Number")
    bp10MB <- getBPnum(CN_data,chrlen)
    print("Get Oscillations")
    osCN <- getOscilation(CN_data,chrlen)
    print("Get Centromere Distance Counts")
    bpchrarm <- getCentromereDistCounts(CN_data,centromeres,chrlen)
    print("Get Changepoint Copynumbers")
    changepoint <- getChangepointCN(CN_data)
    print("Get Copynumbers")
    copynumber <- getCN(CN_data)

    list(segsize=segsize,bp10MB=bp10MB,osCN=osCN,bpchrarm=bpchrarm,changepoint=changepoint,copynumber=copynumber)
  }

}


# MIXTURE MODELLING ----

# Used to breakdown each feature distribution into mixtures of Gaussian or Poisson distributions using the flexmix package

# two functions mentioned

# fitcomponent
# can be used on each individual feature and all compiled at the end

# fitMixtureModels
# used to do all the individual fitComponents^ within just one function

# Component paramaters will be used to generate sample by component matrix and extract signatures

# fitComponent

fitComponent <- function(CN_feats, dist="norm", seed=77777, model_selection="BIC", min_prior=0.001, niter=1000, nrep=1, min_comp=2, max_comp=10){
  control  <-  new("FLXcontrol")  # creates s4 object of class FLXcontrol
  control@minprior  <-  min_prior # assign values stated in arguments ^ to minprior and iter.max slots of FLXcontrol class
  control@iter.max  <-  niter

  set.seed(seed)
  # seed - start point used for generation of sequence of random numbers
  # will obtain the same numbers because same seed specified
  # psuedorandom numbers are the seeds you run through the mathematical function

  if(dist=="norm")
    # if the distribution is normal
  {
    if(min_comp == max_comp)
      # and min components is the same as max components
    {
      fit <- flexmix::flexmix(CN_feats ~ 1, model=flexmix::FLXMCnorm1(), k=min_comp, control=control)
      # implementing model-based clustering of Gaussian data
      # dat is response/data to be clustered, 1 is the number of predictors, model- FLXMCnorm1() for the components of the model
      # k is Number of clusters we want, control class from above

    }else{
      fit <- flexmix::stepFlexmix(CN_feats ~ 1,model = flexmix::FLXMCnorm1(), k=min_comp:max_comp, nrep=nrep, control=control)
      # Runs flexmix repeatedly for different numbers of components and returns the maximum likelihood solution for each.

      fit <- flexmix::getModel(fit,which=model_selection)
      # Get single model from a collection of models
      # 'BIC' = Bayesian information criterion
    }

  }else if(dist=="pois")
    # if the distribution is poisson
  {
    if(min_comp==max_comp)
      # and min components is the same as max components?
    {
      fit <- flexmix::flexmix(CN_feats ~ 1,model=flexmix::FLXMCmvpois(),k=min_comp,control=control)
      # implementing model-based clustering of Poisson distributed data.
    }else{
      fit <- flexmix::stepFlexmix(CN_feats ~ 1,model = flexmix::FLXMCmvpois(),k=min_comp:max_comp,nrep=nrep,control=control)
      fit <- flexmix::getModel(fit,which=model_selection)
    }
  }
  fit
}


# fitMixtureModels - same as above but all in one function

fitMixtureModels <- function(CN_feats, seed=77777, min_comp=2, max_comp=10, min_prior=0.001, model_selection="BIC", nrep=1, niter=1000, cores = 1, featsToFit = seq(1, 6))
{

  if(cores > 1) {
    require(foreach)

    doMC::registerDoMC(cores)

    temp_list = foreach(i=1:6) %dopar% {

      if(i == 1 & i %in% featsToFit ){

        dat <- as.numeric(CN_feats[["segsize"]][,2])
        list( segsize = fitComponent(dat,seed=seed,model_selection=model_selection,
                                     min_prior=min_prior,niter=niter,nrep=nrep,min_comp=min_comp,max_comp=max_comp) )

      } else if (i == 2 & i %in% featsToFit ) {

        dat <- as.numeric(CN_feats[["bp10MB"]][,2])
        list( bp10MB = fitComponent(dat,dist="pois",seed=seed,model_selection=model_selection,
                                    min_prior=min_prior,niter=niter,nrep=nrep,min_comp=min_comp,max_comp=max_comp) )

      } else if (i == 3 & i %in% featsToFit ) {

        dat <- as.numeric(CN_feats[["osCN"]][,2])
        list( osCN = fitComponent(dat,dist="pois",seed=seed,model_selection=model_selection,
                                  min_prior=min_prior,niter=niter,nrep=nrep,min_comp=min_comp,max_comp=max_comp) )

      } else if (i == 4 & i %in% featsToFit ) {

        dat <- as.numeric(CN_feats[["bpchrarm"]][,2])
        list( bpchrarm = fitComponent(dat,dist="pois",seed=seed,model_selection=model_selection,
                                      min_prior=min_prior,niter=niter,nrep=nrep,min_comp=min_comp,max_comp=max_comp) )

      } else if (i == 5 & i %in% featsToFit ) {

        dat <- as.numeric(CN_feats[["changepoint"]][,2])
        list( changepoint = fitComponent(dat,seed=seed,model_selection=model_selection,
                                         min_prior=min_prior,niter=niter,nrep=nrep,min_comp=min_comp,max_comp=max_comp) )

      } else if (i == 6 & i %in% featsToFit) {

        dat <- as.numeric(CN_feats[["copynumber"]][,2])
        list( copynumber = fitComponent(dat,seed=seed,model_selection=model_selection,
                                        nrep=nrep,min_comp=min_comp,max_comp=max_comp,min_prior=0.005,niter=2000) )

      }

    }
    unlist( temp_list, recursive = FALSE )
  } else {
    print("Working on segsize")
    dat <- as.numeric(CN_feats[["segsize"]][,2])
    segsize_mm <- fitComponent(dat,seed=seed,model_selection=model_selection,
                             min_prior=min_prior,niter=niter,nrep=nrep,min_comp=min_comp,max_comp=max_comp)

    print("Working on bp10MB")
    dat <- as.numeric(CN_feats[["bp10MB"]][,2])
    bp10MB_mm <- fitComponent(dat,dist="pois",seed=seed,model_selection=model_selection,
                            min_prior=min_prior,niter=niter,nrep=nrep,min_comp=min_comp,max_comp=max_comp)
    print("Working on osCN")
    dat <- as.numeric(CN_feats[["osCN"]][,2])
    osCN_mm <- fitComponent(dat,dist="pois",seed=seed,model_selection=model_selection,
                          min_prior=min_prior,niter=niter,nrep=nrep,min_comp=min_comp,max_comp=max_comp)
    print("Working on bpchrarm")
    dat <- as.numeric(CN_feats[["bpchrarm"]][,2])
    bpchrarm_mm <- fitComponent(dat,dist="pois",seed=seed,model_selection=model_selection,
                              min_prior=min_prior,niter=niter,nrep=nrep,min_comp=min_comp,max_comp=max_comp)
    print("Working on changepoint")
    dat <- as.numeric(CN_feats[["changepoint"]][,2])
    changepoint_mm <- fitComponent(dat,seed=seed,model_selection=model_selection,
                                 min_prior=min_prior,niter=niter,nrep=nrep,min_comp=min_comp,max_comp=max_comp)
    print("Working on copynumber")
    dat <- as.numeric(CN_feats[["copynumber"]][,2])
    copynumber_mm <- fitComponent(dat,seed=seed,model_selection=model_selection,
                                nrep=nrep,min_comp=min_comp,max_comp=max_comp,min_prior=0.005,niter=2000)

    list(segsize=segsize_mm,bp10MB=bp10MB_mm,osCN=osCN_mm,bpchrarm=bpchrarm_mm,changepoint=changepoint_mm,copynumber=copynumber_mm)
  }
}


# SAMPLE_BY_COMPONENT ----

# calculateSumOfPosteriors
# takes CN_features list and calculated components
# computes posterior probability of a CN feature property belonging to a component
# posterior probability vectors summed
# used within the generateSampleByComponentMatrix function

calculateSumOfPosteriors  <-  function(CN_feature, components, name, rowIter = 1000, cores = 1)
{

  if(cores > 1){
    require(foreach)
    require(doMC)  # all for parallel processing, windows doesn't have this, maybe try alternative

    len = dim(CN_feature)[1]           # number of rows of data in that feature
    iters = floor( len / rowIter )     # rounds the number - number of rows divided by 1000 iter
    lastiter = iters[length(iters)]    # last iteration

    registerDoMC(cores)                # for parallel processing^
    curr_posterior  <-  foreach( i=0:iters, .combine=rbind) %dopar% {
      start = i*rowIter+1
      if(i != lastiter) { end = (i+1)*rowIter } else { end = len }
      flexmix::posterior(components,data.frame(dat=as.numeric(CN_feature[start:end,2])))
    }

  } else {  # do it slower
    curr_posterior <- flexmix::posterior(components,data.frame(dat=as.numeric(CN_feature[,2])))
  }

  mat <- cbind(CN_feature,curr_posterior)    # create matrix of feature and current posterior values
  posterior_sum <- c()                       # emty vector for sum of posterior probabilities

  for(i in unique(mat$ID))                 # for each unique sample id
  {
    posterior_sum <- rbind(posterior_sum,colSums(mat[mat$ID==i,c(-1,-2)]))  # explained below
  }
  params <- flexmix::parameters(components)  # gets the paramters for a model
  if(!is.null(nrow(params)))               # if there are no rows within the parameters (=NULL)
  {
    posterior_sum <- posterior_sum[,order(params[1,])] # reorder
  }
  else
  {
    posterior_sum <- posterior_sum[,order(params)]
  }
  colnames(posterior_sum) <- paste0(name,1:ncol(posterior_sum))  # assign names to each  of the columns and rows
  rownames(posterior_sum) <- rownames(unique(mat$ID))
  posterior_sum
}

# generateSampleByComponentMatrix
# computes posterior probability of a CN feature property belonging to a component
# posterior probability vectors summed and compiled into a sample by component matrix
# can use the

generateSampleByComponentMatrix <- function(CN_features, comp_params_rds = NULL, all_components=NULL, cores = 1, rowIter = 1000, subcores = 2)
{
  if(is.null(all_components) & is.null(comp_params_rds)){
    print("Must include one of all_components, comp_params_rds")
    break
  }
  if(is.null(all_components))
  {
    all_components  <-  readRDS(comp_params_rds)

  }

  # do not have to specify component parameters, creating components was an optional step, they can provide this data if needed

  if(cores > 1){
    require(foreach) # loads package - foreach - looping construct for parallel processing

    feats = c( "segsize", "bp10MB", "osCN", "changepoint", "copynumber", "bpchrarm" )   # the features
    doMC::registerDoMC(cores)                                                           # not available on windows, but should work fine anyway

    full_mat  <-  foreach(feat = feats, .combine = cbind) %dopar% {
      calculateSumOfPosteriors(CN_features[[feat]],all_components[[feat]],              # calculatesumofposteriors all in one
                               feat, rowIter = rowIter, cores = subcores)               # and columnbind into matrix
    }
  } else {                                                                              # or individually
    full_mat  <-  cbind(
      calculateSumOfPosteriors(CN_features[["segsize"]],all_components[["segsize"]],"segsize"),
      calculateSumOfPosteriors(CN_features[["bp10MB"]],all_components[["bp10MB"]],"bp10MB"),
      calculateSumOfPosteriors(CN_features[["osCN"]],all_components[["osCN"]],"osCN"),
      calculateSumOfPosteriors(CN_features[["changepoint"]],all_components[["changepoint"]],"changepoint"),
      calculateSumOfPosteriors(CN_features[["copynumber"]],all_components[["copynumber"]],"copynumber"),
      calculateSumOfPosteriors(CN_features[["bpchrarm"]],all_components[["bpchrarm"]],"bpchrarm"))
  }

  rownames(full_mat)  <-  unique(CN_features[["segsize"]][,1])
  full_mat[is.na(full_mat)]  <-  0
  full_mat
}



# SIGNATURE EXTRACTION ----

# chooseNumberSignatures
# used to detect molecular patterns within the sample_by_comp matrix data

# signature search interval  <-  [3,12]

# matrix factorisation ran 1000 times for each number of signatures, each with a different random seed

# cophenetic, dispersion, silhouette, and sparseness scores compared for the signature-feature matrix (basis),
# patient-signature matrix (coefficients) and a consensus matrix of patient-by-patient across the 1000 runs.

# 1000 random shuffles of the input matrix performed to get a null estimate of each of the scores
# (Macintyre et al. 2018)

chooseNumberSignatures <- function(sample_by_component, outfile = NULL, min_sig=3, max_sig=12, iter=100, cores=1)
{

  if(is.null(outfile)){
    outfile <- "sigNums.pdf"
  }

  nmfalg <- "brunet"  # specific algorithm
  seed <- 77777       # seed specified

  # estimate rank (signatures) - transposed samp_by_comp, min sigs = 3, max = 12^
  # It defines the number of basis effects used to approximate the target matrix
  estim.r  <-  NMF::nmfEstimateRank(t(sample_by_component), min_sig:max_sig,seed = seed,nrun=iter,
                                  verbose=FALSE, method=nmfalg, .opt=list(shared.memory=FALSE, paste0("p", cores) ) )

  V.random  <-  NMF::randomize(t(sample_by_component))  # used to randomize the data (shuffles?)

  # estimate rank of the random data^
  estim.r.random  <-  NMF::nmfEstimateRank(V.random, min_sig:max_sig, seed =seed,nrun=iter,
                                         verbose=FALSE, method=nmfalg, .opt=list(shared.memory=FALSE, paste0("p", cores) ) )

  # comparing the 4 different scores outlined in explanation above^ in both the non random estimation and random estimation
  p <- NMF::plot(estim.r,estim.r.random,
               what = c("cophenetic", "dispersion","sparseness", "silhouette"),
               xname="Observed",yname="Randomised",main="")

  pdf(file=outfile, width=10, height=10 )  # save to pdf
  p
  dev.off()

  return(p)

}


# generateSignatures
# generates the component by signature matrix
# carries out non-negative matrix factorization using sample_by_comp
# must specify number of signatures - optimal number identified from chooseNumberSignatures output
# coefmap and basismap functions to be used in order to generate heatmaps

generateSignatures <- function(sample_by_component, nsig, seed=77777, nmfalg="brunet", cores=1)
{
  NMF::nmf(t(sample_by_component), nsig, seed=seed, nrun=1000, method=nmfalg, .opt = paste0("p", cores) )
}

# quantifySignatures
# quantifies signature exposures using the LCD function from the YAPSA package,
# returning a normalised signature-by-sample matrix
# also uses the helper function - normalizeMatrix

quantifySignatures <- function(sample_by_component, feat_sig_mat_rds = NULL, component_by_signature=NULL)
{
  if(is.null(component_by_signature))
  {
    ##allow user specify the feat_sig_mat_rds input instead of hardcoding
    ##(my filesystem is different to yours!)
    if(!is.null(feat_sig_mat_rds)){
      #component_by_signature <- readRDS('C:/Users/conor/Google Drive/FYP/Data/Main Functions/quantifySignatures/feat_sig_mat.rds')
      component_by_signature <- readRDS(feat_sig_mat_rds)
    } else {
      print("No feat_sig_mat_rds supplied, exiting!")
      break
    }
  }
  signature_by_sample <- YAPSA::LCD(t(sample_by_component),
                                  YAPSA:::normalize_df_per_dim(component_by_signature,2))
  signature_by_sample <- normaliseMatrix(signature_by_sample)
  signature_by_sample
}


# FORMAT SIGNATURE DATA ----

# normalized signature by sample results formatted for differential expression workflow  
# maximum signature identified per sample
# signature results formatted into long format and
# a format suitable for merging with clinical data

get_DE_data <- function(norm_sig_by_samp){
  
  # display signature results in long format 
  # 7 signature posterior probability values for each sample
  
  rownames(norm_sig_by_samp) <- paste0('s', c(1:nrow(norm_sig_by_samp)))
  norm_sig_by_samp <- as.data.frame(norm_sig_by_samp)
  CN_signatures <- rownames(norm_sig_by_samp)
  rep_signatures <- rep(CN_signatures, length(colnames(norm_sig_by_samp)))
  long_DE_sig_data <- gather(norm_sig_by_samp, Sample_ID, Exposure,
                             c(colnames(norm_sig_by_samp)))
  long_DE_sig_data$CN_Signature <- as.factor(rep_signatures)
  long_DE_sig_data <- long_DE_sig_data[, c("Sample_ID", "CN_Signature", "Exposure")]
  
  # get maximum signature per sample and format to allow merging with clinical data
  DE_sig_data <- as.data.frame(t(norm_sig_by_samp))
  max_sig <- colnames(DE_sig_data)[max.col(DE_sig_data, ties.method = "first")]
  DE_sig_data$max_signature <- as.factor(max_sig)
  
  DE_data <- list(
    long_data = long_DE_sig_data,  # place into list object where both formats are available
    clinical_sigs = DE_sig_data
  )
}


