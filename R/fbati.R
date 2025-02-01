## Get the actual coding of the markers
xcode <- function( m0, m1, model ) {
  if( model==ADDITIVE )
    return( as.integer(m0==2) + as.integer(m1==2) )
  if( model==DOMINANT )
    return( as.integer(m0==2 | m1==2) )
  if( model==RECESSIVE )
    return( as.integer(m0==2 & m1==2) )
  stop( paste( "xcode: model (", model, ") not understood.", sep="" ) )
}

## Get the 3 codings of the markers
xcodes <- function(model) {
  Xc <- c();
  if( model==ADDITIVE ) {
    Xc <- c(0,1,2);
  }else if( model==DOMINANT ){
    Xc <- c(0,1,1);
  }else if( model==RECESSIVE ){
    Xc <- c(0,0,1);
  }else{
    stop( paste("xcodes: model (",model,") is not understood.", sep="") );
  }

  return( Xc );
}

#############################
## calculating group means ##
#############################
groupMean <- function( group, x ) {
  ugroup <- unique(group);
  xbar <- rep( 0, length(group) );
  for( u in ugroup ){
    wh <- group==u;
    umean <- mean( x[wh] );
    xbar[wh] <- umean;
  }

  return(xbar);
}

## Only uses the first affected (traces to dataComputeGroupG C function),
##  so this is true for the LL method, DR could be different
fbati.calc <- function( data, ## what was passed to datamatrix
                        m0pos=7, m1pos=m0pos+1,
                        groupsG,  ## part of the result of datamatrix
                        affectedIndex, ## the other part of the result
                        envCol=m1pos+1,
                        model, iter=1000,
                        debug=FALSE ) {
  ## actually - we really first need to get rid of the data that isn't needed,
  ##  _before_ the groups are computed! (unaffected are also marked by group...)
  ## keep (1) children that (2) are affected

  ## The fix for when we have missing markers in the affected child...
  #keep <- groupsG$groups!=0 & data$AffectionStatus==2
  keep <- groupsG$groups!=0 & data$AffectionStatus==2 & data[,m0pos]!=0 & data[,m1pos]!=0
  keep <- intersect( which(keep), affectedIndex ) ## Modification so only uses first affected
  data <- data[keep,]
  groupsG <- groupsG[keep,]

  if( debug ) {
    cat( "*** data[keep,] ***\n" )
    print( data )
    cat( "*** groupsG ***\n" )
    print( groupsG )
  }

  ## get X(g), then the means
  x <- xcode( data[,m0pos], data[,m1pos], model )
  xbar <- groupMean( groupsG$groups, x )
  zbar <- groupMean( groupsG$groups, data[,envCol] )

  xmxbar <- x-xbar
  zmzbar <- data[,envCol]-zbar

  return( fbati2( xmxbar, zmzbar, groupsG$groups, iter, debug ) )
}

fbati2 <- function( xmxbar, zmzbar, group, iter=1000, debug=FALSE ){
  if( length(xmxbar)<=1 )  ## changed to LE instead of EQ
    return( list( pvalue=1, numInf=0, strataSum=NA ) )


  o <- order( group );  ## NEEDS TO BE SORTED FIRST

  if( debug ) {
    print( data.frame( xmxbar=xmxbar[o], zmzbar=zmzbar[o], group=group[o] ) )
  }

  #print( data.frame( xmxbar=xmxbar, zmzbar=zmzbar, group=group ) )

  ## then call the C function for speed.
  # pvalue <- as.double(0.0);
  # res = .C( "fbati_cpp",
  #     pvalue,
  #
  #     as.integer(length(group)),
  #     as.double(xmxbar[o]),
  #     as.double(zmzbar[o]),
  #     as.integer(group[o]),
  #
  #     as.integer(iter),
  #
  #     DUP=TRUE); #DUP=FALSE );
  # pvalue = res[[1]] ## 03/21/2014 <-- Make sure that this works!

  ## BEGIN: former c++ code, re-implemented as R code
  xmxbar <- xmxbar[o]
  zmzbar <- zmzbar[o]
  group <- group[o]

  # - default test statistic
  stat <- sum(xmxbar * zmzbar)

  # - cache indexes into the groups
  ugroup <- unique(group)
  G <- length(ugroup)
  groupl <- list()
  for(g in 1:G)
    groupl[[g]] <- which(group == ugroup[g])

  statPerm <- c()
  for(i in 1:iter){
    zmzbar_perm <- zmzbar
    for(g in 1:G)
      zmzbar_perm[groupl[[g]]] <- sample(zmzbar[groupl[[g]]])
    statPerm[i] <- sum(xmxbar * zmzbar_perm)
  }

  pvalue <- mean(abs(statPerm) >= abs(stat))

  ## END: former c++ code

  #print( pvalue )

  numInf <- sum( xmxbar!=0 & zmzbar!=0 )
  strataSum <- strataSum( xmxbar, zmzbar, group )

  #return( pvalue );
  return( list( pvalue=pvalue, numInf=numInf, strataSum=strataSum ) )
}

strataSum <- function( xmxbar, zmzbar, group ) {
  o <- order( group )
  xmxbar <- xmxbar[o]
  zmzbar <- zmzbar[o]
  group <- group[o]

  ssum <- NULL
  for( g in unique(group) ) {
    wh <- group==g
    if( any(xmxbar[wh]!=0 & zmzbar[wh]!=0) ) {
      numInf=sum(xmxbar[wh]!=0 & zmzbar[wh]!=0)
      names(numInf) <- g; #paste( "G", g, sep="" )

      if( is.null(ssum) ) {
        ssum <- numInf
      }else{
        ssum <- c( ssum, numInf )
      }
    }
  }

  res <- t(data.frame(ssum))
  row.names(res) <- NULL
  return( res )
}

strataNameDehash <- function( number ) {
  return(pG_group_dehash_r(number))

  # str <- as.character("                                                  ")  ## 50 spaces...
  # res <- .C( "pG_group_dehash", as.integer(number), str )
  # ##print( res )
  # return( res[[2]] )
}













if(FALSE){
  source("fbati.R")  # trigger everything with this one...
  source("ibat.R")
  source("merge.R")
  source("cpp_nuclify.R")
  source("cpp_fbatdist.R")
  source("cpp_datamatrix.R")
  source("cpp_cgFbat.R")
  source("joint.R")
  source("cpp_joint.R")

  ##########################
  ## TAKING FROM fbati.Rd ##
  ##########################

  ## Set the seed so you get the same results
  set.seed(13)

  ## Constants (you can vary these)
  NUM_FAMILIES <- 5 #500
  NUM_FAMILIES <- 500
  AFREQ <- 0.1  ## Allele frequency
  BG <- -0.25   ## main effect of gene
  BE <- 0       ## main effect of environment
  BGE <- 0.75   ## main gene-environment effect
  ENV <- 0.2    ## environmental exposure frequency

  ## (but don't modify this one)
  MAX_PROB <- exp( BG*2 + BE*1 + BGE*2*1 )

  #####################################
  ## Create a random dataset (trios) ##
  #####################################

  ## -- genotypes are set to missing for now,
  ##    everyone will be affected
  ## -- these are usually wrapped in as.ped and as.phe
  ped <-        data.frame( pid = kronecker(1:NUM_FAMILIES,c(1,1,1)),
                            id = kronecker( rep(1,NUM_FAMILIES), c(1,2,3) ),
                            idfath = kronecker( rep(1,NUM_FAMILIES), c(0,0,1) ),
                            idmoth = kronecker( rep(1,NUM_FAMILIES), c(0,0,2) ),
                            sex = rep(0,NUM_FAMILIES*3),
                            AffectionStatus = kronecker( rep(1,NUM_FAMILIES), c(0,0,2) ),
                            m0.a = rep(0,NUM_FAMILIES*3),      ## missing for now
                            m0.b = rep(0,NUM_FAMILIES*3) )   ## missing for now
  ## -- envioronment not generated yet
  phe <-        data.frame( pid = ped$pid,
                            id = ped$id,
                            env = rep(NA,NUM_FAMILIES*3) )   ## missing for now

  ## 50/50 chance of each parents alleles
  mendelTransmission <- function( a, b ) {
    r <- rbinom( length(a), 1, 0.75 )
    return( a*r + b*(1-r) )
  }

  ## Not the most efficient code, but it gets it done;
  ##  takes < 5 sec on pentium M 1.8Ghz
  CUR_FAMILY <- 1
  while( CUR_FAMILY<=NUM_FAMILIES ) {
    ## Indexing: start=father, (start+1)=mother, (start+2)=child
    start <- CUR_FAMILY*3-2

    ## Draw the parental genotypes from the population
    ped$m0.a[start:(start+1)] <- rbinom( 1, 1, AFREQ ) + 1
    ped$m0.b[start:(start+1)] <- rbinom( 1, 1, AFREQ ) + 1

    ## Draw the children's genotype from the parents
    ped$m0.a[start+2] <- mendelTransmission( ped$m0.a[start], ped$m0.b[start] )
    ped$m0.b[start+2] <- mendelTransmission( ped$m0.a[start+1], ped$m0.b[start+1] )

    ## Generate the environment
    phe$env[start+2] <- rbinom( 1, 1, ENV )

    ## and the affection status
    Xg <- as.integer(ped$m0.a[start+2]==2) + as.integer(ped$m0.b[start+2]==2)
    if( rbinom( 1, 1, exp( BG*Xg + BE*phe$env[start+2] + BGE*Xg*phe$env[start+2] ) / MAX_PROB ) == 1 )
      CUR_FAMILY <- CUR_FAMILY + 1
    ## otherwise it wasn't a valid drawn individual
  }


  ##############
  ## Analysis ##
  ##############

  ## Print the first 4 families
  print( head( ped, n=15 ) )
  print( head( phe, n=15 ) )

  ## NOTE: We could have just put all of this info into a single dataframe otherwise,
  ##  that would look like just the results of this
  data <- mergePhePed( ped, phe )
  fbati::nuclifyMerged(data)  ## compare to old version - IDENTICAL
  data <- nuclifyMerged(data)
  data

  fbati::strataReduce(data=data, envCol=9, m0=7, m1=8, maxSib=3)
  strataReduce_r(din=data, envCol=9, m0=7, m1=8, maxSib=2)  ## error with adding in the offspring on strataReduce_r
  strataReduce( data=data, envCol=9, m0=7, m1=8, maxSib=3 ) # calls strataReduce_r, almost directly


  fbati:::ibat2(data=data, marker=7:8, markerNames="m0", envCols=9, envColsNames="env", method="fbati", model=ADDITIVE, iter=10000, seed=13, maxSib=2, debug=!TRUE)


  ibat2(data=data, marker=7:8, markerNames="m0", envCols=9, envColsNames="env", method="fbati", model=ADDITIVE, iter=100, seed=13, maxSib=2, debug=!TRUE)

  print( head( data, n=12 ) )

  nuclify_r(ped)

  ## And run the analysis on all the markers
  fbati( ped=ped, phe=phe, env="env" )

  ## Or do it via the merged data.frame object
  ##  7 and 8 correspond to the marker columns
  fbati( data=data, env="env", marker=c(7,8) )

  # this is broken now? fix...
  fbatj( ped=ped, phe=phe, env="env", verbose=TRUE)
  fbatj( data=data, env="env", marker=c(7,8), verbose=TRUE)

  fbatj( data=data, env="env", marker=c(7,8), verbose=TRUE)


  # for hopping inside the fbatj routine
  env <- "env"
  marker <- c(7,8)
  verbose <- TRUE
  trait <- "AffectionStatus"
  model <- "additive"
  mode <- "univariate"
  fixNames <- TRUE


  ###############################
  ## END: TAKING FROM fbati.Rd ##
  ###############################
}