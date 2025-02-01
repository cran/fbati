REXP_joint_r <- function(din, marker, trait, env, model) {
  # Extract dimensions and data from ddata
  R <- nrow(din)
  C <- ncol(din)
 
  # Initialize output vectors and matrices
  NN <- length(marker)
  a <- rep( as.double(0), NN )
  b <- matrix( as.double(0), nrow=NN, ncol=NN )
  numInf <- rep( as.double(0), 1 )   ## 

  RET_a <- a
  RET_b <- b
  RET_numInf <- numInf


  # Find family indices
  start_end <- c(-1, -1)
  start_end <- getNextFamily(din, start_end)
  while(length(start_end) > 0) {  
    cat("start_end\n")
    print(start_end)

    # Get family indices
    start <- start_end[1]
    end <- start_end[2]

    # Loop through markers
    for (m in seq(1, length(marker), 2)) {
      m0 <- marker[m]
      m1 <- marker[m + 1]

      # first find all the children and parental alleles
      p1 <- p2 <- c(0, 0)
      curP <- 0
      ca <- cb <- c()
      cEnv <- c()
      curC <- 0
      y <- c()

      for(i in start:end){
        if(din[i, C_FATH]==0 & din[i, C_MOTH]==0) {
          # it's a parent
          if(curP > 1){
            stop(paste0("Too many parents in family ", din[i, C_PID], "\n"))
          }else{
            if(curP == 0){
              p1 <- c(din[i, m0], din[i, m1])
            }else if(curP == 1){
              p2 <- c(din[i, m0], din[i, m1])
            }
            curP <- curP + 1
          }
        }else{
          # assume a child then
          curC <- curC + 1
          ca[curC] <- din[i, m0]
          cb[curC] <- din[i, m0]
          cEnv[curC] <- din[i, env]
          y[curC] <- din[i, trait]
        }
      } # i

      # LEFT OFF HERE WITH RE-WRITING THIS ONE... HOW MUCH MORE OF THIS DO WE HAVE TO DO???

      # Calculate Si0, Si1, Vi00, Vi01, Vi10, Vi11
      # (Replace fbat_Si_joint_G_GE with your actual R implementation)
      # Assuming fbat_Si_joint_G_GE is an R function with appropriate arguments
      ####fbat_Si_joint_G_GE_r <- function(n, p1, p2, ca, cb, y, z, model, offset, nPhenotyped) {
      cat("just before fbat_si_joint_g_ge_r\n")
      print(p1); print(p2); print(ca); print(cb); print(y); print(cEnv); print(model);
      print(start_end)
      results <- fbat_Si_joint_G_GE_r(curC, p1, p2, ca, cb, y, cEnv, model, 0.0, curC)
      cat("just after fbat_si_joint_g_ge_r\n")
      print(results)

      # Si0 <- results[1]
      # Si1 <- results[2]
      # Vi00 <- results[3]
      # Vi01 <- results[4]
      # Vi10 <- results[5]
      # Vi11 <- results[6]
      Si0 <- results['Si0']
      Si1 <- results['Si1']
      Vi00 <- results['Vi00']
      Vi01 <- results['Vi01']
      Vi10 <- results['Vi10']
      Vi11 <- results['Vi11']

      # Update RET_a and RET_b
      m1idx <- floor(m/2) + 1
      m2idx <- floor(m/2) + 2
      RET_a[m1idx] <- RET_a[m1idx] + Si0
      RET_a[m2idx] <- RET_a[m2idx] + Si1
      RET_b[m1idx, 1] <- RET_b[m1idx, 1] + Vi00
      RET_b[m1idx, 2] <- RET_b[m1idx, 2] + Vi01
      RET_b[m2idx, 1] <- RET_b[m2idx, 1] + Vi10
      RET_b[m2idx, 2] <- RET_b[m2idx, 2] + Vi11

      # Update numInf (replace with your actual condition)
      RET_numInf <- RET_numInf + (abs(Si0) > 0.00001 || abs(Si1) > 0.00001)
    }

    start_end <- getNextFamily(din, start_end)
  }

  print("about to return")

  return(list(a = RET_a, b = RET_b, numInf = RET_numInf))
}

### TODO - SEE WHERE THIS WAS CALLED, AND INSERT THIS IN