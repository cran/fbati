computeGroupG_r <- function(din, c_m0, c_m1) {
  start <- -1
  end <- -1

  #print(din)   ## DEBUG

  prevInformativePid <- -1  # For when firstAffectedOnly=true
  ddata_num_families <- 0

  g0 <- c()
  g1 <- c()
  g2 <- c()
  groups <- c()

  affected_index <- c()
  affected_index_size <- 0

  start_end <- c(start, end)
  start_end <- getNextFamily(din, start_end)

  while (length(start_end) > 0) { 
    start <- start_end[1]
    end <- start_end[2]

    ddata_num_families <- ddata_num_families + 1

    # Get alleles for this family
    numChild <- 0
    ca <- c()
    cb <- c()
    childi <- c()
    curParent <- 1
    p <- rep(0, 4)

    for (i in start:end) {
      if (din[i, C_FATH] == 0 && din[i, C_MOTH] == 0) {
        # Parent
        if (curParent > 2) {
          print(din[start:end, ])
          print(c_m0)
          print(c_m1)
          print(i)
          stop("More than two parents in a pedigree. Current code can only handle nuclear pedigrees where the parents have no parents.")
        }
        p[curParent * 2 - 1] <- din[i, c_m0]
        p[curParent * 2] <- din[i, c_m1]
        curParent <- curParent + 1
      } else {
        # Child
        numChild <- numChild + 1
        ca[numChild] <- din[i, c_m0]
        cb[numChild] <- din[i, c_m1]
        childi[numChild] <- i
      }
    }

    if (numChild == 0) {
      warning("No children in pedigree.")
      next
    }

    group_g <- pG_group_r(numChild, p[1:2], p[3:4], ca, cb)
    group <- group_g$hash
    g <- group_g$pg

    if (group != -1) {
      # Informative family
      for (c in 1:numChild) {
        groups[childi[c]] <- group
        g0[childi[c]] <- g[gBB+1]  # Note: g0, g1, g2 are mapped backwards
        g1[childi[c]] <- g[gAB+1]
        g2[childi[c]] <- g[gAA+1]  # i don't think we should add to this index... <- no, it should

        if (din[childi[c], C_AFF] == 2) {
          # Only use the first affected in a family
          if (prevInformativePid != din[childi[c], C_PID]) {
            affected_index_size <- affected_index_size + 1
            affected_index[affected_index_size] <- childi[c]
            prevInformativePid <- din[childi[c], C_PID]
          }
        }
      }
    }

    start_end = getNextFamily(din, start_end)
  }

  # cat("groups\n")
  # print(groups)
  # cat("g0\n")
  # print(g0)
  # cat("g1\n")
  # print(g1)
  # cat("g2\n")
  # print(g2)
  # cat("affected_index\n")
  # print(affected_index)

  return( list(
    groupsG=data.frame( groups=groups, g0=g0, g1=g1, g2=g2 ),
    affectedIndex=affected_index ) )   # i don't think i need to add to the affected index anymore...

}