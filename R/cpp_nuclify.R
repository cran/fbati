# source("cpp_nuclify.R", verbose=TRUE)

C_PID <- 0 + 1
C_ID <- 1 + 1
C_FATH <- 2 + 1
C_MOTH <- 3 + 1
C_SEX <- 4 + 1
C_AFF <- 5 + 1

DEBUG_NUCLIFY_CPP_R <- !TRUE

SEX_MALE <- 1
SEX_FEMALE <- 2


sample_one <- function(x){
  if(length(x) == 1)  # special case we needed this function for...
    return(x)
  return(sample(x, 1))
}

if(FALSE){
  sample_one(5)
  sample_one(1:5)
}

# data is a matrix, sorted on C_PID
# start_end is c(start, end)
# returns new c(start, end) following passed in one, or NULL if none
getNextFamily <- function(data, start_end) {
  R <- nrow(data)

  start <- start_end[1]
  end <- start_end[2]

  if (start == -1) {
    start <- 1
  } else {
    start <- end + 1
  }

  if (start >= R) {
    return(NULL)
  }

  pid <- data[start, C_PID]
  for (i in start:R) {
    if (data[i, C_PID] == pid) {
      end <- i
    } else {
      return(c(start, end)) # reached the end (assumes sorted)
    }
  }
  return(c(start, end))
}

if (DEBUG_NUCLIFY_CPP_R) {
  ped <- data.frame(
    pid = c(1, 1, 1, 2, 2, 2, 3, 3),
    id = 1:8,
    idfath = 0,
    idmoth = 0,
    sex = 1,
    AffectionStatus = 2,
    m0.a = 1,
    m0.b = 1,
    stringsAsFactors = FALSE
  )

  start_end <- c(-1, -1)
  start_end <- getNextFamily(ped, start_end)
  start_end
  start_end <- getNextFamily(ped, start_end)
  start_end
  start_end <- getNextFamily(ped, start_end)
  start_end
  start_end <- getNextFamily(ped, start_end)
  start_end

  # Should get
  # 1 3
  # 4 6
  # 7 8
  # NULL
}


pushDataRow <- function(din, dinI, newPid, clearParents = FALSE, fullyClearParents = FALSE, clearAffection = FALSE, envCol = -1) {
  # "copy" it in
  dout <- din[dinI, ]

  # reset the pid
  dout[, C_PID] <- newPid

  # should the parents of this one be erased?
  if (clearParents) {
    dout[, C_FATH] <- 0
    dout[, C_MOTH] <- 0
  }
  # fully clear the parents? (alleles as well)
  if (fullyClearParents) {
    dout[, (C_AFF + 1):ncol(dout)] <- 0
  }

  # and clear the affection (for the other sib)
  if (clearAffection) {
    dout[, C_AFF] <- 0
  }

  if (envCol != -1) {
    dout[, envCol] <- NA
  }

  return(dout)
}

if (DEBUG_NUCLIFY_CPP_R) {
  ped[1, ]
  # - to test if erases...
  ped[1, C_FATH] <- 6
  ped[1, C_MOTH] <- 7
  ped[1, ]
  pushDataRow(ped, 1, 666, clearParents = FALSE, fullyClearParents = FALSE, clearAffection = FALSE, envCol = -1)
  pushDataRow(ped, 1, 666, clearParents = TRUE, fullyClearParents = FALSE, clearAffection = FALSE, envCol = -1)
  pushDataRow(ped, 1, 666, clearParents = FALSE, fullyClearParents = TRUE, clearAffection = FALSE, envCol = -1)
  pushDataRow(ped, 1, 666, clearParents = FALSE, fullyClearParents = FALSE, clearAffection = TRUE, envCol = -1)
}

pushEmptyRow <- function(din, newPid, id, sex) {
  dout <- din[1, ]
  # zero it out
  dout[, 1:ncol(dout)] <- 0

  # then setup the other stuff
  dout[, C_PID] <- newPid
  dout[, C_ID] <- id
  # C_FATH - already zeroed out
  # C_MOTH - already zeroed out
  dout[, C_SEX] <- sex
  # C_AFFECTION - already zeroed out

  return(dout)
}

if (DEBUG_NUCLIFY_CPP_R) {
  pushEmptyRow(ped, 777, 13, 1)
}

# replacement for nuclify_cpp
nuclify_r <- function(din) {
  dout <- list()

  # start by going across all families
  start_end <- c(-1, -1)
  start_end <- getNextFamily(din, start_end)
  while (length(start_end) > 0) {
    start <- start_end[1]
    end <- start_end[2]

    newPid <- din[start, C_PID] * 100

    # pull out each unique mother-father pair
    for (i in start:end) {

      idfath <- din[i, C_FATH]
      idmoth <- din[i, C_MOTH]

      idfathRow <- -1
      idmothRow <- -1

      numSibs <- 0
      sibs <- c()
      validFamily <- TRUE
      for (ii in start:end) {
        if (validFamily) {
          if (idfath == din[ii, C_FATH] && idmoth == din[ii, C_MOTH]) {
            # we know it's a new/unique sub-family if it hasn't been found before...
            if (ii < i) { # not the first member of the family
              validFamily <- FALSE
            } else {
              numSibs <- numSibs + 1
              sibs[numSibs] <- ii
            }
          } else if (din[ii, C_ID] == idfath) {
            idfathRow <- ii
          } else if (din[ii, C_ID] == idmoth) {
            idmothRow <- ii
          }
        } # if validFamily
      } # for ii

      if (idfath == 0 && idmoth == 0) {
        validFamily <- FALSE
      }

      if (validFamily) {
        # push it onto the output

        # the father, inserting one if not in the pedigree
        if (idfathRow != -1) {
          dout[[length(dout) + 1]] <- pushDataRow(din = din, dinI = idfathRow, newPid = newPid, clearParents = TRUE)
        } else {
          dout[[length(dout) + 1]] <- pushEmptyRow(din = din, newPid = newPid, id = idfath, sex = SEX_MALE)
        } ## if idfathRow

        # the mother, inserting one if not in the pedigree
        if (idmothRow != -1) {
          dout[[length(dout) + 1]] <- pushDataRow(din = din, dinI = idmothRow, newPid = newPid, clearParents = TRUE)
        } else {
          dout[[length(dout) + 1]] <- pushEmptyRow(din = din, newPid = newPid, id = idmoth, sex = SEX_FEMALE)
        } ## if idmothRow

        # and all of their children
        for (ch in 1:numSibs) {
          dout[[length(dout) + 1]] <- pushDataRow(din = din, dinI = sibs[ch], newPid = newPid)
        } # ch

        # increment the newPid
        newPid <- newPid + 1
        #print(paste("newPid", newPid))
      } # if validFamily
    } # ii

    # setup for next round
    start_end <- getNextFamily(din, start_end)
  }

  return(do.call("rbind", dout))
}

if (DEBUG_NUCLIFY_CPP_R) {
  # Taken from the example help file
  ped <- data.frame(
    pid = rep(1, 8),
    id = c(100, 101, 201, 202, 301, 302, 303, 304),
    idfath = c(0, 0, 100, 0, 201, 201, 201, 201),
    idmoth = c(0, 0, 101, 0, 202, 202, 202, 202),
    sex = c(1, 2, 1, 2, 2, 2, 2, 2),
    AffectionStatus = rep(0, 8),
    m0.a = rep(2, 8),
    m0.b = rep(2, 8)
  )

  din <- ped ## for stepping through...

  ped
  getNextFamily(ped, c(-1, -1)) # yup, this works as expected
  nuclify_r(ped)

  start_end <- c(-1, -1)
  start_end <- getNextFamily(ped, start_end)
  start_end
}

strataReduceRemove <- function(array, elt) {
  return(unique(array)) # I think that is all this code does

  # # find elt, and replace with last elt
  # arraySize <- length(array)
  # for (i in 1:arraySize) {
  #   if (i <= arraySize) {
  #     if (array[i] == elt) {
  #       array[i] <- array[arraySize - 1]
  #       arraySize <- arraySize - 1
  #     }
  #   }
  # }
  # return(array[1:arraySize])
}

strataReduce_r <- function(din, envCol, m0, m1, maxSib) {
  dout <- list()

  # do the processing
  start_end <- c(-1, 1)
  start_end <- getNextFamily(din, start_end)
  while (length(start_end) > 0) {
    start <- start_end[1]
    end <- start_end[2]

    pid <- din[start, C_PID]

    # find the mother/father ids
    idmoth <- -1
    idfath <- -1
    for (i in start:end) {
      if (din[i, C_MOTH] != 0) {
        idmoth <- din[i, C_MOTH]
      }
      if (din[i, C_FATH] != 0) {
        idfath <- din[i, C_FATH]
      }
    } # i

    # find all children and parents locations for current family
    # - basic children/parents location
    parentRow <- c()
    parentRowSize <- 0
    # - affected, environmental info, genotyped
    childRowAffected <- c()
    childRowAffectedSize <- 0
    # - genotyped (for sibpair)
    childRowGeno <- c()
    childRowGenoSize <- 0
    # -- so the childRowAffected set is contained in childRowGeno set
    for (i in start:end) {
      if (din[i, C_ID] == idfath || din[i, C_ID] == idmoth) {
        if (parentRowSize == 2) {
        } else {
          parentRowSize <- parentRowSize + 1
          parentRow[parentRowSize] <- i
        }
      } else {
        # childRowSize <- childRowSize + 1

        if (din[i, C_AFF] == 2 && ## affected
          !is.nan(din[i,envCol]) && ## environment
          din[i, m0] != 0 & din[i, m1] != 0) { # genotyped
          childRowAffectedSize <- childRowAffectedSize + 1
          childRowAffected[childRowAffectedSize] <- i
        } # fi

        if (din[i, m0] != 0 && din[i, m1] != 0) { # genotyped
          childRowGenoSize <- childRowGenoSize + 1
          childRowGeno[childRowGenoSize] <- i
        }
      } # else
    } # for i start:end

    if (parentRowSize == 2 &&
      din[parentRow[1], m0] != 0 && din[parentRow[1], m1] != 0 &&
      din[parentRow[2], m0] != 0 && din[parentRow[2], m1] != 0) {
      # yes we've got both parents!

      # have at least a usable affected?
      if (childRowAffectedSize > 0) {
        # push the parents on
        for (p in 1:2) {
          dout[[length(dout) + 1]] <- pushDataRow(din = din, dinI = parentRow[p], newPid = pid, clearParents = TRUE)
        }

        # choose a random affected
        dout[[length(dout) + 1]] <- pushDataRow(din = din, dinI = sample_one(childRowAffected), newPid = pid)
      }
    } else {
      # No, we don't have both parents

      # have at least a usable affected, and another genotyped?
      if (childRowAffectedSize > 0 && childRowGenoSize > 1) {
        affectedRow <- sample_one(childRowAffected)
        dout[[length(dout) + 1]] <- pushDataRow(din = din, dinI = affectedRow, newPid = pid)
        ##
        childRowGeno <- strataReduceRemove(childRowGeno, affectedRow)
        childRowGenoSize <- length(childRowGeno)

        if ((childRowGenoSize + 1) < maxSib) {
          # push all offspring on
          dout[[length(dout) + 1]] <- pushDataRow(din = din, dinI = childRowGeno[c], newPid = pid, clearAffection = TRUE, envCol = envCol)
        } else {
          # random # of offspring to push on
          for (c in 2:maxSib) { # +1 for affected
            unaffectedRow <- sample_one(childRowGeno)
            dout[[length(dout) + 1]] <- pushDataRow(din = din, dinI = unaffectedRow, newPid = pid, clearAffection = TRUE)

            childRowGeno <- strataReduceRemove(childRowGeno, unaffectedRow)
            childRowGenoSize <- length(childRowGeno)
          }
        }
      }
    } # fi parentRowSize

    # setup for next loop
    start_end <- getNextFamily(din, start_end)
  } # while

  return(do.call("rbind", dout))
}

if(DEBUG_NUCLIFY_CPP_R){
  nped <- nuclify_r(ped)
  set.seed(13)
  nped$env <- rnorm(nrow(nped))
  nped$AffectionStatus <- 1
  nped$AffectionStatus[c(3,6)] <- 2
  nped$m0.a <- 1
  nped
  
  din <- nped
  envCol <- "env"
  m0 <- "m0.a"
  m1 <- "m0.b"
  maxSib <- 20
  
  strataReduce_r(din = nped, envCol = "env", m0 = "m0.a", m1 = "m0.b", 20)
}

