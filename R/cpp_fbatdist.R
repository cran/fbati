# Define constants
gMiss <- -1
gAA <- 0
gAB <- 1
gBB <- 2

# PROBABLY NEEDS A PLUS ONE WHEN THESE ARE USED
gAAgAA <- 0
gAAgAB <- 1
gAAgBB <- 2
gABgAA <- 3
gABgAB <- 4
gABgBB <- 5
gBBgAA <- 6
gBBgAB <- 7
gBBgBB <- 8

MODEL_ADDITIVE <- 0
MODEL_DOMINANT <- 1
MODEL_RECESSIVE <- 2
MODEL_FBAT_CHARS <- c('a', 'd', 'r')

DA_MA <- 0
DA_MB <- 1
DB_MA <- 2
DB_MB <- 3

ALLELE_A <- 2
ALLELE_B <- 1


# Define helper functions
first <- function(a, b, c) {
  return(a > 0 && b == 0 && c == 0)
}

second <- function(b, a, c) {
  return(a > 0 && b == 0 && c == 0)
}

third <- function(c, b, a) {
  return(a > 0 && b == 0 && c == 0)
}

first_second <- function(a, b, c) {
  return(a > 0 && b > 0 && c == 0)
}

first_third <- function(a, c, b) {
  return(a > 0 && b > 0 && c == 0)
}

second_third <- function(c, b, a) {
  return(a > 0 && b > 0 && c == 0)
}

fbatdist_all <- function(a, b, c) {
  return(a > 0 && b > 0 && c > 0)
}

printFamily <- function(p1, p2,
                        ca, cb,
                        numSibs ) {
  cat(sprintf("P: %d %d, %d %d\nC: ", p1[1+0], p1[1+1], p2[1+0], p2[1+1]), sep="")
  for(i in 1:numSibs)
    cat(sprintf("%d %d, ", ca[i], cb[i]))
  cat("\n")
}

# NEW 2025-01-31
pG_cacb <- function(gP1, gP2, ca, cb) {
  nG <- c(0, 0, 0)
  for(i in 1:length(ca)) {
    gcodefme <- gCode(ca[i], cb[i])
    if(gcodefme != -1)
      nG[gcodefme] <- nG[gcodefme] + 1
  }

  return(pG(gCode(gP1[1], gP1[2]), gCode(gP2[1], gP2[2]), nG[1], nG[2], nG[3]))
}

# Define the main function
pG <- function(gP1, gP2, nAA, nAB, nBB) {
  # Swap gP1 and gP2 if gP1 is missing
  if (gP1 == gMiss) {
    temp <- gP1
    gP1 <- gP2
    gP2 <- temp
  }

  nSum <- nAA + nAB + nBB

  pg <- c(0, 0, 0)

  # Neither parent is missing
  if (gP1 != gMiss) {
    if (gP1 == gAA) {
      if (gP2 == gAA) {
        pg <- c(1,0,0)
        return(pg)
      } else if (gP2 == gAB) {
        pg <- c(0.5,0.5,0)
        return(pg)
      } else if (gP2 == gBB) {
        pg <- c(0,1,0)
        return(pg)
      }
    } else if (gP1 == gAB) {
      if (gP2 == gAA) {
        pg <- c(0.5,0.5,0)
        return(pg)
      } else if (gP2 == gAB) {
        pg <- c(0.25,0.5,0.25)
        return(pg)
      } else if (gP2 == gBB) {
        pg <- c(0,0.5,0.5)
        return(pg)
      }
    } else if (gP1 == gBB) {
      if (gP2 == gAA) {
        pg <- c(0,1,0)
        return(pg)
      } else if (gP2 == gAB) {
        pg <- c(0,0.5,0.5)
        return(pg)
      } else if (gP2 == gBB) {
        pg <- c(0,0,1)
        return(pg)
      }
    }
  }

  # One parent is missing
  if (gP1 == gAA || gP1 == gBB) {
    # One homozygous parent
    if (first(nAA, nAB, nBB)) {
      pg <- c(1,0,0)
      return(pg)
    } else if (second(nAA, nAB, nBB)) {
      pg <- c(0,1,0)
      return(pg)
    } else if (third(nAA, nAB, nBB)) {
      pg <- c(0,0,1)
      return(pg)
    } else if (first_second(nAA, nAB, nBB)) {
      pg <- c(0.5,0.5,0)
      return(pg)
    } else if (first_third(nAA, nAB, nBB)) {
      pg <- c(0.5,0,0.5)
      return(pg)
    } else if (second_third(nAA, nAB, nBB)) {
      pg <- c(0,0.5,0.5)
      return(pg)
    }
  }

  # Otherwise, both parents are missing, or the parent is heterozygous
  if (gP1 == gAB || gP1 == gMiss) {
    # One heterozygous parent
    if (first(nAA, nAB, nBB)) {
      pg <- c(1,0,0)
      return(pg)
    } else if (second(nAA, nAB, nBB)) {
      pg <- c(0,1,0)
      return(pg)
    } else if (third(nAA, nAB, nBB)) {
      pg <- c(0,0,1)
      return(pg)
    } else if (!first_third(nAA, nAB, nBB) && !all(nAA, nAB, nBB)) {
      pg <- c(
        nAA / nSum,
        nAB / nSum,
        nBB / nSum)
      return(pg)
    } else {
      pg <- c(
        (4^(nSum - 1) - 3^(nSum - 1)) / (4^nSum - 2 * 3^nSum + 2^nSum),
        pg[gAA],
        1.0 - pg[gAA] - pg[gBB])
      return(pg)
    }
  }

  # Handle special cases
  if ((gP1 == gAA || gP1 == gBB) && gP2 == gMiss && all(nAA, nAB, nBB)) {
    # Handle potential error or warning based on your specific needs
    # For example:
    # warning("WARNING: Impossible genotype in file.")
    return(pg)
  }

  # This really only happens when genotypes of all the children are missing
  # In which case they will be tossed anyway, so no need to warn anymore.
  return(pg)
}


## NOT FULLY CONFIDENT ON BELOW CODE, from fbatdist.cpp



# Helper function to extract a digit from the right-hand side
extractDigitRHS <- function(number, digit) {
  number <- as.numeric(number)   # ??? really shouldn't need this???
  #cat("extractDigitRHS\n")
  #print(number)
  for (i in 1:digit) {
    number <- floor(number / 10)
  }
  temp <- floor(number / 10)
  number <- number - temp * 10
  return(number)
}

if(FALSE){
  extractDigitRHS(525600, 3)
}


pG_group_dehash_gstr_r <- function(g) {
  if (g == (gAA + 1)) {
    str <- "AA"
  } else if (g == (gAB + 1)) {
    str <- "AB"
  } else if (g == (gBB + 1)) {
    str <- "BB"
  } else {
    str <- "?"
  }
  
  return(str)
}


pG_group_dehash_r <- function(number) { 
  p1 <- extractDigitRHS(number, 7)
  p2 <- extractDigitRHS(number, 6)
  
  n <- c()
  for (j in 0:2)
    n[j+1] <- extractDigitRHS(number, j * 2) + 10 * extractDigitRHS(number, j * 2 + 1)
  
  # Assuming pG_group_dehash_gstr is a defined function in R
  p1str <- pG_group_dehash_gstr_r(p1)
  p2str <- pG_group_dehash_gstr_r(p2)
  
  str = ""
  if (p1 != 0 && p2 != 0) {
    str <- sprintf("%s,%s", p1str, p2str)
  } else {
    str <- sprintf("%s,%s - AA%i AB%i BB%i", p1str, p2str, n[gAA], n[gAB], n[gBB])
  }
  
  return(str)
}

# Helper function to encode genotypes (AA, AB, BB)
gCode <- function(a, b) {
  if (a == 0 || b == 0) {
    return(gMiss)
  } else if (a == ALLELE_A && b == ALLELE_A) {
    return(gAA)
  } else if (a == ALLELE_B && b == ALLELE_B) {
    return(gBB)
  } else {
    return(gAB)
  }
}

# Function to calculate genotype probabilities
pG_group_r <- function(
  n,
  p1, p2, # parental alleles
  ca, cb) { # children alleles
  # Calculate genotypes of children
  nG <- c(0, 0, 0)
  for (i in 1:n) {
    nG[gCode(ca[i], cb[i]) + 1] <- nG[gCode(ca[i], cb[i]) + 1] + 1
  }

  # Calculate genotypes of parents
  gP1 <- gCode(p1[1], p1[2])
  gP2 <- gCode(p2[1], p2[2])

  # Calculate genotype probabilities (pG function not provided, assume it's available)
  pg <- pG(gP1, gP2, nG[1], nG[2], nG[3])

  #cat("pg\n")
  #print(pg)

  # Calculate parent code
  pCode <- ifelse(gP1 > gP2, 
                 (gP1 + 1) * 1e7 + (gP2 + 1) * 1e6, 
                 (gP2 + 1) * 1e6 + (gP1 + 1) * 1e7)

  # Return appropriate value
  hash = NA
  if (gP2 != -1 && gP1 != -1) {
    hash <- pCode
  } else {
    hash <- (nG[1] * 1e0 + nG[2] * 1e2 + nG[3] * 1e4) + pCode
  }

  return(list(hash=hash, pg=pg))
}












# Translate C++ functions xCode and ggConvert to R

# Define genotype codes (replace with your actual definitions if different)
ALLELE_A <- 1
ALLELE_B <- 0 
gAA <- 0
gAB <- 1
gBB <- 2
MODEL_ADDITIVE <- 0
MODEL_DOMINANT <- 1
MODEL_RECESSIVE <- 2

# xCode function (version 1)
xCode1 <- function(a, b, MODEL) {
  if(MODEL == MODEL_ADDITIVE)
    return((a == ALLELE_A) + (b == ALLELE_A))
  if(MODEL == MODEL_DOMINANT)
    return(as.integer(a == ALLELE_A | b == ALLELE_A))
  if(MODEL == MODEL_RECESSIVE)
    return(as.integer(a == ALLELE_A & b == ALLELE_A))
  stop(paste0("xCode (1) out of bounds! a = ", a, ", b = ", b))
}
#xCode1(1, 1, 0) # doesn't work with switch - fixed

# xCode function (version 2)
xCode2 <- function(g, MODEL) {
  if(MODEL == gAA)
    return(xCode1(ALLELE_A, ALLELE_A, MODEL))
  if(MODEL == gAB)
    return(xCode1(ALLELE_A, ALLELE_B, MODEL))
  if(MODEL == gBB)
    return(xCode1(ALLELE_B, ALLELE_B, MODEL))
  stop("xCode (2) out of bounds! g = ", g)
}

# ggConvert function
ggConvert <- function(g1, g2) {
  g1 * 3 + g2
}


pGG <- function(gP1, gP2, nAA, nAB, nBB) {
  # Define genotype codes (replace with your actual definitions if different)
  gMiss <- -1
  gAA <- 0
  gAB <- 1
  gBB <- 2
  gAAgAA <- 0
  gAAgAB <- 1
  gAAgBB <- 2
  gABgAA <- 3
  gABgAB <- 4
  gABgBB <- 5
  gBBgAA <- 6
  gBBgAB <- 7
  gBBgBB <- 8

  # Swap gP1 and gP2 if gP1 is missing
  if (gP1 == gMiss) {
    temp <- gP1
    gP1 <- gP2
    gP2 <- temp
  }

  n <- nAA + nAB + nBB

  # Initialize pgg
  pgg <- rep(0, 9)

  # Case 1: Neither parent is missing
  if (gP2 != gMiss) {
    pg <- pG(gP1, gP2, nAA, nAB, nBB)  # Assuming pG function is defined
    pgg[1:9] <- c(pg[1] * pg[1], pg[1] * pg[2], pg[1] * pg[3], 
                 pg[2] * pg[1], pg[2] * pg[2], pg[2] * pg[3], 
                 pg[3] * pg[1], pg[3] * pg[2], pg[3] * pg[3])
    return(pgg)
  }

  # Case 2: Only one parent is missing
  if (gP1 != gMiss) {
    if (gP1 == gBB || gP1 == gAA) {
      if (first(nAA, nAB, nBB)) {
        pgg[gAAgAA] <- 1
      } else if (second(nAA, nAB, nBB)) {
        pgg[gABgAB] <- 1
      } else if (third(nAA, nAB, nBB)) {
        pgg[gBBgBB] <- 1
      } else if (first_second(nAA, nAB, nBB)) {
        pgg[gAAgAB] <- 2^(n - 2) / (2^n - 2)
        pgg[gAAgAA] <- (2^(n - 2) - 1) / (2^n - 2)
        pgg[gABgAB] <- pgg[gAAgAA]
        pgg[gABgAA] <- pgg[gAAgAB]
      } else if (first_third(nAA, nAB, nBB)) {
        stop("Impossible genotypes, 1 missing parent.")
      } else if (second_third(nAA, nAB, nBB)) {
        pgg[gBBgAB] <- 2^(n - 2) / (2^n - 2)
        pgg[gBBgBB] <- (2^(n - 2) - 1) / (2^n - 2)
        pgg[gABgAB] <- pgg[gBBgBB]
        pgg[gABgBB] <- pgg[gBBgAB]
      }
    }
  }

  # Case 3: Both parents are missing or gP1 is heterozygous
  if (gP1 == gAB || gP1 == gMiss) {
    if (first(nAA, nAB, nBB)) {
      pgg[gAAgAA] <- 1
    } else if (second(nAA, nAB, nBB)) {
      pgg[gABgAB] <- 1
    } else if (third(nAA, nAB, nBB)) {
      pgg[gBBgBB] <- 1
    } else if (first_second(nAA, nAB, nBB)) {
      pgg[gAAgAA] <- (nAA * (nAA - 1)) / (n * (n - 1))
      pgg[gABgAB] <- (nAB * (nAB - 1)) / (n * (n - 1))
      pgg[gAAgAB] <- (nAA * nAB) / (n * (n - 1))
      pgg[gABgAA] <- pgg[gAAgAB]
    } else if (first_third(nAA, nAB, nBB) || all(nAA, nAB, nBB)) {
      # ... (Implement the remaining calculations for this case)
    } else if (second_third(nAA, nAB, nBB)) {
      pgg[gBBgBB] <- (nBB * (nBB - 1)) / (n * (n - 1))
      pgg[gABgAB] <- (nAB * (nAB - 1)) / (n * (n - 1))
      pgg[gBBgAB] <- (nBB * nAB) / (n * (n - 1))
      pgg[gABgBB] <- pgg[gBBgAB]
    }
  }

  return(pgg)
}


# updated in R export

fbat_Si_joint_G_GE_r <- function(n, p1, p2, ca, cb, y, z, model, offset, nPhenotyped) {
  cat("hello?\n")
  
  Si0 = 0
  Si1 = 0
  Vi00 = 0
  Vi01 = 0
  Vi10 = 0
  Vi11 = 0.0

  # Get rid of individuals with missing trait or bad genotypes
  ngenopheno <- 0
  y_genopheno <- c()
  z_genopheno <- c()
  ca_genopheno <- c()
  cb_genopheno <- c()
  ca_geno <- c()
  cb_geno <- c()
  ngeno <- 0
  for(j in 1:n){
    if(ca[j]!=0 && cb[j]!=0){
      if(!is.nan(y[j]) & !is.nan(z[j])){
        ngenopheno <- ngenopheno + 1
        y_genopheno[ngenopheno] <- y[j]
        z_genopheno[ngenopheno] <- z[j]
        ca_genopheno[ngenopheno] <- ca[j]
        cb_genopheno[ngenopheno] <- cb[j]
      }else{
        ngeno <- ngeno + 1
        ca_geno[ngeno] <- ca[j]
        cb_geno[ngeno] <- cb[j]
      }
    }
  }
  
  if (ngenopheno == 0) {
    return(c(Si0 = 0, Si1 = 0, Vi00 = 0, Vi01 = 0, Vi10 = 0, Vi11 = 0))
  }
  
  y <- c(y_genopheno, rep(-99999, ngeno))
  z <- z_genopheno
  ca <- c(ca_genopheno, ca_geno)
  cb <- c(cb_genopheno, cb_geno)

  # Calculate pg (probabilities of genotypes)
  cat("is it dieing here?\n")
  print(n)
  print(nPhenotyped)
  print("evil")
  print(p1)
  print(p2)
  print(ca)
  print(cb)
  # pG <- function(gP1, gP2, nAA, nAB, nBB) {
  #pg <- pG(n, p1, p2, ca, cb) 
  pg <- pG_cacb(p1, p2, ca, cb)
  cat("is it dieing here...?\n")
  
  # Calculate pgg (probabilities of genotype pairs)
  if (n > 1 && nPhenotyped > 1) {
    pgg <- pGG(n, p1, p2, ca, cb)
  } else {
    pgg <- NULL
  }

  cat("pgg\n")
  print(pgg)

  cat("or maybe here?\n")

  # Calculate exj (expected genotype)
  #exj <- sum(pg * xCode2(0:2, model)) 
  exj <- 0
  for(g in 1:3)
    exj <- exj + pg[g] * xCode2(g, model)

  cat("Si0[ahhhhh]=", Si0, "\n")

  for(j in 1:min(n, nPhenotyped)){
    cat("ca\n"); print(ca); cat("cb\n"); print(cb); cat("model", model, "\n")
    temp <- xCode1(ca[j], cb[j], model)
    cat("j", j, "temp", temp, "\n")
    Si0 <- Si0 + (y[j]-offset) * temp;
    Si1 <- Si1 + (y[j]-offset) * z[j] * temp;
  }

  cat("Si0=", Si0, "\n")

  # now the variance
  # we've got a special case if there is only one child... (power hack)
  if( n==1 || nPhenotyped==1 ) {
    Vi = 0.0
    # E(X^2)
    for(g in 1:3){
      x = xCode2(g, model)
      Vi <- Vi + x*x*pg[g]
    }
    # - E(X)^2
    Vi <- Vi - exj*exj;
    Vi00 = Vi * (y[1]-offset)*(y[1]-offset)
    Vi01 = Vi * (y[1]-offset)*(y[1]-offset) * z[1]
    Vi10 = Vi01
    Vi11 = Vi * (y[1]-offset)*(y[1]-offset) * z[1] * z[1]
  }else{
    sumTj0 = 0.0;
    sumTj1 = 0.0;
    for(j in 1:min(n, nPhenotyped)){
      sumTj0 <- sumTj0 + (y[j]-offset)
      sumTj1 <- sumTj1 + (y[j]-offset) * z[j]
    }
    Vi = 0.0
    for(g1 in 1:3)
      for(g2 in 1:3)
        Vi = Vi + xCode2(g1,model)*xCode2(g2,model)*( pgg[ggConvert(g1,g2)] - pg[g1]*pg[g2] )
    Vi00 = Vi * sumTj0 * sumTj0
    Vi01 = Vi * sumTj0 * sumTj1
    Vi10 = Vi * sumTj1 * sumTj0
    Vi11 = Vi * sumTj1 * sumTj1

    # the second term
    for(j in 1:min(n, nPhenotyped)){
      sum = 0.0; # jth term sum over g, g~ -- really x-exs
      for(g1 in 1:3){
        sum = sum + (xCode2(g1,model)^2)*pg[g1]
        for(g2 in 1:3){
          sum = sum - xCode2(g1,model)*xCode2(g2,model) * pgg[ggConvert(g1,g2)]
        }
      }

      temp = (y[j]-offset) * (y[j]-offset)
      sum00 = sum * temp
      sum01 = sum * temp * z[j]
      sum10 = sum * temp * z[j]
      sum11 = sum * temp * z[j] * z[j]

      Vi00 = Vi00 + sum00
      Vi01 = Vi01 + sum01
      Vi10 = Vi10 + sum10
      Vi11 = Vi11 + sum11
    }
  }

  # Return results
  return(c(Si0 = Si0, Si1 = Si1, Vi00 = Vi00, Vi01 = Vi01, Vi10 = Vi10, Vi11 = Vi11))
}
