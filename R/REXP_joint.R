
REXP_joint <- function(  data, marker, trait, env, model, RET_a, RET_b ) {
  rexport_result <- .C( 'eREXP_joint',  as.double(data), as.integer(dim(data)), as.double(marker), as.integer(length(marker)), as.double(trait), as.double(env), as.double(model), as.double(RET_a), as.integer(length(RET_a)), as.double(RET_b), as.integer(dim(RET_b)), dup=FALSE, NAOK=TRUE )
  list(   a = rexport_result[[8]] ,   b = rexport_result[[10]] )
}

