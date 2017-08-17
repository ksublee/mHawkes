# m1 <- m2 <- m3 <- m4 <- matrix(1:4, nrow=2, dimnames=list(c("a", "b"), c("d", "e")))
# fn <- function(x) setNames(data.frame(.=paste("", rownames(x)), x, check.names=F, row.names=NULL),c(" ", colnames(x)))
# matrix.names <- Filter( function(x) 'matrix' %in% class( get(x) ), ls(pattern = "m") )
# matrix.list <- lapply(matrix.names, get)
# matrix.chain <- do.call(cbind, lapply(matrix.list, fn))
# cat(" ", paste0(matrix.names, collapse = "     "), "\n");
# print(matrix.chain, row.names = FALSE)


