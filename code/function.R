

# arc2d8 ------------------------------------------------------------------

arc2d8 <- function(x) {
  ## whitebox uses flow pointer D8
  # 64, 128,  1
  # 32,   0,  2 
  # 16,   8,  4
  
  ## ArcGIS uses D8 algorithm
  # 32, 64, 128
  # 16,  0,   1 
  #  8,  4,   2
  
  # check class()
  if(class(x) != "SpatRaster") stop("Provide data in class 'SpatRaster'")
  
  # begin with northeast through north
  fdir_arc <- c(0, 2^(0:7), 247, 255)
  fdir_d8 <- c(0, fdir_arc[3:9], fdir_arc[2], NA, NA)
  y <- terra::subst(x, from = fdir_arc, to = fdir_d8)
  
  return(y)  
}

# from matrix to vector ---------------------------------------------------

m2v <- function(x) {
  
  if (is.null(colnames(x))) {
    index_from <- rep(1:dim(x)[1], times = dim(x)[2])
    index_to <- rep(1:dim(x)[2], each = dim(x)[1])
  } else {
    index_from <- rep(colnames(x), times = dim(x)[2])
    index_to <- rep(rownames(x), each = dim(x)[1])
  }
  
  return(dplyr::tibble(distance = c(x),
                       from = index_from,
                       to = index_to))
  
}

