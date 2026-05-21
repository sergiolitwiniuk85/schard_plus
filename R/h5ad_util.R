loadRequiredPackages <- function(pkgs) {
  for (pkg in pkgs) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      stop(
        "Package '", pkg, "' is required. Please install it:\n",
        "  install.packages('", pkg, "')\n",
        "  # or for Bioconductor packages:\n",
        "  BiocManager::install('", pkg, "')"
      )
    }
  }
  invisible(TRUE)
}

#' Extracts data.frame from h5ad file
#'
#' @param filename file name or H5IdComponent to read data.frame from
#' @param name name of group that contains data.frame to be loaded
#' @param keep.rownames.as.column whether to keep rownames (index) as column in output. TRUE by default
#'
#' @return a data.frame
#' @export
#' @examples
#' obs = h5ad2data.frame('adata.h5ad','obs')
h5ad2data.frame = function(filename,name,keep.rownames.as.column=TRUE){
  collist = rhdf5::h5read(filename,name,read.attributes = TRUE)
  attr = attributes(collist)

  # HDF5 group flattening: slashes in column names create nested HDF5 groups
  # that rhdf5 returns as sub-lists. Flatten them into top-level entries
  # with "/" in the name so downstream code can handle them uniformly.
  repeat {
    ll <- vapply(collist, length, integer(1))
    max_ll <- max(ll)
    candidates <- which(ll != max_ll & names(collist) != "__categories")

    # exclude factors (categories/codes) and masked arrays (mask/values)
    irregular <- c()
    for (i in candidates) {
      val <- collist[[i]]
      if (!is.list(val) || length(val) != 2 ||
          (!all(names(val) %in% c("categories", "codes")) &&
           !all(names(val) %in% c("mask", "values")))) {
        irregular <- c(irregular, i)
      }
    }

    if (length(irregular) == 0) break

    idx <- irregular[1]
    children <- collist[[idx]]
    for (child_name in names(children)) {
      flat_name <- paste0(names(collist)[idx], "/", child_name)
      collist[[flat_name]] <- children[[child_name]]
    }
    collist[[idx]] <- NULL
  }

  # parse regular factors, vectors and masked arrays
  collist = lapply(collist,parseVector)

  # some other ways to store factors:
  # first way to store factors (all levels are in collist[['__categories']])
  for(fn in names(collist[['__categories']])){
    codes = collist[[fn]]+1
    codes[codes==0] = NA
    collist[[fn]] = as.vector(collist[['__categories']][[fn]][codes])
  }

  collist[['__categories']] = NULL

  # another way to store factors with names stored in /uns
  uns_names = rhdf5::h5ls(filename)
  uns_names = uns_names[uns_names$group=='/uns' & endsWith(uns_names$name,'_categories'),]
  uns_names$var_name = sub('_categories$','',uns_names$name)
  uns_names = uns_names[uns_names$var_name %in% names(collist),]
  for(i in seq_len(nrow(uns_names))){
    fn = uns_names$var_name[i]
    cats = rhdf5::h5read(filename,paste0('/uns/',uns_names$name[i]))
    collist[[fn]] = cats[collist[[fn]]+1]
  }

  ll = sapply(collist,length)
  if (any(ll != max(ll)))
    warning("h5ad2data.frame: unexpected data.frame format for group '", name, "' in file '", filename, "', some columns may be missing")
  res = as.data.frame(collist[ll==max(ll)],check.names=FALSE)

  index.col = NULL
  if('index' %in% colnames(res))
    index.col = 'index'
  if('_index' %in% names(attr) && attr$`_index` %in% colnames(res))
    index.col = attr$`_index`

  if(!is.null(index.col))
    rownames(res) = make.unique(res[,index.col])

  ord =NULL
  if(!is.null(attr$`column-order`) & all(attr$`column-order` %in% colnames(res))){
    ord = attr$`column-order`
    if(keep.rownames.as.column && !is.null(index.col))
      ord = c(index.col,ord)
    res = res[,ord,drop=FALSE]
  }
  res
}

# pasres masked arrays, factors and ordinary vectors into vectors
parseVector = function(x){
  if(is(x, "array")){
    return(as.vector(x))
  }
  if(is(x, "list")){
    # factors
    if(all(names(x) %in% c("categories","codes"))){
      codes = as.vector(x$codes+1)
      codes[codes==0] = NA
      return(parseVector(x$categories)[codes])
    }
    # masked arrays
    if(all(names(x) %in% c("mask","values"))){
      values = x$values
      values[x$mask=='TRUE'] = NA # don't know why, by mask if factor...
      return(as.vector(values))
    }
  }
  # return as is if it is something else
  return(x)
}



#' Parse matrix from h5ad file
#'
#' @param filename file name or H5IdComponent to read data.frame from
#' @param name name of group that contains matrix to be loaded
#' @param use_spam logical, whether to use spam instead of Matrix. Can be used if matrix has more than 2^31 -1 non-zero elements. Keep in mind that spam is not compatible with Seurat
#'
#' @return R matrix, dense or sparse - in dependence on input
#' @export
#' @examples
#' obs = h5ad2data.frame('adata.h5ad','X')
h5ad2Matrix = function(filename,name,use_spam = FALSE){
  if(!startsWith(name,'/'))
    name = paste0('/',name)
  attr = rhdf5::h5readAttributes(filename,name)
  # load as dataframe if it appers to be a dataframe, but then convert to matrix as the matrix was requested
  if(!is.null(attr$`encoding-type`) && attr$`encoding-type` =='dataframe'){
    mtx = h5ad2data.frame(filename,name,keep.rownames.as.column=FALSE)
    mtx = as.matrix(mtx)
    return(mtx)
  }

  ls = rhdf5::h5ls(filename)
  nvalues = as.numeric(ls[ls$group==name & ls$name=='data','dim'])
  if(use_spam){
    require(spam)
    require(spam64)
    m = rhdf5::h5read(filename, name,bit64conversion='double')

    mtx=new('spam',entries=as.numeric(m$data),
             colindices=as.numeric(m$indices) + 1,
             rowpointers=as.numeric(m$indptr) + 1,
             dimension=as.numeric(attr$shape))
    mtx = t(mtx)
    return(mtx)
  }
  # if there is data subgroup then it should be sparse. We can only load sparse if it has less than 2^32-1 values
  if(!is.null(nvalues) && length(nvalues)==1 && nvalues >= 2^31){
    stop("h5ad2Matrix: The object '", name, "' in file '", filename, "' is too large for Seurat and R in general: it has more than (2^31 -1) non-zero values. Consider setting use_spam=TRUE or use python.\n For more information please check:\n 1. https://github.com/cellgeni/schard/issues/1\n 2. https://github.com/chanzuckerberg/cellxgene-census/issues/1095")
  }
  m = rhdf5::h5read(filename,name)
  format = attr$`encoding-type`
  if(is.null(format))
    format = attr$h5sparse_format
  shape = attr$shape
  if(is.null(shape))
    shape = attr$h5sparse_shape
  if(is.array(m)){
    mtx = m
  }else{
    if(startsWith(format,'csr')){
      mtx = Matrix::sparseMatrix(j=m$indices+1, p=m$indptr,x = as.numeric(m$data), repr="R")
      if(!inherits(mtx, "dgCMatrix")) mtx = as(mtx, "CsparseMatrix")
      if(identical(as.integer(dim(mtx)), as.integer(shape))) mtx = Matrix::t(mtx)
    }else{
      mtx = Matrix::sparseMatrix(i=m$indices+1, p=m$indptr,x = as.numeric(m$data),dims = shape)
      mtx = Matrix::t(mtx)
    }
  }
  mtx
}



#' Load Visium images from h5ad
#'
#' @param filename path to h5ad file with Visium spatial data
#'
#' @return list of slides. Each slide has scale.factors, hires, and lowres images
#'
#' @keywords internal
#' @examples
#' imgs = h5ad2images('adata.h5ad')
h5ad2images = function(filename){
  h5struct = rhdf5::h5ls(filename)
  library_ids = h5struct$name[h5struct$group=='/uns/spatial']
  library_ids = setdiff(library_ids,'is_single') # to accommodate for cellxgene schema
  result = list()

  for(lid in library_ids){
    result[[lid]] = list()
    result[[lid]]$scale.factors = rhdf5::h5read(filename,paste0('/uns/spatial/',lid,'/scalefactors'))
    for(res in h5struct$name[h5struct$group==paste0('/uns/spatial/',lid,'/images')]){
      img = rhdf5::h5read(filename,paste0('/uns/spatial/',lid,'/images/',res))
      result[[lid]][[res]] = aperm(img,rev(seq_len(length(dim(img)))))
      mode(result[[lid]][[res]]) = 'numeric'
      if(max(result[[lid]][[res]])>1)
        result[[lid]][[res]] = result[[lid]][[res]] / 255
    }
  }
  result
}
