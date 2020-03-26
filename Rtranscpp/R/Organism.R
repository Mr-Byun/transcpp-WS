require(Rcpp)
require(colorspace)

ap2lin <- function(x) { return(x+8191.5) }
lin2ap <- function(x) { return(x-8192.5) }

set_dim <- function(n)
{
  nx <- 0
  ny <- 0
  square_root <- sqrt(n)
  int_part    <- floor(square_root)
  remainder   <- square_root %% 1
  if (!remainder)
  {
    nx <- square_root
    ny <- square_root
    dimension <- n
    dif       <- 0
  } 
  else if (remainder < 0.5)
  {
    nx <- int_part
    ny <- int_part  + 1
    dimension <- nx*ny
    dif       <- dimension - n
  }
  else
  {
    nx <- int_part + 1
    ny <- int_part + 1
    dimension <- nx*ny
    dif       <- dimension - n
  }
  return(c(n,nx,ny,dif))
}

is.graph_param <- function(x)
{
  p <- par()
  if (is.null(p[[x]]))
    return(FALSE)
  else
    return(TRUE)
}

graph_params <- function(x)
{
  for (i in names(x))
  {
    if (!is.graph_param(i))
      x[[i]] <- NULL
  }
  return(x)
}

heat.ramp <- colorRampPalette(colors=c('black',"blue2", "#007FFF", "cyan2","#7FFF7F", "yellow2", "#FF7F00", "red2","#550000"))(100)
rwb <- colorRampPalette(colors=c('black','blue','red','white'))(100)

# While we have to set s4 methods in zzz.R after library load, we can do s3 methods
# with ease! Here I overload the plot command for organism.

plot.Rcpp_Organism <- function(object, type="rate", setdim=TRUE, ...)
{
  require(colorspace)
  # get all the options that were provided
  opts <- as.list(substitute(list(...)))[-1L]

  # setup par
  newpar <- list()
  
  for (i in ls(object$par_defaults))
  {
    if (is.null(opts[[i]]) & i != "initialize")
      newpar[[i]] = object$par_defaults[[i]]
  }
  oldpar <- par(newpar)
    
  # now dispatch to the appropriate plotting function
  if (type == "rate")
    plot_rate(object, setdim, opts)
  else if (type == "R2D")
    plot_R2D(object, setdim, opts)
  else if (type == "N2D")
    plot_N2D(object, setdim, opts)
  else if (type == "T2D")
    plot_T2D(object, setdim, opts)
  else if (type == "RT2D")
    plot_RT2D(object, setdim, opts)
  else if (type == "f" | type=="fa" | type=="F")
    plot_f(object, type, setdim, opts)
  else if (type=="K" | type=="score" | type=="GoF" | type=="adj")
    plot_scores(object, type, setdim, opts)
  #else if (type=="pwm")
  #  plot_pwms()
  #else if (type=="")
  else
    cat(paste("Did not recognize option of type",type,"\n"))
    
  par(oldpar)
}

plot_points <- function(object, setdim, opts)
{
  rate <- as.numeric(as.matrix(object$rate))
  data <- as.numeric(as.matrix(object$scaled_rate_data))

  if (setdim) par(mfrow=c(1,1))
  
  mrate <- max(rate)
  mdata <- max(data)
  ylim <- c(0, mrate*1.06)
  xlim <- c(0, mdata*1.06)
  ylab <- "fit"
  xlab <- "data"
  cols <- c('black')
  
  if (!is.null(opts[["xlim"]]))
    xlim = opts[["xlim"]]
  if (!is.null(opts[["ylim"]]))
    ylim = opts[["ylim"]]
  if (!is.null(opts[["xlab"]]))
    xlab = opts[["xlab"]]
  if (!is.null(opts[["ylab"]]))
    ylab = opts[["ylab"]]
  if (!is.null(opts[["col"]]))
  {
    cols = opts[['col']]
    opts[['col']] <- NULL
  }
}

plot_rate <- function(object, setdim, opts)
{
  require(colorspace)
  gnames <- object$gene_names
  if (!is.null(opts[["gene"]]))
    gnames <- eval(opts[["gene"]]) # now you can set the gene with gene="gname"
  
  ngenes <- length(gnames)
  d <- set_dim(ngenes)
  if (setdim) par(mfrow=c(d[2],d[3]))
  
  rate <- object$rate
  data <- object$scaled_rate_data
  
  is_numeric <- is.numeric(rownames(rate))
  
  m <- max(max(rate),max(data))
  ylim <- c(0, m*1.06)
  xlab <- "nuc"
  ylab <- "rate"
  cols <- c('black','red')
  if (!is.null(opts[["ylim"]]))
    ylim = opts[["ylim"]]
  if (!is.null(opts[["xlab"]]))
    xlab = opts[["xlab"]]
  if (!is.null(opts[["ylab"]]))
    ylab = opts[["ylab"]]
  if (!is.null(opts[["col"]]))
  {
    cols = opts[['col']]
    opts[['col']] <- NULL
  }
  if (length(cols) <2)
    cols <- rep(cols[1],2)

  for (i in gnames)
  {
    y1 <- data[[i]]
    y2 <- rate[[i]]
    call_list <- list()
    if (is_numeric) {
      call_list <- c(call_list, list(x=as.numeric(rownames(rate)), y=y1))
    } else {
      call_list <- c(call_list, list(x=y1))
    }
    call_list <- c(call_list, list(xlab=xlab,ylab=ylab,ylim=ylim, main=i,type='l',col=cols[1]))
    call_list <- c(call_list, graph_params(opts))
    do.call(plot, call_list)
    
    call <- list()
    if (is_numeric) {
      call <- c(call, list(x=as.numeric(rownames(rate)), y=y2))
    } else {
      call <- c(call, list(x=y2))
    }
    call <- c(call, list(col=cols[2]))
    call <- c(call, graph_params(opts))
    do.call(lines, call)
  }
}

plot_R2D <- function(object, setdim, opts)
{
  require(colorspace)
  gnames <- object$gene_names
  if (!is.null(opts[["gene"]]))
    gnames <- opts[["gene"]] # now you can set the gene with gene="gname"
    
  ngenes <- length(gnames)
  d <- set_dim(ngenes)
  if (setdim) par(mfrow=c(d[2],d[3]))
  
  xlab = "bp"
  ylab = "nuc"
  if (!is.null(opts[["xlab"]]))
    xlab = opts[["xlab"]]
  if (!is.null(opts[["ylab"]]))
    ylab = opts[["ylab"]]

  for (i in gnames)
  {
    mat <- object$R2D[[i]]
    z <- t(mat)
    x <- as.numeric(colnames(mat))
    y <- c()
    if (is.numeric(rownames(mat))) {
      y <- as.numeric(rownames(mat))
    } else {
      y <- seq(length.out=length(rownames(mat)))
    }
    call <- list(x,y,z, xlab=xlab, ylab=ylab, main=paste(i,"R"), useRaster=T)
    call <- c(call, graph_params(opts))
    do.call(image, call)
  }
}

plot_N2D <- function(object, setdim, opts)
{
  require(colorspace)
  gnames <- object$gene_names
  if (!is.null(opts[["gene"]]))
    gnames <- opts[["gene"]] # now you can set the gene with gene="gname"
    
  ngenes <- length(gnames)
  d <- set_dim(ngenes)
  if (setdim) par(mfrow=c(d[2],d[3]))
  
  xlab = "bp"
  ylab = "nuc"
  if (!is.null(opts[["xlab"]]))
    xlab = opts[["xlab"]]
  if (!is.null(opts[["ylab"]]))
    ylab = opts[["ylab"]]

  for (i in gnames)
  {
    mat <- object$N2D[[i]]
    z <- t(mat)
    x <- as.numeric(colnames(mat))
    y <- c()
    if (is.numeric(rownames(mat))) {
      y <- as.numeric(rownames(mat))
    } else {
      y <- seq(length.out=length(rownames(mat)))
    }
    call <- list(x,y,z, xlab=xlab, ylab=ylab, main=paste(i,"R"), useRaster=T)
    call <- c(call, graph_params(opts))
    do.call(image, call)
  }
}

plot_T2D <- function(object, setdim, opts)
{
  require(colorspace)
  gnames <- object$gene_names
  if (!is.null(opts[["gene"]]))
    gnames <- opts[["gene"]] # now you can set the gene with gene="gname"
    
  ngenes <- length(gnames)
  d <- set_dim(ngenes)
  if (setdim) par(mfrow=c(d[2],d[3]))
  
  xlab = "bp"
  ylab = "nuc"
  if (!is.null(opts[["xlab"]]))
    xlab = opts[["xlab"]]
  if (!is.null(opts[["ylab"]]))
    ylab = opts[["ylab"]]

  for (i in gnames)
  {
    mat <- object$T2D[[i]]
    z <- t(mat)
    x <- as.numeric(colnames(mat))
    y <- c()
    if (is.numeric(rownames(mat))) {
      y <- as.numeric(rownames(mat))
    } else {
      y <- seq(length.out=length(rownames(mat)))
    }
    call <- list(x,y,z, xlab=xlab, ylab=ylab, main=paste(i,"R"), useRaster=T)
    call <- c(call, graph_params(opts))
    do.call(image, call)
  }
}


plot_RT2D <- function(object, setdim, opts)
{
  require(colorspace)
  gnames <- object$gene_names
  if (!is.null(opts[["gene"]]))
  {
    gnames <- opts[["gene"]] # now you can set the gene with gene="gname"
    opts[["gene"]] <- NULL
  }
    
  ngenes <- length(gnames)
  d <- set_dim(ngenes)
  if (setdim) par(mfrow=c(d[2],d[3]))
  
  xlab = "bp"
  ylab = "nuc"
  if (!is.null(opts[["xlab"]]))
    xlab = opts[["xlab"]]
  if (!is.null(opts[["ylab"]]))
    ylab = opts[["ylab"]]

  for (i in gnames)
  {
    tmat <- object$T2D[[i]] 
    rmat <- object$R2D[[i]]
    mat <- tmat*rmat
    z <- t(mat)
    x <- as.numeric(colnames(mat))
    y <- as.numeric(rownames(mat))
    call <- list(x,y,z, xlab=xlab, ylab=ylab, main=paste(i,"R"), useRaster=T)
    call <- c(call, graph_params(opts))
    do.call(image, call)
  }
}

plot_f <- function(object, type, setdim, opts)
{
  f <- c()
  if (type == "f")  f <- object$f
  if (type == "fa") f <- object$fa
  if (type == "F")  f <- object$F

  require(colorspace)
  gnames   <- object$gene_names
  tfnames  <- object$tf_names
  nucnames <- object$nuc_names
  
  if (!is.null(opts[["gene"]]))
  {
    gnames <- opts[["gene"]] # now you can set the gene with gene="gname"
    opts[["gene"]] <- NULL
  }
  if (!is.null(opts[["tf"]]))
  {
    tfnames <- eval(opts[["tf"]]) # now you can set the tf with tf="tfname"
    opts[["tf"]] <- NULL
  }
    
  ngenes <- length(gnames)
  ntfs   <- length(tfnames)
  nnuc   <- length(nucnames)

  n <- 0
  for (j in tfnames)
  {
    nmodes <- length(object$tfs[[j]]$coefs)
    if (type=="f") nmodes <- 1
    n      <- n + nmodes
  }
  
  d <- set_dim(ngenes*n)
  if (setdim) par(mfrow=c(d[2],d[3]))
  
  xlab = "bp"
  ylab = "nuc"
  if (!is.null(opts[["xlab"]]))
    xlab = opts[["xlab"]]
  if (!is.null(opts[["ylab"]]))
    ylab = opts[["ylab"]]

  
  for (i in gnames)
  {
    bindings <- object$bindings[[i]]
    gene_f   <- f[[i]]
    gene     <- object$genes[[i]]
    l        <- gene$length
    left     <- gene$left_bound
    
    for (j in tfnames)
    {
      nmodes <- length(object$tfs[[j]]$coefs)
      if (type=="f") nmodes <- 1
      tf_bindings <- bindings[which(bindings$tf==j),]
      for (m in 1:nmodes)
      {
        n <- j
        if (m>1)
          n <- paste0(n, m-1)
          
        tf_f <- gene_f[[n]]
        mat <- matrix(nrow=nnuc, ncol=l, data=0)
        ntf_bindings <- nrow(tf_bindings)
        if (ntf_bindings>0)
        {
          for (k in 1:ntf_bindings)
          {
            site <- tf_bindings[k,]
            start <- max(0, site$start - left + 1)
            end   <- min(site$end   - left + 1, l)
            
            occ <- tf_f[,(site$index)]
            mat[,start:end] <- mat[,start:end] + occ
          }
        }
        z <- t(mat)
        x <- seq(from=left, length.out=l)
        ynames <- nucnames
        y <- 1:length(nucnames)
        call <- list(x,y,z, xlab=xlab, ylab=ylab, main=paste(i,n,type), useRaster=T, zlim=c(0,1), yaxt='n')
        call <- c(call, graph_params(opts))
        do.call(image, call)
        axis(2, at=y, labels=nucnames, tck=0)
      }
    }
  }
}

get_f_matrix <- function(object, gene, tf)
{
  f        <- object$f
  gene_f   <- f[[gene]]
  gnames   <- object$gene_names
  tfnames  <- object$tf_names
  nucnames <- object$nuc_names
  nnuc     <- length(nucnames)
  bindings <- object$bindings[[gene]]
  g        <- object$genes[[gene]]
  l        <- g$length
  left     <- g$left_bound
  mat      <- matrix(nrow=nnuc, ncol=l, data=0)
  
  tf_bindings <- bindings[which(bindings$tf==tf),]
  nmodes <- length(object$tfs[[tf]]$coefs)
  
  for (m in 1:nmodes)
  {
    n <- tf
    if (m>1)
      n <- paste0(n, m-1)
      
    tf_f <- gene_f[[n]]
    ntf_bindings <- nrow(tf_bindings)
    if (ntf_bindings>0)
    {
      for (k in 1:ntf_bindings)
      {
        site <- tf_bindings[k,]
        start <- site$start - left + 1;
        end   <- site$end   - left + 1;
        start <- max(c(start, 1))
        end   <- min(c(end, l))
        occ <- tf_f[,(site$index)]
        mat[,start:end] <- mat[,start:end] + occ
      }
    }
  }
  return(mat)
}

get_F_matrix <- function(object, gene, tf)
{
  f        <- object$F
  gene_f   <- f[[gene]]
  gnames   <- object$gene_names
  tfnames  <- object$tf_names
  nucnames <- object$nuc_names
  mat      <- matrix(nrow=nnuc, ncol=l, data=0)
  bindings <- object$bindings[[gene]]
  g        <- object$genes[[gene]]
  l        <- g$length
  left     <- g$left_bound
  
  tf_bindings <- bindings[which(bindings$tf==tf),]
  nmodes <- length(object$tfs[[tf]]$coefs)
  
  for (m in 1:nmodes)
  {
    n <- tf
    if (m>1)
      n <- paste0(n, m-1)
      
    tf_f <- gene_f[[n]]
    ntf_bindings <- nrow(tf_bindings)
    if (ntf_bindings>0)
    {
      for (k in 1:ntf_bindings)
      {
        site <- tf_bindings[k,]
        start <- site$start - left + 1;
        end   <- site$end   - left + 1;
        occ <- tf_f[,(site$index)]
        mat[,start:end] <- mat[,start:end] + occ
      }
    }
  }
  return(mat)
}
  
get_F_matrices <- function(object, gene)
{
  f        <- object$F
  gene_f   <- f[[gene]]
  gnames   <- object$gene_names
  tfnames  <- object$tf_names
  nucnames <- object$nuc_names
  bindings <- object$bindings[[gene]]
  g        <- object$genes[[gene]]
  l        <- g$length
  left     <- g$left_bound
  mat      <- matrix(nrow=nnuc, ncol=l, data=0)
  
  tf_bindings <- bindings[which(bindings$tf==tf),]
  nmodes <- length(object$tfs[[tf]]$coefs)
  
  for (m in 1:nmodes)
  {
    n <- tf
    if (m>1)
      n <- paste0(n, m-1)
      
    tf_f <- gene_f[[n]]
    ntf_bindings <- nrow(tf_bindings)
    if (ntf_bindings>0)
    {
      for (k in 1:ntf_bindings)
      {
        site <- tf_bindings[k,]
        start <- site$start - left + 1;
        end   <- site$end   - left + 1;
        occ <- tf_f[,(site$index)]
        mat[,start:end] <- mat[,start:end] + occ
      }
    }
  }
  return(mat)
}

bcdhbf <- function(object)
{
  bcdf <- get_F_matrix(object, 'MSE2', 'bcd')
  hbf  <- get_F_matrix(object, 'MSE2', 'hb')
  bcdf <- bcdf * object$tfs$bcd$coefs[1]
  hbf  <- hbf *  object$tfs$hb$coefs[2]
  bcdhbf <- bcdf+hbf
  image(t(bcdhbf), col=rwb, useRaster=TRUE)
}
  

get_Q_matrix <- function(object, gene, tf)
{
  f        <- object$fa
  gene_f   <- f[[gene]]
  gnames   <- object$gene_names
  tfnames  <- object$tf_names
  nucnames <- object$nuc_names
  mat      <- matrix(nrow=nnuc, ncol=l, data=1)
  bindings <- object$bindings[[gene]]
  g        <- object$genes[[gene]]
  l        <- g$length
  left     <- g$left_bound
  
  tf_bindings <- bindings[which(bindings$tf==tf),]
  nmodes <- length(object$tfs[[tf]]$coefs)
  
  for (m in 1:nmodes)
  {
    coef <- object$tfs[[tf]]$coefs[m]
    if (coef >= 0) next
    coef <- -coef
        
    tf_f <- gene_f[[tf]]
    ntf_bindings <- nrow(tf_bindings)
    if (ntf_bindings>0)
    {
      for (k in 1:ntf_bindings)
      {
        
        site <- tf_bindings[k,]
        start <- site$start - left + 1
        end   <- site$end   - left + 1
        occ <- tf_f[,(site$index)]
        
        for (i in 1:l)
        {
          sdist <- abs(start-i)
          edist <- abs(end-i)
          d <- min(sdist, edist)
          dfun <- distfunc(d)
          if (i > start && i < end)
            dfun <- 1
          
          mat[,i] <- mat[,i]*(1 - coef*occ*dfun)
        }
      }
    }
  }
  return(mat)
}

get_Q_matrices <- function(object, gene)
{
  f        <- object$fa
  gene_f   <- f[[gene]]
  gnames   <- object$gene_names
  tfnames  <- object$tf_names
  nucnames <- object$nuc_names
  nnuc     <- length(nucnames)
  l        <- object$genes[[gene]]$length
  mat      <- matrix(nrow=nnuc, ncol=l, data=1)
  bindings <- object$bindings[[gene]]
  g        <- object$genes[[gene]]
  l        <- g$length
  left     <- g$left_bound
  
  for (tf in tfnames)
  {
    tf_bindings <- bindings[which(bindings$tf==tf),]
    nmodes <- length(object$tfs[[tf]]$coefs)
    for (m in 1:nmodes)
    {
      coef <- object$tfs[[tf]]$coefs[m]
      if (coef >= 0) next
      coef <- -coef
          
      tf_f <- gene_f[[tf]]
      ntf_bindings <- nrow(tf_bindings)
      if (ntf_bindings>0)
      {
        for (k in 1:ntf_bindings)
        {
          
          site <- tf_bindings[k,]
          start <- site$start - left + 1
          end   <- site$end   - left + 1
          occ <- tf_f[,(site$index)]
          
          for (i in 1:l)
          {
            sdist <- abs(start-i)
            edist <- abs(end-i)
            d <- min(sdist, edist)
            dfun <- distfunc(d)
            if (i > start && i < end)
              dfun <- 1
            
            mat[,i] <- mat[,i]*(1 - coef*occ*dfun)
          }
        }
      }
    }
  }
  return(mat)
}

distfunc <- function(distance, a=100, b=50)
{
  if (distance < 0)
    distance = -distance
  
  if (distance <= a)
    return(1)
  else if (distance < (a+b) )
  {
    x <- distance - a
    return(1-x/b)
  }
  else
    return(0)
}
  
  

# allow type of scoring to be set with scoretype=c("score","adj","gof", "K")
plot_scores <- function(object, type, setdim, opts)
{
  require(colorspace)
  gnames   <- object$gene_names
  tfnames  <- object$tf_names
  
  if (!is.null(opts[["gene"]]))
  {
    gnames <- opts[["gene"]] # now you can set the gene with gene="gname"
    opts[["gene"]] <- NULL
  }
  if (!is.null(opts[["tf"]]))
  {
    tfnames <- opts[["tf"]] # now you can set the tf with tf="tfname"
    opts[["tf"]] <- NULL
  }
  
  ngenes <- length(gnames)
  ntfs   <- length(tfnames)
  
  d <- set_dim(ngenes*ntfs)
  if (setdim) par(mfrow=c(d[2],d[3]))
  
  ylab = "Percent of Max Score"
  if (type=="score" | type=="adj")
    ylab = "log likelihood ratio"
  if (type=="GoF")
    ylab = "Percent of Max Score"
  if (type=="K")
    ylab = "Relative Affinity"
    
  xlab = "bp"
  if (!is.null(opts[["xlab"]]))
    xlab = opts[["xlab"]]
  if (!is.null(opts[["ylab"]]))
    ylab = opts[["ylab"]]
  
  cols <- c('red','blue')
  if (!is.null(opts[["col"]]))
  {
    cols = opts[['col']]
    opts[['col']] <- NULL
  }
  if (length(cols) <2)
    cols <- rep(cols[1],2)
    
  scores <- object$scores
  for (gene in gnames)
  {
    x <- as.numeric(rownames(scores[[gene]]$forward))
    for (tf in tfnames)
    {
      call <- list()
      
      maxscore <- object$tfs[[tf]]$pwm$max_score
      lambda   <- object$tfs[[tf]]$lambda
      
      fscore <- scores[[gene]]$forward[[tf]]
      rscore <- scores[[gene]]$reverse[[tf]]
      
      if (type=="GoF")
      {
        fscore <- fscore / maxscore
        rscore <- rscore / maxscore
        fscore[fscore < 0] <- 0
        rscore[rscore < 0] <- 0
        if (is.null(opts[["ylim"]]))
          call <- c(call, list(ylim=c(0,1)))
      }
      if (type=="K")
      {
        fscore <- exp( (fscore-maxscore)/lambda)
        rscore <- exp( (rscore-maxscore)/lambda)
        
        if (is.null(opts[["ylim"]]))
          call <- c(call, list(ylim=c(0,1)))
      }
        
      call <- c(call, list(x=x, y=fscore, type='l',xlab=xlab,ylab=ylab,main=paste(tf,"on",gene), col=cols[1]))
      call <- c(call, opts)
      do.call(plot, call)
      call <- list(x=x, y=rscore, col=cols[2])
      call <- c(call, opts)
      do.call(lines, call)
    }
  }
}
      
partial_derivative <- function(o, gene, tf, h)
{
  rate      <- o$rate
  row_names <- rownames(rate)
  rate0     <- rate[,gene]
  
  df   <- o$tf_data$data_frame
  
  tfph <- df[row_names, tf] + h
  tfnh <- df[row_names, tf] - h
  
  use_forward <- which(tfnh <= 0)
  tfnh[use_forward] <- 0
  
  dfph <- df
  dfnh <- df
  
  dfph[row_names, tf] <- tfph
  dfnh[row_names, tf] <- tfnh
  
  o$tf_data$set(dfph, "ID", "TF")
  o$recalculate()
  
  rate_p <- o$rate[,gene]
  
  o$tf_data$set(dfnh, "ID", "TF")
  o$recalculate()
  
  rate_n <- o$rate[,gene]
  
  o$tf_data$set(df, "ID", "TF")
  o$recalculate()
  
  pder <- c()
  pder <- (rate_p - rate_n)/(2*h)
  pder[use_forward] <- (rate_p[use_forward] - rate0[use_forward])/h
 
  return(pder)
}

deltas <- function(o, gene, tf)
{
  rate      <- o$rate
  row_names <- rownames(rate)
  
  df   <- o$tf_data$data_frame
  tf   <- df[row_names, tf]
  
  n <- length(row_names)
  
  delta <- c()
  i <- 1
  delta[i] = tf[i] - tf[i+1]
  delta[n] = tf[n-1] - tf[n]
  for (i in 2:(n-1)) delta[i] = (tf[i-1] - tf[i+1])/2
  
  return(delta)
}

plot_border <- function(o, gene)
{
  cols <- rev(brewer.pal(9, "Paired"))
  row_names <- o$tf_names
  ntfs      <- length(row_names)
  
  rate      <- o$rate
  col_names <- rownames(rate)
  nnuc      <- length(col_names)

  m <- matrix(nrow=ntfs, ncol=nnuc)
  for (i in 1:ntfs)
    m[i,] <- partial_derivative(o, gene, row_names[i], 1) * deltas(o, gene, row_names[i])
    
  mpos <- m
  mneg <- m
  mpos[mpos<0] <- 0
  mneg[mneg>0] <- 0
  mpos <- mpos/2.55
  mneg <- mneg/2.55
  
  ymin <- min(colSums(mneg))
  ymax <- max(colSums(mpos))
  
  lim <- max(abs(ymin),ymax)
  lim <- lim+0.04*lim
  
  plot.new()
  plot.window(xlim=c(35.5,92.5),ylim=c(-lim,lim), xaxs='i', yaxs='i', bty='o')
  axis(1, tck=0.02)
  axis(2, tck=0.02)
  axis(3, tck=0.02, labels=FALSE)
  axis(4, tck=0.02, labels=FALSE)
  box(lwd=2)
  #plot(mpos[1,], type='l', ylim=c(ymin, ymax))
  #lines(mneg[1,])
  polygon(x=c(35.5:92.5,92.5:35.5),y=c(mpos[1,],rep(0,nnuc)),border=NA, col=cols[1])
  polygon(x=c(35.5:92.5,92.5:35.5),y=c(mneg[1,],rep(0,nnuc)),border=NA, col=cols[1])
  for (i in 2:ntfs)
  {
    #lines(colSums(mpos[1:i,]))
    #lines(colSums(mneg[1:i,]))
    polygon(x=c(35.5:92.5,92.5:35.5),y=c(colSums(mpos[1:i,]),rev(colSums(matrix(mpos[1:(i-1),], ncol=nnuc)))),border=NA, col=cols[i])
    polygon(x=c(35.5:92.5,92.5:35.5),y=c(colSums(mneg[1:i,]),rev(colSums(matrix(mneg[1:(i-1),], ncol=nnuc)))),border=NA, col=cols[i])

  }
  legend(x=70, y=5, fill=cols, legend=row_names) 
}

plot_activation <- function(o, gene)
{
  require(RColorBrewer)
  # which parameters are activation parameters
  coefs <- which(o$parameter_table$name=='coef')
  
  coef_acts <- c()
  tfs       <- c()
  
  for (i in coefs)
  {
    param <- o$parameters[[i]]
    if (param$value >= 0)
    {
      coef_acts <- c(coef_acts,i)
      tfs       <- c(tfs, param$tf)
    }
  }
  n <- length(coef_acts)
  
  col_names <- rownames(o$rate)
  nnuc      <- length(col_names)
  
  m <- matrix(nrow=n, ncol=nnuc)
  t <- o$T2D$eve_mel
  N <- as.numeric(rowSums(o$N2D$eve_mel * t))
  R <- o$rate$eve_mel/2.55
  # are we using competition equations?
  for (i in 1:length(coef_acts))
  {
    index <- coef_acts[[i]]
    param <- o$parameters[[index]]
    param$value <- 0
    m[i,] <- (N-as.numeric(rowSums(o$N2D$eve_mel * t))) / N
    m[i,] <- m[i,]*R
    param$restore()
  }
  
  cols <- rev(brewer.pal(n, "Paired"))
  
  ymax <- max(colSums(m))
  plot.new()
  plot.window(xlim=c(35.5,92.5),ylim=c(0,ymax), xaxs='i', yaxs='i', bty='o')
  axis(1, tck=0.02)
  axis(2, tck=0.02)
  axis(3, tck=0.02, labels=FALSE)
  axis(4, tck=0.02, labels=FALSE)
  box(lwd=2)
  polygon(x=c(35.5:92.5,92.5:35.5),y=c(m[1,],rep(0,nnuc)),border=NA, col=cols[1])
  for (i in 2:n)
    polygon(x=c(35.5:92.5,92.5:35.5),y=c(colSums(m[1:i,]),rev(colSums(matrix(m[1:(i-1),], ncol=nnuc)))),border=NA, col=cols[i])
    
    
}

    
  



random_sequence <- function(n, gc=0.5)
{
  at    <- 1-gc
  gc2   <- gc/2
  at2   <- at/2
  bases <- c("A","C","G","T")
  probs <- c(at2,gc2,gc2,at2)
  s     <- sample(bases, size=n, prob=probs, replace=TRUE)
  return(s)
}
          
get_rate_from_seq <- function(seq, model_file)
{
  require(Rtranscpp)
  
  m <- new(Organism, model_file)
  s <- m$genes[[1]]
  s$sequence <- seq
  m$recalculate()
  return(m$rate[,1])
}

get_rate_from_seqs <- function(seqs, model_file)
{
  require(Rtranscpp)
  
  m <- new(Organism, model_file)
  s <- m$genes[[1]]
  
  out <- c()
  for (i in seqs)
  {
    s$sequence <- i
    m$recalculate()
    out <- rbind(out, m$rate[,1])
  }
  return(out)
}
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
    
