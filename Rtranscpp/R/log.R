files <- 'fits/locus_7stripe/21.xml.log'
update=0
logplot=F
variable="mu"
n=1000


# the log functiont takes the xml file and then decides which extention to use
plot_log <- function(files, update=0, logplot=F, variable="mu", n=1000, schedule="lam", ...)
{
  lam_headers <- c("steps", "1/T", "dS", "mu", "sigma", "E(mu)", "E(sigma)", "acceptance ratio", "alpha")
  exp_headers <- c("steps", "1/T", "mu")
  opts <- as.list(substitute(list(...)))[-1L]
  
  nfiles <- length(files)
  d <- set_dim(nfiles)
  
  last <- list()
  
  repeat
  {
    par(mfrow=c(d[2],d[3]))
    for (i in files)
    {
      m <- n
      cmd <- paste('wc -l', i)
      nlines <- read.table(pipe(cmd))$V1[1]
      m <- min(c(m, nlines))
      cmd <- paste('tail','-n',formatC(m, format="d"), i)
      tmp <- try(read.table(pipe(cmd)),silent=T)
      if (!inherits(tmp, 'try-error'))
      {
        if (schedule=="lam") {
          colnames(tmp) <- lam_headers
        } else if (schedule=="exp") {
          colnames(tmp) <- exp_headers
        }
        
        nrow = nrow(tmp)
        x = tmp[,"steps"]
        y = tmp[,variable]
        
        best <- min(y);
        
        ylab=variable
        
        if (logplot) 
        {
          y    <- log(y, 10)
          ylab <- paste0('log_10(',variable,')')
        }
        
        ylim <- c(0, max(y[is.finite(y)]))
        if (!is.null(opts[["ylim"]]))
          ylim <- opts[["ylim"]]
        call <- list(x=x,y=y,xlab="steps",ylab=ylab,type='l',ylim=ylim,main=i)
        call <- c(call, graph_params(opts))
        do.call(plot, call)
        text(x=(max(x)+min(x))/2, y=best/2, labels=best)
        steps <- tmp[nrow,"steps"]
        if (!is.null(last[[i]]))
        {
          steps_per_sec <- (steps - last[[i]])/update
          title(sub=paste(steps_per_sec,"/ s"))
        }
        last[[i]] <- steps
      }
    }
    if (!update) break
    Sys.sleep(update);
  } 
}

get_active_runs <- function(dir='./', wait=10, suffix=".log", recursive=FALSE)
{
  # get all the log files in the dir
  fnames <- list.files(dir,pattern=suffix, recursive=recursive)
  files <- paste0(dir,'/',fnames)
  steps_start <- list()
  steps_end   <- list()
  for (i in files)
  {
    tmp <- try(read.table(i),silent=T)
    if (inherits(tmp, 'try-error'))
      steps_start[[i]] <- 0
    else
      steps_start[[i]] <- tmp[nrow(tmp),1]
  }
  Sys.sleep(wait)
  for (i in files)
  {
    tmp <- try(read.table(i),silent=T)
    if (inherits(tmp, 'try-error'))
      steps_end[[i]] <- 0
    else
      steps_end[[i]] <- tmp[nrow(tmp),1]
  }
  out <- c()
  for (i in files)
  {
    if (steps_end[[i]] > steps_start[[i]])
      out <- c(out, substring(i, 1, nchar(i)-4))
  }
  return(out)
}

    
      
f <- function(...)
{
  opts <- as.list(substitute(list(...)))[-1L]
  print(opts)
}
