
## Up until R 2.15.0, the require("methods") is needed but (now)
## triggers an warning from R CMD check
#.onLoad <- function(libname, pkgname){
#    #require("methods")  ## needed with R <= 2.15.0
#    loadRcppModules()
#}


## For R 2.15.1 and later this also works. Note that calling loadModule() triggers
## a load action, so this does not have to be placed in .onLoad() or evalqOnLoad().
loadModule("mod_organism", TRUE)
loadModule("mod_datatable", TRUE)
loadModule("mod_mode", TRUE)
loadModule("mod_par_defaults", TRUE)
loadModule("mod_gene", TRUE)
loadModule("mod_tf", TRUE)
loadModule("mod_pwm", TRUE)
loadModule("mod_param", TRUE)

evalqOnLoad({
  
setMethod("show", Organism, function(object)
{
  cat(paste0("\n\tOrganism Object\n\n"))
  
  nspaces <- max(nchar(ls(object)))
  for (i in ls(object))
    cat(paste0("\t$",formatC(i, width=-nspaces), "\n"))
  cat("\n")
})

setMethod("show", PWMParameter, function(object)
{
  cat(paste0("PWM Parameter Object\n"))
  cat("$name\n")
  print(object$name)
  cat("$type\n")
  print(object$type)
  cat("$tf\n")
  print(object$tf)
  cat("$value\n")
  print(object$value)
})

setMethod("show", DoubleParameter, function(object)
{
  cat(paste0("Double Parameter Object\n"))
  cat("$name\n")
  print(object$name)
  cat("$type\n")
  print(object$type)
  cat("$tf\n")
  print(object$tf)
  cat("$value\n")
  print(object$value)
})

setMethod("show", SeqParameter, function(object)
{
  cat(paste0("Sequence Parameter Object\n"))
  cat("$name\n")
  print(object$name)
  cat("$type\n")
  print(object$type)
  cat("$tf\n")
  print(object$tf)
  cat("$value\n")
  print(object$value)
})

})
