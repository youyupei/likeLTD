## Test unit 'relatedness'
library(svUnit)

###############################################################
# The new two functions are to set up the unit test environment
###############################################################

.setUp <-
function () {
	## Specific actions for svUnit: prepare context
	if ("package:svUnit" %in% search()) {
		.Log <- Log() ## Make sure .Log is created
		.Log$..Unit <- "inst/unitTests/runit_relatedness.R"
		.Log$..File <- ""
		.Log$..Obj <- ""
		.Log$..Tag <- ""
		.Log$..Msg <- ""
		rm(..Test, envir = .Log)
	}
}

.tearDown <-
function () {
	## Specific actions for svUnit: clean up context
	if ("package:svUnit" %in% search()) {
		.Log$..Unit <- ""
		.Log$..File <- ""
		.Log$..Obj <- ""
		.Log$..Tag <- ""
		.Log$..Msg <- ""
		rm(..Test, envir = .Log)
	}
}


test_relatedness <- svTest(function()  {

  if(! "relatedness" %in% ls(.GlobalEnv))
    relatedness = getFromNamespace("relatedness", "likeLTD")
  if(! "all.genotypes.per.locus" %in% ls(.GlobalEnv))
    all.genotypes.per.locus = getFromNamespace("all.genotypes.per.locus", "likeLTD")

  A = all.genotypes.per.locus(5, 2)
  result = relatedness(5, 2, A, c(2, 3))
  checkTrue( all(apply(result[[1]], 2, function(n) 2 %in% n[1])) )
  checkEquals( ncol(result[[1]]), 75)
  checkTrue( all(apply(result[[2]], 2, function(n) 3 %in% n[1])) )
  checkEquals( ncol(result[[2]]), 75)
  checkTrue( result[[3]][1, ] == 2)
  checkTrue( result[[3]][2, ] == 3)
  checkEquals( ncol(result[[3]]), 15)
})
