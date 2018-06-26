# LoadPackage

LoadPackage = function(package.name){
	# LoadPackage 
	#
	# This function loads a package into the current workspace.
	# If the package is not installed, then LoadPackage will attempt to install it. 
	# The installed package is then silently loaded

	if(!require(package.name, character.only = TRUE, quietly = TRUE, warn.conflicts = FALSE)){
		install.packages(package.name)
		library(package.name, character.only = TRUE, quietly = TRUE, warn.conflicts = FALSE)
	}
	return()
}

AutoloadPackages = function(vector.of.package.names){
	# AutoloadPackages
	#
	# This function calls LoadPackage on a vector of package names.
	#
	# Args:
	#	auto.loads: a vector of package names, e.g. c("ggplot2", "MASS") 
	#
	# Output: NULL

	invisible(sapply(vector.of.package.names, LoadPackage))
	return()
}
