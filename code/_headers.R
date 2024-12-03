# Install missing packages and load necessary packages
options(repos = c(CRAN = "http://cran.us.r-project.org"))
install.packages("xfun")
xfun::pkg_attach2(required_packages, message = FALSE)

