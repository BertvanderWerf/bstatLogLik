#install.packages("usethis")
#library(usethis)
usethis::use_readme_rmd()       # Creates a README.Rmd
usethis::use_mit_license()      # Or another license
usethis::use_roxygen_md()       # For documentation
usethis::use_testthat()         # For testing package functions
usethis::use_git()              # If Git is not already initialized

# already done
#usethis::use_github()          # Connects and pushes to GitHub

usethis::use_package_doc()
usethis::use_vignette("getting_started.Rmd")

usethis::use_package_doc()
usethis::use_package("Deriv")
usethis::use_package("MASS")

usethis::use_dev_package("bstatErr")
usethis::use_dev_package("bstatUtils")

usethis::use_r("create_parameter_info")
usethis::use_test("create_parameter_info")

usethis::use_r("evaluate_function")
usethis::use_test("evaluate_function")

usethis::use_r("loglik")
usethis::use_test("loglik")

usethis::use_r("loglik_control")
usethis::use_test("loglik_control")

usethis::use_r("loglik_derivatives")
usethis::use_test("loglik_derivatives")

usethis::use_r("loglik_utilities")
usethis::use_test("loglik_utilities")

usethis::use_r("loglik_sandwich")
usethis::use_test("loglik_sandwich")

usethis::use_r("loglik_methods")
usethis::use_test("loglik_methods")

usethis::use_r("distributions")

usethis::use_r("distributions_weibull")
usethis::use_test("distributions_weibull")

usethis::use_r("distributions_gamma")
usethis::use_test("distributions_gamma")

usethis::use_r("distributions_normal")
usethis::use_test("distributions_normal")

usethis::use_r("distributions_lognormal")
usethis::use_test("distributions_lognormal")

usethis::use_r("loglik_inference")
usethis::use_test("loglik_inference")

usethis::use_r("parameter_info_methods")
usethis::use_test("parameter_methods")


#usethis::rename_files()
