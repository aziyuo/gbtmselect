% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/hello.R
\name{gbtmselect}
\alias{gbtmselect}
\title{gbtmselect}
\usage{
gbtmselect(files_oe, files_of, files_os)
}
\arguments{
\item{files_oe}{string vectors containing paths of GBTM outcome "oe" files.}

\item{files_of}{string vectors containing paths of GBTM outcome "of" files.}

\item{files_os}{string vectors containing paths of GBTM outcome "os" files.}
}
\description{
As we all know, SAS "Traj" program can produce lots of files since we need to iterate our models. You can use this package to easily address the outcomes can converge them.
What you should do is nameing you outcome files as "oe/of/os_group_polynomial degree". Groups may be 1 to 5, and polynomial degree depends on the group numbers such as 33.
"oe" means parameter and covariance matrix estimates. "of" means group assignments and membership probabilities. "os" means parameter estimates used by TRAJPLOT macro.
More information can be found in "https://www.andrew.cmu.edu/user/bjones/index.htm".
}
\examples{
gbtmselect(files_oe,files_of,files_os)
files_oe=c("oe//oe_1_1.sas7bdat","oe//oe_1_2.sas7bdat","oe//oe_1_3.sas7bdat","oe//oe_2_11.sas7bdat","oe//oe_2_22.sas7bdat","oe//oe_2_33.sas7bdat","oe//oe_3_011.sas7bdat","oe//oe_3_221.sas7bdat","oe//oe_3_333.sas7bdat","oe//oe_4_0111.sas7bdat","oe//oe_4_2221.sas7bdat","oe//oe_4_2331.sas7bdat")
files_of=c("of//of_1_1.sas7bdat","of//of_1_2.sas7bdat","of//of_1_3.sas7bdat","of//of_2_11.sas7bdat","of//of_2_22.sas7bdat","of//of_2_33.sas7bdat","of//of_3_011.sas7bdat","of//of_3_221.sas7bdat","of//of_3_333.sas7bdat","of//of_4_0111.sas7bdat","of//of_4_2221.sas7bdat","of//of_4_2331.sas7bdat")
files_os=c("os//os_1_1.sas7bdat","os//os_1_2.sas7bdat","os//os_1_3.sas7bdat","os//os_2_11.sas7bdat","os//os_2_22.sas7bdat","os//os_2_33.sas7bdat","os//os_3_011.sas7bdat","os//os_3_221.sas7bdat","os//os_3_333.sas7bdat","os//os_4_0111.sas7bdat","os//os_4_2221.sas7bdat","os//os_4_2331.sas7bdat")
}
