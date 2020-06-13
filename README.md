# gpirt-replication

This repository contains the code and raw data necessary to reproduce the results in "GPIRT: A Gaussian Process Model for Item Response Theory" by JBrandon Duck-Mayr, Roman Garnett, and Jacob M. Montgomery.

The "data" folder should contain three files:
  - "npi.csv", which contains responses to the 40-item Narcisistic Personality Inventory (Raskin and Terry 1988):
    + 10,440 responses are from the Open Source Psychometrics Project (Open-Source Psychometrics Project 2012)
    +  2,945 responses are from a convenience sample from the Qualtrics panel
  - "H116_votes.csv," which contains members' roll call votes in the 116th House of Representatives, obtained from Voteview.com (Lewis et al. 2020) on January 4, 2020
  - "H116_members.csv," which contains information about members of the 116th House of Representatives, again obtained from Voteview.com (Lewis et al. 2020) on January 4, 2020

The "code" folder should contain eleven files:
  - "covSEiso.R": This script contains an R implementation of the squared exponential covariance function.
    + This function was simply used in plotting example IRFs in Figure 1.
    + The script "figure1-example-IRFs.R" relies on this file, but sources it in; there is no need for the user to run this script manually.
  - "figure1-example-IRFs.R": This script contains the code necessary to reproduce Figure 1.
  - "fix_directories.R": This script contains an R function to ensure subdirectories exist to save plots and model output to.
    + All the application reproduction scripts rely on this file, but source it in themselves; there is no need for the user to run this script manually.
    + Note this means that running our analyses will create subdirectories in the project's folder.
  - "prepare-data.R": This script processes the raw data into the formats analyzed in the paper.
    + Namely, the Congress roll call responses need to be reshaped into a matrix where the rows represent legislators and the columns represent roll call votes, and votes and members that do not meet certain criteria need to be removed, as described in the paper.
    + We also split the NPI data into training and adaptive testing subsamples in this file.
  - "plot_irf.R": This script contains an R function to plot IRFs with a particular aesthetic.
    + The script "congress-application.R" relies on this file, but sources it in; there is no need for the user to run this script manually.
  - "congress-application.R": This script reproduces the figures presented in the main body of the paper for the U.S. Congress application (Section 5.1).
  - "gpirt-cat.cpp": This script contains the main C++ code to implement the active learning algorithm presented in Section 3.2.
    + The script "npi-application.R" relies on this file, but a command there compiles this code and loads the resulting shared object; there is no need for the user to compile this script manually.
    + The script "supplemental-results.R" also relies on this file, but compiles and loads it within the script itself.
  - "cat-helpers.R": This script contains some R helper functions for the active learning experiment.
    + The scripts "npi-application.R" and "supplemental-results.R" rely on this file, but source it in; there is no need for the user to run this script manually.
  - "gplvm.py": This script trains and obtains predictions for the GPLVM model for the held-out experiment in Section 5.2 (see Table 1).
  - "npi-application.R": This script reproduces the results reported in Tables 1 and 2 (relying in part on the results of the "gplvm.py" python script).
  - "supplemental-results.R": This script reproduces the results reported in Table A1 and the text of Appendix A in the Supplementary Material.

To reproduce all results presented in the paper, you must:
  - Run the R script "code/prepare-data.R". This script *must* be run before anything else.
  - Run the R script "code/figure1-example-IRFs.R"
  - Run the python script "code/gplvm.py". This script at the least must be run before "code/npi-application.R".
  - Run the R script "code/npi-application.R"
  - Run the R script "code/congress-application.R"
  - Run the R script "code/supplemental-results.R"

In addition to R (versions 3.6.3 and 4.0.0 were used by the authors) and python (versions 3.8.2 and 3.8.3 were used by the authors), you will need:
  - Sufficient tools to compile C++ code (g++ versions 9.3.0 and 10.1.0 were used by the authors)
  - The python package "numpy" (version 1.18.2 was used by the authors)
  - The python package "GPy" (version 1.9.9 was used by the authors)
  - The R package "bggum" (version 1.0.2 was used by the authors)
  - The R package "ltm" (version 1.1-1 was used by the authors)
  - The R package "KernSmoothIRT" (version 6.4 was used by the authors)
  - The R package "Rcpp" (version 1.0.4.6 was used by the authors)
  - The R package "RcppArmadillo" (version 0.9.900.1.0 was used by the authors)
  - The R package "devtools" (version 2.3.0 was used by the authors)
  - The R package "mvtnorm" (version 1.1-0 was used by the authors)
  - The R package "dplyr" (version 0.8.5 was used by the authors)
  - The R package "tidyr" (version 1.0.3 was used by the authors)
  - The R package "pROC" (version 1.16.2 was used by the authors)
  - The R package "emIRT" (version 0.0.11 was used by the authors)

Please note that additional R or python packages, or system dependencies, may be needed to satisfy the dependencies of the above listed packages directly relied on by the authors.
To the extent that exact results may be affected by operating system differences, you may note the authors conducted analysis on Linux systems (Manjaro Linux and Ubuntu Linux).


References

Lewis, J. B., Poole, K., Rosenthal, H., Boche, A., Rudkin, A., and Sonnet, L. Voteview: Congressional roll-call votes database, 2020. URL: https://voteview.com/.

Open-Source Psychometrics Project. Open psychology data: Raw data from online personality tests, 2012. URL: https://openpsychometrics.org/_rawdata/.

Raskin, R. and Terry, H. A principal-components analysis of the narcissistic personality inventory and further evidence of its construct validity. Journal of Personality and Social Psychology, 54(5):890–902, 1988.
