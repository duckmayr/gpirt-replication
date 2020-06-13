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
