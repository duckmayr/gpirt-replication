##### Setup -----
## Load required packages
library(dplyr)
library(tidyr)
## If you do not have gpirt installed, you can do so via
# library(devtools)
# install_github("duckmayr/gpirt")
library(gpirt)


##### Prepare 116th House of Representatives Data -----
## Read in the Voteview data (retrieved January 4, 2020)
votes   <- read.csv("data/H116_votes.csv")
members <- read.csv("data/H116_members.csv")
## Make a response matrix from the Voteview data
responses  <- votes %>%
    select(icpsr, rollnumber, cast_code) %>% ## trim to only variables we need
    mutate(cast_code = ifelse(cast_code == 6, -1, ## Put votes in {NA, -1, 1}
                              ifelse(cast_code == 1, 1, NA))) %>%
    spread(icpsr, cast_code) %>% ## spread the data out
    as.matrix() ## Convert to matrix
rownames(responses) <- responses[ , "rollnumber"] ## Deal with dimnames
responses <- t(responses[ , -1]) ## Put votes as columns, MCs as rows
party_codes <- members$party_code[match(rownames(responses), members$icpsr)]
## Get rid of Beaman, who wasn't in any non-unanimous votes
no_votes  <- apply(responses, 1, function(x) all(is.na(x)))
responses <- responses[!no_votes, ]
## Get rid of Amash
amash_rows <- which(grepl("21143|91143", rownames(responses)))
responses <- responses[-amash_rows, ]
## Remove overly lopsided votes
is_lopsided <- function(x, cutoff = 0.025) {
    tab  <- table(x)            ## Tabulate the votes
    prop <- min(tab) / sum(tab) ## Get proportion of minority votes
    return(prop < cutoff)       ## Return TRUE if that proporion is < cutoff
}
lopsided_votes <- which(apply(responses, 2, is_lopsided, cutoff = 0.01))
responses <- responses[ , -lopsided_votes]
## Convert to response_matrix class and save to disk
vote_codes <- list(yea = 1, nay = -1, missing = NA)
responses  <- as.response_matrix(responses, response_codes = vote_codes)
saveRDS(responses, file = "data/H116-responses.rds")


##### Prepare NPI Data -----
## Read in full data, eliminating sample ID
## (We don't care if the response is from Open Psychometrics or Qualtrics)
full_responses <- read.csv("data/npi.csv") %>% select(-sample)
## We'll randomly take 2,000 of them (setting the seed for reproducibility)
set.seed(314)
train_idx <- sample(x = 1:nrow(full_responses), size = 2000)
responses <- full_responses[train_idx, ]
## Convert to response_matrix class and save to disk
codes     <- list(yea = 1, nay = 0, missing = NA)
responses <- as.response_matrix(x = responses, response_codes = codes)
saveRDS(responses, file = "data/npi-responses.rds") ## Use RDS for R analyses
write.table(responses, file = "data/npi-responses.csv", row.names = FALSE,
          col.names = FALSE, sep = ",") ## CSV for python
## Now we randomly take another 1,000 for the CAT experiment
test_idx  <- sample(x = setdiff(1:nrow(full_responses), train_idx), size = 1000)
responses <- full_responses[test_idx, ]
## Convert to response_matrix class and save to disk
responses <- as.response_matrix(x = responses, response_codes = codes)
saveRDS(responses, file = "data/npi-cat-responses.rds")
