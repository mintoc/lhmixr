##--------------------
## make the outputted data for the package
## from the data Rui sent
## CM: Tue Sep 27 2016
##--------------------
library(gdata)

##----------
## E.SPINAX
##----------
Espinax <- read.xls("Etmopterus_spinax_dataset.xlsx", stringsAsFactors = FALSE)

## change variable names
names(Espinax)[names(Espinax) == "Species"] <- "species"
names(Espinax)[names(Espinax) == "Sex"] <- "sex"
names(Espinax)[names(Espinax) == "Final_age"] <- "age"
names(Espinax)[names(Espinax) == "TL_cm"] <- "length"
names(Espinax)[names(Espinax) == "Maturity"] <- "maturity"

Espinax$maturity[Espinax$maturity == "Juvenile"] <- "immature"
Espinax$maturity[Espinax$maturity == "Adult"] <- "mature"

vars2keep <- c("species", "sex", "age", "length", "maturity")

Espinax <- Espinax[, vars2keep]

save("Espinax", file = "../data/Espinax.RData")


##----------
## E.PUSILLUS
##----------
Epusillus <- read.xls("Etmopterus_pusillus_dataset.xlsx", stringsAsFactors = FALSE)

## change variable names
names(Epusillus)[names(Epusillus) == "Species"] <- "species"
names(Epusillus)[names(Epusillus) == "Sex"] <- "sex"
names(Epusillus)[names(Epusillus) == "Final_age"] <- "age"
names(Epusillus)[names(Epusillus) == "TL_cm"] <- "length"
names(Epusillus)[names(Epusillus) == "Maturity"] <- "maturity"

Epusillus$maturity[Epusillus$maturity == "Juvenile"] <- "immature"
Epusillus$maturity[Epusillus$maturity == "Adult"] <- "mature"

Epusillus$length  <- round(Epusillus$length, 1)

vars2keep <- c("species", "sex", "age", "length", "maturity")

Epusillus <- Epusillus[, vars2keep]

save("Epusillus", file = "../data/Epusillus.RData")


