library(ggplot2)
library(plyr)
head(slr_all)

setwd("~/Documents/projects/slr_pipeline")

scratch_predictions <- read.table("data/scratch_all.tab", sep="\t", header=T, stringsAsFactors=F)
scratch_predictions$ens_pos <- scratch_predictions$human_idx
head(scratch_predictions)

# Match with the PDB structures to see what the accuracy is
str_predictions <- join(pdb_master_globular, scratch_predictions)
head(str_predictions)
str_predictions <- subset(str_predictions, !is.na(human_idx))
str_predictions$buried[str_predictions$buried == TRUE] <- "buried"
str_predictions$buried[str_predictions$buried == FALSE] <- "exposed"
1 - sum(str_predictions$buriedness != str_predictions$buried) /
sum(str_predictions$buriedness == str_predictions$buried)


str_predictions$ss3 <- NA
str_predictions$ss3[str_predictions$sec_simple == "Helix"] <- "H"
str_predictions$ss3[str_predictions$sec_simple == "Beta sheet"] <- "E"
str_predictions$ss3[str_predictions$sec_simple %in% c("Turn", "Coil") ] <- "C"
str_predictions$ss3 <- factor(str_predictions$ss3, levels=c("H", "E", "C"))
1 - sum(str_predictions$ss_pred != str_predictions$ss3) /
  sum(str_predictions$ss_pred == str_predictions$ss3)

slr_predictions <- join(slr_all, scratch_predictions)
slr_predictions <- subset(slr_predictions, !is.na(lower))
head(slr_predictions)

pred_fractions <- make.fractions(slr_predictions, c("ss_pred", "buriedness_pred"))
pred_fractions_lower <- make.fractions.lower(slr_predictions, c("ss_pred", "buriedness_pred"))

ggplot(pred_fractions_lower, aes(x=buriedness_pred, y=Fraction, fill=ss_pred)) +
  geom_bar(stat="identity", position=position_dodge()) +
  ggtitle("Positive selection and structure") +
  theme_bw(base_size=18)

# Get the domain and transmembrane annotations
head(subset(slr_predictions, is.na(ss_pred)))


