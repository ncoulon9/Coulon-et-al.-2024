#################################
## CALCULATE PROCESS INTENSITY ##
#################################

tropicalization <- subset(species_change, species_change$change > 0 & species_change$sst_median > CTI_mean)
tropicalization$sst_median_diff <- tropicalization$sst_median - CTI_mean
tropicalization$sst_median_wtd_change <- tropicalization$change * tropicalization$sst_median_diff
tropicalization$process <- "tropicalization"
trop_intensity <- sum(abs(tropicalization$sst_median_wtd_change))

deborealization <- subset(species_change, species_change$change < 0 & species_change$sst_median < CTI_mean)
deborealization$sst_median_diff <- deborealization$sst_median - CTI_mean
deborealization$sst_median_wtd_change <- deborealization$change * deborealization$sst_median_diff
deborealization$process <- "deborealization"
deb_intensity <- sum(abs(deborealization$sst_median_wtd_change))

borealization <- subset(species_change, species_change$change > 0 & species_change$sst_median < CTI_mean)
borealization$sst_median_diff <- borealization$sst_median - CTI_mean
borealization$sst_median_wtd_change <- borealization$change * borealization$sst_median_diff
borealization$process <- "borealization"
bor_intensity <- sum(abs(borealization$sst_median_wtd_change))

detropicalization <- subset(species_change, species_change$change < 0 & species_change$sst_median > CTI_mean)
detropicalization$sst_median_diff <- detropicalization$sst_median - CTI_mean
detropicalization$sst_median_wtd_change <- detropicalization$change * detropicalization$sst_median_diff
detropicalization$process <- "detropicalization"
detrop_intensity <- sum(abs(detropicalization$sst_median_wtd_change))

