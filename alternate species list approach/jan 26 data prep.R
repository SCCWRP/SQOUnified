# sand box to play around with tuning the alt-sqo code

output_path="C:/Users/davidg/SCCWRP/Bight 2023 - Documents/Infauna/index code sub committee/B23 Index Coding/alternate species list approach/alt sqo outputs"
file_id="jan-26 test runs"
infauna_path="C:/Users/davidg/SCCWRP/Bight 2023 - Documents/Infauna/index code sub committee/B23 Index Coding/example SQO data/example SQO data.csv"
EG_Scheme="Hybrid"


output_path="C:/Users/GisUser/SCCWRP/Bight 2023 - Infauna/index code sub committee/B23 Index Coding/alternate species list approach/alt sqo outputs"
file_id="jan-26 test runs"
infauna_path="C:/Users/GisUSer/SCCWRP/Bight 2023 - Infauna/index code sub committee/B23 Index Coding/example SQO data/example SQO data.csv"
EG_Scheme="Hybrid"
#orig.socal_sqo<-read.csv("C:/Users/davidg/SCCWRP/David Gillett's Marine-Estuarine Benthic Taxonomy Database - General/Marine Taxonomy Database/individual files for R stuff/Ref - SoCal SQO Species List.csv")

#save(file="Reference Files/SoCal SQO LU.RData", orig.socal_sqo) 

#load("alternate species list approach/Reference Files/SoCal SQO LU.RData")


BenthicData=BenthicData.3


setwd("C:/Users/davidg/SCCWRP/Bight 2023 - Documents/Infauna/index code sub committee/B23 Index Coding/alternate species list approach")



#####

third.time<- alt.SQO_BLOE.generic(file_id="jan-26 test runs", 
                                            infauna_path="C:/Users/davidg/SCCWRP/Bight 2023 - Documents/Infauna/index code sub committee/B23 Index Coding/example SQO data/example SQO data.csv", 
                                            output_path="C:/Users/davidg/SCCWRP/Bight 2023 - Documents/Infauna/index code sub committee/B23 Index Coding/alternate species list approach/alt sqo outputs")


third.5.time<- alt.SQO_BLOE.generic(file_id="jan-26 test runs", 
                                  infauna_path="C:/Users/GisUSer/SCCWRP/Bight 2023 - Infauna/index code sub committee/B23 Index Coding/example SQO data/example SQO data.csv", 
                                  output_path="C:/Users/GisUser/SCCWRP/Bight 2023 - Infauna/index code sub committee/B23 Index Coding/alternate species list approach/alt sqo outputs")


ed14.rollups<-read.csv("alternate species list approach/SQO-ed14 roll up taxa.csv")
SoCal.SQO.ed14.link<-read.csv("alternate species list approach/Bight 23 SQO Appendix - djg.csv")
ed.14.complex<-read.csv("alternate species list approach/SQO-ed 14 complex change taxa.csv")  

save(file="Reference Files/ed14_to_SQO.RData", ed14.rollups,SoCal.SQO.ed14.link,ed.14.complex)
