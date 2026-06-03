#### David's example script to run the SQO and BRI indices over different data sets #####

#BRI calculations

#following the instructions included in the function, i specify the four pieces of information
# (the things in quotations) needed to run the calculator

# this one runs the whole example data

full.bri<-offshore_bri_calc_ed14_v2(file_id="djg test",
                                 infauna_path="example BRI data/example infauna.csv",
                                 station_path="example BRI data/example stations.csv",
                                 output_path = "BRI dry run")

# this one runs the unchanged taxa example data

uct.bri<-offshore_bri_calc_ed14_v2(file_id="djg test - unchanged",
                                 infauna_path="example BRI data/example BRI infauna - unchanged taxa.csv",
                                 station_path="example BRI data/example BRI stations - unchanged taxa.csv",
                                 output_path = "BRI dry run")

#SQO Calculations

# this one runs the whole data example

full.sqo<-SQO_BLOE.generic_v2(file_id = "djg test",
                               infauna_path = "example SQO data/example SQO data.csv",
                               output_path = "SQO dry run")

# this one runs the unchanged taxa example data

uct.sqo<-SQO_BLOE.generic_v2(file_id = "djg test - unchanged",
                               infauna_path = "example SQO data/example SQO data - unchanged taxa.csv",
                               output_path = "SQO dry run")


