#' Compute the multivariate AMBI (M-AMBI) index score.
#'
#' @description
#'     This is a function to calculate multivariate AMBI (M-AMBI) index scores following Pelletier et al. 2018
#'     which is in turn built upon the work of Sigovini et al. 2013 and Muxica et al. 2007.
#'     This is an alternate version that allows for manipulation of the data within R before
#'     submitting it to the function in lieu of directly reading in an excel file to the function.
#'
#'     The function is designed for use in US estuarine waters and requires three arguments:
#'     BenthicData, EG_Ref_values, and EG_Scheme. More details are given below.
#'
#'     Two additional dataframes are also needed to run the script: Saline and Tidal freshwater
#'     good-bad standards for the M-AMBI that are in TidalFresh_Standards.RData and Saline_Standards.RData
#'
#'     For the function to run, the following packages NEED to be installed:  tidyverse, reshape2, vegan,
#'     and readxl. Additionally the EQR.R function must also be installed and is included with this code.
#'
#' @param \strong{BenthicData} a data frame with the following columns:
#'
#'    \strong{\code{StationID}} - an alpha-numeric identifier of the location;
#'
#'    \strong{\code{Replicate}} - a numeric identifying the replicate number of samples taken at the location;
#'
#'    \strong{\code{SampleDate}} - the date of sample collection;
#'
#'    \strong{\code{Latitude}} - latitude in decimal degrees;
#'
#'    \strong{\code{Longitude}} - longitude in decimal degrees. Make sure there is a negative sign for the Western coordinates;
#'
#'    \strong{\code{Taxon}} - name of the fauna, ideally in SCAMIT ed12 format, do not use sp. or spp.,
#'        use sp only or just the Genus. If no animals were present in the sample use
#'        NoOrganismsPresent with 0 abundance;
#'
#'    \strong{\code{Abundance}} - the number of each Species observed in a sample;
#'
#'    \strong{\code{Salinity}} - the salinity observed at the location in PSU, ideally at time of sampling.
#'
#' @param
#'     \strong{EG_Ref_values} -  A data frame with the suite of US Ecological Groups assigned
#'     initially in Gillett et al. 2015. This EG Ref values has multiple versions of the EG values and a Yes/No designation
#'     if the fauna are Oligochaetes or not. The default dataframe is one called EG_Ref which was originally
#'     read in from a csv called "Ref - EG Values 2018.csv."
#'
#'     Replace with other data as you see fit, but make sure the data you use is in a similar format and uses
#'     the same column names. Additionally, new taxa can be added at the bottom of the list with the EG values the user
#'     feels appropriate, THOUGH THIS IS NOT RECOMMENDED
#' @param
#'     \strong{EG_Scheme} A quoted string with the name of the EG Scheme to be used in the AMBI scoring.
#'     The default is Hybrid, though one could use US (all coasts),
#'     Standard (Values from Angel Borja and colleagues),
#'     US_East (US East Coast), US_Gulf (US Gulf of Mexico Coast), or US_West (US West Coast).
#'
#' @usage MAMBI(benthic_data, EG_Ref_values = EG_Ref, EG_Scheme = "Hybrid")
#'
#' @examples
#' data(benthic_sampledata) # load sample dataset to environment
#' data(TidalFresh_Standards) # Load tidal fresh standards if you want to look at them
#' data(Saline_Standards) # Load Saline standards if you want to look at them
#' data(EG_Ref) # load the default EG_Ref values
#' MAMBI(benthic_sampledata, EG_Ref, "Hybrid")
#' MAMBI(benthic_sampledata, EG_Ref, "US Gulf")
#' MAMBI(benthic_sampledata) # uses default values for the last two args
#'
#' @return
#'     The output of the function will be a dataframe with
#'
#'     \code{StationID}, \code{Replicate}, \code{SampleDate},
#'
#'     \code{Latitude}, \code{Longitude},
#'
#'     \code{SalZone} (The Salinity Zone assigned by M-AMBI),
#'
#'     \code{AMBI_Score},
#'
#'     \code{S} (Species Richness),
#'
#'     \code{H} (Species Diversity),
#'
#'     \code{Oligo_pct} (Relative Abundance of Oligochaetes),
#'
#'     \code{MAMBI_Score}, \code{Orig_MAMBI_Condition}, \code{New_MAMBI_Condition},
#'
#'     \code{Use_MAMBI} (Can M-AMBI be applied?), \code{Use_AMBI} (Can AMBI be applied?),
#'
#'     \code{YesEG} (percentage of Abundance with a EG value)
#'
#' @author David Gillett \email{davidg@@sccwrp.org}
#'
#' @import stringr
#' @import reshape2
#' @import vegan
#' @importFrom dplyr mutate filter select group_by summarise left_join bind_rows bind_cols
#' @importFrom tidyr replace_na
#' @importFrom purrr map_dfr
#' @export

MAMBI<-function(BenthicData, EG_Ref_values = NULL, EG_Scheme="Hybrid", logfile = file.path(getwd(), 'logs', paste0(format(Sys.time(), "%Y-%m-%d_%H:%M:%S"), '-log.txt') ), verbose = T)
{

  if (is.null(EG_Ref_values)) {
    EG_Ref_values <- EG_Ref
  }

  #Saline_Standards <- saline_standards
  #TidalFresh_Standards <- tidalFresh_Standards
  #Saline_Standards<-read_xlsx("data/Pelletier2018_Standards.xlsx", sheet = "Saline Sites")# Good-Bad Benchmarks following Pelletier et al. 2018
  #TidalFresh_Standards<-read_xlsx("data/Pelletier2018_Standards.xlsx", sheet = "Tidal Fresh Sites") #Good-Bad Benchmarks following Pelletier et al. 2018

  Input_File.0 <- BenthicData %>%
    mutate(
      Species_ended_in_sp = (str_detect(Taxon," sp$")),
      Taxon=(str_replace(Taxon, " sp$",""))
    ) %>%
    mutate(
      Coast = (ifelse(Longitude<=-115,"West","Gulf-East"))
    )%>%
    mutate(
      SalZone = case_when(
        Salinity > 30 & Salinity <= 40 & Coast == "Gulf-East" ~ "EH",
        Salinity > 18 & Salinity <= 30 & Coast == "Gulf-East" ~ "PH",
        Salinity > 5 & Salinity <= 18 ~ "MH",
        Salinity > 0.2 & Salinity <= 5 ~ "OH",
        Salinity >= 0 & Salinity <= 0.2 ~ "TF",
        Salinity > 40 ~ "HH",
        Salinity > 30 & Salinity <= 40 & Coast == "West" ~ "WEH",
        Salinity > 18 & Salinity <= 30 & Coast == "West" ~ "WPH"
      )
    )



  EG_Ref_values <- EG_Ref_values %>%
    select(Taxon, Exclude, EG=EG_Scheme) %>%
    mutate(
      EG = ifelse(
        Taxon=="Oligochaeta", "V", EG
      )
    )
  #EG_Ref_values <- EG_Ref_values %>% select(.,Taxon, Exclude, EG=EG_Scheme) %>% mutate(EG=(ifelse(Taxon=="Oligochaeta", "V", EG)))

  #TODO: Need to fix these lines!!!
  azoic.samples<-Input_File.0 %>%
    filter(Taxon=="No Organisms Present") %>%
    select(StationID, Replicate, SampleDate, Latitude, Longitude, SalZone, Stratum)


  #azoic.samples <- if(dim(azoic.samples)[1] == 0){StationID = 1 & Replicate = 1}
  #  dplyr::mutate(StationID = case_when(dim(StationID) == dim(NA) ~ NA))
  #  dplyr::mutate(if(dim(azoic.samples)[1] == 0){AMBI_Score = 7  & S=0 & H=0 & Oligo_pct=0 & MAMBI_Score=0 & Orig_MAMBI_Condition="Bad" & New_MAMBI_Condition="High Disturbance" & Use_MAMBI="Yes" & Use_AMBI="Yes - Azoic" & YesEG=NA})

  Input_File<-Input_File.0 %>% filter(Taxon != "No Organisms Present")


  total.abundance<-Input_File %>%
    group_by(StationID, Replicate, SampleDate) %>%
    summarise(Tot_abun=sum(Abundance))

  Sample.info<-Input_File %>%
    select(StationID, Replicate, SampleDate, Latitude, Longitude, Salinity, Coast, SalZone, Stratum) %>%
    unique()

  Input_File2<-Input_File %>%
    filter(!is.na(SalZone))

  EG.Assignment<-Input_File %>%
    left_join(., EG_Ref_values, by="Taxon") %>% #filter(Exclude!="Yes") #%>%
    left_join(.,total.abundance, by=c("StationID", "Replicate", "SampleDate")) %>%
    mutate(
      Rel_abun=((Abundance/Tot_abun)*100)
    )

  EG.Assignment.cast<-data.frame(NoEG=numeric(),
                                 YesEG=numeric())

  AMBI.applicability<-EG.Assignment %>%
    mutate(EG_Test=ifelse(is.na(EG),"NoEG", "YesEG")) %>%
    reshape2::dcast(.,StationID+Replicate+SampleDate~EG_Test, value.var = "Rel_abun", fun.aggregate = sum) %>%
    left_join(.,EG.Assignment.cast) %>%
    mutate(
      Use_AMBI = case_when(
        NoEG <= 20 ~ "Yes",
        NoEG > 20 & NoEG <= 50 ~ "With Care",
        NoEG > 50 ~ "Not Recommended",
        is.na(NoEG) ~ "Yes"
      )
    )

  MAMBI.applicability <- Sample.info %>%
    mutate(
      Use_MAMBI = ifelse(is.na(SalZone),"No - No Salinity Value","Yes")
    ) %>%
    select(StationID, Replicate, SampleDate, Stratum, Use_MAMBI)

  Sal_range.dataset<-unique(Input_File2$SalZone)


  ######Saline calcs ################

  AMBI.Scores<-EG.Assignment %>%
    group_by(StationID, Replicate, SampleDate,Tot_abun,EG) %>%
    summarise(Sum_Rel=sum(Rel_abun)) %>%
    replace_na(list(EG="NoEG")) %>%
    mutate(
      EG_Score = case_when(
        EG == "I" ~ Sum_Rel*0,
        EG == "II" ~ Sum_Rel*1.5,
        EG == "III" ~ Sum_Rel*3,
        EG == "IV" ~ Sum_Rel*4.5,
        EG == "V" ~ Sum_Rel*6,
        EG == "NoEG" ~ 0
      )
    ) %>%
    mutate(EG_Score=ifelse(Tot_abun==0,700,EG_Score)) %>%
    group_by(StationID, Replicate, SampleDate) %>%
    summarise(AMBI_Score=(sum(EG_Score, na.rm=TRUE)/100))

  Rich<-Input_File %>%
    group_by(StationID, Replicate, SampleDate) %>%
    summarise(S=length(Taxon))

  Rich$S<-as.numeric(Rich$S)

  Divy<-Input_File %>%
    reshape2::dcast(StationID+Replicate+SampleDate~Taxon, value.var = "Abundance", fill=0) %>%
    mutate(H=vegan::diversity((select(.,4:(ncol(.)))),index = "shannon", base = 2)) %>%
    select(.,StationID, Replicate, SampleDate, H)



  metrics<-AMBI.Scores %>%
    left_join(.,Rich, by=c("StationID", "Replicate", "SampleDate")) %>%
    left_join(.,Divy, by=c("StationID", "Replicate", "SampleDate")) %>%
    mutate(
      S = (ifelse(AMBI_Score==7,0,S)),H=(ifelse(AMBI_Score==7,0,H))
    )

  metrics.1 <- Sample.info %>%
    left_join(., metrics, by=c("StationID", "Replicate", "SampleDate")) %>%
    select(StationID, Replicate, SampleDate,AMBI_Score, S, H, SalZone)

  no.SalZone.data<-filter(metrics.1, is.na(SalZone)) %>%
    left_join(.,Sample.info, by=c("StationID", "Replicate", "SalZone", "SampleDate")) %>%
    select(1:9)

  metrics.2<- metrics.1 %>%
    filter(!is.na(SalZone)) %>%
    bind_rows(.,Saline_Standards)


  saline.mambi<-map_dfr(Sal_range.dataset, function(sal){
    sal.df<- filter(metrics.2, SalZone == sal)
    METRICS.tot<-sal.df[,c(4:6)]


    options(warn = -1)
    METRICS.fa2 <- princomp(METRICS.tot, cor = T, covmat = cov(METRICS.tot))
    options(warn = 0)
    METRICS.fa2.load <- loadings(METRICS.fa2) %*% diag(METRICS.fa2$sdev)
    METRICS.fa2.load.varimax <- loadings(varimax(METRICS.fa2.load))
    METRICS.scores2 <- scale(METRICS.tot) %*% METRICS.fa2.load.varimax
    colnames(METRICS.scores2) <- c("x", "y", "z")
    METRICS.tr <- METRICS.scores2



    eqr <-EQR(METRICS.tr)
    colnames(eqr)<-c("MAMBI_Score")
    eqr<-data.frame(eqr)

    results<-sal.df %>%
      bind_cols(.,eqr) %>%
      left_join(.,Sample.info, by=c("StationID", "Replicate", "SampleDate", "SalZone")) %>%
      select(1,2,3,9,10,7,4:6,8) %>%
      filter(!StationID%in%Saline_Standards$StationID, SalZone!="TF") %>%
      mutate(
        Orig_MAMBI_Condition = case_when(
          MAMBI_Score < 0.2 ~ "Bad",
          MAMBI_Score >= 0.2 & MAMBI_Score < 0.39 ~ "Poor",
          MAMBI_Score >= 0.39 & MAMBI_Score < 0.53 ~ "Moderate",
          MAMBI_Score >= 0.53 & MAMBI_Score < 0.77 ~ "Good",
          MAMBI_Score >= 0.77 ~ "High"
        ),
        New_MAMBI_Condition = case_when(
          MAMBI_Score <= 0.387 ~ "High Disturbance",
          MAMBI_Score > 0.387 & MAMBI_Score < 0.483 ~ "Moderate Disturbance",
          MAMBI_Score >= 0.483 & MAMBI_Score < 0.578 ~ "Low Disturbance",
          MAMBI_Score >= 0.578 ~ "Reference")
        )
    saline.mambi<-results
  })

  ###################

  if(any(Sal_range.dataset=="TF"))
  {

    TF.EG.Assignment <- EG.Assignment %>% filter(SalZone=="TF")
    TF.EG_Ref_values <- EG_Ref %>% select(.,Taxon, Exclude, EG=EG_Scheme, Oligochaeta)

    TF.AMBI.Scores <- TF.EG.Assignment %>%
      group_by(StationID, Replicate, SampleDate,Tot_abun,EG) %>%
      summarise(Sum_Rel=sum(Rel_abun)) %>%
      replace_na(list(EG="NoEG")) %>%
      mutate(
        EG_Score = case_when(
          EG == "I" ~ Sum_Rel*0,
          EG == "II" ~ Sum_Rel*1.5,
          EG == "III" ~ Sum_Rel*3,
          EG == "IV" ~ Sum_Rel*4.5,
          EG == "V" ~ Sum_Rel*6,
          EG == "NoEG" ~ 0)
        ) %>%
      mutate(
        EG_Score = ifelse(Tot_abun == 0,700,EG_Score)
      ) %>%
      group_by(StationID, Replicate, SampleDate) %>%
      summarise(
        AMBI_Score = sum(EG_Score)/100
      )

    TF.Oligos <- Input_File %>%
      left_join(., total.abundance, by =c("StationID", "Replicate", "SampleDate")) %>%
      left_join(., TF.EG_Ref_values, by="Taxon" ) %>%
      filter(Oligochaeta=="Yes", SalZone=="TF") %>%
      group_by(StationID, Replicate, SampleDate) %>%
      summarise(
        Oligo_pct = sum(Abundance/Tot_abun) * 100
      )

    TF.Divy <- Input_File %>%
      filter(SalZone=="TF") %>%
      reshape2::dcast(StationID+Replicate+SampleDate~Taxon, value.var = "Abundance", fill=0) %>%
      mutate(
        H = diversity((select(.,4:(ncol(.)))), index = "shannon", base = 2)
      ) %>%
      select(.,StationID, Replicate, SampleDate,H)


    TF.metrics <- TF.AMBI.Scores %>%
      left_join(.,TF.Divy, by=c("StationID", "Replicate", "SampleDate")) %>%
      left_join(.,TF.Oligos, by=c("StationID", "Replicate", "SampleDate")) %>%
      mutate(
        Oligo_pct = (ifelse(AMBI_Score==7,0,Oligo_pct)),
        H=(ifelse(AMBI_Score==7,0,H)),
        Oligo_pct=(ifelse(is.na(Oligo_pct),0,Oligo_pct))
      )

    TF.metrics.1<-Sample.info %>%
      left_join(., TF.metrics, by=c("StationID", "Replicate", "SampleDate")) %>%
      select(StationID, Replicate, SampleDate,AMBI_Score, H, Oligo_pct, SalZone) %>%
      filter(SalZone=="TF")

    TF.metrics.2<-bind_rows(TF.metrics.1,TidalFresh_Standards)

    TF.METRICS.tot<-TF.metrics.2[,c(4:6)]

    options(warn = -1)
    TF.METRICS.fa2 <- princomp (TF.METRICS.tot, cor = T, covmat = cov(TF.METRICS.tot))
    options(warn = 0)
    TF.METRICS.fa2.load <- loadings(TF.METRICS.fa2) %*% diag(TF.METRICS.fa2$sdev)
    TF.METRICS.fa2.load.varimax <- loadings(varimax(TF.METRICS.fa2.load))
    TF.METRICS.scores2 <- scale(TF.METRICS.tot) %*% TF.METRICS.fa2.load.varimax
    colnames(TF.METRICS.scores2) <- c("x", "y", "z")
    TF.METRICS.tr <- TF.METRICS.scores2

    TF.eqr <-EQR(TF.METRICS.tr)
    colnames(TF.eqr) <- c("MAMBI_Score")
    TF.eqr <- data.frame(TF.eqr)

    TF.mambi <- TF.metrics.2 %>%
      bind_cols(.,TF.eqr) %>%
      left_join(.,Sample.info, by=c("StationID", "Replicate", "SalZone", "SampleDate")) %>%
      select(1,2,3,9,10,7,4:6,8) %>%
      filter(!StationID%in%TidalFresh_Standards$StationID) %>%
      mutate(
        Orig_MAMBI_Condition = case_when(
          MAMBI_Score < 0.2 ~ "Bad",
          MAMBI_Score >= 0.2 & MAMBI_Score < 0.39 ~ "Poor",
          MAMBI_Score >= 0.39 & MAMBI_Score < 0.53 ~ "Moderate",
          MAMBI_Score >= 0.53 & MAMBI_Score < 0.77 ~ "Good",
          MAMBI_Score >= 0.77 ~ "High"
        ),
        New_MAMBI_Condition = case_when(
          MAMBI_Score <= 0.387 ~ "High Disturbance",
          MAMBI_Score > 0.387 & MAMBI_Score < 0.483 ~ "Moderate Disturbance",
          MAMBI_Score >= 0.483 & MAMBI_Score < 0.578 ~ "Low Disturbance",
          MAMBI_Score >= 0.578 ~ "Reference")
        )

    saline.mambi<- saline.mambi %>% mutate(Oligo_pct=NA) %>% select(1:9,13,10,11,12)
    TF.mambi<-TF.mambi %>% mutate(S=NA) %>% select(1:7,13,8:12)

    Overall.Results <- bind_rows(saline.mambi, TF.mambi) %>%
      bind_rows(.,no.SalZone.data) %>%
      left_join(.,AMBI.applicability[,c(1,2,3,5,6)], by=c("StationID", "Replicate", "SampleDate")) %>%
      left_join(MAMBI.applicability,.,by=c("StationID", "Replicate", "SampleDate")) %>%
      #rename(B13_Stratum = Stratum) %>%
      select(1:3,5:14,4,16,15) %>%
      bind_rows(.,azoic.samples) %>%
      mutate(Index = "M-AMBI")
  }


  else
  {
    Overall.Results<-saline.mambi %>%
      bind_rows(.,no.SalZone.data) %>%
      mutate(Oligo_pct=NA) %>%
      left_join(.,AMBI.applicability[,c(1,2,3,5,6)], by=c("StationID", "Replicate", "SampleDate")) %>%
      left_join(MAMBI.applicability,.,by=c("StationID", "Replicate", "SampleDate")) %>%
      #select(1:3,5:10,14,11:13,4,16,15) %>%
      select(StationID, SampleDate, Replicate, Stratum, Use_MAMBI, Latitude, Longitude, SalZone, AMBI_Score, S, H, MAMBI_Score, Orig_MAMBI_Condition, New_MAMBI_Condition, Oligo_pct, YesEG, Use_AMBI) %>%
      #rename(B13_Stratum = Stratum) %>%
      bind_rows(.,azoic.samples) %>%
      mutate(Index = "M-AMBI")
  }

  return(Overall.Results)

}

# Can run this to test
#test = MAMBI(benthic_data, EG_File_Name="data/Ref - EG Values 2018.csv", EG_Scheme="Hybrid")
#write.csv(test, file = "data/output-overall-results.csv", row.names = FALSE)

