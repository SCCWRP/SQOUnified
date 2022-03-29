library(tidyverse)
require(DBI) # needed to connect to database
require(dbplyr) # needed to connect to database
require(RPostgreSQL) # needed to connect to our database
require(rstudioapi) # just so we can type the password as we run the script, so it is not written in the clear
con <- DBI::dbConnect(
  PostgreSQL(),
  host = rstudioapi::showPrompt("username", "Please enter the hostname of IP address of the database"),
  dbname = rstudioapi::showPrompt("dbname", "Please enter the name of the database"),
  user = rstudioapi::showPrompt("username", "Please enter the username for the database"),
  rstudioapi::askForPassword()
)

chem <- dbGetQuery(con, "SELECT stationid, analytename, result, rl, mdl FROM tbl_chemresults")


prep_SQO_03_08_df <- function(chem){
  # Here chemdata consists of data in the same format as our database, with the columns
  # stationid, analytename, result, rl, mdl

  # result, rl, mdl should be numeric fields
  chem <- chem %>% mutate(
    result = as.numeric(result),
    rl = as.numeric(rl),
    mdl = as.numeric(mdl),
  )

  # lowercase column names
  names(chem) <- names(chem) %>% tolower()

  # Analytes that are not grouped in any particular category
  single_analytes <- c('Cadmium','Copper','Lead','Mercury','Zinc',
                       'alpha-Chlordane','gamma-Chlordane','trans-Nonachlor',"4,4'-DDT")

  # High PAH
  hpah <- c('Benz(a)anthracene', 'Benzo(a)pyrene', 'Benzo(e)pyrene',
                'Chrysene', 'Dibenz(a,h)anthracene', 'Fluoranthene', 'Perylene','Pyrene')

  # Low PAH
  lpah <- c('1-Methylnaphthalene', '1-Methylphenanthrene', '2,6-Dimethylnaphthalene',
                '2-Methylnaphthalene', 'Acenaphthene', 'Anthracene',
                'Biphenyl', 'Fluorene', 'Naphthalene', 'Phenanthrene')

  # The PCB's that we care about
  relevant_pcbs <- paste(
    'PCB', c(008, 018, 028, 044, 052, 066, 101, 105, 110, 118, 128, 138, 153, 180, 187, 195), sep = "-"
  )

  # The other group is the "DD's" DDT, DDE, DDD
  # Can be done with grepl

  # its hard to think of useful, good variable names. Cut me some slack.
  chemdata <- chem %>%
    # create a new column called compound. This is what we will group by,
    mutate(
      compound = case_when(
        analytename %in% single_analytes ~ analytename,
        analytename %in% relevant_pcbs ~ "PCBs_total",
        analytename %in% hpah ~ "HPAH",
        analytename %in% lpah ~ "LPAH",
        grepl("DDD",analytename) ~ "DDDs_total",
        grepl("DDE",analytename) ~ "DDEs_total",
        # The DDT total will be handled separately
        TRUE ~ NA_character_
      )
    ) %>%
    filter(
      !is.na(compound)
    ) %>%
    mutate(
      result = case_when(
        result == -88 & compound %in% single_analytes ~ as.numeric(1/2*mdl),
        result == -88 & !(compound %in% single_analytes) ~ 0,
        TRUE ~ result
      )
    ) %>%
    group_by(
      stationid, compound
    ) %>%
    summarize(
      # if the sum of the results is zero, assign it the max of the RL's
      result = if_else(
        sum(result, na.rm = T) != 0, sum(result, na.rm = T), max(rl)
      )
    ) %>%
    ungroup()


  ddts_total <- chem %>%
    filter(grepl("DDT",analytename)) %>%
    mutate(
      compound = "DDTs_total",
      result = if_else(
        result == -88, 0, result
      )
    ) %>%
    group_by(stationid,compound) %>%
    summarize(
      result = if_else(
        sum(result, na.rm = T) != 0, sum(result, na.rm = T), max(rl)
      )
    ) %>% ungroup()

  chemdata <- chemdata %>%
    bind_rows(ddts_total) %>%
    arrange(stationid, compound)


  # single_compound <- chemdata %>%
  #   filter(analytename %in% single_analytes) %>%
  #   mutate(
  #     result = case_when(
  #       result == -88 ~ as.double(1/2*mdl),
  #       TRUE ~ result
  #     )
  #   )






  # singles <- single_compound %>%
  #   select(stationid, analytename, result) %>%
  #   unique() %>%
  #   group_by(stationid, analytename) %>%
  #   mutate(grouped_id = row_number()) %>%
  #   spread(., key = analytename, value = result ) %>%
  #   select(-grouped_id)


  # data for total compounds ------------------------------------------------
  totals_match <- c('DDD', 'DDT', 'DDE', 'PCB')

  totals_compound <- df_select %>%
    filter(
      grepl(paste(totals_match, collapse = '|'), analytename)
    ) %>%
    separate(analytename, into = c('cat1', 'cat2'), sep = '-')


  # total DDDs --------------------------------------------------------------

  total_DDD <- totals_compound %>%
    filter(cat2 == 'DDD') %>%
    group_by(stationid) %>%
    mutate(result = case_when(
      result == -99 ~ 0,
      TRUE ~ result
    )
    ) %>%
    summarise(
      DDDs_total = sum(result),
      RL = max(rl)
    ) %>%
    mutate(
      DDDs_total = case_when(
        DDDs_total == 0L ~ RL,
        TRUE ~ DDDs_total
      ),
      RL = NULL
    )


  # total DDE ---------------------------------------------------------------

  total_DDE <- totals_compound %>%
    filter(cat2 == 'DDE') %>%
    group_by(stationid) %>%
    mutate(result = case_when(
      result == -99 ~ 0,
      TRUE ~ result
    )
    ) %>%
    summarise(
      DDEs_total = sum(result),
      RL = max(rl)
    ) %>%
    mutate(
      DDEs_total = case_when(
        DDEs_total == 0 ~ RL,
        TRUE ~ DDEs_total
      ),
      RL = NULL
    )


  # total DDT ---------------------------------------------------------------

  total_DDT <- totals_compound %>%
    filter(cat2 == 'DDT') %>%
    group_by(stationid) %>%
    mutate(result = case_when(
      result == -99 ~ 0,
      TRUE ~ result
    )
    ) %>%
    summarise(
      DDTs_total = sum(result),
      RL = max(rl)
    ) %>%
    mutate(
      DDTs_total = case_when(
        DDTs_total == 0 ~ RL,
        TRUE ~ DDTs_total
      ),
      RL = NULL
    )



  # data for PCB ------------------------------------------------------------

  # ACCORDING TO CASQO Manual page 30, looks like the sum of the results must be multiplied by
  # a correction factor of 1.72
  total_PCB <- totals_compound %>%
    filter(cat1 == 'PCB',
           cat2 %in% c(008, 018, 028, 044,
                       052, 066, 101, 105,
                       110, 118, 128, 138,
                       153, 180, 187, 195)
    ) %>%
    mutate(result = case_when(
      result == -99 ~ 0,
      TRUE ~ result
    )
    ) %>%
    group_by(stationid) %>%
    summarise(
      PCBs_total = sum(result),
      RL = max(rl)
    ) %>%
    mutate(
      PCBs_total = case_when(
        PCBs_total == 0L ~ RL,
        TRUE ~ PCBs_total
      ),
      RL = NULL
    )


  # HPAH --------------------------------------------------------------------

  hpah_cal <- c('Benz(a)anthracene', 'Benzo(a)pyrene', 'Benzo(e)pyrene',
                'Chrysene', 'Dibenz(a,h)anthracene', 'Fluoranthene', 'Perylene','Pyrene')

  HPAH_df <- df_select %>%
    select(stationid, analytename, result, rl) %>%
    filter(
      analytename %in% hpah_cal
    ) %>%
    mutate(result = case_when(
      result == -99 ~ 0,
      TRUE ~ result
    )
    ) %>%
    group_by(stationid) %>%
    summarise(
      HPAH = sum(result),
      RL = max(rl)
    ) %>%
    mutate(
      HPAH = case_when(
        HPAH == 0L ~ RL,
        TRUE ~ HPAH
      ),
      RL = NULL
    )


  # LPAH --------------------------------------------------------------------

  lpah_cal <- c('1-Methylnaphthalene', '1-Methylphenanthrene', '2,6-Dimethylnaphthalene',
                '2-Methylnaphthalene', 'Acenaphthene', 'Anthracene',
                'Biphenyl', 'Fluorene', 'Naphthalene', 'Phenanthrene')

  LPAH_df <- df_select %>%
    select(stationid, analytename, result, rl) %>%
    filter(
      analytename %in% lpah_cal
    ) %>%
    mutate(result = case_when(
      result == -99 ~ 0,
      TRUE ~ result
    )
    ) %>%
    group_by(stationid) %>%
    summarise(
      LPAH = sum(result),
      RL = max(rl)
    ) %>%
    mutate(
      LPAH = case_when(
        LPAH == 0L ~ RL,
        TRUE ~ LPAH
      ),
      RL = NULL
    )


  # total DD's --------------------------------------------------------------
  DDT_match <- c('DDT', 'DDD', 'DDE', 'DDMU')
  total_DDs <- df_select %>%
    filter(
      grepl(paste(DDT_match, collapse = '|'), analytename)
    ) %>%
    mutate(result = case_when(
      result == -99 ~ 0,
      TRUE ~ result
    )
    )%>%
    group_by(stationid) %>%
    summarise(
      total_DDs = sum(result)
    )

  # total PAH ---------------------------------------------------------------
  # total_PAH <- df_select %>%
  #   filter(
  #     analyteclass == 'PAH'
  #   ) %>%
  #   mutate(result = case_when(
  #     result == -99 ~ 0L,
  #     TRUE ~ result
  #   )
  #   )%>%
  #   group_by(stationid) %>%
  #   summarise(
  #     total_PAH = sum(result)
  #   )

  # PCB all -----------------------------------------------------------------
  # total_PCB_all <- df_select %>%
  #   filter(
  #     analyteclass == 'PCB'
  #   ) %>%
  #   mutate(result = case_when(
  #     result == -99 ~ 0,
  #     TRUE ~ result
  #   )
  #   )%>%
  #   group_by(stationid) %>%
  #   summarise(
  #     total_PCB_all = sum(result)
  #   )

  # total Phi ---------------------------------------------------------------
  # total_PHI <- df_select %>%
  #   filter(
  #     analytename %in% c(paste('Phi', seq(4.5, 20, by = 0.5)))
  #   ) %>%
  #   mutate(result = case_when(
  #     result == -99 ~ 0,
  #     TRUE ~ result
  #     )
  #   )%>%
  #   group_by(stationid) %>%
  #   summarise(
  #     total_Phi = sum(result)
  #   )

  # total Fipronils ---------------------------------------------------------
  # total_FIP <- df_select %>%
  #   filter(
  #     analyteclass == 'FIPRONIL'
  #   ) %>%
  #   mutate(result = case_when(
  #     result == -99 ~ 0,
  #     TRUE ~ result
  #   )
  #   ) %>%
  #   group_by(stationid) %>%
  #   summarise(
  #     total_Fipronil = sum(result)
  #   )

  # total PBDE --------------------------------------------------------------
  # total_PBDE <- df_select %>%
  #   filter(
  #     analyteclass == 'PBDE'
  #   ) %>%
  #   mutate(result = case_when(
  #     result == -99 ~ 0,
  #     TRUE ~ result
  #   )
  #   ) %>%
  #   group_by(stationid) %>%
  #   summarise(
  #     total_PBDE = sum(result)
  #   )

  # total Pyrethroid --------------------------------------------------------
  # total_PYR <- df_select %>%
  #   filter(
  #     analyteclass == 'Pyrethroid'
  #   ) %>%
  #   mutate(result = case_when(
  #     result == -99 ~ 0,
  #     TRUE ~ result
  #   )
  #   ) %>%
  #   group_by(stationid) %>%
  #   summarise(
  #     total_Pyrethroid = sum(result)
  #   )

  # total Chlordanes --------------------------------------------------------
  # total_CHLOR <- df_select %>%
  #   filter(
  #     !grepl(paste(DDT_match, collapse = '|'), analytename),
  #     analyteclass == 'Chlorinated Hydrocarbons'
  #   ) %>%
  #   mutate(result = case_when(
  #     result == -99 ~ 0,
  #     TRUE ~ result
  #   )
  #   ) %>%
  #   group_by(stationid) %>%
  #   summarise(
  #     total_Chloridanes = sum(result)
  #   )





  # combine all calculations ------------------------------------------------
  formated_df <- singles %>%
    inner_join(., total_DDD, by = 'stationid') %>%
    inner_join(., total_DDE, by = 'stationid') %>%
    inner_join(., total_DDT, by = 'stationid') %>%
    inner_join(., total_PCB, by = 'stationid') %>%
    rename(`4_4_DDT` = `4,4'-DDT`,
           Gamma_Chlordane = `gamma-Chlordane`,
           Alpha_Chlordane = `alpha-Chlordane`,
           # Trans_Nonachlor = `trans-Nonachlor`
    ) %>%
    left_join(., HPAH_df, by = 'stationid') %>%
    left_join(., LPAH_df, by = 'stationid') # %>%
    # left_join(., total_DDs, by = 'stationid') %>%
    # left_join(., total_PAH, by = 'stationid') %>%
    # left_join(., total_PCB_all, by = 'stationid') %>%
    # inner_join(., total_PHI, by = 'stationid') %>%
    # left_join(., total_FIP, by = 'stationid') %>%
    # left_join(., total_PBDE, by = 'stationid') %>%
    # left_join(., total_PYR, by = 'stationid') %>%
    # left_join(., total_CHLOR, by = 'stationid')

  return(formated_df)
}



