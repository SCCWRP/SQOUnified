# Script to generate updated sample data for SQOUnified package
# Creates: chem_sampledata, tox_sampledata, benthic_sampledata, offshore_bri_sampledata
# (the offshore infauna and station info are combined into a single offshore_bri_sampledata frame)

library(dplyr)

# ============================================================================
# 1. CHEMISTRY SAMPLE DATA
# ============================================================================
# Stations represent a range of contamination levels:
#   "0000"     = control (all non-detects)
#   "S24-1001" = clean/reference site
#   "S24-1002" = low contamination
#   "S24-1003" = moderate contamination
#   "S24-1004" = high contamination
#   "S24-1005" = mixed contamination

# Bay/estuary stations (shared across chem, tox, benthic)
# Offshore stations (shared across chem, tox, offshore benthic/stations — NOT in benthic)
chem_stations <- c("0000", "S24-1001", "S24-1002", "S24-1003", "S24-1004", "S24-1005",
                   "OFS-01", "OFS-02", "OFS-03", "OFS-04")

# -- Metals (mg/dry kg) --
metals <- c("Aluminum", "Antimony", "Arsenic", "Barium", "Beryllium",
            "Cadmium", "Chromium", "Copper", "Iron", "Lead",
            "Mercury", "Nickel", "Selenium", "Silver", "Zinc")

metal_rl  <- c(5.0, 0.05, 0.05, 0.05, 0.05,
               0.005, 0.005, 0.005, 5.0, 0.005,
               0.03, 0.02, 0.05, 0.02, 0.05)
metal_mdl <- c(1.0, 0.025, 0.025, 0.025, 0.025,
               0.0025, 0.0025, 0.0025, 1.0, 0.0025,
               0.003, 0.01, 0.025, 0.01, 0.025)

# Values per station (rows = analytes in order of metals vector)
# 0000 = control, all non-detect (-88)
metal_vals <- list(
  "0000"     = rep(-88, length(metals)),
  "S24-1001" = c(12000, 0.3, 5.2, 40.1, 0.5,  0.08, 30.5, 15.0, 18000, 8.0,   0.04, 10.2, 0.2, 0.03, 55.0),
  "S24-1002" = c(15000, 0.8, 7.5, 55.3, 0.7,  0.15, 42.0, 55.0, 22000, 28.0,  0.10, 15.8, 0.4, 0.06, 120.0),
  "S24-1003" = c(20000, 1.5, 12.0, 80.2, 1.2, 0.50, 65.0, 110.0, 28000, 70.0, 0.50, 25.3, 0.8, 0.15, 220.0),
  "S24-1004" = c(25000, 3.0, 18.5, 120.0, 2.0, 1.80, 95.0, 450.0, 35000, 180.0, 2.50, 40.0, 1.5, 0.50, 700.0),
  "S24-1005" = c(18000, 1.0, 9.0, 60.0, 0.9,  0.30, 50.0, 80.0, 25000, 45.0, 0.25, 18.0, 0.5, 0.08, 160.0),
  "OFS-01"   = c(11000, 0.2, 4.8, 35.0, 0.4,  0.06, 28.0, 12.0, 16000, 6.5,  0.03, 9.0, 0.15, 0.02, 48.0),
  "OFS-02"   = c(14000, 0.5, 6.8, 48.0, 0.6,  0.12, 38.0, 42.0, 20000, 22.0, 0.08, 13.5, 0.3, 0.05, 95.0),
  "OFS-03"   = c(16000, 0.9, 8.5, 62.0, 0.8,  0.20, 48.0, 65.0, 24000, 35.0, 0.15, 17.0, 0.45, 0.07, 140.0),
  "OFS-04"   = c(22000, 2.0, 15.0, 95.0, 1.5, 1.00, 78.0, 200.0, 30000, 100.0, 1.00, 30.0, 1.0, 0.25, 350.0)
)

chem_metals <- do.call(rbind, lapply(chem_stations, function(stn) {
  data.frame(
    stationid      = stn,
    analytename    = metals,
    result         = as.character(metal_vals[[stn]]),
    rl             = metal_rl,
    mdl            = metal_mdl,
    units          = "mg/kg dw",
    fieldrep       = 1,
    labrep         = 1,
    sampletypecode = "Result",
    stringsAsFactors = FALSE
  )
}))

# -- PAHs (ug/dry kg) --
lpahs <- c("Acenaphthene", "Acenaphthylene", "Anthracene", "Biphenyl",
            "Fluorene", "Naphthalene", "Phenanthrene",
            "1-Methylnaphthalene", "2-Methylnaphthalene",
            "2,6-Dimethylnaphthalene", "1,6,7-Trimethylnaphthalene",
            "1-Methylphenanthrene")
hpahs <- c("Benz(a)anthracene", "Benzo(a)pyrene", "Benzo(b)fluoranthene",
            "Benzo(e)pyrene", "Benzo(g,h,i)perylene", "Benzo(k)fluoranthene",
            "Chrysene", "Dibenz(a,h)anthracene", "Fluoranthene",
            "Indeno(1,2,3-cd)pyrene", "Perylene", "Pyrene")

pah_rl <- 0.5
pah_mdl <- 0.1

# LPAH values per station
lpah_vals <- list(
  "0000"     = rep(-88, length(lpahs)),
  "S24-1001" = c(2.0, 1.5, 3.0, 1.0, 2.5, 3.5, 8.0, 1.5, 2.0, 0.8, 0.5, 1.2),
  "S24-1002" = c(5.0, 3.5, 8.0, 2.5, 6.0, 8.5, 22.0, 4.0, 5.5, 2.0, 1.5, 3.0),
  "S24-1003" = c(15.0, 10.0, 25.0, 8.0, 18.0, 22.0, 65.0, 12.0, 16.0, 6.0, 4.5, 9.0),
  "S24-1004" = c(80.0, 55.0, 130.0, 40.0, 90.0, 120.0, 350.0, 60.0, 85.0, 30.0, 22.0, 45.0),
  "S24-1005" = c(10.0, 7.0, 15.0, 5.0, 12.0, 16.0, 45.0, 8.0, 11.0, 4.0, 3.0, 6.0),
  "OFS-01"   = c(1.5, 1.0, 2.0, 0.8, 2.0, 2.5, 6.0, 1.0, 1.5, 0.6, 0.4, 0.8),
  "OFS-02"   = c(4.0, 3.0, 6.5, 2.0, 5.0, 7.0, 18.0, 3.5, 4.5, 1.8, 1.2, 2.5),
  "OFS-03"   = c(8.0, 5.5, 12.0, 4.0, 9.0, 13.0, 38.0, 6.5, 9.0, 3.5, 2.5, 5.0),
  "OFS-04"   = c(40.0, 28.0, 65.0, 20.0, 45.0, 60.0, 180.0, 30.0, 42.0, 15.0, 11.0, 22.0)
)

# HPAH values per station
hpah_vals <- list(
  "0000"     = rep(-88, length(hpahs)),
  "S24-1001" = c(5.0, 4.0, 5.5, 3.5, 3.0, 2.5, 6.0, 1.0, 10.0, 2.5, 2.0, 8.0),
  "S24-1002" = c(18.0, 15.0, 20.0, 12.0, 10.0, 8.0, 22.0, 3.5, 35.0, 8.5, 6.0, 30.0),
  "S24-1003" = c(55.0, 48.0, 60.0, 38.0, 32.0, 25.0, 68.0, 10.0, 110.0, 28.0, 18.0, 95.0),
  "S24-1004" = c(280.0, 250.0, 320.0, 200.0, 170.0, 130.0, 350.0, 55.0, 580.0, 145.0, 90.0, 490.0),
  "S24-1005" = c(38.0, 32.0, 42.0, 26.0, 22.0, 17.0, 48.0, 7.0, 75.0, 19.0, 12.0, 65.0),
  "OFS-01"   = c(3.5, 3.0, 4.0, 2.5, 2.0, 1.8, 4.5, 0.7, 7.0, 1.8, 1.5, 6.0),
  "OFS-02"   = c(14.0, 12.0, 16.0, 10.0, 8.0, 6.5, 18.0, 2.8, 28.0, 7.0, 5.0, 24.0),
  "OFS-03"   = c(30.0, 26.0, 34.0, 22.0, 18.0, 14.0, 38.0, 6.0, 60.0, 15.0, 10.0, 52.0),
  "OFS-04"   = c(140.0, 125.0, 160.0, 100.0, 85.0, 65.0, 175.0, 28.0, 290.0, 72.0, 45.0, 245.0)
)

chem_lpahs <- do.call(rbind, lapply(chem_stations, function(stn) {
  data.frame(
    stationid      = stn,
    analytename    = lpahs,
    result         = as.character(lpah_vals[[stn]]),
    rl             = pah_rl,
    mdl            = pah_mdl,
    units          = "ug/kg dw",
    fieldrep       = 1,
    labrep         = 1,
    sampletypecode = "Result",
    stringsAsFactors = FALSE
  )
}))

chem_hpahs <- do.call(rbind, lapply(chem_stations, function(stn) {
  data.frame(
    stationid      = stn,
    analytename    = hpahs,
    result         = as.character(hpah_vals[[stn]]),
    rl             = pah_rl,
    mdl            = pah_mdl,
    units          = "ug/kg dw",
    fieldrep       = 1,
    labrep         = 1,
    sampletypecode = "Result",
    stringsAsFactors = FALSE
  )
}))

# -- Pesticides (ug/dry kg) --
pesticides <- c("alpha-Chlordane", "gamma-Chlordane", "trans-Nonachlor",
                "cis-Nonachlor", "Oxychlordane",
                "4,4'-DDD", "4,4'-DDE", "4,4'-DDT", "4,4'-DDMU",
                "2,4'-DDD", "2,4'-DDE", "2,4'-DDT")

pest_rl <- 0.5
pest_mdl <- 0.1

pest_vals <- list(
  "0000"     = rep(-88, length(pesticides)),
  "S24-1001" = c(0.1, 0.1, -88, -88, -88, 0.2, 0.3, -88, -88, -88, -88, -88),
  "S24-1002" = c(0.5, 0.6, 0.2, 0.1, -88, 0.8, 1.2, 0.3, 0.2, 0.1, -88, -88),
  "S24-1003" = c(1.5, 1.8, 0.8, 0.4, 0.2, 3.0, 5.5, 1.5, 0.8, 0.5, 0.2, 0.3),
  "S24-1004" = c(12.0, 15.0, 5.0, 2.5, 1.5, 28.0, 50.0, 8.0, 5.0, 3.0, 1.5, 2.0),
  "S24-1005" = c(0.8, 1.0, 0.4, 0.2, -88, 2.0, 3.5, 0.8, 0.5, 0.3, -88, 0.1),
  "OFS-01"   = rep(-88, length(pesticides)),
  "OFS-02"   = c(0.3, 0.4, -88, -88, -88, 0.5, 0.8, 0.2, -88, -88, -88, -88),
  "OFS-03"   = c(0.8, 1.0, 0.3, 0.2, -88, 1.5, 2.8, 0.6, 0.3, 0.2, -88, -88),
  "OFS-04"   = c(6.0, 7.5, 2.5, 1.2, 0.8, 14.0, 25.0, 4.0, 2.5, 1.5, 0.8, 1.0)
)

chem_pesticides <- do.call(rbind, lapply(chem_stations, function(stn) {
  data.frame(
    stationid      = stn,
    analytename    = pesticides,
    result         = as.character(pest_vals[[stn]]),
    rl             = pest_rl,
    mdl            = pest_mdl,
    units          = "ug/kg dw",
    fieldrep       = 1,
    labrep         = 1,
    sampletypecode = "Result",
    stringsAsFactors = FALSE
  )
}))

# -- PCBs (ug/dry kg) --
pcbs <- c("PCB-008", "PCB-018", "PCB-028", "PCB-037", "PCB-044",
          "PCB-049", "PCB-052", "PCB-066", "PCB-070", "PCB-074",
          "PCB-077", "PCB-081", "PCB-087", "PCB-099", "PCB-101",
          "PCB-105", "PCB-110", "PCB-114", "PCB-118", "PCB-119",
          "PCB-123", "PCB-126", "PCB-128", "PCB-138", "PCB-149",
          "PCB-151", "PCB-153", "PCB-156", "PCB-157", "PCB-158",
          "PCB-167", "PCB-168", "PCB-169", "PCB-170", "PCB-177",
          "PCB-180", "PCB-183", "PCB-187", "PCB-189", "PCB-194",
          "PCB-195", "PCB-201", "PCB-206")

pcb_rl <- 0.5
pcb_mdl <- 0.1

# Generate PCB values - mostly non-detect for clean, increasing for contaminated
make_pcb_vals <- function(total_target) {
  # Distribute across congeners with realistic pattern
  n <- length(pcbs)
  if (total_target <= 0) return(rep(-88, n))
  # Approximate distribution - heavier congeners get more
  wts <- c(rep(0.01, 10), rep(0.02, 10), rep(0.03, 10), rep(0.04, 8), rep(0.03, 5))
  wts <- wts / sum(wts)
  vals <- round(total_target * wts, 2)
  # Some below detection
  vals[vals < 0.1] <- -88
  vals
}

pcb_vals <- list(
  "0000"     = rep(-88, length(pcbs)),
  "S24-1001" = make_pcb_vals(5),
  "S24-1002" = make_pcb_vals(15),
  "S24-1003" = make_pcb_vals(30),
  "S24-1004" = make_pcb_vals(300),
  "S24-1005" = make_pcb_vals(20),
  "OFS-01"   = make_pcb_vals(3),
  "OFS-02"   = make_pcb_vals(10),
  "OFS-03"   = make_pcb_vals(18),
  "OFS-04"   = make_pcb_vals(150)
)

chem_pcbs <- do.call(rbind, lapply(chem_stations, function(stn) {
  data.frame(
    stationid      = stn,
    analytename    = pcbs,
    result         = as.character(pcb_vals[[stn]]),
    rl             = pcb_rl,
    mdl            = pcb_mdl,
    units          = "ug/kg dw",
    fieldrep       = 1,
    labrep         = 1,
    sampletypecode = "Result",
    stringsAsFactors = FALSE
  )
}))

# -- Pyrethroids (ug/dry kg) --
pyrethroids <- c("Bifenthrin", "Cyfluthrin", "Cyhalothrin, lambda",
                 "Cypermethrin", "Deltamethrin", "Esfenvalerate",
                 "Permethrin, cis", "Permethrin, trans")

pyr_rl <- 0.5
pyr_mdl <- 0.1

pyr_vals <- list(
  "0000"     = rep(-88, length(pyrethroids)),
  "S24-1001" = rep(-88, length(pyrethroids)),
  "S24-1002" = c(0.5, -88, -88, -88, -88, -88, 0.3, 0.2),
  "S24-1003" = c(2.0, 0.5, 0.3, 0.8, -88, -88, 1.5, 1.0),
  "S24-1004" = c(8.0, 2.5, 1.5, 3.0, 1.0, 0.5, 6.0, 4.0),
  "S24-1005" = c(1.0, -88, -88, 0.3, -88, -88, 0.8, 0.5),
  "OFS-01"   = rep(-88, length(pyrethroids)),
  "OFS-02"   = rep(-88, length(pyrethroids)),
  "OFS-03"   = c(0.3, -88, -88, -88, -88, -88, 0.2, -88),
  "OFS-04"   = c(4.0, 1.2, 0.8, 1.5, 0.5, 0.3, 3.0, 2.0)
)

chem_pyrethroids <- do.call(rbind, lapply(chem_stations, function(stn) {
  data.frame(
    stationid      = stn,
    analytename    = pyrethroids,
    result         = as.character(pyr_vals[[stn]]),
    rl             = pyr_rl,
    mdl            = pyr_mdl,
    units          = "ug/kg dw",
    fieldrep       = 1,
    labrep         = 1,
    sampletypecode = "Result",
    stringsAsFactors = FALSE
  )
}))

# -- Fipronils (ug/dry kg) --
fipronils <- c("Fipronil", "Fipronil Desulfinyl", "Fipronil Sulfide", "Fipronil Sulfone")
fip_rl <- 0.5
fip_mdl <- 0.1

fip_vals <- list(
  "0000"     = rep(-88, 4),
  "S24-1001" = rep(-88, 4),
  "S24-1002" = rep(-88, 4),
  "S24-1003" = c(0.3, 0.1, 0.2, 0.15),
  "S24-1004" = c(1.5, 0.8, 1.0, 0.6),
  "S24-1005" = rep(-88, 4),
  "OFS-01"   = rep(-88, 4),
  "OFS-02"   = rep(-88, 4),
  "OFS-03"   = rep(-88, 4),
  "OFS-04"   = c(0.8, 0.4, 0.5, 0.3)
)

chem_fipronils <- do.call(rbind, lapply(chem_stations, function(stn) {
  data.frame(
    stationid      = stn,
    analytename    = fipronils,
    result         = as.character(fip_vals[[stn]]),
    rl             = fip_rl,
    mdl            = fip_mdl,
    units          = "ug/kg dw",
    fieldrep       = 1,
    labrep         = 1,
    sampletypecode = "Result",
    stringsAsFactors = FALSE
  )
}))

# -- PBDEs (ug/dry kg) --
pbdes <- c("PBDE-017", "PBDE-028", "PBDE-047", "PBDE-049", "PBDE-066",
           "PBDE-085", "PBDE-099", "PBDE-100", "PBDE-138", "PBDE-153",
           "PBDE-154", "PBDE-183", "PBDE-190")

pbde_rl <- 0.5
pbde_mdl <- 0.1

pbde_vals <- list(
  "0000"     = rep(-88, length(pbdes)),
  "S24-1001" = rep(-88, length(pbdes)),
  "S24-1002" = c(-88, -88, 0.5, -88, -88, -88, 0.3, 0.1, -88, -88, -88, -88, -88),
  "S24-1003" = c(-88, 0.2, 2.0, 0.3, 0.1, -88, 1.5, 0.5, 0.1, 0.2, 0.1, -88, -88),
  "S24-1004" = c(0.3, 0.8, 8.0, 1.2, 0.5, 0.3, 6.0, 2.0, 0.5, 0.8, 0.5, 0.2, 0.1),
  "S24-1005" = c(-88, 0.1, 1.0, 0.1, -88, -88, 0.8, 0.3, -88, 0.1, -88, -88, -88),
  "OFS-01"   = rep(-88, length(pbdes)),
  "OFS-02"   = c(-88, -88, 0.3, -88, -88, -88, 0.2, -88, -88, -88, -88, -88, -88),
  "OFS-03"   = c(-88, 0.1, 1.2, 0.2, -88, -88, 0.9, 0.3, -88, 0.1, -88, -88, -88),
  "OFS-04"   = c(0.2, 0.5, 5.0, 0.8, 0.3, 0.2, 3.5, 1.2, 0.3, 0.5, 0.3, 0.1, -88)
)

chem_pbdes <- do.call(rbind, lapply(chem_stations, function(stn) {
  data.frame(
    stationid      = stn,
    analytename    = pbdes,
    result         = as.character(pbde_vals[[stn]]),
    rl             = pbde_rl,
    mdl            = pbde_mdl,
    units          = "ug/kg dw",
    fieldrep       = 1,
    labrep         = 1,
    sampletypecode = "Result",
    stringsAsFactors = FALSE
  )
}))

# -- Sediment parameters --
sed_params <- c("TOC", "TN", "Lipids")
sed_rl <- c(0.01, 0.01, 0.1)
sed_mdl <- c(0.005, 0.005, 0.05)

sed_vals <- list(
  "0000"     = c(-88, -88, -88),
  "S24-1001" = c(0.5, 0.04, 0.1),
  "S24-1002" = c(1.2, 0.08, 0.2),
  "S24-1003" = c(2.0, 0.15, 0.4),
  "S24-1004" = c(3.5, 0.28, 0.8),
  "S24-1005" = c(1.5, 0.10, 0.3),
  "OFS-01"   = c(0.4, 0.03, 0.08),
  "OFS-02"   = c(0.9, 0.06, 0.15),
  "OFS-03"   = c(1.5, 0.10, 0.3),
  "OFS-04"   = c(2.8, 0.22, 0.6)
)

chem_sed <- do.call(rbind, lapply(chem_stations, function(stn) {
  data.frame(
    stationid      = stn,
    analytename    = sed_params,
    result         = as.character(sed_vals[[stn]]),
    rl             = sed_rl,
    mdl            = sed_mdl,
    units          = "%",
    fieldrep       = 1,
    labrep         = 1,
    sampletypecode = "Result",
    stringsAsFactors = FALSE
  )
}))

# -- Grain size (Phi values, percent) --
phi_names <- c("Phi -1.5", "Phi -1.0", "Phi -0.5", "Phi 0.0", "Phi 0.5",
               "Phi 1.0", "Phi 1.5", "Phi 2.0", "Phi 2.5", "Phi 3.0",
               "Phi 3.5", "Phi 4.0", "Phi 4.5", "Phi 5.0", "Phi 5.5",
               "Phi 6.0", "Phi 6.5", "Phi 7.0", "Phi 7.5", "Phi 8.0",
               "Phi 8.5", "Phi 9.0", "Phi 9.5", "Phi 10.0", "Phi 10.5",
               "Phi 11.0", "Phi 11.5", "Phi 12.0", "Phi 12.5", "Phi 13.0",
               "Phi 13.5", "Phi 14.0", "Phi 14.5")

phi_rl <- 0.01
phi_mdl <- 0.01

# Sandy site
phi_sandy <- c(0.5, 1.0, 2.0, 5.0, 10.0, 15.0, 20.0, 18.0, 12.0, 8.0,
               3.0, 2.0, 1.0, 0.5, 0.5, 0.3, 0.2, 0.2, 0.1, 0.1,
               0.1, 0.1, 0.1, 0.1, 0.05, 0.05, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0)
# Silty site
phi_silty <- c(0.0, 0.0, 0.5, 1.0, 2.0, 3.0, 4.0, 5.0, 5.0, 5.0,
               6.0, 7.0, 8.0, 8.0, 7.0, 6.0, 5.0, 5.0, 4.0, 4.0,
               3.0, 3.0, 2.5, 2.0, 1.5, 1.0, 0.5, 0.5, 0.3, 0.1, 0.05, 0.05, 0.0)

phi_vals <- list(
  "0000"     = rep(-88, length(phi_names)),
  "S24-1001" = phi_sandy,
  "S24-1002" = round((phi_sandy + phi_silty) / 2, 2),
  "S24-1003" = phi_silty,
  "S24-1004" = phi_silty * 1.02,
  "S24-1005" = round((phi_sandy * 0.7 + phi_silty * 0.3), 2),
  "OFS-01"   = phi_sandy * 0.95,
  "OFS-02"   = round((phi_sandy * 0.6 + phi_silty * 0.4), 2),
  "OFS-03"   = round((phi_sandy * 0.3 + phi_silty * 0.7), 2),
  "OFS-04"   = phi_silty * 0.98
)

chem_phi <- do.call(rbind, lapply(chem_stations, function(stn) {
  data.frame(
    stationid      = stn,
    analytename    = phi_names,
    result         = as.character(phi_vals[[stn]]),
    rl             = phi_rl,
    mdl            = phi_mdl,
    units          = "%",
    fieldrep       = 1,
    labrep         = 1,
    sampletypecode = "Result",
    stringsAsFactors = FALSE
  )
}))

# Combine all chemistry data
chem_sampledata <- rbind(
  chem_metals, chem_lpahs, chem_hpahs, chem_pesticides,
  chem_pcbs, chem_pyrethroids, chem_fipronils, chem_pbdes,
  chem_sed, chem_phi
)

cat("Chemistry sample data:", nrow(chem_sampledata), "rows x", ncol(chem_sampledata), "cols\n")
cat("  Stations:", length(unique(chem_sampledata$stationid)), "\n")
cat("  Analytes:", length(unique(chem_sampledata$analytename)), "\n")


# ============================================================================
# 2. TOXICITY SAMPLE DATA
# ============================================================================
# Two species:
#   Eohaustorius estuarius     -> Whole Sediment
#   Mytilus galloprovincialis   -> Sediment Water Interface
#
# Station "0000" = control (CNEG), others = Grab samples
# 5 lab replicates per station-species pair
# Stations match chemistry stations (minus control as standalone)

tox_stations <- c("S24-1001", "S24-1002", "S24-1003", "S24-1004", "S24-1005")

# Toxbatch links controls to their samples
tox_batches_ee <- c("06/24S24eea", "07/24S24eea", "08/24S24eea")
tox_batches_mg <- c("NOAAS24SedTest1", "NOAAS24SedTest2")

# Control results (high survival/normal development)
set.seed(42)

make_tox_rows <- function(stationid, toxbatch, species, sampletypecode, matrix, results, qacodes = rep(NA_character_, 5)) {
  data.frame(
    stationid = stationid,
    toxbatch = toxbatch,
    species = species,
    sampletypecode = sampletypecode,
    matrix = matrix,
    labrep = 1:5,
    result = results,
    qacode = qacodes,
    stringsAsFactors = FALSE
  )
}

tox_rows <- list()

# ---- Eohaustorius estuarius (Whole Sediment) ----
# Batch 1: controls + S24-1001, S24-1002
tox_rows[[1]] <- make_tox_rows("0000", tox_batches_ee[1], "Eohaustorius estuarius", "CNEG", "Whole Sediment",
                                c(100, 95, 100, 100, 95))
tox_rows[[2]] <- make_tox_rows("S24-1001", tox_batches_ee[1], "Eohaustorius estuarius", "Grab", "Whole Sediment",
                                c(95, 90, 95, 100, 90))
tox_rows[[3]] <- make_tox_rows("S24-1002", tox_batches_ee[1], "Eohaustorius estuarius", "Grab", "Whole Sediment",
                                c(85, 80, 75, 85, 80))

# Batch 2: controls + S24-1003, S24-1004
tox_rows[[4]] <- make_tox_rows("0000", tox_batches_ee[2], "Eohaustorius estuarius", "CNEG", "Whole Sediment",
                                c(100, 100, 95, 100, 100))
tox_rows[[5]] <- make_tox_rows("S24-1003", tox_batches_ee[2], "Eohaustorius estuarius", "Grab", "Whole Sediment",
                                c(50, 55, 45, 60, 40),
                                c(NA, NA, "X", NA, NA))
tox_rows[[6]] <- make_tox_rows("S24-1004", tox_batches_ee[2], "Eohaustorius estuarius", "Grab", "Whole Sediment",
                                c(10, 15, 5, 20, 0))

# Batch 3: controls + S24-1005
tox_rows[[7]] <- make_tox_rows("0000", tox_batches_ee[3], "Eohaustorius estuarius", "CNEG", "Whole Sediment",
                                c(95, 100, 100, 95, 100))
tox_rows[[8]] <- make_tox_rows("S24-1005", tox_batches_ee[3], "Eohaustorius estuarius", "Grab", "Whole Sediment",
                                c(70, 75, 65, 70, 75))

# ---- Mytilus galloprovincialis (Sediment Water Interface) ----
# Batch 1: controls + S24-1001, S24-1002, S24-1003
tox_rows[[9]] <- make_tox_rows("0000", tox_batches_mg[1], "Mytilus galloprovincialis", "CNEG", "Sediment Water Interface",
                                c(90, 85, 90, 95, 88))
tox_rows[[10]] <- make_tox_rows("S24-1001", tox_batches_mg[1], "Mytilus galloprovincialis", "Grab", "Sediment Water Interface",
                                 c(88, 85, 90, 82, 86))
tox_rows[[11]] <- make_tox_rows("S24-1002", tox_batches_mg[1], "Mytilus galloprovincialis", "Grab", "Sediment Water Interface",
                                 c(75, 70, 78, 72, 68),
                                 c(NA, "Q", NA, NA, NA))
tox_rows[[12]] <- make_tox_rows("S24-1003", tox_batches_mg[1], "Mytilus galloprovincialis", "Grab", "Sediment Water Interface",
                                 c(40, 35, 45, 38, 42))

# Batch 2: controls + S24-1004, S24-1005
tox_rows[[13]] <- make_tox_rows("0000", tox_batches_mg[2], "Mytilus galloprovincialis", "CNEG", "Sediment Water Interface",
                                 c(92, 88, 90, 94, 90))
tox_rows[[14]] <- make_tox_rows("S24-1004", tox_batches_mg[2], "Mytilus galloprovincialis", "Grab", "Sediment Water Interface",
                                 c(8, 12, 5, 15, 10),
                                 c(NA, NA, NA, "X", NA))
tox_rows[[15]] <- make_tox_rows("S24-1005", tox_batches_mg[2], "Mytilus galloprovincialis", "Grab", "Sediment Water Interface",
                                 c(60, 55, 65, 58, 62))

# ---- Offshore stations: Eohaustorius estuarius (Whole Sediment) ----
# Batch 4: controls + OFS-01, OFS-02
tox_batches_ee_ofs <- "09/24S24eea"
tox_rows[[16]] <- make_tox_rows("0000", tox_batches_ee_ofs, "Eohaustorius estuarius", "CNEG", "Whole Sediment",
                                 c(100, 95, 100, 95, 100))
tox_rows[[17]] <- make_tox_rows("OFS-01", tox_batches_ee_ofs, "Eohaustorius estuarius", "Grab", "Whole Sediment",
                                 c(90, 95, 85, 90, 95))
tox_rows[[18]] <- make_tox_rows("OFS-02", tox_batches_ee_ofs, "Eohaustorius estuarius", "Grab", "Whole Sediment",
                                 c(75, 80, 70, 78, 72))

# Batch 5: controls + OFS-03, OFS-04
tox_batches_ee_ofs2 <- "10/24S24eea"
tox_rows[[19]] <- make_tox_rows("0000", tox_batches_ee_ofs2, "Eohaustorius estuarius", "CNEG", "Whole Sediment",
                                 c(100, 100, 95, 100, 95))
tox_rows[[20]] <- make_tox_rows("OFS-03", tox_batches_ee_ofs2, "Eohaustorius estuarius", "Grab", "Whole Sediment",
                                 c(80, 85, 75, 82, 78))
tox_rows[[21]] <- make_tox_rows("OFS-04", tox_batches_ee_ofs2, "Eohaustorius estuarius", "Grab", "Whole Sediment",
                                 c(20, 25, 15, 30, 10))

# ---- Offshore stations: Mytilus galloprovincialis (Sediment Water Interface) ----
# Batch 3: controls + OFS-01, OFS-02, OFS-03, OFS-04
tox_batches_mg_ofs <- "NOAAS24SedTest3"
tox_rows[[22]] <- make_tox_rows("0000", tox_batches_mg_ofs, "Mytilus galloprovincialis", "CNEG", "Sediment Water Interface",
                                 c(90, 92, 88, 90, 94))
tox_rows[[23]] <- make_tox_rows("OFS-01", tox_batches_mg_ofs, "Mytilus galloprovincialis", "Grab", "Sediment Water Interface",
                                 c(85, 88, 82, 86, 84))
tox_rows[[24]] <- make_tox_rows("OFS-02", tox_batches_mg_ofs, "Mytilus galloprovincialis", "Grab", "Sediment Water Interface",
                                 c(68, 72, 65, 70, 66))
tox_rows[[25]] <- make_tox_rows("OFS-03", tox_batches_mg_ofs, "Mytilus galloprovincialis", "Grab", "Sediment Water Interface",
                                 c(72, 75, 68, 70, 74))
tox_rows[[26]] <- make_tox_rows("OFS-04", tox_batches_mg_ofs, "Mytilus galloprovincialis", "Grab", "Sediment Water Interface",
                                 c(15, 20, 10, 18, 12))

tox_sampledata <- do.call(rbind, tox_rows)

cat("\nToxicity sample data:", nrow(tox_sampledata), "rows x", ncol(tox_sampledata), "cols\n")
cat("  Stations:", length(unique(tox_sampledata$stationid)), "\n")
cat("  Species:", unique(tox_sampledata$species), "\n")


# ============================================================================
# 3. BENTHIC SAMPLE DATA
# ============================================================================
# Uses taxa from the SQO lookup list (xl_tool.SoCalLUList)
# Stations match chemistry/tox stations
# Column name: SampleDepth (for backward compat with existing benthic data)
# The benthic.sqo function lowercases all column names and expects "depth"

# Common SoCal embayment taxa with a range of tolerance scores
# Mix of sensitive, intermediate, and tolerant taxa
# Include molluscs, crustaceans, polychaetes, and others

# Reference/clean community
ref_taxa <- data.frame(
  Taxon = c(
    "Ampelisca cristata", "Ampelisca brevisimulata", "Amphideutopus oculatus",
    "Amphiodia sp", "Amphicteis scaphobranchiata", "Apoprionospio pygmaea",
    "Notomastus sp", "Leitoscoloplos pugettensis", "Heteromastus filiforme",
    "Lumbrineris californiensis", "Prionospio malmgreni",
    "Alia carinata", "Acteocina harpa", "Barleeia sp",
    "Photis californica", "Rhepoxynius abronius",
    "Edwardsia californica", "Fabricinuda limnicola",
    "Glycera americana", "Pherusa neopapillata",
    "Spiophanes bombyx", "Spiophanes duplex",
    "Mediomastus sp", "Cossura candida"
  ),
  stringsAsFactors = FALSE
)

# Disturbed community - more tolerant taxa
disturbed_taxa <- data.frame(
  Taxon = c(
    "Capitella capitata Cmplx", "Polydora nuchalis",
    "Dorvillea (Schistomeringos) annulata", "Exogone lourei",
    "Aphelochaeta glandaria Cmplx", "Aphelochaeta petersenae",
    "Anoplodactylus viridintestinalis", "Boccardiella hamata",
    "Grandidierella japonica", "Tubificidae",
    "Streblospio benedicti", "Tubificoides brownae",
    "Euchone limnicola", "Armandia brevis"
  ),
  stringsAsFactors = FALSE
)

make_benthic_station <- function(stationid, lat, lon, depth, sampledate,
                                  salinity, stratum, taxa_df, abundances) {
  n <- nrow(taxa_df)
  data.frame(
    StationID   = stationid,
    Replicate   = 1L,
    SampleDate  = as.POSIXct(sampledate),
    Latitude    = lat,
    Longitude   = lon,
    Depth       = depth,
    Taxon       = taxa_df$Taxon[1:n],
    Abundance   = as.integer(abundances[1:n]),
    Salinity    = salinity,
    Stratum     = stratum,
    Exclude     = "No",
    stringsAsFactors = FALSE
  )
}

benthic_rows <- list()

# S24-1001: Reference site in a bay (clean, diverse, many sensitive taxa)
benthic_rows[[1]] <- make_benthic_station(
  "S24-1001", 33.7512, -118.2165, 8.0, "2024-08-15", 34.2, "Bays",
  ref_taxa,
  c(15, 12, 8, 5, 3, 10, 20, 6, 4, 8, 12, 3, 2, 4, 6, 3, 5, 8, 2, 3, 7, 5, 10, 6)
)

# S24-1002: Low disturbance in an estuary (mostly reference taxa, some tolerant)
mixed_1002 <- rbind(ref_taxa[1:16, , drop = FALSE], disturbed_taxa[1:4, , drop = FALSE])
benthic_rows[[2]] <- make_benthic_station(
  "S24-1002", 33.7680, -118.1950, 5.0, "2024-08-20", 30.5, "Estuaries",
  mixed_1002,
  c(10, 8, 5, 3, 2, 8, 15, 4, 6, 5, 8, 2, 1, 3, 4, 2, 25, 10, 8, 5)
)

# S24-1003: Moderate disturbance in a port (mixed community)
mixed_1003 <- rbind(ref_taxa[c(7, 9, 10, 11, 19, 23, 24), , drop = FALSE],
                     disturbed_taxa)
benthic_rows[[3]] <- make_benthic_station(
  "S24-1003", 33.7450, -118.2750, 10.0, "2024-08-22", NA, "Ports",
  mixed_1003,
  c(8, 15, 12, 10, 3, 5, 4, 45, 20, 30, 35, 18, 12, 8, 15, 10, 6, 22, 8, 5, 3)
)

# S24-1004: High disturbance in a port (dominated by tolerant taxa)
high_dist <- rbind(disturbed_taxa, ref_taxa[c(7, 9), , drop = FALSE])
benthic_rows[[4]] <- make_benthic_station(
  "S24-1004", 33.7380, -118.2820, 6.0, "2024-08-25", NA, "Ports",
  high_dist,
  c(120, 80, 45, 95, 60, 30, 15, 25, 40, 55, 35, 20, 18, 10, 3, 5)
)

# S24-1005: Low-moderate disturbance in a marina
mixed_1005 <- rbind(ref_taxa[c(1:8, 11, 12, 18, 20:24), , drop = FALSE],
                     disturbed_taxa[c(1, 3, 4, 11), , drop = FALSE])
benthic_rows[[5]] <- make_benthic_station(
  "S24-1005", 33.8560, -118.3990, 4.0, "2024-08-28", 33.8, "Marinas",
  mixed_1005,
  c(8, 6, 4, 3, 2, 7, 12, 3, 6, 2, 5, 2, 3, 8, 5, 4, 15, 8, 12, 6)
)

benthic_sampledata <- do.call(rbind, benthic_rows)

cat("\nBenthic sample data:", nrow(benthic_sampledata), "rows x", ncol(benthic_sampledata), "cols\n")
cat("  Stations:", length(unique(benthic_sampledata$StationID)), "\n")
cat("  Taxa:", length(unique(benthic_sampledata$Taxon)), "\n")


# ============================================================================
# 4. OFFSHORE BENTHIC SAMPLE DATA (INFAUNA)
# ============================================================================
# Uses taxa from ed.14.ptaxa (the offshore p-code lookup)
# Various depths to test shallow/mid/deep zone calculations

offshore_taxa_shallow <- c(
  "Ampelisca cristata", "Rhepoxynius abronius",
  "Spiophanes bombyx", "Prionospio malmgreni",
  "Glycera americana", "Lumbrineris californiensis",
  "Mediomastus sp", "Heteromastus filiforme",
  "Pectinaria californiensis", "Notomastus sp"
)

offshore_taxa_mid <- c(
  "Adontorhina cyclia", "Ampelisca pacifica",
  "Amphicteis scaphobranchiata", "Axinopsida serricata",
  "Paradiopatra parva", "Pectinaria californiensis",
  "Scalibregma californicum", "Rhodine bitorquata",
  "Melinna oculata", "Eranno lagunae"
)

offshore_taxa_deep <- c(
  "Adontorhina cyclia", "Acila castrensis",
  "Paradiopatra parva", "Scalibregma californicum",
  "Melinna oculata", "Amphicteis scaphobranchiata",
  "Cirrophorus branchiatus", "Cossura candida",
  "Rhodine bitorquata", "Nicippe tumida"
)

make_offshore_infauna <- function(stationid, sampledate, taxa, abundances) {
  data.frame(
    stationid  = stationid,
    replicate  = 1L,
    sampledate = sampledate,
    taxon      = taxa,
    abundance  = as.integer(abundances),
    stringsAsFactors = FALSE
  )
}

off_rows <- list()

# OFS-01: Shallow site (24m)
off_rows[[1]] <- make_offshore_infauna(
  "OFS-01", "8/15/2024", offshore_taxa_shallow,
  c(12, 8, 15, 10, 3, 5, 20, 8, 4, 6)
)

# OFS-02: Mid-depth site (85m)
off_rows[[2]] <- make_offshore_infauna(
  "OFS-02", "8/18/2024", offshore_taxa_mid,
  c(5, 3, 2, 1, 10, 2, 2, 1, 1, 1)
)

# OFS-03: Deep site (180m)
off_rows[[3]] <- make_offshore_infauna(
  "OFS-03", "8/20/2024", offshore_taxa_deep,
  c(3, 2, 8, 1, 1, 2, 1, 4, 1, 1)
)

# OFS-04: Shallow site (28m) - more disturbed
off_rows[[4]] <- make_offshore_infauna(
  "OFS-04", "8/22/2024",
  c("Capitella capitata Cmplx", "Polydora nuchalis", "Streblospio benedicti",
    "Mediomastus sp", "Heteromastus filiforme", "Spiophanes bombyx",
    "Glycera americana", "Lumbrineris californiensis"),
  c(45, 30, 25, 15, 20, 8, 3, 5)
)

# OFS-05: Mid-depth site (60m) - moderate
off_rows[[5]] <- make_offshore_infauna(
  "OFS-05", "8/25/2024",
  c("Ampelisca pacifica", "Paradiopatra parva", "Adontorhina cyclia",
    "Pectinaria californiensis", "Eranno lagunae", "Rhodine bitorquata",
    "Capitella capitata Cmplx", "Polydora nuchalis"),
  c(4, 6, 3, 2, 2, 1, 10, 8)
)

# OFS-06: Deep site (200m) - reference
off_rows[[6]] <- make_offshore_infauna(
  "OFS-06", "8/28/2024",
  c("Acila castrensis", "Adontorhina cyclia", "Paradiopatra parva",
    "Scalibregma californicum", "Melinna oculata", "Cossura candida",
    "Cirrophorus branchiatus", "Amphicteis scaphobranchiata",
    "Nicippe tumida", "Rhodine bitorquata", "Axinopsida serricata",
    "Campylaspis rubromaculata"),
  c(4, 6, 12, 3, 2, 5, 2, 3, 1, 2, 1, 1)
)

offshore_infauna <- do.call(rbind, off_rows)

cat("\nOffshore infauna sample data:", nrow(offshore_infauna), "rows x", ncol(offshore_infauna), "cols\n")
cat("  Stations:", length(unique(offshore_infauna$stationid)), "\n")


# ============================================================================
# 5. OFFSHORE STATION DATA
# ============================================================================
offshore_station_info <- data.frame(
  stationid  = c("OFS-01", "OFS-02", "OFS-03", "OFS-04", "OFS-05", "OFS-06"),
  sampledate = c("8/15/2024", "8/18/2024", "8/20/2024", "8/22/2024", "8/25/2024", "8/28/2024"),
  latitude   = c(34.0237, 32.8075, 32.5859, 33.6954, 34.3870, 34.2053),
  longitude  = c(-118.5931, -117.3441, -117.3411, -118.2960, -119.8142, -119.5064),
  depth      = c(24.0, 85.0, 180.0, 28.0, 60.0, 200.0),
  stringsAsFactors = FALSE
)

cat("\nOffshore station info:", nrow(offshore_station_info), "rows x", ncol(offshore_station_info), "cols\n")


# Combine infauna with the per-station depth/lat/long into the single frame the package ships.
# BRI.Offshore (and SQOUnified, via offshore_benthic) expect one combined table, so the infauna and
# station info are joined here on stationid + sampledate.
offshore_bri_sampledata <- offshore_infauna %>%
  dplyr::left_join(offshore_station_info, by = c("stationid", "sampledate")) %>%
  dplyr::select(stationid, sampledate, replicate, taxon, abundance, latitude, longitude, depth)

cat("\nOffshore BRI sample data (combined):", nrow(offshore_bri_sampledata), "rows x", ncol(offshore_bri_sampledata), "cols\n")
cat("  Stations:", length(unique(offshore_bri_sampledata$stationid)), "\n")


# ============================================================================
# SAVE ALL DATASETS
# ============================================================================
save(chem_sampledata, file = "data/chem_sampledata.RData")
save(tox_sampledata, file = "data/tox_sampledata.RData")
save(benthic_sampledata, file = "data/benthic_sampledata.RData")
save(offshore_bri_sampledata, file = "data/offshore_bri_sampledata.RData")

cat("\nAll sample data saved to data/ directory.\n")
cat("\nSummary:\n")
cat("  chem_sampledata:           ", nrow(chem_sampledata), "rows\n")
cat("  tox_sampledata:            ", nrow(tox_sampledata), "rows\n")
cat("  benthic_sampledata:        ", nrow(benthic_sampledata), "rows\n")
cat("  offshore_bri_sampledata:   ", nrow(offshore_bri_sampledata), "rows\n")
