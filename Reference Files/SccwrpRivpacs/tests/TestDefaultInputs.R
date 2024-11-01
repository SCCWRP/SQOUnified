library(SccwrpRivpacs)

# Run using default user data.

socal <- SoCalRivpacs()
sfbay <- SFBayRivpacs()

HtmlOutput(rivpacs = socal, timestamp = "1359585702", user.filename = "test", path = "")
HtmlOutput(rivpacs = sfbay, timestamp = "1359585702", user.filename = "test", path = "")

# Checks of O/E output. Checks verified agaist first 4 digits of Utah website output.

# From dput(socal$oe.table$O.over.E)
socal.check <- c(0.760857895425563, 0.729367649061569, 0.915839501477777, 0.788327686332406, 
                 0.096393342963495, 0.827434126368914, 1.1449592266525, 0.701677608697632, 
                 0.960783736873942, 0.975227016442242, 0.693104263431818, 0.923388779236679, 
                 1.49924606553788, 0.896408353250279, 0.863144286004965)

stopifnot(round(socal$oe.table$O.over.E, digits = 4) == round(socal.check, digits = 4))

# From dput(sfbay$oe.table$O.over.E)
sfbay.check <- c(0.533834398961792, 0.246269609096255, 0.559903922144785, 0.394449800590905, 
                 0.608493467223863, 0.569864067418964, 0.395852343703887, 0.625397099646503, 
                 0.507117435125047, 0.761272753945334, 1.48277657855199, 0.394415191748079, 
                 0.586576813392697, 0.59771701859227, 0.507326979997167)

stopifnot(round(sfbay$oe.table$O.over.E, digits = 4) == round(sfbay.check, digits = 4))
