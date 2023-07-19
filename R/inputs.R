
#' clean up catch data
#'
#' @param year  year of assessment
#' @param species species of interest e.g., "SABL", "DUSK"
#' @param fishery identify the fishery default is "fsh"
#' @param TAC last three TAC in form: c(year-3, year-2, year-1)
#' @param discard if summarizing catch by discard/retained is desired change to TRUE
#' @param gear if summarizing catch by gear type is desired change to TRUE
#' @param fixed_catch if early catch is frozen place the file in user_input folder (format: year, catch)
#' @param alt alternate folder to save to - will be placed in "year/alt/data" folder
#' @param save
#'
#' @return
#' @export clean_catch
#'
#' @examples
#' \dontrun{
#' clean_catch(year, TAC = c(2874, 2756, 3100))
#' }
#'
clean_catch <- function(year, species, TAC = c(3333, 2222, 1111), discard = FALSE, gear = FALSE, fixed_catch = NULL, alt=NULL, save = TRUE){
  years = (year-3):(year-1)
  yr = year
  if(sum(TAC == c(3333, 2222, 1111)) == 3) {
    stop("check your TAC!")
  }

  if(!is.null(fixed_catch)){
    fc = vroom::vroom(here::here(year, "data", "user_input", fixed_catch))
  } else if(is.null(fixed_catch)){
    if(species == "NORK"){
      fc = afscdata::goa_nork_catch_1961_1992
    }
    if(species == "SABL"){
      fc = afscdata::sabl_fixed_abundance %>%
        tidytable::filter(variable == "catch")
    }
    if(species %in% c("REBS", "REYE")){
      fc = afscdata::goa_rebs_catch_1977_2004
    }
    if(species == "DUSK"){
      fc = afscdata::goa_dusk_catch_1977_1990
    }
    if(species %in% c("POP", "POPA")){
      fc = afscdata::goa_pop_catch_1960_1990
    }
    if(species %in% c("ATF", "ARTH")){
      fc = afscdata::goa_atf_catch_1961_1990
    }
    if(species == "FSOL"){
      fc = afscdata::goa_fhs_catch_1978_1990
    }
  }

  names(fc) <- c("year", "catch")

  # Fishery catch data ----
  vroom::vroom(here::here(year, "data", "raw", "fsh_catch_data.csv")) -> catch_data
  vroom::vroom(here::here(year, "data", "raw", "fsh_obs_data.txt"),
               delim = ",",
               col_type = c(join_key = "c", haul_join = "c")) -> obs_data

  # Estimate catch ratio in final year to end of year
  obs_data %>%
    tidytable::filter(year %in% years) %>%
    tidytable::mutate(tot_catch = sum(extrapolated_weight),
                      test_date = lubridate::`year<-`(max(catch_data$week_end_date), year), .by = year) %>%
    tidytable::filter(haul_date <= test_date) %>%
    tidytable::summarise(oct_catch = round(sum(extrapolated_weight)),
                         tot_catch = round(mean(tot_catch)), .by = year) %>%
    tidytable::summarise(ratio = 1 + (sum(tot_catch) - sum(oct_catch)) / sum(oct_catch)) -> rat

  # Compute catch
  if(nrow(fc)>=1){
    catch_data %>%
      tidytable::select(year, catch = weight_posted) %>%
      tidytable::filter(year > max(fc$year)) %>%
      tidytable::summarise(catch = round(sum(catch), 4), .by = year) %>%
      tidytable::bind_rows(fc) %>%
      tidytable::arrange(year) -> catch
  } else {
    catch_data %>%
      tidytable::select(year, catch = weight_posted) %>%
      tidytable::summarise(catch = round(sum(catch), 4), .by = year) %>%
      tidytable::arrange(year) -> catch
  }


  # estimate yield ratio of previous 3 years relative to TAC
  catch %>%
    tidytable::filter(year %in% years) %>%
    tidytable::bind_cols(tac = TAC) %>%
    tidytable::summarise(yld = mean(catch / tac)) -> yield

  # estimate catch through end of the year
  catch %>%
    tidytable::filter(year==yr) %>%
    tidytable::mutate(proj_catch = catch * rat$ratio) %>%
    tidytable::bind_cols(rat, yield) -> yld

  # replace current year's catch with projected catch
  catch$catch[length(catch$catch)] <- yld$proj_catch

  if(!(is.null(alt))) {
    vroom::vroom_write(catch, here::here(year, alt, "output",  "fsh_catch.csv"), delim = ",")
    vroom::vroom_write(yld, here::here(year, alt, "output", "yld_rat.csv"), delim = ",")
    catch
  }else if(isTRUE(save)){
    vroom::vroom_write(catch, here::here(year, "data", "output",  "fsh_catch.csv"), delim = ",")
    vroom::vroom_write(yld, here::here(year, "data", "output", "yld_rat.csv"), delim = ",")
    catch
  } else {
    list(catch = catch, yld_rat = yld)
  }

}

#' survey biomass
#'
#' @param year of interest
#' @param area options are bs, bsslope, nbs, ai, goa, old_bs - can only call a single area
#' @param type "depth", "stratum", "area", "total", "inpfc", "inpfc_depth" - only available for goa/ai (default: "total") - can only use a single switch
#' @param file if not using the design-based abundance, the file name must be stated (e.g. "GAP_VAST.csv")
#' @param rmv_yrs any survey years to exclude
#' @param alt alternate folder to save to - will be placed in "year/alt/data" folder
#' @param save save to default location
#'
#' @return
#' @export bts_biomass
#'
#' @examples
#'
bts_biomass <- function(year, area = "goa", type = "total", file = NULL, rmv_yrs = NULL, alt=NULL, save = TRUE){

  area = tolower(area)
  type = tolower(type)

  if(is.null(file)){

    vroom::vroom(here::here(year, "data", "raw", paste0(area, "_", type, "_bts_biomass_data.csv"))) %>%
      dplyr::rename_with(tolower)  -> df

    # sablefish are different...
    if("summary_depth" %in% names(df)){
      df %>%
        dplyr::filter(summary_depth < 995, year != 2001) %>%
        dplyr::group_by(year) %>%
        dplyr::summarise(biom = sum(area_biomass) / 1000,
                         se = sqrt(sum(biomass_var)) / 1000) %>%
        dplyr::mutate(lci = biom - 1.96 * se,
                      uci = biom + 1.96 * se) -> sb
    } else {
      df %>%
        tidytable::summarise(biomass = sum(total_biomass),
                             se = sqrt(sum(biomass_var)),
                             .by = year) %>%
        tidytable::mutate(lci = biomass - 1.96 * se,
                          uci = biomass + 1.96 * se,
                          lci = ifelse(lci < 0, 0, lci)) %>%
        tidytable::mutate(tidytable::across(tidytable::where(is.double), round)) %>%
        tidytable::filter(biomass > 0) -> sb
    }

  } else {
    vroom::vroom(here::here(year, "data", "user_input", file)) -> sb
  }

  if(!is.null(rmv_yrs)){
    sb <- dplyr::filter(sb, !(year %in% rmv_yrs))
  }

  if(!(is.null(alt))) {
    vroom::vroom_write(sb, here::here(year, alt, "data", paste(area, type, "bts_biomass.csv", sep="_")), ",")
    sb
  } else if(isTRUE(save)){
    vroom::vroom_write(sb, here::here(year, "data", "output", paste(area, type, "bts_biomass.csv", sep="_")), delim=",")
    sb
  }

}


#' aging error analysis
#'
#' @param read_tester = looks for a file in the user_input folder, e.g., :"reader_tester.csv"
#' @param species = "NORK"
#' @param year = year of the assessment
#' @param admb_home = location admb exists on your computer - if is "c:/admb" can leave NULL
#' @param area = "GOA" (BSAI not currently setup)
#' @param rec_age = recruitment age
#' @param plus_age = max age for modeling
#' @param max_age = max age for age error analysis - default = 100
#' @param alt alternate folder to save to - will be placed in "year/alt/data" folder
#' @param save = default is TRUE, FALSE outputs a list, names outputs to a data folder within a specific folder (e.g., save = "age_plus")
#'
#' @return
#' @export age_error
#'
#' @examples ageage(species = "NORK", year = 2020, admb_home = NULL, area = "GOA", rec_age = 2, plus_age = 45, max_age = 100)

age_error <- function(reader_tester, species, year, admb_home = NULL, area = "GOA", rec_age = 2, plus_age = 45, max_age = 100, alt = NULL, save = TRUE){


  rt = vroom::vroom(here::here(year, "data", "user_input", reader_tester))
  area = toupper(area)

  if(is.null(admb_home)){
    R2admb::setup_admb()
  } else {
    R2admb::setup_admb(admb_home)
  }

  nages = length(rec_age:plus_age)

  rt %>%
    tidytable::filter(Species %in% sp_switch(species),
                      Region == area,
                      Read_Age > 0,
                      Test_Age > 0,
                      Final_Age > 0) %>%
    tidytable::rename_with(tolower) %>%
    tidytable::summarise(freq = tidytable::n(),
                         .by = c(test_age, read_age)) -> dat


  tidytable::left_join(expand.grid(test_age = unique(dat$test_age),
                               read_age = unique(dat$read_age)),
                   dat) %>%
    tidytable::replace_na(list(freq = 0)) %>%
    tidytable::mutate(num = tidytable::case_when(test_age != read_age ~ freq),
                      .by = test_age) %>%
    tidytable::summarise(num = sum(num, na.rm = TRUE),
                     den = sum(freq),
                     .by = test_age) %>%
    tidytable::mutate(ape = 1 - (num / den),
                      ape = ifelse(is.nan(ape), 0, ape)) %>%
    tidytable::select(age = test_age, ape, ss = den) %>%
    tidytable::left_join(data.frame(age = min(dat$test_age):max(dat$test_age)), .) %>%
    tidytable::replace_na(list(ape = -9, ss = -9)) -> dats

  c("# Number of obs", nrow(dats),
    "# Age vector", dats$age,
    "# Percent agreement vector", dats$ape,
    "# Sample size vector", dats$ss) %>%
    write.table(here::here(year, "data", "models", "ageage", "ageage.dat"),
                sep="", quote=F, row.names=F, col.names=F)

  setwd(here::here(year, "data", "models", "ageage"))
  R2admb::compile_admb("ageage", verbose = TRUE)
  R2admb::run_admb("ageage", verbose=TRUE)

  setwd(here::here())
  read.delim(here::here(year, "data", "models", "ageage", "ageage.std"), sep="") %>%
    dplyr::filter(grepl("_a", name)) %>%
    dplyr::bind_cols(dats) %>%
    dplyr::select(age, value) -> sds

  fit = lm(value ~ age, data = sds)

  # fit out to age 100 (aka: max_age)
  data.frame(age = rec_age:max_age) %>%
    dplyr::mutate(ages_sd = predict(fit, .)) -> fits

  ages = fits$age

  mtx100 = matrix(nrow = nrow(fits), ncol = nrow(fits))
  colnames(mtx100) = rownames(mtx100) = ages

  for(j in 1:nrow(fits)){
    mtx100[j,1] = pnorm(ages[1] + 0.5,
                        ages[j],
                        fits[which(fits[,1] == ages[j]), 2])


    for(i in 2:(nrow(fits) - 1)){
      mtx100[j,i] = pnorm(ages[i] + 0.5,
                          ages[j],
                          fits[which(fits[,1] == ages[j]), 2]) -
        pnorm(ages[i-1] + 0.5,
              ages[j],
              fits[which(fits[,1] == ages[j]), 2])

    }
    mtx100[j,nrow(fits)] = 1 - sum(mtx100[j, 1:(nrow(fits) - 1)])
  }

  # Compute ageing error matrix for model
  ae_Mdl = matrix(nrow=length(ages), ncol=nages)
  ae_Mdl[, 1:(nages-1)] = as.matrix(mtx100[, 1:(nages-1)])
  ae_Mdl[, nages] = rowSums(mtx100[, nages:length(ages)])
  ae_Mdl = round(ae_Mdl, digits=4)
  r = which(ae_Mdl[, nages]>=0.999)[1]
  ae_Mdl = ae_Mdl[1:r,]


  if(!is.null(alt)){
    write.csv(mtx100, here::here(year, alt, "data", paste0("ae_mtx_", max_age, ".csv")), row.names = F)
    vroom::vroom_write(fits,  here::here(year,alt, "data", "ae_SD.csv"), ",")
    write.csv(ae_Mdl,  here::here(year, alt, "data", "ae_model.csv"), row.names = F)
    ae_Mdl
  } else if(isTRUE(save)){
    write.csv(mtx100, here::here(year, "data", "output", paste0("ae_mtx_", max_age, ".csv")), row.names = F)
    vroom::vroom_write(fits,  here::here(year,"data", "output", "ae_SD.csv"), ",")
    write.csv(ae_Mdl,  here::here(year, "data", "output", "ae_model.csv"), row.names = F)
    ae_Mdl
  } else {
    list(mtx100 = mtx100, ae_sd = fits, ae_model = ae_Mdl)
  }

}

#' size at age analysis
#'
#' @param year analysis year
#' @param area area data are from e.g., "goa"
#' @param admb_home location admb exists on your computer
#' @param rec_age recruitment age
#' @param max_age max age for age error analysis - default = 100
#' @param lenbins length bin file
#' @param alt #' @param alt alternate folder to save to - will be placed in "year/alt/data" folder
#' @param save save in the default location

#' @return
#' @export size_at_age
#'
#' @examples
size_at_age <- function(year, area, admb_home = NULL, rec_age, lenbins = NULL, alt=NULL, save = TRUE){


  if(is.null(admb_home)){
    R2admb::setup_admb()
  } else {
    R2admb::setup_admb(admb_home)
  }

  if (!file.exists(here::here(year, "data", "output", "ae_model.csv"))){
    stop("You must first run the age-error function 'ageage()")
  } else {
    nages_m = nrow(read.csv(here::here(year, "data", "output", "ae_model.csv")))
  }


  ages_m = rec_age:(rec_age + nages_m - 1)

  if(is.null(lenbins)){
    stop("Please provide a vector of length buns or the file that is in the user_input folder e.g.,('lengthbins.csv') with a column names 'len_bins'")
  }

  if(!is.vector(lenbins)){
    lenbins =  vroom::vroom(here::here(year, "data", "user_input", lenbins), delim = ",")$len_bins
  }

  vroom::vroom(here::here(year, "data", "raw", "bts_specimen_data.csv")) %>%
    tidytable::rename_with(tolower) %>%
    tidytable::select(year, age, length) %>%
    tidytable::filter(year>=1990, !is.na(age))  %>%
    tidytable::select(-year) %>%
    tidytable::filter(tidytable::n()>1, .by = age) %>%
    tidytable::mutate(n_l = tidytable::n(), .by = length) %>%
    tidytable::arrange(age, length) %>%
    tidytable::mutate(sample_size = tidytable::n(), .by = age) -> inter

  vroom::vroom(here::here(year, "data", "raw", paste0("bts_length_data.csv"))) %>%
    tidytable::rename_with(tolower) %>%
    tidytable::filter(year>=1990, !is.na(length)) -> dat

  if(is.null(dat$frequency)){
    dat %>%
      tidytable::summarise(tot = tidytable::n(), .by = length) -> dat
  } else {
    dat %>%
      tidytable::select(frequency, length) %>%
      tidytable::summarise(tot = sum(frequency), .by = length) -> dat
  }

  dat %>%
    tidytable::left_join(inter, .) %>%
    tidytable::mutate(prop =  dplyr::n() / n_l * tot, .by = c(age, length)) %>%
    tidytable::distinct() %>%
    tidytable::mutate(sample_size = mean(sample_size),
                      Lbar = sum(prop * length) / sum(prop) * 0.1,
                      .by=age) %>%
    tidytable::summarise(SD_Lbar = sqrt(1 / (sum(prop) - 1) * sum(prop * (length / 10 - Lbar)^2)),
                         sample_size=mean(sample_size),
                         Lbar = mean(Lbar),
                         .by = age) %>%
    tidytable::filter(SD_Lbar>=0.01) -> laa_stats

  if(!is.null(alt)){
    vroom::vroom_write(laa_stats, here::here(year, alt, "data", "laa_stats.csv"))
  } else {
    vroom::vroom_write(laa_stats, here::here(year, "data", "output", "laa_stats.csv"))
  }

  laa_stats


  # run models ----

  setwd(here::here(year, "data", "models", "vonb"))
  # Estimate mean length
  c("# Data file for LVB model of mean length",
    "# Number of ages (nages)",
    nrow(laa_stats),
    "# Ages with observed mean length (ages)",
    paste(laa_stats$age, collapse=" "),
    "# Observed mean length (Lbar_obs)",
    paste(laa_stats$Lbar, collapse=" "),
    "# SD in Observed mean length (Lbar_obs)",
    paste(laa_stats$SD_Lbar, collapse=" ")) %>%
    write.table("vbl.dat", quote=FALSE, row.names=FALSE, col.names=FALSE)

  R2admb::compile_admb("vbl", verbose = TRUE)
  R2admb::run_admb("vbl", verbose = TRUE)

  # retrieve output

  REP = readLines("vbl.rep", warn=FALSE)
  Linf = as.numeric(sub(".*? ", "", REP[1]))
  k = as.numeric(sub(".*? ", "", REP[2]))
  t0 = as.numeric(sub(".*? ", "", REP[3]))


  # run model 2
  setwd(here::here(year, "data", "models", "length_sd"))

  c("# Data file for LVB model of mean length",
    "# Number of ages (nages)",
    nrow(laa_stats),
    "# Ages with observed mean length (ages)",
    paste(laa_stats$age, collapse=" "),
    "# SD in Observed mean length (Lbar_obs)",
    paste(laa_stats$SD_Lbar, collapse=" "),
    "# Sample size vector",
    paste(laa_stats$sample_size, collapse=" ")) %>%
    write.table("lengthsd.dat", quote=FALSE, row.names=FALSE, col.names=FALSE)


  R2admb::compile_admb("lengthsd", verbose = TRUE)
  R2admb::run_admb("lengthsd", verbose = TRUE)
  STD <- read.delim("lengthsd.std", sep="")
  a <- STD$value[1]
  b <- STD$value[2]
  (params <- cbind(Linf, k, t0, a, b))

  if(!is.null(alt)) {
    write.csv(params, here::here(year, alt, "data", "lbar_params.csv"), row.names=FALSE)
  } else {
    write.csv(params, here::here(year, "data", "output", "lbar_params.csv"), row.names=FALSE)
  }



  # Compute Sz@A transition matrix

  expand.grid(age = ages_m,
              length = lenbins) %>%
    dplyr::mutate(Lbar = Linf * (1 - exp(-k * (age - t0))),
                  Lbar = ifelse(age == max(ages_m), 0.5 * (Lbar + Linf), Lbar),
                  SD_Lbar = a * log(age) + b,
                  prob = ifelse(length == min(length),
                                pnorm(length + 0.5, Lbar, SD_Lbar),
                                pnorm(length + 0.5, Lbar, SD_Lbar) -
                                  pnorm(length -0.5, Lbar, SD_Lbar)),
                  prob = round(prob, digits = 4)) %>%
    dplyr::select(age, length, prob) %>%
    tidyr::pivot_wider(names_from = length, values_from = prob) %>%
    dplyr::mutate(!!rev(names(.))[1] := 1 - rowSums(.[2:(ncol(.) - 1)])) %>%
    dplyr::mutate_at(2:ncol(.), round, 4) -> saa

  if(!is.null(alt)){
    vroom::vroom_write(saa, here::here(year, alt, "data", "saa.csv"), ",")
    saa
  } else if(isTRUE(save)){
    vroom::vroom_write(saa, here::here(year, "data", "output", "saa.csv"), ",")
    saa
  } else {
    saa
  }


}

#' size at age for 1960s POP
#'
#' @param year analysis year
#' @param rec_age recruitment age
#' @param lenbins length bin file - this should be inyour user_input folder
#' @return
#' @export size_at_age_pop_60
#'
#' @examples
size_at_age_pop_60 <- function(year, rec_age, lenbins = NULL){

  if (!file.exists(here::here(year,"data", "output", "ae_model.csv"))){
    stop("You must first run the age-error function 'ageage()")
  } else {
    nages_m = nrow(read.csv(here::here(year, "data", "output", "ae_model.csv")))
    ages_m = rec_age:(rec_age + nages_m - 1)
  }

  pars = vroom::vroom(here::here(year, "data", "user_input", "saa_pop_60.csv"), delim = "\t")

  expand.grid(age = ages_m,
              length = lenbins) %>%
    dplyr::mutate(lbar = pars$linf * (1 - exp(-pars$k * (age - pars$t0))),
                  lbar = ifelse(age == max(ages_m), 0.5 * (lbar + pars$linf), lbar),
                  sd_lbar = pars$a * log(age) + pars$b,
                  prob = ifelse(length == min(length),
                                pnorm(length + 0.5, lbar, sd_lbar),
                                pnorm(length + 0.5, lbar, sd_lbar) -
                                  pnorm(length -0.5, lbar, sd_lbar)),
                  prob = round(prob, digits = 4)) %>%
    dplyr::select(age, length, prob) %>%
    tidyr::pivot_wider(names_from = length, values_from = prob) %>%
    dplyr::mutate(!!rev(names(.))[1] := 1 - rowSums(.[2:(ncol(.) - 1)])) %>%
    dplyr::mutate_at(2:ncol(.), round, 4) -> saa

  write.csv(saa, here::here(year, "data", "output", "saa_60.csv"), row.names = FALSE)
  saa

}

#' fishery age composition analysis
#'
#' @param year assessment year
#' @param fishery default is fsh, change if age comps from multiple fisheries (e.g., fsh1, fsh2)
#' @param exp_meth expansion method: marg - use marginal ages, marg_len - expand by marginal lengths thru alk, exp_len - expand by expanded lengths thru alk
#' @param rec_age recruitment age
#' @param plus_age plus age group
#' @param lenbins length bins in comp data
#' @param rmv_yrs any years to remove form the age comp e.g. c(1987, 1989)
#' @param id id a specific comp name - will be placed at end of file name e.g., id='use' will create 'fsh_age_comp-use.csv' in the data/output folder
#' @param save whether to save the file - wll be placed in "year/data/output" folder
#'
#' @return
#' @export  fish_age_comp
#'
#' @examples
#' \dontrun{
#' fish_age_comp(year, fishery = "fsh", rec_age, plus_age)
#' }
fish_age_comp <- function(year, fishery = "fsh", exp_meth, rec_age, plus_age, lenbins = NULL, rmv_yrs = NULL, id = NULL, save = TRUE){

  # compute age comps with marginal ages
  if(exp_meth == 'marg'){
    vroom::vroom(here::here(year, "data", "raw", paste0(fishery, "_specimen_data.txt")),
                 delim = ",",
                 col_type = c(join_key = "c", haul_join = "c", port_join = "c")) %>%
      tidytable::filter(age >= rec_age,
                        !(year %in% rmv_yrs),
                        !is.na(length),
                        !is.na(performance)) %>%
      tidytable::mutate(age = ifelse(age > plus_age, plus_age, age)) %>%
      tidytable::mutate(tot = tidytable::n(), .by = year) %>%
      tidytable::filter(tot>49) %>%
      tidytable::mutate(n_h = length(unique(na.omit(haul_join))) +
                          length(unique(na.omit(port_join))),
                        .by = year) %>%
      tidytable::summarise(n_s = mean(tot),
                           n_h = mean(n_h),
                           age_tot = tidytable::n(),
                           .by = c(year, age)) %>%
      tidytable::mutate(prop = age_tot / n_s) %>%
      tidytable::left_join(expand.grid(year = unique(.$year),
                                       age = rec_age:plus_age), .) %>%
      tidytable::replace_na(list(prop = 0)) %>%
      tidytable::mutate(AA_Index = 1,
                        n_s = mean(n_s, na.rm = T),
                        n_h = mean(n_h, na.rm = T),
                        .by = year) %>%
      tidytable::select(-age_tot) %>%
      tidytable::pivot_wider(names_from = age, values_from = prop) -> fac
  }

  # expand age comps with marginal lengths and alk
  if(exp_meth == 'marg_len'){

    # get marginal length comp
    vroom::vroom(here::here(year, "data", "raw", paste0(fishery, "_length_data.txt")),
                 delim = ",",
                 col_type = c(haul_join = "c", port_join = "c")) %>%
      tidytable::filter(!(year %in% rmv_yrs),
                        !is.na(length),
                        !is.na(performance)) %>%
      tidytable::mutate(length = ifelse(length >= max(lenbins), max(lenbins), length),
                        .by = year) %>%
      tidytable::summarise(length_tot = sum(frequency),
                           .by = c(year, length)) -> mar_len1

    vroom::vroom(here::here(year, "data", "raw", paste0(fishery, "_specimen_data.txt")),
                 delim = ",",
                 col_type = c(join_key = "c", haul_join = "c", port_join = "c")) %>%
      tidytable::filter(age >= rec_age,
                        !(year %in% rmv_yrs),
                        !is.na(length),
                        !is.na(performance)) %>%
      tidytable::distinct(year, length) %>%
      tidytable::left_join(mar_len1) %>%
      tidytable::drop_na() %>%
      tidytable::mutate(prop_l = length_tot / sum(length_tot), .by = c(year)) -> mar_len

    # get age comp expanded by marginal length comp and alk
    vroom::vroom(here::here(year, "data", "raw", paste0(fishery, "_specimen_data.txt")),
                 delim = ",",
                 col_type = c(join_key = "c", haul_join = "c", port_join = "c")) %>%
      tidytable::filter(age >= rec_age,
                        !(year %in% rmv_yrs),
                        !is.na(length),
                        !is.na(performance)) %>%
      tidytable::mutate(n_l = tidytable::n(), .by = c(year, age, length)) %>%
      tidytable::select(year, age, length, n_l) %>%
      tidytable::arrange(age, length) %>%
      tidytable::summarise(n_l = mean(n_l), .by = c(year, age, length)) %>%
      tidytable::mutate(N_l = sum(n_l), .by = c(year, length)) %>%
      tidytable::mutate(prop_al = n_l / N_l) %>%
      tidytable::left_join(mar_len) %>%
      tidytable::drop_na() %>%
      tidytable::mutate(prop = prop_al * prop_l) %>%
      tidytable::mutate(age = ifelse(age > plus_age, plus_age, age)) %>%
      tidytable::summarise(prop = sum(prop), .by = c(year, age)) -> age_comp

    # put it all together
    vroom::vroom(here::here(year, "data", "raw", paste0(fishery, "_specimen_data.txt")),
                 delim = ",",
                 col_type = c(join_key="c", haul_join="c", port_join="c")) %>%
      tidytable::mutate(tot = tidytable::n(), .by = year) %>%
      tidytable::filter(tot > 49,
                        !(year %in% rmv_yrs)) %>%
      tidytable::mutate(n_h = length(unique(na.omit(haul_join))) +
                          length(unique(na.omit(port_join))),
                        .by = year) %>%
      tidytable::summarise(n_s = mean(tot),
                           n_h = mean(n_h),
                           .by = c(year)) %>%
      tidytable::mutate(AA_Index = 1) %>%
      tidytable::left_join.(expand.grid(year = unique(.$year),
                                        age = rec_age:plus_age), .) %>%
      tidytable::left_join(age_comp) %>%
      tidytable::replace_na(list(prop = 0)) %>%
      tidytable::pivot_wider(names_from = age, values_from = prop) -> fac
  }

  # expand age comps with expanded lengths and alk
  if(exp_meth == 'exp_len'){

    # get expanded length comp (weighted by observer catch)
    vroom::vroom(here::here(year, "data", "raw", paste0(fishery, "_length_data.txt")),
                 delim = ",",
                 col_type = c(haul_join = "c", port_join = "c")) %>%
      tidytable::filter(!(year %in% rmv_yrs),
                        !is.na(length),
                        !is.na(performance)) %>%
      tidytable::left_join(vroom::vroom(here::here(year, "data", "raw", "fsh_obs_data.txt"),
                                        delim = ",",
                                        col_type = c(join_key = "c", haul_join = "c")) %>%
                             tidytable::select(haul_join, extrapolated_number)) %>%
      tidytable::select(year, haul_join, length, frequency, extrapolated_number) %>%
      tidytable::drop_na() %>%
      tidytable::distinct(year, haul_join, extrapolated_number) %>%
      tidytable::mutate(p_haul = extrapolated_number / sum(extrapolated_number),
                        .by = c(year)) -> p_haul # proportion of catch across hauls sampled for length

    vroom::vroom(here::here(year, "data", "raw", paste0(fishery, "_length_data.txt")),
                 delim = ",",
                 col_type = c(haul_join = "c", port_join = "c")) %>%
      tidytable::filter(!(year %in% rmv_yrs),
                        !is.na(length),
                        !is.na(performance)) %>%
      tidytable::left_join(p_haul) %>%
      tidytable::select(year, haul_join, length, frequency, p_haul) %>%
      tidytable::drop_na() %>%
      tidytable::mutate(length = ifelse(length >= max(lenbins), max(lenbins), length),
                        .by = year) %>%
      tidytable::mutate(p_hlen = frequency / sum(frequency), .by = c (year, haul_join)) %>% # compute haul length comps
      tidytable::mutate(n_len = sum(frequency), .by = c(year, length)) %>% # compute number of samples by length bin
      tidytable::mutate(wtd_freq = p_haul * p_hlen * n_len) %>% # compute weighted length frequencies per haul
      tidytable::summarise(length_tot = sum(wtd_freq), .by = c(year, length)) -> exp_len # expanded length frequencies weighted by haul catch

    exp_len %>%
      tidytable::distinct(year) -> exp_yrs

    # compute marginal length frequencies in years where observer catch data not available for expansion
    vroom::vroom(here::here(year, "data", "raw", paste0(fishery, "_length_data.txt")),
                 delim = ",",
                 col_type = c(haul_join = "c", port_join = "c")) %>%
      tidytable::filter(!(year %in% rmv_yrs),
                        !(year %in% exp_yrs),
                        !is.na(length),
                        !is.na(performance)) %>%
      tidytable::mutate(length = ifelse(length >= max(lenbins), max(lenbins), length),
                        .by = year) %>%
      tidytable::summarise(length_tot = sum(frequency),
                           .by = c(year, length)) -> mar_len

    # combine marginal and expanded frequencies
    mar_len %>%
      tidytable::bind_rows(exp_len) -> comb_len

    # compute length comps
    vroom::vroom(here::here(year, "data", "raw", paste0(fishery, "_specimen_data.txt")),
                 delim = ",",
                 col_type = c(join_key = "c", haul_join = "c", port_join = "c")) %>%
      tidytable::filter(age >= rec_age,
                        !(year %in% rmv_yrs),
                        !is.na(length),
                        !is.na(performance)) %>%
      tidytable::distinct(year, length) %>%
      tidytable::left_join(comb_len) %>%
      tidytable::drop_na() %>%
      tidytable::mutate(prop_l = length_tot / sum(length_tot), .by = c(year)) -> len_comp

    # get age comp expanded by length comp and alk
    vroom::vroom(here::here(year, "data", "raw", paste0(fishery, "_specimen_data.txt")),
                 delim = ",",
                 col_type = c(join_key = "c", haul_join = "c", port_join = "c")) %>%
      tidytable::filter(age >= rec_age,
                        !(year %in% rmv_yrs),
                        !is.na(length),
                        !is.na(performance)) %>%
      tidytable::mutate(n_l = tidytable::n(), .by = c(year, age, length)) %>%
      tidytable::select(year, age, length, n_l) %>%
      tidytable::arrange(age, length) %>%
      tidytable::summarise(n_l = mean(n_l), .by = c(year, age, length)) %>%
      tidytable::mutate(N_l = sum(n_l), .by = c(year, length)) %>%
      tidytable::mutate(prop_al = n_l / N_l) %>%
      tidytable::left_join(len_comp) %>%
      tidytable::drop_na() %>%
      tidytable::mutate(prop = prop_al * prop_l) %>%
      tidytable::mutate(age = ifelse(age > plus_age, plus_age, age)) %>%
      tidytable::summarise(prop = sum(prop), .by = c(year, age)) -> age_comp

    # put it all together
    vroom::vroom(here::here(year, "data", "raw", paste0(fishery, "_specimen_data.txt")),
                 delim = ",",
                 col_type = c(join_key="c", haul_join="c", port_join="c")) %>%
      tidytable::mutate(tot = tidytable::n(), .by = year) %>%
      tidytable::filter(tot > 49,
                        !(year %in% rmv_yrs)) %>%
      tidytable::mutate(n_h = length(unique(na.omit(haul_join))) +
                          length(unique(na.omit(port_join))),
                        .by = year) %>%
      tidytable::summarise(n_s = mean(tot),
                           n_h = mean(n_h),
                           .by = c(year)) %>%
      tidytable::mutate(AA_Index = 1) %>%
      tidytable::left_join(expand.grid(year = unique(.$year),
                                       age = rec_age:plus_age), .) %>%
      tidytable::left_join(age_comp) %>%
      tidytable::replace_na(list(prop = 0)) %>%
      tidytable::pivot_wider(names_from = age, values_from = prop) -> fac
  }

  if(!is.null(id)) {
    vroom::vroom_write(fac, here::here(year, "data", "output", paste0(fishery, "_age_comp-", id, ".csv")), ",")
    fac
  } else if(isTRUE(save)) {
    vroom::vroom_write(fac, here::here(year, "data", "output", paste0(fishery, "_age_comp.csv")), ",")
    fac
  } else {
    fac
  }

}

#' trawl survey age comp analysis
#'
#' @param year assessment year
#' @param area default is "goa"
#' @param rec_age recruitment age
#' @param plus_age plus group age
#' @param rmv_yrs any survey years to exclude
#' @param alt alternate folder to save to - will be placed in "year/alt/data" folder
#' @param save save in the default location
#'
#' @return
#' @export bts_age_comp
#'
#' @examples bts_age_comp(year = 2020, rec_age = 2, plus_age = 45)
bts_age_comp <- function(year, area = "goa", rec_age, plus_age, rmv_yrs = NULL, alt = NULL, save = TRUE){

  read.csv(here::here(year, "data", "raw", "bts_specimen_data.csv")) %>%
    dplyr::filter(!is.na(age)) %>%
    dplyr::group_by(year) %>%
    dplyr::summarise(n_s = dplyr::n(),
                     n_h = length(unique(hauljoin))) -> dat1


  read.csv(here::here(year, "data", "raw", paste0(area, "_total_bts_agecomp_data.csv"))) %>%
    dplyr::rename_with(tolower) %>%
    tidytable::filter.(age >= rec_age) %>%
    tidytable::mutate.(tot = sum(round(population_count, digits = 0)),
                       age = ifelse(age < plus_age, age, plus_age),
                       .by = year) %>%
    tidytable::summarise.(prop = sum(population_count) / mean(tot),
                          .by = c(age, year)) %>%
    tidytable::left_join.(dat1) %>%
    tidytable::left_join.(expand.grid(year = unique(.$year),
                                      age = rec_age:plus_age), .) %>%
    tidytable::replace_na.(list(prop = 0)) %>%
    tidytable::mutate.(AA_Index = 1,
                       n_s = mean(n_s, na.rm = T),
                       n_h = mean(n_h, na.rm = T),
                       .by = year) %>%
    tidytable::pivot_wider.(names_from = age, values_from = prop) %>%
    tidytable::arrange.(year) -> age_comp


  if(!is.null(rmv_yrs)){
    age_comp  %>%
      tidytable::filter.(!(year %in% rmv_yrs)) -> age_comp
  }

  if(!is.null(alt)) {
    vroom::vroom_write(age_comp, here::here(year, alt, "data", paste0(area, "_bts_age_comp.csv")), ",")
    age_comp
  } else if(isTRUE(save)){
    vroom::vroom_write(age_comp, here::here(year, "data", "output", paste0(area, "_bts_age_comp.csv")), ",")
    age_comp
  }

  age_comp

}

#' fishery length composition analysis
#'
#' @param year assessment year
#' @param fishery default is "fsh"
#' @param lenbins lenbin file if left NULL it looks for here::here(year, "data", "user_input", "len_bin_labels.csv")
#' @param rec_age recruitment age
#' @param lenbins length bins for composition data
#' @param fmv_yrs years to remove from comp data
#' @param rmv_yrs any survey years to exclude
#' @param id id a specific comp name - will be placed at end of file name e.g., id='use' will create 'fsh_length_comp-use.csv' in the data/output folder
#' @param save
#'
#' @return
#' @export fish_length_comp
#'
#' @examples
fish_length_comp <- function(year, fishery = "fsh", rec_age, lenbins = NULL, rmv_yrs = NULL, id = NULL, save = TRUE){

  if(is.null(lenbins)){
    stop("Please provide a vector of length buns or the file that is in the user_input folder e.g.,('lengthbins.csv') with a column names 'len_bins'")
  }

  if(!is.vector(lenbins)){
    lenbins =  vroom::vroom(here::here(year, "data", "user_input", lenbins), delim = ",")$len_bins
  }

  yr = year
  vroom::vroom(here::here(year, "data", "raw", paste0(fishery, "_specimen_data.txt")),
               delim = ",",
               col_type = c(join_key="c", haul_join="c", port_join="c")) %>%
    tidytable::filter(!is.na(age), age>=rec_age) %>%
    tidytable::group_by(year) %>%
    tidytable::tally(name = "age") %>%
    tidytable::filter(age >= 50) %>%
    tidytable::ungroup() -> ages

  vroom::vroom(here::here(year, "data", "raw", paste0(fishery,"_length_data.txt")),
               delim = ",",
               col_type = c(haul_join="c", port_join="c")) %>%
    tidytable::filter(!(year %in% unique(ages$year))) %>%
    tidytable::mutate(tot = sum(frequency),
                      length = ifelse(length >= max(lenbins), max(lenbins), length),
                      n_h = length(unique(na.omit(haul_join))) + length(unique(na.omit(port_join))),
                      .by = year) %>%
    tidytable::summarise(n_s = mean(tot),
                         n_h = mean(n_h),
                         length_tot = sum(frequency),
                         .by = c(year, length)) %>%
    tidytable::mutate(prop = length_tot / n_s) %>%
    tidytable::left_join(expand.grid(year = unique(.$year), length = lenbins), .) %>%
    tidytable::replace_na(list(prop = 0)) %>%
    tidytable::mutate(SA_Index = 1,
                      n_s = mean(n_s, na.rm = T),
                      n_h = mean(n_h, na.rm = T),
                      .by = year) %>%
    tidytable::select(-length_tot) %>%
    tidytable::pivot_wider(names_from = length, values_from = prop) -> flc

  if(!is.null(rmv_yrs)){
    flc  %>%
      tidytable::filter(!(year %in% rmv_yrs)) -> flc
  }


  if(!is.null(id)) {
    vroom::vroom_write(flc, here::here(year, "data", "output", paste0(fishery, "_length_comp-", id, ".csv")), ",")
    flc
  } else if(isTRUE(save)){
    vroom::vroom_write(flc, here::here(year, "data", "output", paste0(fishery, "_length_comp.csv")), ",")
    flc
  }
  flc


}

#' fishery length composition analysis for pop
#'
#' @param year assessment year
#' @param fishery default is "fsh"
#' @param lenbins lenbin file if left NULL it looks for here::here(year, "data", "user_input", "len_bin_labels.csv")
#' @param rec_age recruitment age
#' @param lenbins length bins for composition data
#' @param fmv_yrs years to remove from comp data
#' @param id id a specific comp name - will be placed at end of file name e.g., id='use' will create 'fsh_length_comp-use.csv' in the data/output folder
#' @param save
#'
#' @return
#' @export fish_length_comp_pop
#'
#' @examples
fish_length_comp_pop <- function(year, fishery = "fsh", rec_age, lenbins = NULL, rmv_yrs = NULL, id = NULL, save = TRUE){

  if(is.null(lenbins)){
    stop("Please provide a vector of length buns or the file that is in the user_input folder e.g.,('lengthbins.csv') with a column names 'len_bins'")
  }

  if(!is.vector(lenbins)){
    lenbins =  vroom::vroom(here::here(year, "data", "user_input", lenbins), delim = ",")$len_bins
  }

  # get current data
  yr = year
  vroom::vroom(here::here(year, "data", "raw", paste0(fishery, "_specimen_data.txt")),
               delim = ",",
               col_type = c(join_key="c", haul_join="c", port_join="c")) %>%
    tidytable::filter(!is.na(age), age>=rec_age) %>%
    tidytable::group_by(year) %>%
    tidytable::tally(name = "age") %>%
    tidytable::filter(age >= 50) %>%
    tidytable::ungroup() -> ages

  vroom::vroom(here::here(year, "data", "raw", paste0(fishery,"_length_data.txt")),
               delim = ",",
               col_type = c(haul_join="c", port_join="c")) %>%
    tidytable::filter(!(year %in% unique(ages$year)),
                      length >= 11) %>%
    tidytable::drop_na(haul_join) %>%
    tidytable::mutate(tot = sum(frequency),
                      length = ifelse(length >= max(lenbins), max(lenbins), length),
                      length = ifelse(length <= min(lenbins), min(lenbins), length),
                      n_h = length(unique(na.omit(haul_join))) + length(unique(na.omit(port_join))),
                      .by = year) %>%
    tidytable::summarise(n_s = mean(tot),
                         n_h = mean(n_h),
                         length_tot = sum(frequency),
                         .by = c(year, length)) %>%
    tidytable::mutate(prop = length_tot / n_s) %>%
    tidytable::left_join(expand.grid(year = unique(.$year), length = lenbins), .) %>%
    tidytable::replace_na(list(prop = 0)) %>%
    tidytable::mutate(SA_Index = 2,
                      n_s = mean(n_s, na.rm = T),
                      n_h = mean(n_h, na.rm = T),
                      .by = year) %>%
    tidytable::select(-length_tot) %>%
    tidytable::pivot_wider(names_from = length, values_from = prop) -> flc

  if(!is.null(rmv_yrs)){
    flc  %>%
      tidytable::filter(!(year %in% rmv_yrs)) -> flc
  }


  # get historical data
  vroom::vroom(here::here(year, "data", "user_input", "goa_pop_fixed_fish_length_comp.csv"), delim = "\t") %>%
    dplyr::rename_with(tolower) %>%
    tidytable::mutate(tot = sum(value),
                      length = ifelse(length >= max(lenbins), max(lenbins), length),
                      length = ifelse(length <= min(lenbins), min(lenbins), length),
                      .by = year) %>%
    tidytable::summarise(n_s = mean(tot),
                         n_h = mean(tot),
                         length_tot = sum(value),
                         .by = c(year, length)) %>%
    tidytable::mutate(prop = length_tot / n_s) %>%
    tidytable::left_join(expand.grid(year = unique(.$year), length = lenbins), .)  %>%
    tidytable::replace_na(list(prop = 0)) %>%
    tidytable::mutate(SA_Index = 1,
                      n_s = mean(n_s, na.rm = T),
                      n_h = mean(n_h, na.rm = T), # for now, set hauls as sample size, don't have haul-level data
                      .by = year) %>%
    tidytable::select(-length_tot) %>%
    tidytable::pivot_wider(names_from = length, values_from = prop) -> flc_hist

  # get historical input sample size
  fsc_iss <- vroom::vroom(here::here(year, "data", "user_input", "goa_pop_fsc_iss.csv"), delim = ",")$iss

  # put 'em together
  flc_hist %>%
    tidytable::bind_rows(flc) %>%
    tidytable::mutate(n_s = fsc_iss,
                      n_h = fsc_iss) -> flc

  if(!is.null(id)) {
    vroom::vroom_write(flc, here::here(year, "data", "output", paste0(fishery, "_length_comp-", id, ".csv")), ",")
    flc
  } else if(isTRUE(save)){
    vroom::vroom_write(flc, here::here(year, "data", "output", paste0(fishery, "_length_comp.csv")), ",")
    flc
  }
  flc


}

#' trawl survey length composition analysis
#'
#' @param year assessment year
#' @param area survey area default = "goa"
#' @param lenbins lenbin file if left NULL it looks for (year/data/user_input/len_bins.csv")
#' @param bysex should the length comp be calculated by sex - default is null (not differentiated)
#' @param alt alternate folder to save to - will be placed in "year/alt/data" folder
#' @param save
#' @return
#' @export bts_length_comp
#'
#' @examples
#'
bts_length_comp <- function(year, area = "goa", lenbins = NULL, bysex = NULL, alt=NULL, save = TRUE){


  read.csv(here::here(year, "data", "raw", "bts_length_data.csv")) %>%
    dplyr::rename_with(tolower) %>%
    tidytable::summarize(n_s = sum(frequency),
                         n_h = length(unique(hauljoin)),
                         .by = c(year))-> df

  if(is.null(lenbins)){
    stop("Please provide a vector of length buns or the file that is in the user_input folder e.g.,('lengthbins.csv') with a column names 'len_bins'")
  }

  if(!is.vector(lenbins)){
    lenbins =  vroom::vroom(here::here(year, "data", "user_input", lenbins), delim = ",")$len_bins
  }

  vroom::vroom(here::here(year, "data", "raw", paste0(area, "_total_bts_sizecomp_data.csv"))) %>%
    dplyr::rename_with(tolower) %>%
    dplyr::filter(!is.na(length_mm)) %>%
    dplyr::mutate(length = length_mm / 10) %>%
    tidytable::select(-length_mm) -> dat

  if(!is.null(bysex)){
    # note that this code needs to still be changed to have column names consistent with gap_products
    df %>%
      dplyr::rename_with(tolower) %>%
      dplyr::filter(summary_depth < 995, year != 2001) %>%
      tidyr::pivot_longer(cols = c(males, females, unsexed)) %>%
      dplyr:: mutate(bin = round((length / 10 - 0.5) / 20, 1) * 20 + 1) %>%
      dplyr::filter(bin %in% lenbins) %>%
      dplyr::group_by(year, name, bin) %>%
      dplyr::summarise(value = sum(value)) %>%
      dplyr::ungroup() %>%
      tidyr::complete(bin, tidyr::nesting(year, name), fill = list(value = 0)) %>%
      dplyr::group_by(year, name) %>%
      dplyr::mutate(prop = value / sum(value)) %>%
      dplyr::select(-value) %>%
      dplyr::left_join(dat) %>%
      dplyr::group_by(year) %>%
      dplyr::mutate(SA_Index = 1,
                    n_s = mean(n_s, na.rm = T),
                    n_h = mean(n_h, na.rm = T)) %>%
      tidyr::pivot_wider(names_from = bin, values_from = prop) -> size_comp
  } else {
    dat %>%
      tidytable::mutate(length = ifelse(length >= max(lenbins), max(lenbins), length)) %>%
      tidytable::filter(length %in% lenbins) %>%
      tidytable::mutate(tot = sum(population_count), .by = c(year)) %>%
      tidytable::summarise(prop = sum(population_count) / mean(tot), .by = c(year, length)) %>%
      tidyr::replace_na(list(prop = 0)) %>%
      tidytable::left_join(df) %>%
      tidytable::mutate(SA_Index = 1) %>%
      tidyr::pivot_wider(names_from = length, values_from = prop) -> size_comp
  }

  if(!is.null(alt)) {
    write.csv(size_comp, here::here(year, alt, "data", paste0(area, "_ts_length_comp.csv")), row.names = F)
  } else if(isTRUE(save)){
    write.csv(size_comp, here::here(year, "data", "output", paste0(area, "_ts_length_comp.csv")), row.names = F)
  }
  size_comp

}

#' estimate allometric relationship and weight-at-age
#'
#' @param year model year
#' @param admb_home = location admb exists on your computer
#' @param rec_age recruitment age
#' @param area currently fixed at "goa"
#' @param alt alternate folder to save to - will be placed in "year/alt/data" folder
#' @param save save in the default location
#' @return
#' @export weight_at_age
#'
#' @examples weight_at_age(year = 2020, admb_home = "C:/Program Files (x86)/ADMB-12.1", rec_age = 2)
weight_at_age <- function(year, admb_home = NULL, rec_age, area = "goa", alt=NULL, save = TRUE){

  area = tolower(area)
  if(is.null(admb_home)){
    R2admb::setup_admb()
  } else {
    R2admb::setup_admb(admb_home)
  }

  if (!file.exists(here::here(year,"data", "output", "ae_model.csv"))){
    stop("You must first run the age-error function 'ageage()")
  } else {
    nages_m = nrow(vroom::vroom(here::here(year, "data", "output", "ae_model.csv")))
  }
  ages_m = rec_age:(rec_age + nages_m - 1)



  # data ----
  vroom::vroom(here::here(year, "data", "raw", "bts_length_data.csv")) %>%
    dplyr::rename_with(tolower) %>%
    dplyr::filter(year >= 1990, !is.na(length)) -> length_data_raw

  if(!("frequency" %in% colnames(length_data_raw))){
    length_data_raw %>%
      dplyr::select(age, length) %>%
      dplyr::group_by(age, length) %>%
      dplyr::summarise(frequency = dplyr::n()) -> length_data_raw
  }


  vroom::vroom(here::here(year, "data", "raw", "bts_specimen_data.csv")) %>%
    dplyr::rename_with(tolower) %>%
    dplyr::select(year, age, length, weight) %>%
    dplyr::filter(year >= 1990, !is.na(age))  %>%
    dplyr::select(-year) -> age_data_raw


  # Get parameters
  ages = sort(unique(age_data_raw$age))
  nages = length(ages)
  lengths = sort(unique(age_data_raw$length))
  nlengths = length(lengths)

  # Subset to ages with >1 obs
  n_a = table(age_data_raw$age)
  r = which(n_a<2)
  if(length(r)>0){
    n_a = n_a[-r]
  }
  ages = as.numeric(names(n_a))
  nages = length(ages)
  age_data_1 = NULL
  for(a in 1:nages){
    t = subset(age_data_raw, age_data_raw$age==ages[a])
    age_data_1 = rbind(age_data_1,t)
  }
  # Get Age-length key together
  n_al = table(age_data_1$age, age_data_1$length)
  n_l = colSums(n_al)
  r = which(n_l<2)

  if(length(r)>0){
    n_l = n_l[-r]
    n_al = n_al[,-r]
  }

  lengths = as.numeric(names(n_l))
  nlengths = length(lengths)
  N_l = matrix(nrow=nlengths)
  rownames(N_l) = lengths

  for(l in 1:nlengths){
    N_l[l,1] = sum(subset(length_data_raw$frequency, length_data_raw$length==lengths[l]))
  }

  N_al = matrix(0, nrow=nages, ncol=nlengths)
  rownames(N_al) = ages
  colnames(N_al) = lengths

  for(l in 1:nlengths){
    N_al[,l] = n_al[,l] / n_l[l] * N_l[l]
  }

  # Get mean weight and r age-length key together
  Wbar_la = r_la = V_Wbar_la = V_r_la = theta_la = matrix(NA,nrow=nages,ncol=nlengths)
  rownames(Wbar_la) = rownames(r_la) = rownames(V_Wbar_la) = rownames(V_r_la) = rownames(theta_la) = ages
  colnames(Wbar_la) = colnames(r_la) = colnames(V_Wbar_la) = colnames(V_r_la) = colnames(theta_la) = lengths

  theta_a = vector(length=nages)
  alpha_l = vector(length=nlengths)

  for(a in 1:nages){
    for(l in 1:nlengths){
      awl_data = subset(age_data_1,
                        age_data_1$age == ages[a] &
                          age_data_1$length == lengths[l])

      if(length(awl_data$weight) > 0){
        Wbar_la[a,l] = mean(awl_data$weight, na.rm=TRUE)

        if(length(awl_data$weight)>1){
          V_Wbar_la[a,l] = var(awl_data$weight, na.rm=TRUE) / length(awl_data$weight)
        }
      }

      alpha_l[l] = N_l[l] / sum(N_l)
      theta_la[a,l] = n_al[a,l] / sum(n_al[,l])
      r_la[a,l] = alpha_l[l] * theta_la[a,l]
    }
    theta_a[a] = sum(r_la[a,])
  }

  L = sum(N_l)
  A_l = colSums(n_al)

  for(a in 1:nages){
    for(l in 1:nlengths){
      V_r_la[a,l] = alpha_l[l]^2 * theta_la[a,l] *
        (1 - theta_la[a,l]) / (A_l[l] - 1) + alpha_l[l] *
        (theta_la[a,l] - theta_a[a])^2 / L
    }
  }

  # Get/compile weight-at-age statistics
  Age = ages
  SS = vector(length = nages)
  Wbar = vector(length = nages)
  SD_Wbar = vector(length = nages)
  for(a in 1:nages){
    SS[a] = length(subset(age_data_1$weight, age_data_1$age == ages[a]))
    Wbar[a] = sum(r_la[a,] * Wbar_la[a,], na.rm=TRUE) / sum(r_la[a,])

    SD_Wbar[a] = sqrt(sum(r_la[a,]^2 * V_Wbar_la[a,] +
                            (Wbar_la[a,] - Wbar[a])^2 * V_r_la[a,], na.rm=TRUE) /
                        theta_a[a]^2) *
      sqrt(length(subset(age_data_1$weight, age_data_1$age == ages[a])))
  }
  WaA_stats = as.data.frame(cbind(Age,SS,Wbar,SD_Wbar))
  r = which(WaA_stats$SD_Wbar == 0)
  WaA_stats = WaA_stats[-r,]
  r = which(WaA_stats$SS < 30)
  WaA_stats = WaA_stats[-r,]

  # Write data

  if(!is.null(alt)) {
    vroom::vroom_write(WaA_stats,
                       here::here(year, alt, "data", "waa_stats.csv"), ",")
  } else {
    vroom::vroom_write(WaA_stats,
                       here::here(year, "data", "output", "waa_stats.csv"), ",")
  }


  # Get/compile weight-at-length statistics
  age_data_raw %>%
    dplyr::select(length, weight) %>%
    dplyr::filter(length > 0, !is.na(weight)) -> lw_data

  lw_data %>%
    dplyr::group_by(length) %>%
    dplyr::summarise(Wbar = mean(weight, na.rm = T),
                     SD_Wbar = sd(weight, na.rm = T)) %>%
    tidyr::drop_na() -> lw_mdl_data


  # Write data
  if(!is.null(alt)) {
    vroom::vroom_write(lw_mdl_data, here::here(year, alt, "data", "wal_stats.csv"), ",")
  } else {
    vroom::vroom_write(lw_mdl_data,
                       here::here(year, "data", "output", "wal_stats.csv"), ",")
  }


  # Run allometric model ----
  DAT = c("# Data file for allometric model of mean weight by length",
          "# Number of lengths (nlengths)",
          nrow(lw_mdl_data),
          "# Lengths with observed mean weight (lengths)",
          paste(lw_mdl_data$length, collapse=" "),
          "# Observed mean weight (Wbar_obs)",
          paste(lw_mdl_data$Wbar, collapse=" "),
          "# SD in Observed mean weight (SD_Wbar)",
          paste(lw_mdl_data$SD_Wbar, collapse=" "))

  setwd(here::here(year, "data", "models", "allometric"))
  write.table(DAT, "allometric.dat", quote=FALSE,
              row.names=FALSE, col.names=FALSE)

  # run model

  R2admb::compile_admb("allometric", verbose = TRUE)
  R2admb::run_admb("allometric", verbose = TRUE)
  par = readLines("allometric.par", warn = FALSE)
  alpha_lw = as.numeric(strsplit(par[grep("alpha", par) + 1]," ")[[1]])
  beta_lw = as.numeric(strsplit(par[grep("beta", par) + 1]," ")[[1]])

  setwd(here::here())

  allo = data.frame(alpha_lw = alpha_lw, beta_lw = beta_lw)

  if(!is.null(alt)) {
    vroom::vroom_write(allo, here::here(year, alt, "data", "alpha_beta_lw.csv"), ",")
  } else {
    vroom::vroom_write(allo, here::here(year, "data", "output", "alpha_beta_lw.csv"), ",")
  }

  setwd(here::here(year, "data", "models", "wvonb"))

  # Run LVBmodel and estimate mean weight
  PIN <- c("# Parameter starting values for LVB model of mean weight",
           "# Winf", "800",
           "# k", "0.1",
           "# t0", "0",
           "# beta", as.character(beta_lw))

  write.table(PIN, "wVBL.PIN", quote=FALSE, row.names=FALSE, col.names=FALSE)

  WaA_stats = data.frame(WaA_stats)
  DAT<-c("# Data file for LVB model of mean weight",
         "# Number of ages (nages)",
         length(WaA_stats$Age),
         "# Ages with observed mean weight (ages)",
         paste(WaA_stats$Age, collapse=" "),
         "# Observed mean weight (Wbar_obs)",
         paste(WaA_stats$Wbar, collapse=" "),
         "# SD in Observed mean weight (Wbar_obs)",
         paste(WaA_stats$SD_Wbar, collapse=" "))

  write.table(DAT, file="wvbl.dat", quote=FALSE, row.names=FALSE, col.names=FALSE)

  # run model

  R2admb::compile_admb("wvbl", verbose = TRUE)
  R2admb::run_admb("wvbl", verbose = TRUE)

  REP <- readLines("wvbl.rep", warn=FALSE)

  setwd(here::here())

  Winf = as.numeric(strsplit(REP[grep("Winf", REP)[1]], " ")[[1]][2])
  k = as.numeric(strsplit(REP[grep("k", REP)[1]], " ")[[1]][2])
  t0 = as.numeric(strsplit(REP[grep("t0", REP)[1]], " ")[[1]][2])

  Wbar = Winf * (1 - exp(-k * (ages_m - t0)))^beta_lw
  Wbar[nages_m] = 0.5 * (Wbar[nages_m] + Winf)
  Wbar = round(Wbar, digits=1)
  Wbar_params = cbind(Winf, k, t0, beta_lw)

  if(!is.null(alt)) {
    write.csv(Wbar_params, here::here(year, alt, "data", "Wbar_params.csv"))
    write.csv(Wbar, here::here(year, alt, "data", "waa.csv"))
  } else {
    write.csv(Wbar_params, here::here(year, "data", "output", "Wbar_params.csv"))
    write.csv(Wbar, here::here(year, "data", "output", "waa.csv"))
  }

  Wbar
}


#' Concatenate a .dat file
#'
#' @param year assessment year
#' @param species "NORK", "REBS", "SABL"
#' @param area "goa", "bsai", "everywhere"
#' @param folder folder that the `.tpl` will be in
#' @param dat_name what to call the .dat file - ".dat" will be appended to the name
#' @param rec_age recruitment age
#' @param plus_age plus age group
#' @param spawn_mo spawning month
#' @param maturity if maturity is from outside the model (should be placed in user_input folder)
#' @param n_ageage number of age error transmission matrices default is 1
#' @param n_sizeage number of size at age transmission matrices default is 1
#' @param retro not yet implemented
#' @param n_fleets number of fishing fleets e.g., ll fleet, trawl fleet, default is 1
#' @param n_ts number of trawl surveys dafault is 1
#' @param n_lls number of longline surveys e.g., domestic and japanes, default is 1
#' @export concat_dat
#'
#' @examples concat_dat(year = 2020, species = "NORK",  area = "goa", folder = "base", dat_name = "goa_nr", rec_age = 2, plus_age = 45)
#'
concat_dat <- function(year, species, area = "goa", folder, dat_name, rec_age, plus_age, spawn_mo = 5,
                       maturity = NULL, alt = NULL, n_ageage = 1, n_sizeage = 1,
                       retro = NULL, n_fleets = 1, n_ts = NULL, n_lls = NULL){

  # create directory
  if (!dir.exists(here::here(year, folder))){
    dir.create(here::here(year, folder), recursive=TRUE)
  }

  if(is.null(alt)) {

    if(length(grep(paste0(area,"_lls"),
                   list.files(here::here(year, "data", "output")), value=TRUE)) > 0){
      llslc = read.csv(here::here(year, "data", "output", paste0(area, "_lls_length_comp.csv")))
      llsb = read.csv(here::here(year, "data", "output", paste0(area, "_lls_biomass.csv")))
    }

    if(!is.null(maturity)){
      mature = as.vector(read.csv(paste0(here::here(year, "data", "user_input", maturity))) %>%
                           dplyr::rename_with(tolower) %>%
                           dplyr::select(-age))
    }

    fishery = grep("fsh", list.files(here::here(year, 'data', "output")), value=TRUE)
    survey = grep("ts_", list.files(here::here(year, 'data', "output")), value=TRUE)
    ll_survey = grep("lls_", list.files(here::here(year, 'data', "output")), value=TRUE)

    catch = read.csv(here::here(year, "data", "output", grep("catch", fishery, value=TRUE)))
    waa = read.csv(here::here(year, "data", "output", "waa.csv"))
    saa = read.csv(here::here(year, "data", "output", "saa.csv"))
    ae = read.csv(here::here(year, "data", "output", "ae_model.csv"))
    fishac = read.csv(here::here(year, "data", "output", grep("age", fishery, value=TRUE)))
    fishlc = read.csv(here::here(year, "data", "output", grep("length", fishery, value=TRUE)))
    tsac = read.csv(here::here(year, "data", "output", grep("age", survey, value=TRUE)))
    tslc = read.csv(here::here(year, "data", "output", grep("length", survey, value=TRUE)))
    tsb = read.csv(here::here(year, "data", "output", grep("biomass", survey, value=TRUE)))

    if(length(ll_survey) > 0){
      llsrpw = read.csv(here::here(year, "data", "output", grep("biomass", ll_survey, value=TRUE)))
      llsslc = read.csv(here::here(year, "data", "output", grep("length", ll_survey, value=TRUE)))
      llsrpn = read.csv(here::here(year, "data", "output", grep("numbers", ll_survey, value=TRUE)))
    }
  } else {
    if(length(grep(paste0(area,"_lls"),
                   list.files(here::here(year, alt, "data")), value=TRUE)) > 0){
      llslc = read.csv(here::here(year, alt, "data", paste0(area, "_lls_length_comp.csv")))
      llsb = read.csv(here::here(year, alt, "data", paste0(area, "_lls_biomass.csv")))
    }

    if(!is.null(maturity)){
      mature = as.vector(read.csv(paste0(here::here(year, "data", "user_input", maturity))) %>%
                           dplyr::rename_with(tolower) %>%
                           dplyr::select(-age))
    }

    fishery = grep("fsh", list.files(here::here(year, alt, 'data')), value=TRUE)
    survey = grep("ts_", list.files(here::here(year, alt, 'data')), value=TRUE)
    ll_survey = grep("lls_", list.files(here::here(year, alt, 'data')), value=TRUE)

    catch = read.csv(here::here(year, alt, "data", grep("catch", fishery, value=TRUE)))
    waa = read.csv(here::here(year, alt, "data", "waa.csv"))
    saa = read.csv(here::here(year, alt, "data", "saa.csv"))
    ae = read.csv(here::here(year, alt, "data", "ae_model.csv"))
    fishac = read.csv(here::here(year, alt, "data", grep("age", fishery, value=TRUE)))
    fishlc = read.csv(here::here(year, alt, "data", grep("length", fishery, value=TRUE)))
    tsac = read.csv(here::here(year, alt, "data", grep("age", survey, value=TRUE)))
    tslc = read.csv(here::here(year, alt, "data", grep("length", survey, value=TRUE)))
    tsb = read.csv(here::here(year, alt, "data", grep("biomass", survey, value=TRUE)))

    if(length(ll_survey) > 0){
      llsrpw = read.csv(here::here(year, alt, "data", grep("biomass", ll_survey, value=TRUE)))
      llsslc = read.csv(here::here(year, alt, "data", grep("length", ll_survey, value=TRUE)))
      llsrpn = read.csv(here::here(year, alt, "data", grep("numbers", ll_survey, value=TRUE)))
    }
  }

  names(tsb) <- c("year", "biomass", "se", "lci", "uci")
  m_nages = nrow(ae)
  nages = length(rec_age:plus_age)

  # get length bin info
  lbin = as.numeric(gsub("[^0-9.]", "",  colnames(tslc)))
  lbin = lbin[!is.na(lbin)]
  nlenbins = length(lbin)

  if(is.null(n_ageage)){
    n_ageage = 1
  }

  if(is.null(n_sizeage)){
    n_sizeage = 1
  }

  sep = "# -------------------------------------------------------------------"

  # header ----
  header = c(sep,
             paste0("#", area, " ", species, " Rockfish .dat file for ADMB optimization"),
             paste ("# New data provided on:", read.table(file = here::here(year, "data/raw/data_called.txt"),
                                                          sep = "\t")[2,1]),
             "# Notes:",
             "#   ~ Weight-at-age and length-age transition matrix automatically updated",
             "#   ~ Formatted to conduct automated retrospective analysis",
             "#   ~ Does not use most recent years fishery size data",
             "#   ~ Does not use fishery size data in years when ages are expected",
             sep,
             "#",
             "#")

  # model inputs ----

  if(is.null(maturity)){
    mipv <- c(sep,
              "# Model input parameters/vectors",
              sep,
              "# Start and end years, recruitment age, number of age and length bins",
              "# Model start year (styr):",
              as.character(min(catch$Year)),
              "# Model end year (endyr): #!",
              as.character(year),
              "# Age at recruitment (rec_age): #-",
              as.character(rec_age),
              "# Number of ages in data (nages_D):",
              as.character(nages),
              "# Number of ages in model (nages_M):",
              as.character(m_nages),
              "# Number of length bins (nlenbins):",
              as.character(nlenbins),
              "# Number of age-age transition matrices (n_ageage_mat):",
              as.character(n_ageage),
              "# Number of size-age transition matrices (n_sizeage_mat):",
              as.character(n_sizeage),
              "# Length bin labels (len_bin_labels):",
              paste(lbin, collapse=" "),
              "# Spawn month (spawn_fract):",
              as.character(spawn_mo),
              "#",
              "#")

  } else {
    mipv <- c(sep,
              "# Model input parameters/vectors",
              sep,
              "# Start and end years, recruitment age, number of age and length bins",
              "# Model start year (styr):",
              as.character(min(catch$Year)),
              "# Model end year (endyr): #!",
              as.character(year),
              "# Age at recruitment (rec_age): #-",
              as.character(rec_age),
              "# Number of ages in data (nages_D):",
              as.character(nages),
              "# Number of ages in model (nages_M):",
              as.character(m_nages),
              "# Number of length bins (nlenbins):",
              as.character(nlenbins),
              "# Number of age-age transition matrices (n_ageage_mat):",
              as.character(n_ageage),
              "# Number of size-age transition matrices (n_sizeage_mat):",
              as.character(n_sizeage),
              "# Length bin labels (len_bin_labels):",
              paste(lbin, collapse=" "),
              "# Spawn month (spawn_fract):",
              as.character(spawn_mo),
              "#",
              "#")

    mat = c(sep,
            "# Proportion mature at age (p_mature):",
            paste0("#! ",
                   paste(mature$mature, collapse = " ")),
            "#-",
            "",
            "")
  }

  waa = c(sep,
          "# Weight-at-age (wt):",
          paste0("#! ",
                 paste(waa$x, collapse=" ")),
          "#-",
          "#",
          "#")

  # fishery catch ----
  fishery_catch = c(sep,
                    "# Fishery catch (mt): obs_catch(styr,endyr)",
                    sep,
                    paste0("#! ", paste(min(catch$Year):year, collapse=" ")),
                    paste(catch$Catch, collapse=" "),
                    "#-",
                    "",
                    "")
  # cpue ----
  # not currently used for northern rockfish
  cpue = c(sep,
           "# CPUE Data",
           sep,
           "# Number of CPUE years",
           "0",
           "# CPUE observations (leave blank if 0)",
           "",
           "")

  # trawl biomass ----
  trawl_biomass = c(sep,
                    "# Trawl Survey Biomass",
                    sep,
                    "#! Number of trawl surveys: nyrs_srv1",
                    as.character(nrow(tsb)),
                    "#- Trawl survey years: yrs_srv1(1,nyrs_srv1) #!",
                    paste(tsb$year, collapse=" "),
                    "#- Observed trawl survey biomass (mt): obs_srv1_biom(1,nyrs_srv1) #!",
                    paste(tsb$biomass, collapse=" "),
                    "#- SE of observed trawl survey biomass: obs_srv1_se(1,nyrs_srv1) #!",
                    paste(tsb$se, collapse=" "),
                    "#- Lower CI, 1.96*SE #!",
                    paste(tsb$lci, collapse=" "),
                    "#- Upper CI, 1.96*SE #!",
                    paste(tsb$uci, collapse=" "),
                    "#-",
                    "",
                    "")
  # long line survey biomass ----

  if(exists("llsrpw")){
    ll_biomass = c(
      sep,
      "# Longline Survey Biomass",
      sep,
      "# Number of longline surveys: nyrs_srv2",
      as.character(nrow(llsb)),
      "# Longline survey years: yrs_srv2(1,nyrs_srv2)",
      paste(llsb$year, collapse=" "),
      "# Observed longline survey biomass (mt): obs_srv2_biom(1,nyrs_srv2)",
      paste(llsb$rpw, collapse=" "),
      "# SE of observed longline survey biomass: obs_srv2_se(1,nyrs_srv2)",
      paste(llsb$sd, collapse=" "),
      "# Lower CI, 1.96*SE",
      paste(llsb$lci, collapse=" "),
      "# Upper CI, 1.96*SE",
      paste(llsb$uci, collapse=" "),
      "",
      "")
  } else {
    ll_biomass = c(
      sep,
      "# Longline Survey Biomass",
      sep,
      "# Number of longline surveys: nyrs_srv2",
      "1",
      "# Longline survey years: yrs_srv2(1,nyrs_srv2)",
      "1999",
      "# Observed longline survey biomass (mt): obs_srv2_biom(1,nyrs_srv2)",
      "1000",
      "# SE of observed longline survey biomass: obs_srv2_se(1,nyrs_srv2)",
      "100",
      "# Lower CI, 1.96*SE",
      "10",
      "# Upper CI, 1.96*SE",
      "10000",
      "",
      "")
  }

  # fishery age comp ----
  fac <- c(
    sep,
    "# Fishery Age Composition",
    sep,
    "#! Number of years: nyrs_fish_age",
    as.character(nrow(fishac)),
    "#- Fishery age comp years: yrs_fish_age #!",
    paste(fishac$year, collapse=" "),
    "#- Number of samples: nsamples_fish_age(1,nyrs_fish_age) #!",
    paste(fishac$n_s, collapse=" "),
    "#- Number of hauls: nhauls_fish_age(1,nyrs_fish_age) #!",
    paste(fishac$n_h, collapse=" "),
    "#- Index for age-age error matrix #!",
    paste(fishac$AA_Index, collapse=" "),
    "#- Observed fishery age compositions (proportions at age): oac_fish(1,nyrs_fish_age,1,nages) #!",
    collapse_row(dplyr::select(fishac, -year, -n_s, -n_h, -AA_Index)),
    "#-",
    "",
    "")

  # trawl survey age comp ----

  tsac <- c(sep,
            "# Trawl Survey Age Composition",
            sep,
            "#! Number of years: nyrs_srv1_age",
            as.character(nrow(tsac)),
            "#- Trawl Survey age comp years: yrs_srv1_age #!",
            paste(tsac$year, collapse=" "),
            "#- Number of samples: nsamples_srv1_age(1,nyrs_srv1_age) #!",
            paste(tsac$n_s, collapse=" "),
            "#- Number of hauls: nhauls_srv1_age(1,nyrs_srv1_age) #!",
            paste(tsac$n_h, collapse=" "),
            "#- Index for age-age error matrix #!",
            paste(tsac$AA_Index, collapse=" "),
            "#- Observed trawl survey age compositions (proportions at age): oac_srv1(1,nyrs_srv1_age,1,nages) #!",
            collapse_row(dplyr::select(tsac, -year, -n_s, -n_h, -AA_Index)),
            "#-",
            "",
            "")

  # fishery length comp ----
  flc <- c(
    sep,
    "# Fishery Size Composition",
    sep,
    "#! Number of years:",
    as.character(nrow(fishlc)),
    "#- Fishery size comp years: #!",
    paste(fishlc$year, collapse=" "),
    "#- Number of samples:  #!",
    paste(fishlc$n_s, collapse=" "),
    "#- Number of hauls:  #!",
    paste(fishlc$n_h, collapse=" "),
    "#- Index for size-age error matrix #!",
    paste(fishlc$SA_Index, collapse=" "),
    "#- Observed fishery size compositions (proportions at age)#!",
    collapse_row(dplyr::select(fishlc, -year, -n_s, -n_h, -SA_Index)),
    "#-",
    "",
    "")

  # trawl survey size comp ----
  tslc <- c(
    sep,
    "# Trawl Survey Size Composition",
    sep,
    "#! Number of years:",
    as.character(nrow(tslc)),
    "#- Survey Years: #!",
    paste(tslc$year, collapse=" "),
    "#- Number of samples:#!",
    paste(tslc$n_s, collapse=" "),
    "#- Number of hauls: #!",
    paste(tslc$n_h, collapse=" "),
    "#- Index for size-age error matrix #!",
    paste(tslc$SA_Index, collapse=" "),
    "#- Observed survey size compositions (proportions at age): oac_fish(1,nyrs_fish_age,1,nages) #!",
    collapse_row(dplyr::select(tslc, -year, -n_s, -n_h, -SA_Index)),
    "#-",
    "",
    "")

  # longline survey size comp ----
  if(exists("llslc")){

    llsc <- c(sep,
              "# Longline Survey Size Composition",
              sep,
              "# Number of years: nyrs_srv2_size",
              as.character(nrow(llslc)),
              "# Longline Survey size comp years: yrs_srv1_size",
              paste(llslc$year, collapse=" "),
              "# Number of samples: nsamples_srv2_size(1,nyrs_srv2_size)",
              paste(llslc$n_s, collapse=" "),
              "# Number of hauls: nhauls_srv2_size(1,nyrs_srv2_size)",
              paste(llslc$n_h, collapse=" "),
              "# Index for size-age error matrix",
              paste(llslc$SA_Index, collapse=" "),
              "# Observed longline survey size compositions (proportions at length): osc_srv2(1,nyrs_srv2_size,1,nlenbins)",
              collapse_row(dplyr::select(llslc, -year, -n_s, -n_h, -SA_Index)),
              "",
              "")
  } else {
    llsc <- c(sep,
              "# Longline Survey Size Composition, NOT USED IN MODEL, include one year of fake data",
              sep,
              "# Number of years: nyrs_srv2_size",
              "1",
              "# Longline Survey size comp years: yrs_srv1_size",
              "1999",
              "# Number of samples: nsamples_srv2_size(1,nyrs_srv2_size)",
              "99",
              "# Number of hauls: nhauls_srv2_size(1,nyrs_srv2_size)",
              "99",
              "# Index for size-age error matrix",
              "1",
              "# Observed longline survey size compositions (proportions at length): osc_srv2(1,nyrs_srv2_size,1,nlenbins)",
              paste(seq(1/nlenbins, 1/nlenbins, length.out=nlenbins), collapse=" "),
              "",
              "")
  }

  # size-age transition matrix ----
  sizeage <- c(sep,
               "# Size-age transition matrix: proportion at size given age: ",
               sep,
               collapse_row(dplyr::select(saa, -age)),
               "#",
               "",
               "")

  # age error matrix ----
  aa <- c(sep,
          "# age error transition matrix: ",
          sep,
          collapse_row(ae),
          "#",
          "",
          "")

  # eof ----
  eof <- c(sep,
           "# end of file marker",
           sep,
           "42",
           "#!")

  # Compile DAT file for ADMB ----

  if(is.null(maturity)){
    dat <- c(header,
             mipv,
             waa,
             fishery_catch,
             cpue,
             trawl_biomass,
             ll_biomass,
             fac,
             tsac,
             flc,
             tslc,
             llsc,
             sizeage,
             aa,
             eof)
  } else {
    dat <- c(header,
             mipv,
             mat,
             waa,
             fishery_catch,
             cpue,
             trawl_biomass,
             ll_biomass,
             fac,
             tsac,
             flc,
             tslc,
             llsc,
             sizeage,
             aa,
             eof)
  }


  write.table(dat, file = here::here(year, folder, paste0(dat_name, ".dat")) ,
              quote=FALSE, row.names=FALSE, col.names=FALSE)
}

#' Concatenate a .dat file for goa pop
#'
#' @param year assessment year
#' @param species "NORK", "REBS", "SABL"
#' @param area "goa", "bsai", "everywhere"
#' @param folder folder that the `.tpl` will be in
#' @param dat_name what to call the .dat file - ".dat" will be appended to the name
#' @param rec_age recruitment age
#' @param plus_age plus age group
#' @param spawn_mo spawning month
#' @param maturity if maturity is from outside the model (should be placed in user_input folder)
#' @param n_ageage number of age error transmission matrices default is 1
#' @param n_sizeage number of size at age transmission matrices default is 1
#' @param retro not yet implemented
#' @param n_fleets number of fishing fleets e.g., ll fleet, trawl fleet, default is 1
#' @param n_ts number of trawl surveys dafault is 1
#' @param n_lls number of longline surveys e.g., domestic and japanes, default is 1
#' @export concat_dat_pop
#'
#' @examples concat_dat(year = 2020, species = "NORK",  area = "goa", folder = "base", dat_name = "goa_nr", rec_age = 2, plus_age = 45)
#'
concat_dat_pop <- function(year, species, area = "goa", folder, dat_name, rec_age, plus_age,
                           spawn_mo = 5, maturity = NULL, alt = NULL, n_ageage = 1, n_sizeage = 2,
                           retro = NULL, n_fleets = 1, n_ts = NULL, n_lls = NULL){

  # create directory
  if (!dir.exists(here::here(year, folder))){
    dir.create(here::here(year, folder), recursive=TRUE)
  }

  # get maturity data if fit outside model
  if(!is.null(maturity)){
    mature = as.vector(read.csv(paste0(here::here(year, "data", "user_input", maturity))) %>%
                         dplyr::rename_with(tolower) %>%
                         dplyr::select(-age))
  }

  # read in data parts
  fishery = grep("fsh", list.files(here::here(year, 'data', "output")), value=TRUE)
  survey = grep("ts_", list.files(here::here(year, 'data', "output")), value=TRUE)

  catch = read.csv(here::here(year, "data", "output", grep("catch", fishery, value=TRUE)))
  waa = read.csv(here::here(year, "data", "output", "waa.csv"))
  saa = read.csv(here::here(year, "data", "output", "saa.csv"))
  saa60 = read.csv(here::here(year, "data", "output", "saa_60.csv"))
  ae = read.csv(here::here(year, "data", "output", "ae_model.csv"))
  fishac = read.csv(here::here(year, "data", "output", grep("age", fishery, value=TRUE)))
  fishlc = read.csv(here::here(year, "data", "output", grep("length", fishery, value=TRUE)))
  tsac = read.csv(here::here(year, "data", "output", grep("age", survey, value=TRUE)))
  tslc = read.csv(here::here(year, "data", "output", grep("length", survey, value=TRUE)))
  tsb = read.csv(here::here(year, "data", "output", grep("biomass", survey, value=TRUE)))

  names(tsb) <- c("year", "biomass", "se", "lci", "uci")
  m_nages = nrow(ae)
  nages = length(rec_age:plus_age)

  # get length bin info
  lbin = as.numeric(gsub("[^0-9.]", "",  colnames(tslc)))
  lbin = lbin[!is.na(lbin)]
  nlenbins = length(lbin)

  if(is.null(n_ageage)){
    n_ageage = 1
  }

  if(is.null(n_sizeage)){
    n_sizeage = 2
  }

  sep = "# -------------------------------------------------------------------"

  # header ----
  header = c(sep,
             paste0("# ",toupper(area), " Pacific ocean perch .dat file for ADMB optimization"),
             paste ("# New data provided on:", read.table(file = here::here(year, "data/raw/data_called.txt"),
                                                          sep = "\t")[2,1]),
             "# Notes:",
             "#   ~ Total catch prior to 1990 frozen",
             "#   ~ Total catch from 1991 on uses catch downloaded from AKFIN",
             "#   ~ Weight-at-age and current length-age transition matrix automatically updated",
             "#   ~ Formatted to conduct automated retrospective analysis",
             "#   ~ This data file provides ABL kluge for 2001 survey biomass for Eastern Gulf",
             "#   ~ Does not use most recent years fishery size data",
             "#   ~ Does not use fishery size data in years when ages are expected",
             sep,
             "",
             "")

  # model inputs ----

  if(is.null(maturity)){
    mipv <- c(sep,
              "# Model input parameters/vectors",
              sep,
              "# Start and end years, recruitment age, number of age and length bins",
              "# Model start year (styr):",
              as.character(min(catch$year)),
              "# Model end year (endyr): #!",
              as.character(year),
              "# Age at recruitment (rec_age): #-",
              as.character(rec_age),
              "# Number of ages in data (nages_D):",
              as.character(nages),
              "# Number of ages in model (nages_M):",
              as.character(m_nages),
              "# Number of length bins (nlenbins):",
              as.character(nlenbins),
              "# Number of age-age transition matrices (n_ageage_mat):",
              as.character(n_ageage),
              "# Number of size-age transition matrices (n_sizeage_mat):",
              as.character(n_sizeage),
              "# Length bin labels (len_bin_labels):",
              paste(lbin, collapse=" "),
              "# Spawn month (spawn_fract):",
              as.character(spawn_mo),
              "# Weight-at-age (wt):",
              paste(waa$x, collapse=" "),
              "",
              "")

  } else {
    mipv <- c(sep,
              "# Model input parameters/vectors",
              sep,
              "# Start and end years, recruitment age, number of age and length bins",
              "# Model start year (styr):",
              as.character(min(catch$Year)),
              "# Model end year (endyr): #!",
              as.character(year),
              "# Age at recruitment (rec_age): #-",
              as.character(rec_age),
              "# Number of ages in data (nages_D):",
              as.character(nages),
              "# Number of ages in model (nages_M):",
              as.character(m_nages),
              "# Number of length bins (nlenbins):",
              as.character(nlenbins),
              "# Number of age-age transition matrices (n_ageage_mat):",
              as.character(n_ageage),
              "# Number of size-age transition matrices (n_sizeage_mat):",
              as.character(n_sizeage),
              "# Length bin labels (len_bin_labels):",
              paste(lbin, collapse=" "),
              "# Spawn month (spawn_fract):",
              as.character(spawn_mo),
              "# Weight-at-age (wt):",
              paste(waa$x, collapse=" "),
              "# Proportion mature at age (p_mature):",
              paste(mature$mature, collapse = " "),
              "",
              "")

  }


  # fishery catch ----
  fishery_catch = c(sep,
                    "# Fishery catch (mt): obs_catch(styr,endyr)",
                    sep,
                    paste0("#! ", paste(min(catch$year):year, collapse=" ")),
                    paste(catch$catch, collapse=" "),
                    "#-",
                    "",
                    "")

  # cpue ----
  cpue = c(sep,
           "# CPUE Data, NOT USED IN MODEL",
           sep,
           "# Number of CPUE years",
           "0",
           "# CPUE observations (leave blank if 0)",
           "",
           "")

  # trawl biomass ----
  trawl_biomass = c(sep,
                    "# Trawl Survey Biomass",
                    sep,
                    "#! Number of trawl surveys: nyrs_srv1",
                    as.character(nrow(tsb)),
                    "#- Trawl survey years: yrs_srv1(1,nyrs_srv1) #!",
                    paste(tsb$year, collapse=" "),
                    "#- Observed trawl survey biomass (mt): obs_srv1_biom(1,nyrs_srv1) #!",
                    paste(tsb$biomass, collapse=" "),
                    "#- SE of observed trawl survey biomass: obs_srv1_se(1,nyrs_srv1) #!",
                    paste(tsb$se, collapse=" "),
                    "#- Lower CI, 1.96*SE #!",
                    paste(tsb$lci, collapse=" "),
                    "#- Upper CI, 1.96*SE #!",
                    paste(tsb$uci, collapse=" "),
                    "#-",
                    "",
                    "")

  # long line survey biomass ----
  ll_biomass = c(
    sep,
    "# Longline Survey Biomass, NOT USED IN MODEL",
    sep,
    "# Number of longline surveys: nyrs_srv2",
    "1",
    "# Longline survey years: yrs_srv2(1,nyrs_srv2)",
    "1999",
    "# Observed longline survey biomass (mt): obs_srv2_biom(1,nyrs_srv2)",
    "1000",
    "# SE of observed longline survey biomass: obs_srv2_se(1,nyrs_srv2)",
    "100",
    "# Lower CI, 1.96*SE",
    "10",
    "# Upper CI, 1.96*SE",
    "10000",
    "",
    "")

  # fishery age comp ----
  fac <- c(
    sep,
    "# Fishery Age Composition",
    sep,
    "#! Number of years: nyrs_fish_age",
    as.character(nrow(fishac)),
    "#- Fishery age comp years: yrs_fish_age #!",
    paste(fishac$year, collapse=" "),
    "#- Number of samples: nsamples_fish_age(1,nyrs_fish_age) #!",
    paste(fishac$n_s, collapse=" "),
    "#- Number of hauls: nhauls_fish_age(1,nyrs_fish_age) #!",
    paste(fishac$n_h, collapse=" "),
    "#- Index for age-age error matrix #!",
    paste(fishac$AA_Index, collapse=" "),
    "#- Observed fishery age compositions (proportions at age): oac_fish(1,nyrs_fish_age,1,nages) #!",
    collapse_row(dplyr::select(fishac, -year, -n_s, -n_h, -AA_Index)),
    "#-",
    "",
    "")

  # trawl survey age comp ----
  tsac <- c(sep,
            "# Trawl Survey Age Composition",
            sep,
            "#! Number of years: nyrs_srv1_age",
            as.character(nrow(tsac)),
            "#- Trawl Survey age comp years: yrs_srv1_age #!",
            paste(tsac$year, collapse=" "),
            "#- Number of samples: nsamples_srv1_age(1,nyrs_srv1_age) #!",
            paste(tsac$n_s, collapse=" "),
            "#- Number of hauls: nhauls_srv1_age(1,nyrs_srv1_age) #!",
            paste(tsac$n_h, collapse=" "),
            "#- Index for age-age error matrix #!",
            paste(tsac$AA_Index, collapse=" "),
            "#- Observed trawl survey age compositions (proportions at age): oac_srv1(1,nyrs_srv1_age,1,nages) #!",
            collapse_row(dplyr::select(tsac, -year, -n_s, -n_h, -AA_Index)),
            "#-",
            "",
            "")

  # fishery length comp ----
  flc <- c(
    sep,
    "# Fishery Size Composition",
    sep,
    "#! Number of years:",
    as.character(nrow(fishlc)),
    "#- Fishery size comp years: #!",
    paste(fishlc$year, collapse=" "),
    "#- Number of samples:  #!",
    paste(fishlc$n_s, collapse=" "),
    "#- Number of hauls:  #!",
    paste(fishlc$n_h, collapse=" "),
    "#- Index for size-age error matrix #!",
    paste(fishlc$SA_Index, collapse=" "),
    "#- Observed fishery size compositions (proportions at age)#!",
    collapse_row(dplyr::select(fishlc, -year, -n_s, -n_h, -SA_Index)),
    "#-",
    "",
    "")

  # trawl survey size comp ----
  tslc <- c(
    sep,
    "# Trawl Survey Size Composition, NOT USED IN MODEL",
    sep,
    "#! Number of years:",
    as.character(nrow(tslc)),
    "#- Survey Years: #!",
    paste(tslc$year, collapse=" "),
    "#- Number of samples:#!",
    paste(tslc$n_s, collapse=" "),
    "#- Number of hauls: #!",
    paste(tslc$n_h, collapse=" "),
    "#- Index for size-age error matrix #!",
    paste(tslc$SA_Index, collapse=" "),
    "#- Observed survey size compositions (proportions at age): oac_fish(1,nyrs_fish_age,1,nages) #!",
    collapse_row(dplyr::select(tslc, -year, -n_s, -n_h, -SA_Index)),
    "#-",
    "",
    "")

  # longline survey size comp ----
  if(exists("llslc")){

    llsc <- c(sep,
              "# Longline Survey Size Composition",
              sep,
              "# Number of years: nyrs_srv2_size",
              as.character(nrow(llslc)),
              "# Longline Survey size comp years: yrs_srv1_size",
              paste(llslc$year, collapse=" "),
              "# Number of samples: nsamples_srv2_size(1,nyrs_srv2_size)",
              paste(llslc$n_s, collapse=" "),
              "# Number of hauls: nhauls_srv2_size(1,nyrs_srv2_size)",
              paste(llslc$n_h, collapse=" "),
              "# Index for size-age error matrix",
              paste(llslc$SA_Index, collapse=" "),
              "# Observed longline survey size compositions (proportions at length): osc_srv2(1,nyrs_srv2_size,1,nlenbins)",
              collapse_row(dplyr::select(llslc, -year, -n_s, -n_h, -SA_Index)),
              "",
              "")
  } else {
    llsc <- c(sep,
              "# Longline Survey Size Composition, NOT USED IN MODEL, include one year of fake data",
              sep,
              "# Number of years: nyrs_srv2_size",
              "1",
              "# Longline Survey size comp years: yrs_srv1_size",
              "1999",
              "# Number of samples: nsamples_srv2_size(1,nyrs_srv2_size)",
              "99",
              "# Number of hauls: nhauls_srv2_size(1,nyrs_srv2_size)",
              "99",
              "# Index for size-age error matrix",
              "1",
              "# Observed longline survey size compositions (proportions at length): osc_srv2(1,nyrs_srv2_size,1,nlenbins)",
              paste(seq(1/nlenbins, 1/nlenbins, length.out=nlenbins), collapse=" "),
              "",
              "")
  }

  # size-age transition matrix ----

  sizeage <- c(sep,
               "# Size-age transition matrix (proportion at size given age): sizeage(1,nages,1,nlenbins)",
               sep,
               "# Size-age matrix representing slower growth in the 1960s by 6% in the 60's and 70's",
               collapse_row(dplyr::select(saa60, -age)),
               "",
               "# Size-age matrix using mean estimation for 1980s-2000s",
               collapse_row(dplyr::select(saa, -age)),
               "",
               "")

  # age error matrix ----
  aa <- c(sep,
          "# Ageing error matrix (proportion at reader age given tester age): ageage(1,nages_D,1,nages_M)",
          sep,
          collapse_row(ae),
          "",
          "")

  # eof ----
  eof <- c(sep,
           "# end of file marker",
           sep,
           "42",
           "#!")

  # Compile DAT file for ADMB ----

  if(is.null(maturity)){
    dat <- c(header,
             mipv,
             fishery_catch,
             cpue,
             trawl_biomass,
             ll_biomass,
             fac,
             tsac,
             flc,
             tslc,
             llsc,
             sizeage,
             aa,
             eof)
  } else {
    dat <- c(header,
             mipv,
             fishery_catch,
             cpue,
             trawl_biomass,
             ll_biomass,
             fac,
             tsac,
             flc,
             tslc,
             llsc,
             sizeage,
             aa,
             eof)
  }


  write.table(dat, file = here::here(year, folder, paste0(dat_name, "_", year, ".dat")) ,
              quote=FALSE, row.names=FALSE, col.names=FALSE)
}

#' write a .ctl file for goa pop
#'
#' @param year assessment year
#' @param base_mdl_fldr name of folder in which last full assessment model is stored (should be in 'models' folder)
#' @param curr_mdl_fldr name of folder that stores current model
#' @param mdl_name optional name for model run (default = "Model_1)
#' @param ctl_name what you want to name your ctl file
#' @param dat_name name of data file
#' @param folder where you want ctl file written
#' @param mcmc toggle to freeze sigr in mcmc estimation
#' @export write_ctl_pop
#'
#' @examples
#'
write_ctl_pop <- function(year, base_mdl_fldr, curr_mdl_fldr, mdl_name = "Model_1", ctl_name = "goa_pop", dat_name = "goa_pop", folder, mcmc = FALSE){

  ctl_orig = grep("ctl", list.files(here::here(year, 'mgmt', base_mdl_fldr)), value=TRUE)

  ctl_base = read.delim(here::here(year, 'mgmt', base_mdl_fldr, ctl_orig), sep = "", header = F)

  ctl_base[1,1] <- mdl_name # to change model name if desired
  ctl_base[2,1] <- paste0(dat_name, "_", year, ".dat")
  ctl_base[4,1] <- year
  ctl_base[length(ctl_base[,1]),1] <- read.csv(here::here(year, "data", "output", "yld_rat.csv"))$yld

  write.table(ctl_base, file = here::here(year, folder, paste0(dat_name, "_", year, ".ctl")) ,
              quote = FALSE, row.names = FALSE, col.names = FALSE)

  if(mcmc == TRUE){
    ctl_curr = grep("ctl", list.files(here::here(year, 'mgmt', curr_mdl_fldr)), value=TRUE)
    std_curr = grep("std", list.files(here::here(year, 'mgmt', curr_mdl_fldr)), value=TRUE)


    ctl_base = read.delim(here::here(year, 'mgmt', curr_mdl_fldr, ctl_curr), sep = "", header = F)
    std_base = read.delim(here::here(year, 'mgmt', curr_mdl_fldr, std_curr), sep = "", header = F)


    sigr <- as.numeric(std_base[which(std_base[,2] == "sigr"),3])

    ctl_base[which(ctl_base[,3] == "sigrprior"),1] <- sigr
    ctl_base[which(ctl_base[,3] == "ph_sigr"),1] <- -1

    write.table(ctl_base, file = here::here(year, folder, "mcmc", paste0(dat_name, "_", year, ".ctl")) ,
                quote = FALSE, row.names = FALSE, col.names = FALSE)
  }

}

