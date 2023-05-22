
#' clean up catch data
#'
#' @param year  year of assessment
#' @param species species of interest e.g., "SABL", "DUSK"
#' @param fishery identify the fishery default is "fsh"
#' @param TAC last three TAC in form: c(year-3, year-2, year-1)
#' @param discard if summarizing catch by discard/retained is desired change to TRUE
#' @param gear if summarizing catch by gear type is desired change to TRUE
#' @param fixed_catch if early catch is frozen place the file in user_input folder (format: year, catch)
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
clean_catch <- function(year, species, TAC = c(3333, 2222, 1111), discard = FALSE, gear = FALSE, fixed_catch = NULL, save = TRUE){
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
        dplyr::filter(variable == "catch")
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
  vroom::vroom(here::here(year, "data", "raw",  "fsh_obs_data.csv")) -> obs_data

  # Estimate catch ratio in final year to end of year
  obs_data %>%
    dplyr::filter(year %in% years) %>%
    dplyr::group_by(year) %>%
    dplyr::mutate(tot_catch = sum(extrapolated_weight),
                      test_date = lubridate::`year<-`(max(catch_data$week_end_date), year)) %>%
    dplyr::filter(haul_date <= test_date) %>%
    dplyr::summarise(oct_catch = round(sum(extrapolated_weight)),
                         tot_catch = round(mean(tot_catch))) %>%
    dplyr::ungroup() %>%
    dplyr::summarise(ratio = 1 + (sum(tot_catch) - sum(oct_catch)) / sum(oct_catch)) -> rat

  # Compute catch
  if(nrow(fc)>=1){
    catch_data %>%
      dplyr::select(year, catch = weight_posted) %>%
      dplyr::filter(year > max(fc$year)) %>%
      dplyr::group_by(year) %>%
      dplyr::summarise(catch = round(sum(catch), 4)) %>%
      dplyr::ungroup() %>%
      dplyr::bind_rows(fc) %>%
      dplyr::arrange(year) -> catch
  } else {
    catch_data %>%
      dplyr::select(year, catch = weight_posted) %>%
      dplyr::group_by(year) %>%
      dplyr::summarise(catch = round(sum(catch), 4)) %>%
      dplyr::ungroup() %>%
      dplyr::arrange(year) -> catch
  }


  # estimate yield ratio of previous 3 years relative to TAC
  catch %>%
    dplyr::filter(year %in% years) %>%
    dplyr::bind_cols(tac = TAC) %>%
    dplyr::summarise(yld = mean(catch / tac)) -> yield

  # estimate catch through end of the year
  catch %>%
    dplyr::filter(year==yr) %>%
    dplyr::mutate(proj_catch = catch * rat$ratio) %>%
    dplyr::bind_cols(rat, yield) -> yld

  if(isTRUE(save)){
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
#' @param by "depth", "stratum", "area", "total", "inpfc", "inpfc_depth" - only available for goa/ai (default: "total") - can only use a single switch
#' @param file if not using the design-based abundance, the file name must be stated (e.g. "GAP_VAST.csv")
#' @param rmv_yrs any survey years to exclude
#' @param save save to default location
#'
#' @return
#' @export bts_biomass
#'
#' @examples
#'
bts_biomass <- function(year, area = "goa", by = "total", file = NULL, rmv_yrs = NULL, save = TRUE){

  area = tolower(area)
  by = tolower(by)

  if(is.null(file)){

    vroom::vroom(here::here(year, "data", "raw", paste0(area, "_", by, "_bts_biomass_data.csv")) %>%
      dplyr::rename_all(tolower)  -> df

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
                             .by=year) %>%
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
    sb <- tidytable::filter(sb, !(year %in% rmv_yrs))
  }

  if(isTRUE(save)){
    vroom::vroom_write(sb, here::here(year, "data", "output", paste0(area, "_bts_biomass.csv")), delim=",")
  }

    sb


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
#' @param save = default is TRUE, FALSE outputs a list, names outputs to a data folder within a specific folder (e.g., save = "age_plus")
#'
#' @return
#' @export age_error
#'
#' @examples ageage(species = "NORK", year = 2020, admb_home = NULL, area = "GOA", rec_age = 2, plus_age = 45, max_age = 100)

age_error <- function(reader_tester, species, year, admb_home = NULL, area = "GOA", rec_age = 2, plus_age = 45, max_age = 100, save = TRUE){


  rt = vroom::vroom(here::here(year, "data", "user_input", reader_tester))


  if(is.null(admb_home)){
    R2admb::setup_admb()
  } else {
    R2admb::setup_admb(admb_home)
  }

  nages = length(rec_age:plus_age)

  rt %>%
    dplyr::filter(Species %in% sp_switch(species),
                  Region == area,
                  Read_Age > 0,
                  Test_Age > 0,
                  Final_Age > 0) %>%
    dplyr::group_by(Test_Age, Read_Age) %>%
    dplyr::rename_all(tolower) %>%
    dplyr::summarise(freq = dplyr::n()) -> dat


  dplyr::left_join(expand.grid(test_age = unique(dat$test_age),
                               read_age = unique(dat$read_age)),
                   dat) %>%
    tidyr::replace_na(list(freq = 0)) %>%
    dplyr::group_by(test_age) %>%
    dplyr::mutate(num = dplyr::case_when(test_age != read_age ~ freq)) %>%
    dplyr::summarise(num = sum(num, na.rm = TRUE),
                     den = sum(freq)) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(ape = 1 - (num / den),
                  ape = ifelse(is.nan(ape), 0, ape)) %>%
    dplyr::select(age = test_age, ape, ss = den) %>%
    dplyr::left_join(data.frame(age = min(dat$test_age):max(dat$test_age)), .) %>%
    tidyr::replace_na(list(ape = -9, ss = -9)) -> dats

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

  if(isFALSE(save) | is.null(save)) {
    list(mtx100 = mtx100, ae_sd = fits, ae_model = ae_Mdl)
  } else if(isTRUE(save)) {
    write.csv(mtx100, here::here(year, "data", "output", paste0("ae_mtx_", max_age, ".csv")), row.names = F)
    vroom::vroom_write(fits,  here::here(year,"data", "output", "ae_SD.csv"), ",")
    write.csv(ae_Mdl,  here::here(year, "data", "output", "ae_model.csv"), row.names = F)
    ae_Mdl
  } else {
    if(dir.exists(here::here(year, save, "data")) == FALSE){
      dir.create(here::here(year, save, "data"), recursive=TRUE)
    }
    write.csv(mtx100, here::here(year, save, "data", paste0("ae_mtx_", max_age, ".csv")), row.names = F)
    vroom::vroom_write(fits,  here::here(year,save, "data", "ae_SD.csv"), ",")
    write.csv(ae_Mdl,  here::here(year, save, "data", "ae_model.csv"), row.names = F)
    ae_Mdl
  }

}

#' size at age analysis
#'
#' @param year analysis year
#' @param admb_home location admb exists on your computer
#' @param rec_age recruitment age
#' @param max_age max age for age error analysis - default = 100
#' @param lenbins length bin file
#' @param save save in the default location

#' @return
#' @export size_at_age
#'
#' @examples
size_at_age <- function(year, admb_home = NULL, rec_age, lenbins = NULL, save = TRUE){

  if(!(isTRUE(save)) | !(isFALSE(save)) | !(is.null(save))) {

    if(!dir.exists(here::here(year, save, "data"))){
      dir.create(here::here(year, save, "data"), recursive=TRUE)
    }
  }

  if(is.null(admb_home)){
    R2admb::setup_admb()
  } else {
    R2admb::setup_admb(admb_home)
  }

  if(!(isTRUE(save)) | !(isFALSE(save)) | !(is.null(save))) {
        if (!file.exists(here::here(year, save, "data", "ae_model.csv"))){
          stop("You must first run the age-error function 'ageage()")
            } else {
              nages_m = nrow(read.csv(here::here(year, save, "data", "ae_model.csv")))
            }

    } else {
        if (!file.exists(here::here(year, "data", "output", "ae_model.csv"))){
          stop("You must first run the age-error function 'ageage()")
        } else {
          nages_m = nrow(read.csv(here::here(year, "data", "output", "ae_model.csv")))
      }
    }

    ages_m = rec_age:(rec_age + nages_m - 1)

  if(is.null(lenbins)){
    stop("Please provide the length bin file that is in the user_input folder e.g.,('lengthbins.csv')")
  } else {
    lenbins = read.csv(here::here(year, "data", "user_input", lenbins))$len_bins
  }

  vroom::vroom(here::here(year, "data", "raw", "goa_ts_saa_age_data.csv")) %>%
    dplyr::rename_all(tolower) %>%
    dplyr::select(year, age, length) %>%
    dplyr::filter(year>=1990, !is.na(age))  %>%
    dplyr::select(-year) %>%
    dplyr::group_by(age) %>%
    dplyr::filter(dplyr::n()>1) %>%
    dplyr::group_by(length) %>%
    dplyr::mutate(n_l = dplyr::n()) %>%
    dplyr::arrange(age, length) %>%
    dplyr::group_by(age) %>%
    dplyr::mutate(sample_size =  dplyr::n()) -> inter

  vroom::vroom(here::here(year, "data", "raw", "goa_ts_saa_length_data.csv")) %>%
    dplyr::rename_all(tolower) %>%
    dplyr::filter(year>=1990, !is.na(length)) -> dat

  if(is.null(dat$frequency)){
    dat %>%
      dplyr::group_by(length) %>%
      dplyr::summarise(tot = dplyr::n()) -> dat
  } else {
    dat %>%
      dplyr::select(frequency, length) %>%
      dplyr::group_by(length) %>%
      dplyr::summarise(tot = sum(frequency)) -> dat
  }

  dat %>%
    dplyr::left_join(inter, .) %>%
    dplyr::group_by(age, length) %>%
    dplyr::mutate(prop =  dplyr::n() / n_l * tot) %>%
    dplyr::distinct() %>%
    dplyr::group_by(age) %>%
    dplyr::summarise(sample_size = mean(sample_size),
                     Lbar = sum(prop * length) / sum(prop) * 0.1,
                     SD_Lbar = sqrt(1 / (sum(prop) - 1) * sum(prop * (length / 10 - Lbar)^2))) %>%
    dplyr::filter(SD_Lbar>=0.01) -> laa_stats

  if(!(isTRUE(save)) | !(isFALSE(save)) | !(is.null(save))){
    vroom::vroom_write(laa_stats, here::here(year, save, "data", "laa_stats.csv"))
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

  if(!(isTRUE(save)) | !(isFALSE(save)) | !(is.null(save))) {
      write.csv(params, here::here(year, save, "data", "lbar_params.csv"), row.names=FALSE)
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

  if(isFALSE(save) | is.null(save)) {
    saa
  } else if(isTRUE(save)) {
    readr::write_csv(saa, here::here(year, "data", "output", "saa.csv"))
    saa
  } else {
    if(dir.exists(here::here(year, save, "data")) == FALSE){
      dir.create(here::here(year, save, "data"), recursive=TRUE)
    }
    readr::write_csv(saa, here::here(year, save, "data", "saa.csv"))
    saa
  }
}

#' fishery age composition analysis
#'
#' @param year assessment year
#' @param fishery default is fsh1, change if age comps from multiple fisheries (e.g., fsh2)
#' @param rec_age recruitment age
#' @param plus_age plus age group
#' @param save whether to save the file
#'
#' @return
#' @export  fish_age_comp
#'
#' @examples
#' \dontrun{
#' fish_age_comp(year, fishery = "fsh1", rec_age, plus_age)
#' }
fish_age_comp <- function(year, fishery = "fsh", rec_age, plus_age, save = TRUE){

  vroom::vroom(here::here(year, "data", "raw", paste0(fishery, "_age_comp_data.csv")),
               col_types = list(HAUL_JOIN = "c",
                                PORT_JOIN = "c")) %>%
    dplyr::rename_all(tolower) %>%
    dplyr::filter(age>=rec_age) %>%
    dplyr::mutate(age = ifelse(age>plus_age, plus_age, age)) %>%
    dplyr::group_by(year) %>%
    dplyr::mutate(tot = dplyr::n()) %>%
    dplyr::filter(tot>49) %>%
    dplyr::mutate(n_h = length(unique(na.omit(haul_join))) +
                    length(unique(na.omit(port_join)))) %>%
    dplyr::group_by(year, age) %>%
    dplyr::summarise(n_s = mean(tot),
                     n_h = mean(n_h),
                     age_tot = dplyr::n()) %>%
    dplyr::mutate(prop = age_tot / n_s) %>%
    dplyr::left_join(expand.grid(year = unique(.$year),
                                 age = rec_age:plus_age), .) %>%
    tidyr::replace_na(list(prop = 0)) %>%
    dplyr::group_by(year) %>%
    dplyr::mutate(AA_Index = 1,
                  n_s = mean(n_s, na.rm = T),
                  n_h = mean(n_h, na.rm = T)) %>%
    dplyr::select(-age_tot) %>%
    tidyr::pivot_wider(names_from = age, values_from = prop) -> fac

  if(!(isTRUE(save)) | !(isFALSE(save)) | (!is.null(save))) {
    readr::write_csv(fac, here::here(year, alt, "data", paste0(fishery, "_age_comp.csv")))
    fac
  } else if(isTRUE(save)) {
    readr::write_csv(fac, here::here(year, "data", "output", paste0(fishery, "_age_comp.csv")))
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
#' @param save save in the default location
#'
#' @return
#' @export ts_age_comp
#'
#' @examples ts_age_comp(year = 2020, rec_age = 2, plus_age = 45)
ts_age_comp <- function(year, area = "goa", rec_age, plus_age, rmv_yrs = NULL, save = TRUE){

  read.csv(here::here(year, "data", "raw", paste0(area, "_ts_age_specimen_data.csv"))) %>%
    dplyr::rename_all(tolower) %>%
    dplyr::filter(!is.na(age)) %>%
    dplyr::group_by(year) %>%
    dplyr::summarise(n_s = dplyr::n(),
                     n_h = length(unique(hauljoin))) -> dat1


  read.csv(here::here(year, "data", "raw", paste0(area, "_ts_age_data.csv"))) %>%
    dplyr::rename_all(tolower) %>%
    tidytable::rename.(year = survey_year) %>%
    tidytable::filter.(age >= rec_age) %>%
    tidytable::mutate.(tot = sum(agepop),
                      age = ifelse(age < plus_age, age, plus_age),
                      .by = year) %>%
    tidytable::summarise.(prop = sum(agepop) / mean(tot),
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
    age_comp |>
      tidytable::filter.(!(year %in% rmv_yrs)) -> age_comp
  }

  if(!(isTRUE(save)) | !(isFALSE(save)) | (!is.null(save))) {
    readr::write_csv(age_comp, here::here(year, save, "data", paste0(area, "_ts_age_comp.csv")))
  } else if(isTRUE(save)){
      readr::write_csv(age_comp, here::here(year, "data", "output", paste0(area, "_ts_age_comp.csv")))
  }

  age_comp

}

#' fishery length composition analysis
#'
#' @param year assessment year
#' @param fishery default is "fsh1"
#' @param lenbins lenbin file if left NULL it looks for here::here(year, "data", "user_input", "len_bin_labels.csv")
#' @param rec_age recruitment age
#' @param save
#'
#' @return
#' @export fish_length_comp
#'
#' @examples
fish_length_comp <- function(year, fishery = "fsh", rec_age, lenbins = NULL, save = TRUE){

  if(is.null(lenbins)){
    stop("Please provide the length bin file that is in the user_input folder e.g.,('lengthbins.csv')")
  } else {
    lenbins =  vroom::vroom(here::here(year, "data", "user_input", lenbins), delim = ",")$len_bins
  }

  Y = year
  read.csv(here::here(year, "data", "raw", paste0(fishery, "_age_comp_data.csv"))) %>%
    dplyr::rename_all(tolower) %>%
    dplyr::filter(!is.na(age), age>=rec_age, year>1990 & year<Y) %>%
    dplyr::group_by(year) %>%
    dplyr::tally(name = "age") %>%
    dplyr::filter(age >= 50) %>%
    dplyr::ungroup() -> ages

  vroom::vroom(here::here(year, "data", "raw", paste0(fishery,"_length_comp_data.csv")),
               col_types = list(HAUL_JOIN = "c",
                                PORT_JOIN = "c")) %>%
    dplyr::rename_all(tolower) %>%
    dplyr::filter(!(year %in% unique(ages$year)), year>1990 & year<Y) %>%
    dplyr::group_by(year) %>%
    dplyr::mutate(tot = sum(frequency),
                  length = ifelse(length >= max(lenbins), max(lenbins), length),
                  n_h = length(unique(na.omit(haul_join))) + length(unique(na.omit(port_join)))) %>%
    dplyr::group_by(year, length) %>%
    dplyr::summarise(n_s = mean(tot),
                     n_h = mean(n_h),
                     length_tot = sum(frequency)) %>%
    dplyr::mutate(prop = length_tot / n_s) %>%
    dplyr::left_join(expand.grid(year = unique(.$year), length = lenbins), .) %>%
    tidyr::replace_na(list(prop = 0)) %>%
    dplyr::group_by(year) %>%
    dplyr::mutate(SA_Index = 1,
                  n_s = mean(n_s, na.rm = T),
                  n_h = mean(n_h, na.rm = T)) %>%
    dplyr::select(-length_tot) %>%
    tidyr::pivot_wider(names_from = length, values_from = prop) -> flc

  if(!(isTRUE(save)) | !(isFALSE(save)) | (!is.null(save))) {
    vroom::vroom_write(flc, here::here(year, save, "data", paste0(fishery, "_length_comp.csv")), ",")
  } else if(isTRUE(save)){
    vroom::vroom_write(flc, here::here(year, "data", "output", paste0(fishery, "_length_comp.csv")), ",")

  }
    flc


}

#' trawl survey length composition analysis
#'
#' @param year assessment year
#' @param area survey area default = "goa"
#' @param lenbins lenbin file if left NULL it looks for (year/data/user_input/len_bins.csv")
#' @param bysex should the length comp be calculated by sex - default is null (not differentiated)
#' @param save
#' @return
#' @export ts_length_comp
#'
#' @examples
#'
ts_length_comp <- function(year, area = "goa", lenbins = NULL, bysex = NULL, save = TRUE){


  read.csv(here::here(year, "data", "raw", paste0(area, "_ts_length_specimen_data.csv"))) %>%
    dplyr::rename_all(tolower) -> df

  if(is.null(lenbins)){
    stop("Please provide the length bin file that is in the user_input folder e.g.,('lengthbins.csv')")
  } else {
    lenbins =  vroom::vroom(here::here(year, "data", "user_input", lenbins), delim = ",")$len_bins
  }

  vroom::vroom(here::here(year, "data", "raw", paste0(area, "_ts_length_data.csv"))) %>%
    dplyr::rename_all(tolower) %>%
    dplyr::filter(!is.na(length)) %>%
    dplyr::mutate(length = length / 10) -> dat

  if("frequency" %in% colnames(dat)){
    dat %>%
      dplyr::group_by(year) %>%
      dplyr::summarise(n_s = sum(frequency),
                       n_h = length(unique(hauljoin))) %>%
      dplyr::ungroup() -> dat
  } else {
    dat %>%
      dplyr::group_by(year) %>%
      dplyr::summarise(n_s = dplyr::n(),
                       n_h = length(unique(hauljoin))) %>%
      dplyr::ungroup() -> dat
  }

  if(!is.null(bysex)){
    df %>%
      dplyr::rename_all(tolower) %>%
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
    df %>%
      dplyr::rename_all(tolower) %>%
      dplyr::mutate(length = length / 10,
                    length = ifelse(length >= max(lenbins), max(lenbins), length)) %>%
      dplyr::filter(length %in% lenbins) %>%
      dplyr::group_by(year) %>%
      dplyr::mutate(tot = sum(total)) %>%
      dplyr::group_by(year, length) %>%
      dplyr::summarise(prop = sum(total) / mean(tot)) %>%
      dplyr::ungroup() %>%
      dplyr::left_join(expand.grid(year = unique(.$year), length = lenbins), .) %>%
      tidyr::replace_na(list(prop = 0)) %>%
      dplyr::left_join(dat) %>%
      dplyr::group_by(year) %>%
      dplyr::mutate(SA_Index = 1,
                    n_s = mean(n_s, na.rm = T),
                    n_h = mean(n_h, na.rm = T)) %>%
      tidyr::pivot_wider(names_from = length, values_from = prop) -> size_comp
  }

  if(!(isTRUE(save)) | !(isFALSE(save)) | (!is.null(save))) {
    write.csv(size_comp, here::here(year, save, "data", paste0(area, "_ts_length_comp.csv")), row.names = F)
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
#' @param save save in the default location
#' @return
#' @export weight_at_age
#'
#' @examples weight_at_age(year = 2020, admb_home = "C:/Program Files (x86)/ADMB-12.1", rec_age = 2)
weight_at_age <- function(year, admb_home = NULL, rec_age, area = "goa", save = TRUE){

  if(is.null(admb_home)){
    R2admb::setup_admb()
  } else {
    R2admb::setup_admb(admb_home)
  }

  if (!file.exists(here::here(year,"data", "output", "ae_model.csv"))){
    stop("You must first run the age-error function 'ageage()")
  } else {
    if(!(isTRUE(save)) | !(isFALSE(save)) | (!is.null(save))) {
    nages_m = nrow(vroom::vroom(here::here(year, save, "data", "ae_model.csv")))
    } else {
      nages_m = nrow(vroom::vroom(here::here(year, "data", "output", "ae_model.csv")))
    }
    ages_m = rec_age:(rec_age + nages_m - 1)
  }


  # data ----
  vroom::vroom(here::here(year, "data", "raw", paste0(area, "_ts_saa_length_data.csv"))) %>%
    dplyr::rename_all(tolower) %>%
    dplyr::filter(year >= 1990, !is.na(length)) -> length_data_raw

  if(!("frequency" %in% colnames(length_data_raw))){
    length_data_raw %>%
      dplyr::select(age, length) %>%
      dplyr::group_by(age, length) %>%
      dplyr::summarise(frequency = dplyr::n()) -> length_data_raw
  }


  vroom::vroom(here::here(year, "data", "raw", paste0(area, "_ts_saa_age_data.csv"))) %>%
    dplyr::rename_all(tolower) %>%
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

  if(!(isTRUE(save)) | !(isFALSE(save)) | (!is.null(save))) {
    vroom::vroom_write(WaA_stats,
                       here::here(year, save, "data", "waa_stats.csv"), ",")
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
  if(!(isTRUE(save)) | !(isFALSE(save)) | (!is.null(save))) {
    vroom::vroom_write(lw_mdl_data, here::here(year, save, "data", "wal_stats.csv"), ",")
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

  if(!(isTRUE(save)) | !(isFALSE(save)) | (!is.null(save))) {
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

  if(!(isTRUE(save)) | !(isFALSE(save)) | (!is.null(save))) {
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
#' @param model folder that the `.tpl` will be in
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
#' @examples concat_dat(year = 2020, species = "NORK",  area = "goa", model = "base", dat_name = "goa_nr", rec_age = 2, plus_age = 45)
#'
concat_dat <- function(year, species, area = "goa", model, dat_name, rec_age, plus_age, spawn_mo = 5,
                       maturity = NULL, alt = NULL, n_ageage = 1, n_sizeage = 1,
                       retro = NULL, n_fleets = 1, n_ts = NULL, n_lls = NULL){

  # create directory
  if (!dir.exists(here::here(year, model))){
    dir.create(here::here(year, model), recursive=TRUE)
  }

  if(is.null(alt)) {

    if(length(grep(paste0(area,"_lls"),
                   list.files(here::here(year, "data", "output")), value=TRUE)) > 0){
      llslc = read.csv(here::here(year, "data", "output", paste0(area, "_lls_length_comp.csv")))
      llsb = read.csv(here::here(year, "data", "output", paste0(area, "_lls_biomass.csv")))
    }

    if(!is.null(maturity)){
      mature = as.vector(read.csv(paste0(here::here(year, "data", "user_input", maturity))) %>%
                           dplyr::rename_all(tolower) %>%
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
                           dplyr::rename_all(tolower) %>%
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


  write.table(dat, file = here::here(year, model, paste0(dat_name, ".dat")) ,
              quote=FALSE, row.names=FALSE, col.names=FALSE)
}

