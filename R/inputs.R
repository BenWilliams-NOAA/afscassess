run_model <- function(func, pars, data) {
  obj = RTMB::MakeADFun(func, pars, data = data, silent = TRUE)
  fit = nlminb(obj$par, obj$fn, obj$gr)
  list(
    report = obj$report(fit$par),
    fit = fit,
    sd = RTMB::sdreport(obj)
  )
}

#' get last 3 year's TAC values
#' @param year year of assessment
#' @param area "GOA" or "BSAI"
#' @export
get_tac <- function(year, area) {
  yr = year
  area = toupper(area)
  vroom::vroom(here::here(year, "data", "raw", "specs.csv")) %>%
    tidytable::filter(area_label == area) %>%
    tidytable::select(year, tac = total_allowable_catch) %>%
    tidytable::filter(year %in% (yr-3):(yr-1)) %>%
    tidytable::pull(tac)
}

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
  vroom::vroom(here::here(year, "data", "raw", "fish_catch_data.csv")) -> catch_data
  vroom::vroom(here::here(year, "data", "raw", "fish_obs_data.txt"),
               delim = ",",
               col_type = c(join_key = "c", haul_join = "c")) -> obs_data

  # Estimate catch ratio in final year to end of year
  obs_data %>%
    tidytable::filter(year %in% years) %>%
    tidytable::mutate(tot_catch = sum(extrapolated_weight, na.rm=T),
                      test_date = lubridate::`year<-`(max(catch_data$week_end_date), year), .by = year) %>%
    tidytable::filter(haul_date <= test_date) %>%
    tidytable::summarise(oct_catch = round(sum(extrapolated_weight, na.rm=T)),
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

  if(!(is.null(alt))) {
    vroom::vroom_write(catch, here::here(year, alt, "output",  "fish_catch.csv"), delim = ",")
    vroom::vroom_write(yld, here::here(year, alt, "output", "yld_rat.csv"), delim = ",")
    catch
  }else if(isTRUE(save)){
    vroom::vroom_write(catch, here::here(year, "data", "output",  "fish_catch.csv"), delim = ",")
    vroom::vroom_write(yld, here::here(year, "data", "output", "yld_rat.csv"), delim = ",")
    catch
  } else {
    list(catch = catch, yld_rat = yld)
  }

}

#' survey biomass
#'
#' @param year model year
#' @param area options are bs, bsslope, nbs, ai, goa, old_bs - can only call a single area
#' @param type "depth", "stratum", "area", "total", "inpfc", "inpfc_depth" - only available for goa/ai (default: "total") - can only use a single switch
#' @param file if not using the design-based abundance, the file name must be stated (e.g. "GAP_VAST.csv") which is stored in "data/user_input"
#' @param rmv_yrs any survey years to exclude
#' @param save save to default location, default: TRUE
#' @param id identifier will be appended to file name - "bts_biomass_vast, default:NULL
#'

#' @export bts_biomass
#'
#'
bts_biomass <- function(year,area = "goa", type = "total", file = NULL, rmv_yrs = NULL, save = TRUE, id = NULL){

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

    if(!is.null(id) & isTRUE(save)){
    vroom::vroom_write(sb, here::here(year, "data", "output", paste0(area, "_", type, "_bts_biomass_", id,".csv")), ",")
    } else if(is.null(id) & isTRUE(save)){
      vroom::vroom_write(sb, here::here(year, "data", "output", paste(area, type, "bts_biomass.csv", sep="_")), ",")
      sb
    } else
    sb

}

#' GAP survey biomass
#'
#' @param year model year
#' @param area options are bs, bsslope, nbs, ai, goa, old_bs - can only call a single area
#' @param type "depth", "stratum", "area", "total", "inpfc", "inpfc_depth" - only available for goa/ai (default: "total") - can only use a single switch
#' @param file if not using the design-based abundance, the file name must be stated (e.g. "GAP_VAST.csv") which is stored in "data/user_input"
#' @param rmv_yrs any survey years to exclude
#' @param save save to default location, default: TRUE
#' @param id identifier will be appended to file name - "bts_biomass_vast, default:NULL
#'

#' @export
#'
#'
bts_gap_biomass <- function(year, area = "goa", type = "region", file = NULL, rmv_yrs = NULL, save = TRUE, id = NULL){

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
        tidytable::summarise(biomass = sum(biomass_mt),
                             se = sqrt(sum(biomass_var)),
                             .by = year) %>%
        tidytable::mutate(lci = biomass - 1.96 * se,
                          uci = biomass + 1.96 * se,
                          lci = ifelse(lci < 0, 0, lci)) %>%
        tidytable::mutate(tidytable::across(tidytable::where(is.double), round)) %>%
        tidytable::filter(biomass > 0) -> sb
    }

  } else {
    if(grepl(".rds", file, ignore.case = TRUE)) {
      readRDS(here::here(year, "data", "user_input", file)) %>%
        tidytable::summarise(biomass = est / 1000,
                             se = sqrt((exp(se^2) - 1) * exp(2 * log_est + se^2)) / 1000,
                             lci = lwr / 1000,
                             uci = upr / 1000, .by = year) -> sb
    } else {
    vroom::vroom(here::here(year, "data", "user_input", file)) -> sb
    }
  }

  if(!is.null(rmv_yrs)){
    sb <- dplyr::filter(sb, !(year %in% rmv_yrs))
  }

  if(!is.null(id) & isTRUE(save)){
    vroom::vroom_write(sb, here::here(year, "data", "output", paste0(area, "_", type, "_bts_biomass_", id,".csv")), ",")
  } else if(is.null(id) & isTRUE(save)){
    vroom::vroom_write(sb, here::here(year, "data", "output", paste(area, type, "bts_biomass.csv", sep="_")), ",")
    sb
  } else
    sb

}
#' aging error analysis
#'
#' @param read_tester = looks for a file in the user_input folder, e.g., :"reader_tester.csv"
#' @param species = "NORK"
#' @param year = year of the assessment
#' @param area = "GOA" (BSAI not currently setup)
#' @param rec_age = recruitment age
#' @param plus_age = max age for modeling
#' @param max_age = max age for age error analysis, default 100
#' @param alt alternate folder to save to - will be placed in "year/alt/data" folder
#' @param save = default is TRUE, FALSE outputs a list, names outputs to a data folder within a specific folder (e.g., save = "age_plus")
#'
#' @export
#'
#' @examples ageage(species = "NORK", year = 2020, admb_home = NULL, area = "GOA", rec_age = 2, plus_age = 45, max_age = 100)

age_error <- function(reader_tester, species, year, area = "goa", rec_age = 2, plus_age = 45, max_age = 100, alt = NULL, save = TRUE){
  f_age_error <- function(pars, data) {
    RTMB::getAll(pars, data)
    s1 = exp(log_sigma1)
    s2 = exp(log_sigma2)

    # linear interpolation of sd across age range
    sigma_a = s1 + (age - min_age) * ((s2 - s1) / (max_age - min_age))

    # percent agreement likelihood
    p_corr  = RTMB::pnorm(0.5, 0, sigma_a) - RTMB::pnorm(-0.5, 0, sigma_a)
    p_corr1 = RTMB::pnorm(-0.5, 0, sigma_a) - RTMB::pnorm(-1.5, 0, sigma_a)
    p_corr2 = RTMB::pnorm(-1.5, 0, sigma_a) - RTMB::pnorm(-2.5, 0, sigma_a)

    phat = p_corr^2 + 2 * p_corr1^2 + 2 * p_corr2^2

    RTMB::REPORT(s1); RTMB::REPORT(s2); RTMB::REPORT(sigma_a)
    return(sum(sqrt(ss) * (phat - ape)^2))
  }
  target_species = sp_switch(species)
  # data prep
  rt = vroom::vroom(here::here(year, "data", "user_input", reader_tester))

  rt %>%
    tidytable::rename_with(tolower) %>%
    tidytable::filter(species %in% target_species,
                      region == toupper(area),
                      read_age > 0, test_age > 0, final_age > 0) %>%
    tidytable::summarise(freq = .N, .by = c(test_age, read_age)) %>%
    # Use tidytable's complete logic to fill missing age combinations
    tidytable::complete(test_age = full_seq(test_age, 1),
                        read_age = full_seq(read_age, 1),
                        fill = list(freq = 0)) %>%
    tidytable::summarise(num = sum(freq[test_age != read_age]),
                         den = sum(freq), .by = test_age) %>%
    tidytable::mutate(ape = tidytable::replace_na(1 - (num / den), 0)) %>%
    tidytable::select(age = test_age, ape, ss = den) %>%
    tidytable::filter(ss > 0) -> dats

  # reference original data range to maintain the correct slope
  min_dat_age = min(dats$age)
  max_dat_age = max(dats$age)

  # RTMB estimation
  dat_age = list(age = dats$age, ape = dats$ape, ss = dats$ss, min_age = min_dat_age, max_age = max_dat_age)
  pars = list(log_sigma1 = log(0.5), log_sigma2 = log(5.0))

  res = RTMButils::run_model(f_age_error,
                  pars = pars,
                  data = dat_age,
                  proj = FALSE)$rpt

  # results
  s1 <- res$s1
  s2 <- res$s2


  # calculate sd's for all ages
  tidytable::tidytable(age = rec_age:max_age) %>%
    tidytable::mutate(sd = s1 + (age - min_dat_age) * ((s2 - s1) / (max_dat_age - min_dat_age))) -> fits

  # outer() to build the probability matrix without a loop
  # rows = true age (j), cols = observed age (i)
  mtx100 = outer(fits$age, fits$age, function(true_a, obs_a) {
    # match SDs to the true age
    s = fits$sd[match(true_a, fits$age)]
    # standard interval probability calculation
    pnorm(obs_a + 0.5, true_a, s) - pnorm(obs_a - 0.5, true_a, s)
  })

  # Correct the edges (0 to first bin and last bin to Infinity)
  mtx100[, 1] = pnorm(fits$age[1] + 0.5, fits$age, fits$sd)
  mtx100[, ncol(mtx100)] = 1 - rowSums(mtx100[, 1:(ncol(mtx100)-1)])

  colnames(mtx100) = rownames(mtx100) = fits$age

  # collapse to assessment dimensions
  n_model_ages = length(rec_age:plus_age)
  ae_model = mtx100[, 1:n_model_ages]
  ae_model[, n_model_ages] = rowSums(mtx100[, n_model_ages:ncol(mtx100)])

  # trim matrix where the plus group is fully absorbed
  ae_model = ae_model[1:which(ae_model[, n_model_ages] >= 0.999)[1], ]

  result = list(mtx100 = mtx100, ae_sd = fits, ae_model = round(ae_model, 4))

  if(isTRUE(save) | !is.null(alt)) {
    path = if(!is.null(alt)) here::here(year, alt, "data") else here::here(year, "data", "output")
    if(!dir.exists(path)) dir.create(path, recursive = TRUE)
    vroom::vroom_write(as.data.frame(mtx100), file.path(path, "ae_mtx.csv"), ",")
    vroom::vroom_write(fits, file.path(path, "ae_sd.csv"), ",")
    vroom::vroom_write(as.data.frame(ae_model), file.path(path, "ae_model.csv"), ",")
  }

  return(result)
}


prep_alw_data <- function(age_data, length_data, model_ages, len_bins, rec_age) {

   # 1. Setup metadata
  max_age_err <- rec_age + model_ages - 1
  age_range <- rec_age:max_age_err

  # 2. Clean and summarize length frequency (dat)
  len_clean <- length_data %>%
    tidytable::rename_with(tolower) %>%
    tidytable::filter(year >= 1990, !is.na(length))

  if(!("frequency" %in% names(len_clean))){
    dat <- len_clean %>% tidytable::summarise(tot = .N, .by = length)
  } else {
    dat <- len_clean %>% tidytable::summarise(tot = sum(frequency), .by = length)
  }

  L_total <- sum(dat$tot)

  # 3. Size at Age (SAA) 
  laa_stats <- age_data %>%
    tidytable::rename_with(tolower) %>%
    tidytable::filter(year >= 1990, !is.na(age)) %>%
    tidytable::select(age, length) %>%
    tidytable::filter(.N > 1, .by = age) %>%
    tidytable::mutate(n_l = .N, .by = length) %>%
    tidytable::mutate(sample_size = .N, .by = age) %>%
    tidytable::left_join(dat, by = "length") %>%
    tidytable::mutate(prop = .N / n_l * first(tot), .by = c(age, length)) %>%
    tidytable::distinct(age, length, .keep_all = TRUE) %>%
    tidytable::mutate(lbar = sum(prop * length) / sum(prop) * 0.1, .by = age) %>%
    tidytable::summarise(
      sample_size = mean(sample_size),
      lbar = first(lbar),
      sd = sqrt(1 / (sum(prop) - 1) * sum(prop * (length / 10 - lbar)^2)),
      .by = age
    ) %>%
    tidytable::filter(sd >= 0.01)

  # 4. Length-Weight (LW) -
  lw_data <- age_data %>%
    tidytable::rename_with(tolower) %>%
    tidytable::filter(year >= 1990, length > 0, !is.na(weight)) %>%
    tidytable::summarise(wbar = mean(weight, na.rm = TRUE),
                         sd = sd(weight, na.rm = TRUE), .by = length) %>%
    tidytable::drop_na()

  # 5. weight at age 
  
  age_step1 <- age_data %>%
    tidytable::rename_with(tolower) %>%
    tidytable::filter(year >= 1990, !is.na(age)) %>%
    tidytable::mutate(SS = .N, .by = age) %>% 
    tidytable::filter(SS >= 2)
  
  valid_lengths <- age_step1 %>%
    tidytable::summarise(n_l_aged = .N, .by = length) %>%
    tidytable::filter(n_l_aged >= 2) %>%
    tidytable::pull(length)
  
  dat_legacy <- dat %>% tidytable::filter(length %in% valid_lengths)
  L_legacy <- sum(dat_legacy$tot)
  
  waa_stats <- age_step1 %>%
    tidytable::filter(length %in% valid_lengths) %>%
    tidytable::mutate(n_al = .N, .by = c(age, length)) %>%
    tidytable::mutate(n_l = .N, .by = length) %>%
    tidytable::left_join(dat_legacy %>%
                           tidytable::mutate(alpha_l = tot / L_legacy), by = "length") %>%
    tidytable::summarise(
      wbar_la = mean(weight, na.rm = TRUE),
      v_wbar_la = var(weight, na.rm = TRUE) / .N,
      theta_la = .N / first(n_l),
      r_la = (first(tot) / L_legacy) * (.N / first(n_l)),
      SS = first(SS), 
      alpha_l = first(alpha_l),
      n_l = first(n_l),
      .by = c(age, length)
    ) %>%
    tidytable::mutate(
      theta_a = sum(r_la), 
      wbar = sum(r_la * wbar_la, na.rm = TRUE) / sum(r_la), 
      v_r_la = alpha_l^2 * theta_la * (1 - theta_la) / (n_l - 1) + alpha_l * (theta_la - theta_a)^2 / L_legacy,
      .by = age 
    ) %>%
    tidytable::summarise(
      wbar = first(wbar),
      sd = sqrt(sum(r_la^2 * v_wbar_la + (wbar_la - wbar)^2 * v_r_la, na.rm = TRUE) / (first(theta_a)^2)) * sqrt(first(SS)),
      sample_size = first(SS),
      .by = age 
    ) %>%
    tidytable::filter(!is.na(sd), sd > 0, sample_size >= 30) %>%
    tidytable::arrange(age)

  list(dat_saa = list(age = laa_stats$age,
                      n = laa_stats$sample_size,
                      lbar = laa_stats$lbar,
                      sd = laa_stats$sd),
       dat_lw = list(length = lw_data$length,
                     wbar = lw_data$wbar,
                     sd = lw_data$sd),
       dat_waa = list(age = waa_stats$age,
                      wbar = waa_stats$wbar,
                      sd = waa_stats$sd),
       laa_stats = laa_stats,
       lw_data = lw_data,
       waa_stats = waa_stats)
}

#' size at age analysis
#'
#' @param year analysis year
#' @param age_data gap survey specimen data
#' @param length_data gap survey length data
#' @param lenbins length bins used e.g., 15:45
#' @param rec_age recruitment age
#' @param alt alternate folder to save to - will be placed in "year/alt/data" folder
#' @param save save in the default location


#' @export
#'

saa_waa <- function(year, age_data, length_data, len_bins, rec_age, alt=NULL, save = TRUE) {
  cmb <- function(f, d) function(p) f(p, d)
  if (!file.exists(here::here(year, "data", "output", "ae_model.csv"))){
    stop("You must first run the age-error function 'ageage()")
  } else {
    nages_m = nrow(read.csv(here::here(year, "data", "output", "ae_model.csv")))
  }
  ages = rec_age:(rec_age + nages_m - 1)

  prep = prep_alw_data(age_data, length_data, nages_m, len_bins, rec_age)

  # model: saa
  f_saa_mean <- function(pars, data) {
    RTMB::getAll(pars, data)
    linf = exp(log_linf)
    k = exp(log_k)
    pred = linf * (1 - exp(-k * (age - t0)))
    yvar = log(1 + sd^2 / lbar^2)
    rss = sum(0.5 * (log(2 * pi * yvar * lbar^2) + (log(pred) - log(lbar))^2 / yvar))
    
    RTMB::REPORT(linf); RTMB::REPORT(k); RTMB::REPORT(t0)
    return(rss)
  }

 f_saa_sd <- function(pars, data) {
    RTMB::getAll(pars, data)
    
    # 1. standard linear space (no exp() bounds) to match admb
    lpred_sd = alpha * log(age) + beta 
    
    # 2. ADMB Likelihood
    lnll = sum(sqrt(n) * (log(lpred_sd) - log(sd))^2) 
    
    RTMB::REPORT(alpha); RTMB::REPORT(beta)
    return(lnll)
  }

  rep_saa_sd = RTMButils::run_model(model = f_saa_sd,
                                    # Use exact ADMB initialization values!
                                    pars = list(alpha = 0.005, beta = 0.5), 
                                    data = prep$dat_saa, proj = FALSE)$rpt

  # rep_saa_sd = RTMButils::run_model(model = f_saa_sd,
  #                                     pars = list(log_alpha = log(5), beta = -2), 
  #                                     data = prep$dat_saa, proj = FALSE)$rpt
  
 
  rep_saa_mean = RTMButils::run_model(model = f_saa_mean,
                                 pars = list(log_linf = log(50), log_k = log(0.2), t0 = -0.1),
                                 data = prep$dat_saa, proj = FALSE)$rpt

  # Combine them into a single list so the rest of your script still works!
  rep_saa = list(linf = rep_saa_mean$linf, 
                 k = rep_saa_mean$k, 
                 t0 = rep_saa_mean$t0, 
                 alpha = rep_saa_sd$alpha, 
                 beta = rep_saa_sd$beta)

  # model: length-weight
  f_lw <- function(pars, data) {
    RTMB::getAll(pars, data)
    alpha = exp(log_alpha)
    beta = exp(log_beta)
    pred = alpha * length^beta
    yvar = log(1 + sd^2 / wbar^2)
    rss = sum(0.5 * (log(2 * pi * yvar * wbar^2) + (log(pred) - log(wbar))^2 / yvar))
    RTMB::REPORT(alpha); RTMB::REPORT(beta)
    return(rss)
  }

  rep_lw = RTMButils::run_model(f_lw,
                                pars = list(log_alpha=log(0.01),
                                            log_beta=log(3)),
                                data = prep$dat_lw,
                                proj = FALSE)$rpt

  # model: waa
  f_waa <- function(pars, data) {
    RTMB::getAll(pars, data)
    winf = exp(log_winf)
    k = exp(log_k)
    beta = exp(log_beta)
    pred = winf * (1 - exp(-k * (age - t0)))^beta
    yvar = log(1 + sd^2 / wbar^2)
    rss = sum(0.5 * (log(2 * pi * yvar * wbar^2) + (log(pred) - log(wbar))^2 / yvar))
    RTMB::REPORT(winf); RTMB::REPORT(k); RTMB::REPORT(t0); RTMB::REPORT(beta)
    return(rss)
  }

  
  rep_waa = RTMButils::run_model(f_waa,
                                 pars = list(log_winf=log(max(prep$dat_waa$wbar)), log_k=log(0.1), t0=0, log_beta=log(rep_lw$beta)),
                                 data = prep$dat_waa,
                                 map = list(log_beta = factor(NA)),
                                 proj = FALSE)$rpt

  # final matrices

  saa_matrix = tidytable::expand_grid(age = ages, length = len_bins) %>%
    tidytable::mutate(
      Lbar = rep_saa$linf * (1 - exp(-rep_saa$k * (age - rep_saa$t0))),
      Lbar = ifelse(age == max(age), 0.5 * (Lbar + rep_saa$linf), Lbar),
      sd_L = rep_saa$alpha * log(age) + rep_saa$beta,
      prob = ifelse(length == min(length),
                    pnorm(length + 0.5, Lbar, sd_L),
                    pnorm(length + 0.5, Lbar, sd_L) - pnorm(length - 0.5, Lbar, sd_L))
    ) %>%
    tidytable::select(age, length, prob) %>%
    tidytable::pivot_wider(names_from = length, values_from = prob) 
  
  # deal with plus group
  nc = ncol(saa_matrix)
  last_col_name = names(saa_matrix)[nc]
  
  # coerce to data.frame explicitly for the column slice so it behaves predictably
  saa_matrix[[last_col_name]] = 1 - rowSums(as.data.frame(saa_matrix)[, 2:(nc - 1)])
  
  # apply the rounding matrix-wide
  saa_matrix = saa_matrix %>% 
    tidytable::mutate(dplyr::across(2:nc, ~ round(.x, 4)))

  waa_table = tidytable::tidytable(age = ages) %>%
    tidytable::mutate(
      wbar = rep_waa$winf * (1 - exp(-rep_waa$k * (age - rep_waa$t0)))^rep_waa$beta,
      wbar = ifelse(age == max(age), 0.5 * (wbar + rep_waa$winf), wbar),
      wbar = round(wbar, 1)
    )

  if(isTRUE(save) | !is.null(alt)) {
    path = if(!is.null(alt)) here::here(year, alt, "data") else here::here(year, "data", "output")
    if(!dir.exists(path)) dir.create(path, recursive = TRUE)
    vroom::vroom_write(as.data.frame(saa_matrix), file.path(path, "saa.csv"), ",")
    vroom::vroom_write(waa_table, file.path(path, "waa.csv"), ",")
    vroom::vroom_write(as.data.frame(rep_saa), file.path(path, "params_saa.csv"), ",")
    vroom::vroom_write(as.data.frame(rep_waa), file.path(path, "params_waa.csv"), ",")
  }

  return(list(saa = saa_matrix, waa = waa_table, params_saa = rep_saa, params_waa = rep_waa))
}

#' fishery age composition analysis
#'
#' @param year assessment year
#' @param fishery default is fsh, change if age comps from multiple fisheries (e.g., fsh1, fsh2)
#' @param rec_age recruitment age
#' @param plus_age plus age group
#' @param rmv_yrs any years to remove form the age comp e.g. c(1987, 1989)
#' @param id id a specific comp name - will be placed at end of file name e.g., id='use' will create 'fsh_age_comp-use.csv' in the data/output folder
#' @param save whether to save the file - wll be placed in "year/data/output" folder
#'

#' @export  fish_age_comp
#'
#' @examples
#' \dontrun{
#' fish_age_comp(year, fishery = "fsh", rec_age, plus_age)
#' }
fish_age_comp <- function(year, fishery = "fish", rec_age, plus_age, rmv_yrs = NULL, id=NULL, save = TRUE){

  vroom::vroom(here::here(year, "data", "raw", paste0(fishery, "_specimen_data.txt")),
               delim = ",",
               col_type = c(join_key="c", haul_join="c", port_join="c")) %>%
    tidytable::filter(age>=rec_age, !(year %in% rmv_yrs), !is.na(length), specimen_type==1) %>%
    tidytable::mutate(age = ifelse(age>plus_age, plus_age, age)) %>%
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

  if(!is.null(id)) {
    vroom::vroom_write(fac, here::here(year, "data", 'output', paste0(fishery, "_age_comp-", id, ".csv")), ",")
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
#' @param id id a specific comp name - will be placed at end of file name e.g., id='use' will create 'bts_age_comp-use.csv' in the data/output folder
#' @param save save in the default location
#'

#' @export bts_age_comp
#'
#' @examples bts_age_comp(year = 2020, rec_age = 2, plus_age = 45)
bts_age_comp <- function(year, area = "goa", rec_age, plus_age, rmv_yrs = NULL, id=NULL, save = TRUE){

  area = tolower(area)
  read.csv(here::here(year, "data", "raw", paste0(area, "_bts_specimen_data.csv"))) %>%
    dplyr::filter(!is.na(age)) %>%
    dplyr::group_by(year) %>%
    dplyr::summarise(n_s = dplyr::n(),
                     n_h = length(unique(hauljoin))) -> dat1


  read.csv(here::here(year, "data", "raw", paste0(area, "_bts_agecomp_data.csv"))) %>%
    dplyr::rename_with(tolower) %>%
    tidytable::rename(year = survey_year) %>%
    tidytable::filter(age >= rec_age) %>%
    tidytable::mutate(tot = sum(agepop),
                      age = ifelse(age < plus_age, age, plus_age),
                      .by = year) %>%
    tidytable::summarise(prop = sum(agepop) / mean(tot),
                         .by = c(age, year)) %>%
    tidytable::left_join(dat1) %>%
    tidytable::left_join(expand.grid(year = unique(.$year),
                                     age = rec_age:plus_age), .) %>%
    tidytable::replace_na(list(prop = 0)) %>%
    tidytable::mutate(AA_Index = 1,
                      n_s = mean(n_s, na.rm = T),
                      n_h = mean(n_h, na.rm = T),
                      .by = year) %>%
    tidytable::pivot_wider(names_from = age, values_from = prop) %>%
    tidytable::arrange(year) -> age_comp


  if(!is.null(rmv_yrs)){
    age_comp  %>%
      tidytable::filter(!(year %in% rmv_yrs)) -> age_comp
  }

  if(!is.null(id)) {
    vroom::vroom_write(age_comp, here::here(year, "data", "output", paste0(area, "_bts_age_comp", id, ".csv")), ",")
    age_comp
  } else if(isTRUE(save)){
    vroom::vroom_write(age_comp, here::here(year, "data", "output", paste0(area, "_bts_age_comp.csv")), ",")
    age_comp
  }

  age_comp

}

#' GAP trawl survey age comp analysis
#'
#' @param year assessment year
#' @param area default is "goa"
#' @param rec_age recruitment age
#' @param plus_age plus group age
#' @param rmv_yrs any survey years to exclude
#' @param id id a specific comp name - will be placed at end of file name e.g., id='use' will create 'bts_age_comp-use.csv' in the data/output folder
#' @param save save in the default location
#'

#' @export
#'
#' @examples bts_gap_age_comp(year = 2020, rec_age = 2, plus_age = 45)
bts_gap_age_comp <- function(year, area = "goa", rec_age, plus_age, rmv_yrs = NULL, id=NULL, save = TRUE){

    area = tolower(area)
    # read.csv(here::here(year, "data", "raw", paste0(area, "_bts_gap_specimen_data.csv"))) %>%
    #   dplyr::filter(!is.na(age)) %>%
    #   dplyr::group_by(year) %>%
    #   dplyr::summarise(n_s = dplyr::n(),
    #                    n_h = length(unique(hauljoin))) -> dat1


    read.csv(here::here(year, "data", "raw", paste0(area, "_bts_gap_agecomp_data.csv"))) %>%
      dplyr::filter(age >= rec_age) %>%
      dplyr::mutate(tot = sum(population_count),
                        age = ifelse(age < plus_age, age, plus_age),
                        .by = year) %>%
      dplyr::summarise(prop = sum(population_count) / mean(tot),
                           .by = c(age, year)) %>%
      # tidytable::left_join(dat1) %>%
      dplyr::left_join(expand.grid(year = unique(.$year),
                                       age = rec_age:plus_age), .) %>%
      tidyr::replace_na(list(prop = 0)) %>%
      # tidytable::mutate(AA_Index = 1,
      #                   n_s = mean(n_s, na.rm = T),
      #                   n_h = mean(n_h, na.rm = T),
      #                   .by = year) %>%
      dplyr::pivot_wider(names_from = age, values_from = prop) %>%
      dplyr::arrange(year) -> age_comp


    if(!is.null(rmv_yrs)){
      age_comp  %>%
        tidytable::filter(!(year %in% rmv_yrs)) -> age_comp
    }

    if(!is.null(id)) {
      vroom::vroom_write(age_comp, here::here(year, "data", "output", paste0(area, "_bts_age_comp", id, ".csv")), ",")
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
#' @param fishery default is "fsh1"
#' @param lenbins lenbin file if left NULL it looks for here::here(year, "data", "user_input", "len_bin_labels.csv")
#' @param rec_age recruitment age
#' @param rmv_yrs any survey years to exclude
#' @param id id a specific comp name - will be placed at end of file name e.g., id='use' will create 'fsh_length_comp-use.csv' in the data/output folder
#' @param save default TRUE
#'

#' @export fish_length_comp
#'

fish_length_comp <- function(year, fishery = "fish", rec_age, lenbins = NULL, rmv_yrs = NULL, id=NULL, save = TRUE){

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
    tidytable::filter(!is.na(age), age>=rec_age, year<yr) %>%
    tidytable::group_by(year) %>%
    tidytable::tally(name = "age") %>%
    tidytable::filter(age >= 50) %>%
    tidytable::ungroup() -> ages

  vroom::vroom(here::here(year, "data", "raw", paste0(fishery,"_length_data.txt")),
               delim = "\t",
               col_type = c(haul_join="c", port_join="c")) %>%
    tidytable::filter(!(year %in% c(unique(ages$year), yr))) %>%
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

#' trawl survey length composition analysis
#'
#' @param year assessment year
#' @param area survey area default = "goa"
#' @param lenbins lenbin file if left NULL it looks for (year/data/user_input/len_bins.csv")
#' @param bysex should the length comp be calculated by sex - default is null (not differentiated)
#' @param rmv_yrs any survey years to exclude
#' @param alt alternate folder to save to - will be placed in "year/alt/data" folder
#' @param save default TRUE

#' @export bts_length_comp
#'
#'
bts_length_comp <- function(year, area = "goa", lenbins = NULL, bysex = NULL, rmv_yrs = NULL, alt=NULL, save = TRUE){


  area = tolower(area)
  read.csv(here::here(year, "data", "raw", paste0(area, "_bts_sizecomp_data.csv"))) %>%
    dplyr::rename_with(tolower) -> df

  if(is.null(lenbins)){
    stop("Please provide a vector of length buns or the file that is in the user_input folder e.g.,('lengthbins.csv') with a column names 'len_bins'")
  }

  if(!is.vector(lenbins)){
    lenbins =  vroom::vroom(here::here(year, "data", "user_input", lenbins), delim = ",")$len_bins
  }

  vroom::vroom(here::here(year, "data", "raw", paste0(area, "_bts_length_data.csv"))) %>%
    dplyr::rename_with(tolower) %>%
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
    df %>%
      dplyr::rename_with(tolower) %>%
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

  if(!is.null(rmv_yrs)){
    size_comp  %>%
      tidytable::filter(!(year %in% rmv_yrs)) -> size_comp
  }

  if(!is.null(alt)) {
    write.csv(size_comp, here::here(year, alt, "data", paste0(area, "_bts_sizecomp.csv")), row.names = F)
  } else if(isTRUE(save)){
    write.csv(size_comp, here::here(year, "data", "output", paste0(area, "_bts_sizecomp.csv")), row.names = F)
  }
    size_comp

}

#' GAP trawl survey length composition analysis
#'
#' @param year assessment year
#' @param area survey area default = "goa"
#' @param lenbins lenbin file if left NULL it looks for (year/data/user_input/len_bins.csv")
#' @param bysex should the length comp be calculated by sex - default is null (not differentiated)
#' @param rmv_yrs any survey years to exclude
#' @param alt alternate folder to save to - will be placed in "year/alt/data" folder
#' @param save default TRUE

#' @export
#'
#'
bts_gap_length_comp <- function(year, area = "goa", lenbins = NULL, bysex = NULL, rmv_yrs = NULL, alt=NULL, save = TRUE){


  area = tolower(area)
  read.csv(here::here(year, "data", "raw", paste0(area, "_bts_gap_sizecomp_data.csv"))) %>%
    dplyr::rename_with(tolower)  %>%
    dplyr::rename(length = length_mm, total = population_count) -> df

  if(is.null(lenbins)){
    stop("Please provide a vector of length buns or the file that is in the user_input folder e.g.,('lengthbins.csv') with a column names 'len_bins'")
  }

  if(!is.vector(lenbins)){
    lenbins =  vroom::vroom(here::here(year, "data", "user_input", lenbins), delim = ",")$len_bins
  }

  vroom::vroom(here::here(year, "data", "raw", paste0(area, "_bts_gap_length_data.csv"))) %>%
    dplyr::rename_with(tolower) %>%
    dplyr::filter(!is.na(length)) %>%
    dplyr::mutate(length = length / 10) -> dat

  if("frequency" %in% colnames(dat)){
    dat %>%
      tidytable::summarise(n_s = sum(frequency),
                       n_h = length(unique(hauljoin)),
                      .by = year) -> dat
  } else {
    dat %>%
      dplyr::group_by(year) %>%
      dplyr::summarise(n_s = dplyr::n(),
                       n_h = length(unique(hauljoin))) %>%
      dplyr::ungroup() -> dat
  }

  if(!is.null(bysex)){
    df %>%
      dplyr::rename_with(tolower) %>%
      tidytable::filter(summary_depth < 995, year != 2001) %>%
      tidyr::pivot_longer(cols = c(males, females, unsexed)) %>%
      tidytable:: mutate(bin = round((length / 10 - 0.5) / 20, 1) * 20 + 1) %>%
      tidytable::filter(bin %in% lenbins) %>%
      tidytable::group_by(year, name, bin) %>%
      tidytable::summarise(value = sum(value)) %>%
      tidytable::ungroup() %>%
      tidyr::complete(bin, tidyr::nesting(year, name), fill = list(value = 0)) %>%
      tidytable::group_by(year, name) %>%
      tidytable::mutate(prop = value / sum(value)) %>%
      dptidytablelyr::select(-value) %>%
      tidytable::left_join(dat) %>%
      tidytable::group_by(year) %>%
      tidytable::mutate(SA_Index = 1,
                    n_s = mean(n_s, na.rm = T),
                    n_h = mean(n_h, na.rm = T)) %>%
      tidyr::pivot_wider(names_from = bin, values_from = prop) -> size_comp
  } else {
    df %>%
      tidytable::rename_with(tolower) %>%
      tidytable::mutate(length = length / 10,
                        length = ifelse(length >= max(lenbins), max(lenbins), length)) %>%
      tidytable::filter(length %in% lenbins) %>%
      tidytable::mutate(tot = sum(total), .by = year) %>%
      tidytable::summarise(prop = sum(total) / mean(tot),
                        .by = c(year, length)) %>%
      tidytable::left_join(expand.grid(year = unique(.$year), length = lenbins), .) %>%
      tidytable::replace_na(list(prop = 0)) %>%
      tidytable::left_join(dat) %>%
      tidytable::mutate(SA_Index = 1,
                    n_s = mean(n_s, na.rm = T),
                    n_h = mean(n_h, na.rm = T), .by = year) %>%
      tidyr::pivot_wider(names_from = length, values_from = prop) -> size_comp
  }

  if(!is.null(rmv_yrs)){
    size_comp  %>%
      tidytable::filter(!(year %in% rmv_yrs)) -> size_comp
  }

  if(!is.null(alt)) {
    write.csv(size_comp, here::here(year, alt, "data", paste0(area, "_bts_sizecomp.csv")), row.names = F)
  } else if(isTRUE(save)){
    write.csv(size_comp, here::here(year, "data", "output", paste0(area, "_bts_sizecomp.csv")), row.names = F)
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
  vroom::vroom(here::here(year, "data", "raw", paste0(area, "_bts_length_data.csv"))) %>%
    dplyr::rename_with(tolower) %>%
    dplyr::filter(year >= 1990, !is.na(length)) -> length_data_raw

  if(!("frequency" %in% colnames(length_data_raw))){
    length_data_raw %>%
      dplyr::select(age, length) %>%
      dplyr::group_by(age, length) %>%
      dplyr::summarise(frequency = dplyr::n()) -> length_data_raw
  }


  vroom::vroom(here::here(year, "data", "raw", paste0(area, "_bts_specimen_data.csv"))) %>%
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

