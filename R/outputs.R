#' Process model results for tables and figs
#'
#' @param year  assessment year
#' @param model   model being evaluated (folder name)
#' @param model_name   name of the model e.g., updated_nr
#' @param dat_name name of dat file e.g., goa_nr_2020
#' @param rec_age recruitment age
#' @param plus_age plus age group
#' @param mcmc number of mcmcs run
#' @param mcsave the number of mcmcs saved
#' @param len_bins file name (stored in the "user_input" folder) that has length bins
#' @param ... future functions
#'
#' @return
#' @export process_results
#'
#' @examples process_results (year = 2020, model = m18.2, model_name = "goa_nr", dat_name = "goa_nr_2020", rec_age = 2, plus_age = 45, mcmc = 1e+07, mcsave = 2000, len_bins = "lbins.csv")
#'
process_results <- function(year, model, model_name, dat_name,
                            rec_age, plus_age, mcmc, mcsave, len_bins, ...){

  # setup
  if (!dir.exists(here::here(year, folder, "processed"))){
    dir.create(here::here(year, folder, "processed"), recursive=TRUE)
  }

  if (!dir.exists(here::here(year, folder, "figs"))){
    dir.create(here::here(year, folder, "figs"), recursive=TRUE)
  }

  if (!dir.exists(here::here(year, folder, "tables"))){
    dir.create(here::here(year, folder, "tables"), recursive=TRUE)
  }

  # helper functions
  rep_item <- function(name){
    t <- strsplit(REP[grep(name, REP)]," ")
    t <- subset(t[[1]], t[[1]]!="")
    if(t[[1]][1] == "TWL"){
      as.numeric(t[3:length(t)])
    } else {
      as.numeric(t[2:length(t)])
    }
  }


  # read in rep and ctl files
  REP <- readLines(here::here(year, folder, paste0(model_name, ".rep")))
  CTL <- readLines(here::here(year, folder, paste0(dat_name, ".ctl")))
  PSV <- file(here::here(year, folder, paste0(model_name, ".psv")), "rb")
  STD <- read.delim(here::here(year, model, paste0(model_name, ".std")), sep="", header = TRUE)
  mceval <- read.delim(here::here(year, folder, "evalout.prj"), sep="", header=FALSE)

  # clean rep file
  suppressWarnings(data.frame(year = unlist(base::strsplit(REP[grep("Year", REP)[1]]," "))) %>%
                     tidytable::mutate.(year = as.numeric(year)) %>%
                     tidytable::drop_na.() %>%
                     tidytable::pull.(year)) -> yrs

  suppressWarnings(data.frame(age = unlist(base::strsplit(REP[grep("Age", REP)[1]]," "))) %>%
                     tidytable::mutate.(age = as.numeric(age)) %>%
                     tidytable::drop_na.() %>%
                     tidytable::pull.(age)) -> ages

  styr_rec <- yrs[1] - length(ages) + rec_age

  suppressWarnings(as.data.frame(cbind(yrs = yrs, ages = ages, styr_rec = styr_rec)) %>%
                     tidytable::mutate.(ages = replace(ages, duplicated(ages), NA),
                                        styr_rec = replace(styr_rec, duplicated(styr_rec), NA))) %>%
    write.csv(here::here(year, folder, "processed", "ages_yrs.csv"), row.names = FALSE)

  # MCMC parameters ----

  npar = readBin(PSV, what = integer(), n=1)
  mcmcs = readBin(PSV, what = numeric(), n = (npar * mcmc / mcsave))
  close(PSV)
  mcmc_params = matrix(mcmcs, byrow=TRUE, ncol=npar)
  # thin the string
  mcmc_params = mcmc_params[501:nrow(mcmc_params),]
  colnames(mcmc_params) = STD$name[1:ncol(mcmc_params)]
  write.csv(mcmc_params, here::here(year, folder, "processed", "mcmc.csv"), row.names = FALSE)

  # mceval phase output ----

  #Curry's Change
  mceval = mceval[501:nrow(mceval),]

  #Length colnames = 286
  # columns mcmc_other = 271

  #1-8: Through objective function value

  colnames(mceval) = c("sigr", "q_srv1", "q_srv2", "F40", "natmort", "spawn_biom_proj",
                       "ABC", "obj_fun",
                       paste0("tot_biom_", yrs),
                       paste0("log_rec_dev_", seq(styr_rec, yrs[length(yrs)])),
                       paste0("spawn_biom_", yrs),
                       "log_mean_rec",
                       paste0("spawn_biom_proj_", max(yrs) + 1:15),
                       paste0("pred_catch_proj_", max(yrs) + 1:15),
                       paste0("rec_proj_", max(yrs) + 1:10),
                       paste0("tot_biom_proj_", max(yrs)))
  write.csv(mceval, here::here(year, folder, "processed", "mceval.csv"), row.names = FALSE)

  # catch data ----

  pred = base::strsplit(REP[grep("Pred_Catch", REP)], " ")
  r1 = which(pred[[1]] == "Pred_Catch")
  r2 = which(pred[[1]] == "Pred_catch_later")
  r3 = which(pred[[1]] == "")
  pred = as.numeric(pred[[1]][-c(r1, r2, r3)])

  obs = base::strsplit(REP[grep("Obs_Catch", REP)], " ")
  r1 = which(obs[[1]] == "Obs_Catch")
  r2 = which(obs[[1]] == "Obs_Catch_Later")
  r3 = which(obs[[1]] == "")
  obs = as.numeric(obs[[1]][-c(r1, r2, r3)])

  data.frame(year = yrs, obs = obs, pred = pred) %>%
    write.csv(here::here(year, model, "processed", "catch.csv"))



  # survey data ----
  syr = REP[grep("Survey Biomass",REP)[1]:(grep("Survey Biomass",REP)[2]-2)][2]
  syr = base::strsplit(syr," ")
  syr = subset(syr[[1]], syr[[1]]!="")
  syr = as.numeric(syr[2:length(syr)])

  obs = REP[grep("Survey Biomass",REP)[1]:(grep("Survey Biomass",REP)[2]-2)][4]
  obs = base::strsplit(obs," ")
  obs = subset(obs[[1]], obs[[1]]!="")
  obs = as.numeric(obs[2:length(obs)])

  se = REP[grep("Survey Biomass",REP)[1]:(grep("Survey Biomass",REP)[2]-2)][6]
  se = base::strsplit(se," ")
  se = subset(se[[1]], se[[1]]!="")
  se = as.numeric(se[2:length(se)])

  pred = REP[grep("Survey Biomass",REP)[1]:(grep("Survey Biomass",REP)[2]-2)][3]
  pred = base::strsplit(pred," ")
  pred = subset(pred[[1]], pred[[1]]!="")
  pred = as.numeric(pred[2:length(pred)])


  data.frame(year = syr, biomass = obs, pred = pred, se = se) %>%
    tidytable::mutate.(lci = biomass - 1.96 * se,
                       uci = biomass + 1.96 * se) %>%
    write.csv(here::here(year, folder, "processed", "survey.csv"), row.names = FALSE)


  # recruitment ----

  N = REP[grep("Numbers", REP):(grep("Obs_P_fish_age", REP)-2)]
  t = NA
  for(i in 1:length(yrs)){
    ts = as.numeric(base::strsplit(N[i+1]," ")[[1]][3])
    t = c(t, ts)
  }
  pred_rec = t[!is.na(t)]

  # biomass & F & recruits ----
  data.frame(year = yrs,
             tot_biom = afscassess::rep_item("Tot_biom"),
             sp_biom = afscassess::rep_item("SpBiom"),
             F = afscassess::rep_item("Fully_selected_F"),
             recruits = pred_rec) %>%
    write.csv(here::here(year, folder, "processed", "bio_rec_f.csv"), row.names = FALSE)


  # selectivity ----
  data.frame(age = ages,
             fish = afscassess::rep_item("Fishery_Selectivity"),
             srv1 = afscassess::rep_item("TWL Survey_Selectivity"),
             maturity = afscassess::rep_item("Maturity")) %>%
    write.csv(here::here(year, folder, "processed", "selex.csv"), row.names = FALSE)

  # yield ratio B40 & B35----

  data.frame(B40 = STD$value[which(STD$name=="B40")],
             B35 = as.numeric(REP[(grep("B_35",REP)+1):(grep("F_40",REP)[1]-1)]),
             yld_rat = as.numeric(unlist(base::strsplit(CTL[grep("yieldratio", CTL)], "\t"))[1])) %>%
    write.csv(here::here(year, folder, "processed", "b35_b40_yld.csv"), row.names = FALSE)

  # size comps ----

  #! this will need a switch for multiple surveys

  obs = REP[grep("Obs_P_fish_age",REP):(grep("Pred_P_fish_age",REP)-2)]
  pred = REP[grep("Pred_P_fish_age",REP):(grep("Obs_P_fish_size",REP)-2)]

  obs_l = REP[grep("Obs_P_fish_size",REP):(grep("Pred_P_fish_size",REP)-2)]
  pred_l = REP[grep("Pred_P_fish_size",REP):(grep("Obs_P_srv1_age",REP)-2)]

  s_obs = REP[grep("Obs_P_srv1_age",REP):(grep("Pred_P_srv1_age",REP)-2)]
  s_pred = REP[grep("Pred_P_srv1_age",REP):(grep("Obs_P_srv1_size",REP)-2)]

  s_obs_l = REP[grep("Obs_P_srv1_size",REP):(grep("Pred_P_srv1_size",REP)-2)]

  afscassess::purrit(obs, pred, rec_age, plus_age, comp = "age", lenbins = lenbins) %>%
    write.csv(here::here(year, folder, "processed", "fac.csv"))

  afscassess::purrit(obs_l, pred_l, rec_age, plus_age, comp = "length", lenbins = lenbins) %>%
    write.csv(here::here(year, folder, "processed", "fsc.csv"))

  afscassess::purrit(s_obs, s_pred, rec_age, plus_age, comp = "age", lenbins = lenbins) %>%
    write.csv(here::here(year, folder, "processed", "sac.csv"))

  afscassess::purrit(s_obs_l, pred = NULL, rec_age, plus_age, comp = "length", lenbins = lenbins) %>%
    write.csv(here::here(year, folder, "processed", "ssc.csv"))

}

#' best F
#'
#'  Use this year's model to recalculate the F for the previous year that would have produced the specified OFL for the previous year
#'
#'
#' @param data data.frame in form (age, naa, waa, saa)
#' @param m natural mortality (for both sexes(\code{type=1}) or for sex 1 (\code{type=2} or \code{type=3} estimates))
#' @param last_ofl OFL from the previous assessment
#' @param type 1 = 1 sex - 1 gear; 2 = 1 sex - 2 gear; 3 = 2 sex - 1 gear - 2m; 4 = 2 sex - 2 gear - 2m
#' @param f_ratio  is the final ratio of catch between gears in the previous year, as estimated by this year's model - only needed if using multiple gear types
#' @param m2 is the second estimate of natural mortality
#' @param last_f F value from the previous assessment
#'
#' @return The estimated F_OFL
#' @export best_f
#'
#' @examples
#'\dontrun{
#' oneone <- data.frame(age = 1:9,
#'                      naa = c(15.27, 10.85, 0.24, 8.72, 2.49, 10.66, 3.92, 3.71, 2.6),
#'                      waa = c(0.69, 1.11, 1.54, 1.95, 2.32, 2.64, 2.91, 3.41, 3.33),
#'                      saa = c(0.01,	0.1	, 0.57, 0.94, 0.99, 1, 1,	1, 1))
#'
#'
#' onetwo <- data.frame(age = 1:9,
#'                     naa = c(15.27, 10.85, 0.24, 8.72, 2.49, 10.66, 3.92, 3.71, 2.6),
#'                     waa = c(0.69, 1.11, 1.54, 1.95, 2.32, 2.64, 2.91, 3.41, 3.33),
#'                     saa01 = c(0.01,	0.1	, 0.57, 0.94, 0.99, 1, 1,	1, 1),
#'                     saa02 = c(0.01, 0.06, 0.34, 0.82, 0.97, 1, 1, 1, 1))
#'
#' twotwo = structure(list(age = 2:10,
#'                        naa1 = c(7.64, 5.43, 0.12, 4.36, 1.24,5.33, 1.96, 1.85, 1.3),
#'                        naa2 = c(7.64, 5.42, 0.12, 4.4, 1.27,5.44, 2, 1.89, 1.32),
#'                        waa1 = c(0.97, 1.46, 1.88, 2.22, 2.48, 2.67, 2.81, 2.91, 2.99),
#'                        waa2 = c(0.92, 1.48, 2.05, 2.6, 3.09, 3.52, 3.89, 4.19, 4.44),
#'                        saa11 = c(0.01, 0.1, 0.57, 0.94, 0.99, 1, 1, 1, 1),
#'                        saa21 = c(0.01, 0.06, 0.34, 0.82, 0.97, 1, 1, 1, 1),
#'                        saa12 = c(0.12, 0.37, 0.62, 0.82, 0.94, 1, 1, 0.96, 0.89),
#'                        saa22 = c(0.01, 0.13, 0.38, 0.66, 0.89, 1, 1, 0.92, 0.79)),
#'                   class = "data.frame", row.names = c(NA, -9L))
#'
#' type 1
#' best_f(m = 0.1, last_ofl = 15, data = oneone, type = 1)
#' type 2
#' best_f(m = 0.1, last_ofl = 15, data = onetwo, type = 2, f_ratio = 0.1)
#' type 3
#' best_f(m = 0.12, last_ofl = 15, data = twotwo, type = 3, f_ratio = 0.1, m2 = 0.08)
#' type 4
#' best_f(m = 0.12, last_ofl = 15, data = twotwo, type = 4, f_ratio = 0.1, m2 = 0.08)
#' }

#'

best_f <- function(data, m, last_ofl, type = 1, f_ratio = NULL, m2 = NULL, last_f){

  # helper function
  catch = function(naa, waa, saa, m, f_ofl, f_ratio = NULL, gear = NULL){


    if(!is.null(gear)) f_ratio = 1 - f_ratio
    if(!is.null(f_ratio)){
      naa * waa * saa * f_ofl * f_ratio /
        (f_ofl * f_ratio * saa + m) * (1 - exp(-f_ofl * f_ratio * saa - m))

    } else {

      naa * waa * saa * f_ofl /
        (f_ofl * saa + m) * (1 - exp(-f_ofl * saa - m))
    }
  }


  if(type == 1){
    # 1 sex - 1 gear
    names(data) <- c("age", "naa", "waa", "saa")

    naa = data$naa
    waa = data$waa
    saa = data$saa

    g = function(f_ofl){
      (last_ofl - sum(catch(naa, waa, saa, m, f_ofl)))^2
    }
  }

  if(type == 2){
    # 1 sex - 2 gear

    if(is.null(f_ratio)){
      stop("you must provide the final ratio of catch between gears in the
             previous year, as estimated by this year's model ")
    }

    names(data) <- c("age", "naa", "waa", "saa_01", "saa_02")

    naa = data$naa
    waa = data$waa
    saa_01 = data$saa_01
    saa_02 = data$saa_02

    g = function(f_ofl){
      (last_ofl - sum(catch(naa, waa, saa_01, m, f_ofl, f_ratio, gear = 1),
                      catch(naa, waa, saa_02, m, f_ofl, f_ratio)))^2
    }
  }

  if(type == 3){
    # 2 sex - 1 gear - 2 M

    if(is.null(m2)){
      stop("you must provide two natural mortality values")
    }

    names(data) <- c("age", "naa_1", "naa_2", "waa_1", "waa_2", "saa_1", "saa_2")

    naa_1 = data$naa_1
    naa_2 = data$naa_2
    waa_1 = data$waa_1
    waa_2 = data$waa_2
    saa_1 = data$saa_1
    saa_2 = data$saa_2



    g = function(f_ofl){
      (last_ofl - sum(catch(naa_1, waa_1, saa_1, m, f_ofl),
                      catch(naa_2, waa_2, saa_2, m2, f_ofl)))^2
    }
  }

  if(type == 4){
    # 2 sex - 2 gear - 2 M

    if(is.null(f_ratio)){
      stop("you must provide the final ratio of catch between gears in the
             previous year, as estimated by this year's model ")
    }

    if(is.null(m2)){
      stop("you must provide two natural mortality values")
    }

    names(data) <- c("age", "naa_1", "naa_2", "waa_1", "waa_2", "saa_11", "saa_21", "saa_12", "saa_22")

    naa_1 = data$naa_1
    naa_2 = data$naa_2
    waa_1 = data$waa_1
    waa_2 = data$waa_2
    saa_11 = data$saa_11
    saa_21 = data$saa_21
    saa_12 = data$saa_12
    saa_22 = data$saa_22


    g = function(f_ofl){
      (last_ofl - sum(catch(naa_1, waa_1, saa_11, m, f_ofl, f_ratio, gear = 1),
                      catch(naa_2, waa_2, saa_21, m2, f_ofl, f_ratio, gear = 1),
                      catch(naa_1, waa_1, saa_12, m, f_ofl, f_ratio),
                      catch(naa_2, waa_2, saa_22, m2, f_ofl, f_ratio)))^2
    }
  }

  f_ofl = optimize(g, c(0, 1))$minimum

  if(f_ofl < last_f){
    stop("The F from the previous assessment is greater than F_OFL: something is wrong!!")
  } else {
    paste0("Type ", type, " corresponding F_OFL = ", f_ofl)
  }

}

run_retro <- function(year, model, model_name, dat_name, mcmc = 10000000, mcsave=2000) {

  # setup
  if (!dir.exists(here::here(year, model, "retro", "model"))){
    dir.create(here::here(year, model, "retro", "model"), recursive=TRUE)
  }

  if (!dir.exists(here::here(year, model, "retro", "results"))){
    dir.create(here::here(year, model, "retro", "results"), recursive=TRUE)
  }


  file.copy(here::here(year, model, paste0(model_name, ".exe")),
            here::here(year, model, "retro", "model"))

  if(file.exists(here::here(year, model, "mat.dat"))){
    file.copy(here::here(year, model, "mat.dat"),
              here::here(year, model, "retro", "model"))
  }


  ctl <- read.delim(here::here(year, model, paste0(dat_name, ".ctl")), header=F)
  dat <- readLines(here::here(year, model, paste0(dat_name, ".dat")))


  Sec_st = grep("#-",dat)
  Sec_end = grep("#!",dat)

  st_end<-matrix(NA,nrow=length(Sec_st),ncol=2)
  st_end[,1] = Sec_st
  st_end[,2] = Sec_end

  mcmcon = "YES"
  mcmcruns = mcmc# Could change these, but I like 5000 as a manageable number to deal with
  mcmcsave = mcsave

  styr = as.numeric(dat[Sec_st[2]-3]) # start of model (example 1961 for POP)
  nages = as.numeric(dat[Sec_st[2]+3]) # number of age bins
  nlens = as.numeric(dat[Sec_st[2]+5]) # number of length bins

  # Set up some results files
  RES_SB <-matrix(nrow=length(seq(styr, year)),ncol=n_retro)
  rownames(RES_SB) <-seq(styr, year)
  colnames(RES_SB) <-seq( year-n_retro+1, year)
  RES_Rec<-matrix(nrow=length(seq(styr, year)),ncol=n_retro)
  rownames(RES_Rec)<-seq(styr, year)
  colnames(RES_Rec)<-seq( year-n_retro+1, year)

  T_start<-Sys.time() #Timer start

  for(y in 1:n_retro){

    # Set endyr
    yrs_retro<-seq(year -n_retro+1, year )
    endyr<-yrs_retro[y]
    nyrs<-endyr-styr+1
    DAT_retro<-c(dat[st_end[1,1]:st_end[1,2]],as.character(endyr), dat[st_end[2,1]:st_end[2,2]])

    # Fishery catch
    DAT_retro<-c(DAT_retro,paste(scan(text=dat [Sec_st[3]-1])[1:nyrs],collapse=" "),dat [st_end[3,1]:st_end[3,2]])


    # Trawl survey biomass
    BTSb_yrs<-length(which(scan(text=dat [Sec_st[5]-1])<=endyr))
    DAT_retro<-c(
      DAT_retro,
      as.character(BTSb_yrs),
      dat [st_end[4,1]:st_end[4,2]],
      paste(scan(text=dat [Sec_st[5]-1])[1:BTSb_yrs],collapse=" "),
      dat [st_end[5,1]:st_end[5,2]],
      paste(scan(text=dat [Sec_st[6]-1])[1:BTSb_yrs],collapse=" "),
      dat [st_end[6,1]:st_end[6,2]],
      paste(scan(text=dat [Sec_st[7]-1])[1:BTSb_yrs],collapse=" "),
      dat [st_end[7,1]:st_end[7,2]],
      paste(scan(text=dat [Sec_st[8]-1])[1:BTSb_yrs],collapse=" "),
      dat[st_end[8,1]:st_end[8,2]],
      paste(scan(text=dat[Sec_st[9]-1])[1:BTSb_yrs],collapse=" "),
      dat[st_end[9,1]:st_end[9,2]])

    # Fish age comp
    FAC_yrs<-length(which(scan(text=dat [Sec_st[11]-1])<(endyr-1)))
    DAT_retro<-c(DAT_retro,
                 as.character(FAC_yrs),
                 dat [st_end[10,1]:st_end[10,2]],
                 paste(scan(text=dat [Sec_st[11]-1])[1:FAC_yrs],collapse=" "),
                 dat [st_end[11,1]:st_end[11,2]],
                 paste(scan(text= dat[Sec_st[12]-1])[1:FAC_yrs],collapse=" "),
                 dat [st_end[12,1]:st_end[12,2]],
                 paste(scan(text= dat[Sec_st[13]-1])[1:FAC_yrs],collapse=" "),
                 dat [st_end[13,1]:st_end[13,2]],
                 paste(scan(text= dat[Sec_st[14]-1])[1:FAC_yrs],collapse=" "),
                 dat [st_end[14,1]:st_end[14,2]])
    for(i in 1:FAC_yrs)   DAT_retro<-c(DAT_retro,paste(scan(text=dat[Sec_st[15]-FAC_yrs-1+i]) ,collapse = " "))
    DAT_retro<-c(DAT_retro,dat[st_end[15,1]:st_end[15,2]])

    # Survey age comp
    SAC_yrs<-length(which(scan(text=dat[Sec_st[17]-1])<=(endyr-1)))
    DAT_retro<-c(DAT_retro,
                 as.character(SAC_yrs),
                 dat[st_end[16,1]:st_end[16,2]],
                 paste(scan(text=dat[Sec_st[17]-1])[1:SAC_yrs],collapse=" "),
                 dat[st_end[17,1]:st_end[17,2]],
                 paste(scan(text=dat[Sec_st[18]-1])[1:SAC_yrs],collapse=" "),
                 dat[st_end[18,1]:st_end[18,2]],
                 paste(scan(text=dat[Sec_st[19]-1])[1:SAC_yrs],collapse=" "),
                 dat[st_end[19,1]:st_end[19,2]],
                 paste(scan(text=dat[Sec_st[20]-1])[1:SAC_yrs],collapse=" "),
                 dat[st_end[20,1]:st_end[20,2]])
    for(i in 1:SAC_yrs)   DAT_retro<-c(DAT_retro,paste(scan(text=dat[Sec_st[21]-SAC_yrs-1+i]) ,collapse = " "))
    DAT_retro<-c(DAT_retro,dat[st_end[21,1]:st_end[21,2]])

    # Fish size comp
    FSC_yrs<-length(which(scan(text=dat[Sec_st[23]-1])<=(endyr-1)))
    DAT_retro<-c(DAT_retro,
                 as.character(FSC_yrs),
                 dat[st_end[22,1]:st_end[22,2]],
                 paste(scan(text=dat[Sec_st[23]-1])[1:FSC_yrs],collapse=" "),
                 dat[st_end[23,1]:st_end[23,2]],
                 paste(scan(text=dat[Sec_st[24]-1])[1:FSC_yrs],collapse=" "),
                 dat[st_end[24,1]:st_end[24,2]],
                 paste(scan(text=dat[Sec_st[25]-1])[1:FSC_yrs],collapse=" "),
                 dat[st_end[25,1]:st_end[25,2]],
                 paste(scan(text=dat[Sec_st[26]-1])[1:FSC_yrs],collapse=" "),
                 dat[st_end[26,1]:st_end[26,2]])
    for(i in 1:FSC_yrs)   DAT_retro<-c(DAT_retro,paste(scan(text=dat[Sec_st[27]-FSC_yrs-1+i]) ,collapse = " "))
    DAT_retro<-c(DAT_retro,dat[st_end[27,1]:st_end[27,2]])

    # Survey size comp
    SSC_yrs<-length(which(scan(text=dat[Sec_st[29]-1])<=endyr))
    DAT_retro<-c(DAT_retro,
                 as.character(SSC_yrs),
                 dat[st_end[28,1]:st_end[28,2]],
                 paste(scan(text=dat[Sec_st[29]-1])[1:SSC_yrs],collapse=" "),
                 dat[st_end[29,1]:st_end[29,2]],
                 paste(scan(text=dat[Sec_st[30]-1])[1:SSC_yrs],collapse=" "),
                 dat[st_end[30,1]:st_end[30,2]],
                 paste(scan(text=dat[Sec_st[31]-1])[1:SSC_yrs],collapse=" "),
                 dat[st_end[31,1]:st_end[31,2]],
                 paste(scan(text=dat[Sec_st[32]-1])[1:SSC_yrs],collapse=" "),
                 dat[st_end[32,1]:st_end[32,2]])
    for(i in 1:SSC_yrs)   DAT_retro<-c(DAT_retro,paste(scan(text=dat[Sec_st[33]-SSC_yrs-1+i]) ,collapse = " "))
    DAT_retro<-c(DAT_retro,dat[st_end[33,1]:st_end[33,2]])

    # Write data and control file
    write.table(DAT_retro, file=here::here(year, model, "retro", 'model',paste0("goa_nr_",endyr,".dat")),
                quote=FALSE,row.names=FALSE,col.names=FALSE)

    ctl[2,1] = paste0("goa_nr_", endyr, ".dat")
    ctl[5,1] = as.character(endyr)
    write.table(ctl, file=here::here(year, model, "retro", 'model', paste0(dat_name, ".ctl")),
                quote=FALSE,row.names=FALSE,col.names=FALSE)

    #Updated to account for fact that .tpl is looking for 2018.
    # write.table(CTL_retro,file=paste(pathM,"/goa_",species,"_",modelyear,".ctl",sep=""),quote=FALSE,row.names=FALSE,col.names=FALSE)


    #/\/\/\/\/\/\/\/\ Run model
    #/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\

    ## set your number of MCMC runs at the top of the program...
    setwd(here::here(year, model, "retro", "model"))

    #Compile the Model
    # compile_admb(MDL_name)

    #Determine Operating system

    if(mcmcon=="YES") {
      system(paste0(model_name, '.exe',' -mcmc ',mcmcruns,' -mcsave ',mcmcsave))
    }else {
      # system(paste(MDL_name,'.exe ','-nox'))
      run_admb(model_name, verbose=TRUE)
    }


    #/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\
    #/\/\/\/\/\/\/\/\ Get/write results
    #/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\

    if(mcmcon=="YES") {
      system(paste0(model_name, '.exe',' -mceval'))
      file.copy(from=here::here(year, model, "retro", "model", "evalout.prj"),
                to=here::here(year, model, "retro", "results", paste0("mcmc_", endyr,".std")),
                overwrite=TRUE)
      file.copy(from=here::here(year, model, "retro", "model", paste0(model_name, ".std")),
                to=  here::here(year, model, "retro", "results", paste0("std_", endyr,".std")),
                overwrite=TRUE)
    }

    # Compile SSB/recruitment results

    STD <- read.delim(here::here(year, model, "retro", "model", paste0(model_name, ".std")),
                      header = T, sep = "")
    RES_SB[1:nyrs,y] <- STD$value[which(STD$name=="spawn_biom")]
    RES_Rec[1:nyrs,y] <- STD$value[which(STD$name=="pred_rec")]

    #---------------------------------------------
    # End of retrospective model running loop
    #---------------------------------------------
  }

  write.csv(RES_SB, here::here(year, model, "retro", "results", "RES_SB.csv"))
  write.csv(RES_Rec, here::here(year, model, "retro", "results", "RES_Rec.csv"))

  T_end<-Sys.time()

  T_end-T_start

}

#' @param year  assessment year
#' @param model_dir   directory of model being evaluated (folder name)
#' @export fac_table
fac_table <- function(year, model_dir){

  options(scipen = 999)
  fsc = read.csv(here::here(year, "data", "output", "fsh_age_comp.csv"))

  fsc %>%
    dplyr::select(n_s, n_h) %>%
    t(.) %>%
    as.data.frame() %>%
    tibble::rownames_to_column("name") -> samps

  fsc %>%
    dplyr::select(-n_s, -n_h, -AA_Index) %>%
    tidyr::pivot_longer(-c(year)) %>%
    tidyr::pivot_wider(names_from = year, values_from = value, names_prefix = "y") %>%
    as.data.frame() %>%
    dplyr::mutate_if(is.numeric, round, digits = 4) %>%
    dplyr::mutate(name = gsub("X", "", name),
                  name = ifelse(dplyr::row_number() == dplyr::n(), paste0(name, "+"), name )) %>%
    dplyr::rename_all(~stringr::str_replace(., "y", "")) -> comp

  names(samps) <- names(comp)

  dplyr::bind_rows(comp, samps) %>%
    write.csv(here::here(model_dir, "tables", "fac.csv"), row.names = FALSE)

}

#' @param year  assessment year
#' @param model_dir   directory of model being evaluated (folder name)
#' @export fsc_table
fsc_table <- function(year, model_dir){

  options(scipen = 999)
  fsc = read.csv(here::here(year, "data", "output", "fsh_length_comp.csv"))

  fsc %>%
    dplyr::select(n_s, n_h) %>%
    t(.) %>%
    as.data.frame() %>%
    tibble::rownames_to_column("name") -> samps

  fsc %>%
    dplyr::select(-n_s, -n_h, -SA_Index) %>%
    tidyr::pivot_longer(-c(year)) %>%
    tidyr::pivot_wider(names_from = year, values_from = value, names_prefix = "y") %>%
    as.data.frame() %>%
    dplyr::mutate_if(is.numeric, round, digits = 4) %>%
    dplyr::mutate(name = gsub("X", "", name),
                  name = ifelse(dplyr::row_number() == dplyr::n(), paste0(name, "+"), name )) %>%
    dplyr::rename_all(~stringr::str_replace(., "y", "")) -> comp

  names(samps) <- names(comp)

  dplyr::bind_rows(comp, samps) %>%
    write.csv(here::here(model_dir, "tables", "fsc.csv"), row.names = FALSE)

}

#' @param year  assessment year
#' @param model_dir   directory of model being evaluated (folder name)
#' @export sac_table
sac_table <- function(year, model_dir){

  options(scipen = 999)
  fsc = read.csv(here::here(year, "data", "output", "goa_bts_age_comp.csv"))

  fsc %>%
    dplyr::select(n_s, n_h) %>%
    t(.) %>%
    as.data.frame() %>%
    tibble::rownames_to_column("name") -> samps

  fsc %>%
    dplyr::select(-n_s, -n_h, -AA_Index) %>%
    tidyr::pivot_longer(-c(year)) %>%
    tidyr::pivot_wider(names_from = year, values_from = value, names_prefix = "y") %>%
    as.data.frame() %>%
    dplyr::mutate_if(is.numeric, round, digits = 4) %>%
    dplyr::mutate(name = gsub("X", "", name),
                  name = ifelse(dplyr::row_number() == dplyr::n(), paste0(name, "+"), name )) %>%
    dplyr::rename_all(~stringr::str_replace(., "y", "")) -> comp

  names(samps) <- names(comp)

  dplyr::bind_rows(comp, samps) %>%
    write.csv(here::here(model_dir, "tables", "sac.csv"), row.names = FALSE)

}

#' @param year  assessment year
#' @param model_dir   directory of model being evaluated (folder name)
#' @export ssc_table
ssc_table <- function(year, model_dir){

  options(scipen = 999)
  fsc = read.csv(here::here(year, "data", "output", "goa_ts_length_comp.csv"))

  fsc %>%
    dplyr::select(n_s, n_h) %>%
    t(.) %>%
    as.data.frame() %>%
    tibble::rownames_to_column("name") -> samps

  fsc %>%
    dplyr::select(-n_s, -n_h, -SA_Index) %>%
    tidyr::pivot_longer(-c(year)) %>%
    tidyr::pivot_wider(names_from = year, values_from = value, names_prefix = "y") %>%
    as.data.frame() %>%
    dplyr::mutate_if(is.numeric, round, digits = 4) %>%
    dplyr::mutate(name = gsub("X", "", name),
                  name = ifelse(dplyr::row_number() == dplyr::n(), paste0(name, "+"), name )) %>%
    dplyr::rename_all(~stringr::str_replace(., "y", "")) -> comp

  names(samps) <- names(comp)

  dplyr::bind_rows(comp, samps) %>%
    write.csv(here::here(model_dir, "tables", "ssc.csv"), row.names = FALSE)

}

recruit_tbl <- function(year, model, model_name, rec_age){

  # hard to filter year with year so change the name
  mod_year = year

  if (!dir.exists(here::here(year, model, "processed"))){
    stop("must run 'process_results' before creating tables")
  }

  # read in data
  REP <- readLines(here::here(year, model, paste0(model_name, ".rep")))
  STD <- read.delim(here::here(year, model, paste0(model_name, ".std")), sep="", header = TRUE)

  suppressWarnings(data.frame(year = unlist(strsplit(REP[grep("Year", REP)[1]]," "))) %>%
                     dplyr::mutate(year = as.numeric(year)) %>%
                     tidyr::drop_na() %>%
                     dplyr::pull(year)) -> yrs

  bio_rec = read.csv(here::here(year, model, "processed", "bio_rec_f.csv")) %>%
    dplyr::select(-F) %>%
    dplyr::mutate(recruits = round(recruits * 1000))

  bio_rec %>%
    dplyr::filter(year %in% (1977 + rec_age):(mod_year - rec_age)) %>%
    dplyr::summarise(recruits = round(mean(recruits))) -> pred_rec

  data.frame(year = mod_year + 1:2,
             tot_biom = STD$value[which(STD$name=="tot_biom_proj")][1:2],
             sp_biom = STD$value[which(STD$name=="spawn_biom_proj")][1:2],
             recruits = pred_rec$recruits) -> std_data

  values = dplyr::bind_rows(bio_rec, std_data) %>%
    dplyr::mutate_all(dplyr::funs(round(.)))



  # get mcmc data - clean it, calculate annual uci and lci
  read.csv(here::here(year, model, "processed", "mceval.csv")) %>%
    dplyr::select(dplyr::starts_with(c( "tot_biom", "spawn_biom", "log_rec_dev", "rec_proj"))) %>%
    dplyr::mutate(id = 1:dplyr::n()) %>%
    dplyr::mutate_if(is.character, dplyr::funs(as.numeric(gsub(",", "", .)))) %>%
    tidyr::pivot_longer(-id) %>%
    dplyr::mutate(value = ifelse(grepl("log", name), exp(value), value),
                  year =  as.numeric(stringr::str_extract(name, "[[:digit:]]+")),
                  name = dplyr::case_when(stringr::str_detect(name, "tot_biom") ~ "tot_biom",
                                          stringr::str_detect(name, "spawn_biom") ~ "sp_biom",
                                          stringr::str_detect(name, "log_rec") ~ "recruits",
                                          stringr::str_detect(name, "rec_proj") ~ "recruits")) %>%
    group_by(year, name) %>%
    dplyr::summarise(lci = quantile(value, 0.025),
                     uci = quantile(value, 0.975)) %>%
    dplyr::left_join(values, .) %>%
    dplyr::filter(year >= 1977 & year <= mod_year + rec_age) %>%
    dplyr::mutate(tot_lci = ifelse(name == 'tot_biom', lci, NA),
                  tot_uci = ifelse(name == 'tot_biom', uci, NA),
                  sp_lci = ifelse(name == 'sp_biom', lci, NA),
                  sp_uci = ifelse(name == 'sp_biom', uci, NA),
                  rec_lci = ifelse(name == 'recruits', lci * 1000, NA),
                  rec_uci = ifelse(name == 'recruits', uci * 1000, NA)) %>%
    group_by(year) %>%
    dplyr::summarise(recruits = mean(recruits, na.rm = T),
                     rec_lci = mean(rec_lci, na.rm = T),
                     rec_uci = mean(rec_uci, na.rm = T),
                     tot_biom = mean(tot_biom, na.rm = T),
                     tot_lci = mean(tot_lci, na.rm = T),
                     tot_uci = mean(tot_uci, na.rm = T),
                     sp_biom = mean(sp_biom, na.rm = T),
                     sp_lci = mean(sp_lci, na.rm = T),
                     sp_uci = mean(sp_uci, na.rm = T)) %>%
    write.csv(here::here(year, model, "tables", "ts.csv"), row.names = FALSE)

}

#' @param year  assessment year
#' @param model_dir  full path of model being evaluated
#' @param MCMC = logical, does this run include MCMC evaluations to be processed?
#' @param no_mcmc = number of mcmc runs
#' @param rec_age recruitment age
#' @param plus_age plus age group
#' @param mcsave the number of mcmcs saved
#' @param ... future functions

process_results_pop <- function(year = 2023,
                                model_dir = NULL,
                                rec_age = 2,
                                plus_age = 25,
                                size_bins = NULL,
                                MCMC=FALSE,
                                no_mcmc = 100000,
                                mcsave = 100, ...){

  # setup

  if (!dir.exists(paste0(model_dir, "/processed"))){
    dir.create(paste0(model_dir, "/processed"), recursive=TRUE)
  }

  if (!dir.exists(paste0(model_dir, "/figs"))){
    dir.create(paste0(model_dir, "/figs"), recursive=TRUE)
  }
  if (!dir.exists(paste0(model_dir, "/tables"))){
    dir.create(paste0(model_dir, "/tables"), recursive=TRUE)
  }


  # helper functions
  rep_item <- function(name){
    t <- strsplit(REP[grep(name, REP)]," ")
    t <- subset(t[[1]], t[[1]]!="")
    if(t[[1]][1] == "TWL"){
      as.numeric(t[3:length(t)])
    } else {
      as.numeric(t[2:length(t)])
    }
  }


  # read in report and ctl files
  dats <- list.files(model_dir, full.names = TRUE)
  DAT <- readLines(dats[grepl("*.dat",dats) & !grepl('proj',dats)])
  REP <- readLines(list.files(model_dir, pattern="*.rep", full.names = TRUE))
  modname <- gsub(".rep","",basename(list.files(model_dir, pattern="*.rep", full.names = TRUE))) ## strip model name from rep file
  CTL <- readLines(list.files(model_dir, pattern="*.ctl", full.names = TRUE))
  STD <- read.delim(list.files(model_dir, pattern="*.std", full.names = TRUE), sep="", header = TRUE)


  # clean rep file
  suppressWarnings(data.frame(year = unlist(base::strsplit(REP[grep("Year", REP)[1]]," "))) %>%
                     tidytable::mutate.(year = as.numeric(year)) %>%
                     tidytable::drop_na.() %>%
                     tidytable::pull.(year)) -> yrs

  suppressWarnings(data.frame(age = unlist(base::strsplit(REP[grep("Age", REP)[1]]," "))) %>%
                     tidytable::mutate.(age = as.numeric(age)) %>%
                     tidytable::drop_na.() %>%
                     tidytable::pull.(age)) -> ages

  styr_rec <- yrs[1] - length(ages) + rec_age

  suppressWarnings(as.data.frame(cbind(yrs = yrs, ages = ages, styr_rec = styr_rec)) %>%
                     tidytable::mutate.(ages = replace(ages, duplicated(ages), NA),
                                        styr_rec = replace(styr_rec, duplicated(styr_rec), NA))) %>%
    write.csv(paste0(model_dir, "/processed/ages_yrs.csv"), row.names = FALSE)

  ## pull out likelihoods ----
  ## this function extracts & binds them rowwise
  do.call(rbind,
          lapply(unlist(base::strsplit(REP[grep('Likelihood|Priors|Penalty|Objective', REP)][2:19],"\n")),
                 FUN = function(x){
                   tmpl <- unlist(strsplit(x," "))
                   tempr <- matrix(c(tmpl[1], tmpl[2], paste0(tmpl[3:length(tmpl)], collapse = ' ')), ncol = 3)
                   return(tempr)
                 })) %>%
    data.frame() %>%
    tidytable::mutate(value = as.numeric(X2),
                      weight = as.numeric(X1)) %>%
    tidytable::select(weight, value , variable = X3) %>%
    tidytable::mutate(model = basename(model_dir),
                      weight = ifelse(is.na(weight) == TRUE, 1, weight)) -> likes
  write.csv(likes, paste0(model_dir, "/processed/likelihoods.csv"), row.names = FALSE)


  # # MCMC parameters ----
  # if(MCMC){
  #   mceval <- read.delim(list.files(model_dir, pattern="*evalout.prj", full.names = TRUE), sep="", header = FALSE)
  #   PSV <- file(paste0(model_dir,"/",modname,".psv"), "rb")
  #
  #   npar = readBin(PSV, what = integer(), n=1)
  #   mcmcs = readBin(PSV, what = numeric(), n = (npar * no_mcmc / mcsave))
  #   close(PSV)
  #   mcmc_params = matrix(mcmcs, byrow=TRUE, ncol=npar)
  #   # thin the string
  #   mcmc_params = mcmc_params[501:nrow(mcmc_params),]
  #   colnames(mcmc_params) = STD$name[1:ncol(mcmc_params)]
  #   write.csv(mcmc_params,paste0(model_dir, "/processed/mcmc.csv"), row.names = FALSE)
  #
  #   # mceval phase output ----
  #
  #   #Curry's Change
  #   mceval = mceval[501:nrow(mceval),]
  #
  #   #Length colnames = 286
  #   # columns mcmc_other = 271
  #
  #   #1-8: Through objective function value
  #
  #   colnames(mceval) = c("sigr", "q_srv1", "q_srv2", "F40", "natmort", "spawn_biom_proj",
  #                        "ABC", "obj_fun",
  #                        paste0("tot_biom_", yrs),
  #                        paste0("log_rec_dev_", seq(styr_rec, yrs[length(yrs)])),
  #                        paste0("spawn_biom_", yrs),
  #                        "log_mean_rec",
  #                        paste0("spawn_biom_proj_", max(yrs) + 1:15),
  #                        paste0("pred_catch_proj_", max(yrs) + 1:15),
  #                        paste0("rec_proj_", max(yrs) + 1:10),
  #                        paste0("tot_biom_proj_", max(yrs)))
  #   write.csv(mceval, paste0(model_dir, "/processed/mceval.csv"), row.names = FALSE)
  #
  # } ## end if MCMC == T

  # catch data ----
  pred = base::strsplit(REP[grep("Pred_Catch", REP)], " ")
  r1 = which(pred[[1]] == "Pred_Catch")
  r2 = which(pred[[1]] == "Pred_catch_later")
  r3 = which(pred[[1]] == "")
  pred = as.numeric(pred[[1]][-c(r1, r2, r3)])

  obs = base::strsplit(REP[grep("Obs_Catch", REP)], " ")
  r1 = which(obs[[1]] == "Obs_Catch")
  r2 = which(obs[[1]] == "Obs_Catch_Later")
  r3 = which(obs[[1]] == "")
  obs = as.numeric(obs[[1]][-c(r1, r2, r3)])

  data.frame(year = yrs, obs = obs, pred = pred) -> catch
  write.csv(catch, paste0(model_dir, "/processed/catch.csv"), row.names = FALSE)

  # survey data ----
  syr = REP[grep("Survey Biomass",REP)[1]:(grep("Survey Biomass",REP)[2]-2)][2]
  syr = base::strsplit(syr," ")
  syr = subset(syr[[1]], syr[[1]]!="")
  syr = as.numeric(syr[2:length(syr)])

  obs = REP[grep("Survey Biomass",REP)[1]:(grep("Survey Biomass",REP)[2]-2)][4]
  obs = base::strsplit(obs," ")
  obs = subset(obs[[1]], obs[[1]]!="")
  obs = as.numeric(obs[2:length(obs)])

  se = REP[grep("Survey Biomass",REP)[1]:(grep("Survey Biomass",REP)[2]-2)][5]
  se = base::strsplit(se," ")
  se = subset(se[[1]], se[[1]]!="")
  se = as.numeric(se[2:length(se)])

  pred = REP[grep("Survey Biomass",REP)[1]:(grep("Survey Biomass",REP)[2]-2)][3]
  pred = base::strsplit(pred," ")
  pred = subset(pred[[1]], pred[[1]]!="")
  pred = as.numeric(pred[2:length(pred)])


  data.frame(year = syr, biomass = obs, pred = pred, se = se) %>%
    tidytable::mutate.(lci = biomass - 1.96 * se,
                       uci = biomass + 1.96 * se,
                       lci = ifelse(lci < 0, 0 ,lci)) -> srv
  write.csv(srv, paste0(model_dir, "/processed/survey.csv"), row.names = FALSE)

  # recruitment ----
  N = REP[grep("Numbers", REP):(grep("Obs_P_fish_age", REP)-2)]
  t = NA
  for(i in 1:length(yrs)){
    ts = as.numeric(base::strsplit(N[i+1]," ")[[1]][3])
    t = c(t, ts)
  }
  pred_rec = t[!is.na(t)]

  # biomass & F & recruits ----
  data.frame(year = yrs,
             tot_biom = rep_item("Tot_biom"),
             sp_biom = rep_item("SpBiom"),
             F = rep_item("Fully_selected_F"),
             recruits = pred_rec) -> bio_rec_f
  write.csv(bio_rec_f, paste0(model_dir, "/processed/bio_rec_f.csv"), row.names = FALSE)

  # selectivity ----
  data.frame(age = ages,
             fish1 = rep_item("Fishery_Selectivity_1967-1976"),
             fish2 = rep_item("Fishery_Selectivity_1977-1995"),
             fish3 = rep_item("Fishery_Selectivity_1996-2006"),
             fish4 = rep_item(paste0("Fishery_Selectivity_2007-", year)),
             srv1 = rep_item("Trawl_Survey_Selectivity")) -> selex
  write.csv(selex, paste0(model_dir, "/processed/selex.csv"), row.names = FALSE)

  # weight-at-age, maturity ----
  data.frame(age = ages,
             srv1 = rep_item("Weight"),
             maturity = rep_item("Maturity")) -> waa_mat
  write.csv(waa_mat, paste0(model_dir, "/processed/selex.csv"), row.names = FALSE)

  # number of parameters ----

  data.frame(num_param = as.numeric(REP[(grep("Number parameters estimated",REP)+1)])) -> num_param

  # key parameters ----

  data.frame(q_trawl = as.numeric(REP[(grep("q_trawl",REP)+1):(grep("nat_mort",REP)[1]-1)]),
             nat_mort = as.numeric(REP[(grep("nat_mort",REP)+1):(grep("sigr",REP)[1]-1)]),
             sigr = as.numeric(REP[(grep("sigr",REP)+1):(grep("log_mean_rec",REP)[1]-1)]),
             log_mean_rec = as.numeric(REP[(grep("log_mean_rec",REP)+1):(grep("q_alt",REP)[1]-1)])) -> key_param
  write.csv(key_param, paste0(model_dir, "/processed/key_param.csv"), row.names = FALSE)

  # fishery age comp ----

  obs = REP[grep("Obs_P_fish_age",REP):(grep("Pred_P_fish_age",REP)-2)]
  pred = REP[grep("Pred_P_fish_age",REP):(grep("yrs_fish_age",REP)-2)]

  fac <- afscassess::purrit(obs, pred, rec_age, plus_age, comp = "age", lenbins = size_bins)

  afscassess::fac_table(year, model_dir)

  # fishery size comp ----

  obs = REP[grep("Obs_P_fish_size",REP):(grep("Pred_P_fish_size",REP)-2)]
  pred = REP[grep("Pred_P_fish_size",REP):(grep("yrs_fish_size",REP)-2)]

  fsc <- afscassess::purrit(obs, pred, rec_age, plus_age, comp = "length", lenbins = size_bins)

  afscassess::fsc_table(year, model_dir)

  # survey age comp ----

  obs = REP[grep("Obs_P_srv1_age",REP):(grep("Pred_P_srv1_age",REP)-2)]
  pred = REP[grep("Pred_P_srv1_age",REP):(grep("yrs_srv1_age",REP)-2)]

  sac <- afscassess::purrit(obs, pred, rec_age, plus_age, comp = "age", lenbins = size_bins)

  afscassess::sac_table(year, model_dir)
  afscassess::ssc_table(year, model_dir)

  # put all the results into a list ----
  # yrs = years of model
  # ages is ages
  # styr_rec is the start year for recruitment (further back to populate initial abundance)
  # likes is likehood values
  # catch is obs and pred catch
  # srv is obs, pred, and se, lci, uci of obs
  # bio_rec_f is time series of biomass (sp and tot), F, and recruitment
  # selex is selectivity from survey and fishery
  # waa_mat is weight-at-age and maturity
  # num_param is number of parameters in model
  # key_param are parameter estimates for survey q, M, sigr, and mean recruitment (log-scale)
  # fac is obs and pred fishery age comp
  # fsc is obs and pred fishery size comp
  # sac is obs and pred survey age comp



  proc_res <- list(yrs = yrs,
                   ages = ages,
                   styr_rec = styr_rec,
                   likes = likes,
                   catch = catch,
                   srv = srv,
                   bio_rec_f = bio_rec_f,
                   selex = selex,
                   waa_mat = waa_mat,
                   num_param = num_param,
                   key_param = key_param,
                   fac = fac,
                   fsc = fsc,
                   sac = sac)

  proc_res
}
