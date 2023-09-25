#' Plot catch
#'
#' @param year model year
#' @param folder folder name model is in
#' @param save default is TRUE, saves fig to the folder the model is in
#'
#' @return
#' @export plot_catch
#'
#' @examples
plot_catch <- function(year, folder, save=TRUE){

  if (!dir.exists(here::here(year, folder, "processed"))){
    stop("must run 'process_results' before creating figures")
  }
  # set view
  ggplot2::theme_set(afscassess::theme_report())

  vroom::vroom(here::here(year, folder, "processed", "ages_yrs.csv"))$yrs -> yrs
  vroom::vroom(here::here(year, folder, "processed", "catch.csv")) %>%
    tidytable::bind_cols.(year = yrs) %>%
    tidytable::summarise.(Observed = obs / 1000,
                          Estimated = pred / 1000,
                          years = "All years") -> dat

  tidytable::filter.(dat, year %in% (max(yrs) - 20):max(yrs)) %>%
    tidytable::mutate.(years = "Recent years") %>%
    tidytable::bind_rows.(dat) %>%
    tidytable::pivot_longer.(c(-year, -years)) %>%
    ggplot2::ggplot(ggplot2::aes(year, value, color = name, lty = name)) +
    ggplot2::geom_line() +
    scico::scale_color_scico_d(name = "", palette = "roma") +
    ggplot2::scale_linetype_manual(name = "",
                                   values = c(1,1)) +
    ggplot2::facet_wrap(~years, scales = "free",
                        dir = "v") +
    ggplot2::ylab("Catch (kt)") +
    ggplot2::xlab("Year") +
    ggplot2::expand_limits(y = 0) +
    afscassess::scale_x_tickr(data=dat, var=year, start=1960) +
    ggplot2::theme(legend.justification=c(1,0),
                   legend.position=c(0.98,0.8))

  if(isTRUE(save)) {
    ggplot2::ggsave(here::here(year, folder, "figs", "catch.png"),
                    width = 6.5, height = 6.5, units = "in", dpi = 200)
  }

}

#' Plot catch
#'
#' @param year model year
#' @param folder folder name model is in
#' @param save default is TRUE, saves fig to the folder the model is in
#'
#' @return
#' @export plot_biomass
#'
plot_biomass <- function(year, folder, save=TRUE) {

  if (!dir.exists(here::here(year, folder, "processed"))){
    stop("must run 'process_results' before creating figures")
  }
  # set view
  ggplot2::theme_set(afscassess::theme_report())
    bio = read.csv(here::here(year, folder, "processed", "bio_rec_f.csv"))
    yrs = bio$year

    vroom::vroom(here::here(year, folder, "processed", "mceval.csv"))  %>%
      tidytable::select.(paste0("tot_biom_", yrs)) %>%
      tidytable::mutate.(group = 1:tidytable::n.()) %>%
      tidytable::pivot_longer.(-group) %>%
      tidytable::mutate.(year = as.numeric(gsub("tot_biom_", "", name)),
                         name = "Total biomass") %>%
      tidytable::bind_rows.( vroom(here::here(year, folder, "processed", "mceval.csv")) %>%
                               tidytable::select(paste0("spawn_biom_", yrs)) %>%
                               tidytable::mutate.(group = 1) %>%
                               tidytable::pivot_longer.(-group) %>%
                               tidytable::mutate.(year = as.numeric(gsub("spawn_biom_", "", name)),
                                                  name = "Spawning biomass")) %>%
      tidytable::mutate.(name = factor(name, levels = c("Total biomass", "Spawning biomass"))) %>%
      tidytable::summarise.(median = median(value) / 1000,
                            lci = quantile(value, 0.025) / 1000,
                            uci = quantile(value, 0.975) / 1000,
                            .by = c(year, name)) %>%
      tidytable::left_join.(data.frame(year = yrs,
                                       tot = bio$tot_biom / 1000,
                                       bio = bio$sp_biom / 1000)) %>%
      tidytable::mutate.(biomass = ifelse(name == "Total biomass", tot, bio)) %>%
      tidytable::select.(-tot, -bio) -> dat

    dat %>%
      ggplot2::ggplot(ggplot2::aes(year, biomass)) +
      ggplot2::geom_line() +
      ggplot2::geom_ribbon(ggplot2::aes(ymin = lci, ymax = uci), alpha = 0.1) +
      ggplot2::facet_wrap(~name, dir = "v", scales = "free_y") +
      ggplot2::scale_y_continuous(name = "Biomass (kt)", labels = scales::comma) +
      ggplot2::expand_limits(y = 0) +
      afscassess::scale_x_tickr(name = "Year", data=dat, var=year, to=10, start = 1960)

    if(isTRUE(save)) {
        ggplot2::ggsave(here::here(year, folder, "figs", "est_biomass.png"),
                        width = 6.5, height = 6.5, units = "in", dpi = 200)
    }
}

#' Plot age and length compositions
#'
#' @param year = assessment year
#' @param folder = folder the model lives in
#' @return
#' @export plot_comps
#' @param save default is TRUE, saves fig to the folder the model is in
#'
#' @examples
#' plot_comps(year, folder)
plot_comps <- function(year, folder, save = TRUE){


  # SAVE NOT IMPLEMETED


  # set view
  ggplot2::theme_set(afscassess::theme_report())

  fac = read.csv(here::here(year, folder, "processed", "fac.csv"))
  fsc = read.csv(here::here(year, folder, "processed", "fsc.csv"))
  sac = read.csv(here::here(year, folder, "processed", "sac.csv"))
  ssc = read.csv(here::here(year, folder, "processed", "ssc.csv"))

  # function to create plots
  # have to use dplyr, tidytable doesn't play well with enquo
  comp <- function(data, variable, fleet) {
    var <- rlang::enquo(variable)

    data %>%
      dplyr::mutate(level = factor(.data[[var]]),
                    Year = factor(Year)) %>%
      dplyr::filter(groups == "obs") %>%
      ggplot2::ggplot(ggplot2::aes(.data[[var]], value)) +
      ggplot2::geom_col(ggplot2::aes(fill = level), width = 1, color = 1) +
      ggplot2::facet_wrap(~Year, strip.position="right",
                          dir = "v",
                          ncol = 1) +
      ggplot2::geom_point(data = dplyr::filter(data, groups == "pred"),
                          ggplot2::aes(.data[[var]], value, group = 1)) +
      ggplot2::theme(panel.spacing.y = grid::unit(0, "mm")) +
      ggplot2::theme(axis.text.y = ggplot2::element_blank(),
                     axis.ticks.y = ggplot2::element_blank()) +
      ggplot2::xlab(Hmisc::capitalize(variable)) +
      ggplot2::ylab(paste0(Hmisc::capitalize(fleet)," ", variable,  " composition")) +
      afscassess::scale_x_tickr(data=data, var=.data[[var]], start = 0, min = 5) +
      ggplot2::theme(legend.position = "none")
  }

  comp(data=fac, var="age", fleet="fishery")

  ggplot2::ggsave(here::here(year, folder, "figs", "fsh_age_comp.png"),
                  width = 6.5, height = 6.5, units = "in", dpi = 200)

  comp(data=sac, var="age", fleet="survey")

  ggplot2::ggsave(here::here(year, folder, "figs", "srv_age_comp.png"),
                  width = 6.5, height = 6.5, units = "in", dpi = 200)

  comp(data=fsc, var="length", fleet="fishery")

  ggplot2::ggsave(here::here(year, folder, "figs", "fsh_length_comp.png"),
                  width = 6.5, height = 6.5, units = "in", dpi = 200)

  comp(data=ssc, var="length", fleet="survey")

  ggplot2::ggsave(here::here(year, folder, "figs", "srv_length_comp.png"),
                  width = 6.5, height = 6.5, units = "in", dpi = 200)

}


#' Phase plot
#'
#' @param year model year
#' @param folder  folder model is in
#' @param model_name e.g. goa_nr_2020
#' @param save default is TRUE, saves fig to the folder the model is in
#'
#' @return
#' @export plot_phase
#'
#' @examples
plot_phase <- function(year, folder, model_name, save = TRUE){

  if (!dir.exists(here::here(year, folder, "processed"))){
    stop("must run 'process_results' before creating figures")
  }

  STD <- read.delim(here::here(year, folder, paste0(model_name, ".std")), sep="", header = TRUE)
  bio = read.csv(here::here(year, folder, "processed", "bio_rec_f.csv"))
  bby = read.csv(here::here(year, folder, "processed", "b35_b40_yld.csv"))
  yrs = read.csv(here::here(year, folder, "processed", "ages_yrs.csv"))$yrs
  B35 = bby$B35
  yield = bby$yld_rat

  Fabc <- 0.35/0.4

  segs = data.frame(x1 = rep(c(0.05, 0.4/0.35), 2),
                    x2 = rep(c(0.4/0.35, 2.8), 2),
                    y1 = c(0, 1, 0, Fabc),
                    y2 = c(1, 1, Fabc, Fabc),
                    group = factor(c("ofl", "ofl", "abc", "abc"),
                                   levels = c("ofl", "abc")))

  # set view
  ggplot2::theme_set(afscassess::theme_report())

  data.frame(year = min(yrs):(max(yrs) + 2),
             x = c(bio$sp_biom, STD$value[which(STD$name=="spawn_biom_proj")][1],
                   STD$value[which(STD$name=="spawn_biom_proj")][2]) / B35,
             y = c(bio$F, STD$value[which(STD$name=="F40")] * yield,
                   STD$value[which(STD$name=="F40")] * yield) /
               STD$value[which(STD$name=="F35")]) %>%
    tidytable::mutate.(label = stringr::str_sub(year, 3),
                       decade = (floor(year / 10) * 10))  %>%
    ggplot2::ggplot(ggplot2::aes(x, y)) +
    ggplot2::geom_path(color = "darkgray", show.legend = FALSE) +
    ggplot2::geom_label(ggplot2::aes(label=label, color = decade), label.size = 0,
                        show.legend = FALSE, size = 3, family="Times", alpha = 0.5) +
    ggplot2::geom_segment(data = segs, ggplot2::aes(x1, y1, xend = x2, yend = y2, linetype = group)) +
    ggplot2::scale_linetype_manual(values = c(1, 3),
                                   labels = c(expression(italic(F[OFL])),
                                              expression(italic(F[ABC]))),
                                   name = "") +
    scico::scale_color_scico(palette = "roma") +
    ggplot2::ylab(expression(italic(F/F["35%"]))) +
    ggplot2::xlab(expression(italic(SSB/B["35%"]))) +
    ggplot2::theme(legend.justification=c(1,0),
                   legend.position=c(0.9,0.85))

  # Can create a zoomed in figure with this code
  # library(ggforce)
  # facet_zoom(xlim = c(0, 2.8),
  #            ylim = c(0, 1),
  #            zoom.data = z,
  #            horizontal = FALSE) +
  # ggplot2::theme(zoom.y = element_blank(), validate = FALSE)

if(isTRUE(save)) {
  ggplot2::ggsave(here::here(year, folder, "figs", "phase_plane.png"),
         width = 6.5, height = 6.5, units = "in", dpi = 200)
}

}

#' plot random effects model by GOA area
#'
#' @param file  where the files live e.g., "Y:/ABL_MESA/SAFES2020/Apportionment/Plot data - NR"
#' @param save default is TRUE, saves fig to the folder the model is in
#'
#'
#' @return
#' @export plot_re
#'
#' @examples
#' \dontrun{
#' plot_re(location)
#' }
#'
plot_re <- function(file, save=TRUE){

  # IGNORE THIS _ DOESN"T WORK WITH REMA




  list.files(file,
             pattern="_re", full.names = TRUE) %>%
    purrr::map_df( ~ read.csv(.), .id = "region") %>%
    dplyr::mutate(region = dplyr::case_when(region==1 ~ "CGOA",
                                            region==2 ~ "EGOA",
                                            region==3 ~ "WGOA")) %>%
    dplyr::left_join(list.files(file,
                                pattern="_srv", full.names = TRUE) %>%
                       purrr::map_df( ~ read.csv(.), .id = "region") %>%
                       dplyr::mutate(region = dplyr::case_when(region==1 ~ "CGOA",
                                                               region==2 ~ "EGOA",
                                                               region==3 ~ "WGOA")) %>%
                       dplyr::mutate(lci = srv_est - srv_est * srv_sd * 1.96,
                                     uci = srv_est + srv_est * srv_sd * 1.96,
                                     yrs = yrs_srv,
                                     lci = ifelse(lci <0, 0, lci)),
                     by = c("region", "yrs")) %>%
    dplyr::mutate(region = factor(region,
                                  levels = c("WGOA", "CGOA", "EGOA"))) -> dat

  # set view
  ggplot2::theme_set(afscassess::theme_report())
  dat %>%
    ggplot2::ggplot(ggplot2::aes(yrs, biomA / 1000)) +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = LCI / 1000, ymax = UCI / 1000), alpha = 0.15) +
    ggplot2::geom_errorbar(ggplot2::aes(ymin = lci / 1000, ymax = uci / 1000), color = "darkgray") +
    ggplot2::geom_point(ggplot2::aes(y = srv_est / 1000), color = "darkgray") +
    ggplot2::geom_line() +
    ggplot2::facet_wrap(~region, scales = "free_y", nrow = 3, drop=TRUE) +
    ggplot2::ylab("Survey biomass (kt)") +
    ggplot2::xlab("Year") +
    ggplot2::scale_x_continuous(breaks = funcr::tickr(dat, year, 5, start = 1980)$breaks,
                                labels = funcr::tickr(dat, year, 5, start = 1980)$labels)

  if(isTRUE(save)) {
    ggplot2::ggsave(here::here(year, model, "figs", "random_effect.png"),
                  width = 6.5, height = 6.5, units = "in", dpi = 200)
  }
}



#' Plot recruitment estimates
#'
#' @param year assessment year
#' @param folder folder the model is in
#' @param rec_age age at first recruitment
#' @param save default is TRUE, saves fig to the folder the model is in
#'
#' @return
#' @export plot_rec
#'
#' @examples
plot_rec <- function(year, folder, rec_age, save=TRUE){

  recs = paste0(rec_age, "+")
  if (!dir.exists(here::here(year, folder, "processed"))){
    stop("must run 'process_results' before creating figures")
  }

  vroom::vroom(here::here(year, folder, "processed", "ages_yrs.csv"))$yrs -> yrs
  vroom::vroom(here::here(year, folder, "processed", "bio_rec_f.csv")) %>%
    tidytable::select(recruits) %>%
    tidytable::bind_cols(vroom::vroom(here::here(year, folder, "processed", "mceval.csv")) %>%
                           tidytable::select.(log_mean_rec, paste0("log_rec_dev_", yrs)) %>%
                           tidytable::mutate.(group = 1:tidytable::n()) %>%
                           tidytable::pivot_longer.(c(-group, -log_mean_rec), values_to = "rec_dev") %>%
                           tidytable::mutate.(year = as.numeric(gsub("log_rec_dev_", "", name)),
                                              rec_dev = exp(log_mean_rec + rec_dev)) %>%
                           tidytable::summarise.(lci = quantile(rec_dev, 0.025),
                                                 uci = quantile(rec_dev, 0.975),
                                                 .by = year)) -> dat
  # set view
  ggplot2::theme_set(afscassess::theme_report())

  dat %>%
    ggplot2::ggplot(ggplot2::aes(year, recruits)) +
    ggplot2::geom_col(width = 0.8, alpha = 0.6) +
    ggplot2::geom_errorbar(ggplot2::aes(ymin = lci, ymax = uci), width = 0.5) +
    afscassess::scale_x_tickr(name = "Year", data= dat, var = year, to = 10, start = 1960) +
    ggplot2::ylab(paste0("Age-", recs, " Recruitment (millions)"))

  if(isTRUE(save)) {
    ggplot2::ggsave(here::here(year, folder, "figs", "recruits.png"),
                    width = 6.5, height = 5.5, units = "in", dpi = 200)
  }
}

#' Recruitment/SSB plot
#'
#' @param year model year
#' @param folder folder model is in
#' @param rec_age age at first recruitment
#' @param save default is TRUE, saves fig to the folder the model is in
#'
#' @return
#' @export plot_rec_ssb
#'
#' @examples plot_rec_ssb(year, folder, rec_age)
plot_rec_ssb <- function(year, folder, rec_age, save=TRUE){

  if (!dir.exists(here::here(year, folder, "processed"))){
    stop("must run 'process_results' before creating figures")
  }

  vroom::vroom(here::here(year, folder, "processed", "ages_yrs.csv"))$yrs -> yrs
  vroom::vroom(here::here(year, folder, "processed", "bio_rec_f.csv")) %>%
    tidytable::select.(sp_biom, recruits) %>%
    tidytable::bind_cols.(year = yrs) -> dat

  # set view
  ggplot2::theme_set(afscassess::theme_report())

  data.frame(ssb = dat$sp_biom[1:(length(yrs) - rec_age)] / 1000,
             rec = dat$recruits[(rec_age + 1):length(yrs)],
             year = yrs[1:(length(yrs)-rec_age)]) %>%
    tidytable::mutate.(label = stringr::str_sub(year, 3),
                       decade = (floor(year / 10) * 10)) %>%
    ggplot2::ggplot(ggplot2::aes(ssb, rec)) +
    ggplot2::geom_label(ggplot2::aes(label=label, color = decade),
                        label.size = 0, show.legend = FALSE, size = 3, family="Times", alpha = 0.85) +
    ggplot2::expand_limits(x = 0, y = 0) +
    scico::scale_color_scico(palette = "romaO") +
    ggplot2::xlab("SSB (kt)") +
    ggplot2::ylab("Recruitment (millions)")

  if(isTRUE(save)) {
    ggplot2::ggsave(here::here(year, folder, "figs", "recr-ssb.png"),
                    width = 6.5, height = 5.5, units = "in", dpi = 200)
  }
}

#' Plot selectivity at age
#'
#' @param year model year
#' @param folder folder model is in
#' @param save default is TRUE, saves fig to the folder the model is in
#' @return
#' @export plot_selex
#'
#' @examples plot_selex(year, folder)
plot_selex <- function(year, folder, save=TRUE){

  if (!dir.exists(here::here(year, folder, "processed"))){
    stop("must run 'process_results' before creating figures")
  }


  vroom::vroom(here::here(year, folder, "processed", "selex.csv")) %>%
    tidytable::select.(age, Fishery = fish, Survey = srv1, Maturity = maturity) -> dat

  # set view
  ggplot2::theme_set(afscassess::theme_report())

  dat %>%
    tidytable::pivot_longer.(-age) %>%
    ggplot2::ggplot(ggplot2::aes(age, value, color = name)) +
    ggplot2::geom_line() +
    scico::scale_color_scico_d(name = "", palette = 'roma', begin = 0.2)+
    ggplot2::ylab("Selectivity/Maturity\n") +
    afscassess::scale_x_tickr(name = "Year", data=dat, var=age, to=10, start = 0) +
    ggplot2::theme(legend.justification=c(1,0),
                   legend.position=c(0.9,0.2))

  if(isTRUE(save)) {
    ggplot2::ggsave(here::here(year, folder, "figs", "selex.png"),
                    width = 6.5, height = 5.5, units = "in", dpi = 200)
  }
}

#' Plot survey biomass
#'
#' @param year model year
#' @param folder  folder model is in
#' @param save default is TRUE, saves fig to the folder the model is in
#'
#' @return
#' @export plot_survey
#'
#' @examples plot_survey(year, folder)
plot_survey <- function(year, folder, save=TRUE){

  if (!dir.exists(here::here(year, folder, "processed"))){
    stop("must run 'process_results' before creating figures")
  }

  vroom::vroom(here::here(year, folder, "processed", "survey.csv")) %>%
    tidytable::rename_with.(tolower) %>%
    tidytable::select.(year = starts_with("y"),
                       Observed = starts_with("bio"),
                       Predicted = pred,
                       se, lci, uci) %>%
    tidytable::pivot_longer.(-c(year, se, uci, lci)) %>%
    tidytable::mutate.(value = value / 1000,
                       uci = uci / 1000,
                       lci = lci / 1000) -> dat
  # set view
  ggplot2::theme_set(afscassess::theme_report())

  dat %>%
    ggplot2::ggplot(ggplot2::aes(year, value, color = name, linetype = name)) +
    ggplot2::geom_point() +
    ggplot2::geom_line() +
    ggplot2::geom_errorbar(ggplot2::aes(ymin = lci, ymax = uci), color = "gray", width = 0.4) +
    ggplot2::scale_color_grey(name = "",
                              start = 0.8, end = 0.2) +
    ggplot2::scale_linetype_manual(name = "",
                                   values = c(0,1)) +
    afscassess::scale_x_tickr(data=dat, var=year) +
    ggplot2::scale_y_continuous(labels = scales::comma) +
    ggplot2::xlab("Year") +
    ggplot2::ylab("Survey biomass (kt)") +
    ggplot2::expand_limits(y = 0) +
    ggplot2::guides(color = ggplot2::guide_legend(override.aes = list(linetype = 0))) +
    ggplot2::theme(legend.justification=c(1,0),
                   legend.position=c(0.98,0.8))

  if(isTRUE(save)) {
    ggplot2::ggsave(here::here(year, folder, "figs", "srv1_biomass.png"),
           width = 6.5, height = 6.5, units = "in", dpi = 200)
  }
}

' Swath plot
#'
#' @param year assessment year
#' @param folder   the folder with the model in it
#' @param save default is TRUE, saves fig to the folder the model is in
#'
#' @return
#' @export plot_swath
#'
#' @examples plot_swath(year, folder)
plot_swath <- function(year, folder, save=TRUE){

  yr = year
  if (!dir.exists(here::here(year, folder, "processed"))){
    stop("must run 'process_results' before creating figures")
  }

  # set view
  ggplot2::theme_set(afscassess::theme_report())

  # establish quantiles
  q_name <- tidytable::map_chr.(seq(.025,.975,.05), ~ paste0("q", .x*100))
  q_fun <- tidytable::map.(seq(.025,.975,.05), ~ purrr::partial(quantile, probs = .x, na.rm = TRUE)) %>%
    purrr::set_names(nm = q_name)

  # read in data calculate quantiles/median and plot
  vroom::vroom(here::here(year, folder, "processed", "ages_yrs.csv"))$yrs -> yrs
  vroom::vroom(here::here(year, folder, "processed", "b35_b40_yld.csv")) -> bby

  vroom::vroom(here::here(year, folder, "processed", "mceval.csv")) %>%
    tidytable::select.(paste0("spawn_biom_", yrs),
                       paste0("spawn_biom_proj_", (max(yrs)+1):(max(yrs+15)))) %>%
    tidytable::mutate.(group = 1:tidytable::n.()) %>%
    tidytable::pivot_longer.(c(-group), values_to = "biomass") %>%
    tidytable::mutate.(year = as.numeric(stringr::str_extract(name, "[[:digit:]]+")),
                       biomass = biomass / 1000) %>%
    dplyr::group_by(year) %>%
    dplyr::summarise_at(dplyr::vars(biomass), tibble::lst(!!!q_fun, median)) %>%
    tidytable::pivot_longer.(-c(year, median)) %>%
    tidytable::mutate.(grouping = tidytable::case_when.(name == q_name[1] | name == q_name[20] ~ 1,
                                                        name == q_name[2] | name == q_name[19] ~ 2,
                                                        name == q_name[3] | name == q_name[18] ~ 3,
                                                        name == q_name[4] | name == q_name[17] ~ 4,
                                                        name == q_name[5] | name == q_name[16] ~ 5,
                                                        name == q_name[6] | name == q_name[15] ~ 6,
                                                        name == q_name[7] | name == q_name[14] ~ 7,
                                                        name == q_name[8] | name == q_name[13] ~ 8,
                                                        name == q_name[9] | name == q_name[12] ~ 9,
                                                        name == q_name[10] | name == q_name[11] ~ 10)) %>%
    tidytable::mutate.(min = min(value),
                       max = max(value),
                       fill = as.character(ifelse(year>yr, 0, 1)),
                       .by = c(year, grouping)) -> dat


  dat %>%
    ggplot2::ggplot(ggplot2::aes(year, group = grouping)) +
    ggplot2::geom_ribbon(data = dplyr::filter(dat, year<yr),
                         ggplot2::aes(ymin = min, ymax = max), alpha = 0.07) +
    ggplot2::geom_line(ggplot2::aes(y = median)) +
    ggplot2::guides(linetype = "none", fill = "none") +
    ggplot2::geom_ribbon(data = dat, ggplot2::aes(ymin = min, ymax = max), alpha = 0.07) +
    afscassess::scale_x_tickr(data=dat, var=year, to=10, start = 1960) +
    ggplot2::ylab("Spawning biomass (kt)") +
    ggplot2::xlab("Year") +
    ggplot2::expand_limits(y = 0) +
    ggplot2::geom_hline(yintercept = c(bby$B40/1000, bby$B35/1000),
                        lty = c(3, 1))

  if(isTRUE(save)) {
    ggplot2::ggsave(here::here(year, folder, "figs", "swath.png"),
                    width = 6.5, height = 5.5, units = "in", dpi = 200)
  }
}

' Parameter histogram plots
#'
#' @param year assessment year
#' @param folder   the folder with the model in it
#' @param model_name the name of your .tpl file
#' @param pars parameter names
#' @param save default is TRUE, saves fig to the folder the model is in
#'
#' @return
#' @export plot_params
#'
#' @examples plot_swath(year, folder)
plot_params <- function(year, folder, model_name, pars = c("q_srv1", "ABC", "nattymort", "tot_biom", "F40","spawn_biom"), save=TRUE) {
  if (!dir.exists(here::here(year, folder, "processed"))){
    stop("must run 'process_results()' before creating figures")
  }

  # set view
  ggplot2::theme_set(afscassess::theme_report())

  vroom::vroom(here::here(year, folder, "processed", "ages_yrs.csv"))$yrs -> yrs

  read.delim(here::here(year, folder, paste0(model_name, ".std")), sep="", header = TRUE) %>%
    tidytable::filter(name %in% pars) %>%
    tidytable::slice(tidytable::n.(),
                     .by = name) %>%
    tidytable::mutate(name = tidytable::case_when(name =="q_srv1" ~ "q_srv",
                                                  name =="nattymort" ~ "natmort",
                                                  name == "F40" ~ "F",
                                                  name == "spawn_biom" ~ "spawn_biom_",
                                                  name == "tot_biom" ~ "tot_biom_",
                                                  TRUE ~ name),
                      value = tidytable::case_when(name == "ABC" ~ value / 1000,
                                                   name == "spawn_biom_" ~ value / 1000,
                                                   name == "tot_biom_" ~ value / 1000,
                                                   TRUE ~ value),
                      name = factor(name, levels = c("q_srv", "natmort", "F", "ABC", "tot_biom_", "spawn_biom_"))) -> fits

  vroom::vroom(here::here(year, folder, "processed", "mceval.csv"))  %>%
    tidytable::select.(q_srv1, ABC, natmort, paste0("tot_biom_", yrs),
                       F40, paste0("spawn_biom_", yrs)) %>%
    tidytable::mutate.(group = 1:dplyr::n()) %>%
    tidytable::pivot_longer.(-group) %>%
    tidytable::mutate.(years = as.numeric(gsub('\\D+','', name)),
                       name = gsub('[[:digit:]]+', '', name),
                       value = tidytable::case_when.(name=="spawn_biom_" ~ value / 1000,
                                                     name=="tot_biom_" ~ value / 1000,
                                                     name=="ABC" ~ value / 1000,
                                                     TRUE ~ value),
                       name = factor(name, levels = c("q_srv", "natmort", "F", "ABC", "tot_biom_", "spawn_biom_"))) -> dat

  p1 = dat %>%
    tidytable::filter.(name == "q_srv") %>%
    ggplot2::ggplot(ggplot2::aes(value)) +
    # facet_wrap(~name, scales = "free", dir = "v") +
    ggplot2::geom_histogram(ggplot2::aes(y=(..count..)/tapply(..count..,..PANEL..,sum)[..PANEL..]),
                            fill = "lightgray", color = "darkgray", bins = 50) +
    # scico::scale_fill_scico(palette = "grayC", direction = -1) +
    # ggplot2::scale_x_continuous(breaks = seq(0,2.5,0.5)) +
    ggplot2::geom_segment(data = tidytable::filter.(fits, name == "q_srv"),
                          mapping = ggplot2::aes(x = value, xend = value, y = 0, yend = 0.05),
                          size = 2, color = "black") +
    ggplot2::theme(axis.title.y = ggplot2::element_blank()) +
    ggplot2::xlab(expression("Trawl survey catchability ("*italic(q)*")"))

  p2 = dat %>%
    tidytable::filter.(name == "natmort") %>%
    ggplot2::ggplot(ggplot2::aes(value)) +
    # facet_wrap(~name, scales = "free", dir = "v") +
    ggplot2::geom_histogram(ggplot2::aes(y=(..count..)/tapply(..count..,..PANEL..,sum)[..PANEL..]),
                            fill = "lightgray", color = "darkgray", bins = 50) +
    # scico::scale_fill_scico(palette = "grayC", direction = -1) +
    # ggplot2::scale_x_continuous(breaks = seq(0,2.5,0.5)) +
    ggplot2::geom_segment(data = tidytable::filter.(fits, name == "natmort"),
                          mapping = ggplot2::aes(x = value, xend = value, y = 0, yend = 0.05), size = 2, color = "black") +
    ggplot2::ylab("Probability density") +
    ggplot2::xlab(expression("Natural mortality ("*italic(M)*")")) +
    ggplot2::theme(legend.position = "none")

  p3 = dat %>%
    tidytable::filter.(name == "F") %>%
    ggplot2::ggplot(ggplot2::aes(value)) +
    # facet_wrap(~name, scales = "free", dir = "v") +
    ggplot2::geom_histogram(ggplot2::aes(y=(..count..)/tapply(..count..,..PANEL..,sum)[..PANEL..]),
                            fill = "lightgray", color = "darkgray", bins = 50) +
    # scico::scale_fill_scico(palette = "grayC", direction = -1) +
    # ggplot2::scale_x_continuous(breaks = seq(0,2.5,0.5)) +
    ggplot2::geom_segment(data = tidytable::filter.(fits, name == "F"),
                          mapping = ggplot2::aes(x = value, xend = value, y = 0, yend = 0.05),
                          size = 2, color = "black") +
    ggplot2::theme(axis.title.y = ggplot2::element_blank()) +
    ggplot2::xlab(expression(italic(F)["40%"])) +
    ggplot2::theme(legend.position = "none")

  p4 = dat %>%
    tidytable::filter.(name == "ABC") %>%
    ggplot2::ggplot(ggplot2::aes(value)) +
    # facet_wrap(~name, scales = "free", dir = "v") +
    ggplot2::geom_histogram(ggplot2::aes(y=(..count..)/tapply(..count..,..PANEL..,sum)[..PANEL..]),
                            fill = "lightgray", color = "darkgray", bins = 50) +
    # scico::scale_fill_scico(palette = "grayC", direction = -1) +
    # ggplot2::scale_x_continuous(breaks = seq(0,2.5,0.5)) +
    ggplot2::geom_segment(data = tidytable::filter.(fits, name == "ABC"),
                          mapping = ggplot2::aes(x = value, xend = value, y = 0, yend = 0.05),
                          size = 2, color = "black") +
    ggplot2::theme(axis.title.y = ggplot2::element_blank()) +
    ggplot2::xlab("ABC (kt)") +
    ggplot2::theme(legend.position = "none")


  p5 = dat %>%
    tidytable::filter.(name == "tot_biom_", years == year) %>%
    ggplot2::ggplot(ggplot2::aes(value)) +
    # facet_wrap(~name, scales = "free", dir = "v") +
    ggplot2::geom_histogram(ggplot2::aes(y=(..count..)/tapply(..count..,..PANEL..,sum)[..PANEL..]),
                            fill = "lightgray", color = "darkgray", bins = 50) +
    # scico::scale_fill_scico(palette = "grayC", direction = -1) +
    # ggplot2::scale_x_continuous(breaks = seq(0,2.5,0.5)) +
    ggplot2::geom_segment(data = tidytable::filter.(fits, name == "tot_biom_"),
                          mapping = ggplot2::aes(x = value, xend = value, y = 0, yend = 0.05),
                          size = 2, color = "black") +
    ggplot2::theme(axis.title.y = ggplot2::element_blank()) +
    ggplot2::xlab("Current total biomass (kt)") +
    ggplot2::theme(legend.position = "none")

  p6 = dat %>%
    tidytable::filter.(name == "spawn_biom_", years == year) %>%
    ggplot2::ggplot(ggplot2::aes(value)) +
    # facet_wrap(~name, scales = "free", dir = "v") +
    ggplot2::geom_histogram(ggplot2::aes(y=(..count..)/tapply(..count..,..PANEL..,sum)[..PANEL..]),
                            fill = "lightgray", color = "darkgray", bins = 50) +
    # scico::scale_fill_scico(palette = "grayC", direction = -1) +
    # ggplot2::scale_x_continuous(breaks = seq(0,2.5,0.5)) +
    ggplot2::geom_segment(data = tidytable::filter.(fits, name == "spawn_biom_"),
                          mapping = ggplot2::aes(x = value, xend = value, y = 0, yend = 0.05),
                          size = 2, color = "black") +
    ggplot2::theme(axis.title.y = ggplot2::element_blank()) +
    ggplot2::xlab("Current spawning biomass (kt)") +
    ggplot2::theme(legend.position = "none")

  print((cowplot::plot_grid(p1, p4, p2,  p5, p3, p6, align = "v", ncol = 2, rel_heights = c(0.5, 0.5))))

  if(isTRUE(save)) {
    png(filename=here::here(year, folder, "figs", "hists.png"), width = 6.5, height = 6.5,
        units = "in", type ="cairo", res = 200)

    print((cowplot::plot_grid(p1, p4, p2,  p5, p3, p6, align = "v", ncol = 2, rel_heights = c(0.5, 0.5))))

    dev.off()

  }

}


' retrospective plots
#'
#' @param year assessment year
#' @param folder   the folder with the model in it
#' @param n_retro = 10 number of years to retro
#' @param save default is TRUE, saves fig to the folder the model is in
#'
#' @return
#' @export plot_retro
#'
#' @examples plot_retro(year, folder)
plot_retro <- function(year, folder, n_retro=10, save=TRUE) {
  peels = n_retro - 1
  max_year = year
  # loop through mcmc output
  age_yr = read.csv(here::here(year, folder, "processed", "ages_yrs.csv"))
  yrs = age_yr %>%
    dplyr::select(yrs) %>%
    tidyr::drop_na() %>%
    dplyr::pull(yrs)
  styr_rec = age_yr[1,3]
  retro_yrs = (year - n_retro + 1):year

  dat = list()

  for(i in 1:n_retro) {

    read.delim(here::here(year, folder, "retro", "results",
                          paste0("mcmc_", retro_yrs[i], ".std")),
               sep = "",  header = FALSE) -> df

    df = df[floor((0.2 * nrow(df))):nrow(df),] # drop burn in

    nyrs_tot_biom_proj = ncol(df)- (2*length(yrs[1]:retro_yrs[i]) + 7 +
      length(styr_rec:retro_yrs[i]) +  2*length(max(retro_yrs[i]) + 1:15) +
      length(max(retro_yrs[i]) + 1:10))-1

    colnames(df) = c("sigr", "q_srv1", "q_srv2", "F40", "natmort",
                     "ABC", "obj_fun",
                     paste0("tot_biom_", yrs[1]:retro_yrs[i]),
                     paste0("log_rec_dev_", styr_rec:retro_yrs[i]),
                     paste0("spawn_biom_", yrs[1]:retro_yrs[i]),
                     "log_mean_rec",
                     # "spawn_biom_proj",
                     paste0("spawn_biom_proj_", max(retro_yrs[i]) + 1:15),
                     paste0("pred_catch_proj_", max(retro_yrs[i]) + 1:15),
                     paste0("rec_proj_", max(retro_yrs[i]) + 1:10),
                     paste0("tot_biom_proj_", max(retro_yrs[i]) + 1:nyrs_tot_biom_proj))

    dat[[i]] = df %>% dplyr::mutate(retro_year = retro_yrs[i])

  }

  # clean up columns so can bind all together
  col <- unique(unlist(sapply(dat, names)))
  dat <- lapply(dat, function(df) {
    df[, setdiff(col, names(df))] <- NA
    df
  })

  do.call(rbind, dat)  -> retro_mc

  # save output
  write.csv(retro_mc, here::here(year, folder, "processed", "retro_mcmc.csv"), row.names = FALSE)

  # functions for quantiles
  q_name <- purrr::map_chr(c(.025,.975), ~ paste0("q", .x*100))
  q_fun <- purrr::map(c(.025,.975), ~ purrr::partial(quantile, probs = .x, na.rm = TRUE)) %>%
    purrr::set_names(nm = q_name)

  retro_mc %>%
    dplyr::select(paste0("spawn_biom_", yrs), retro_year) %>%
    tidyr::pivot_longer(c(-retro_year), values_to = "biomass") %>%
    dplyr::mutate(year = as.numeric(stringr::str_extract(name, "[[:digit:]]+")),
                  biomass = biomass / 1000) %>%
    dplyr::group_by(year, retro_year) %>%
    dplyr::summarise_at(dplyr::vars(biomass), tibble::lst(!!!q_fun, median)) %>%
    dplyr::mutate(Retro = factor(retro_year)) %>%
    dplyr::ungroup() -> dat


  dat %>%
    dplyr::select(year, retro_year, median) %>%
    dplyr::group_by(year) %>%
    dplyr::mutate(pdiff = (median - median[retro_year==max_year]) /
                    median[retro_year==max_year]) %>%
    tidyr::drop_na() %>%
    dplyr::filter(year %in% (max_year-peels):max_year) %>%
    dplyr::ungroup() %>%
    dplyr::filter(year == retro_year, year !=max_year) %>%
    dplyr::summarise(rho = mean(pdiff)) %>%
    dplyr::pull(rho) -> ssb_rho



  # set view
  ggplot2::theme_set(afscassess::theme_report())

  dat %>%
    # filter(retro_year==2022) %>%
    ggplot2::ggplot(ggplot2::aes(year, median, color = Retro, group = Retro)) +
    ggplot2::geom_line() +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = q2.5, ymax = q97.5, fill = Retro),
                         alpha = .05, color = NA) +
    ggplot2::ylab("Spawning biomass (kt)\n") +
    ggplot2::xlab("\nYear") +
    ggplot2::expand_limits(y = 0) +
    scico::scale_fill_scico_d(palette = "roma") +
    scico::scale_color_scico_d(palette = "roma") +
    funcr::theme_report() +
    ggplot2::scale_x_continuous(breaks = afscassess::tickr(dat, year, 10, start = 1960)$breaks,
                                labels = afscassess::tickr(dat, year, 10, start = 1960)$labels) +
    ggplot2::annotate(geom = "text", x=1963, y=Inf, hjust = -0.05, vjust = 2,
                      label = paste0("Mohn's rho = ", round(ssb_rho, 3)),
                      family = "Times") +
    ggplot2::theme(legend.position = "none") -> p1


  retro_mc %>%
    dplyr::select(paste0("spawn_biom_", yrs), retro_year) %>%
    tidyr::pivot_longer(c(-retro_year), values_to = "biomass") %>%
    dplyr::mutate(year = as.numeric(stringr::str_extract(name, "[[:digit:]]+")),
                  biomass = biomass / 1000,
                  pdiff = (biomass - biomass[retro_year==max_year]) /
                    biomass[retro_year==max_year])  %>%
    dplyr::group_by(year, retro_year) %>%
    dplyr::summarise_at(dplyr::vars(pdiff), tibble::lst(!!!q_fun, median)) %>%
    dplyr::mutate(Retro = factor(retro_year)) %>%
    dplyr::ungroup() -> dat2


  dat2 %>%
    ggplot2::ggplot(ggplot2::aes(year, median, color = Retro, group = Retro)) +
    ggplot2::geom_line() +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = q2.5, ymax = q97.5, fill = Retro),
                         alpha = .05, color = NA) +
    ggplot2::ylab(paste0("Percent change from ", year, "\n")) +
    ggplot2::xlab("\nYear") +
    ggplot2::expand_limits(y = 0) +
    scico::scale_fill_scico_d(palette = "roma") +
    scico::scale_color_scico_d(palette = "roma") +
    funcr::theme_report() +
    ggplot2::scale_x_continuous(breaks = funcr::tickr(dat, year, 10, start = 1960)$breaks,
                                labels = funcr::tickr(dat, year, 10, start = 1960)$labels) +
    ggplot2::theme(legend.position = "none") -> p2

  png(filename=here::here(year, folder, "figs", "retro.png"), width = 6.5, height = 6.5,
      units = "in", type ="cairo", res = 200)

  print((cowplot::plot_grid(p1, p2, ncol = 1)))

  dev.off()
}


' Base plots
#'
#' @param year assessment year
#' @param folder   the folder with the model in it
#' @param model_name the name of your .tpl file
#' @param rec_age recruitment age
#' @export base_plots
base_plots <- function(year, folder, model_name, rec_age) {
  if (!dir.exists(here::here(year, folder, "processed"))){
    stop("must run 'process_results()' before creating figures")
  }
  plot_catch(year, folder)
  plot_survey(year, folder)
  plot_selex(year, folder)
  plot_params(year, folder, model_name)
  plot_phase(year, folder, model_name)
  plot_rec(year, folder)
  plot_rec_ssb(year, folder, rec_age)
  plot_swath(year, folder)
  plot_comps(year, folder)
}

plot_F <- function(year, model) {


  read.csv(here::here(year, model, "processed", "bio_rec_f.csv")) %>%
    dplyr::select(year, F) -> dat

  png(filename=here::here(year, model, "figs", "fF.png"), width = 6.5, height = 6.5,
      units = "in", type ="cairo", res = 200)

  dat %>%
    ggplot2::ggplot(ggplot2::aes(year, F)) +
    ggplot2::geom_line() +
    ggplot2::expand_limits(y=0) +
    ggplot2::ylab("Fishing mortality rate (F)\n") +
    afscassess::scale_x_tickr(data=dat, var=year, to=10, start = 0)

  dev.off()
}


comp <- function(data, variable, fleet) {
  var <- rlang::enquo(variable)

  data %>%
    dplyr::mutate(level = factor(.data[[var]]),
                  Year = factor(Year)) %>%
    dplyr::filter(groups == "obs") %>%
    ggplot2::ggplot(ggplot2::aes(.data[[var]], value)) +
    ggplot2::geom_col(ggplot2::aes(fill = level), width = 1, color = 1) +
    ggplot2::facet_wrap(~Year, strip.position="right",
                        dir = "v",
                        ncol = 1) +
    ggplot2::geom_point(data = dplyr::filter(data, groups == "pred"),
                        ggplot2::aes(.data[[var]], value, group = 1)) +
    ggplot2::theme(panel.spacing.y = grid::unit(0, "mm")) +
    ggplot2::theme(axis.text.y = ggplot2::element_blank(),
                   axis.ticks.y = ggplot2::element_blank()) +
    ggplot2::xlab(Hmisc::capitalize(variable)) +
    ggplot2::ylab(paste0(Hmisc::capitalize(fleet)," ", variable,  " composition")) +
    afscassess::scale_x_tickr(data=data, var=.data[[var]], start = 0, min = 5) +
    ggplot2::theme(legend.position = "none")
}

