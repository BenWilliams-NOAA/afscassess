# code for generating data input plots for fishery stock assessments
#' Plot cumulative annual catch
#'
#' @param year model year
#' @param save default is TRUE, saves fig to the folder the model is in
#'
#' @export
#'
#' @examples
#' \dontrun{
#' plot_cum_catch(year=2026)
#' }
plot_cum_catch <- function(year, save = TRUE) {
  # set view
  ggplot2::theme_set(afscassess::theme_report())

  vroom::vroom(here::here(year, 'data', 'raw', 'fish_catch_data.csv')) %>%
    tidytable::mutate(day = lubridate::yday(week_end_date)) %>%
    tidytable::summarise(catch = sum(weight_posted),
                         .by = c(year, day)) %>%
    tidytable::mutate(cumsum = cumsum(catch),
                      .by = year) -> df

  df %>%
    tidytable::filter(year==max(year)) %>%
    tidytable::filter(day==max(day))  %>%
    tidytable::mutate(date = as.Date(paste(year, "01", "01", sep = "-")) + day) -> pt

  df %>%
    ggplot2::ggplot(ggplot2::aes(day, cumsum, color = year, group = year)) +
    ggplot2::geom_line() +
    ggplot2::geom_point(data=pt, size = 2) +
    scico::scale_color_scico("Year", palette = 'roma') +
    ggplot2::scale_y_continuous(labels=scales::comma) +
    ggplot2::xlab('Julian Day') +
    ggplot2::ylab('Cumulative catch (t)') +
    ggplot2::theme(legend.position.inside = c(x=0.2, y=0.8)) +
    ggplot2::ggtitle(paste("Catch as of", pt$date)) -> fig

  if(isTRUE(save)) {
    if (!dir.exists(here::here(year, "figs"))){
      dir.create(here::here(year, "figs"))
    }
    ggplot2::ggsave(plot = fig, filename = here::here(year, "figs", "cum_catch.png"),
                    width = 6.5, height = 6.5, units = "in", dpi = 200)
  }
  fig
}



#' Plot OFL, ABC, TAC and annual catch
#'
#' @param year model year
#' @param save default is TRUE, saves fig to the folder the model is in
#'
#' @export
#'
#' @examples
#' \dontrun{
#' plot_ofl_abc(year=2026)
#' }
plot_ofl_abc <- function(year, save = TRUE) {
  # set view
  ggplot2::theme_set(afscassess::theme_report())

  catch = vroom::vroom(here::here(year, "data", "output", "fish_catch.csv"))

  vroom::vroom(here::here(year, "data", "raw", "specs.csv")) %>%
    tidytable::filter(area_label == "GOA") %>%
    tidytable::select(year, OFL = overfishing_level, ABC = acceptable_biological_catch,  TAC = total_allowable_catch) %>%
    tidytable::left_join(catch) %>%
    tidytable::rename(Catch = catch) %>% 
    tidytable::pivot_longer(-year) %>%
    tidytable::mutate(name = factor(name, levels = c("OFL", "ABC", "TAC", "Catch"))) %>%
    ggplot2::ggplot(ggplot2::aes(year, value, color = name)) +
    ggplot2::geom_line() +
    geom_point() +
    ggplot2::scale_y_continuous(labels = scales::comma) +
    scico::scale_color_scico_d("Metric", palette='roma') +
    tickr::scale_x_tickr(data=catch, var=year) +
    ggplot2::xlab("Year") +
    ggplot2::ylab("Metric tons") +
    expand_limits(y = 0) -> fig

  if(isTRUE(save)) {
    if (!dir.exists(here::here(year, "figs"))){
      dir.create(here::here(year, "figs"))
    }
    ggplot2::ggsave(plot = fig, filename = here::here(year, "figs", "ofl_abc_tac_catch.png"),
                    width = 6.5, height = 6.5, units = "in", dpi = 200)
  }
  fig
}


#' Plot age comps for the fishery and survey that are input into the assessment
#'
#' @param year model year
#' @param save default is TRUE, saves fig to the folder the model is in
#'
#' @export
#'
#' @examples
#' \dontrun{
#' plot_age_comps_in(year=2026)
#' }
plot_age_comps_in <- function(year, save = TRUE) {
  ggplot2::theme_set(afscassess::theme_report())
  vroom::vroom(here::here(year, "data", "output", "fish_age_comp.csv"))[,-c(2:4)] %>% 
  as.data.frame() %>% 
  tidyr::pivot_longer(-year) %>% 
  mutate(source = "fishery") %>% 
bind_rows(
vroom::vroom(here::here(year, "data", "output", "goa_bts_age_comp.csv")) %>% 
  as.data.frame() %>% 
  tidyr::pivot_longer(-year) %>% 
  mutate(source = "survey")
) %>% 
  mutate(age = as.numeric(name)) -> dat

dat %>% 
  ggplot(aes(year, age, size = value, color = source)) +
  geom_point(alpha = 0.8) +
  scale_size_area() +
  scico::scale_color_scico_d(palette = 'roma', begin = 0.1, end = 0.8) +
  tickr::scale_x_tickr(data=dat, var = year) +
  tickr::scale_y_tickr(data=dat, var = age) +
  expand_limits(y = 0) +
  afscassess::theme_report() +
  xlab("\nYear") +
  ylab("Age\n") -> fig 
  
  if(isTRUE(save)) {
    if (!dir.exists(here::here(year, "figs"))){
      dir.create(here::here(year, "figs"))
    }
    ggplot2::ggsave(plot = fig, filename = here::here(year, "figs", "age_comp_in.png"),
                    width = 6.5, height = 6.5, units = "in", dpi = 200)
  }
  fig

}


#' Plot size comps for the fishery and survey that are input into the assessment
#'
#' @param year model year
#' @param save default is TRUE, saves fig to the folder the model is in
#'
#' @export
#'
#' @examples
#' \dontrun{
#' plot_size_comps_in(year=2026)
#' }
#' 
plot_size_comps_in <- function(year, save = TRUE) {
  ggplot2::theme_set(afscassess::theme_report())
  vroom::vroom(here::here(year, "data", "output", "fish_length_comp.csv"))[,-c(2:4)] %>% 
  as.data.frame() %>% 
  tidyr::pivot_longer(-year) %>% 
  mutate(source = "fishery") %>% 
  bind_rows(
    vroom::vroom(here::here(year, "data", "output", "goa_bts_sizecomp.csv"))[,-c(2:4)]  %>% 
    as.data.frame() %>% 
    tidyr::pivot_longer(-year) %>% 
    mutate(source = "survey")
  ) %>% 
    mutate(length = as.numeric(name)) -> dat

  dat %>% 
    ggplot(aes(year, length, size = value, color = source)) +
    geom_point(alpha = 0.7) +
    scale_size_area() +
    scico::scale_color_scico_d(palette = 'roma', begin = 0.1, end = 0.8) +
    tickr::scale_x_tickr(data=dat, var = year) +
    tickr::scale_y_tickr(data=dat, var = length) +
    afscassess::theme_report() +
    xlab("\nYear") +
    ylab("Length (cm)\n") -> fig 
  
  if(isTRUE(save)) {
    if (!dir.exists(here::here(year, "figs"))){
      dir.create(here::here(year, "figs"))
    }
    ggplot2::ggsave(plot = fig, filename = here::here(year, "figs", "size_comp_in.png"),
                    width = 6.5, height = 6.5, units = "in", dpi = 200)
  }
  fig

}