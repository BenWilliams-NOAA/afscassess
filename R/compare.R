#' Plot MCMC estimated biomass for multiple models
#'
#' @param year model year
#' @param models year and folder the models are in e.g., c('2020, plus_age', '2020, plus_length')
#' @param save default is TRUE, saves fig to the folder the model is in
#'
#' @return
#' @export plot_compare_biomass
plot_compare_biomass <- function(year, models = c('2020, plus_age', '2020, plus_length'), save=TRUE) {

  if (!dir.exists(here::here(year, "compare_models"))){
    dir.create(here::here(year, "compare_models"))
  }

  dat = data.frame()
  m = list(rep(NA, length(models)))
  for(i in 1:length(models)) {
    m[[i]] = scan(text = models[i], sep = ",", what = "")
    m[[i]][2] = gsub(" ", "", m[[i]][2])

    year = m[[i]][1]
    model = m[[i]][2]

    yrs = vroom::vroom(here::here(year, model, "processed", "ages_yrs.csv"))$yrs
    bio = vroom::vroom(here::here(year, model, "processed", "bio_rec_f.csv"))

    dat %>%
      tidytable::bind_rows.(
        vroom::vroom(here::here(year, model, "processed", "mceval.csv"))  %>%
          tidytable::select.(paste0("tot_biom_", yrs)) %>%
          tidytable::mutate.(group = 1:tidytable::n.()) %>%
          tidytable::pivot_longer.(-group) %>%
          tidytable::mutate.(year = as.numeric(gsub("tot_biom_", "", name)),
                             name = "Total biomass") %>%
          tidytable::bind_rows.( vroom::vroom(here::here(year, model, "processed", "mceval.csv")) %>%
                                   tidytable::select.(paste0("spawn_biom_", yrs)) %>%
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
          tidytable::mutate.(biomass = ifelse(name == "Total biomass", tot, bio),
                             model = model) %>%
          tidytable::select.(-tot, -bio)) -> dat
  }

  dat %>%
    ggplot2::ggplot(ggplot2::aes(year, biomass, color = model, fill = model)) +
    ggplot2::geom_line() +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = lci, ymax = uci), alpha = 0.1, color = NA) +
    ggplot2::facet_wrap(~name, dir = "v", scales = "free_y") +
    ggplot2::scale_y_continuous(name = "Biomass (kt)\n", labels = scales::comma) +
    ggplot2::expand_limits(y = 0) +
    afscassess::scale_x_tickr(data=dat, var=year, to=10, start = 1900) +
    scico::scale_color_scico_d(palette = 'roma') +
    scico::scale_fill_scico_d(palette = 'roma')

  if(isTRUE(save)) {
    ggplot2::ggsave(here::here(year, "compare_models", "figs", "compare_est_biomass.png"),
                    width = 6.5, height = 6.5, units = "in", dpi = 200)
  }
}
