#' Pipe operator
#'
#' See \code{magrittr::\link[magrittr:pipe]{\%>\%}} for details.
#'
#' @name %>%
#' @rdname pipe
#' @keywords internal
#' @export
#' @importFrom magrittr %>%
#' @usage lhs \%>\% rhs
NULL

#' helper functions for species codes
#' @export
#'
sp_switch <- function(species) {

  if(species == "NORK"){
    sp = 303
  }
  if(species == "DUSK"){
    sp = 330
  }
  if(species == "POP" | species == "POPA"){
    sp = 301
  }
  if(species == "REBS"){
    sp = c(307, 357)
  }
  sp
}

#' Setup folder structure and add .tpl's for a suite of biological functions
#'
#' Creates a common folder structure for assessment data
#'
#' @param year assessment year
#' @param dirs directories to write
#' @return creates a designated/named folder structure
#' @export setup_tpl
#'
#' @examples
#' \dontrun{
#' setup_tpl(2022)
#'}
setup_tpl <- function(year, dirs = "models"){

    folders = c("ageage", "allometric", "vonb", "wvonb", "length_sd")

    for(i in 1:length(dirs)){
      if(dir.exists(here::here(year, "data", dirs[i])) == FALSE){
        dir.create(here::here(year, "data", dirs[i]), recursive=TRUE)
      }
    }
    for(i in 1:length(folders)){
      dir.create(here::here(year, "data", "models", folders[i]))
    }

    file.copy(system.file("models", "ageage.tpl", package = "afscassess"),
              here::here(year, "data", "models", "ageage"))

    file.copy(system.file("models", "allometric.tpl", package = "afscassess"),
              here::here(year, "data", "models", "allometric"))

    file.copy(system.file("models", "vbl.tpl", package = "afscassess"),
              here::here(year, "data", "models", "vonb"))

    file.copy(system.file("models", "wvbl.tpl", package = "afscassess"),
              here::here(year, "data", "models", "wvonb"))

    file.copy(system.file("models", "wvbl.ctl", package = "afscassess"),
              here::here(year, "data", "models", "wvonb"))

    file.copy(system.file("models", "lengthsd.tpl", package = "afscassess"),
              here::here(year, "data", "models", "length_sd"))

  # dir.create(here::here(year, "proj"))
  # dir.create(here::here(year, "apportionment"))

}

# make data.table work
.datatable.aware <- TRUE

#' helper function for comp data
#'
#' obs  observed data from .rep file
#' pred predicted data from .rep file (if used)
#' rec_age  recruitement age
#' plus_age plus age group
#' comp age or length - default is length
#' lenbins set to base unless using alt in which case the file should be in the user_input folder and the name needs to be provided e.g., lengthbins.csv - the column must be named len_bin

#' @export rep_item
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

#' @export purrit

purrit <- function(obs, pred = NULL, rec_age, plus_age, comp = "length", lenbins){

  obs = stringr::str_split(obs, " ")

  purrr::map_if(obs[-1], is.character, as.numeric) %>%
    purrr::map_df(., ~as.data.frame(t(.))) %>%
    dplyr::select_if(~sum(!is.na(.)) > 0) -> obs


  if(comp == "age" & !is.null(pred)){
    obs %>%
      dplyr::select(1:(plus_age - rec_age + 2)) -> obs

    pred = stringr::str_split(pred, " ")
    purrr::map_if(pred[-1], is.character, as.numeric) %>%
      purrr::map_df(., ~as.data.frame(t(.))) %>%
      dplyr::select_if(~sum(!is.na(.)) > 0) %>%
      dplyr::select(1:(plus_age - rec_age + 2)) -> pred

    names(pred) <- names(obs) <- c("year", rec_age:plus_age)

    obs %>%
      tidyr::pivot_longer(-year, names_to = "age") %>%
      dplyr::mutate(groups = "obs") %>%
      dplyr::bind_rows(pred %>%
                         tidyr::pivot_longer(-year, names_to = "age") %>%
                         dplyr::mutate(groups = "pred")) %>%
      dplyr::mutate(age = as.integer(age),
                    Age = factor(age),
                    Year = factor(year)) -> dat


  } else if(comp != "age" & !is.null(pred)){

    obs %>%
      dplyr::select(1:(length(lenbins) + 1)) -> obs

    pred = stringr::str_split(pred, " ")
    purrr::map_if(pred[-1], is.character, as.numeric) %>%
      purrr::map_df(., ~as.data.frame(t(.))) %>%
      dplyr::select_if(~sum(!is.na(.)) > 0) %>%
      dplyr::select(1:(length(lenbins) + 1)) -> pred

    names(pred) <- names(obs) <- c("year", lenbins)

    obs %>%
      tidyr::pivot_longer(-year, names_to = "length") %>%
      dplyr::mutate(groups = "obs") %>%
      dplyr::bind_rows(pred %>%
                         tidyr::pivot_longer(-year, names_to = "length") %>%
                         dplyr::mutate(groups = "pred")) %>%
      dplyr::mutate(length = as.integer(length),
                    Length = factor(length),
                    Year = factor(year)) -> dat

  } else if(comp == "age" & is.null(pred)){
    obs %>%
      dplyr::select(1:(plus_age - rec_age + 2)) -> obs

    names(obs) <- c("year", rec_age:plus_age)

    obs %>%
      tidyr::pivot_longer(-year, names_to = "age") %>%
      dplyr::mutate(groups = "obs") %>%
      dplyr::mutate(age = as.integer(age),
                    Age = factor(age),
                    Year = factor(year)) -> dat

  } else if(comp != "age" & is.null(pred)){

    obs %>%
      dplyr::select(1:(length(lenbins) + 1)) -> obs

    names(obs) <- c("year", lenbins)

    obs %>%
      tidyr::pivot_longer(-year, names_to = "length") %>%
      dplyr::mutate(groups = "obs") %>%
      dplyr::mutate(length = as.integer(length),
                    Length = factor(length),
                    Year = factor(year)) -> dat
  }

  dat
}

#' helper function for creating dat file
#' @export collapse_row
#'
collapse_row <- function(data){

  l1 = paste(as.vector(data[1,]), collapse = " ")

  for(i in 2:nrow(data)){
    l2 = paste(as.vector(data[i,]), collapse = " ")
    l1 = c(l1, l2)
  }
  l1
}


#' Set figure theme for reports
#'
#' @param base_size font size
#' @param base_family font family
#'
#' @export
#'
#' @importFrom ggplot2 element_blank
#' @importFrom ggplot2 element_line
#' @importFrom ggplot2 element_rect
#' @importFrom ggplot2 element_text
#'
#'
#' @examples
#' \dontrun{
#'theme_report(base_size = 11, base_family = "Times")
#'
#'Other fonts are available, though sans font is
#'the easiest to implement using the following.
#'
#'theme_report(base_family = "")
#'
#'Updating font size is accomplished by changing the base_size.
#'
#'theme_report(base_size = 20, base_family = "")
#'}
theme_report <- function(base_size = 11, base_family = "Times") {

  windowsFonts(Times=windowsFont("TT Times New Roman"))

  half_line <- base_size/2

  ggplot2::theme_light(base_size = base_size, base_family = base_family) +
    ggplot2::theme(
      panel.grid.major = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
      axis.ticks.length = grid::unit(half_line / 2.2, "pt"),
      strip.background = ggplot2::element_rect(fill = NA, colour = NA),
      strip.text.x = ggplot2::element_text(colour = "black"),
      strip.text.y = ggplot2::element_text(colour = "black"),
      panel.border = ggplot2::element_rect(fill = NA),
      legend.key.size = grid::unit(0.9, "lines"),
      legend.key = ggplot2::element_rect(colour = NA, fill = NA),
      legend.background = ggplot2::element_rect(colour = NA, fill = NA)
    )
}
