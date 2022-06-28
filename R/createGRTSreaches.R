#' @title Create GRTS reaches
#'
#' @description Based on GRTS points and a stream layer, determines the boundaries of a reach that will be sampled as part of the study design. Strives to keep that reach on the same stream, and extend for a specified distance upstream and downstream of the GRTS point.
#'
#' @author Kevin See
#'
#' @param grts_pts `sf` object of GRTS points. Could be the output of `spsurvey::grts`.
#' @param strm_sf `sf` object containing the line layer of the GRTS frame. Should be in a projection with feet as the distant unit, and contain a column `llid` with unique codes for each stream.
#' @param length_buffer_mi distance (in miles) that the reach should extend on either side of the GRTS point.
#'
#' @import dplyr units sf lwgeom
#' @export
#' @return sf
#' @examples createGRTSreaches()

createGRTSreaches <- function(grts_pts = NULL,
                              strm_sf = NULL,
                              length_buffer_mi = 0.5) {

  stopifnot(!is.null(grts_pts),
            !is.null(strm_sf))

  pt_buffer_mi = length_buffer_mi * 2
  units(pt_buffer_mi) <- units::as_units("mi")


  # generate all pt buffers
  pt_buffs <- grts_pts %>%
    sf::st_buffer(dist = pt_buffer_mi) %>%
    select(siteID,
           wtrbdy_d,
           wtrbdy_n,
           llid)

  # get the stream within that buffer
  strm_nr_pt <- strm_sf %>%
    sf::st_intersection(pt_buffs) %>%
    filter(llid == llid.1) %>%
    group_by(siteID) %>%
    summarize() %>%
    mutate(lngth = sf::st_length(geometry)) %>%
    rowwise() %>%
    mutate(geom_type = sf::st_geometry_type(geometry)) %>%
    ungroup()

  if(sum(strm_nr_pt$geom_type == "MULTILINESTRING") > 0) {
    strm_nr_pt %>%
      filter(geom_type == "MULTILINESTRING") %>%
      mutate(across(geometry,
                    sf::st_line_merge)) %>%
      bind_rows(strm_nr_pt %>%
                  filter(geom_type != "MULTILINESTRING")) %>%
      ungroup() %>%
      arrange(siteID) %>%
      rowwise() %>%
      mutate(geom_type = sf::st_geometry_type(geometry)) %>%
      ungroup() -> strm_nr_pt
  }

  samp_rch <- strm_nr_pt %>%
    group_by(siteID) %>%
    nest() %>%
    rename(line_sf = data) %>%
    left_join(grts_pts %>%
                group_by(siteID) %>%
                nest() %>%
                rename(pt_sf = data),
              by = "siteID") %>%
    mutate(strm_split = map2(line_sf,
                             pt_sf,
                             .f = function(x, y) {
                               strm_split = x %>%
                                 select(geometry) %>%
                                 lwgeom::st_split(y %>%
                                                    sf::st_buffer(dist = .1)) %>%
                                 sf::st_collection_extract(type = "LINESTRING") %>%
                                 mutate(lngth = sf::st_length(geometry)) %>%
                                 filter(lngth > min(lngth)) %>%
                                 arrange(desc(lngth)) %>%
                                 slice(c(1:2)) %>%
                                 mutate(across(lngth,
                                               as.numeric)) %>%
                                 mutate(across(lngth,
                                               conv_unit,
                                               from = "ft",
                                               to = "mi")) %>%
                                 mutate(split_num = 1:n()) %>%
                                 mutate(meet_req = if_else(lngth >= length_buffer_mi, T, F),
                                        lngth_keep = if_else(sum(meet_req) == length(meet_req),
                                                             length_buffer_mi,
                                                             NA_real_)) %>%
                                 select(split_num,
                                        everything(),
                                        geometry)

                               if(sum(is.na(strm_split$lngth_keep)) > 0) {
                                 strm_split %<>%
                                   mutate(lngth_keep = if_else(!meet_req,
                                                               lngth,
                                                               length_buffer_mi + (length_buffer_mi - lngth[!meet_req])))
                               }
                               return(strm_split)
                             }),
           samp_rchs = map2(strm_split,
                            pt_sf,
                            .f = function(x, pt) {

                              x %>%
                                group_by(split_num) %>%
                                nest() %>%
                                ungroup() %>%
                                mutate(start_pt = list(pt)) %>%
                                mutate(rch_keep = map2(data,
                                                       start_pt,
                                                       .f = function(y, y_pt) {
                                                         subsampleReach(rch_sf = y,
                                                                        start_pt = y_pt,
                                                                        needed_length = y$lngth_keep)
                                                       }))
                            })) %>%
    select(siteID, samp_rchs) %>%
    unnest(samp_rchs) %>%
    select(siteID,
           split_num,
           rch_keep) %>%
    unnest(rch_keep) %>%
    sf::st_sf()

  grts_rchs <- samp_rch %>%
    group_by(siteID) %>%
    summarize() %>%
    left_join(grts_pts %>%
                sf::st_drop_geometry(),
              by = "siteID") %>%
    select(siteID:stratum,
           wtrbdy_d,
           wtrbdy_n,
           llid) %>%
    mutate(rch_lngth_mi = sf::st_length(geometry),
           across(rch_lngth_mi,
                  as.numeric),
           across(rch_lngth_mi,
                  conv_unit,
                  from = "ft",
                  to = "mi"))

  return(grts_rchs)
}
