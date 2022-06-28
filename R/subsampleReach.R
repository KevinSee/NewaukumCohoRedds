#' @title Subsample Reach
#'
#' @description Starting from a fixed point, extracts a subsample of a reach in one direction from that point, extending it for a specified distance.
#'
#' @author Kevin See
#'
#' @param rch_sf `sf` object containing the line layer of a possible GRTS reach.
#' @param start_pt `sf` object containing the starting point (e.g. GRTS point) where reach will be built around
#' @param needed_length distance (in miles) that the reach should extend on one side of the GRTS point.
#'
#' @import dplyr units sf lwgeom sfnetworks purrr
#' @export
#' @return sf
#' @examples subsampleReach()

subsampleReach <- function(rch_sf = NULL,
                           start_pt = NULL,
                           needed_length = 0.5) {

  stopifnot(!is.null(rch_sf),
            !is.null(start_pt))

  # convert to a network object
  rch_net <- rch_sf %>%
    select(geometry) %>%
    lwgeom::st_subdivide(5) %>%
    sf::st_collection_extract(type = "LINESTRING") %>%
    sfnetworks::as_sfnetwork(directed = F) %>%
    sfnetworks::activate("edges") %>%
    mutate(edge_id = 1:n(),
           w = sfnetworks::edge_length()) %>%
    sfnetworks::activate("nodes") %>%
    mutate(node_id = 1:n())

  # which is the start node?
  start_node = rch_net %>%
    sfnetworks::activate("nodes") %>%
    sf::st_as_sf() %>%
    mutate(start_dist = sf::st_distance(., start_pt)[,1]) %>%
    filter(start_dist == min(start_dist))

  # calculate network paths
  paths <- sfnetworks::st_network_paths(rch_net,
                                        from = start_node,
                                        weights = "w") %>%
    mutate(end_node = purrr::map_dbl(node_paths,
                                     .f = function(x) {
                                       x[length(x)]
                                     }),
           path_length = purrr::map_dbl(edge_paths,
                                        .f = function(x) {
                                          rch_net %>%
                                            sfnetworks::activate("edges") %>%
                                            slice(x) %>%
                                            pull(w) %>%
                                            sum()
                                        })) %>%
    arrange(path_length) %>%
    mutate(across(path_length,
                  conv_unit,
                  from = "ft",
                  to = "mi")) %>%
    mutate(keep_seg = if_else(path_length <= needed_length,
                              T, F))
  # keep_seg = if_else(path_length == min(path_length[path_length > needed_length]),
  #                    T, keep_seg))

  # which nodes should be kept?
  keep_nodes <- paths %>%
    filter(keep_seg) %>%
    pull(end_node)

  # what edges are connected to those nodes?
  rch_keep <- rch_net %>%
    sfnetworks::activate("edges") %>%
    sf::st_as_sf() %>%
    filter(from %in% keep_nodes |
             to %in% keep_nodes) %>%
    select(geometry) %>%
    # combine back to single multilinestring object
    sf::st_combine() %>%
    sf::st_sf() %>%
    mutate(lngth = sf::st_length(geometry)) %>%
    mutate(across(lngth,
                  as.numeric),
           across(lngth,
                  units::conv_unit,
                  from = "ft",
                  to = "mi")) %>%
    select(lngth,
           geometry)

  return(rch_keep)
}
