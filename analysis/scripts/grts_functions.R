# Author: Kevin See
# Purpose: functions to test GRTS draws in the Newaukum
# Created: 6/27/22
# Last Modified: 6/27/2022
# Notes:

# #-----------------------------------------------------------------
# # load needed libraries
# library(sf)
# library(dplyr)
# library(sfnetworks)

#-----------------------------------------------------------------
# extract a sub-reach of a certain distance from a starting point
reach_subsample_length <- function(rch_sf = NULL,
                                   needed_length = 0.5,
                                   start_pt = NULL) {

  stopifnot(!is.null(rch_sf),
            !is.null(start_pt))

  # convert to a network object
  rch_net <- rch_sf %>%
    select(geometry) %>%
    st_subdivide(5) %>%
    st_collection_extract(type = "LINESTRING") %>%
    as_sfnetwork(directed = F) %>%
    activate("edges") %>%
    mutate(edge_id = 1:n(),
           w = edge_length()) %>%
    activate("nodes") %>%
    mutate(node_id = 1:n())

  # which is the start node?
  start_node = rch_net %>%
    activate("nodes") %>%
    st_as_sf() %>%
    mutate(start_dist = st_distance(., start_pt)[,1]) %>%
    filter(start_dist == min(start_dist))

  # calculate network paths
  paths <- st_network_paths(rch_net,
                            from = start_node,
                            weights = "w") %>%
    mutate(end_node = map_dbl(node_paths,
                              .f = function(x) {
                                x[length(x)]
                              }),
           path_length = map_dbl(edge_paths,
                                 .f = function(x) {
                                   rch_net %>%
                                     activate("edges") %>%
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
    activate("edges") %>%
    st_as_sf() %>%
    filter(from %in% keep_nodes |
             to %in% keep_nodes) %>%
    select(geometry) %>%
    # combine back to single multilinestring object
    st_combine() %>%
    st_sf() %>%
    mutate(lngth = st_length(geometry)) %>%
    mutate(across(lngth,
                  as.numeric),
           across(lngth,
                  conv_unit,
                  from = "ft",
                  to = "mi")) %>%
    select(lngth,
           geometry)

  return(rch_keep)
}

#-----------------------------------------------------------------
# take a sample of GRTS points and create a reach associated
# with each point (of a certain length)

create_grts_reaches_old <- function(grts_pts = NULL,
                                strm_sf = NULL,
                                length_buffer_mi = 0.5,
                                pt_buffer_mi = 1,
                                verbose = T) {

  stopifnot(!is.null(grts_pts),
            !is.null(strm_sf))

  units(pt_buffer_mi) <- units::as_units("mi")

  grts_rchs <- NULL
  for(i in 1:nrow(grts_pts)) {
    if(verbose) {
      cat(paste("Working on point", i, "of", nrow(grts_pts), "\n"))
    }

    # pull out single GRTS point
    ex_pt = grts_pts %>%
      slice(i) %>%
      select(siteID:stratum,
             wgt,
             wtrbdy_d,
             wtrbdy_n,
             llid,
             geometry)

    pt_buff <- ex_pt %>%
      st_buffer(dist = pt_buffer_mi) %>%
      select(wtrbdy_d,
             wtrbdy_n,
             llid)

    # get the stream within that buffer
    strm_nr_pt <- strm_sf %>%
      st_intersection(pt_buff) %>%
      filter(llid == llid.1) %>%
      group_by(llid) %>%
      summarize() %>%
      # st_cast(to = "LINESTRING") %>%
      # st_sf() %>%
      mutate(lngth = st_length(geometry))

    if(st_geometry_type(strm_nr_pt) == "MULTILINESTRING") {
      strm_nr_pt <- strm_nr_pt %>%
        st_line_merge()
    }

    # strm_nr_pt <- strm_sf %>%
    #   st_intersection(pt_buff) %>%
    #   filter(llid == llid.1) %>%
    #   st_cast(to = "LINESTRING") %>%
    #   mutate(lngth = st_length(geometry)) %>%
    #   filter(lngth == max(lngth)) %>%
    #   slice(1)

    # split the stream at the GRTS point
    strm_split <- strm_nr_pt %>%
      select(geometry) %>%
      st_split(ex_pt %>%
                 st_buffer(dist = .1)) %>%
      st_collection_extract(type = "LINESTRING") %>%
      mutate(lngth = st_length(geometry)) %>%
      filter(lngth > min(lngth)) %>%
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

    samp_rchs = vector("list",
                       length = nrow(strm_split))
    for(j in 1:nrow(strm_split)) {
      samp_rchs[[j]] = reach_subsample_length(rch_sf = strm_split %>%
                                                slice(j),
                                              needed_length = strm_split$lngth_keep[j],
                                              start_pt = ex_pt)
    }

    # put the upstream and downstream portions back together
    samp_rch = st_union(samp_rchs[[1]],
                        samp_rchs[[2]]) %>%
      select(geometry) #%>%
      # st_line_merge()
      # st_cast(to = "LINESTRING")

    # add some info from the GRTS point to the reach
    ex_rch = samp_rch %>%
      bind_cols(ex_pt %>%
                  st_drop_geometry() %>%
                  as_tibble()) %>%
      select(everything(),
             geometry)


    if(i == 1) {
      grts_rchs <- ex_rch
    } else {
      grts_rchs %<>%
        rbind(ex_rch)
    }

    rm(ex_pt,
       pt_buff,
       strm_nr_pt,
       strm_split,
       samp_rchs,
       samp_rch,
       ex_rch)
  }

  grts_rchs %<>%
    mutate(rch_lngth_mi = st_length(geometry),
           across(rch_lngth_mi,
                  as.numeric),
           across(rch_lngth_mi,
                  conv_unit,
                  from = "ft",
                  to = "mi"))

  return(grts_rchs)
}

#-----------------------------------------------------------------
# take a sample of GRTS points and create a reach associated
# with each point (of a certain length)
# trying to vectorize it more

create_grts_reaches <- function(grts_pts = NULL,
                                strm_sf = NULL,
                                length_buffer_mi = 0.5,
                                pt_buffer_mi = 1) {

  stopifnot(!is.null(grts_pts),
            !is.null(strm_sf))

  units(pt_buffer_mi) <- units::as_units("mi")


  # generate all pt buffers
  pt_buffs <- grts_pts %>%
    st_buffer(dist = pt_buffer_mi) %>%
    select(siteID,
           wtrbdy_d,
           wtrbdy_n,
           llid)

  # get the stream within that buffer
  strm_nr_pt <- strm_sf %>%
    st_intersection(pt_buffs) %>%
    filter(llid == llid.1) %>%
    group_by(siteID) %>%
    summarize() %>%
    mutate(lngth = st_length(geometry)) %>%
    rowwise() %>%
    mutate(geom_type = st_geometry_type(geometry)) %>%
    ungroup()

  if(sum(strm_nr_pt$geom_type == "MULTILINESTRING") > 0) {
    strm_nr_pt %>%
      filter(geom_type == "MULTILINESTRING") %>%
      mutate(across(geometry,
                    st_line_merge)) %>%
      bind_rows(strm_nr_pt %>%
                  filter(geom_type != "MULTILINESTRING")) %>%
      ungroup() %>%
      arrange(siteID) %>%
      rowwise() %>%
      mutate(geom_type = st_geometry_type(geometry)) %>%
      ungroup() -> strm_nr_pt
  }

  #
  # strm_nr_pt %>%
  #   filter(geom_type == "MULTILINESTRING") %>%
  #   group_by(siteID) %>%
  #   mutate(n_str = map_dbl(geometry,
  #                          .f = function(x) {
  #                            x %>%
  #                              st_cast(to = "LINESTRING") %>%
  #                              length()
  #                          }))


  # strm_nr_pt %>%
  #   select(geometry) %>%
  #   st_split(grts_pts %>%
  #              st_buffer(dist = .1))

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
                                st_split(y %>%
                                           st_buffer(dist = .1)) %>%
                                st_collection_extract(type = "LINESTRING") %>%
                                mutate(lngth = st_length(geometry)) %>%
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
                                                         reach_subsample_length(rch_sf = y,
                                                                                needed_length = y$lngth_keep,
                                                                                start_pt = y_pt)
                                                       }))
                            })) %>%
    select(siteID, samp_rchs) %>%
    unnest(samp_rchs) %>%
    select(siteID,
           split_num,
           rch_keep) %>%
    unnest(rch_keep) %>%
    st_sf()

  grts_rchs <- samp_rch %>%
    group_by(siteID) %>%
    summarize() %>%
    left_join(grts_pts %>%
                st_drop_geometry(),
              by = "siteID") %>%
    select(siteID:stratum,
           wtrbdy_d,
           wtrbdy_n,
           llid) %>%
    mutate(rch_lngth_mi = st_length(geometry),
           across(rch_lngth_mi,
                  as.numeric),
           across(rch_lngth_mi,
                  conv_unit,
                  from = "ft",
                  to = "mi"))

  return(grts_rchs)
}
