# The file is used to generate polygons informed by genetic or outbreak data 

require( readr )
require( sf )
require( lubridate )
require( ggplot2 )
require( dplyr )
require( ks )
require( grDevices )


# FUNCTIONS ###

df_to_sf = function( point_df )
{
  # df with lon + lat
  # default crs = "EPSG:4326" 
  
  require( sf )
  
  
  lls = list()
  for( i in 1: nrow( point_df ) )
  {
    point.sfg = sf::st_point( x = c( point_df$lon[i], point_df$lat[i] ), dim = "XY" )
    
    lls[[i]] = point.sfg
  }
  
  points.sfc = sf::st_sfc( lls, crs = "EPSG:4326" )
  points_sf = sf::st_sf( data.frame( id           = seq(1, length(points.sfc) ), 
                                     geometry     = points.sfc ) ) 
  
  return( points_sf )
}

test_pos = function( in_polygon, points_sf )
{
  sf::sf_use_s2(use_s2 = FALSE)
  in_area_i = sf::st_contains( in_polygon, points_sf )[[1]]
  sf::sf_use_s2(use_s2 = TRUE)
  
  .df = data.frame( new_ = points_sf$new_outbreak,
                    pos  = FALSE,
                    year = lubridate::year(points_sf$sample_time) )
  
  .df$pos[ in_area_i ] = TRUE
  
  x.1 = length( which( .df$pos & .df$new_ ) )
  n.1 = length( which( .df$new_ ) )
  
  n.2 = .df %>% group_by( year ) %>%  dplyr::count() %>% dplyr::select(n) %>% pull(n)
  x.2 = .df %>% group_by( year ) %>%  dplyr::filter( pos == TRUE ) %>% dplyr::count() %>% pull(n)
  
  sf::sf_use_s2(use_s2 = FALSE)
  x.3 = units::set_units( sf::st_area( in_polygon ), "km^2" )  
  sf::sf_use_s2(use_s2 = TRUE)
  
  out_vec = c( x.1, n.1, x.2, n.2, x.3 )
  names(out_vec) = c( "new_pos", "new_sum", 
                      paste0("pos_", seq(2015,2022) ), as.character( seq(2015,2022) ),
                      "area")
  return(out_vec)
}
test_pos2 = function( in_polygon, points_sf )
{
  if( length(in_polygon) ==0 )
  {
    return(0)
  }else
  {
    sf::sf_use_s2(use_s2 = FALSE)
    in_area_i = sf::st_contains( in_polygon, points_sf )[[1]]
    sf::sf_use_s2(use_s2 = TRUE)
    
    .df = data.frame( new_ = points_sf$new_outbreak,
                      pos  = FALSE,
                      year = lubridate::year(points_sf$sample_time) )
    
    .df$pos[ in_area_i ] = TRUE
    
    x.1 = length( which( .df$pos & .df$new_ ) )
    n.1 = length( which( .df$new_ ) )
    
    
    null_df = data.frame(year=2015:2022)
    
    #n.2 = .df %>% group_by( year ) %>%  dplyr::count() %>% dplyr::select(n) %>% pull(n)
    
    n.2a = .df %>% group_by( year ) %>%  dplyr::count() 
    n.2b = left_join( null_df, n.2a, "year" )
    n.2b[ is.na(n.2b) ] = 0 
    
    n.2 = n.2b$n
    
    #x.2 = .df %>% group_by( year ) %>%  dplyr::filter( pos == TRUE ) %>% dplyr::count() %>% pull(n)
    if( !TRUE %in% .df$pos ){ x.2 = rep(0,length(2015:2022) ) }else
    {
      x.2a = .df %>% group_by( year ) %>% filter( pos ) %>% dplyr::count() 
      x.2b = left_join( null_df, x.2a, "year" )
      x.2b[ is.na(x.2b) ] = 0
      
      x.2 = x.2b$n
    }
    
    sf::sf_use_s2(use_s2 = FALSE)
    x.3 = units::set_units( sf::st_area( sf::st_intersection( in_polygon, X_county.t )  ), "km^2" )  
    sf::sf_use_s2(use_s2 = TRUE)
    
    out_vec = c( x.1, n.1, x.2, n.2, x.3 )
    names(out_vec) = c( "new_pos", "new_sum", 
                        paste0("pos_", seq(2015,2022) ), as.character( seq(2015,2022) ),
                        "area")
    return(out_vec)  
  }
  
  #20231012
}

trees_uniPGs = function( folder, treesN, phaseCutoff=2016.667, X_county, areaD=1000, silent=T )
{
  ls_node_sf = list()
  
  for( i in 1: treesN )
  {
    tab = read.csv( paste0( folder,"/TreeExtractions_", i ,".csv" ), header=T )
    
    startingNodeID = which( !tab[,"node1"] %in% tab[,"node2"] ) # root 
    
    tab_nootless = tab[ -startingNodeID, ]
    
    tab.1 = tab_nootless[ ,c("startLon", "startLat", "startYear") ]; names(tab.1) = c( "lon", "lat", "time" )
    tab.2 = tab_nootless[ ,c("endLon", "endLat", "endYear") ]; names(tab.2) = c( "lon", "lat", "time" )
    
    
    # 1 load all nodes 
    all_nodes = rbind( tab.1, tab.2 )
    
    # 1.1 filtering 
    all_nodes = subset(all_nodes, time >= phaseCutoff )
    
    # 2 unique nodes
    uniq_nodes = unique( all_nodes[ ,c("lon","lat") ] )
    
    # 3 nodes in X 
    nodes.sf = df_to_sf( uniq_nodes )
    
    nodes_X_i = st_contains( X_county, nodes.sf )[[1]]
    nodes_X.sf = nodes.sf[ nodes_X_i, ]
    
    
    # 4 area with diameter D
    nodes_X.sf.t = sf::st_transform( nodes_X.sf, crs = "EPSG:7801" )
    
    buffer_nodes_X.sf.t = sf::st_buffer( nodes_X.sf.t, dist = areaD )
    
    Union_nodes_X.sfc.t = sf::st_union( buffer_nodes_X.sf.t )
    
    ls_node_sf[[i]] = Union_nodes_X.sfc.t
    
    if( !silent ){ cat("tree", i, "\n") }
  }
  
  return( ls_node_sf )
}

get_node_num = function( folder, treesN, phaseCutoff=2016.667, X_county, silent=T )
{
  node_num = c()
  
  for( i in 1: treesN )
  {
    tab = read.csv( paste0( folder,"/TreeExtractions_", i ,".csv" ), header=T )
    
    startingNodeID = which( !tab[,"node1"] %in% tab[,"node2"] ) # root 
    
    tab_nootless = tab[ -startingNodeID, ]
    
    tab.1 = tab_nootless[ ,c("startLon", "startLat", "startYear") ]; names(tab.1) = c( "lon", "lat", "time" )
    tab.2 = tab_nootless[ ,c("endLon", "endLat", "endYear") ]; names(tab.2) = c( "lon", "lat", "time" )
    
    
    # 1 load all nodes 
    all_nodes = rbind( tab.1, tab.2 )
    
    # 1.1 filtering 
    all_nodes = subset(all_nodes, time >= phaseCutoff )
    
    # 2 unique nodes
    uniq_nodes = unique( all_nodes[ ,c("lon","lat") ] )
    
    # 3 nodes in X 
    nodes.sf = df_to_sf( uniq_nodes )
    
    nodes_X_i = st_contains( X_county, nodes.sf )[[1]]
    nodes_X.sf = nodes.sf[ nodes_X_i, ]
    
    
    # 4 area with diameter D
    # nodes_X.sf.t = sf::st_transform( nodes_X.sf, crs = "EPSG:7801" )
    # 
    # buffer_nodes_X.sf.t = sf::st_buffer( nodes_X.sf.t, dist = areaD )
    # 
    # Union_nodes_X.sfc.t = sf::st_union( buffer_nodes_X.sf.t )
    
    node_num = c( node_num, dim(nodes_X.sf)[1] )
    
    if( !silent ){ cat("tree", i, "\n") }
  }
  
  return( node_num )
  
}

site_uniPGs = function( vec_node, points_X_farm_t2.t, Nseed=666, areaD=1000, silent=T )
{
  set.seed( Nseed )
  
  ls_site_sf = list()
  
  Nsite_X_t2 = dim(points_X_farm_t2.t)[1]
  
  for( i in 1: length(vec_node) )
  {
    # tab = read.csv( paste0( folder,"/TreeExtractions_", i ,".csv" ), header=T )
    # startingNodeID = which( !tab[,"node1"] %in% tab[,"node2"] ) # root 
    # tab.1 = tab_nootless[ ,c("startLon", "startLat", "startYear") ]; names(tab.1) = c( "lon", "lat", "time" )
    # tab.2 = tab_nootless[ ,c("endLon", "endLat", "endYear") ]; names(tab.2) = c( "lon", "lat", "time" )
    
    # 1 get the nodes
    
    sub_sites_i = sample( x = seq(1, Nsite_X_t2), size = vec_node[i], replace = TRUE )
    sub_sites_i = unique( sub_sites_i )
    
    sub_sites.t = points_X_farm_t2.t[ sub_sites_i, ]
    
    
    # 2 area with diameter D
    
    buffer_sub_sites.t = sf::st_buffer( sub_sites.t, dist = areaD )
    
    Union_sites.sfc.t = sf::st_union( buffer_sub_sites.t )
    
    ls_site_sf[[i]] = Union_sites.sfc.t
    
    if( !silent ){ cat("tree", i, "\n") }
  }
  
  return( ls_site_sf )
}

null_uniPGs = function( vec_node, X_county.t, Nseed=666, areaD=1000, silent=T )
{
  set.seed( Nseed )
  
  ls_site_sf = list()
  
  for( i in 1: length(vec_node) )
  {
    # tab = read.csv( paste0( folder,"/TreeExtractions_", i ,".csv" ), header=T )
    # startingNodeID = which( !tab[,"node1"] %in% tab[,"node2"] ) # root 
    # tab.1 = tab_nootless[ ,c("startLon", "startLat", "startYear") ]; names(tab.1) = c( "lon", "lat", "time" )
    # tab.2 = tab_nootless[ ,c("endLon", "endLat", "endYear") ]; names(tab.2) = c( "lon", "lat", "time" )
    
    # 1 get the nodes

    sub_sites_i = sf::st_sample( x = X_county.t, size = vec_node[i] )
    
    # 2 area with diameter D
    
    buffer_sub_sites.t = sf::st_buffer( sub_sites_i, dist = areaD )
    
    Union_sites.sfc.t = sf::st_union( buffer_sub_sites.t )
    
    ls_site_sf[[i]] = Union_sites.sfc.t
    
    if( !silent ){ cat("tree", i, "\n") }
  }
  
  return(ls_site_sf)
}



# 0 BACKGROUND DATA ####

# 1. BAPHIQ website records
# 2. Taiwan county shape file

f_twai_page  = "data/baphiq_website.tsv"
f_county_shp = "data/predict/COUNTY_MOI_1090820.shp"


phase_cutoff    = 2016.667
latest_sample   = 2019.227
latest_sample_s = "2019-03-25"

# BAPHIQ website 
twai_page  = read_tsv( f_twai_page, show_col_types = F )


# shp file 
tw_county0 = sf::read_sf( f_county_shp )

tw_county   = sf::st_transform( tw_county0, crs = "EPSG:4326" )
tw_county.t = sf::st_transform( tw_county0, crs = "EPSG:7801" )



# 1 OUTBREAK SITES ####

lls = list()
for( i in 1: nrow( twai_page ) )
{
  point.sfg = sf::st_point( x = c( twai_page$Long[i], twai_page$Lat[i] ), dim = "XY" )
  
  lls[[i]] = point.sfg
}

points.sfc = sf::st_sfc( lls, crs = "EPSG:4326" )

points_sf = sf::st_sf( data.frame( id           = twai_page$no, 
                                   sample_time  = twai_page$Sampling_date, 
                                   confirm_time = twai_page$Confirming_date,
                                   species      = twai_page$Species,
                                   prev_loc     = twai_page$location3,
                                   geometry     = points.sfc ) ) 

points_map = sf::st_join( points_sf, tw_county, join=st_within )


non_poultry_i   = setdiff( grep( "野", points_map$species ), grep( "鹿野雞", points_map$species ) ) 
points_map_farm = points_map[ -non_poultry_i, ]  # n=1513 


f_fd_trees1 = "data/predict/cTW.contin_extract_trees/"
f_fd_trees2 = "data/predict/cTW.contin_extract_trees-2/"


target_areas = c( "Yunlin", "Changhua", "Chiayi", "Tainan", "Pingtung", "Kaohsiung" )

for( y in target_areas )
{
  X_county = tw_county[ grep( y, tw_county$COUNTYENG ), ] 
  
  if( dim(X_county)[1] > 1 )
  {
    X_county$geometry[1] = st_union(X_county)
    X_county = X_county[1,]
  }
  
  X_county.t = tw_county.t[ grep( y, tw_county.t$COUNTYENG ), ]  
  
  if( dim(X_county.t)[1] > 1 )
  {
    X_county.t$geometry[1] = st_union(X_county.t)
    X_county.t = X_county.t[1,]
  }
  

  
  points_X_farm  = points_map_farm[ grep( y, points_map_farm$COUNTYENG ), ]
  
  new_oubreaks_i = which( as_date(points_X_farm$sample_time) > as_date( latest_sample_s ) )
  points_X_farm$new_outbreak                   = FALSE
  points_X_farm$new_outbreak[ new_oubreaks_i ] = TRUE
  
  points_map_farm_t2 = points_map_farm[ which( decimal_date(points_map_farm$sample_time) >= phase_cutoff ),  ]
  points_map_farm_t2 = points_map_farm_t2[ which( decimal_date(points_map_farm_t2$sample_time) <= latest_sample ),  ]
  
  
  points_X_farm_t2  = points_map_farm_t2[ grep( y, points_map_farm_t2$COUNTYENG ), ]
  
  points_X_farm.t    = sf::st_transform( points_X_farm, crs = "EPSG:7801" )
  points_X_farm_t2.t = sf::st_transform( points_X_farm_t2, crs = "EPSG:7801" )
  
  
  # buffer area for outbreak sites 
  buffer0.5k_X_farm_t2.t = sf::st_buffer( points_X_farm_t2.t, dist = 500 )
  buffer1k_X_farm_t2.t   = sf::st_buffer( points_X_farm_t2.t, dist = 1000 )
  buffer2k_X_farm_t2.t   = sf::st_buffer( points_X_farm_t2.t, dist = 2000 )
  
  UniBuffer0.5k_X_farm_t2.t = sf::st_union( buffer0.5k_X_farm_t2.t )
  UniBuffer1k_X_farm_t2.t   = sf::st_union( buffer1k_X_farm_t2.t )
  UniBuffer2k_X_farm_t2.t   = sf::st_union( buffer2k_X_farm_t2.t )
  
  
  
  # 2 POLYGONS 
  
  # source other/random.seq.R
  

  
  # list_pg_rand1_0.5k = trees_uniPGs( folder = f_fd_rand1, treesN = 1000, areaD = 500, X_county=X_county )
  # list_pg_rand1_1k   = trees_uniPGs( folder = f_fd_rand1, treesN = 1000, areaD = 1000, X_county=X_county )
  # list_pg_rand1_2k   = trees_uniPGs( folder = f_fd_rand1, treesN = 1000, areaD = 2000, X_county=X_county )
  
  
  # circle areas by tree nodes 
  list_pg_trees1_0.5k = trees_uniPGs( folder = f_fd_trees1, treesN = 1000, areaD = 500, X_county=X_county )
  list_pg_trees1_1k   = trees_uniPGs( folder = f_fd_trees1, treesN = 1000, areaD = 1000, X_county=X_county )
  list_pg_trees1_2k   = trees_uniPGs( folder = f_fd_trees1, treesN = 1000, areaD = 2000, X_county=X_county )
  
  
  # tree2
  list_pg_trees2_0.5k = trees_uniPGs( folder = f_fd_trees2, treesN = 1000, areaD = 500, X_county=X_county )
  list_pg_trees2_1k   = trees_uniPGs( folder = f_fd_trees2, treesN = 1000, areaD = 1000, X_county=X_county )
  list_pg_trees2_2k   = trees_uniPGs( folder = f_fd_trees2, treesN = 1000, areaD = 2000, X_county=X_county )
  
  
  
  X_t2_num         = get_node_num( folder = f_fd_trees1, treesN = 1000, X_county=X_county )
  # list_pg_site_0.5k = site_uniPGs( X_t2_num, points_X_farm_t2.t, Nseed = 666, areaD = 500 )
  # list_pg_site_1k   = site_uniPGs( X_t2_num, points_X_farm_t2.t, Nseed = 666, areaD = 1000 )
  # list_pg_site_2k   = site_uniPGs( X_t2_num, points_X_farm_t2.t, Nseed = 666, areaD = 2000 )
  
  
  # circle areas by random sampling (n=node)
  list_null1_0.5k = null_uniPGs( vec_node = X_t2_num, areaD = 500, X_county.t = X_county.t )
  list_null1_1k   = null_uniPGs( vec_node = X_t2_num, areaD = 1000, X_county.t = X_county.t )
  list_null1_2k   = null_uniPGs( vec_node = X_t2_num, areaD = 2000, X_county.t = X_county.t )
  
  
  # circle areas by tree nodes + outbreak farms 
  list_pg_merge1_0.5k = lapply( as.list( seq_along(list_pg_trees1_0.5k) ), function(x) sf::st_union( list_pg_trees1_0.5k[[x]], UniBuffer0.5k_X_farm_t2.t ) )
  list_pg_merge1_1k   = lapply( as.list( seq_along(list_pg_trees1_1k) ), function(x) sf::st_union( list_pg_trees1_1k[[x]], UniBuffer1k_X_farm_t2.t ) )
  list_pg_merge1_2k   = lapply( as.list( seq_along(list_pg_trees1_2k) ), function(x) sf::st_union( list_pg_trees1_2k[[x]], UniBuffer2k_X_farm_t2.t ) )
  
  
  # tree2
  list_pg_merge2_0.5k = lapply( as.list( seq_along(list_pg_trees2_0.5k) ), function(x) sf::st_union( list_pg_trees2_0.5k[[x]], UniBuffer0.5k_X_farm_t2.t ) )
  list_pg_merge2_1k   = lapply( as.list( seq_along(list_pg_trees2_1k) ), function(x) sf::st_union( list_pg_trees2_1k[[x]], UniBuffer1k_X_farm_t2.t ) )
  list_pg_merge2_2k   = lapply( as.list( seq_along(list_pg_trees2_2k) ), function(x) sf::st_union( list_pg_trees2_2k[[x]], UniBuffer2k_X_farm_t2.t ) )
  
  
  
  # circle areas by random sampling (n=node + outbreak farms)
  Npoints_X_farm_t2 = dim(points_X_farm_t2.t)[1]
  list_null2_0.5k = null_uniPGs( vec_node = X_t2_num+Npoints_X_farm_t2, areaD = 500, X_county.t = X_county.t )
  list_null2_1k   = null_uniPGs( vec_node = X_t2_num+Npoints_X_farm_t2, areaD = 1000, X_county.t = X_county.t )
  list_null2_2k   = null_uniPGs( vec_node = X_t2_num+Npoints_X_farm_t2, areaD = 2000, X_county.t = X_county.t )
  
  
  
  
  # 3 EVALUATION
  
  # eval_rand1_0.5k = sapply( list_pg_rand1_0.5k, function(x) suppressMessages( test_pos( x, points_X_farm.t ) ) )
  # eval_rand1_1k   = sapply( list_pg_rand1_1k, function(x) suppressMessages( test_pos( x, points_X_farm.t ) ) )
  # eval_rand1_2k   = sapply( list_pg_rand1_2k, function(x) suppressMessages( test_pos( x, points_X_farm.t ) ) )
  
  eval_trees1_0.5k = sapply( list_pg_trees1_0.5k, function(x) suppressMessages( test_pos2( x, points_X_farm.t ) ) )
  eval_trees1_1k   = sapply( list_pg_trees1_1k, function(x) suppressMessages( test_pos2( x, points_X_farm.t ) ) )
  eval_trees1_2k   = sapply( list_pg_trees1_2k, function(x) suppressMessages( test_pos2( x, points_X_farm.t ) ) )
  
  eval_null1_0.5k = sapply( list_null1_0.5k, function(x) suppressMessages( test_pos2( x, points_X_farm.t ) ) )
  eval_null1_1k   = sapply( list_null1_1k, function(x) suppressMessages( test_pos2( x, points_X_farm.t ) ) )
  eval_null1_2k   = sapply( list_null1_2k, function(x) suppressMessages( test_pos2( x, points_X_farm.t ) ) )
  
  # eval_site_0.5k = sapply( list_pg_site_0.5k, function(x) suppressMessages( test_pos( x, points_X_farm.t ) ) )
  # eval_site_1k   = sapply( list_pg_site_1k, function(x) suppressMessages( test_pos( x, points_X_farm.t ) ) )
  # eval_site_2k   = sapply( list_pg_site_2k, function(x) suppressMessages( test_pos( x, points_X_farm.t ) ) )
  
  eval_merge1_0.5k = sapply( list_pg_merge1_0.5k, function(x) suppressMessages( test_pos2( x, points_X_farm.t ) ) )
  eval_merge1_1k   = sapply( list_pg_merge1_1k, function(x) suppressMessages( test_pos2( x, points_X_farm.t ) ) )
  eval_merge1_2k   = sapply( list_pg_merge1_2k, function(x) suppressMessages( test_pos2( x, points_X_farm.t ) ) )
  
  eval_null2_0.5k = sapply( list_null2_0.5k, function(x) suppressMessages( test_pos2( x, points_X_farm.t ) ) )
  eval_null2_1k   = sapply( list_null2_1k, function(x) suppressMessages( test_pos2( x, points_X_farm.t ) ) )
  eval_null2_2k   = sapply( list_null2_2k, function(x) suppressMessages( test_pos2( x, points_X_farm.t ) ) )
  
  
  
  # tree2
  
  eval_trees2_0.5k = sapply( list_pg_trees2_0.5k, function(x) suppressMessages( test_pos2( x, points_X_farm.t ) ) )
  eval_trees2_1k   = sapply( list_pg_trees2_1k, function(x) suppressMessages( test_pos2( x, points_X_farm.t ) ) )
  eval_trees2_2k   = sapply( list_pg_trees2_2k, function(x) suppressMessages( test_pos2( x, points_X_farm.t ) ) )
  
  eval_merge2_0.5k = sapply( list_pg_merge2_0.5k, function(x) suppressMessages( test_pos2( x, points_X_farm.t ) ) )
  eval_merge2_1k   = sapply( list_pg_merge2_1k, function(x) suppressMessages( test_pos2( x, points_X_farm.t ) ) )
  eval_merge2_2k   = sapply( list_pg_merge2_2k, function(x) suppressMessages( test_pos2( x, points_X_farm.t ) ) )
  
  
  ## out 
  
  ls_all  = ls()
  X_evals = grep( "^eval", ls_all )

  for( j in 1: length(X_evals) )
  {
    assign( paste0( y, "__", ls_all[ X_evals[j] ] ), get( ls_all[ X_evals[j] ] ) )
  }
  
  
  
}






output_names = paste0( "polygonX", ".", gsub( "-", "", date(Sys.time()) ), ".rda" )

save.image( file= output_names )








