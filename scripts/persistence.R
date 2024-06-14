### parse_persistence ----

parse_persistence = function( treeFolder, cutoff_time, nCore=2 )
{
  # the default state name is `geo`
  
  require( readr )
  require( sets )
  require( purrr )
  require( doMC )
  
  registerDoMC( cores = nCore )
  
  dir_tsv = list.files( path = treeFolder, pattern = "tree", full.names = TRUE )
  n_tsv   = length( dir_tsv )
  
  # states = c()
  # n_states = 0
  
  pooled_df = 
    foreach( i = 1:n_tsv, .combine = "rbind" ) %dopar% 
    {
      tree_tab = suppressMessages( read_tsv( dir_tsv[i] ) )
      
      min_time = min( tree_tab$time_s )  
      max_time = max( tree_tab$time )
      
      tree_tab$geo_s = tree_tab$geo[ tree_tab$parent ]
      tree_tab_local = tree_tab[ which( tree_tab$geo == tree_tab$geo_s ), ]
      
      tree_df = tree_tab_local
      
      states = sort( unique(tree_df$geo) )
      c_state = c()
      c_t1_p  = c()
      c_t2_p  = c()
      
      t1_span = cutoff_time - min_time 
      t2_span = max_time - cutoff_time
      
      for( j in 1: length(states) )
      {
        state_i = states[j]
        
        df_i = subset( tree_df, geo==state_i )
        
        int0 = apply( cbind( df_i$time_s, df_i$time ),
                      1, 
                      function(x) sets::interval( x[1], x[2], "[]") )
        
        int1 = interval_union( lapply( int0, function(x){ sets::interval_intersection( x, sets::interval(-Inf, cutoff_time) ) } ) )
        int2 = interval_union( lapply( int0, function(x){ sets::interval_intersection( x, sets::interval(cutoff_time, Inf) ) } ) )
        
        t1_p = sum( unlist( sapply( int1, function(x) base::range(x)[2]-base::range(x)[1] ) ) ) / t1_span
        t2_p = sum( unlist( sapply( int2, function(x) base::range(x)[2]-base::range(x)[1] ) ) ) / t2_span
        
        c_state = c( c_state, state_i )
        c_t1_p  = c( c_t1_p, t1_p )
        c_t2_p  = c( c_t2_p, t2_p )
      }
      
      out_df = data.frame( tree = i, state = c_state, t1_p = c_t1_p, t2_p = c_t2_p )
      
      return(out_df)
    }
  
  return( pooled_df )
  cat( "done\n" )
  
  #v20231006
}