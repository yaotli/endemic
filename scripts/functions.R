### fastEx ----

fastaEx <- function( filedir = file.choose() )
{
  suppressWarnings(suppressMessages( require(seqinr) ) )
  
  file     <- seqinr::read.fasta( filedir, forceDNAtolower = FALSE )
  file_seq <- seqinr::getSequence( file )
  file_id  <- attributes( file )$names
  
  return( list( seq = file_seq, 
                id  = file_id ) )
  #v20220708
}



### tagExtra ----

tagExtra <- function( filedir = file.choose() )
{
  require(stringr)
  
  anno.tre <- read.csv( filedir, stringsAsFactors = FALSE )
  taxa.s   <- grep( x = anno.tre[,1], pattern = "taxlabels" ) + 1
  ntax     <- as.numeric( str_match( grep( x       = anno.tre[,1],
                                           pattern = "ntax", 
                                           value   = TRUE) , 
                                     "(ntax=)([0-9]+)")[,3] )
  
  taxa.e <- taxa.s + ntax - 1
  
  tre_name <- gsub( "\\[.*\\]$", "", gsub( "\t|'", "", anno.tre[, 1][taxa.s: taxa.e] ) )
  
  tre_tag  <- str_match( string  = gsub( "\t|'", "", anno.tre[, 1][taxa.s: taxa.e] ),
                         pattern = "&!color=#([A-Z0-9a-z]+)" )[,2]
  
  if( TRUE %in% is.na(tre_name) ){ warning( paste0( "no. NA in name= ", 
                                                    length( which( is.na(tre_name) ) ) 
                                                   ) 
                                           ) }
  
  return( df = data.frame( id  = tre_name,
                           tag = tre_tag, stringsAsFactors = FALSE) )
  
  #v20220908
}



### jumpMx ----

jumpMx <- function( states = c("cnC", "cnE", "cnN", "cnS", "cnSW"), prefix = "" )
{
  states = sort( unique( states ) )
  lth <- length( states )
  mx  <- matrix( rep(0,lth^2), nrow = lth)
  ord <- combn( lth, 2)
  
  # two-by-two transition
  dirname <- c() 
  string  <- c()
  for( x in 1: ( dim(ord)[2] ) )
  {
    dirname = c( dirname, paste0( states[ ord[,x][1] ], "_to_", states[ ord[,x][2] ] ) )
    dirname = c( dirname, paste0( states[ rev(ord[,x])[1] ], "_to_", states[ rev(ord[,x])[2] ] ) )
  }
  
  for( x in 1: ( dim(ord)[2] ) )
  {
    mx1  <- matrix( rep(0,lth^2), nrow = lth)
    mx1[ ord[,x][1], ord[,x][2]  ] <- 1
    string  <- c( string, paste( paste( as.vector( t(mx1) ), ".0", sep = ""), collapse = " ") )
    
    mx2  <- matrix( rep(0,lth^2), nrow = lth)
    mx2[ ord[,x][2], ord[,x][1] ] <- 1
    string  <- c( string, paste( paste( as.vector( t(mx2) ), ".0", sep = ""), collapse = " ") )
  }
  
  out <- c()
  for( y in 1: length(dirname) )
  {
    out <- c(out, paste("<parameter id=", dirname[y], " value=", string[y], "/>",  sep = '\"' ) )
  }
  
  # all_in and all_out
  
  dirname2 <- as.vector( sapply( as.list( states ), 
                                 function(x){ y = c( paste0( "to_", x ), paste0( "from_", x ) )  } )  )
  
  string2 <- c()
  for( j in 1: lth )
  {
    in.mx = matrix( rep(0,lth^2), nrow = lth)
    in.mx[ ,j ]  = 1
    in.mx[ j,j ] = 0
    string2  <- c( string2, paste( paste( as.vector( t(in.mx) ), ".0", sep = ""), collapse = " ") )
    
    out.mx = matrix( rep(0,lth^2), nrow = lth)
    out.mx[ j, ]  = 1
    out.mx[ j,j ] = 0
    string2  <- c( string2, paste( paste( as.vector( t(out.mx) ), ".0", sep = ""), collapse = " ") )
  }
  out2 <- sapply( as.list( seq(1, length(dirname2)) ), 
                  function(x) paste("<parameter id=", dirname2[x], " value=", string2[x], "/>",  sep = '\"' )  )
  
  
  # reward 
  blls = c()
  for( r in 1: length(states) )
  {
    bll  <- rep(0, length(states) )
    bll[ r ] = 1
    blls = c( blls, paste0( bll, collapse = " ") )
  }
  
  rw <- paste0( " <parameter value=\"", blls, "\" id=\"reward_", states, "\"/>")
  rw <- c( "<rewards>", rw, "</rewards>")
  
  write.table( rw,   file = paste0( prefix, "_reward.mx.txt" ), sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
  write.table( out,  file = paste0( prefix, "_trans.mx.txt" ), sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
  write.table( out2, file = paste0( prefix, "_trans2.mx.txt" ), sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
  
  #v20221206
}





### para_mx ----

para_mx = function( df, std = TRUE, one_way_input = TRUE, para_name = "", col_data = 2 )
{
  .log_std = function( x )
  {
    log_x = log(x)
    
    mean_log_x = mean(log_x)
    sd_log_x   = sd(log_x)
    
    out = (log_x - mean_log_x)/sd_log_x
    
    return( out )
  }
  
  if( !(is.data.frame( df ) | is.matrix( df ) ) ){ stop( "only datadrame/matrix as input" ) }
  
  
  if( one_way_input )
  {
    df_data = df[, col_data]
    
    dim_file = dim( df )
    
    n = length( df_data )
    
    null_mx_st1 = c()
    null_mx_st2 = c()
    
    for( x in 1:(n-1) )
    {
      for( y in (x+1): n )
      {
        null_mx_st1 = c( null_mx_st1, paste0( x, "_", y ) ) 
      }
    }
    
    for( y in 1:(n-1) )
    {
      for( x in (y+1):n )
      {
        null_mx_st2 = c( null_mx_st2, paste0( x, "_", y ) ) 
      }
    }
    
    null_mx = c( null_mx_st1, null_mx_st2 )
    
    origin_i = as.numeric( gsub( "_[0-9]+$", "", null_mx ) )
    dest_i   = as.numeric( gsub( "^[0-9]+_", "", null_mx ) )
    
    origin_raw = df_data[ origin_i ]
    dest_raw   = df_data[ dest_i ]
    
    
    if( std )
    {
      origin_out = paste0( .log_std( origin_raw ), collapse = " " )
      dest_out   = paste0( .log_std( dest_raw ), collapse = " " )
      
    }else
    {
      origin_out = paste0( origin_raw, collapse = " " )
      dest_out   = paste0( dest_raw, collapse = " " )
    }
    
    
    origin_st = paste("<parameter id=", paste0(para_name, "_origin"), " value=", origin_out, "/>",  sep = '\"' )
    dest_st   = paste("<parameter id=", paste0(para_name, "_destination"), " value=", dest_out, "/>",  sep = '\"' )
    
    write.table( origin_st,  
                 file = paste0(para_name, "_origin.txt"), 
                 sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE )
    write.table( dest_st,  
                 file = paste0(para_name, "_dest.txt"), 
                 sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE )
    
  }else
  {
    dim_file = dim( df )
    
    n = min( dim_file )
    
    max_file = df
    
    out_st = c()
    
    for( x in 1:(n-1) )
    {
      for( y in (x+1): n )
      {
        out_st = c( out_st, max_file[[x,y]] ) 
      }
    }
    
    for( y in 1:(n-1) )
    {
      for( x in (y+1):n )
      {
        out_st = c( out_st, max_file[[x,y]] ) 
      }
    }
    
    if( std )
    {
      out_value = paste0( .log_std( out_st ), collapse = " " )
      
    }else
    {
      out_value = paste0( out_st, collapse = " " )
    }
    
    out = paste("<parameter id=", paste0(para_name, "_mx"), " value=", out_value, "/>",  sep = '\"' )
    
    write.table( out,  
                 file = paste0(para_name, "_mx.txt"), 
                 sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE )
  }
  
  
  #v20230605
}





### jumpMx2 ----

jumpMx2 <- function( states = c("cnC", "cnE", "cnN", "cnS", "cnSW"), prefix = "" )
{
  states = sort( unique( states ) )
  
  lth    = length( states )
  mx_lth = lth*(lth-1)
  
  
  # module from function:para_mx
  null_mx_st1 = c()
  null_mx_st2 = c()
  for( x in 1:(lth-1) )
  {
    for( y in (x+1): lth )
    {
      null_mx_st1 = c( null_mx_st1, paste0( x, "_", y ) ) 
    }
  }
  
  for( y in 1:(lth-1) )
  {
    for( x in (y+1):lth )
    {
      null_mx_st2 = c( null_mx_st2, paste0( x, "_", y ) ) 
    }
  }
  
  null_mx = c( null_mx_st1, null_mx_st2 )
  
  dirname2 <- as.vector( sapply( as.list( states ), 
                                 function(x){ y = c( paste0( "from_", x ), paste0( "to_", x ) )  } )  )
  
  string2 <- c()
  for( i in 1: lth )
  {
    # from
    mx_i = grep( paste0( "^", i, "_" ), null_mx )
    st_i = rep(0, mx_lth)
    st_i[mx_i] = 1
    
    # to
    mx_j = grep( paste0( "_", i, "$"  ), null_mx )
    st_j = rep(0, mx_lth)
    st_j[mx_j] = 1
    
    string2 = c( string2,  paste( paste( st_i, ".0", sep = ""), collapse = " " ) )
    string2 = c( string2,  paste( paste( st_j, ".0", sep = ""), collapse = " " ) )
  }
  
  out2 <- sapply( as.list( seq(1, length(dirname2)) ), 
                  function(x) paste("<parameter id=", dirname2[x], " value=", string2[x], "/>",  sep = '\"' )  )
  
  
  write.table( out2, file = paste0( prefix, "_trans2.mx.txt" ), sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
  
  #v20230609
}




### pyCol ----

pyCol <- function( name = c( "red", "blue", "green" ) )
{
  Col.code <- c( "#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd", 
                 "#8c564b", "#e3777c2", "#7f7f7f", "#bcbd22", "#17becf" )
  Col.name <- c( "blue", "orange", "green", "red", "purple",
                 "brown", "pink", "gray", "yellow", "cyan")
  
  
  return( Col.code[ match( name, Col.name ) ] )
  
  #v20180110
  
}






### BF ----

BF = function( mean_indicator = NULL, n_location = 15, asym = TRUE )
{
  K   = n_location
  
  if( asym )
  {
    eta = K - 1 # prior 
    
    qk = (eta + 0)/( K*(K-1) ) # prior expection 
    
  }else
  {
    eta = log(2) # prior 
    
    qk = (eta + K - 1)/( K*(K-1)/2 ) # prior expection   
  }
  
  prior_odds = qk/(1-qk)
  post_odds  = mean_indicator/(1-mean_indicator)
  
  bf = post_odds / prior_odds
  
  return(bf)
  
  
  #v20230909
}



### read_beast_log ----

read_beast_log = function( log_file, n_burin = 10^7, n_logEvery = 10^4, sampleRow = 1 )
{
  library( tracerer )
  
  in_file = read.table( log_file, header = TRUE, comment.char = "#" )
  
  n_log = dim( in_file )[1]
  r_s   = n_burin / n_logEvery
  r_e   = n_log
  
  if( sampleRow != 1 )
  {
    rr = sort( sample( seq( r_s, r_e ), floor( sampleRow*n_log ) ) )
    
  }else
  {
    rr = seq( r_s+1, r_e )  
  }
  
  data_log = in_file[ rr, ]
  summ_log = apply( data_log, 2, function(x) calc_summary_stats_trace(x,1) )
  
  summ_df = do.call( rbind, summ_log )
  
  return( summ_df )
}





### sum_mj ----

sum_mj = function( history_txt, cutoff_time )
{
  # history data frame should include columns of
  # 1 treeId, 2 startLocation, 3 endLocation, 4 time
  
  require( tracerer )
  require( dplyr )
  suppressMessages( require( tidyverse ) )
  
  if( !is.data.frame( history_txt ) ){ stop( "input must be a dataframe" ) }
  
  tab_in = history_txt
  
  tab_in$state = paste( tab_in$startLocation, tab_in$endLocation, sep = "_" )
  
  #unique_dir = sort( unique( tab_in$state ) )
  uniq_state = sort( unique( c( tab_in$startLocation, tab_in$endLocation ) ) )
  
  epoch0_i = which( tab_in$time > cutoff_time )
  epoch1_i = which( tab_in$time <= cutoff_time )
  
  tab_in$state[ epoch0_i ] = paste0( tab_in$state[ epoch0_i ], "_0" )
  tab_in$state[ epoch1_i ] = paste0( tab_in$state[ epoch1_i ], "_1" )
  
  
  tab_sum = 
    tab_in %>% 
    dplyr::group_by( treeId ) %>%
    dplyr::count(state ) %>%
    tidyr::pivot_wider( names_from = state, values_from = n  ) 
  
  
  tab_sum[ is.na(tab_sum) ] = 0
  
  e0_c = grep( "_0$", names(tab_sum ) )
  e1_c = grep( "_1$", names(tab_sum ) )
  
  tab_sum$sum0 = apply( tab_sum, 1, function(x) sum(x[e0_c]) )
  tab_sum$sum1 = apply( tab_sum, 1, function(x) sum(x[e1_c]) )
  
  tab_sum$sum = apply( tab_sum, 1, function(x) x[length(x)] + (x[length(x)-1]) )
  
  print( "summarizing stat..." )
  sum_ls = apply( tab_sum, 2, function(x) calc_summary_stats_trace(x,1))
  sum_df = do.call( rbind, sum_ls )
  
  
  from_  = c()
  to_    = c()
  epoch_ = c()
  n_     = c()
  
  for( f in uniq_state )
  {
    for( g in  uniq_state )
    {
      for( e in c(0,1) )
      {
        dir_ = paste( f, g, e, sep = "_" )
        
        from_  = c( from_, f )
        to_    = c( to_, g )
        epoch_ = c( epoch_, e )
        
        rr = which( rownames( sum_df ) == dir_ )
        
        if( any( rr ) )
        {
          
          n_ = c( n_, sum_df$mean[ rr ] )
          
        }else
        {
          n_ = c( n_, 0 )  
        }
      }
    }
  }
  
  
  mx_df = data.frame( from = from_, to = to_, epoch = epoch_, n = n_ )
  
  mx_df_prop = mx_df
  mx_df_prop$n[ which( mx_df_prop$epoch==0 ) ] = mx_df_prop$n[ which( mx_df_prop$epoch==0 ) ]/ sum_df["sum0",1]
  mx_df_prop$n[ which( mx_df_prop$epoch==1 ) ] = mx_df_prop$n[ which( mx_df_prop$epoch==1 ) ]/ sum_df["sum1",1]
  
  
  out_ls = list( sum_df, mx_df, mx_df_prop )
  
  return( out_ls )
}




randIdx =  function( x, n.sample=1 )
{
  dup.i = which( duplicated( x ) )
  
  Ndup.i = seq_along( x )[ -dup.i ]
  
  dup.ls = as.list( Ndup.i )
  
  for( i in 1: length( dup.ls ) )
  {
    dup.ls[[i]] = which( x == x[ dup.ls[[i]] ] )
  }
  
  s.ls = sapply( dup.ls, function(x) x[ sample( seq_along(x),size = n.sample ) ] )
  
  out = sort( s.ls )
  
  return( out )
  
  #v20240304
}
  