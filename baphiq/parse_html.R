require( xml2 )
require( stringr )
require( readr )
require( lubridate )

require(googleway)

yao_key=""

## 0 read-in the html ####

# Obsolete - 
# f_baphiq_2015Jan         = "output/twai_baphiq-2015Jan.sim.html"
# f_baphiq_2015Feb_2022Sep = "output/twai_baphiq-2015Feb-2022Sep.sim.html"
# baphiq_2015Jan         = read_lines( f_baphiq_2015Jan )
# baphiq_2015Feb_2022Sep = read_lines( f_baphiq_2015Feb_2022Sep )


f_baphiq_2015Jan         = "twai_baphiq-2015Jan.html"
f_baphiq_2015Feb_2022Dec = "twai_baphiq-2015Feb-2022Dec.html"

for( f in 1:2 )
{
  html_i = scan( file = c( f_baphiq_2015Jan, f_baphiq_2015Feb_2022Dec )[f], 
                 what = "", sep = "\n", blank.lines.skip = FALSE )
  
  for( i in 1: length(html_i) )
  {
    if( html_i[i] == '            <div id="drag_tbldata_list" class="drag-content">' )
    {
      listtime_l = i+2
      break
    }
  }
  
  assign( c( "baphiq_2015Jan", "baphiq_2015Feb_2022Dec" )[f], html_i[ listtime_l ] )
}



# function for text parsing 

parse_listtime = function( x )
{
  popup_unit = strsplit( x, "listtime" )[[1]]
  popup_unit = popup_unit[-1]
  
  no = as.numeric( str_extract( popup_unit, "^[0-9]{1,3}" ) )
  
  if( TRUE %in% duplicated(no) ){ stop( "duplicated no" ) }
  
  details = str_match( popup_unit, "getDetail\\((.*[0-9]+)\\)" )[,2]
  
  ls_details = strsplit( gsub( "'","", details ), "," )
  
  long_details = lapply( ls_details, 
                         function(x)
                         {
                           if( length(x) > 10 ){ x[9] = paste( x[9: length(x)], collapse = "_" ) }
                           
                           df_x = data.frame( t(x) )
                           names( df_x ) = c( "Lat", "Long", "Species", "Sampling_date", "Confirming_date", "Subtype", 
                                              "Patho_type", "Report_type", "Info9", "Info10" )
                           
                           return(df_x)
                         } )
  
  out_df = do.call( rbind, long_details )
  out_df = cbind( no, out_df )
  
  
  return(out_df)
}



# 1 combine the data ####

df_baphiq_2015Jan         = parse_listtime(baphiq_2015Jan)
df_baphiq_2015Feb_2022Dec = parse_listtime(baphiq_2015Feb_2022Dec)


df_baphiq_2015Jan$no         = paste0( df_baphiq_2015Jan$no, "a" )
df_baphiq_2015Feb_2022Dec$no = paste0( df_baphiq_2015Feb_2022Dec$no, "b" )


df_baphiq = rbind( df_baphiq_2015Jan, df_baphiq_2015Feb_2022Dec )

df_baphiq$Sampling_date   = as_date(df_baphiq$Sampling_date)
df_baphiq$Confirming_date = as_date(df_baphiq$Confirming_date)

df_baphiq$Subtype = str_replace_all( df_baphiq$Subtype, c( "H5N2、H5N3"="H5N2_H5N3",
                                                           "H5N2、H5N5"="H5N2_H5N5",
                                                           "H5N2、H5N6"="H5N2_H5N6",
                                                           "H5N2、H5N8"="H5N2_H5N8",
                                                           "H5N3、H5N8"="H5N3_H5N8",
                                                           "H5高病原"="H5Nx") )


write_tsv( df_baphiq, "baphiq_map.tsv",  na = "" )



# 2 assign the city/county by searching the Google map ####

# Inactivated by default -

tem_address = c()
for( n in 1: nrow(df_baphiq) )
{
  lat_long = paste( df_baphiq$Lat[n], df_baphiq$Long[n], sep = " " )

  g_output = google_places( search_string = lat_long, key = yao_key )

  tem_address = c( tem_address, g_output$results$formatted_address )

  print( tem_address[ length(tem_address) ] )
}

write.table( tem_address, "googlemap_results.20230222.txt", quote = FALSE )



st_city_county = "[A-Za-z ]+ City|[A-Za-z ]+ County"

city_county = str_extract_all( tem_address, st_city_county )
city_county = sapply( city_county, 
                      function(x)
                      {
                        if( length(x) > 1 )
                        {
                          x = x[ grep( "County", x ) ]
                          
                        }else if( length(x) == 0 )
                        {
                          x = NA
                        }
                        
                        return(x)
                      } )



na_j = which( is.na(city_county) ) # tem_address[na_j]

city_county[ na_j ] = c( "Yunlin County", "Yunlin County", "Pingtung County",
                         "Kaohsiung City" )


# 3 compile data with searching results and export  ####

df_baphiq$location3 = city_county
df_baphiq$location3 = gsub( "^ ", "", df_baphiq$location3 )

write_tsv( df_baphiq, "baphiq_website.20230222.tsv",  na = "" )



