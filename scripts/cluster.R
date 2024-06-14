library( treeio )
library( ggtree )
library( ggplot2 )
library( stringr )

source( "scripts/functions.R" )


# METHODS ####

getDes = function(node)
{
  phytools::getDescendants( tree =  merged_tree@phylo, node = node)
}





# 0 Input ####

in_time_tree = "data/cluster/cluster_timetree.nexus"

in_anno_tree  = "data/cluster/cluster_annotated_tree.nexus"
in_anno_tree2 = "data/cluster/cluster_annotated_tree2.nexus"



# 1 parsing annotated trees ####

time_tree = read.beast( in_time_tree )

# tree1 
anno_tree = read.beast( in_anno_tree )

merged_tree = merge_tree( time_tree, anno_tree )

treedf = fortify( merged_tree )
treedf = treedf[ ,-which(names(treedf) == "mutations") ]

treedf$location2 = treedf$location
treedf$location2[ which( treedf$location2 == "Hong_Kong" ) ] = "China"

treedf$location_p = as.character( treedf$location2[ treedf$parent ] )

treedf$date = as.numeric( treedf$date )



# tree2: adjusted freq rate
anno_tree2 = read.beast( in_anno_tree2 )

merged_tree2 = merge_tree( time_tree, anno_tree2 )

treedf2 = fortify( merged_tree2 )
treedf2 = treedf2[ ,-which(names(treedf2) == "mutations") ]

treedf2$location2 = treedf2$location
treedf2$location2[ which( treedf2$location2 == "Hong_Kong" ) ] = "China"

treedf2$location_p = as.character( treedf2$location2[ treedf2$parent ] )

treedf2$date = as.numeric( treedf2$date )




# 2 defining clusters ####

treedf$diffusion          = treedf$location2 != treedf$location_p
treedf$diffusion_sing     = treedf$diffusion & treedf$isTip
treedf$diffusion_internal = treedf$diffusion & (!treedf$isTip) 

treedf2$diffusion          = treedf2$location2 != treedf2$location_p
treedf2$diffusion_sing     = treedf2$diffusion & treedf2$isTip
treedf2$diffusion_internal = treedf2$diffusion & (!treedf2$isTip) 




## 2.2 method ####


findCluster = function( fort_tab )
{
  # fort_tab = fortified table
  out_list = list()
  
  fort_tab$cluster = NA
  
  ls_cluster = list()
  for( i in 1: length( which( fort_tab$diffusion_internal ) ) )
  {
    diffusion_internal_i = which( fort_tab$diffusion_internal )[i]
    
    diffusion_internal_i.dec = getDes( diffusion_internal_i )
    
    diffusion_internal_i.dec_diffusion = diffusion_internal_i.dec[ which( diffusion_internal_i.dec %in% which( fort_tab$diffusion ) ) ]
    
    
    if( length( diffusion_internal_i.dec_diffusion ) == 0 )
    {
      diffusion_internal_i.dec_local = c( diffusion_internal_i,  diffusion_internal_i.dec )
      
    }else
    {
      diffusion_internal_i.dec_diffusion.dec = unique( unlist( sapply( as.list( diffusion_internal_i.dec_diffusion ),
                                                                       function(x)
                                                                       {
                                                                         y = c(x, getDes( x ) )
                                                                         return(y)
                                                                       } ) ) )
      
      diffusion_internal_i.dec_local = diffusion_internal_i.dec[ -which( diffusion_internal_i.dec %in% 
                                                                           diffusion_internal_i.dec_diffusion.dec ) ]  
      
      diffusion_internal_i.dec_local = c( diffusion_internal_i, diffusion_internal_i.dec_local )
    }
    
    ls_cluster[[ length(ls_cluster)+1 ]] = diffusion_internal_i.dec_local
    
    fort_tab$cluster[ diffusion_internal_i.dec_local ] = paste0( "cluster_", i )
  }
  
  out_list[[1]] = fort_tab
  out_list[[2]] = ls_cluster
  
  return( out_list )
}

sumCluster = function( outlist )
{
  
  .tredata     = outlist[[1]]
  .clusterList = outlist[[2]]
    
  tip.n = length( which( .tredata$isTip ) )
  
  
  temp_out = 
  lapply( .clusterList, 
          function(x)
          {
            .no = .tredata$cluster[ x[1] ] 
              
            .loc.d = .tredata$location2[ x[1] ] 
            
            .loc.o = .tredata$location_p[ x ] 
            .loc.o = .loc.o[ which( !.loc.o %in% .loc.d ) ]
            
            .span  = range( .tredata$date[ x ])[2] - range( .tredata$date[ x ])[1] 
            .span2 = range( .tredata$date[ x[x<= tip.n ] ] )[2] - 
              range( .tredata$date[ x[x<=tip.n] ] )[1]
              
            .intro_date  = .tredata$date[ x[1] ]
            .intro_date2 = min( .tredata$date[ x[x<=tip.n ] ] )
              
            .rep_strain = .tredata$label[ x[x<=tip.n][1] ]
            
            out = data.frame( cluster     = .no,
                              loc.d       = .loc.d,
                              loc.o       = .loc.o,
                              span        = .span,
                              span2       = .span2,
                              intro_date  = .intro_date,
                              intro_date2 = .intro_date2,
                              rep_strain  = .rep_strain )
            return( out )           
          })
  
  summ_out = do.call( rbind, temp_out )
   
  return( summ_out )  
     
}



# run 
cluster_treedf = findCluster( treedf )

summ_df = sumCluster( cluster_treedf )


cluster_treedf2 = findCluster( treedf2 )

summ_df2 = sumCluster( cluster_treedf2 )



# examining with trees
merged_tree@extraInfo$cluster = NA
merged_tree@extraInfo$cluster = cluster_treedf[[1]]$cluster[ match( merged_tree@extraInfo$node,
                                                                    cluster_treedf[[1]]$node ) ]
merged_tree2@extraInfo$cluster = NA
merged_tree2@extraInfo$cluster = cluster_treedf2[[1]]$cluster[ match( merged_tree2@extraInfo$node,
                                                                    cluster_treedf2[[1]]$node ) ]

# write.beast( merged_tree, "cluster_tree1.tre" )
# write.beast( merged_tree2, "cluster_tree2.tre" )



# 3 summarizing clusters ####

# mean/median
# 95% quantile 
# to/from proportion
# >95% quantile

summary( summ_df$span2 )
summary( summ_df2$span2 )

quantile( summ_df$span2, 0.95 ); length( which( summ_df$span2>3 ) )
quantile( summ_df2$span2, 0.95 ); length( which( summ_df2$span2>3 ) )


cluster3yr = rbind( summ_df[ which( summ_df$span2 > 3 ), ],
                    summ_df2[ which( summ_df2$span2 > 3 ), ] )

cluster3yr$method = c( rep(1, length( which( summ_df$span2>3 ) ) ), 
                       rep(2, length( which( summ_df2$span2>3 ) ) ) )



# 4 Visualization ####

summ_df$area  = NA
summ_df2$area = NA 

## 4.1 Classification ####

gCN  = c( "China", "Hong_Kong" )
gSEA = c( "Lao", "Malaysia", "Myanmar", "Thailand", "Vietnam", "Indonesia", "Viet_Nam", "Lao_s_Democratic_Republic", 
          "Cambodia", "Afghanistan" )
gWA  = c( "Bangladesh", "Bhutan", "India", "Lebanon", "Nepal", "United_Arab_Emirates", "Egypt", "Pakistan", "Kuwait", "Kazakhstan",
          "Israel", "Iraq", "Iran", "Gaza_Strip", "Saudi_Arabia", 
          "Palestinian")
gNA  = c( "Japan", "Mongolia", "Russia", "South_Korea", "Korea",
          "Taiwan" )
gA   = c( "Burkina_Faso", "Cote_dIvoire", "Niger", "Nigeria", "Togo", "Ghana", "Zimbabwe", "Uganda", "South_Africa",
          "Democratic_Republic_of_the_Congo", "Cameroon", "Benin", "Sudan", 
          "Congo", "Senegal", "Lesotho" )
gE   = c( "Austria", "Bulgaria", "Belgium", "Denmark", "France", "Germany", "Hungary", "Italy", "Netherlands", "Poland", 
          "Sweden", "Switzerland", "United_Kingdom", "Croatia", "Czech_Republic", "Ukraine", "Turkey", "Spain", "Slovenia", "Slovakia",
          "Romania", "Macedonia", "Bosnia_and_Herzegovina", "Azerbaijan", "Greece", "Georgia",
          "Finland", "Estonia", "Norway", "Lithuania", "Albania", "Kosovo", "Ireland", "Luxembourg" )
gAm  = c( "USA", "Canada", "United_States" )

g.ls = list( gCN, gSEA, gWA, gNA, gA, gE, gAm )
g.ls = unlist( lapply( g.ls, function(x) paste0( x, collapse = "|" ) ) )

for( i in 1: length(g.ls) )
{
  summ_df$area[ grep( g.ls[i], summ_df$loc.d, ignore.case = TRUE ) ] = c( "gCN", "gSEA", "gWA", "gNA", "gA", "gE", "gAm" )[i]
  
  summ_df2$area[ grep( g.ls[i], summ_df2$loc.d, ignore.case = TRUE ) ] = c( "gCN", "gSEA", "gWA", "gNA", "gA", "gE", "gAm" )[i]
}




