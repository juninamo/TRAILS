createEmptyDf = function( nrow, ncol, colnames = c() ){
  if( missing( ncol ) && length( colnames ) > 0 ){
    ncol = length( colnames )
  }
  data.frame( matrix( vector(), nrow, ncol, dimnames = list( c(), colnames ) ) )
}