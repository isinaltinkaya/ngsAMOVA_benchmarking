check_na_col <- function(df){
  for(i in colnames(df)){
    if(sum(is.na(df))!=0){
      warning(i, " contains NA")
    }else{
      cat(i,": OK\n")
    }
  }
}
