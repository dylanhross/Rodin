#' lipid.list.from.msdial()
#'
#' Clean and prepare a character vector prior to input it in the lipid.miner function.
#'
#' @param xlsx path to MS-DIAL output file ("Metabolite name" column used as lipid list) 
#'
#' @details 
#' This function is meant to be used prior to clean.lipid.list(), it loads a lipid
#' list from the specified XLSX file (MS-DIAL output file, uses the Metabolites 
#' column), makes various substitutions to the lipid names to make them compatible
#' with the parser (in 02_lipid_miner.R -> lipid.miner()), and returns the lipid
#' list as a character vector
#'
#' @return 1 column df with lipid names
#'
#' @author Dylan H. Ross
#' @export
#'

lipid.list.from.msdial<-function(xlsx){
  
  # load lipids from "Metabolite name" column from MS-DIAL results
  # as character vector
  msd.lipid.list<-pull(read_xlsx(xlsx, .name_repair = "unique_quiet"), "Metabolite name")
  
  # convert lines with multiple values (delimited by |) to only the rightmost value
  msd.lipid.list.cleaned<-unlist(lapply(strsplit(msd.lipid.list, split = "|", fixed = TRUE), tail, n = 1))
  
  # if there are still some lines with "<some_lipid>;<some_lipid>"
  # split those into separate lines (i.e. decombine)
  # need to use this complicated setup to distinguish between cases where
  # the line contains "<some_lipid>;<some_lipid>" or a single lipid that
  # has the oxygen annotations for a FA e.g. "Cer 18:1;O2/15:0"
  msd.lipid.list.decombined<-c()
  n.pre.decombine<-length(msd.lipid.list.cleaned)
  for (i in 1:n.pre.decombine) {
    l<-msd.lipid.list.cleaned[i]
    # ! does not match O after ; !
    mtch<-regexpr("(?<delim>;[ ]?)[ABCDEFGHIJKLMNPQRSTUVWXYZ]", l, perl = TRUE)
    cap.len<-attributes(mtch)$capture.length
    if (cap.len > 0) {
      split.start.pos<-attributes(mtch)$capture.start
      split.end.pos<-split.start.pos + cap.len
      l1<-substr(l, 1, split.start.pos - 1)
      l2<-substr(l, split.end.pos, nchar(l))
      # add each of the two split entries to the list
      msd.lipid.list.decombined<-append(msd.lipid.list.decombined, l1)
      msd.lipid.list.decombined<-append(msd.lipid.list.decombined, l2)
    } else {
      msd.lipid.list.decombined<-append(msd.lipid.list.decombined, l)
    }
  }
  n.post.decombine<-length(msd.lipid.list.decombined)
  # emit a warning if the length of the list changed due to decombining
  if (n.pre.decombine != n.post.decombine) {
    warning("lipid list length has changed due to decombining")
  }
  
  # remove ;O2 from sphingos?
  
  # convert "X C:U" to "X(C:U)"
  msd.lipid.list.paren<-gsub(" ", "(", msd.lipid.list.decombined)
  msd.lipid.list.paren[grepl("\\(", msd.lipid.list.paren)]<-unlist(lapply(msd.lipid.list.paren[grepl("\\(", msd.lipid.list.paren)], paste0, ")"))
  
  # make some statically defined substitutions
  msd.lipid.list.subbed<-gsub("CAR", "carnitine", msd.lipid.list.paren)  # change names of carnitines
  msd.lipid.list.subbed<-gsub("_", "/", msd.lipid.list.subbed)  # change _ chain separators to /
  msd.lipid.list.subbed<-gsub("(?<front>^LP[A-Z][(][0-9]+:[0-9]+)[)]", "\\1/0:0)", msd.lipid.list.subbed, perl = TRUE)  # ensure LPX lipids have explicit 0:0 chain
  msd.lipid.list.subbed<-gsub("^LP", "P", msd.lipid.list.subbed)  # change LPX lipids to PX
  msd.lipid.list.subbed<-gsub(";O[0-9]?", "", msd.lipid.list.subbed)  # remove oxygen annotations like ;O2 from sphingos
  msd.lipid.list.subbed<-gsub("(?<cls>(Cer|SM)[(])", "\\1d", msd.lipid.list.subbed, perl = TRUE)  # add d to sphingos
  msd.lipid.list.subbed<-gsub("FA[(][0-9]+:[0-9]+[)]", "", msd.lipid.list.subbed)  # get rid of weird FA(C:U) inside of some sphingos
  # others?
  
  # return as a 1-col data frame
  data.frame(msd.lipid.list.subbed)
}


#' clean.lipid.list()
#'
#' Clean and prepare a character vector prior to input it in the lipid.miner function.
#'
#' @param X is a character vector containing the lipids names.
#'
#' @details We recommend to use the clean.lipid.list() prior to use the lipid.miner()
#' @details to maximize its compatibility.
#' @details This function removes missing values, usless space characters, and redundancy.
#' @details clean.lipid.list() informs you using warnings on the steps it has achieved.
#'
#' @return The R object outputed by this function is a character vector
#'
#' @examples clean.lipid.list(QueryExample)
#' @examples clean.lipid.list(UniverseExample)
#'
#' @author Geremy Clair
#' @export
#'

clean.lipid.list<-function(X){
  
  # if the data had more than one column -> remove the data in the other columns
  cleaned<- X
  if(ncol(cleaned)>0){
    cleaned<-as.character(unlist(cleaned[,1]))
    warning("Your lipid list format was inappropriate (it had more than one column)")
    warning("only the first column was kept as lipid identifier for the analysis")
  }

  #convert in character vector (if not already)
  if (typeof(cleaned)!="character"){
    cleaned<- as.character(unlist(cleaned[,1]))
  }

  #remove m
  if(sum(cleaned=="")>0){
    cleaned[cleaned==""]<-NA
    cleaned<- na.omit(cleaned)
    warning("Your lipid list was containing empty ids")
    warning("These ids were removed")
  }

  #Remove duplicates (original duplicates)
  if(sum(duplicated(cleaned))>0){
    n.duplicates<-sum(duplicated(cleaned))
    cleaned<-unique(cleaned)
    warning(paste("Your list was containing", n.duplicates, "duplicated lipid identifiers"))
    warning("the duplicated values were removed")
  }

  #Convert return delimited list into vector
  if(length(cleaned)==1&&grepl("\n",cleaned)){
    cleaned<- unlist(strsplit(cleaned,"\n"))
    warning("Your list was return delimited and was converted in an appropriate format")
    if(length(cleaned)>1&&grepl(";",cleaned)){
      warning("Your list contain some elements with multiple identifiers separated by a semicolon or contains both semicolons and returns. Please ensure to use only one identifier for each lipid")
    }
    if(length(cleaned)>1&&grepl(",",cleaned)){
      warning("Your list contain some elements with multiple identifiers separated by a comma or contains both comma and returns. Please ensure to use only one identifier for each lipid")
    }
  }

  #Convert return delimited list into vector
  if(length(cleaned)==1&&grepl("\t",cleaned)){
    cleaned<- unlist(strsplit(cleaned,"\t"))
    warning("Your list was tab delimited and was converted in an appropriate format")
    if(length(cleaned)>1&&grepl(";",cleaned)){
      warning("Your list contain some elements with multiple identifiers separated by a semicolon or contains both semicolons and tab. Please ensure to use only one identifier for each lipid")
    }
    if(length(cleaned)>1&&grepl(",",cleaned)){
      warning("Your list contain some elements with multiple identifiers separated by a comma or contains both comma and tab. Please ensure to use only one identifier for each lipid")
    }
  }

  #Convert comma delimited list into vector
  if(length(cleaned)==1&&grepl(",",cleaned)){
    cleaned<- unlist(strsplit(cleaned,","))
    warning("Your list was comma delimited and was converted in an appropriate format")
    if(length(cleaned)>1&&grepl(";",cleaned)){
      warning("Your list contain some elements with multiple identifiers separated by a semicolon or contains both semicolons and commas. Please ensure to use only one identifier for each lipid")
    }
  }

  #Convert semicolon delimited list into vector
  if(length(cleaned)==1&&grepl(";",cleaned)){
    cleaned<- unlist(strsplit(cleaned,";"))
    warning("Your list was semicolon delimited and was converted in an appropriate format")
  }

  #remove _CoE, _NEG, _POS, and spaces
  rm.list<- c("_CoE", "_Coe", "_coe", "_NEG", "_Neg", "_neg", "_POS", "_Pos", "_pos", " ")
  for (i in 1:length(rm.list)){
    cleaned<-gsub(rm.list[i],"",cleaned)
  }

  #Remove duplicates (final duplicates)
  if(sum(duplicated(cleaned))>0){
    n.duplicates<-sum(duplicated(cleaned))
    cleaned<-unique(cleaned)
    warning(paste("Your list was containing", n.duplicates, "duplicates"))
    warning("the duplicated values were removed")
  }

  #count identifiers
  warning(paste("After reformating, your", "lipid" , "list comprise", length(cleaned), "lipid(s)"))
  cleaned
}

#' clean.RankingTable()
#'
#' Clean and prepare a character vector prior to input it in the lipid.miner function.
#'
#' @param X is a ranking table
#'
#' @details We recommend to use the clean.RankingTable() prior to use the lipid.miner() on the first column of the rankingTable (containing the lipid list)
#' @details to maximize its compatibility.
#' @details This function removes missing values, usless space characters, and redundancy.
#' @details clean.lipid.list() informs you using warnings on the steps it has achieved.
#'
#' @return The R object outputed by this function is a data.frame with two columns the first one contains the lipids, the second one the ranked weights
#'
#' @examples cleaned.rankingTable<-clean.rankingTable(RankingTableExample)
#'
#' @author Geremy Clair
#' @export
#'

clean.rankingTable<-function(df){
  # if the data had more than one column -> remove the data in the other columns
    if(ncol(df)>2){
    df<-df[,1:2]
    warning("Your lipid list format was inappropriate (it had more than two columns)")
    warning("only the first column was kept as lipid identifier and the second columns as ranking values")
  }

  #ensure that the first column of df is of character type
  df[,1]<-as.character(t(df[,1]))

  #ensure that the secon column of df is of numeric type
  df[,2]<-as.numeric(t(df[,2]))

  #remove the rows containing NAs or missing values
  df[df==""]<-NA
  df[df=="NA"]<-NA
  if(sum(is.na(df))>0){
    warning("Your lipid list was containing empty ids or ranking values, these were removed from the list")
    df<-df[!rowSums(is.na(df))>0,]
  }

  #Remove ID duplicates
  if(sum(duplicated(df[,1]))>0){
    n.duplicates<-sum(duplicated(df))
    df<-unique(df)
    warning(paste("Your list was containing", n.duplicates, "only the first of the duplicated ID(s) were kept"))
  }

  #remove _CoE, _NEG, _POS, and spaces from the lipid names
  rm.list<- c("_CoE", "_Coe", "_coe", "_NEG", "_Neg", "_neg", "_POS", "_Pos", "_pos", " ")
  for (i in 1:length(rm.list)){
    cleaned<-gsub(rm.list[i],"",df)
  }

  #count identifiers
  warning(paste("After reformating, your", "lipid" , "list comprise", nrow(df), "lipid(s)"))

  #output the table
  return(df)

}

