#' Calculate metric values
#'
#' This function calculates metric values for bugs, fish, and algae.
#' Inputs are a data frame of taxa, SampleID, and community.
#' NonUnique = taxa to be excluded from taxa metrics (i.e., ambiguous taxa in the sample)
#' Benthic Macroinvertebrates
#' Required fields: SampleID (provided by user), Count, NonUnique, FinalID, Order, Family, Tribe, Genus, FFG, TolVal, Habit, Voltinism, BCG_Atr
#' FFG = CF, CG, SC, SH, PR
#' Habit = BU, CB, CN, SP, SW
#' Voltinism = multivoltine, semivoltine, univoltine
#' Fish
#' Required fields: SampleID (provided by user), Count, NonUnique, FinalID, Genus, Trophic, Lithophil, StWidAvg, StLength, MassTotal
#' Trophic =
#' Lithophil =
#' Stratum =
#' Tolerance =
#' Metric names have a prefix and suffix.  The prefix follows the naming convention listed below.  The suffix is a short name of the group being measured.  Only phylogenetic names are first letter caps.
#' ni = number of individuals
#' nt = number of taxa
#' pi = percent of individuals
#' pt = percent of taxa
#' x = index
#'
#
#' @param fun.DF Data frame of taxa (list required fields)
#' @param fun.SampID Sample ID column name
#' @param fun.Community Community name for which to calculate metric values (bugs, fish, or algae)
#' @param fun.MetricNames Optional vector of metric names to be returned.  If none are supplied then all will be returned.
#' @param boo.Adjust Optional boolean value on whether to perform adjustments of values prior to scoring.  Default = FALSE.
#' @return data frame of SampleID and metric values
#' @examples
#' # Index, Benthic Macroinvertebrates, genus
#' myIndex <- "MBSS.2005.Bugs"
#' # Thresholds
#' thresh <- metrics_scoring
#' # get metric names for myIndex
#' (myMetrics.Bugs <- as.character(droplevels(unique(thresh[thresh[,"Index.Name"]==myIndex,"Metric"]))))
#' # Taxa Data
#' myDF.Bugs <- taxa_bugs_genus
#' myMetric.Values.Bugs <- metric.values(myDF.Bugs, "SampleID", "bugs", myMetrics.Bugs)
#' View(myMetric.Values.Bugs)
#'
#' # Index, Benthic Macroinvertebrates, family
#' myIndex <- "MSW.1999.Bugs"
#' # Thresholds
#' thresh <- metrics_scoring
#' # get metric names for myIndex
#' (myMetrics.Bugs <- as.character(droplevels(unique(thresh[thresh[,"Index.Name"]==myIndex,"Metric"]))))
#' # Taxa Data
#' myDF.Bugs <- taxa_bugs_family
#' myMetric.Values.Bugs <- metric.values(myDF.Bugs, "SampleID", "bugs", myMetrics.Bugs)
#' View(myMetric.Values.Bugs)
#'
#' # Index, Fish
#' myIndex <- "MBSS.2005.Fish"
#' # Thresholds
#' thresh <- metrics_scoring
#' # get metric names for myIndex
#' (myMetrics.Fish <- as.character(droplevels(unique(thresh[thresh[,"Index.Name"]==myIndex,"Metric"]))))
#' # Taxa Data
#' myDF.Fish <- taxa_fish
#' myMetric.Values.Fish <- metric.values(myDF.Fish, "SampleID", "fish", myMetrics.Fish, TRUE)
#' View(myMetric.Values.Fish)
#
#' @export
metric.values <- function(fun.DF, fun.SampID, fun.Community, fun.MetricNames=NULL, boo.Adjust=FALSE){##FUNCTION.metric.values.START
  # convert community to lowercase
  fun.Community <- tolower(fun.Community)
  # run the proper sub function
  if (fun.Community=="bugs") {##IF.START
    metric.values.bugs(fun.DF, fun.SampID, fun.MetricNames, boo.Adjust)
  } else if(fun.Community=="fish"){
    metric.values.fish(fun.DF, fun.SampID, fun.MetricNames, boo.Adjust)
  } else if(fun.Community=="algae"){
    metric.values.algae(fun.DF, fun.SampID, fun.MetricNames, boo.Adjust)
  }##IF.END
}##FUNCTION.metric.values.START
#
#
#' @export
metric.values.bugs <- function(myDF, SampleID, MetricNames=NULL, boo.Adjust){##FUNCTION.metric.values.bugs.START
  # Remove Non-Target Taxa
  myDF <- myDF[myDF[,"NonTarget"]==0,]
  # Add extra columns for FFG and Habit (need unique values for functions in summarise)
  # each will be TRUE or FALSE
    myDF["Habit_BU"] <- grepl("BU",toupper(myDF[,"Habit"]))
    myDF["Habit_CB"] <- grepl("CB",toupper(myDF[,"Habit"]))
    myDF["Habit_CN"] <- grepl("CN",toupper(myDF[,"Habit"]))
    myDF["Habit_SP"] <- grepl("SP",toupper(myDF[,"Habit"]))
    myDF["Habit_SW"] <- grepl("SW",toupper(myDF[,"Habit"]))
    myDF["FFG_col"] <- grepl("col",tolower(myDF[,"FFG"]))
    myDF["FFG_fil"] <- grepl("fil",tolower(myDF[,"FFG"]))
    myDF["FFG_pre"] <- grepl("pre",tolower(myDF[,"FFG"]))
    myDF["FFG_scr"] <- grepl("scr",tolower(myDF[,"FFG"]))
    myDF["FFG_shr"] <- grepl("shr",tolower(myDF[,"FFG"]))
  # Calculate Metrics (could have used pipe, %>%)
  met.val <- dplyr::summarise(dplyr::group_by(myDF, SampleID, Index.Name, Index.Region)
             #
             # individuals, total
             ,ni_total=sum(Count)
             #
             # number of individuals
             ,ni_Ephem=sum(Count[Order=="Ephemeroptera"])
             ,ni_Trich=sum(Count[Order=="Trichoptera"])
             ,ni_Pleco=sum(Count[Order=="Plecoptera"])
             ,ni_EPT=sum(Count[Order=="Ephemeroptera" | Order=="Trichoptera" | Order=="Plecoptera"])
              #
             # percent individuals
             ,pi_Amph=sum(Count[Order=="Amphipoda"]) / ni_total
             ,pi_Bival=sum(Count[Class=="Bivalvia"]) / ni_total
             ,pi_Caen=sum(Count[Family=="Caenidae"]) / ni_total
             ,pi_Coleo=sum(Count[Order=="Coleoptera"]) / ni_total
             # Cole2Odon,
             # Colesensitive
             ,pi_Corb=sum(Count[Genus=="Corbicula"]) / ni_total
             #CruMol
             #Crus
             ,pi_Deca=sum(Count[Order=="Decapoda"]) / ni_total
             ,pi_Dipt=sum(Count[Order=="Diptera"]) / ni_total
             ,pi_Ephem=sum(Count[Order=="Ephemeroptera"]) / ni_total
             #EphemNoCaen
             #EPTsenstive
             ,pi_EPT=sum(Count[Order=="Ephemeroptera" | Order=="Trichoptera" | Order=="Plecoptera"]) / ni_total
             ,pi_Gast=sum(Count[Class=="Gastropoda"]) / ni_total
             ,pi_Iso=sum(Count[Order=="Isopoda"]) / ni_total
             #Moll
             ,pi_NonIns=sum(Count[Order!="Insecta" | is.na(Class)]) / ni_total
             ,pi_Odon=sum(Count[Order=="Odonata"]) / ni_total
             #oligo
             ,pi_Pleco=sum(Count[Order=="Plecoptera"]) / ni_total
             ,pi_Trich=sum(Count[Order=="Trichoptera"]) / ni_total
             ,pi_Tubif=sum(Count[Family=="Tubificidae"]) / ni_total
             #
             # number of taxa
              ,nt_total=dplyr::n_distinct(FinalID[NonUnique!=1])
              ,nt_Coleo=dplyr::n_distinct(FinalID[NonUnique!=1 & Order=="Coleoptera"])
             # ,nt_CruMol=dplyr::n_distinct(FinalID[NonUnique!=1 & (Phylum=="Mollusca" | SubPhylum="Crustacea")])
             ,nt_Dipt=dplyr::n_distinct(FinalID[NonUnique!=1 & Order=="Diptera"])
             ,nt_Ephem=dplyr::n_distinct(FinalID[NonUnique!=1 & Order=="Ephemeroptera"])
             ,nt_EPT=dplyr::n_distinct(FinalID[NonUnique!=1 & (Order=="Ephemeroptera"| Order=="Trichoptera" | Order=="Plecoptera")])
             ,nt_Oligo=dplyr::n_distinct(FinalID[NonUnique!=1 & Class=="Oligochaeta"])
             ,nt_Pleco=dplyr::n_distinct(FinalID[NonUnique!=1 & Order=="Plecoptera"])
             ,nt_Ptero=dplyr::n_distinct(FinalID[NonUnique!=1 & Genus=="Pteronarcys"])
             ,nt_Trich=dplyr::n_distinct(FinalID[NonUnique!=1 & Order=="Trichoptera"])
             # Amph, Bival, Gast, Deca, Insect, Isopod, intolMol, Oligo, POET, Tubif
             # intol
             #
             # Midges
             ,nt_Chiro=n_distinct(FinalID[NonUnique!=1 & Family=="Chironomidae"])
             ,pi_Chiro=sum(Count[Family=="Chironomidae"]) / ni_total
             #,pi_CrCh2Chi
             #,pi_Orth2Chi
             #,nt_Ortho
             #MB_pi_OrthocladiinaeCricotopusChironomus2Chironomidae
             ,pi_Tanyt=sum(Count[Tribe=="Tanytarsini"]) / ni_total
             #,pi-Tnyt2Chi,
             # COC2Chi
             # tanyp
             # tanyp2Chir
             #
             # percent of taxa
             # Amph, POET, Bival, Chiro, Deca, Dip, Gast, Iso, NonIns, Toler
             # / nt_total
             #
             # tolerance
             ,nt_tv_intol=n_distinct(FinalID[NonUnique!=1 & TolVal>=7])
             ,nt_tv_toler=n_distinct(FinalID[NonUnique!=1 & TolVal<=3])
             ,pi_tv_intolurb=sum(Count[TolVal>=7])/sum(Count[!is.na(TolVal)])
             # pi_Baet2Eph, pi_Hyd2EPT, pi_Hyd2Tri, pi_intol, pi_toler, , nt_intMol,
             # pt toler
             #
             # ffg
             ,nt_ffg_col=dplyr::n_distinct(FinalID[NonUnique!=1 & FFG_col==TRUE])
             ,nt_ffg_filt=dplyr::n_distinct(FinalID[NonUnique!=1 & FFG_fil==TRUE])
             ,nt_ffg_pred=dplyr::n_distinct(FinalID[NonUnique!=1 & FFG_pre==TRUE])
             ,nt_ffg_scrap=dplyr::n_distinct(FinalID[NonUnique!=1 & FFG_scr==TRUE])
             ,nt_ffg_shred=dplyr::n_distinct(FinalID[NonUnique!=1 & FFG_shr==TRUE])
             ,pi_ffg_col=dplyr::n_distinct(Count[NonUnique!=1 & FFG_col==TRUE]) / ni_total
             ,pi_ffg_filt=dplyr::n_distinct(Count[NonUnique!=1 & FFG_fil==TRUE]) / ni_total
             ,pi_ffg_pred=dplyr::n_distinct(Count[NonUnique!=1 & FFG_pre==TRUE]) / ni_total
             ,pi_ffg_scrap=dplyr::n_distinct(Count[NonUnique!=1 & FFG_scr==TRUE]) / ni_total
             ,pi_ffg_shred=dplyr::n_distinct(Count[NonUnique!=1 & FFG_shr==TRUE]) / ni_total
             # pt for cllct, filtr, pred, scrap, shred
              #
             # habit (need to be wild card)
             ,pi_habit_burrow=sum(Count[Habit_BU==TRUE]) / ni_total
             ,pi_habit_clmbrs=sum(Count[Habit_CB==TRUE]) / ni_total
             ,pi_habit_clngrs=sum(Count[Habit_CN==TRUE]) / ni_total
             ,pi_habit_sprawl=sum(Count[Habit_SP==TRUE]) / ni_total
             ,pi_habit_swmmrs=sum(Count[Habit_SW==TRUE]) / ni_total
             ,nt_habit_burrow=dplyr::n_distinct(FinalID[NonUnique!=1 & Habit_BU==TRUE])
             ,nt_habit_clmbrs=dplyr::n_distinct(FinalID[NonUnique!=1 & Habit_CB==TRUE])
             ,nt_habit_clngrs=dplyr::n_distinct(FinalID[NonUnique!=1 & Habit_CN==TRUE])
             ,nt_habit_sprawl=dplyr::n_distinct(FinalID[NonUnique!=1 & Habit_SP==TRUE])
             ,nt_habit_swmmrs=dplyr::n_distinct(FinalID[NonUnique!=1 & Habit_SW==TRUE])
             # pt for each
              #
             # voltinism
             # pi and nt for mltvol, semvol, univol
             ,pi_volt_multi=sum(Count[Voltinism=="multivoltine"]) / ni_total
             ,pi_volt_semi=sum(Count[Voltinism=="semivoltine"]) / ni_total
             ,pi_volt_uni=sum(Count[Voltinism=="univoltine"]) / ni_total
             ,nt_volt_multi=dplyr::n_distinct(FinalID[NonUnique!=1 & Voltinism=="multivoltine"])
             ,nt_volt_semi=dplyr::n_distinct(FinalID[NonUnique!=1 & Voltinism=="semivoltine"])
             ,nt_volt_uni=dplyr::n_distinct(FinalID[NonUnique!=1 & Voltinism=="univoltine"])
             #
             # indices
             ,pi_dom01=max(Count)/ni_total
             #,x_Becks.Class1=n_distinct(Count[NonUnique!=1 & TolVal>=0 & TolVal<=2.5])
             #,x_Becks.Class2=n_distinct(Count[NonUnique!=1 & TolVal>=2.5 & TolVal<=4])
             ,x_Becks=(2*dplyr::n_distinct(FinalID[NonUnique!=1 & TolVal>=0 & TolVal<=2.5]))+(1*dplyr::n_distinct(FinalID[NonUnique!=1 & TolVal>=2.5 & TolVal<=4]))
             #,x_HBI_num=sum(Count*TolVal)
             #,x_HBI_denom=sum(Count[!is.na(TolVal) & TolVal>0])
             ,x_HBI=sum(Count*TolVal)/sum(Count[!is.na(TolVal) & TolVal>0])
           #  ,x_Shan_Num=log(3.14)
          #   ,x_Shan_e=x_Shan_Num/log(exp(1))
          #   ,x_Shan_2=x_Shan_Num/log(2)
          #   ,x_Shan_10=x_Shan_Num/log(10)
             #, x_D Simpson
             #, x_Hbe
             #, x_D_Mg Margalef
             #, x_H
             # Pielou
              # H / Hmax  Hmax is log(nt_total)
             #
             # BCG
             ,nt_BCG_att123=dplyr::n_distinct(FinalID[NonUnique!=1 & (BCG_Atr=="1" | BCG_Atr=="2" | BCG_Atr=="3")])
             # nt_att 12, 123, 2, 23, 234, 4, 5, 5, 56
             # nt_EPT_att123
             # pi_att 12, 123, 23, 45, 5, 56
             # pi_dom01_att 4, 5, 56
             # pi_dom05_att 123, not 456
             # pi_EPT_att123
             # pt_att 12, 123, 23, 234, 5, 56
             # pt_EPT_att 123
             #
          )## met.val.END
  # replace NA with 0
  met.val[is.na(met.val)] <- 0
  # # subset to only metrics specified by user
  if (!is.null(MetricNames)){
    met.val <- met.val[,c(SampleID,"Index.Name","Index.Region", "ni_total", MetricNames)]
  }
  # df to report back
  return(met.val)
}##FUNCTION.metric.values.bugs.END
#
#
#' @export
metric.values.fish <- function(myDF, SampleID, MetricNames=NULL, boo.Adjust){##FUNCTION.metric.values.fish.START
  # Remove Non-Target Taxa
  myDF <- myDF[myDF[,"NonTarget"]==0,]
  # set case on fields
  myFlds <- c("FinalID", "Genus", "Trophic", "Stratum", "Guild", "Tolerance")
  for (i in myFlds) {
    myDF[,i] <- tolower(myDF[,i])
  }
  # Calculate Metrics (could have used pipe, %>%)
  met.val <- dplyr::summarise(dplyr::group_by(myDF, SampleID, Index.Name, Index.Region, Catchment)
                       #
                       # MBSS 2005, 11 metrics
                       #
                       # individuals, total
                       ,ni_total=sum(Count)
                       #
                       # number of individuals
                       # Num Benthic Species

                       #,ni_Ephem=sum(Count[Order=="Ephemeroptera"])
                       #,ni_Trich=sum(Count[Order=="Trichoptera"])
                       #,ni_Pleco=sum(Count[Order=="Plecoptera"])
                       #,ni_EPT=sum(Count[Order=="Ephemeroptera" | Order=="Trichoptera" | Order=="Plecoptera"])
                       #
                       # percent individuals

                       # % RBS
                       ,pi_rbs=sum(Count[Genus=="hypentelium"|Genus=="moxostoma"|Genus=="minytrema"|Genus=="erimyzon"])/ni_total
                       # Pct Brook Trout
                       ,pi_brooktrout=sum(Count[FinalID=="brook trout"])/ni_total
                       # Pct Sculpins
                       ,pi_sculpin=sum(Count[Genus=="cottus"|Genus=="myoxocephalus"])/ni_total


                       # ,pi_Amph=sum(Count[Order=="Amphipoda"]) / ni_total
                       # ,pi_Bival=sum(Count[Class=="Bivalvia"]) / ni_total
                       # ,pi_Caen=sum(Count[Family=="Caenidae"]) / ni_total
                       # ,pi_Coleo=sum(Count[Order=="Coleoptera"]) / ni_total
                       #
                       #
                       # ,pi_Corb=sum(Count[Genus=="Corbicula"]) / ni_total
                       #
                       #
                       # ,pi_Deca=sum(Count[Order=="Decapoda"]) / ni_total
                       # ,pi_Dipt=sum(Count[Order=="Diptera"]) / ni_total
                       # ,pi_Ephem=sum(Count[Order=="Ephemeroptera"]) / ni_total
                       #
                       #
                       # ,pi_EPT=sum(Count[Order=="Ephemeroptera" | Order=="Trichoptera" | Order=="Plecoptera"]) / ni_total
                       # ,pi_Gast=sum(Count[Class=="Gastropoda"]) / ni_total
                       # ,pi_Iso=sum(Count[Order=="Isopoda"]) / ni_total
                       #
                       # ,pi_NonIns=sum(Count[Order!="Insecta" | is.na(Class)]) / ni_total
                       # ,pi_Odon=sum(Count[Order=="Odonata"]) / ni_total
                       #
                       # ,pi_Pleco=sum(Count[Order=="Plecoptera"]) / ni_total
                       # ,pi_Trich=sum(Count[Order=="Trichoptera"]) / ni_total
                       # ,pi_Tubif=sum(Count[Family=="Tubificidae"]) / ni_total
                       #
                       #
                       # number of taxa
                       ,nt_total=dplyr::n_distinct(FinalID[NonUnique!=1])
                       ,nt_benthic=dplyr::n_distinct(FinalID[NonUnique!=1 & Stratum=="b"])
                       # ,nt_Coleo=n_distinct(Count[NonUnique!=1 & Order=="Coleoptera"])
                       # #,nt_CruMol=n_distinct(Count[NonUnique!=1 & (Phylum=="Mollusca" | Subphylum="Crustacea")])
                       # ,nt_Ephem=n_distinct(Count[NonUnique!=1 & Order=="Ephemeroptera"])
                       # ,nt_EPT=n_distinct(Count[NonUnique!=1 & (Order=="Ephemeroptera"| Order=="Trichoptera" | Order=="Plecoptera")])
                       # ,nt_Oligo=n_distinct(Count[NonUnique!=1 & Class=="Oligochaeta"])
                       # ,nt_Pleco=n_distinct(Count[NonUnique!=1 & Order=="Plecoptera"])
                       # ,nt_Ptero=n_distinct(Count[NonUnique!=1 & Genus=="Pteronarcys"])
                       # ,nt_Trich=n_distinct(Count[NonUnique!=1 & Order=="Trichoptera"])
                       # # Amph, Bival, Gast, Deca, Insect, Isopod, intolMol, Oligo, POET, Tubif
                       # intol
                       # #
                       # # Midges
                       # ,nt_Chiro=n_distinct(Count[NonUnique!=1 & Family=="Chironomidae"])
                       # ,pi_Chiro=sum(Count[Family=="Chironomidae"]) / ni_total
                       # #MB_pi_OrthocladiinaeCricotopusChironomus2Chironomidae
                       # ,pi_Tanyt=sum(Count[Tribe=="Tanytarsini"]) / ni_total
                       #
                       # percent of taxa
                       # Amph, POET, Bival, Chiro, Deca, Dip, Gast, Iso, NonIns, Toler
                       # / nt_total
                       #
                       # tolerance
                       # % Tolerant
                       ,pi_tv_toler=sum(Count[Tolerance=="tol"])/ni_total


                       #,pi_tv_intolurb=sum(Count[TolVal>=7])/sum(Count[!is.na(TolVal)])
                       # pi_Baet2Eph, pi_Hyd2EPT, pi_Hyd2Tri, pi_intol, pi_toler, nt_intol, nt_intMol, nt_toler
                       # pt toler
                       #
                       # ffg

                       # Feeding
                       # % Lithophilic spawners
                       ,pi_lithophil=sum(Count[Lithophil=="Yes"])/ni_total
                       # % gen, omn, invert
                       ,pi_genomninvrt=sum(Count[Trophic=="generalist" | Trophic=="omnivore" | Trophic=="invertivore"])/ ni_total
                       # % insectivore
                       ,pi_insectivore=sum(Count[Trophic=="insectivore"])/ ni_total


                       #
                       # ,nt_ffg_scrap=n_distinct(Count[NonUnique!=1 & FFG=="SC"])
                       # ,pi_ffg_scrap=n_distinct(Count[NonUnique!=1 & FFG=="SC"]) / ni_total
                       # # pi and nt for cllct, filtr, pred,scrap, shred
                       #
                       # # habit (need to be wild card)
                       # ,pi_habit_burrow=sum(Count[Habit=="BU"]) / ni_total
                       # ,pi_habit_clngrs=sum(Count[Habit=="CN"]) / ni_total
                       # ,pi_habit_clmbrs=sum(Count[Habit=="CB"]) / ni_total
                       # ,pi_habit_sprawl=sum(Count[Habit=="SP"]) / ni_total
                       # ,pi_habit_swmmrs=sum(Count[Habit=="SW"]) / ni_total
                       # ,nt_habit_burrow=n_distinct(Count[NonUnique!=1 & Habit=="BU"])
                       # ,nt_habit_clngrs=n_distinct(Count[NonUnique!=1 & Habit=="CN"])
                       # ,nt_habit_clmbrs=n_distinct(Count[NonUnique!=1 & Habit=="CB"])
                       # ,nt_habit_sprawl=n_distinct(Count[NonUnique!=1 & Habit=="SP"])
                       # ,nt_habit_swmmrs=n_distinct(Count[NonUnique!=1 & Habit=="SW"])
                       # # pt Swmmr
                       # #
                       # # voltinism
                       # # pi and nt for mltvol, semvol, univol
                       # ,pi_volt_multi=sum(Count[Voltinism=="multivoltine"]) / ni_total
                       # ,pi_volt_semi=sum(Count[Voltinism=="semivoltine"]) / ni_total
                       # ,pi_volt_uni=sum(Count[Voltinism=="univoltine"]) / ni_total
                       # ,nt_volt_multi=n_distinct(Count[NonUnique!=1 & Voltinism=="multivoltine"])
                       # ,nt_volt_semi=n_distinct(Count[NonUnique!=1 & Voltinism=="semivoltine"])
                       # ,nt_volt_uni=n_distinct(Count[NonUnique!=1 & Voltinism=="univoltine"])
                       #
                       # indices
                       #,pi_dom01/2/3/5 #last? or nth
                       ,pi_dom01=max(Count)/ni_total

                       # Other
                       ,area=max(StWidAvg)*max(StLength)
                       # Abund / sq meter
                       ,ni_m2=ni_total/area #/(StWidAvg*StLength)
                       # biomass per square meter
                       ,x_biomass_m2=max(MassTotal)/area #/(StWidAvg*StLength)
                       # #
                       # # BCG
                       # ,nt_BCG_att123=n_distinct(Count[NonUnique!=1 & (BCG_Atr=="1" | BCG_Atr=="2" | BCG_Atr=="3")])
                       #
  )## met.val.END
  # replace NA with 0
  met.val[is.na(met.val)] <- 0
  # subset to only metrics specified by user
  if (!is.null(MetricNames)){
    met.val <- met.val[,c(SampleID, "Index.Name", "Index.Region", "Catchment", "ni_total", MetricNames)]
  }
  # Adjust metrics
  if (boo.Adjust==TRUE) {##IF.boo.Ajust.START
    # MBSS.2005.Fish
    # nt_benthic
      met.val[,"nt_benthic_Obs"] <- met.val[,"nt_benthic"]
      # Expected constants
      ## m
      met.val[,"nt_benthic_m"] <- NA
      met.val[,"nt_benthic_m"][met.val[,"Index.Region"]=="COASTAL"]   <- 1.69
      met.val[,"nt_benthic_m"][met.val[,"Index.Region"]=="EPIEDMONT"] <- 1.25
      met.val[,"nt_benthic_m"][met.val[,"Index.Region"]=="HIGHLAND"]  <- 1.23
      ## b
      met.val[,"nt_benthic_b"] <- NA
      met.val[,"nt_benthic_b"][met.val[,"Index.Region"]=="COASTAL"]   <- -3.33
      met.val[,"nt_benthic_b"][met.val[,"Index.Region"]=="EPIEDMONT"] <- -2.36
      met.val[,"nt_benthic_b"][met.val[,"Index.Region"]=="HIGHLAND"] <- -2.35
      # Calc Expected
      met.val[,"nt_benthic_Exp"] <- (met.val[,"nt_benthic_m"] * log10(met.val[,"Catchment"])) + met.val[,"nt_benthic_b"]
      # Calc Adjusted
      met.val[,"nt_benthic_Adj"] <- met.val[,"nt_benthic_Obs"] / met.val[,"nt_benthic_Exp"]
      # Rename base metric with adjusted value
      met.val[,"nt_benthic"] <- met.val[,"nt_benthic_Adj"]

  }##IF.boo.Ajust.END
  #
  # df to report back
  return(met.val)
}##FUNCTION.metric.values.fish.END
#
#
#' @export
metric.values.algae <- function(myDF, SampleID, MetricNames=NULL, boo.Adjust){##FUNCTION.metric.values.algae.START
  # Calculate Metrics (could have used pipe, %>%)
    met.val <- dplyr::summarise(dplyr::group_by(myDF, SampleID, "Index.Name", "Index.Region")
                #
                # individuals, total
                ,ni_total=sum(Count)
                #
    )##met.val.END
    # replace NA with 0
    met.val[is.na(met.val)] <- 0
    # subset to only metrics specified by user
    if (!is.null(MetricNames)){
      met.val <- met.val[,c(SampleID, "Index.Name", "Index.Region", "ni_total", MetricNames)]
    }
    # df to report back
    return(met.val)
}##FUNCTION.metric.values.algae.END
