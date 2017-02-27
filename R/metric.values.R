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
#
#' @return data frame of SampleID and metric values
#
#' @examples
#' myDF.Bugs <- read.delim("bugs.agg.tab",stringsAsFactors = FALSE)
#' myMetric.Values.Bugs <- metric.values(myDF,SampleID,"bugs")
#' View(myMetric.Values.Bugs)
#
#
#' @export
metric.values <- function(fun.DF,fun.SampID,fun.Community,fun.MetricNames=NULL){##FUNCTION.metric.values.START
  # convert community to lowercase
  fun.Community <- tolower(fun.Community)
  # run the proper sub function
  if (fun.Community=="bugs") {##IF.START
    metric.values.bugs(fun.DF,fun.SampID,fun.MetricNames)
  } else if(fun.Community=="fish"){
    metric.values.fish(fun.DF,fun.SampID,fun.MetricNames)
  } else if(fun.Community=="algae"){
    metric.values.algae(fun.DF,fun.SampID,fun.MetricNames)
  }##IF.END
}##FUNCTION.metric.values.START
#
#
#' @export
metric.values.bugs <- function(myDF,SampleID,MetricNames=NULL){##FUNCTION.metric.values.bugs.START
  # Remove Non-Target Taxa
  myDF <- myDF[myDF[,"NonTarget"]==0,]
  # Calculate Metrics
  met.val <- dplyr::summarise(dplyr::group_by(myDF,SampleID,IndexName,IndexRegion)
             #
             # individuals, total
             ,ni_total=sum(Count)
             #
             # number of individuals
             #,ni_Ephem=sum(Count[Order=="Ephemeroptera"])
             #,ni_Trich=sum(Count[Order=="Trichoptera"])
             #,ni_Pleco=sum(Count[Order=="Plecoptera"])
             #,ni_EPT=sum(Count[Order=="Ephemeroptera" | Order=="Trichoptera" | Order=="Plecoptera"])
             #
             # percent individuals
             ,pi_Amph=sum(Count[Order=="Amphipoda"]) / ni_total
             ,pi_Bival=sum(Count[Class=="Bivalvia"]) / ni_total
             ,pi_Caen=sum(Count[Family=="Caenidae"]) / ni_total
             ,pi_Coleo=sum(Count[Order=="Coleoptera"]) / ni_total
             
             
             ,pi_Corb=sum(Count[Genus=="Corbicula"]) / ni_total
             
             
             ,pi_Deca=sum(Count[Order=="Decapoda"]) / ni_total
             ,pi_Dipt=sum(Count[Order=="Diptera"]) / ni_total
             ,pi_Ephem=sum(Count[Order=="Ephemeroptera"]) / ni_total
             
             
             ,pi_EPT=sum(Count[Order=="Ephemeroptera" | Order=="Trichoptera" | Order=="Plecoptera"]) / ni_total
             ,pi_Gast=sum(Count[Class=="Gastropoda"]) / ni_total
             ,pi_Iso=sum(Count[Order=="Isopoda"]) / ni_total
             
             ,pi_NonIns=sum(Count[Order!="Insecta" | is.na(Class)]) / ni_total
             ,pi_Odon=sum(Count[Order=="Odonata"]) / ni_total
             
             ,pi_Pleco=sum(Count[Order=="Plecoptera"]) / ni_total
             ,pi_Trich=sum(Count[Order=="Trichoptera"]) / ni_total
             ,pi_Tubif=sum(Count[Family=="Tubificidae"]) / ni_total
             # Cole2Odon, Colesensitive, CruMol, Crus, EphemNoCaen, EPTsenstive,
             # intol, Moll, oligo
             # 
             #
             # number of taxa
             ,nt_total=dplyr::n_distinct(FinalID[NonUnique!=1])
             ,nt_Coleo=dplyr::n_distinct(FinalID[NonUnique!=1 & Order=="Coleoptera"])
             #,nt_CruMol=n_distinct(FinalID[NonUnique!=1 & (Phylum=="Mollusca" | Subphylum="Crustacea")])
             ,nt_Ephem=dplyr::n_distinct(FinalID[NonUnique!=1 & Order=="Ephemeroptera"])
             ,nt_EPT=dplyr::n_distinct(FinalID[NonUnique!=1 & (Order=="Ephemeroptera"| Order=="Trichoptera" | Order=="Plecoptera")])
             ,nt_Oligo=dplyr::n_distinct(FinalID[NonUnique!=1 & Class=="Oligochaeta"])
             ,nt_Pleco=dplyr::n_distinct(FinalID[NonUnique!=1 & Order=="Plecoptera"])
             ,nt_Ptero=dplyr::n_distinct(FinalID[NonUnique!=1 & Genus=="Pteronarcys"])
             ,nt_Trich=dplyr::n_distinct(FinalID[NonUnique!=1 & Order=="Trichoptera"])
             # Amph, Bival, Gast, Deca, Insect, Isopod, intolMol, Oligo, POET, Tubif
             # intol
             # #
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
             ,pi_tv_intolurb=sum(Count[TolVal>=7])/sum(Count[!is.na(TolVal)])
             # pi_Baet2Eph, pi_Hyd2EPT, pi_Hyd2Tri, pi_intol, pi_toler, nt_intol, nt_intMol, nt_toler
             # pt toler
             #
             # ffg
             ,nt_ffg_scrap=dplyr::n_distinct(FinalID[NonUnique!=1 & FFG=="SC"])
             ,pi_ffg_scrap=dplyr::n_distinct(Count[NonUnique!=1 & FFG=="SC"]) / ni_total
             # pi and nt for cllct, filtr, pred,scrap, shred
             #
             # habit (need to be wild card)
             ,pi_habit_burrow=sum(Count[Habit=="BU"]) / ni_total
             ,pi_habit_clngrs=sum(Count[Habit=="CN"]) / ni_total
             ,pi_habit_clmbrs=sum(Count[Habit=="CB"]) / ni_total
             ,pi_habit_sprawl=sum(Count[Habit=="SP"]) / ni_total
             ,pi_habit_swmmrs=sum(Count[Habit=="SW"]) / ni_total
             ,nt_habit_burrow=dplyr::n_distinct(FinalID[NonUnique!=1 & Habit=="BU"])
             ,nt_habit_clngrs=dplyr::n_distinct(FinalID[NonUnique!=1 & Habit=="CN"])
             ,nt_habit_clmbrs=dplyr::n_distinct(FinalID[NonUnique!=1 & Habit=="CB"])
             ,nt_habit_sprawl=dplyr::n_distinct(FinalID[NonUnique!=1 & Habit=="SP"])
             ,nt_habit_swmmrs=dplyr::n_distinct(FinalID[NonUnique!=1 & Habit=="SW"])
             # pt Swmmr
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
             ,x_Shan_Num=log(3.14)
             ,x_Shan_e=x_Shan_Num/log(exp(1))
             ,x_Shan_2=x_Shan_Num/log(2)
             ,x_Shan_10=x_Shan_Num/log(10)
             #, x_D Simpson
             #, x_Hbe
             #, x_D_Mg Margalef
             #, x_H
             # Pielou
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
  # subset to only metrics specified by user
  if (!is.null(MetricNames)){
    met.val <- met.val[,c(SampleID,"IndexName","IndexRegion",MetricNames)]
  }
  # df to report back
  return(met.val)
}##FUNCTION.metric.values.bugs.END
#
#
#' @export
metric.values.fish <- function(myDF,SampleID,MetricNames=NULL){##FUNCTION.metric.values.fish.START
  # Remove Non-Target Taxa
  myDF <- myDF[myDF[,"NonTarget"]==0,]
  # Calculate Metrics
  met.val <- dplyr::summarise(dplyr::group_by(myDF,SampleID)
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
                       ,pi_rbs=sum(Count[Genus=="Hypentelium"|Genus=="Moxostoma"|Genus=="Minytrema"|Genus=="Erimyzon"])/ni_total
                       # Pct Brook Trout
                       ,pi_BrkTrt=sum(Count[FinalID=="Brook Trout"])/ni_total
                       # Pct Sculpins
                       ,pi_sculpin=sum(Count[Genus=="Cottus"|Genus=="Myoxocephalus"])/ni_total
                       
                       
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
                       # nt_Benthic
                       ,nt_benthic=dplyr::n_distinct(FinalID[Stratum=="B"])
                       ,nt_total=dplyr::n_distinct(FinalID[NonUnique!=1])
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
                       ,pi_toler=sum(Count[Tolerance=="Tol"])/ni_total
                       
                       
                       #,pi_tv_intolurb=sum(Count[TolVal>=7])/sum(Count[!is.na(TolVal)])
                       # pi_Baet2Eph, pi_Hyd2EPT, pi_Hyd2Tri, pi_intol, pi_toler, nt_intol, nt_intMol, nt_toler
                       # pt toler
                       #
                       # ffg
                       
                       # Feeding
                       # % Lithophilic spawners
                       ,pi_lith=sum(Count[Lithophil=="Yes"])/ni_total
                       # % gen, omn, invert
                       ,pi_goi=sum(Count[Trophic=="Generalist" | Trophic=="Omnivore" | Trophic=="Invertivore"])/ ni_total
                       # % insectivore
                       ,pi_insectivore=sum(Count[Trophic=="Insectivore"])/ ni_total
                       
                       
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
                       ,ind_m2=ni_total/area #/(StWidAvg*StLength)
                       # biomass per square meter
                       ,biomass_m2=max(MassTotal)/area #/(StWidAvg*StLength)
                       # #
                       # # BCG
                       # ,nt_BCG_att123=n_distinct(Count[NonUnique!=1 & (BCG_Atr=="1" | BCG_Atr=="2" | BCG_Atr=="3")])
                       #
  )## met.val.END
  # replace NA with 0
  met.val[is.na(met.val)] <- 0
  # subset to only metrics specified by user
  if (!is.null(MetricNames)){
    met.val <- met.val[,c(SampleID,MetricNames)]
  }
  # df to report back
  return(met.val)
}##FUNCTION.metric.values.fish.END
# 
#
#' @export
metric.values.algae <- function(myDF,SampleID){##FUNCTION.metric.values.algae.START
    met.val <- dplyr::summarise(dplyr::group_by(myDF,SampleID)
                #
                # individuals, total
                ,ni_total=sum(Count)
                #
    )##met.val.END
    # replace NA with 0
    met.val[is.na(met.val)] <- 0
    # subset to only metrics specified by user
    if (!is.null(MetricNames)){
      met.val <- met.val[,c(SampleID,MetricNames)]
    }
    # df to report back
    return(met.val)
}##FUNCTION.metric.values.algae.END