rm(list=ls())
setwd('/Users/edmondlowjy/datasciencecoursera/Reproducible-Research/Week\ 4')

download.file(url='https://d396qusza40orc.cloudfront.net/repdata%2Fdata%2FStormData.csv.bz2',destfile='./stormdata.csv.bz2',method='curl')

library(dplyr); library(ggplot2)
stormdata=read.csv(file='stormdata.csv.bz2',stringsAsFactors=F); beepr::beep()
stormdata<-tbl_df(stormdata) %>%
  mutate(EVTYPE=trimws(toupper(EVTYPE)))
evnames=unique(stormdata$EVTYPE)

event_types=c("ASTRONOMICAL LOW TIDE","AVALANCHE","BLIZZARD","COASTAL FLOOD",
              "COLD/WIND CHILL","DEBRIS FLOW","DENSE FOG","DENSE SMOKE","DROUGHT",
              "DUST DEVIL","DUST STORM","EXCESSIVE HEAT","EXTREME COLD/WIND CHILL",
              "FLASH FLOOD","FLOOD","FROST/FREEZE","FUNNEL CLOUD","FREEZING FOG","HAIL",
              "HEAT","HEAVY RAIN","HEAVY SNOW","HIGH SURF","HIGH WIND","HURRICANE (TYPHOON)",
              "ICE STORM","LAKE-EFFECT SNOW","LAKESHORE FLOOD","LIGHTNING","MARINE HAIL",
              "MARINE HIGH WIND","MARINE STRONG WIND","MARINE THUNDERSTORM WIND",
              "RIP CURRENT","SEICHE","SLEET","STORM SURGE/TIDE","STRONG WIND",
              "THUNDERSTORM WIND","TORNADO","TROPICAL DEPRESSION","TROPICAL STORM","TSUNAMI",
              "VOLCANIC ASH","WATERSPOUT","WILDFIRE","WINTER STORM","WINTER WEATHER")

### First attempt to find exact matches using grep - finding EVTYPE within event_types

textmatcher=function(textvec,refvec){ #returns event_type of min char length that matches exactly 
  temp=lapply(textvec,FUN=function(text){
    matchvec=grep(text,refvec,fixed = T,value = T)
    indexselect=which.min(nchar(matchvec))
    matchvec[indexselect]
    })
  return(temp)
}

# A list of event_type matches for each unique stormdata$EVTYPE after 1st attempt
textmatched1=textmatcher(evnames,event_types); names(textmatched1)<-evnames

# A vector of (unique) unmatched elements from stormdata$EVTYPE after 1st matching attempt
notmatched=evnames[sapply(textmatched1,FUN=function(textlist){length(unlist(textlist))==0})]

### Second attempt to find exact matches using grep - finding event_types within EVTYPE

textmatcher2=function(textvec,refvec){
  templist=lapply(textvec,FUN=function(text){grep(text,refvec,fixed=T)})
  reslist=as.list(rep(NA,times=length(refvec)))
  for(i in 1:length(templist)){ # 1:48
    if(length(templist[[i]]!=0)){ # if some EVTYPES are found for selected event_type
      for(j in templist[[i]]){ # for element in vector of matches for selected event_type
        reslist[[j]]=c(reslist[[j]],textvec[i])
      }
    }
  }
  names(reslist)=refvec
  reslist=lapply(reslist,FUN=function(vec){vec[which.max(nchar(vec))]})
}

# A list of event_type matches for each unique stormdata$EVTYPE after 2nd attempt
textmatched2=textmatcher2(event_types,notmatched)

# A vector of (unique) unmatched elements from stormdata$EVTYPE after 2nd matching attempt
notmatched=notmatched[sapply(textmatched2,FUN=function(textlist){length(unlist(textlist))==0})]

# Extracting the matches from the 1st and 2nd lists
EVTYPE_match=NULL
Event_match=NULL
for(i in 1:length(textmatched1)){
  if(length(textmatched1[[i]]!=0)){
    EVTYPE_match=c(EVTYPE_match,names(textmatched1[i]))
    Event_match=c(Event_match,textmatched1[[i]])
  }
}
for(i in 1:length(textmatched2)){
  if(length(textmatched2[[i]]!=0)){
    EVTYPE_match=c(EVTYPE_match,names(textmatched2[i]))
    Event_match=c(Event_match,textmatched2[[i]])
  }
}

# Replacing the EVTYPES with the 'official' event_types

stormdata<-mutate(stormdata,EVTYPE_mod=EVTYPE)
for(i in 1:length(Event_match)){
  #this is wrong because grep and gsub takes partial matches e.g. WIND in TSTM WIND throws a match
  #stormdata<-mutate(stormdata,EVTYPE_mod=replace(EVTYPE_mod,grep(EVTYPE_match[i],EVTYPE_mod,fixed=T),Event_match[i]))
  stormdata<-mutate(stormdata,EVTYPE_mod=replace(EVTYPE_mod,which(EVTYPE_mod==EVTYPE_match[i]),Event_match[i]))
}

#which/how many EVTYPES are found in official event_types?
sum(stormdata$EVTYPE %in% event_types)#before work above
sum(stormdata$EVTYPE_mod %in% event_types)#after work above
sum(stormdata$EVTYPE_mod %in% event_types)/nrow(stormdata) #percentage coverage

# Takes an input vector of names and a reference vector, outputs a
# vector of names that needs to be matched in order to achieve x% coverage
topx=function(textvec,refvec,coverage=0.9){
  total=length(textvec); curr=sum(textvec %in% refvec)/total
  textvec=textvec[which(!(textvec %in% refvec))] #taking what does not match
  EVTYPE_counts=sort(table(textvec),decreasing = T) #table frequency of unmatches names
  i=0
  while(curr<coverage){
    i=i+1
    curr=curr+(EVTYPE_counts[i]/total)
  }
  print(paste0('Percentage Coverage: ',round(curr*100,1),'%'))
  return(EVTYPE_counts[1:i])
}
topx(stormdata$EVTYPE_mod,event_types,0.995)

#Manual Replacements to match at least 99.5% of EVTYPES
stormdata<-stormdata %>%
  mutate(EVTYPE_mod=replace(EVTYPE_mod,which(EVTYPE_mod=='TSTM WIND'),'THUNDERSTORM WIND')) %>%
  mutate(EVTYPE_mod=replace(EVTYPE_mod,which(EVTYPE_mod=='MARINE TSTM WIND'),'MARINE THUNDERSTORM WIND')) %>%
  mutate(EVTYPE_mod=replace(EVTYPE_mod,which(EVTYPE_mod=='URBAN/SML STREAM FLD'),'HEAVY RAIN')) %>%
  mutate(EVTYPE_mod=replace(EVTYPE_mod,which(EVTYPE_mod=='WILD/FOREST FIRE'),'WILDFIRE'))
sum(stormdata$EVTYPE_mod %in% event_types)/nrow(stormdata) #percentage coverage

EVTYPE_match=c(EVTYPE_match,'TSTM WIND','MARINE TSTM WIND','URBAN/SML STREAM FLD','WILD/FOREST FIRE')
Event_match=c(Event_match,'THUNDERSTORM WIND','MARINE THUNDERSTORM WIND','HEAVY RAIN','WILDFIRE')
match_table=data.frame(EVTYPE_match,Event_match)
notmatched=notmatched[!(notmatched %in% EVTYPE_match)]

rm(list=setdiff(ls(), c("stormdata","event_types","match_table","notmatched")))

### Analyzing health impact of different EVTYPES

health_impact<-stormdata %>%
  select(EVTYPE_mod,FATALITIES,INJURIES) %>%
  group_by(EVTYPE_mod) %>%
  summarize(FATALITIES=sum(FATALITIES),INJURIES=sum(INJURIES),TOTAL=FATALITIES+INJURIES) %>%
  arrange(desc(TOTAL),desc(FATALITIES),desc(INJURIES))

# Categorize the different EVTYPES according to the official event_types to ensure 
# 99.5% coverage in terms of health impact

cum_health=sum(health_impact$TOTAL)
curr=0; i=0; coverage=0.995
while(curr<coverage){
  i=i+1
  curr=curr+health_impact$TOTAL[i]/cum_health
}
reqvec=health_impact$EVTYPE_mod[1:i]
reqvec[which(!(reqvec %in% event_types))]

stormdata<-stormdata %>%
  mutate(EVTYPE_mod=replace(EVTYPE_mod,grep('HURRICANE',EVTYPE_mod),'HURRICANE (TYPHOON)')) %>%
  mutate(EVTYPE_mod=replace(EVTYPE_mod,grep('GLAZE',EVTYPE_mod),'FROST/FREEZE')) %>%
  mutate(EVTYPE_mod=replace(EVTYPE_mod,grep('WILD',EVTYPE_mod),'WILDFIRE')) %>%
  mutate(EVTYPE_mod=replace(EVTYPE_mod,grep('LANDSLIDE|LANDSLUMP',EVTYPE_mod),'DEBRIS FLOW')) %>%
  mutate(EVTYPE_mod=replace(EVTYPE_mod,grep('WIN(.*) MIX',EVTYPE_mod),'WINTER WEATHER'))

sum(stormdata$EVTYPE_mod %in% event_types)/nrow(stormdata) #percentage coverage

# Adding to the table of matched EVTYPES

EVTYPE_match=NULL; Event_match=NULL
regexpvec=c('HURRICANE','GLAZE','WILD','LANDSLIDE|LANDSLUMP','WIN(.*) MIX')
replacevec=c('HURRICANE (TYPHOON)','FROST/FREEZE','WILDFIRE','DEBRIS FLOW','WINTER WEATHER')
for(i in 1:5){
  tempvec=health_impact$EVTYPE_mod[grep(regexpvec[i],health_impact$EVTYPE_mod)]
  tempvec=tempvec[(!tempvec %in% event_types)]
  EVTYPE_match=c(EVTYPE_match,tempvec)
  Event_match=c(Event_match,rep(replacevec[i],length(tempvec)))
  notmatched=notmatched[!notmatched %in% tempvec]
}
match_table=rbind(match_table,data.frame(EVTYPE_match,Event_match))

rm(list=setdiff(ls(), c("stormdata","event_types","match_table","notmatched")))

### Analyzing health impact

health_impact<-stormdata %>%
  select(EVTYPE_mod,FATALITIES,INJURIES) %>%
  filter(EVTYPE_mod %in% event_types) %>%
  mutate(TOTAL=FATALITIES+INJURIES) %>%
  group_by(EVTYPE_mod) %>%
  summarize(IMPACT_MEAN=mean(TOTAL),FATALITIES=sum(FATALITIES),INJURIES=sum(INJURIES),TOTAL=sum(TOTAL)) %>%
  arrange(desc(TOTAL),desc(FATALITIES),desc(INJURIES))

g<-ggplot(health_impact,aes(EVTYPE_mod,TOTAL))
g+geom_bar(stat='identity')+theme(axis.text.x = element_text(angle = 90, hjust = 1))+labs(x='Event Type',y='Total Health Impact',title='Total Health Impact of each Event Type summed across the years')
qplot(EVTYPE_mod,IMPACT_MEAN,data=health_impact)+theme(axis.text.x = element_text(angle = 90, hjust = 1))+labs(x='Event Type',y='Mean Health Impact',title='Mean Health Impact per instance of each Event Type')


###### Analyzing Economic Damage

economic_impact<-stormdata %>%
  select(EVTYPE_mod,PROPDMG,PROPDMGEXP,CROPDMG,CROPDMGEXP) %>%
  filter(EVTYPE_mod %in% event_types) %>%
  mutate(PROPDMGEXP=trimws(toupper(PROPDMGEXP)),CROPDMGEXP=trimws(toupper(CROPDMGEXP))) %>%
  filter(PROPDMGEXP==''|PROPDMGEXP=='K'|PROPDMGEXP=='M'|PROPDMGEXP=='B') %>%
  filter(CROPDMGEXP==''|CROPDMGEXP=='K'|CROPDMGEXP=='M'|CROPDMGEXP=='B') %>%
  mutate(PROPDMGEXP_mod=as.numeric(sapply(PROPDMGEXP,FUN=function(item){
    if(item=='K'){item=1000}
    else if(item==''){item=1}
    else if(item=='M'){item=1000000}
    else if (item=='B'){item=1000000000}}))) %>%
  mutate(CROPDMGEXP_mod=as.numeric(sapply(CROPDMGEXP,FUN=function(item){
    if(item=='K'){item=1000}
    else if(item==''){item=1}
    else if(item=='M'){item=1000000}
    else if (item=='B'){item=1000000000}}))) %>%
  mutate(PROPDMG_mod=PROPDMG*PROPDMGEXP_mod,CROPDMG_mod=CROPDMG*CROPDMGEXP_mod,TOTALDMG=PROPDMG_mod+CROPDMG_mod) %>%
  group_by(EVTYPE_mod) %>%
  summarize(MEANDMG=mean(TOTALDMG),PROPDMG=sum(PROPDMG_mod),CROPDMG=sum(CROPDMG_mod),TOTALDMG=sum(TOTALDMG)) %>%
  arrange(desc(TOTALDMG),desc(PROPDMG),desc(CROPDMG))

g<-ggplot(economic_impact,aes(EVTYPE_mod,TOTALDMG))
g+geom_bar(stat='identity')+theme(axis.text.x = element_text(angle = 90, hjust = 1))+labs(x='Event Type',y='Economic Impact',title='Total Economic Impact of each Event Type summed across the years')
qplot(EVTYPE_mod,MEANDMG,data=economic_impact)+theme(axis.text.x = element_text(angle = 90, hjust = 1))+labs(x='Event Type',y='Mean Economic Impact',title='Mean Economic Impact per instance of each Event Type')

#rm(list=setdiff(ls(), c("stormdata","event_types","match_table","notmatched","health_impact","economic_impact")))



