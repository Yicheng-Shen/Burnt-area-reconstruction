# is.nan #### 
# aim: replace all nan with na
is.nan.data.frame <- function(x)
  do.call(cbind, lapply(x, is.nan))

# %notin% ####
`%notin%` = function(x,y) !(x %in% y)

# Pollen bin function ####
# aim: bin pollen data
pollen_bin_function<- function(p_data){
  name<- unique(p_data$Entity.name) 
  df_totalbin<- data.frame()
  for (i in 1:length(name)){
    # i=1
    target_p<- p_data %>% dplyr::filter(Entity.name==name[i])
    minage<- min(target_p$median,na.rm = T)
    maxage<- max(target_p$median,na.rm = T)
    x<- maxage %/% 100 
    y<- maxage %% 100 
    if (y==0){
      maxend<- x
    }else{
      maxend<- x+1
    }
    if (minage<= 100){
      age_interval_head<--100 # the first bin (modern bin) is -100~100 (-50~100, almost)
      age_interval_body<- seq(100, 100*maxend, 100)
      age_interval<- c(age_interval_head,age_interval_body)
    }else{
      minstart<- minage %/% 100
      age_interval<- seq(100*minstart, 100*maxend, 100)
    }
    df_bin<- data.frame()
    df_bin[1:(length(age_interval)-1),1]<- age_interval[2:length(age_interval)]
    colnames(df_bin)[1]<- 'age'
    df_bin[1:(length(age_interval)-1), 2:8]<-target_p[1, 1:7] #add basic information
    # j<-3
    for (j in 1:(length(age_interval)-1)){
      target_bin<- target_p %>% filter(median > age_interval[j] & median <= age_interval[j+1]) 
      databin<- target_bin[, which(colnames(target_bin) == "Abies"):length(target_bin)] %>% summarise_each(sum)
      databin<- databin/rowSums(databin,na.rm = T)
      len_df_bin<-8+length(databin) 
      df_bin[j, 9:len_df_bin]<- databin
    }
    #run the is.nan function fist
    df_bin[, 9:len_df_bin][is.nan(df_bin[, 9:len_df_bin])]<-NA
    df_totalbin<- rbind(df_totalbin, df_bin)
  }
  # Remove NA ages
  taxa_start<- which(colnames(df_totalbin)=="reference")+1
  taxa_end<- length(colnames(df_totalbin))
  df_totalbin<- df_totalbin[rowSums(is.na(df_totalbin[ , taxa_start:taxa_end])) != ncol(df_totalbin[ , taxa_start:taxa_end]), ]
  df_totalbin[,taxa_start:taxa_end] <- sapply(df_totalbin[,taxa_start:taxa_end],as.numeric)
  return(df_totalbin)
  # write.csv(p_data, paste0(deparse(substitute(p_data)), "_bin.csv"), row.names = F, fileEncoding = "UTF-8")
}

# Charcoal influx transformation function ####
charcoal_fluxtrans_function<- function(data_query){
  dir.create(file.path(path, '/charcoaltrans'), recursive = TRUE)
  # query label and path to query and query name
  csvpath <- paste0(path,"/charcoaltrans/")
  
  # debug/log file
  debugname <- "RPD_influx_debug.txt"
  # open the debug/log file
  debugfile <- file(paste(csvpath, debugname, sep=""), "w")
  
  data_query <-  c
  # remove negative quantity 
  data_query<- data_query %>% dplyr::filter(quantity>=0) 
  #head(data_query); tail(data_query)
  
  # The main part of the script loops over the individual entities, and does various checks and calculations, including
  # calculation of sedimentation rates and deposition times
  # checking for age or depth reversals
  # setting of various indicator variables and flags
  # calculation of alternative quanties (e.g. influx, given concentrations)
  # further checking for anomalies
  # writing out a .csv file for each site
  
  # number of samples in the file
  nsamp <- length(data_query$ID_ENTITY); nsamp
  
  ##lets check what 'TYPE' codes we have
  
  unique(data_query$TYPE)
  
  ##FOR NOW (17/03/2020): will fill in the TYPE = NA with a character to avoid problems with loop. we foresee that there will be no NAs in this variable in future. so this is a temporary solution. If your data has no NAs under the TYPE variable in entity metadata, dont worry about this. The RPD data has 'NULL' instead of NA which is a character anyway, so if you have NULL for TYPE in your data, sub this in in the loop for 'undefined':
  
  data_query$TYPE[is.na(data_query$TYPE)] <- "undefined" #step not necessary if no NAs in your TYPE variable
  unique(data_query$TYPE)
  
  
  ls1 <- ls() #to check which entities are successfully saved with transformation.
  
  # get names to check if renaming is necessary
  names(data_query)
  
  maxsites <- unique(data_query$ID_ENTITY) 
  class(data_query$ID_ENTITY)
  
  for (j in maxsites) { ##trying looping using entity IDs: changed from j in 1:maxsites
    #j=1701
    nsamp <- 0
    sitedata <- data_query[data_query$ID_ENTITY == j, ]
    nsamp <- length(sitedata$ID_ENTITY)
    
    # Process the data for the j-th site
    # Define some local variables, and replace NA's with a missing values code
    
    if (nsamp > 1) { ##I changed this from 0 to 1, because we have entities that are only one row in v4.
      jchar <- as.character(j)
      nsampchar <- as.character(nsamp)
      writeLines(paste("Site",jchar,nsampchar,"samples", sep=" "), con = debugfile, sep = "\n")
      
      # local variables
      depth <- sitedata$depth; age <- sitedata$median; quant <- sitedata$quantity; xst_level <- sitedata$TYPE #added x1st level here to make consistent with his example he sent me for cygnet (mycsv-to-csv_cygnet...r)
      
      # recode NA's to missing
      miss <- -9999.0
      depth[is.na(depth)] <- miss
      age[is.na(age)] <- miss
      quant[is.na(quant)] <- miss
      
      # define some new variables
      thickness <- rep(miss, nsamp); dep_time <- rep(miss, nsamp); sed_rate <- rep(miss, nsamp)
      unit_dep_time <- rep(miss, nsamp) #removed x1st level here because added it as vector above
      
      
      # sed rate and deposition time
      # first (top) sample
      if (depth[1] != miss && depth[2] != miss) {
        thickness[1] <- (depth[2] - depth[1])  
        dep_time[1] <- age[2] - age[1]
        if (dep_time[1] > 0.0) sed_rate[1] <- thickness[1]/dep_time[1]
        if (sed_rate[1] != miss) unit_dep_time[1] <- 1.0/sed_rate[1] #the 1/sed rate - 1 is a placeholder because it will be *quant
      }
      # samples 2 to nsamp-1
      for (i in 2:(nsamp-1)) {
        if (depth[1] != miss && depth[2] != miss) {
          thickness[i] <- (depth[i+1] - depth[i])  
          dep_time[i] <- ((age[i+1] + age[i])/2.0) - ((age[i] + age[i-1])/2.0)
          if (dep_time[i] > 0.0) sed_rate[i] <- thickness[i]/dep_time[i]
          if (sed_rate[i] != miss) unit_dep_time[i] <- 1.0/sed_rate[i] 
        }
      }
      # last (bottom) sample
      if (depth[nsamp-1] != miss  && depth[nsamp] != miss) {
        thickness[nsamp] <- thickness[nsamp-1] # replicate thickness
        dep_time[nsamp] <- age[nsamp] - age[nsamp-1]
        sed_rate[nsamp] <- sed_rate[nsamp-1] # replicate sed_rate
        unit_dep_time[nsamp] <- unit_dep_time[nsamp-1]
      }
      
      # counts of missing values
      depth_count <- 0; age_count <- 0; quant_count <- 0; sed_rate_count <- 0; sed_rate_flag <- 1
      depth_count <- sum(depth != miss)
      age_count <- sum(age != miss)
      quant_count <- sum(quant != miss)
      sed_rate_count <- sum(sed_rate != miss)
      if (sed_rate_count != nsamp) sed_rateflag = 0
      
      # check for age or depth reversals, and zero or negative sed rates (in nonmissing data)
      depth_reversal <- 0; age_reversal <- 0; sed_rate_zeroneg <- 0         
      for (i in 2:nsamp) {
        if (age[i] != miss && age[i-1] != miss && age[i] <= age[i-1]) age_reversal=1
        if (depth[i] != miss && depth[i-1] != miss) {
          if (depth[i] <= depth[i-1]) depth_reversal=1
        } 
      }
      for (i in 2:nsamp) {
        if (sed_rate[i] != miss && sed_rate[i] <= 0.0) sed_rate_zeroneg=1
      }
      
      # set and write out various flags
      if (depth_count != 0 && depth_count != nsamp) {
        writeLines(paste("**** has a missing depth when some are nonmissing", sep=" "), con = debugfile, sep = "\n")
      }
      if (age_count != 0 && age_count != nsamp) {
        writeLines(paste("**** has a missing age when some are nonmissing", sep=" "), con = debugfile, sep = "\n")
      }
      if (quant_count != 0 && quant_count != nsamp) {
        writeLines(paste("**** has a missing quantity when some are nonmissing", sep=" "), con = debugfile, sep = "\n")
      }
      if (sed_rate_count != 0 && sed_rate_count != nsamp) {
        writeLines(paste("**** has a missing sed rate when some are nonmissing", sep=" "), con = debugfile, sep = "\n")
      }
      if (depth_reversal != 0) {
        writeLines(paste("**** has a depth reversal", sep=" "), con = debugfile, sep = "\n")
      }
      if (age_reversal != 0) {
        writeLines(paste("**** has an age reversal", sep=" "), con = debugfile, sep = "\n")
      }
      if (sed_rate_zeroneg != 0) {
        writeLines(paste("**** has zero or negative sed rates", sep=" "), con = debugfile, sep = "\n")
      }
      
      # alternative quantities
      
      conc <- rep(miss, nsamp); influx <- rep(miss, nsamp)
      influx_source <- rep("none", nsamp) ; conc_source <- rep("none", nsamp)
      
      # select case based on xst_level
      if (xst_level[1] == "INFL")          # adopt influx values as they are, calculate concentration
      {  
        influx <- quant
        influx_source <- "data"
        if (influx != miss && unit_dep_time != miss && sed_rate != 0.0) {
          conc <- influx * unit_dep_time
          conc_source <- "calculated from influx "
        } else {
          conc <- quant
          conc_source <- "copied from quant "
        }
        writeLines("INFL", con = debugfile, sep = "\n")
      }
      
      else if (xst_level[1] == "CONC")     # calculate influx, adopt conc values as they are
      {
        conc <- quant
        conc_source <- "data"
        if (conc != miss && sed_rate != miss && sed_rate != 0.0) {
          influx <- quant * sed_rate
          influx_source <- "calculated from conc "
        } else {
          influx <- quant
          influx_source <- "copied from quant "
        }  
        writeLines("CONC", con = debugfile, sep = "\n")
      }
      
      else if (xst_level[1] == "C0P0")     # assume quantity is concentration like
      {
        conc <- quant
        conc_source <- "C0P0"
        if (sed_rate != miss && sed_rate != 0.0) {
          influx <- quant * sed_rate
          influx_source <- "calculated from C0P0 (conc) "
        } else {
          influx <- quant
          influx_source <- "copied from quant "
        }    
        writeLines("C0P0", con = debugfile, sep = "\n")
      }
      
      
      
      else if (xst_level[1] == "per unit weight")     # for now we are assuming that this is like concentration, but meeting on 19/03 to discuss this. (I changed this from "SOIL") - yes its concentration like according to sandy.
      {
        conc <- quant
        conc_source <- "per unit weight"
        if (sed_rate != miss && sed_rate != 0.0) {
          influx <- quant * sed_rate
          influx_source <- "calculated from per unit weight "
        } else {
          influx <- quant
          influx_source <- "copied from quant "
        }
        writeLines("per unit weight", con = debugfile, sep = "\n")
      }
      
      else if (xst_level[1] == "undefined")     # There are some NAs in the type variable as of 17/03/2020, which I renamed above to "undefined". assuming these are concentration like but ideally one would make sure first and exclude if not (exclude by deleting this block of code)
      {
        conc <- quant
        conc_source <- "data"
        if (conc != miss && sed_rate != miss && sed_rate != 0.0) {
          influx <- quant * sed_rate
          influx_source <- "calculated from conc "
        } else {
          influx <- quant
          influx_source <- "copied from quant "
        }  
        writeLines("undefined", con = debugfile, sep = "\n")
      }
      
      else if (xst_level[1] == "NULL")     # you may have some NULL under type. probably not. assume these are concentration like but ideally one would make sure and exclude if not exclude by deleting this block of code)
      {
        conc <- quant
        conc_source <- "copied from quant "
        influx <- quant
        influx_source <- "copied from quant "
        writeLines("NULL", con = debugfile, sep = "\n")
      }
      else 
      {
        conc <- quant
        conc_source <- "copied from quant "
        influx <- quant
        influx_source <- "copied from quant "
        writeLines("Unknown", con = debugfile, sep = "\n") 
      }
    }
    
    #this script works for the northern extratropics, however may have to add in some other types here for other regions e.g. SOIL and OTHE - doing a unique(entity table$TYPE) will tell you what you need here.Would need to make sure they are concentration like if including them. AREA is another one in the database. Also, there is a field in the database (COUNT), which is deliberately excluded from this calculation because it is not concentration like.
    
    # check for influx == 0.0 everywhere
    nzero <- 0
    nzero <- sum(influx != 0.0)
    if (nzero == 0) {
      writeLines(paste("**** has no non-zero influx values", sep=" "), con = debugfile, sep = "\n")
    }
    
    # .csv out
    if (nsamp > 1 && nzero > 0) { ##Quality control bit after z score calculation
      
      # get siteid string
      entname <- unique(sitedata$entity_name)
      entname <- substr(entname, start = 1, stop = 3)
      entidchar <- as.character(j)
      if (j >= 1) entid <- paste("000", entidchar, entname, sep="") ##this needs to be fixed to be suitable for entity names.
      if (j >= 10) entid <- paste("00", entidchar, entname, sep="")
      if (j >= 100) entid <- paste("0", entidchar, entname, sep="")
      if (j >= 1000) entid <- paste(    entidchar, entname, sep="")
      
      
      # assemble output data and write it out
      # also save a list of entity IDs that get successfully written out, for quality control.
      ls1 <- append(ls1,entidchar,)
      
      siteidvec <- sitedata$ID_SITE
      entidvec <- sitedata$ID_ENTITY
      sampidvec <- sitedata$ID_SAMPLE
      
      outdata <- data.frame(siteidvec, entidvec, sampidvec, depth, age, sed_rate, quant, conc, influx, xst_level, conc_source, influx_source) #removed ID_SAMPLE from this block of code - can add back in if you want
      names(outdata) <- c("ID_SITE", "ID_ENTITY", "ID_SAMPLE", "depth", "est_age", "sed_rate", "quant", "conc",
                          "influx", "xst_level", "conc_source", "influx_source" )
      csvfile <- paste(path,"/charcoaltrans/",entid,"_RPD_influx_17_03_2020.csv", sep="") # add your own desired path
      write.csv(outdata, csvfile, row.names=FALSE)
    } 
    print(j)
  }
  
  close(debugfile)
  
  length(unique(data_query$ID_ENTITY)) 
  
  ##try read all of them in and rbind them into a single dataframe and then find the ID_entities that werent conserved. 
  
  
  checkdat <- outdata[0,]
  for(z in maxsites){
    #compose filename
    RPDent<- data_query
    tmp <- RPDent[RPDent$ID_ENTITY == z,]
    entname <- unique(tmp$entity_name) #make sure this bit is working
    entname <- substr(entname, start = 1, stop = 3)
    entidchar <- as.character(z)
    if (z >= 1) entid <- paste("000", entidchar, entname, sep="") ##this needs to be fixed to be suitable for entity names.
    if (z >= 10) entid <- paste("00", entidchar, entname, sep="")
    if (z >= 100) entid <- paste("0", entidchar, entname, sep="")
    if (z >= 1000) entid <- paste(    entidchar, entname, sep="")
    
    csvin <- paste(path,"/charcoaltrans/",entid,"_RPD_influx_17_03_2020.csv", sep="") # the individual files you wrote out earlier
    a <- read.csv(csvin)
    checkdat <- rbind(checkdat, a)
  }
  
  ent_infl <- RPDent[RPDent$ID_ENTITY %in% checkdat$ID_ENTITY,]
  #save this entity metadata file to the same folder that all the individual entity csv files are saved.
  return(checkdat)
  write.csv(ent_infl, paste0(path,"/charcoaltrans/RPD_infl_mtdata_24_03_2020.csv"),row.names = F)
  write.csv(checkdat, paste0(path,"/charcoaltrans/checkdat.csv"),row.names = F)
}

# Charcoal bin function ####
charcoal_bin_function<- function(c_influx){
  name<- unique(c_influx$entity_name)
  df_totalbin<- data.frame()
  for (i in 1:length(name)){
    target_c<- c_influx %>% filter(entity_name==name[i])
    minage<- min(target_c$est_age,na.rm = T)
    maxage<- max(target_c$est_age,na.rm = T)
    x<- maxage %/% 100
    y<- maxage %% 100 
    if (y==0){
      maxend<- x
    }else{
      maxend<- x+1
    }
    if (minage<=100){
      ##### first part made change 
      age_interval_head<--100 # the first bin (modern bin) is -100~100 (-50~100, almost)
      age_interval_body<- seq(100, 100*maxend, 100)
      age_interval<- c(age_interval_head,age_interval_body)
    }else{
      
      minstart<- minage %/% 100
      age_interval<- seq(100*minstart, 100*maxend, 100)
    }
    df_bin<- data.frame()
    df_bin[1:(length(age_interval)-1),1]<- age_interval[2:length(age_interval)]
    colnames(df_bin)[1]<- 'age'
    target_c_info<- target_c %>% dplyr::select(ID_SITE, ID_ENTITY, entity_name, site_name, latitude, longitude, elevation)
    df_bin[1:(length(age_interval)-1), 2:8]<-target_c_info[1, 1:7] #add basic information
    
    for (i in 1:(length(age_interval)-1)){
      target_bin<- target_c %>% filter(est_age>age_interval[i] & est_age <= age_interval[i+1]) ##second part we make change
      databin<- mean(target_bin$influx)
      df_bin[i,9]<- databin
    }
    colnames(df_bin)[9]<- 'influx'
    # run the is.nan.data.frame function first
    df_bin[,9][is.nan(df_bin[,9])]<-NA
    # max transformation
    max<- df_bin$influx/max(df_bin$influx,na.rm = T)
    df_bin$max<- max
    df_totalbin<- rbind(df_totalbin, df_bin)
  }
  which(colnames(df_totalbin)=="influx") #9
  which(colnames(df_totalbin)=="max") #10
  df_totalbin<- df_totalbin[rowSums(is.na(df_totalbin[ , 9:10])) != ncol(df_totalbin[ , 9:10]), ]
  return(df_totalbin)
}


# Make composite plot of reconstructed burnt area####
Composite_plot<- function(data=pre_core, method="pre_wapls", hw=300, nreps=500){
  library(locfit)
  # target ages for fitted values
  targbeg <- 0
  targend <- 15000
  targstep <- 100
  
  # array sizes
  maxrecs <- 2000
  maxreps <- 1000
  
  # read the list of sites
  sites<- unique(data$Entity.name)
  #sites <- read.csv(sitelistfile)
  head(sites)
  #ns <- length(sites[,1]) #length(sites$ID_SITE)
  ns<- length(sites)
  # arrays for data and fitted values
  age <- matrix(NA, ncol=ns, nrow=maxrecs)
  influx <- matrix(NA, ncol=ns, nrow=maxrecs)
  nsamples <- rep(0, maxrecs)
  targage <- seq(targbeg,targend,targstep)
  targage.df <- data.frame(x=targage)
  lowage <- targage - hw; highage <- targage + hw
  ntarg <- length(targage)
  yfit <- matrix(NA, nrow=length(targage.df$x), ncol=maxreps)
  
  # arrays for sample number and effective window span tracking
  ndec <- matrix(0, ncol=ntarg, nrow=ns)
  ndec_tot <- rep(0, ntarg)
  xspan <- rep(0, ntarg)
  ninwin <- matrix(0, ncol=ntarg, nrow=ns)
  ninwin_tot <- rep(0, ntarg)
  
  # read and store the presample (binned) files as matrices of ages and influx values
  # I didn't change the name of influx. Here influx is our burnt area fraction.
  # This part will change when use different data sets.
  id_method<- which(colnames(data)==method)
  ii <- 0
  for (i in 1:length(data$Entity.name)){
    # i<- 1
    indata<- data %>% dplyr::filter(Entity.name == unique(data$Entity.name)[i])
    nsamp <-  length(indata$age)
    if (nsamp > 0) {
      ii <- ii+1
      age[1:nsamp,ii] <- indata$age 
      influx[1:nsamp,ii] <- as.data.frame(indata[,id_method])[,1] # we only want to extract the target column without geometry
      nsamples[ii] <- nsamp
    }
  }
  nsites <- ii
  
  # number of sites with data
  nsites
  
  # trim samples to age range
  influx[age >= targend+hw] <- NA
  age[age >= targend+hw] <- NA
  
  # censor abs(influx) values > 10
  influx[abs(influx) >= 10] <- NA
  age[abs(influx) >= 10] <- NA
  
  # count number of sites that contributed to each fitted value
  ptm <- proc.time()
  for (i in 1:ntarg) {
    agemax <- -1e32; agemin <- 1e32
    for (j in 1:nsites) {
      for (k in 1:nsamples[j]) {
        if (!is.na(age[k,j])) {
          ii <- (age[k,j]-targage[1])/targstep + 1
          #print (c(i,j,k,ii))
          if (ii > 0 && ii <= ntarg) {ndec[j,ii] = 1}
          if (age[k,j] >= targage[i]-hw && age[k,j] <= targage[i]+hw) {
            ninwin[j,i] = 1
            if (agemax < age[k,j]) {agemax <- age[k,j]}
            if (agemin > age[k,j]) {agemin <- age[k,j]}
          }
        }
      }
    }
    ndec_tot[i] <- sum(ndec[,i])
    ninwin_tot[i] <- sum(ninwin[,i])
    xspan[i] <- agemax - agemin
  }
  proc.time() - ptm
  head(cbind(targage,ndec_tot,xspan,ninwin_tot))
  
  ptm <- proc.time()
  # 1. reshape matrices into vectors 
  x <- as.vector(age)
  y <- as.vector(influx)
  lfdata <- data.frame(x,y)
  lfdata <- na.omit(lfdata)
  x <- lfdata$x; y <- lfdata$y
  
  # 2. locfit
  # initial fit, unresampled (i.e. all) data
  loc01 <- locfit(y ~ lp(x, deg=1, h=hw), maxk=800, family="qrgauss")
  summary(loc01)
  
  # 3. get  fitted values
  pred01 <- predict(loc01, newdata=targage.df, se.fit=TRUE)
  loc01_fit <- data.frame(targage.df$x, pred01$fit)
  fitname <- paste("locfit_",as.character(hw), sep="")
  colnames(loc01_fit) <- c("age", fitname)
  head(loc01_fit)
  proc.time() - ptm
  ptm <- proc.time()
  
  # Bootstrap samples
  
  # Step 1 -- Set up to plot individual replications
  pdf(file = paste0(path, "/rm_rare/", method, "_reconstructions.pdf")) 
  plot(x, y, xlab="Age (BP 1950)", ylab=fitname, xlim=c(15000,0), xaxp  = c(15000, 0, 15),ylim=c(0,0.001), type="n") 
  
  # Step 2 -- Do the bootstrap iterations, and plot each composite curve
  set.seed(10) # do this to get the same sequence of random samples for each run
  
  for (i in 1:nreps) {
    print(i)
    randsitenum <- sample(seq(1:nsites), nsites, replace=TRUE)
    # print(head(randsitenum))
    x <- as.vector(age[,randsitenum])
    y <- as.vector(influx[,randsitenum])
    lfdata <- data.frame(x,y)
    lfdata <- na.omit(lfdata)
    x <- lfdata$x; y <- lfdata$y
    locboot <- locfit(y ~ lp(x, deg=1, h=hw), maxk=800, maxit=20, family="qrgauss")
    predboot <- predict(locboot, newdata=targage.df, se.fit=TRUE)
    yfit[,i] <- predboot$fit
    # note plotting lines is slowww
    lines(targage.df$x, yfit[,i], lwd=2, col=rgb(0.5,0.5,0.5,0.10))
    if (i %% 10 == 0) {print(i)}
  }
  
  # Step 3 -- Plot the unresampled (initial) fit
  fitname <- paste("locfit_",as.character(hw), sep="")
  colnames(loc01_fit) <- c("age", fitname)
  lines(loc01_fit[,1], loc01_fit[,2], lwd=2.5, col="yellow")
  
  # Step 4 -- Find and add bootstrap CIs
  yfit95 <- apply(yfit, 1, function(x) quantile(x,prob=0.975, na.rm=T))
  yfit05 <- apply(yfit, 1, function(x) quantile(x,prob=0.025, na.rm=T))
  lines(targage.df$x, yfit95, lwd=2, col="#FFA500")
  lines(targage.df$x, yfit05, lwd=2, col="#FFA500")
  dev.off()
  curveout <- data.frame(cbind(targage.df$x, pred01$fit, yfit95, yfit05, ndec_tot, xspan, ninwin_tot))
  colnames(curveout) <- c("age", "locfit", "cu95", "cl95", "nsites", "window", "ninwin")
  write.table(curveout, paste0(path, "/rm_rare/", method,"_curveout.csv"), col.names=TRUE, row.names=FALSE, sep=",")
}



# Make composite plot for charcoal####
Composite_plot_char<- function(data=mc_data, method="mc_check", hw=300, nreps=500){
  library(locfit)
  
  # target ages for fitted values
  targbeg <- 0
  targend <- 15000
  targstep <- 100
  
  # array sizes
  maxrecs <- 2000
  maxreps <- 1000
  
  # plot output 
  # plotout <-  "pdf" #"screen" 
  
  # read the list of sites
  sites<- unique(data$entity_name)
  #sites <- read.csv(sitelistfile)
  head(sites)
  #ns <- length(sites[,1]) #length(sites$ID_SITE)
  ns<- length(sites)
  # arrays for data and fitted values
  age <- matrix(NA, ncol=ns, nrow=maxrecs)
  influx <- matrix(NA, ncol=ns, nrow=maxrecs)
  nsamples <- rep(0, maxrecs)
  targage <- seq(targbeg,targend,targstep)
  targage.df <- data.frame(x=targage)
  lowage <- targage - hw; highage <- targage + hw
  ntarg <- length(targage)
  yfit <- matrix(NA, nrow=length(targage.df$x), ncol=maxreps)
  
  # arrays for sample number and effective window span tracking
  ndec <- matrix(0, ncol=ntarg, nrow=ns)
  ndec_tot <- rep(0, ntarg)
  xspan <- rep(0, ntarg)
  ninwin <- matrix(0, ncol=ntarg, nrow=ns)
  ninwin_tot <- rep(0, ntarg)
  
  # read and store the presample (binned) files as matrices of ages and influx values
  # I didn't change the name of influx. Here influx is our burnt area fraction.
  # This part will change when use different data sets.
  ii <- 0
  for (i in 1:length(data$entity_name)){
    # i<- 1
    indata<- data %>% dplyr::filter(entity_name == unique(data$entity_name)[i])
    nsamp <-  length(indata$age)
    if (nsamp > 0) {
      ii <- ii+1
      age[1:nsamp,ii] <- indata$age  
      influx[1:nsamp,ii] <- indata$max
      nsamples[ii] <- nsamp
    }
  }
  nsites <- ii
  
  # number of sites with data
  nsites
  
  # trim samples to age range
  influx[age >= targend+hw] <- NA
  age[age >= targend+hw] <- NA
  
  # censor abs(influx) values > 10
  influx[abs(influx) >= 10] <- NA
  age[abs(influx) >= 10] <- NA
  
  # count number of sites that contributed to each fitted value
  ptm <- proc.time()
  for (i in 1:ntarg) {
    agemax <- -1e32; agemin <- 1e32
    for (j in 1:nsites) {
      for (k in 1:nsamples[j]) {
        if (!is.na(age[k,j])) {
          ii <- (age[k,j]-targage[1])/targstep + 1
          #print (c(i,j,k,ii))
          if (ii > 0 && ii <= ntarg) {ndec[j,ii] = 1}
          if (age[k,j] >= targage[i]-hw && age[k,j] <= targage[i]+hw) {
            ninwin[j,i] = 1
            if (agemax < age[k,j]) {agemax <- age[k,j]}
            if (agemin > age[k,j]) {agemin <- age[k,j]}
          }
        }
      }
    }
    ndec_tot[i] <- sum(ndec[,i])
    ninwin_tot[i] <- sum(ninwin[,i])
    xspan[i] <- agemax - agemin
  }
  proc.time() - ptm
  head(cbind(targage,ndec_tot,xspan,ninwin_tot))
  
  ptm <- proc.time()
  # 1. reshape matrices into vectors 
  x <- as.vector(age)
  y <- as.vector(influx)
  lfdata <- data.frame(x,y)
  lfdata <- na.omit(lfdata)
  x <- lfdata$x; y <- lfdata$y
  
  # 2. locfit
  # initial fit, unresampled (i.e. all) data
  loc01 <- locfit(y ~ lp(x, deg=1, h=hw), maxk=800, family="qrgauss")
  summary(loc01)
  
  # 3. get  fitted values
  pred01 <- predict(loc01, newdata=targage.df, se.fit=TRUE)
  loc01_fit <- data.frame(targage.df$x, pred01$fit)
  fitname <- paste("locfit_",as.character(hw), sep="")
  colnames(loc01_fit) <- c("age", fitname)
  head(loc01_fit)
  proc.time() - ptm
  ptm <- proc.time()
  
  # Bootstrap samples
  
  # Step 1 -- Set up to plot individual replications
  #if (plotout == "pdf") {pdf(file=paste(curvecsvpath,pdffile,sep=""))}
  pdf(file = paste0(path, "/rm_rare/", method, ".pdf")) 
  plot(x, y, xlab="Age (BP 1950)", ylab=fitname, xlim=c(15000,0), xaxp  = c(15000, 0, 15),ylim=c(0,1), type="n") 
  
  # Step 2 -- Do the bootstrap iterations, and plot each composite curve
  set.seed(10) # do this to get the same sequence of random samples for each run
  
  for (i in 1:nreps) {
    print(i)
    randsitenum <- sample(seq(1:nsites), nsites, replace=TRUE)
    # print(head(randsitenum))
    x <- as.vector(age[,randsitenum])
    y <- as.vector(influx[,randsitenum])
    lfdata <- data.frame(x,y)
    lfdata <- na.omit(lfdata)
    x <- lfdata$x; y <- lfdata$y
    locboot <- locfit(y ~ lp(x, deg=1, h=hw), maxk=800, maxit=20, family="qrgauss")
    predboot <- predict(locboot, newdata=targage.df, se.fit=TRUE)
    yfit[,i] <- predboot$fit
    # note plotting lines is slowww
    lines(targage.df$x, yfit[,i], lwd=2, col=rgb(0.5,0.5,0.5,0.10))
    if (i %% 10 == 0) {print(i)}
  }
  
  # Step 3 -- Plot the unresampled (initial) fit
  fitname <- paste("locfit_",as.character(hw), sep="")
  colnames(loc01_fit) <- c("age", fitname)
  lines(loc01_fit[,1], loc01_fit[,2], lwd=2, col="yellow")
  
  # Step 4 -- Find and add bootstrap CIs
  yfit95 <- apply(yfit, 1, function(x) quantile(x,prob=0.975, na.rm=T))
  yfit05 <- apply(yfit, 1, function(x) quantile(x,prob=0.025, na.rm=T))
  lines(targage.df$x, yfit95, lwd=1, col="#FFA500")
  lines(targage.df$x, yfit05, lwd=1, col="#FFA500")
  dev.off()
  curveout <- data.frame(cbind(targage.df$x, pred01$fit, yfit95, yfit05, ndec_tot, xspan, ninwin_tot))
  colnames(curveout) <- c("age", "locfit", "cu95", "cl95", "nsites", "window", "ninwin")
  write.table(curveout, paste0(path, "/rm_rare/", method,"_curveout.csv"), col.names=TRUE, row.names=FALSE, sep=",")
}

# Plot composite curve ####
curveout_plot<- function(curve){
  par(mar = c(5, 5, 3, 5))
  plot(curve$age, curve$locfit*1000, type ="l", 
       ylab = "Reconstructed burnt area fraction",
       xlab = "Age (BP 1950)",
       col = "blue",
       xlim=c(15000,0), xaxp  = c(15000, 0, 15),
       ylim=c(1,10))
  polygon(c(rev(curve$age), curve$age), 
          c(rev(curve$cu95*1000), curve$cl95*1000),
          col = adjustcolor('grey',alpha.f=0.5) , border = NA)
  lines(curve$age, curve$locfit*1000, type ="l", lwd=2.5, col="blue")
}

# Plot char curve
curveout_char_plot<- function(curve){
  par(mar = c(5, 5, 3, 5))
  plot(curve$age, curve$locfit, type ="l", 
       ylab = "Max transformed charcoal",
       xlab = "Age (BP 1950)",
       col = "blue",
       xlim=c(15000,0), xaxp  = c(15000, 0, 15),
       ylim=c(0,1))
  polygon(c(rev(curve$age), curve$age), 
          c(rev(curve$cu95), curve$cl95),
          col = adjustcolor('grey',alpha.f=0.5) , border = NA)
  lines(curve$age, curve$locfit, type ="l", lwd=2.5, col="blue")
}
