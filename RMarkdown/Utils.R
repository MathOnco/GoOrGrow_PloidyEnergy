runAndPlotPDE <- function(this, growBeyondR0 = 0, dorun=T, doplot=T, domov = F, OUTFILE="~/Downloads/ParamSweep.txt", normalize=F, cex=3, cex.lab=1.9, cex.axis=1.6){
  library(animation)
  syso = NULL
  if(dorun){
    this$outFile = "out"
    this$statisticsFile = OUTFILE	
    this$dim = 2
    this$growOutsideR0 = tolower(as.character(as.logical(growBeyondR0)))
    this$saveFiles = "true"
    this$randomCellMotion = 1e-2
    write.table(t(this), file="./goOrGrowTest.i", sep ="\t", quote=F, col.names = F)
    
    # cmd=paste("./go_or_grow",this["tFinal"], this["R"], this["R0"], this["dt"], this["dr"], this["u0"], this["v0"], this["eta"], this["a"], this["xi_u"], this["phi_u"], this["xi_v"], this["phi_v"], "t.txt r.txt E.txt u.txt v.txt", OUTFILE, "1 1", growBeyondR0)
    cmd=paste("./go_or_grow -v 1 -i goOrGrowTest.i")
    print(cmd)
    syso = system(cmd, intern = T)
  }

  t = read.table("out_t.dat", sep="", check.names = F, stringsAsFactors = F)
  r = read.table("out_x.dat", sep="", check.names = F, stringsAsFactors = F)
  u = read.table("out_u.dat", sep="", check.names = F, stringsAsFactors = F)
  v = read.table("out_v.dat", sep="", check.names = F, stringsAsFactors = F)
  # e = read.table("e.dat", sep="", check.names = F, stringsAsFactors = F)
  # ## Number of cells at timepoint of interest (TOI)
  # initCells = sum(u[,1])+sum(v[,1])
  # cells_at_toi = initCells*2^(hours/hcc1954_DT)
  # i_stop=which(apply(v,2,sum) + apply(u,2,sum)>cells_at_toi)[1]
  # if (is.na(i_stop)){ i_stop = ncol(u) }
  i_stop = which(t$V1 > as.numeric(this["tau"]))[1]-1
  if (is.na(i_stop)){ i_stop = nrow(t) }
  
  ##From unitless R to tilde R (measured in um)
  r$V1 = 175*r$V1/as.numeric(this["R0"])
  
  if(!growBeyondR0){
    ## Exclude everything beyond 175 as there are no cells
    ii = which(r$V1 <= 175)
    r=r[ii,,drop=F]
    u=u[ii,,drop=F]
    v=v[ii,,drop=F]
  }
  
  ymax = max(c(0.1,u[,i_stop]))*1.25
  if(normalize){
    u = sweep(u, 2, apply(u, 2, max), FUN="/")
    v = sweep(v, 2, apply(v, 2, max), FUN="/")
    ymax = 1
  }
  ymin = min(c(0.1,u[,i_stop]))
  if(this["v0"]>0){
    ymin = min(v[,i_stop], ymin) 
    ymax = max(v[,i_stop], ymax)*1.25 
  }
  
  if(doplot || domov){
    this_ = unlist(nondimensional2dimensional(eta=this["eta"], a = this["a"], R = this["R"], tau = this["tau"]))
    # Visualize
    vizfun <- function(idx){
      for (i_stop in idx) {
        plot(r$V1, u[,i_stop], pch=20, col="purple", cex.main=0.65, main=paste0(rownames(this),": delta=",round(this_["delta"],1),"; Chi=",round(this_["Chi"]),"; gamma=",round(this_["gamma"],1),"; phi=",round(this["phi_u"],3),"; xi=",round(this["xi_u"],5)), 
             xlab="", ylab="", cex=cex, cex.axis=cex.axis ,ylim=c(ymin,ymax)); #ylim=c(min(c(v[,i_stop],u[,i_stop])), max(c(v[,i_stop]+u[,i_stop],0.3)))
        title(ylab="Cell concentration", xlab=~ "Distance to center " (mu * m), mgp=c(3,3,0), cex.lab=cex.lab)
        if(this["v0"]>0){
          # points(r$V1, (u[,i_stop]+v[,i_stop]), pch=20, col="black")
          points(r$V1, v[,i_stop], pch=20, col="orange", cex=3)
          # points(r$V1, e[,i_stop], pch=20, col="green")
          legend("topleft",c("Grower (v)", "Goer (u)"), fill=c("orange","purple"), cex=max(1, 0.5*cex), bty="n")
        }
      }
    }
    if(domov){
      idx = seq(1, i_stop, by = ceil(i_stop/60));
      # ymax = max(c(0.5,unlist(u)))
      saveGIF({
        vizfun(idx)}, movie.name = paste0("~/Downloads/",fileparts(OUTFILE)$name,"_",rownames(this),".gif"))
    }else{
      vizfun(i_stop)
    }
  }
  return(list(r=r, t=t, u=u, v=v, uv=u+v, i_stop = i_stop, syso = syso))
}

getGreekLetter <- function(v){
  v = strsplit(v,"_")[[1]][1]
  if(v=="eta"){
    return(expression(eta))
  }else  if(v=="xi"){
    return(expression(xi))
  }else  if(v=="k"){
    return(expression(phi))
  }else{
    return(v)
  }
}

hist4ABC <- function(poi, allFits, choices=NULL) {
  o = rep(NA, length(poi))
  names(o) = poi
  for(x in poi){
    uv = unique(allFits[,x])
    if(length(uv)==1){
      o[x] = uv;
      next;
    }
    ho = hist(allFits[,x], 30, xlab=PARAM_ANNO[[x]], col="cyan", border = "white", main="", ylab="# ABC-inferred simulations")
    ii = which.max(ho$counts)
    if(!is.null(choices)){
      ii = which.min(abs(ho$breaks - choices[[x]]))
    }
    arrows(ho$mids[ii], ho$counts[ii], y1 = max(ho$counts)*0.9, col="magenta" )
    o[x] = ho$mids[ii]
  }
  return(o)
}

readAndPrepABCinput <- function(sumStatsFile = "SummaryStats.txt", paramFile = "SimStats.txt", measuredStatsFile = NULL, soi, poi, normalize=T, includeFailed = F){
  # Read simulation output summaries
  allStats = read.table(sumStatsFile, sep="\t", check.names = F, stringsAsFactors = F, row.names = NULL);
  dm = read.table(paramFile, sep="\t", check.names = F, stringsAsFactors = F);
  rownames(dm) = dm$Sim
  param = dm[, !colnames(dm) %in% c("OUTFILE", "Sim")]
  
  ## Prep. sumstats
  sumstats = allStats[grep("u+v",allStats$row.names, fixed = T), ]
  if(!includeFailed){
    sumstats = sumstats[!param$FAIL,]
    param    =    param[!param$FAIL,]
  }
  rownames(sumstats) <- sumstats$simulation
  sumstats = sumstats[,colnames(sumstats) %in% soi, drop=F]
  # param = param[,apply(param, 2, var)>0] 
  
  # ## Transform to dimensional distance to center
  # if(all(!is.na(as.numeric(soi)))){
  #   d2c = as.numeric(colnames(sumstats))
  #   d2c = d2c/max(d2c)
  #   d2c = d2c * (UM_SPOT_DIAM/2)
  #   colnames(sumstats) = as.character(d2c)
  # }
  
  ## COI
  if(!is.null(poi)){
    param = param[,poi];
  }
  if(normalize){
    sumstats = sweep(sumstats, 2, apply(abs(sumstats), 2, max), FUN="/")
  }
  out = list(sumstats = sumstats, param = param)
  
  if(!is.null(measuredStatsFile)){
    ## Read MEP-LINCS summaries and prep for input to ABC
    mepStats = read.table(measuredStatsFile, sep="\t", check.names = F, stringsAsFactors = F)
    mepStats = mepStats[,colnames(mepStats) %in% colnames(sumstats), drop=F]
    # ii = as.character(0:10)
    # # loess.smoot
    # s_mepStats = apply(mepStats[,ii], 1, function(x) loess.smooth(as.numeric(ii), as.numeric(x)))
    # mepStats[,ii] = t(sapply(s_mepStats, function(xy) grpstats(as.matrix(xy$y), round(xy$x), "sum")$sum))
    # s_mepStats = t(sapply(s_mepStats,"[[",2))
    # # # 
    if(normalize){
      mepStats = sweep(mepStats, 2, apply(abs(mepStats), 2, max), FUN="/")
    }
    # # s_mepStats = sweep(s_mepStats, 1, apply(s_mepStats, 1, sum), FUN="/")
    # ##Standardize by a robust estimate of the standard deviation (the median absolute deviation):
    # tmp = sapply(list(mepStats[,colnames(sumstats)], sumstats), function(la) scale(la, center = FALSE, scale = apply(la, 2, mad, na.rm = TRUE)));
    # names(tmp) = c("mepstats", "sumstats")
    # mepStats[,colnames(sumstats)]=tmp$mepstats
    # sumstats=tmp$sumstats
    
    out$measuredStats = mepStats
  }
  return(out)
}


ggplotHist <- function(x, param_2, logx=F, col){
  o = easyGgplot2::ggplot2.histogram(data=param_2, xName = x, groupColors = col, xlim=c(min(param[,x]), min(Inf,max(param[,x]))), groupName="ECM", legendPosition="top", bins=20, alpha=0.5, addDensity=T,
                                     addMeanLine=F, scale = "frequency", meanLineColor="white", meanLineSize=1.5, main = "test") + labs(x=PARAM_ANNO[[x]]) +
    theme(axis.text.x = element_text(face="plain", size = 16),axis.text.y = element_text(face="plain", size = 16),axis.title.x = element_text(face="plain", size = 19.5),axis.title.y = element_text(face="plain", size = 19.5)) +
    theme_bw() + ggtitle(paste(unique(param_2$group), collapse= ";"))
  if(logx){
    o = o+ scale_x_log10()
  }
  return(o)
}


PARAM_ANNOTATION <- function(numberOfSubpopulations = 2){
  param= list(tFinal="Max final time to run the simulation to",
              R="Radius of the dish",
              R_tilde="Radius of the dish",
              R0="Initial seeding radius",
              tau = "Simulation time",
              L = "Characteristic length",
              lambda = "Maximal growth rate",
              dt="Time step size",
              dr="Spatial step size",
              dFUNeta="Diffusion coefficient",
              a="Energy consumption rate",
              Chi = "Chemotactic coefficient",
              delta = "Energy consumption rate",
              gamma = 'Energy diffusion',
              xi_u="Goer's chemotactic inability (xi_u)", # (high = low movement)",  ## Xi_u < Xi_v <=> u begins searching for greener pastures sooner
              xi_v="Grower's chemotactic inability (xi_v)", # (high = low movement)",
              phi_u="Goer's sensitivty to low energy (phi_u)", # (high = sensitive to levels)", ## phi_u > phi_v <=> u seize cell division sooner 
              phi_v="Grower's sensitivty to low energy (phi_v)", # (high = sensitive to levels)",
              u0 = 'u0',
              v0 = 'v0',
              eta = 'Energy diffusion',
              u_max = 'degree of infiltration goers (u)',
              v_max = 'degree of infiltration growers (v)',
              u_tot = 'total # of cells goers (u)',
              v_tot = 'total # of cells growers (v)',
              u_var = "u_var",
              simpson = "simpson",
              tradeoffImbalance = "tradeoffImbalance (goer u - grower v)")
  param[["T"]] = "Simulation time";
  param[["v_rho_u_rho_ratio"]] = "infiltration: growers / goers";
  param[["v_tot_u_tot_ratio"]] = "growers / goers";
  param[["ku_xi_u_ratio"]] = "goer: sensitivity / tendency to stay";
  
  if(numberOfSubpopulations==1){
    param$phi_u = "Sensitivty to low energy"
    param$xi_u = "Chemotactic inability"
    param$u0 = "Confluence at seeding"
  }
  return(param)
}



nondimensional2dimensional <- function(eta, a, R, R_um = 0.5*10*UM_SPOT_DIAM, tau = 1.692957591, anno=F){
  ## characteristic length  
  L = R_um / R
  
  ## Cell growth rate
  lambda = 2 /(2^(hcc1954_DT/24))
  
  days = tau/lambda
  
  ## chemotactic coefficient
  Chi = lambda * L^2
  
  ## Energy diffusion rate
  gamma = Chi * eta
  
  ## Energy consumption rate
  delta = lambda * a;
  
  x = c(lambda, L, days, Chi, gamma, delta)
  names(x) = c("lambda", "L", "T", "Chi", "gamma", "delta")
  if(anno){
    value = c(NA, eta, a, R, tau, NA, NA, R_um, x)
    names(value)[1:8] = c("u0","eta","a", "R", "tau","phi_u", "xi_u", "R_tilde")
    units = c("confluence", rep(NA, 6), "um", "day^-1",NA,"days","um^2 day^-1","um^2 day^-1", "cells^-1 day^-1")
    conversion = c( rep(NA, 8), NA, "R_tilde/R", "tau/lambda", "lambda * L^2", "Chi * eta", "lambda * a")
    return(rbind(value,units, conversion))
  }
  return(x)
}


dimensional2nondimensional <- function(Time, Chi, gamma, delta, R_um, lambda = 2 /(2^(hcc1954_DT/24))){
  a= delta/lambda
  eta = gamma/Chi
  tau = Time*lambda  
  L = sqrt(Chi/lambda)
  R = R_um / L
  
  x = c(eta, a, tau, R);
  names(x) = c("eta", "a", "tau", "R")
  return(x)
}



