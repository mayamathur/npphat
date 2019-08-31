
# analyzes one meta-analytic dataset and adds a row to the results dataframe
# pub.bias can be TRUE or FALSE
# if TRUE, fits the weightr model

analyze_one_meta = function( dat,
                             meta.name,
                             pub.bias,
                             #tail = "above",
                             digits = 2,
                             method) {
  
  # # test only
  # dat = dm
  # meta.name = "Miguel"
  # pub.bias = FALSE
  # digits = 2
  
  ##### Plain Meta-Analysis #####
  if ( pub.bias == FALSE ) {
    meta = rma.uni( yi = dat$yi,
                    vi = dat$vyi,
                    knha = TRUE,
                    method = "REML" )
    
    est = meta$b
    est.se = meta$se
    
    t2 = meta$tau2
    t2.se = meta$se.tau2
    tau.CI = tau_CI(meta)
  }
  
  # ##### With Publication Bias Correction #####
  # if ( pub.bias == TRUE ) {
  #   
  #   # model might not fit
  #   meta = tryCatch( {
  #     
  #     # one-tailed selection model
  #     meta = weightfunct( effect = dat$yi,
  #                         v = dat$vyi,
  #                         steps = c(0.025, 1) )
  #     
  #     # just for fun
  #     cat( "\nEta for ", meta.name, ": ", 1/meta[[2]]$par[3] )
  #     meta
  #     
  #   }, error = function(err) {
  #     # needs to be superassignment because inside fn
  #     NULL
  #   } )
  #   
  #   if ( is.null(meta) ) {
  #     new.row <<- data.frame( Meta = meta.name,
  #                             k = nrow(dat),
  #                             Est = "N/A",
  #                             Tau = "N/A"
  #     )
  #     
  #     Phat.names <<- paste( "Percent above ", unlist(ql), sep = "" )
  #     new.row[, Phat.names ] <<- "N/A"
  #     
  #     
  #     # this should be a global variable
  #     if ( !exists("resE") ) resE <<- new.row
  #     else resE <<- rbind(resE, new.row)
  #     return(resE)
  #   }
  #   
  #   # otherwise we can proceed assuming the meta object exists
  #   H = meta[[2]]$hessian
  #   ses = sqrt( diag( solve(H) ) )
  #   
  #   est = meta[[2]]$par[2]
  #   est.se = ses[2]
  #   
  #   t2 = meta[[2]]$par[1]
  #   t2.se = ses[1]
  #   tau.CI = tau_CI_weightr(meta)
  #   
  # }  # end pub.bias = TRUE loop
  
  
  ##### Parametric (Possibly with Bootstrap) #####
  
  if (method == "parametric") {
    Phat.l = lapply( ql,
                     FUN = function(q) {
                       
                       # set tail based on sign of q
                       if (q >= 0) tail = "above"
                       else tail = "below"
                       
                       MetaUtility::prop_stronger( q = q, 
                                                   M = est,
                                                   t2 = t2, 
                                                   se.M = est.se,
                                                   se.t2 = t2.se,
                                                   tail = tail,
                                                   dat = dat,
                                                   R = boot.reps,
                                                   
                                                   yi.name = "yi",
                                                   vi.name = "vyi" ) } )
  }
  
  if (method == "np.sign") {
    Phat.l = lapply( ql,
                     FUN = function(q) {
                       
                       # set tail based on sign of q
                       if (q >= 0) tail = "above"
                       else tail = "below"
                       
                       # bm
                       prop_stronger_np(q = q,
                                        yi = dat$yi,
                                        vi = dat$vyi,
                                        CI.level = 0.95, 
                                        tail = tail,
                                        R = 2000,
                                        return.vectors = FALSE ) } )
  }
  
  if (method == "calibrated") {
  
    Phat.l = lapply( ql,
                     FUN = function(q) {
                       
                       ens = my_ens( yi = dat$yi, 
                                     sei = sqrt(dat$vyi) )
                       
                       # set tail based on sign of q
                       if (q >= 0) tail = "above"
                       else tail = "below"
                       if ( tail == "above" ) Phat.NP.ens = sum(ens > c(q)) / length(ens)
                       if ( tail == "below" ) Phat.NP.ens = sum(ens < c(q)) / length(ens)
                       
                       Note = NA
                       tryCatch({
                         boot.res.ens = boot( data = dat, 
                                              parallel = "multicore",
                                              R = boot.reps, 
                                              statistic = function(original, indices) {
                                                
                                                b = original[indices,]
                                                
                                                ens.b = my_ens( yi = b$yi, 
                                                                sei = sqrt(b$vyi) )
                                                if ( tail == "above" ) return( sum(ens.b > c(q)) / length(ens.b) )
                                                if ( tail == "below" ) return( sum(ens.b < c(q)) / length(ens.b) )
                                              }
                         )
                         
                         bootCIs.ens = boot.ci(boot.res.ens, type="bca")
                         boot.lo.ens = bootCIs.ens$bca[4]
                         boot.hi.ens = bootCIs.ens$bca[5]
                         
                       }, error = function(err){
                         boot.lo.ens <<- NA
                         boot.hi.ens <<- NA
                         Note <<- err$message
                       } )  # end tryCatch
                       
                       return( data.frame( Est = Phat.NP.ens,
                                           lo = boot.lo.ens,
                                           hi = boot.hi.ens,
                                           boot.note = Note ) )
                  } )  # end lapply

  } # end calibrated method
  
  
  Phat.df = do.call( rbind, 
                     Phat.l )
  Phat.df$string = paste( format_stat( 100*Phat.df$Est,
                                       digits = 0 ),
                          format_CI( 100*Phat.df$lo, 
                                     100*Phat.df$hi,
                                     digits = 0 ),
                          sep = " " )
  

  ##### Put Results in Dataframe #####
  est.string = paste( format_stat( est, digits ),
                      format_CI( est - qnorm(.975) * est.se, 
                                 est + qnorm(.975) * est.se,
                                 digits),
                      sep = " " )
  
  tau.string = paste( format_stat( sqrt(t2), digits),
                      format_CI( tau.CI[1], 
                                 tau.CI[2],
                                 digits ),
                      sep = " " )
  
  new.row = data.frame( Meta = meta.name,
                        k = nrow(dat),
                        Est = est.string,
                        Tau = tau.string,
                        Method = method
  )
  
  # tail is now just for string purposes
  tail = rep("above", length(unlist(ql)))
  tail[unlist(ql) < 0] = "below"
  Phat.names = paste( "Percent ", tail, " ", unlist(ql), sep = "" )
  new.row[, Phat.names ] = NA
  
  new.row[ , Phat.names ] = Phat.df$string
  
  # this should be a global variable
  if ( !exists("resE") ){
    resE <<- new.row
  } else {
    resE <<- rbind(resE, new.row)
  }
} 
