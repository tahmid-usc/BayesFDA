# Tree volume calculation based on Regional Stem Taper Equations for Eleven Conifer Species 
#in the Acadian Region of North America: Development and Assessment (Li et al. 2021)

tree <- read.csv("E:/R Project/BayesFDA/data/tree/tree.csv")


Smalians <- function(d1, d2, leng) {
  L = (d1 / 2)^ 2 * pi
  S = (d2 / 2)^ 2 * pi
  Vol = ((L + S) / 2) * leng
  return(Vol) 
}


KozakExp02Mod <- function(DHT, HT, DBHO) {
  
  reg1 <- 1 # removed reg1 from argument, set it to 1
  p = 1.37 / HT
  Z = DHT / HT
  Xi = (1 - Z ^ (1 / 3)) / (1 - p ^ (1 / 3))
  Qi = 1 - Z ^ (1 / 3)
  KozakExp02Mod = (a0 * (DBHO^a1) * (HT^a2)) * Xi^(b1 * Z^4 + b2 * (exp(-DBHO / HT)) + 
                                                     b3 * Xi^0.1 + b4 * (1 / DBHO) + b5 * HT^Qi + b6 * Xi + b7 * reg1)
  return(KozakExp02Mod)
}


treeVolume <- function(HT, DBHO, spec){
  
  if(spec == "BF"){
    a0 = 0.88075316
    a1 = 1.01488665
    a2 = 0.01958804
    b1 = 0.41951756
    b2 = -0.67232564
    b3 = 0.54329725
    b4 = 1.48181152
    b5 = 0.06470371
    b6 = -0.34684837
    b7 = 0
  } else if(spec == "BS"){
    a0 = 0.80472902
    a1 = 1.00804553
    a2 = 0.05601099
    b1 = 0.35533529
    b2 = -0.41320046
    b3 = 0.41527304
    b4 = 1.11652424
    b5 = 0.0990167
    b6 = -0.40992056
    b7 = 0.11394943
  } else if(spec == "WS"){
    a0 = 0.75826241
    a1 = 0.98481863
    a2 = 0.09956165
    b1 = 0.36505143
    b2 = -0.51501314
    b3 = 0.55913869
    b4 = 0.75846281
    b5 = 0.07011851
    b6 = -0.44928376
    b7 = 0.07830011
  } else if(spec == "RS"){
    a0 = 0.89797987
    a1 = 1.00579742
    a2 = 0.01667313
    b1 = 0.49500865
    b2 = -0.63375155
    b3 = 0.3836274
    b4 = 1.41380994
    b5 = 0.08866994
    b6 = -0.29753964
    b7 = 0.15192029
  } else if(spec == "JP"){
    a0 = 0.931552701
    a1 = 1.008192708
    a2 = -0.004177373
    b1 = 0.431297353
    b2 = -0.863672736
    b3 = 0.511698303
    b4 = 2.232484834
    b5 = 0.059865263
    b6 = -0.331897255
    b7 = 0.039630786
  } else if(spec == "RP"){
    a0 = 0.9717883
    a1 = 1.00113806
    a2 = -0.01597933
    b1 = 0.51143292
    b2 = -0.9739954
    b3 = 0.25844201
    b4 = 4.75315518
    b5 = 0.05846224
    b6 = -0.12372176
    b7 = 0
  } else if(spec == "WP"){
    a0 = 1.04881379
    a1 = 1.00779696
    a2 = -0.04595353
    b1 = 0.38085445
    b2 = -0.85956463
    b3 = 0.34380669
    b4 = 4.60836993
    b5 = 0.111855
    b6 = -0.5523203
    b7 = 0
  } else if(spec == "NS"){
    a0 = 0.9308817
    a1 = 0.97360573
    a2 = 0.03522864
    b1 = 0.65078104
    b2 = -0.30355787
    b3 = 0.37832812
    b4 = 1.18815216
    b5 = 0.03111631
    b6 = -0.03172809
    b7 = 0
  } else if(spec == "EH"){
    a0 = 0.960235102
    a1 = 1.00821143
    a2 = -0.025167937
    b1 = 0.825260258
    b2 = 1.962520834
    b3 = 0.415234319
    b4 = -5.061571874
    b5 = 0.009839526
    b6 = -0.095533007
    b7 = 0
  } else if(spec == "NWC"){
    a0 = 0.86118766
    a1 = 0.98152118
    a2 = 0.0568203
    b1 = 0.40717678
    b2 = -0.05482572
    b3 = 0.47809459
    b4 = -1.32512447
    b5 = 0.1538487
    b6 = -0.53687808
    b7 = 0
  } else if(spec == "TL"){
    a0 = 0.762977581
    a1 = 0.979320526
    a2 = 0.122788251
    b1 = 0.245935863
    b2 = -0.564901858
    b3 = 0.666790795
    b4 = -0.072877893
    b5 = 0.143651488
    b6 = -0.791188037
    b7 = 0
  } else { # temporarily added to account for unspecified species
    a0 = 0.762977581
    a1 = 0.979320526
    a2 = 0.122788251
    b1 = 0.245935863
    b2 = -0.564901858
    b3 = 0.666790795
    b4 = -0.072877893
    b5 = 0.143651488
    b6 = -0.791188037
    b7 = 0
  }
  
  
  
  heights <- seq(0, HT, length.out = 100)[-100]
  L <- HT / 100
  treeVolume <- 0
  
  for(i in heights){
    H1 = i
    H2 = i + L
    
    if(HT - H1 < 0.001) dib1 <- 0 else dib1 <- KozakExp02Mod(H1, HT, DBHO)
    if(HT - H2 < 0.001) dib2 <- 0 else dib2 <- KozakExp02Mod(H2, HT, DBHO)
    
    treeVolume <- treeVolume + Smalians(dib1, dib2, L * 100)/1000000
  }
  
  treeVolume <- ifelse(DBHO == 0, 0, treeVolume) #added
  return(treeVolume)
}

treeVolume <- Vectorize(treeVolume)

bs_data <- read.table(file = "clipboard", sep = "\t", header=TRUE) %>% 
  mutate(spec = 'BS') 
bs_data %>% mutate(volume = treeVolume(Tot.HT..m., DBH..cm., spec))


