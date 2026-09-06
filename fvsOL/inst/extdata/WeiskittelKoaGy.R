#Total Height =f(DBH (in), BAL (ft2/ac), BAPA (ft2/ac)), rain (mm), temp (C)
koa.HT=function(DBH,BAL,BAPA,rain=NA,temp=NA){
  if(is.na(rain) | is.na(temp)){
    b0= 69.78919  #4.410542 4409 15.82327   0e+00
    b1=  0.06144  #0.010558 4409  5.81912   0e+00
    b2=  0.85111  #0.018770 4409 45.34271   0e+00
    b4= -0.05703  #0.011801 4409 -4.83256   0e+00
    b5= -0.11266  #0.017304 4409 -6.51032   0e+00
    b6=  0.02221  #0.005903 4409  3.76238   2e-04 
    HT = 4.5 + (b0) * (1 - exp(-b1 * DBH))^(b2 + b4 * log(BAL + 1) + b5 * log((BAPA/100)+1) 
                                            + b6*log(BAL*BAPA+1))
  }
  else{
    b0= 373.2487 #139.26464 6007  2.68014  0.0074
    b1=  2.112872e-07   # 0.00000 6007  0.06321  0.9496
    b2=   4.1991  # 0.11302 6007 37.15530  0.0000
    b3= -25.4689  #13.08578 6007 -1.94631  0.0517
    b4=  -0.0189  # 0.00190 6007 -9.91404  0.0000
    b5=  -0.0248  # 0.00278 6007 -8.93117  0.0000
    b6=   0.0063  # 0.00093 6007  6.75460  0.0000
    b7=  -0.0325  # 0.01614 6007 -2.01229  0.0442
    b8=   0.5438  # 0.16525 6007  3.29089  0.0010
    HT = 4.5 + (b0+b3*log(rain*temp)) * (1 - exp(-b1* DBH^b2))^(b8+b4 * log(BAL + 1) + b5 * log((BAPA/100)+1) 
                                                                + b6*log(BAL*BAPA+1)+ b7*log(rain*temp+1))
  }
  return(HT=HT)
}

koa.HT(20,0,0,1500,16)
koa.HT(20,0,0)


HT=expand.grid(DBH=seq(0.1,100),BAL=seq(0,100,10),BAPA=seq(5,200,50))
HT$HT=koa.HT(HT$DBH,HT$BAL,HT$BAPA,1500,16)

plot(HT[HT$BAL==10 & HT$BAPA==50,]$DBH,HT[HT$BAL==10 & HT$BAPA==50,]$HT,type='l',
     xlab='DBH',ylab='HT')

#Height to crown base
koa.HCB=function(DBH,HT,BAL,BAPA){
  b0= -1.2329111 #0.4461396 372 -2.763510  0.0060
  b1= -0.2221513 #0.4413884 372 -0.503301  0.6151
  b2=  0.2485774 #0.1229929 372  2.021070  0.0440
  b3=  0.0014994 #0.0008839 372  1.696230  0.0907
  b4=  0.3421674 #0.0992084 372  3.448975  0.0006
  HCB = HT/(1 + exp(b0 + b1 * sqrt(HT/100) + b2*log(HT/DBH) + b3 * sqrt(BAL*BAPA + 1) + b4 * log(BAPA + 1)))
  return(HCB=HCB)
}

koa.HCB(10,40,10,100)

HT$HCB=koa.HCB(HT$DBH,HT$HT,HT$BAL,HT$BAPA)
plot(HT[HT$BAL==10 & HT$BAPA==55,]$DBH,HT[HT$BAL==10 & HT$BAPA==55,]$HT,type='l',
     xlab='DBH',ylab='HT/HCB')
points(HT[HT$BAL==10 & HT$BAPA==55,]$DBH,HT[HT$BAL==10 & HT$BAPA==55,]$HCB,type='l',lty=2)

require(ggplot2)
theme_set(theme_bw())
HT %>% filter(BAL==10, BAPA==55) %>% 
  ggplot()+
  geom_line(aes(x=DBH, y=HCB))+ 
  geom_line(aes(x=DBH, y=HT), linetype='dashed')+
  labs(y='Height')

#Monthly diameter increment
koa.dDBH=function(DBH,BAL,BAPA,rain=NA,temp=NA){
  if(is.na(rain) | is.na(temp)){
    b0 = -1.3223315# 0.20244992 2354  -6.531647  0.0000
    b1 =  0.2260716# 0.11546740 2354   1.957882  0.0504
    b2 = -0.2412006# 0.02857518 2354  -8.440913  0.0000
    b3 = -0.0002227# 0.00003130 2354  -7.114750  0.0000
    b4 = -0.8247849# 0.05729640 2354 -14.395056  0.0000
    b5 =  0.0210812# 0.00539649 2354   3.906472  0.0001
    b6 =  0.2809946# 0.02774182 2354  10.128919  0.0000
    dDBH = exp(b0+b1*log(DBH+1)+b2*DBH+b3*(BAL^2/log(DBH+5))+b4*log(BAL+1)+
                 b5*sqrt(BAPA*DBH)+b6*log(BAL*BAPA+1))
  }
  else{
    b0= -86.87632 #31.315169 3615  -2.774257  0.0056
    b1=   0.46092  #0.105180 3615   4.382214  0.0000
    b2=  -0.31036  #0.025890 3615 -11.987680  0.0000
    b3=  -0.00020  #0.000025 3615  -7.998912  0.0000
    b4=  -0.82359  #0.042371 3615 -19.437708  0.0000
    b5=   0.03270  #0.004141 3615   7.895308  0.0000
    b6=   0.23718  #0.019672 3615  12.056988  0.0000
    b7=   9.27207  #3.408181 3615   2.720534  0.0065
    b8=  -0.33883  #0.127784 3615  -2.651559  0.0080
    dDBH = exp(b0+b1*log(DBH+1)+b2*DBH+b3*(BAL^2/log(DBH+5))+b4*log(BAL+1)+
                 b5*sqrt(BAPA*DBH)+b6*log(BAL*BAPA+1)+b7*log(rain*temp)+b8*(rain*temp/1000))
  }
  return(dDBH=dDBH)
}

koa.dDBH(10,10,100,1500,15)*12

p.dDBH=expand.grid(DBH=seq(0.1,100,1),BAL=seq(0,100,10),BAPA=seq(0,200,10))
p.dDBH$dDBH=koa.dDBH(p.dDBH$DBH,p.dDBH$BAL,p.dDBH$BAPA,1500,15)

plot(p.dDBH[p.dDBH$BAL==10 & p.dDBH$BAPA==100,]$DBH,p.dDBH[p.dDBH$BAL==10 & p.dDBH$BAPA==100,]$dDBH,
     type='l',xlab='DBH (in)',ylab='dDBH (in/mo)',xlim=c(0,25))


#Monthly height increment
koa.dHT=function(DBH,HT,BAL,BAPA,rain=NA,temp=NA){
  if(is.na(rain) | is.na(temp)){
    b0= -1.2984999 #0.10897595 2154 -11.915472  0.0000
    b1=  1.0873012 #0.11443753 2154   9.501264  0.0000
    b2= -0.2133439 #0.02349380 2154  -9.080857  0.0000
    b3= -0.0001268 #0.00002563 2154  -4.946981  0.0000
    b5=  0.0114348 #0.00426040 2154   2.683977  0.0073
    b6= -0.0957501 #0.01055932 2154  -9.067830  0.0000
    b4=0
    dHT = exp(b0+b1*log(DBH+1)+b2*DBH+b3*(BAL^2/log(DBH+5))+b4*log(BAL+1)+
                b5*sqrt(BAPA*DBH)+b6*log(BAL*BAPA+1))
  }
  else{
    b0= 163.54729  #93.92480 3263   1.741258  0.0817
    b1=   2.80342   #0.17512 3263  16.008452  0.0000
    b2=  -0.21923   #0.01181 3263 -18.567832  0.0000
    b3=  -0.00022   #0.00003 3263  -8.327186  0.0000
    b4=  -0.31141   #0.07390 3263  -4.214224  0.0000
    b5=   0.08924   #0.00428 3263  20.873448  0.0000
    b6=   0.10094   #0.03809 3263   2.650235  0.0081
    b7=  13.21554   #7.43355 3263   1.777824  0.0755
    b8= -17.61490   #9.82302 3263  -1.793226  0.0730
    dHT = exp(b0+b1*log(HT+1)+b2*HT+b3*(BAL^2/log(DBH+5))+b4*log(BAL+1)+
                b5*sqrt(BAPA*DBH)+b6*log(BAL*BAPA+1)+b7*sqrt(rain*temp/1000)+b8*log((rain*temp)^2/1000))
  }
  return(dHT=dHT)
}

koa.dHT(10,50,10,100)

p.dDBH$HT=koa.HT(p.dDBH$DBH,p.dDBH$BAL,p.dDBH$BAPA,1500,15)
p.dDBH$dHT=koa.dHT(p.dDBH$DBH,p.dDBH$HT,p.dDBH$BAL,p.dDBH$BAPA,1500,15)

plot(p.dDBH[p.dDBH$BAL==50 & p.dDBH$BAPA==150,]$HT,
     p.dDBH[p.dDBH$BAL==50 & p.dDBH$BAPA==150,]$dHT,
     type='l',xlab='HT',ylab='dHT (ft/mo)')

#Annual % Survival
koa.alive=function(DBH,HT,BAL,BAPA,rain=NA,temp=NA,YIP=1){
  if(is.na(rain) | is.na(temp)){
    b0 =  7.320309     #0.400527   7.637 2.22e-14 ***
    b1 = 0.207024    #0.092536  -3.752 0.000176 ***
    b2 =  -0.685661   #0.004863   4.191 2.78e-05 ***
    b3 = 0.710481  # 0.003332  -4.481 7.45e-06 ***
    b4 = -0.023592   # 0.036881  -2.566 0.010298 *
    b5 = -0.945985 
    lp = b0+b1*(DBH)+b2*log(DBH^2)+b3*log(HT/DBH+1)+b4*((BAL+1)/log(DBH+1))+b5*log(BAPA)
    ps =  (1 / (1 + exp(-lp)))^(1/YIP)
  }
  else{
    b0 =   1.419e+03  #4.166e+01  34.068  < 2e-16 ***
    b1 =   -6.523e-01  #1.050e-01  -6.212 5.22e-10 ***
    b2 =   3.946e-02  #5.703e-03   6.918 4.57e-12 ***
    b3 =  5.664e-01  #3.242e-01   1.747 0.080618 .  
    b4 =  -2.556e-02  #3.027e-03  -8.441  < 2e-16 ***
    b5 =  -1.086e-01  #3.296e-02  -3.294 0.000988 ***
    b6 = -1.711e+02  #5.062e+00 -33.801  < 2e-16 ***
    b7 = 2.019e+00  #6.237e-02  32.377  < 2e-16 ***
    lp = b0+b1*(DBH)+b2*log(DBH^2)+b3*log(HT/DBH+1)+b4*((BAL+1)/log(DBH+1))+b5*sqrt(BAPA)+b6*log(temp*rain)+b7*sqrt(temp*rain)
    ps =  (1 / (1 + exp(-lp)))^(1/YIP)
  }
  return(ps=ps)  
}

koa.alive(10,100,100,150,1500,15,YIP=1)

pmort=expand.grid(DBH=seq(0.1,100,1),BAL=seq(0,100,10),BAPA=seq(0,150,10))
pmort$HT=koa.HT(pmort$DBH,pmort$BAL,pmort$BAPA,1500,15)
pmort$Alive=koa.alive(pmort$DBH,pmort$HT,pmort$BAL,pmort$BAPA,1500,15,1)
hist(pmort$Alive)

plot(pmort$DBH,pmort$Alive,pch=".",xlim=c(0,30))
