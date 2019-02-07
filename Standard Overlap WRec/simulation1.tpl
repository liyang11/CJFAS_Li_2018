//FILE CREATED ON 3/3/2011 BY TRAVIS BRENDEN & KYLE MOLTON & WEIHAI LUI
//DEVELOPED BY YANG LI FOR WHITEFISH MSE SIMULATION PROJECT
//LAST MODIFIED 12/12/2016

GLOBALS_SECTION
  #include "admodel.h"
  #include <time.h>
  #include "qfclib.h"
  //define constant variable 
  const double MathPI = 3.141592654;        // or using MI_PI
  double eps = 1e-30; 
  const double MathE = 2.71828183; 

DATA_SECTION
  !! ad_comm::change_datafile_name("simulation1.dat");
  init_int numSim                           // num sims
  init_int numYear                          // total yrs
  init_int numAge                           // total ages
  init_int numPop                           // num pops
  init_int numProd                          // num of productivity levels
  init_int numMort                          // num of target total mortality levels
  init_int numSce                           // num mixing and productivity scenaros
  init_int numRV
  init_int numDQ
  init_int numMU
  init_int assessYr                         // ass starting yr=21
  init_int assessStepYr                     // ass every 3 year
  init_vector natMort(1,numPop)             // M 
  init_matrix Mottable(1,numMort,1,numPop)  // target A (According to control rule)
  init_vector lvbGrowth(1,3)                // LVB growth paras(Linf, K, to)
  init_vector weightPar(1,2)                // weihgt condition paras(alpha, beta)
  init_number propFemale                    // prop. of fem
  init_number eggsKg                        // eggs/kg   
  init_vector logistMat(1,3)                // logistic maturity paras(maturity at inf length,curv param, length at inflect of logistic model) 
  init_matrix recruitPar(1,numProd,1,2)     //fecundity alpha, density depend. beta
  init_matrix recruitErrPar(1,numRV,1,2)    //recruAut,recruStd
  init_matrix recruitInit(1,numMort,1,numProd)    //intial Rec
  init_vector gammaSel(1,2)                 // selectivity at age
  init_number q                             // catchability for 4 independent pops
  init_matrix move(1,numSce,1,numPop)       // movement parameters_ stay rates for different scenarios
  init_matrix ProdTab(1,numSce,1,numPop)  // productivity levels for different scenarios
  init_vector implSdTAC(1,numPop)           // implement error for TAC(sd)
  init_vector effSampleSizeTab(1,numDQ)        // effective sample size 50 or 200
  init_matrix stdObsCatchTab(1,numDQ,1,numPop)         // std error for observed catch for assessment model
  init_matrix sdEffortTab(1,numDQ,1,numPop)            // std error for fishing effort
  init_adstring outFile                     // output file name
  init_adstring outPerform1                 // output file name for independent assessment model
  init_adstring outPerform2                 // output file name for pooled assessment model
  init_int useAllYr                         // use all year data for assessment, 1 = all,0 != all(provide yr# later) =0
  init_int nYrAssess                        // use how many yr data for assess =20
  init_int assessFAge                       // ass 1st age =3
  init_int assessLAge                       // ass last age =12
  init_int lastNumYrMort                    // how many last yrs mortality as mean mortality for TAC
  init_int lastNumYrRec                     // number of last yrs on first age as mean recruitment for TAC
  init_int forLinux
  init_matrix MixVtable(1,3,1,3) 
  init_int useRnd                           // use fixed or random seed_ we used fixed seed here
  init_int seed                             // if using fixed, then use this seed as provided in dat file,ow using system time as seed
 
  init_int model                            // model 1,2,3 represent three assessment methods, model 4 represent overlap but with catch-at-pop+age
  init_int sceLev
  init_int motLev
  init_int RecLev                           
  init_int DQLev
  init_int MixVLev
  init_vector test(1,3)



  // ***************************************following derived directly from the input data
  vector targetA(1,numPop)            
  vector initAbund(1,numPop)                // initial abundance /population
  number Linf                               // growth para
  number K                                  // growth para
  number t0                                 // growth para
  number wgtAlpha                           // weight condition para
  number wgtBeta                            // weight condition para
  number matInf                             // maturity parameter
  number matK                               // maturity parameter
  number matFlect                           // maturity parameter
  vector targetTotMort(1,numPop)            // vector for Z's
  vector recruAlpha(1,numPop)               // logistic recruit parameter
  vector recruAlpha_Pri(1,numPop)               // logistic recruit parameter  
  vector recruBeta(1,numPop)                // logistic recruit parameter
  vector recruStd(1,numPop)                 // logistic recruit parameter
  vector recruAut(1,numPop)                 // logistic recruit parameter
  vector rec_St_Var(1,numPop)               // logistic recruit parameter
  vector ages(1,numAge)                     // age vector  (12)
  vector years(1,numYear)                   // year vector (100)
  vector len(1,numAge)                      // fish length at age 
  vector weight(1,numAge)                   // weight at length
  vector maturity(1,numAge)                 // maturity at length
  vector select(1,numAge)                   // selectivity at length or age
  vector fishMort(1,numPop)                 // fishing mortality/population
  vector recruLev(1,numPop)                 // translate productivity levels to usable matrix (col=7 rec paras, row=pops)
 


  //following holding the sim result data
  5darray subPop(1,numSim,1,numYear,1,numPop,1,numPop,1,numAge)                   // sub population abundance
  5darray subMixPop(1,numSim,1,numYear,1,numPop,1,numPop,1,numAge)                // sub population abundance only after mixing
  5darray subCatch(1,numSim,1,numYear,1,numPop,1,numPop,1,numAge)                 // harvest number for sub pop
  4darray areaCatch(1,numSim,1,numYear,1,numPop,1,numAge)                         // catch at age for each area
  4darray pop(1,numSim,1,numYear,1,numPop,1,numAge)                               // population abundance at age (at the end of each yr)
  4darray mixStock(1,numSim,1,numYear,1,numPop,1,numAge)                          // area abundance at age
  4darray subSSB(1,numSim,1,numYear,1,numPop,1,numPop)                            // cross all age, (pop,area)
  4darray subSpawn(1,numSim,1,numYear,1,numPop,1,numPop)                          // sub population abundance
  3darray SSB(1,numSim,1,numYear,1,numPop)                                        // cross all pop
  3darray Spawn(1,numSim,1,numYear,1,numPop)                                      // cross all pop
  3darray Recruit(1,numSim,1,numYear,1,numPop)                                    // cross all pop
  3darray effort(1,numSim,1,numYear,1,numPop)                                     // fishing effort for each area
  4darray popBegYr(1,numSim,1,numPop,1,numYear,1,numAge)                          // pop abundance at the beginning of yr

  4darray mixStock1(1,numSim,1,numPop,1,numYear,1,numAge)                         //switch second and 3rd dimension as mixStock 
  4darray areaCatch1(1,numSim,1,numPop,1,numYear,1,numAge)                        //switch second and 3rd dimension as areaCatch    
  4darray areaYield1(1,numSim,1,numPop,1,numYear,1,numAge)                        //switch second and 3rd dimension as areaYield
  4darray pop1(1,numSim,1,numPop,1,numYear,1,numAge)                              //switch second and 3rd dimension as pop
  
  int assessOK                                                                    // control which yr to do assessment
  vector readInTAC(1,numPop)                                                      // from assessment model
  vector realTAC(1,numPop)                                                        // add the implementation error
  vector newF(1,numPop)                                                           // estimate fishing mort for each area, from Newton-raphson on TAC
  vector warnConverg(1,numPop)                                                    // keep track of convergence warning for assessment run
  vector warnHessian(1,numPop)                                                    // keep track of hessian  warning for assessment run
  int totAssess                                                                   // total assessment running
  3darray annualTAC(1,numSim,1,numYear,1,numPop)
  3darray annualYield(1,numSim,1,numPop,1,numYear)
  matrix indicator(1,numSim,1,numPop)
  matrix hess_ind(1,numSim,1,numPop)
  matrix counter(1,numSim,1,numPop)
  matrix hess_count(1,numSim,1,numPop)
 
 
  
  //added by yang

  number catchability      
  3darray F(1,numSim,1,numYear,1,numPop)                                          // fishing mortality for each pop each yr                                              
  //number readInTACpooled                                                          // tmp TACpooled value from ass model         
  //matrix pooledeffort(1,numSim,1,numYear)
  //3darray pooledYield(1,numSim,1,numYear,1,numAge)    
  //3darray pooledBegYr(1,numSim,1,numYear,1,numAge)
  //3darray pooledcatch(1,numSim,1,numYear,1,numAge)
  //number pooledwarnConverg                                                        // keep track of convergence warning for assessment run
  //number pooledwarnHessian                                                        // keep track of hessian  warning for assessment run
  //vector pooledindicator(1,numSim)
  //vector pooledhess_ind(1,numSim)
  //vector pooledcounter(1,numSim)
  //vector pooledhess_count(1,numSim)
  //matrix pooledTAC(1,numSim,1,numYear)
  3darray yieldTolarea(1,numSim,1,numYear,1,numPop)
  //3darray catchTolarea(1,numSim,1,numYear,1,numPop)
  //3darray areaCPE(1,numSim,1,numYear,1,numPop) 
  //matrix TotalCPE(1,numSim,1,numYear)   
  //3darray CPEpar(1,numSim,1,numYear,1,numPop)



  3darray real_ssb(1,numSim,1,numPop,1,numYear)        //area specific abundance *weight
  3darray assess_ssb(1,numSim,1,numPop,1,numYear)                                
  4darray real_rec(1,numSim,1,numPop,1,numYear,1,nYrAssess)     
  4darray assess_rec(1,numSim,1,numPop,1,numYear,1,nYrAssess)   
  4darray real_F(1,numSim,1,numPop,1,numYear,1,nYrAssess)     
  4darray assess_F(1,numSim,1,numPop,1,numYear,1,nYrAssess)   


           
  3darray ssb_mre(1,numSim,1,numPop,1,numYear)         // readin MRE  from assessment model
  3darray ssb_mare(1,numSim,1,numPop,1,numYear)        // readin MARE from assessment model
  4darray rec_mre(1,numSim,1,numPop,1,numYear,1,nYrAssess)         // readin MRE  from assessment model
  //4darray rec_mare(1,numSim,1,numPop,1,numYear,1,nYrAssess)        // readin MARE from assessment model
  4darray F_mre(1,numSim,1,numPop,1,numYear,1,nYrAssess)         // readin MRE  from assessment model

  int totalAssessYrs

  //matrix pooled_real_ssb(1,numSim,1,numYear)          //Used to cal SSB mre, this is abundance *weight ~SSB
  //matrix pooled_assess_ssb(1,numSim,1,numYear)        //used to cal SSB mre, this is abundance *weight ~SSB
  //matrix pooled_ssb_mre(1,numSim,1,numYear)            // total abundance *weight for 4 areas mre
  //matrix pooled_ssb_mare(1,numSim,1,numYear)            // total abundance *weight for 4 areas mare    
  //matrix pooled_real_rec(1,numSim,1,numYear)            //Used to cal rec re
  //matrix pooled_assess_rec(1,numSim,1,numYear)          //used to cal rec re
  //matrix pooled_rec_mre(1,numSim,1,numYear)             // total rec for 4 pops re
  //matrix pooled_rec_mare(1,numSim,1,numYear)            // total rec for 4 pops are    
  //matrix pooledcatch_mre(1,numSim,1,numYear)
  //matrix pooledcatch_mare(1,numSim,1,numYear)  

  //3darray catch_mre(1,numSim,1,numPop,1,numYear)         // readin MRE  from assessment model
  //3darray catch_mare(1,numSim,1,numPop,1,numYear)        // readin MARE from assessment model


  //add by yang for check result
  matrix n_est(1,20,3,12) 
  vector N_rc_est(1,20)  
  vector new_rc_cch_est(1,20)  
  vector rc_s_est(3,12)
  //matrix rc_s_est(1,20,3,12)   
  vector Ftac_est(3,12)
  vector tacAge_est(1,10) 
  matrix movePar(1,numPop,1,numPop)    // movement parameters square matrix
  3darray R_err(1,numSim,1,numYear,1,numPop)     
  5darray areaCatchPop(1,numSim,1,numYear,1,numPop,1,numPop,assessFAge,assessLAge)                      //switch second and 3rd dimension as areaCatch    



  //added by yang on Feb 15 for chapter 4
  int effSampleSize
  vector stdObsCatch(1,numPop)
  vector sdEffort(1,numPop)
  3darray moveParVar(1,numYear,1,numPop,1,numPop)    // movement parameters square matrix when mixing varying annually
  3darray moveOp(1,numYear,1,numPop,1,numPop)    // movement parameters square matrix used in operating model
  matrix moveAss(1,numPop,1,numPop)    // movement parameters square matrix used in Assessment model
  matrix movePartmp(1,numPop,1,numPop);  
  vector mix_mu(1,numPop)  
  vector mix_std(1,numPop)
   
 // ********begin loc_cals to filled data in to targetTotMort,  growth paras,   weight condition paras,  maturity paras, logistic recruit paras... 
 LOC_CALCS

   targetA    = Mottable(motLev); 
   //cout<<"targetA"<<targetA<<endl;
   for(int i=1;i<=numPop;i++)  
   targetTotMort(i)=-log(1.0-targetA(i)); // set the targetTotMort to Z instead of A, ie make it instantaneous
   
   //following derived directly from the input data
   Linf       = lvbGrowth(1);         // growth parameter
   K          = lvbGrowth(2);
   t0         = lvbGrowth(3);
   wgtAlpha   = weightPar(1);         // weight condition parameter
   wgtBeta    = weightPar(2);
   matInf     = logistMat(1);         // maturity parameter
   matK       = logistMat(2);
   matFlect   = logistMat(3);
   
   recruLev   = ProdTab(sceLev);
   //added by yang on Feb 15 for chapter 4
   effSampleSize = effSampleSizeTab(DQLev);
   stdObsCatch = stdObsCatchTab(DQLev);
   sdEffort = sdEffortTab(DQLev);


   
   for(int i=1; i<=numPop; i++) {
    int ind=recruLev[i];  //varying from 1 to 3, correspond to three productivity levels
    recruAlpha_Pri[i] = recruitPar(ind,1);
    recruBeta[i]      = recruitPar(ind,2);
    recruAut[i]       = recruitErrPar(RecLev,1);
    recruStd[i]       = recruitErrPar(RecLev,2);  
    initAbund[i]      = recruitInit(motLev,ind); 
    rec_St_Var[i]     = square(recruStd[i])/(1-square(recruAut[i]));
    recruAlpha[i]     = recruAlpha_Pri[i]*mfexp(-rec_St_Var[i]/2);  // read in value was alpha'= alpha * exp(var/2) for E(R)=alpha'*S*exp(-beta*S)     
   }
   

   ages.fill_seqadd(1,1);                                 // fill age vector
   years.fill_seqadd(1,1);                                // fill year vector
   fishMort  = targetTotMort-natMort;
   len       = Linf*(1.-mfexp(-K*(ages-t0)));             // fish length at age, vector
   weight    = wgtAlpha*pow(len,wgtBeta);                 // weight at length
   maturity  = matInf/(1.+mfexp(-matK*(len- matFlect)));  // maturity at length
   maturity(1,3)=0;
   select.initialize();
   for(int i=assessFAge; i<=numAge; i++)  select(i) = pow(i,gammaSel(1))*mfexp(-gammaSel(2)*i); // selectivity at age/len
   select /= select(10);                               // changed from 8 to 10

   




   //cout<<"targetA"<<targetA<<endl;
   //cout<<"recruLev"<<recruLev<<endl;
   //cout<<"recruAlpha_Pri"<<recruAlpha_Pri<<endl;
   //cout<<"recruAlpha"<<recruAlpha<<endl;
   //cout<<"recruBeta"<<recruBeta<<endl;
   //cout<<"initAbund"<<initAbund<<endl;
   //cout<<"recruAut"<<recruAut<<endl;
   //cout<<"recruAut"<<"recruStd"<<recruStd<<endl;

   //cout<<"test"<<test<<endl;
   //exit(11);

 END_CALCS


  4darray mixBegYr(1,numSim,1,numYear,1,numPop,1,numAge) // mixture abundance for each area at beginning of yr
  !! totalAssessYrs=int(ceil(double(numYear-assessYr+1)/assessStepYr)); // used for defined variables for mre 
  ivector assessYearVec(1,totalAssessYrs)
  !! assessYearVec.fill_seqadd(assessYr,assessStepYr); 
  //added to record grad
  //3darray Grad(1,numSim,1,numPop,1,numYear)        // readin Grad from independent assessment model
  //matrix pooledGrad(1,numSim,1,numYear)            // readin Grad from pooled assessment model
  
  
  
PARAMETER_SECTION
  init_number dummy
  objective_function_value fake


PRELIMINARY_CALCS_SECTION 
  // ******************************************************initialize vectors
  // ******************************************************change random seed based on value from data file
  indicator =0.;
  hess_ind  =0.;
  counter   =0.;
  hess_count=0.;
  if(useRnd == 1) seed=(unsigned)time(0);                 // using system time as seed for total random
  random_number_generator rnd(seed); 
  subPop.initialize();      pop.initialize();       subCatch.initialize();  areaCatch.initialize(); 
  mixStock.initialize();    subSSB.initialize();    SSB.initialize();       warnConverg.initialize();
  warnHessian.initialize(); popBegYr.initialize();  subSpawn.initialize();  Spawn.initialize();
  Recruit.initialize(); F.initialize();yieldTolarea.initialize();
  //pooledBegYr.initialize();pooledTAC.initialize();tacAge_est.initialize(); 
  //catchTolarea.initialize(); areaCPE.initialize();TotalCPE.initialize();CPEpar.initialize();
  mixBegYr.initialize(); subMixPop.initialize();
  real_ssb.initialize(); assess_ssb.initialize(); //pooled_real_ssb.initialize(); pooled_assess_ssb.initialize(); 
  real_rec.initialize(); assess_rec.initialize();
  real_F.initialize(); assess_F.initialize();
  ssb_mre.initialize();ssb_mare.initialize();
  rec_mre.initialize();//rec_mare.initialize();
  F_mre.initialize();
  //assess_rec.initialize(); pooled_real_rec.initialize(); pooled_assess_rec.initialize(); 
  //pooled_ssb_mre.initialize();pooled_ssb_mare.initialize(); catch_mre.initialize(); catch_mare.initialize();


  //pooledcatch_mre.initialize();pooledcatch_mare.initialize();pooled_rec_mre.initialize();pooled_rec_mare.initialize();

  R_err.initialize();

  //pooledindicator   =0.;
  //pooledhess_ind    =0.;
  //pooledcounter     =0.;
  //pooledhess_count  =0.;
  //pooledwarnConverg =0.;
  //pooledwarnHessian =0.;
  


   dvector stay(1,numPop);
   dvector stayVar(1,numPop); 
   for(int i=1; i<=numPop; i++)   {
    stay(i)= 1- move(sceLev,i);
    for(int j=1; j<=3; j++){
     if (stay(i)== MixVtable(j,1)) {       //find the correponding std and mu for a given stay rate, which is defined in move 
         mix_std(i)=MixVtable(j,2);
         mix_mu(i)=MixVtable(j,3);   
     }
    } 
  }
  mixParameter(stay);
  movePar=movePartmp;

 dvector staytmp(1,numPop);
 for(int y=1; y<=numYear; y++){
   for(int i=1; i<=numPop; i++)   {
 	 //stayVar(i)=mfexp(mix_std(i)*randn(rnd)+ mix_mu(i))/(mfexp(mix_std(i)*randn(rnd)+ mix_mu(i))+1);   //stay rate for annual varying mixing rate scenarios
     staytmp(i)=rnorm(mix_mu(i),mix_std(i),rnd);
     stayVar(i)=exp(staytmp(i))/(exp(staytmp(i))+1);   //stay rate for annual varying mixing rate scenarios
   }
   //cout<<"stayVar"<<stayVar<<endl;
   mixParameter(stayVar);
   moveParVar(y)=movePartmp;
 }  


 if (MixVLev==1) {
   for (int y=1; y<=numYear; y++) moveOp(y)=movePar;
   moveAss=movePar;
  }
  if (MixVLev==2) {
   for (int y=1; y<=numYear; y++) moveOp(y)=movePar;
   //moveAss=moveParVar; // defined later in the simulation model
  }
  if (MixVLev==3) {
   for (int y=1; y<=numYear; y++) moveOp(y)=moveParVar(y);
   moveAss=movePar;
  }

  //cout<<"moveTrue"<<endl<<movePar<<endl;
  //cout<<"stayVar"<<stayVar<<endl; 
  //cout<<"model"<<endl<<model<<endl;
  //cout<<"scenarios"<<endl<<MixVLev<<endl;
  //cout<<"moveAss"<<endl<<moveAss<<endl;
  //for (int y=80; y<=numYear; y++) cout<<"moveOp"<<" for Year "<<y<<endl<<moveOp(y)<<endl;
  //cout<<"test"<<test<<endl;
  //exit(11);





   // *************************************data genetation begins********************************
   // year 1
  totAssess=0;
  for(int i=1;i<=numSim;i++){          // **simulation loop  i
    assessOK=assessYr;                 // set assessment start yr 
    for(int j=1;j<=numYear;j++){       // **year loop j
    //*********************************   first year data simulation 
      if(j == 1){ 
        for(int n=1;n<=numPop;n++){    // **population  loop n
          double  tmpPop=initAbund(n); // initial bundance for each pop at first year, first age
          popBegYr(i,n,j,1)=tmpPop;    // keep track of abundance at beginning of yr, first age hold initial abund,other ages as 0
	       
          for(int k=1;k<=numPop;k++){  // **area  loop k
            for(int m=1;m<=numAge;m++){// **age loop m
              if(m == 1){              // age=1                              
                subMixPop(i,j,n,k,m)=tmpPop*moveOp(j,n,k);                                          // only apply mixing on age 0,1st                
                subPop(i,j,n,k,m)=tmpPop*moveOp(j,n,k)*mfexp(-(natMort(k)+select(m)*fishMort(k)));  // apply mortality on this age0, 1st yr ; after fishing: P'= P * exp(-(M+F))
                subCatch(i,j,n,k,m)=(select(m)*fishMort(k)/(natMort(k)+select(m)*fishMort(k)))*tmpPop*moveOp(j,n,k)*(1.-mfexp(-(natMort(k)+select(m)*fishMort(k))));
                                       // catch= P*(1-exp(-(M+F)))*(F/(M+F)) 
              }
              else {                   // age>=2
                                       // initial other ages(based on previous age),1st yr; assume P'[age2, 0yr]=P'[age1,1yr]            
               subPop(i,j,n,k,m)=subPop(i,j,n,k,m-1)*mfexp(-(natMort(k)+select(m)*fishMort(k)));   // eg. P`[age2]= P`[age1] * exp(-(M+F))
               subCatch(i,j,n,k,m)=(select(m)*fishMort(k)/(natMort(k)+select(m)*fishMort(k)))*subPop(i,j,n,k,m-1)*(1.-mfexp(-(natMort(k)+select(m)*fishMort(k)))); 
                                       // eg. catch[age2]= P`[age1]*(1-exp(-(M+F)))*(F/(M+F))
              }
              double tmp=subPop(i,j,n,k,m)*eggsKg*propFemale*maturity(m)*weight(m);  
              subSSB(i,j,n,k)+= tmp;   // accumulate SSB, cross all age 
              subSpawn(i,j,n,k)+=(tmp/eggsKg)/weight(m);
              pop(i,j,n,m)+=subPop(i,j,n,k,m);         // sum up all area to each pop
              if(m>=2) popBegYr(i,n,j,m)=pop(i,j,n,m-1);
              mixBegYr(i,j,k,m)+=subMixPop(i,j,n,k,m); // sum up all population to each area just after mixing
              mixStock(i,j,k,m)+=subPop(i,j,n,k,m);    // sum up all population to each area            
              areaCatch(i,j,k,m)+=subCatch(i,j,n,k,m); // sum up all subgroup catch to each area          
          }                                         // end age m             
            SSB(i,j,n)+=subSSB(i,j,n,k);            // sum up all area to each population SSB
            Spawn(i,j,n)+=subSpawn(i,j,n,k);
            F(i,j,k)= fishMort(k);       
            Recruit(i,j,n)=popBegYr(i,n,j,assessFAge);   
        }                                          // end area k 
        R_err(i,j,n) =0;
      }   
     //added by Yang on Jan11,2017, for chapter4 population composition of catch 
      for(int k=1;k<=numPop;k++)
      	for(int n=1;n<=numPop;n++)
          for (int m=assessFAge;m<=assessLAge;m++)
             areaCatchPop(i,j,k,n,m)=subCatch(i,j,n,k,m);


      for(int k=1;k<=numPop;k++) {
        effort(i,j,k)=1./q*fishMort(k)*mfexp((sdEffort(k)*randn(rnd))-(0.5*square(sdEffort(k))));      // fishing effort:  E= F*exp(sd-0.5*sigma^2) /q
        for(int m=1;m<=numAge;m++)  yieldTolarea(i,j,k) += areaCatch(i,j,k,m)*weight(m);    
      }
       
    } 
    // *************************************************end of first year data simulation 
    //year 2-21
    // **************************************************************after first year &  before starting assessment year data simulation
      if(j > 1 && j<=assessYr){ 
        for(int n=1;n<=numPop;n++){    
          dvector tmpPop(1,numAge);
          R_err(i,j,n)= recruStd(n)*randn(rnd)+ recruAut(n)*R_err(i,j-1,n);
          tmpPop(1)=recruAlpha(n)*SSB(i,j-1,n)*mfexp(-recruBeta(n)*SSB(i,j-1,n))*mfexp(R_err(i,j,n)); //using SSB from previous year to cal recruitment
                                                                  //using SSB from previous year to cal recruitment
     
          

      // if (j==2){
      //     cout<<"s"<<SSB(i,j-1,n)<<endl;
      //	  cout<<"r"<<Recruit(i,j,n)<<endl;
      //	  cout<<"n"<<n<<endl;
      //	  cout<<"alpha"<<recruAlpha(n)<<endl;
      //	  cout<<"beta"<<recruBeta(n)<<endl;
      //	  cout<<"err"<<setprecision(4)<<setw(10)<<R_err(i,j,n)<<endl;
      //       exit(9);
      //    }


          for(int t=2;t<=numAge;t++) tmpPop(t)=pop(i,j-1,n,t-1);  //initialize abun at age: eg. initial tmpPop[age2,yr2]=Pop[age1, end of yr1] (no lastyr class) 
          tmpPop(numAge)+=pop(i,j-1,n,numAge);                    //initialize last yr class: tmp[age12,yr2] =Pop[age11, end of yr1] + Pop[age12, end of yr1]    
          popBegYr(i,n,j)=tmpPop;                                 //keep track of abundance at all age at beginning of yr
          Recruit(i,j,n)=tmpPop(assessFAge);                               // record recruitment
          for(int k=1;k<=numPop;k++){  
            for(int m=1;m<=numAge;m++){ 
              subMixPop(i,j,n,k,m)=tmpPop(m)*moveOp(j,n,k);
              subPop(i,j,n,k,m)=tmpPop(m)*moveOp(j,n,k)*mfexp(-(natMort(k)+select(m)*fishMort(k)));
              subCatch(i,j,n,k,m)=(select(m)*fishMort(k)/(natMort(k)+select(m)*fishMort(k)))*tmpPop(m)*moveOp(j,n,k)*(1.-mfexp(-(natMort(k)+select(m)*fishMort(k))));
              double tmp=subPop(i,j,n,k,m)*eggsKg*propFemale*maturity(m)*weight(m);  
              subSSB(i,j,n,k)+= tmp; 
              subSpawn(i,j,n,k)+=(tmp/eggsKg)/weight(m);
              pop(i,j,n,m)+=subPop(i,j,n,k,m); 
              mixBegYr(i,j,k,m)+=subMixPop(i,j,n,k,m); 
              mixStock(i,j,k,m)+=subPop(i,j,n,k,m);            
              areaCatch(i,j,k,m)+=subCatch(i,j,n,k,m); 
            } 
            SSB(i,j,n)+=subSSB(i,j,n,k);  
            Spawn(i,j,n)+=subSpawn(i,j,n,k);
            F(i,j,k)= fishMort(k);
          }            // end area k                
        }              // end pop n 

     //added by Yang on Jan11,2017, for chapter4 population composition of catch 
      for(int k=1;k<=numPop;k++)
      	for(int n=1;n<=numPop;n++)
          for (int m=assessFAge;m<=assessLAge;m++)
             areaCatchPop(i,j,k,n,m)=subCatch(i,j,n,k,m);


        for(int k=1;k<=numPop;k++)  {
         effort(i,j,k)=1./q*fishMort(k)*mfexp((sdEffort(k)*randn(rnd))-(0.5*square(sdEffort(k))));      // fishing effort:  E= F*exp(sd-0.5*sigma^2) /q
         for(int m=1;m<=numAge;m++)  yieldTolarea(i,j,k) += areaCatch(i,j,k,m)*weight(m); 
        }
      } 
      // end  first year &  before starting assessment year data simulation
      // *********************************************************************
            

      
      
// ****************assessment begin**************
// model=1 :: independent populatation assessment model
// model=2 :: pooled population assessment model (cpe3)
// model=3 :: overlap population assessment model
// model=4 :: overlap population with catch at pop and age
// Begins in year 22
// ***********************************************

  if(j > assessYr){                        // use new F for year after assessment started
        dmatrix tmpPop(1,numPop,1,numAge); // build the population abundance at start of the yr
        dmatrix tmpPopMove(1,numPop,1,numAge);
        for(int n=1;n<=numPop;n++){ 
           R_err(i,j,n)= recruStd(n)*randn(rnd)+ recruAut(n)*R_err(i,j-1,n);
          tmpPop(n,1)=recruAlpha(n)*SSB(i,j-1,n)*mfexp(-recruBeta(n)*SSB(i,j-1,n))*mfexp(R_err(i,j,n)); // using SSB from previous year       
          for(int t=2;t<=numAge;t++) 
          tmpPop(n,t)=pop(i,j-1,n,t-1);    
          tmpPop(n,numAge)+=pop(i,j-1,n,numAge);
          popBegYr(i,n,j)=tmpPop(n);
 	  Recruit(i,j,n)=tmpPop(n,assessFAge);
        }    // end pop n
        dvector tmp(1,numPop);
        if (model!=2){ // begin model 1,3,4: population specific TAC
         for(int k=1;k<=numPop;k++){ 
           tmp(k)=mfexp((implSdTAC(k)*randn(rnd))-(0.5*square(implSdTAC(k))));
           realTAC(k)=readInTAC(k)*tmp(k); // *error
           annualTAC(i,j,k)=readInTAC(k);
         }  // end area k
        }   // end model 1
        

        cout<< "######Sim, Year: " <<i<<" "<<j<<endl;
        //cout<< "TAC for each area(without error): " << annualTAC(i,j)<< "  " <<endl;
        //cout<< "TAC for each area(with error): " << realTAC<< "  " <<endl; 
 
        for(int n=1;n<=numPop;n++){      
         for(int k=1;k<=numPop;k++){   
          for(int m=1;m<=numAge;m++){      
           subMixPop(i,j,n,k,m)=tmpPop(n,m)*moveOp(j,n,k); 
           mixBegYr(i,j,k,m)+=subMixPop(i,j,n,k,m);     
           tmpPopMove(k,m) = mixBegYr(i,j,k,m);     
          }       
         }
        }         
        
        NewtonRaphsonForF(tmpPopMove); // get new F based on annualTAC and using the abundance at two year before                  
        cout<<" updated F: "<< newF <<  "  " <<endl;
        
        for(int n=1;n<=numPop;n++){      
          for(int k=1;k<=numPop;k++){   
            for(int m=1;m<=numAge;m++){ 
              F(i,j,k) = newF(k);  
              subPop(i,j,n,k,m)=tmpPop(n,m)*moveOp(j,n,k)*mfexp(-(natMort(k)+select(m)*newF(k)));
              subCatch(i,j,n,k,m)=(select(m)*newF(k)/(natMort(k)+select(m)*newF(k)))*tmpPop(n,m)*moveOp(j,n,k)*(1.-mfexp(-(natMort(k)+select(m)*newF(k))));
              double tmp=subPop(i,j,n,k,m)*eggsKg*propFemale*maturity(m)*weight(m);  
              subSSB(i,j,n,k)+= tmp; 
              subSpawn(i,j,n,k)+=(tmp/eggsKg)/weight(m);
              pop(i,j,n,m)+=subPop(i,j,n,k,m); 
              mixStock(i,j,k,m)+=subPop(i,j,n,k,m);           
              areaCatch(i,j,k,m)+=subCatch(i,j,n,k,m); 
            } // end age m
            SSB(i,j,n)+=subSSB(i,j,n,k);  
            Spawn(i,j,n)+=subSpawn(i,j,n,k);
          }   // end area k
        }     // end pop n

     //added by Yang on Jan11,2017, for chapter4 population composition of catch 
      for(int k=1;k<=numPop;k++)
      	for(int n=1;n<=numPop;n++)
          for (int m=assessFAge;m<=assessLAge;m++)
             areaCatchPop(i,j,k,n,m)=subCatch(i,j,n,k,m);

        for(int k=1;k<=numPop;k++)  {
          effort(i,j,k)=1./q*newF(k)*mfexp((sdEffort(k)*randn(rnd))-(0.5*square(sdEffort(k))));      // fishing effort:  E= F*exp(sd-0.5*sigma^2) /q
          for(int m=1;m<=numAge;m++)  yieldTolarea(i,j,k) += areaCatch(i,j,k,m)*weight(m); 
        }
     }         // end assessment if branch: j > assessYr
 
 //****************************************************
 //population data simulation ends here
 //****************************************************
  

  //run accessment to get TAC from assessment model
  //the order for following if branch does matters, should not change it 
  //*****assessmentOK ::assessOK+=assessStepYr 
  //**********************assessment on each assessment OK year***********************

      if(j == assessOK){ //start on assessYr, then every year after 
       //*****************independent population model begins here***************
        
       if((model==3)|(model==4))  {//meta
       //check if std file exist or not, if true, then delete it    
        if(model==4){     
        ifstream fexist("sim4assessment4.std");
        if(fexist) {
         fexist.close(); 
         if(forLinux==0) int tmp11=system("del sim4assessment4.std");
         else int tmp11=system("rm -f sim4assessment4.std");
        }
        }
	if(model==3){     
        ifstream fexist("sim4assessment3.std");
        if(fexist) {
         fexist.close(); 
         if(forLinux==0) int tmp11=system("del sim4assessment3.std");
         else int tmp11=system("rm -f sim4assessment3.std");
        }
        }
        totAssess++;
         //cal true mean Rec AND initial Abun, used as starting value for assessment model, to availed hessian warning.
         dvector IntPop(1,numPop);
         dvector RecPop(1,numPop);
         for(int k=1;k<=numPop;k++){ 
          //cout<<endl<<"#==== N  (yr x age) "<<k <<" ===="<<endl; //output catch at age     
          if(useAllYr==1){
            IntPop(k)=mean(popBegYr(i,k,1)(assessFAge+1,assessLAge));
            for(int yr=1;yr<=j-1;yr++) { //use all previous yr data
             dvector TmpPop(1,j-1);
             TmpPop(yr)=popBegYr(i,k,yr,assessFAge);
              //cout<<popBegYr(i,k,yr)<<endl;  //all ages for this yr
             if(yr==j-1)  RecPop(k)=mean(TmpPop);
            }	   
          }else{                        //output fixed num yr data
            int startYr=j-nYrAssess;
            IntPop(k)=mean(popBegYr(i,k,startYr)(assessFAge+1,assessLAge));
            for(int yr=startYr;yr<=j-1;yr++){
              dvector TmpPop(startYr,j-1);
              TmpPop(yr)=popBegYr(i,k,yr,assessFAge);
              //cout<<popBegYr(i,k,yr)<<endl;
              if(yr==j-1) RecPop(k)=mean(TmpPop);
            }
           }                                
         }
         IntPop=log(IntPop);
         RecPop=log(RecPop);
         //cout<<"Ini Mean"<<IntPop<<endl;
         //cout<<"Rec Mean"<<RecPop<<endl;
         //exit(11);


         //catch at age data for first three years would be exactly the same (except last age group) because recruitment in assessment model start from age 3, so first three year age 1 recruitment did not influence the catch at age data for first three years

             //output data file for assessment
          if(model==3) {
          ofstream report("sim4assessment3.dat");
          report<<"# sim,yr "<<i<<", "<<j-1<<endl; //atucally yr data from
          if(nYrAssess>assessYr-1) nYrAssess=assessYr-1;           //in case user input big number than actual assessment starting yr
          report<<"# number of Pops"<<endl;
          report<<numPop<<endl;
          report<<"#intermixing matrix"<<endl;
          if (MixVLev==2) report<<moveParVar(j-1)<<endl;
          else report<<moveAss<<endl;
          report<<"# num years, first Age, last age "<<endl;
          int tmpYr=(useAllYr==1)? j-1:nYrAssess;
          report<<tmpYr <<" "<<assessFAge<<" "<<assessLAge<<endl;
          report<<endl<<"# natural mort, total target mort "<<endl;
          report<<natMort(1) <<" "<<targetTotMort(1)<<endl;
          report<<endl<<"# weight at age "<<endl;
          report<<weight(assessFAge,assessLAge)<<endl;
          report<<endl<<"# selectivity at age "<<endl;
          report<<select(assessFAge,assessLAge)<<endl; 
          report<<endl<<"# maturity at age "<<endl;
          report<<maturity(assessFAge,assessLAge)<<endl;
          report<<endl<<"# female proportion "<<endl;
          report<<propFemale<<endl;     
          report<<endl<<"#Initial"<<endl;
          report<<IntPop<<endl;
          report<<endl<<"#Rec"<<endl;
          report<<RecPop<<endl;                          
          report<<endl<<"#==== catch at age for area 1  (yr x age)  ===="<<endl; //output catch at age     
          if(useAllYr==1){
            for(int yr=1;yr<=j-1;yr++)  //use all previous yr data
              report<<areaCatch(i,yr,1)(assessFAge,assessLAge)<<endl;  //all ages for this yr
          }else{                        //output fixed num yr data
            int startYr=j-nYrAssess;
            for(int yr=startYr;yr<=j-1;yr++)
              report<<areaCatch(i,yr,1)(assessFAge,assessLAge)<<endl;
          }   
          report<<endl<<"#==== catch at age for area 2  (yr x age)  ===="<<endl; //output catch at age     
          if(useAllYr==1){
            for(int yr=1;yr<=j-1;yr++)  //use all previous yr data
              report<<areaCatch(i,yr,2)(assessFAge,assessLAge)<<endl;  //all ages for this yr
          }else{                        //output fixed num yr data
            int startYr=j-nYrAssess;
            for(int yr=startYr;yr<=j-1;yr++)
              report<<areaCatch(i,yr,2)(assessFAge,assessLAge)<<endl;
          }
          report<<endl<<"#==== catch at age for area 3  (yr x age)  ===="<<endl; //output catch at age     
          if(useAllYr==1){
            for(int yr=1;yr<=j-1;yr++)  //use all previous yr data
              report<<areaCatch(i,yr,3)(assessFAge,assessLAge)<<endl;  //all ages for this yr
          }else{                        //output fixed num yr data
            int startYr=j-nYrAssess;
            for(int yr=startYr;yr<=j-1;yr++)
              report<<areaCatch(i,yr,3)(assessFAge,assessLAge)<<endl;
          }
          report<<endl<<"#==== catch at age for area 4  (yr x age)  ===="<<endl; //output catch at age     
          if(useAllYr==1){
            for(int yr=1;yr<=j-1;yr++)  //use all previous yr data
              report<<areaCatch(i,yr,4)(assessFAge,assessLAge)<<endl;  //all ages for this yr
          }else{                        //output fixed num yr data
            int startYr=j-nYrAssess;
            for(int yr=startYr;yr<=j-1;yr++)
              report<<areaCatch(i,yr,4)(assessFAge,assessLAge)<<endl;
          }                              //output fishing effort
          
          report<<endl<<"#==== fishing effort for area 1 (year)  ===="<<endl; 
          if(useAllYr==1){
            for(int yr=1;yr<=j-1;yr++)  //use all previous yr data
              report<<effort(i,yr,1) <<endl;
          }else{                        //output fixed num yr data
            int startYr=j-nYrAssess;
            for(int yr=startYr;yr<=j-1;yr++)
              report<<effort(i,yr,1)<<endl;
          }  
          report<<endl<<"#==== fishing effort for area 2 (year)  ===="<<endl; 
          if(useAllYr==1){
            for(int yr=1;yr<=j-1;yr++)  //use all previous yr data
              report<<effort(i,yr,2) <<endl;
          }else{                        //output fixed num yr data
            int startYr=j-nYrAssess;
            for(int yr=startYr;yr<=j-1;yr++)
              report<<effort(i,yr,2)<<endl;
          } 
          report<<endl<<"#==== fishing effort for area 3 (year)  ===="<<endl; 
          if(useAllYr==1){
            for(int yr=1;yr<=j-1;yr++)  //use all previous yr data
              report<<effort(i,yr,3) <<endl;
          }else{                        //output fixed num yr data
            int startYr=j-nYrAssess;
            for(int yr=startYr;yr<=j-1;yr++)
              report<<effort(i,yr,3)<<endl;
          } 
          report<<endl<<"#==== fishing effort for area 4 (year)  ===="<<endl; 
          if(useAllYr==1){
            for(int yr=1;yr<=j-1;yr++)  //use all previous yr data
              report<<effort(i,yr,4) <<endl;
          }else{                        //output fixed num yr data
            int startYr=j-nYrAssess;
            for(int yr=startYr;yr<=j-1;yr++)
              report<<effort(i,yr,4)<<endl;
          }                               
          report<<endl<<"# number of last yrs mortality  as mean mortality for TAC "<<endl<<lastNumYrMort<<endl;
          report<<endl<<"# number of last yrs on first age as mean recruitment for TAC "<<endl<<lastNumYrRec<<endl;
          report<<endl<<"# lognormal std error for observed catch for each area "<<endl<<stdObsCatch(1)<<endl;
          report<<endl<<"# effective sample size "<<endl<<effSampleSize <<endl;
          report<<endl<<"# testing "<<endl;
          report<<"-9999 -8888 -7777"<<endl;
          report.close();
          //exit(4);



          if(forLinux==0) int tmp12=system("sim4assessment3 -nox > NUL"); //run assessment program
          else int tmp12=system("./sim4assessment3 -nox > /dev/null");    //run assessment program
          
          //read in TAC from assessment report file
          ifstream readin("sim4assessment3.rep");  
 	  readin>>readInTAC(1);				   
          readin>>readInTAC(2);
          readin>>readInTAC(3);
          readin>>readInTAC(4);
          dvariable grad;
          readin>>grad;
          //Grad(i,k,j)= double(grad);
          //readin>>catch_mre(i,1,j);
          //readin>>catch_mre(i,2,j);
          //readin>>catch_mre(i,3,j);
          //readin>>catch_mre(i,4,j);
          //readin>>catch_mare(i,1,j);
          //readin>>catch_mare(i,2,j);
          //readin>>catch_mare(i,3,j);
          //readin>>catch_mare(i,4,j);
          readin>>assess_ssb(i,1,j);
          readin>>assess_ssb(i,2,j);
          readin>>assess_ssb(i,3,j);
          readin>>assess_ssb(i,4,j);
          readin>>assess_rec(i,1,j);
 	  readin>>assess_rec(i,2,j);
 	  readin>>assess_rec(i,3,j);
 	  readin>>assess_rec(i,4,j);
          readin>>assess_F(i,1,j);
 	  readin>>assess_F(i,2,j);
 	  readin>>assess_F(i,3,j);
 	  readin>>assess_F(i,4,j);


          //*******here to check
          //cout<<"sim="<<i<<" yr="<<j<<" area="<<k<<"Fnew"<<newF(k)<<" tac="<<readInTAC(k) <<" N="<<catch_N(i,k,j)<<" mre="<<catch_mre(i,k,j)<<" mare="<<catch_mare(i,k,j)<<"grad"<<grad<<endl;
          readin.close();
           //if max gradient too big
           if(fabs(grad)>0.004) {
            cout<<"****** Warning: Proper convergence could not be reached, gradient >0.004 !"<<endl;
            //exit(5);
            for(int k=1;k<=numPop;k++){ 
             warnConverg(k)++;  //convergence warning counter for each area
             indicator(i,k)=1.0;
             counter(i,k)+=1.0;
             if(annualTAC(i,j,k)!=0) readInTAC(k)=annualTAC(i,j,k);
             }
           }
		       }


          if(model==4) {
	  	       ofstream report("sim4assessment4.dat");
		       report<<"# sim,yr "<<i<<", "<<j-1<<endl; //atucally yr data from
          	       if(nYrAssess>assessYr-1) nYrAssess=assessYr-1;           //in case user input big number than actual assessment starting yr
          	       report<<"# number of Pops"<<endl;
          	       report<<numPop<<endl;
          	       report<<"#intermixing matrix"<<endl;
          	       if (MixVLev==2) report<<moveParVar(j-1)<<endl;
          	       else report<<moveAss<<endl;
          	       report<<"# num years, first Age, last age "<<endl;
          	       int tmpYr=(useAllYr==1)? j-1:nYrAssess;
          	       report<<tmpYr <<" "<<assessFAge<<" "<<assessLAge<<endl;
          	       report<<endl<<"# natural mort, total target mort "<<endl;
          	       report<<natMort(1) <<" "<<targetTotMort(1)<<endl;
          	       report<<endl<<"# weight at age "<<endl;
          	       report<<weight(assessFAge,assessLAge)<<endl;
          	       report<<endl<<"# selectivity at age "<<endl;
          	       report<<select(assessFAge,assessLAge)<<endl; 
          	       report<<endl<<"# maturity at age "<<endl;
          	       report<<maturity(assessFAge,assessLAge)<<endl;
          	       report<<endl<<"# female proportion "<<endl;
          	       report<<propFemale<<endl;     
          	       report<<endl<<"#Initial"<<endl;
          	       report<<IntPop<<endl;
          	       report<<endl<<"#Rec"<<endl;
          	       report<<RecPop<<endl;                          
          	       report<<endl<<"#==== fishing effort for area 1 (year)  ===="<<endl; 
          	       if(useAllYr==1){
            	        for(int yr=1;yr<=j-1;yr++)  //use all previous yr data
              		report<<effort(i,yr,1) <<endl;
          		}else{                        //output fixed num yr data
            		int startYr=j-nYrAssess;
            		for(int yr=startYr;yr<=j-1;yr++)
              		report<<effort(i,yr,1)<<endl;
          		}  
          		report<<endl<<"#==== fishing effort for area 2 (year)  ===="<<endl; 
          		if(useAllYr==1){
            		 for(int yr=1;yr<=j-1;yr++)  //use all previous yr data
              		 report<<effort(i,yr,2) <<endl;
          		}else{                        //output fixed num yr data
            		int startYr=j-nYrAssess;
            		for(int yr=startYr;yr<=j-1;yr++)
			report<<effort(i,yr,2)<<endl;
          		} 
          		report<<endl<<"#==== fishing effort for area 3 (year)  ===="<<endl; 
          if(useAllYr==1){
            for(int yr=1;yr<=j-1;yr++)  //use all previous yr data
              report<<effort(i,yr,3) <<endl;
          }else{                        //output fixed num yr data
            int startYr=j-nYrAssess;
            for(int yr=startYr;yr<=j-1;yr++)
              report<<effort(i,yr,3)<<endl;
          } 
          report<<endl<<"#==== fishing effort for area 4 (year)  ===="<<endl; 
          if(useAllYr==1){
            for(int yr=1;yr<=j-1;yr++)  //use all previous yr data
              report<<effort(i,yr,4) <<endl;
          }else{                        //output fixed num yr data
            int startYr=j-nYrAssess;
            for(int yr=startYr;yr<=j-1;yr++)
              report<<effort(i,yr,4)<<endl;
          }    
          //added by Yang on Jan11,2017, for chapter4 population composition of catch 
 	  report<<endl<<"#==== Catch at pop for area 1  (yr  x pop)  ===="<<endl; 
          if(useAllYr==1){
            for(int yr=1;yr<=j-1;yr++)  //use all previous yr data
              report<<areaCatchPop(i,yr,1)<<endl;  //all ages for this yr
          }else{                        //output fixed num yr data
            int startYr=j-nYrAssess;
            for(int yr=startYr;yr<=j-1;yr++)
              report<<areaCatchPop(i,yr,1)<<endl;
          }  
          report<<endl<<"#==== Catch at pop and age for area 2 (yr  x pop)  ===="<<endl; 
          if(useAllYr==1){
            for(int yr=1;yr<=j-1;yr++)  //use all previous yr data
              report<<areaCatchPop(i,yr,2)<<endl;  //all ages for this yr
          }else{                        //output fixed num yr data
            int startYr=j-nYrAssess;
            for(int yr=startYr;yr<=j-1;yr++)
              report<<areaCatchPop(i,yr,2)<<endl;
          }  
  	  report<<endl<<"#==== Catch at pop and age for area 3  (yr  x pop)  ===="<<endl; 
          if(useAllYr==1){
            for(int yr=1;yr<=j-1;yr++)  //use all previous yr data
              report<<areaCatchPop(i,yr,3)<<endl;  //all ages for this yr
          }else{                        //output fixed num yr data
            int startYr=j-nYrAssess;
            for(int yr=startYr;yr<=j-1;yr++)
              report<<areaCatchPop(i,yr,3)<<endl;
          }  
 	  report<<endl<<"#==== Catch at pop and age for area 4 (yr  x pop)  ===="<<endl; 
          if(useAllYr==1){
            for(int yr=1;yr<=j-1;yr++)  //use all previous yr data
              report<<areaCatchPop(i,yr,4)<<endl;  //all ages for this yr
          }else{                        //output fixed num yr data
            int startYr=j-nYrAssess;
            for(int yr=startYr;yr<=j-1;yr++)
              report<<areaCatchPop(i,yr,4)<<endl;
          }      
             
          report<<endl<<"# number of last yrs mortality  as mean mortality for TAC "<<endl<<lastNumYrMort<<endl;
          report<<endl<<"# number of last yrs on first age as mean recruitment for TAC "<<endl<<lastNumYrRec<<endl;
          report<<endl<<"# lognormal std error for observed catch for each area "<<endl<<stdObsCatch(1)<<endl;
          report<<endl<<"# effective sample size "<<endl<<effSampleSize <<endl;
          report<<endl<<"# testing "<<endl;
          report<<"-9999 -8888 -7777"<<endl;
          report.close();
       
	          if(forLinux==0) int tmp12=system("sim4assessment4 -nox > NUL"); //run assessment program
          else int tmp12=system("./sim4assessment4 -nox > /dev/null");    //run assessment program
          
          //read in TAC from assessment report file
          ifstream readin("sim4assessment4.rep");  
 	  readin>>readInTAC(1);
          readin>>readInTAC(2);
          readin>>readInTAC(3);
          readin>>readInTAC(4);
          dvariable grad;
          readin>>grad;
          //Grad(i,k,j)= double(grad);
          //readin>>catch_mre(i,1,j);
          //readin>>catch_mre(i,2,j);
          //readin>>catch_mre(i,3,j);
          //readin>>catch_mre(i,4,j);
          //readin>>catch_mare(i,1,j);
          //readin>>catch_mare(i,2,j);
          //readin>>catch_mare(i,3,j);
          //readin>>catch_mare(i,4,j);
          readin>>assess_ssb(i,1,j);
          readin>>assess_ssb(i,2,j);
          readin>>assess_ssb(i,3,j);
          readin>>assess_ssb(i,4,j);
          readin>>assess_rec(i,1,j);
 	  readin>>assess_rec(i,2,j);
 	  readin>>assess_rec(i,3,j);
 	  readin>>assess_rec(i,4,j);
          readin>>assess_F(i,1,j);
 	  readin>>assess_F(i,2,j);
 	  readin>>assess_F(i,3,j);
 	  readin>>assess_F(i,4,j);
          //*******here to check
         
          readin.close();
           //if max gradient too big
           if(fabs(grad)>0.004) {
            cout<<"****** Warning: Proper convergence could not be reached, gradient >0.004 !"<<endl;
            //exit(5);
            for(int k=1;k<=numPop;k++){ 
             warnConverg(k)++;  //convergence warning counter for each area
             indicator(i,k)=1.0;
             counter(i,k)+=1.0;
             if(annualTAC(i,j,k)!=0) readInTAC(k)=annualTAC(i,j,k);
 	    
             }
           }
 	  
           }
           /*
            for(int k=1;k<=numPop;k++) {
 	   	   cout<<"sim="<<i<<" yr="<<j<<" area="<<k<<"Fnew"<<F(i,j-1,k)<<" tac="<<readInTAC(k) <<" assess_F="<<assess_F(i,k,j)<<endl;
    		    cout<<" assess_rec="<<assess_rec(i,k,j)<<endl;
     		    // cout<<" assess_ssb="<<assess_ssb(i,k,j)<<endl;
		   }
		   */
           //exit(4);
	   

           if (model==3){
              //if got hessian warning
              ifstream ifs("sim4assessment3.std");
              if(ifs) {
              	      ifs.close();
            	      if(forLinux==0) int tmp14=system("del sim4assessment3.std");
            	      else int tmp14=system("rm -f sim4assessment3.std");
           	      }
              else { 
              	   cout<<"******************************** Warning: The function maximizer failed !" <<endl;
            	   //exit(5);
             	   for(int k=1;k<=numPop;k++){ 
            	   	   warnHessian(k)++;  //convergence warning counter for each area
            	   	   hess_ind(i,k)=1.0;
            	  	    hess_count(i,k)+=1.0;
            	   }
              }
              if(forLinux==0) int tmp13=system("del sim4assessmen3.rep");
              else int tmp13=system("rm -f sim4assessment3.rep");
	   }
 	   if (model==4){
              //if got hessian warning
              ifstream ifs("sim4assessment4.std");
              if(ifs) {
              	      ifs.close();
            	      if(forLinux==0) int tmp14=system("del sim4assessment4.std");
            	      else int tmp14=system("rm -f sim4assessment4.std");
              }
              else { 
              	   cout<<"******************************** Warning: The function maximizer failed !" <<endl;
            	   //exit(5);
             	   for(int k=1;k<=numPop;k++){ 
            	   warnHessian(k)++;  //convergence warning counter for each area
            	   hess_ind(i,k)=1.0;
            	   hess_count(i,k)+=1.0;
            	   }
             }
             if(forLinux==0) int tmp13=system("del sim4assessmen4.rep");
             else int tmp13=system("rm -f sim4assessment4.rep");
	   }

          for(int k=1;k<=numPop;k++){ 
	  	   for(int m=assessFAge;m<=assessLAge;m++) {
            	   real_ssb(i,k,j) +=   popBegYr(i,k,j-1,m)*weight(m)*propFemale*maturity(m);   //number of fish * weight~~~SSB real                      
           	   } //end of m loop  
           	  for(int x=1;x<=nYrAssess;x++){
	   	   int ind=x+j-nYrAssess-1;
  		   real_rec(i,k,j,x)=Recruit(i,ind,k);
		   real_F(i,k,j,x)= F(i,ind,k);
                   //cout<<"ind:"<<ind<<" k:"<<k<<endl;
 		   //cout<<"F:"<<F(i,ind,k)<<endl;
	   }
           //end of m loop                                 
          } //end k area loop        
          for(int k=1;k<=numPop;k++) {
           //pooled_real_ssb(i,j)   +=   real_ssb(i,k,j);
           //pooled_assess_ssb(i,j) +=   assess_ssb(i,k,j);  //number of fish * weight~~~SSB assessment result           
           //pooled_real_rec(i,j) +=  Recruit(i,j-1,k);
           //pooled_assess_rec(i,j) += assess_rec(i,k,j);  
           ssb_mre(i,k,j)=(assess_ssb(i,k,j)- real_ssb(i,k,j))/(real_ssb(i,k,j)+eps); //re for terminal year SSB
 	   ssb_mare(i,k,j)=abs(assess_ssb(i,k,j)- real_ssb(i,k,j))/(real_ssb(i,k,j)+eps); //are for terminal year SSB
           for(int x=1;x<=nYrAssess;x++) {
   	   	   rec_mre(i,k,j,x)=(assess_rec(i,k,j,x)- real_rec(i,k,j,x))/(real_rec(i,k,j,x)+eps); //re for terminal year rec
        
 	   	   //rec_mare(i,k,j,x)=abs(assess_rec(i,k,j,x)- real_rec(i,k,j,x))/(real_rec(i,k,j,x)+eps); //are for terminal year rec
          	   F_mre(i,k,j,x)=(assess_F(i,k,j,x)- real_F(i,k,j,x))/(real_F(i,k,j,x)+eps); //re for terminal year rec

		   }       
		   //cout<<" assess_F="<<assess_F(i,k,j)<<endl;
 		   //cout<<" real_F="<<real_F(i,k,j)<<endl;
                   //cout<<" F_mre="<<F_mre(i,k,j)<<endl;

          }
         //exit(5);
          //pooled_ssb_mre(i,j)=(pooled_assess_ssb(i,j)- pooled_real_ssb(i,j))/(pooled_real_ssb(i,j)+eps); //mre of SSB
          //pooled_ssb_mare(i,j)=abs(pooled_assess_ssb(i,j)- pooled_real_ssb(i,j))/(pooled_real_ssb(i,j)+eps); //mare of SSB
          //pooled_rec_mre(i,j)=(pooled_assess_rec(i,j)- pooled_real_rec(i,j))/(pooled_real_rec(i,j)+eps); //re for terminal year rec
          //pooled_rec_mare(i,j)=abs(pooled_assess_rec(i,j)- pooled_real_rec(i,j))/(pooled_real_rec(i,j)+eps); //are for terminal year rec
         }  //end if model=meta 3
     
       
        assessOK+=assessStepYr; //update to next assess yr     
     }
     //************************end of each time assessment****************************** 
    } //end year loop j
  }//end sim loop i

//**********************************end of data generation*************************************
//*********************************************************************************************
  for(int i=1;i<=numSim;i++){
   for(int j=1;j<=numYear;j++){
    for(int k=1;k<=numPop;k++) {
     for(int m=1;m<=numAge;m++){ 
      mixStock1(i,k,j,m)=mixStock(i,j,k,m);
      pop1(i,k,j,m)=pop(i,j,k,m);
      areaCatch1(i,k,j,m)=areaCatch(i,j,k,m);
      areaYield1(i,k,j,m)=areaCatch1(i,k,j,m)*weight(m);
      annualYield(i,k,j)+=areaYield1(i,k,j,m);
     }    
    }
   }
  }


//**********************************result analysis********************************************* 
 if(model!=2){
  Saveindependent();
  }
 //if(model==2){
 // Savepooled();
 //} 
  SaveSimFile();
  exit(8);
 
 
 
 
PROCEDURE_SECTION
  exit(9);

FUNCTION int mixParameter(const dvector&stayR)
   dvector stayother(1,numPop);
   for(int i=1; i<=numPop; i++){
    stayother(i)= sum(stayR)-stayR(i);
   }
   for(int i=1; i<=numPop; i++){
    for (int j=1; j<=numPop; j++){
     if (j==i)  movePartmp(i,j)=stayR(i);  
     if (j!=i)  movePartmp(i,j)=(1-stayR(i))* stayR(j)/stayother(i);
     if (movePartmp(i,j)<1e-10) movePartmp(i,j)=0;
    }
   } 
  return 0;

// run newton method to estiamte F, based on real TAC
//popN is pop initial abundance at start of yr, with 1st age as recruitment,
// other ages from the N at end of previous yr but excluding the last age(age move up)
FUNCTION int NewtonRaphsonForF(const dmatrix&stockN)
  double criteria=0.0001;
  int maxLoop=100;
  for(int k=1;k<=numPop;k++){ //for each area
    double Fnew=0.1; //starting value, initial as .1
    double holder=0;  //hold Fnew value from last run
    double delta,tmpCatch,tmpCatchH,tmpCatchL;    
    for(int l=1;l<=maxLoop;l++){  //fixed the loop number as 10
      delta=Fnew*0.001; //try make the step small
      tmpCatch=0; tmpCatchH=0;tmpCatchL=0;
      for(int m=1;m<=numAge;m++){ //loop for all ages
       tmpCatch+=(select(m)*Fnew/(natMort(k)+select(m)*Fnew))*stockN(k,m)*(1.-mfexp(-(natMort(k)+select(m)*Fnew)));//use Fnew, fill in catch equation
       tmpCatchH+=(select(m)*(Fnew+delta)/(natMort(k)+select(m)*(Fnew+delta)))*stockN(k,m)*(1.-mfexp(-(natMort(k)+select(m)*(Fnew+delta))));//use Fnew+delta
       tmpCatchL+=(select(m)*(Fnew-delta)/(natMort(k)+select(m)*(Fnew-delta)))*stockN(k,m)*(1.-mfexp(-(natMort(k)+select(m)*(Fnew-delta))));//use Fnew-delta
      }
      //cout<<k<<" "<<l<<"     "<<Fnew<<"     "<<tmpCatch<<"      "<<delta<<endl;
      //following Fratio=f(x0)/f'(x0)
      double Fratio=(tmpCatch-realTAC(k))/((tmpCatchH-tmpCatchL)/(2.*delta));
      Fnew=Fnew-Fratio;  //update Fnew
      if(Fnew<0) Fnew=0.5*(Fnew+Fratio);  //adjust it
      if(Fnew>3) Fnew=3;//At low N, Fnew would sometimes be really high, so limit here
      if(fabs(Fnew-holder)<criteria) break; //if meeting the convergence criteria, then exit
      holder=Fnew;
    }
    newF(k)=Fnew;
  }
  //cout<<endl<<newF<<endl;
  return 0;
 

FUNCTION Saveindependent
  ofstream report(outPerform1);
  report<<"#============= Number of populations ==================" << endl;
  report << "# numPop" << endl; 
  report<<numPop<<endl; 
  report<<"#=============== Number of iterations =======================" << endl;
  report << "# numSim" << endl; 
  report<<numSim<<endl; 
  
  report<<"#================= MRE for pop SSB=================="<<endl;
  for(int i=1;i<=numSim;i++)
  {
    report<<"## Simulation = "<<i<<endl;
    for(int n=1;n<=numPop;n++)
    {
      report<<"### pop_ssb_mre(Year x population) = "<<n<<endl;
      report<<"### Sim"<<i<<"_popssb_mre_pop"<<n<<endl; //that is the name when read in R
      report<<setfixed()<<setprecision(4)<<setw(10)<<ssb_mre(i,n)(assessYearVec)<<endl<<endl;
    }
   }
  report<<"#================= MARE for pop SSB=================="<<endl;
  for(int i=1;i<=numSim;i++)
  {
    report<<"## Simulation = "<<i<<endl;
    for(int n=1;n<=numPop;n++)
    {
      report<<"### pop_ssb_mare(Year x population) = "<<n<<endl;
      report<<"### Sim"<<i<<"_popssb_mare_pop"<<n<<endl; //that is the name when read in R
      report<<setfixed()<<setprecision(4)<<setw(10)<<ssb_mare(i,n)(assessYearVec)<<endl<<endl;
    }
   }
  report<<"#================= MRE for pop rec=================="<<endl;
  for(int i=1;i<=numSim;i++)
  {
    report<<"## Simulation = "<<i<<endl;
    for(int n=1;n<=numPop;n++)
    {
      //for(int j=assessYr;j<=numYear;j++){
      report<<"### pop_rec_mre(Year x population) = "<<n<<endl;
      report<<"### Sim"<<i<<"_poprec_mre_pop"<<n<<endl; //that is the name when read in R
      report<<setfixed()<<setprecision(4)<<setw(10)<<rec_mre(i,n)<<endl<<endl;

      //}
    }
   }

  report<<"#================= MRE for pop F=================="<<endl;
  for(int i=1;i<=numSim;i++)
  {
    report<<"## Simulation = "<<i<<endl;
    for(int n=1;n<=numPop;n++)
    {
      //for(int j=assessYr;j<=numYear;j++){
      report<<"### pop_F_mre(Year x population) = "<<n<<endl;
      report<<"### Sim"<<i<<"_popF_mre_pop"<<n<<endl; //that is the name when read in R
      report<<setfixed()<<setprecision(4)<<setw(10)<<F_mre(i,n)<<endl<<endl;
      //}
    }
   }
 
  /*
   report<<"#================= MARE for pop rec=================="<<endl;
  for(int i=1;i<=numSim;i++)
  {
    report<<"## Simulation = "<<i<<endl;
    for(int n=1;n<=numPop;n++)
    {
      //for(int j=assessYr;j<=numYear;j++){
      report<<"### pop_rec_mare(Year x population) = "<<n<<endl;
      report<<"### Sim"<<i<<"_poprec_mare_pop"<<n<<endl; //that is the name when read in R
      report<<setfixed()<<setprecision(4)<<setw(10)<<rec_mare(i,n,numYear)<<endl<<endl;
      //}
    }
   }
  */
  report<<"#============== Convergence Gradient Indicator (sim x area) ===="<<endl;
  report <<"# assess_conv_warn" << endl;
  report<<setfixed()<<setprecision(0)<<setw(10)<<indicator<<endl<<endl;
  report<<"#============== Hession Warning Indicator (sim x area) ===="<<endl;
  report <<"# assess_hess_warn" << endl;
  report<<setfixed()<<setprecision(0)<<setw(10)<<hess_ind<<endl<<endl;
  report<<"#============== Convergence Gradient Counter (sim x area) ===="<<endl;
  report <<"# counter_conv_warn" << endl;
  report<<setfixed()<<setprecision(0)<<setw(10)<<counter<<endl<<endl;
  report<<"#============== Hession Warning Counter (sim x area) ===="<<endl;
  report <<"# counter_hess_warn" << endl;
  report<<setfixed()<<setprecision(0)<<setw(10)<<hess_count<<endl<<endl;
  report<<"#============== Assessment Convergence Warnings Prop. for each area cross all sim. ===="<<endl;
  report<<"# convergWarn"<<endl; //that is the name when read in R
  report<<setfixed()<<setprecision(6)<<setw(10)<<warnConverg/(totAssess)<<endl<<endl;
  report<<"#============== Assessment Hessian  Warnings Prop. for each area cross all sim. ===="<<endl;
  report<<"# hessianWarn"<<endl; //that is the name when read in R
  report<<setfixed()<<setprecision(6)<<setw(10)<<warnHessian/(totAssess)<<endl<<endl; 
  /*
  report<<"#============== Grad(area x year) ===="<<endl;
  report <<"# Gradient" << endl;
  for(int i=1;i<=numSim;i++)
  {    
   report<<setfixed()<<setprecision(4)<<setw(10)<<Grad(i)<<" ";
   report<<endl;
  }
  */
  report.close();  
  
  
  
  
  
  /*  
FUNCTION Savepooled
  ofstream report(outPerform2);
  report<<"#============= ***For pooled population populations*** ==================" << endl;
  report<<"#============= *************************************** ==================" << endl;
  report<<"#=============== Number of iterations =======================" << endl;
  report << "# numSim" << endl; 
  report<<numSim<<endl;
  
  
  report<<"#======================= MRE for Pooled SSB (sim )  =========================="<<endl;  
  report<<"# Pooled_SSB_mre "<<endl; //that is the name when read in R
  for(int i=1;i<=numSim;i++)
  {    
   report<<"## Simulation ( Year)= "<<i<<endl;
   report<<"### Sim"<<i<<"_poolSSB_mre"<<endl; //that is the name when read in R
   report<<setfixed()<<setprecision(4)<<setw(10)<<pooled_ssb_mre(i)(assessYearVec)<<endl<<endl;
  }
  
  report<<"#======================= MARE for Pooled SSB (sim )  =========================="<<endl;  
  report<<"# Pooled_SSB_mare "<<endl; //that is the name when read in R
  for(int i=1;i<=numSim;i++)
  {    
   report<<"## Simulation ( Year)= "<<i<<endl;
   report<<"### Sim"<<i<<"_poolSSB_mare"<<endl; //that is the name when read in R
   report<<setfixed()<<setprecision(4)<<setw(10)<<pooled_ssb_mare(i)(assessYearVec)<<endl<<endl;
  } 

  report<<"#======================= MRE for Pooled rec (sim )  =========================="<<endl;  
  report<<"# Pooled_rec_mre "<<endl; //that is the name when read in R
  for(int i=1;i<=numSim;i++)
  {    
   report<<"## Simulation ( Year)= "<<i<<endl;
   report<<"### Sim"<<i<<"_poolRec_mre"<<endl; //that is the name when read in R
   report<<setfixed()<<setprecision(4)<<setw(10)<<pooled_rec_mre(i)(assessYearVec)<<endl<<endl;
  }
  
  report<<"#======================= MARE for Pooled rec (sim )  =========================="<<endl;  
  report<<"# Pooled_rec_mare "<<endl; //that is the name when read in R
  for(int i=1;i<=numSim;i++)
  {    
   report<<"## Simulation ( Year)= "<<i<<endl;
   report<<"### Sim"<<i<<"_poolRec_mare"<<endl; //that is the name when read in R
   report<<setfixed()<<setprecision(4)<<setw(10)<<pooled_rec_mare(i)(assessYearVec)<<endl<<endl;
  } 

  report<<"#============== Convergence Gradient Indicator for pooled population (sim ) ===="<<endl;
  report <<"# Pooled_assess_conv_warn" << endl;
  report<<setfixed()<<setprecision(0)<<setw(10)<<pooledindicator<<endl<<endl;
  report<<"#============== Hession Warning Indicator for pooled population(sim ) ===="<<endl;
  report <<"# Pooled_assess_hess_warn" << endl;
  report<<setfixed()<<setprecision(0)<<setw(10)<<pooledhess_ind<<endl<<endl;
  report<<"#============== Convergence Gradient Counter for pooled population(sim ) ===="<<endl;
  report <<"# Pooled_counter_conv_warn" << endl;
  report<<setfixed()<<setprecision(0)<<setw(10)<<pooledcounter<<endl<<endl;
  report<<"#============== Hession Warning Counter for pooled population(sim ) ===="<<endl;
  report <<"# Pooled_counter_hess_warn" << endl;
  report<<setfixed()<<setprecision(0)<<setw(10)<<pooledhess_count<<endl<<endl;
  report<<"#============== Assessment Convergence Warnings Prop. for pooled population cross all sim. ===="<<endl;
  report<<"# Pooled_convergWarn"<<endl; //that is the name when read in R
  report<<setfixed()<<setprecision(6)<<setw(10)<<pooledwarnConverg/totAssess<<endl<<endl;
  report<<"#============== Assessment Hessian  Warnings Prop. for pooled population cross all sim. ===="<<endl;
  report<<"# Pooled_hessianWarn"<<endl; //that is the name when read in R
  report<<setfixed()<<setprecision(6)<<setw(10)<<pooledwarnHessian/totAssess<<endl<<endl;
  report.close();
  
  */

  
// save the simulation  results as file, filename from dat file
FUNCTION SaveSimFile
  ofstream report(outFile);
  report<<"# assessment model"<<endl; //for plot in r
  report<<model<<endl<<endl;
  report<<"# Years"<<endl; //for plot in r
  report<<years<<endl<<endl;

  report<<endl<<"#====================== SSB(Year x Pop) ========================="<<endl;
  for(int i=1;i<=numSim;i++)
  {
    report<<"## Simulation = "<<i<<endl;
    report<<"## Sim"<<i<<"_SSB"<<endl; //that is the name when read in R
    report<<setfixed()<<setprecision(0)<<setw(15)<<SSB(i)<<endl<<endl;  
//    report<<"## Sim"<<i<<"_Spawners"<<endl; //that is the name when read in R
//    report<<setfixed()<<setprecision(0)<<setw(15)<<Spawn(i)<<endl<<endl; 
//    report<<"## Sim"<<i<<"_Recruits"<<endl; //that is the name when read in R
//    report<<setfixed()<<setprecision(0)<<setw(15)<<Recruit(i)<<endl<<endl; 
  }

  report<<"#==================== Fishing Mortality ====================="<<endl;
  for(int i=1;i<=numSim;i++)
  {
    report<<"## Simulation (Year x Area)= "<<i<<endl;
    report<<"### Sim"<<i<<"_F"<<endl; //that is the name when read in R
    report<<setfixed()<<setprecision(5)<<setw(10)<<F(i)<<endl<<endl;
  }
  /*
  report<<"#================= Area TAC=================="<<endl;
  for(int i=1;i<=numSim;i++)
  {
    report<<"## Simulation ( Area x Year)= "<<i<<endl;
    report<<"### Sim"<<i<<"_TAC"<<endl; //that is the name when read in R
    report<<setfixed()<<setprecision(5)<<setw(10)<<annualTAC(i)<<endl<<endl;
  }

  report<<"#================= Area Harvest Yield =================="<<endl;
  for(int i=1;i<=numSim;i++)
  {
    report<<"## Simulation = "<<i<<endl;
    for(int n=1;n<=numPop;n++)
    {
      report<<"### harvest yield (Year x Age) = "<<n<<endl;
      report<<"### Sim"<<i<<"_Yield_Area"<<n<<endl; //that is the name when read in R
      report<<setfixed()<<setprecision(2)<<setw(10)<<areaYield1(i,n)<<endl<<endl;
    }
  }
  */ 
  report<<"#================= Area Yield=================="<<endl;
  for(int i=1;i<=numSim;i++)
  {
    report<<"## Simulation ( Year x Area)= "<<i<<endl;
    report<<"### Sim"<<i<<"_Yield"<<endl; //that is the name when read in R
    report<<setfixed()<<setprecision(0)<<setw(10)<<yieldTolarea(i)<<endl<<endl;
  }


  /*
  report<<"#======================= Stocks Population Abundance (beginning of yr) =========================="<<endl;
  for(int i=1;i<=numSim;i++)
  {
    report<<"## Simulation = "<<i<<endl;
    for(int n=1;n<=numPop;n++)
    {
      report<<"### Pop num (Year x Age) = "<<n<<endl;
      report<<"### Sim"<<i<<"_PopBeg_Pop"<<n<<endl; //that is the name when read in R
      report<<setfixed()<<setprecision(0)<<setw(10)<<popBegYr(i,n)<<endl<<endl;
    }
  }
  */
  report.close();
