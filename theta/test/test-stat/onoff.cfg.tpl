//note: this is not a valid configuration file, but a configuration file template.

//variables to replace: TAU, N_ON, N_OFF

parameters = ("mu_s", "mu_b", "tau");

observables = {
   background-region = {
      nbins = 1;
      range = [0.0, 1.0];
   };
   signal-region = {
      nbins = 1;
      range = [0.0, 1.0];
   };
};

flat-histo-s = {
    type = "fixed_poly";
    observable = "signal-region";
    normalize_to = 1.0;
    coefficients = [1.0];
};

flat-histo-b = {
    type = "fixed_poly";
    observable = "background-region";
    normalize_to = 1.0;
    coefficients = [1.0];
};

//model definition:
onoff = {
   // the background region has only one component: the background, with poisson expectation
   // mu_b times tau
   background-region = {
      background = {
          coefficient-function = { type = "multiply"; factors = ("mu_b", "tau"); };
          histogram = "@flat-histo-b";
      };
   };
   
   // signal region has two components: signal and background, with mu_s and mu_b as expected poisson means
   signal-region = {
       background = {
          coefficient-function = { type = "multiply"; factors = ("mu_b"); };
          histogram = "@flat-histo-s";
       };
       
       signal = {
          coefficient-function = { type = "multiply"; factors = ("mu_s"); };
          histogram = "@flat-histo-s"; 
       };
   };
   
   parameter-distribution = {
      type = "flat_distribution";
      
      mu_s = {
         range = (0.0, "inf");
         fix-sample-value = __N_OFF__;
      };
      
      mu_b = {
         range = (0.0, "inf");
         fix-sample-value = __N_OFF__;
      };
      
      //tau is fixed ...
      tau = {
        range = [__TAU__, __TAU__];  
      };
   };
};

data_source = {
    type = "histo_source";
    name = "data_source";
    signal-region = {
       type = "fixed_poly";
       observable = "signal-region";
       normalize_to = __N_ON__;
       coefficients = [1.0];
    };
    background-region = {
       type = "fixed_poly";
       observable = "background-region";
       normalize_to = __N_OFF__;
       coefficients = [1.0];
    };
};

