//note: this is not a valid configuration file, but a configuration file template.

//variables to replace: MUHAT_B, SIGMA_B, N_ON

parameters = ("mu_s", "mu_b");

observables = {
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

gaussmean = {
    
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
      type = "product_distribution";
      distributions = ("@mu_s-dist", "@mu_b-dist");
   };
};

mu_s-dist = {
    type = "flat_distribution";
    mu_s = {
      range = (0.0, "inf");
      fix-sample-value = __N_ON__;
   }; 
};

mu_b-dist = {
    type = "gauss";
    parameter = "mu_b";
    mean = __MUHAT_B__;
    width = __SIGMA_B__;
    range = (0.0, "inf");
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
};

