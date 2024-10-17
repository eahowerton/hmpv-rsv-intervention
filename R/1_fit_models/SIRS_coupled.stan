data {
  int N;                         // number of weeks
  int cases_mpv[N];              // HMPV detections 
  int cases_rsv[N];              // RSV detections
  real tests_mpv[N];             // moving average of tests scaled [0,1] 
  real tests_rsv[N];             // moving average of tests scaled [0,1] 
  int pop;                       // population size
  real gamma;                    // assumed recovery rate
  real mu;                       // assumed birth/death rate
  real omega;                    // waning rate (assumed same for both)
}

parameters {
  real<lower=0, upper=1> S0_mpv; // initial conditions
  real<lower=0, upper=1> I0_mpv;
  real<lower=0, upper=1> S0_rsv;
  real<lower=0, upper=1> I0_rsv;
  real<lower=0, upper=1> rho_mpv; // HMPV reporting rate
  real<lower=0, upper=1> rho_rsv; // RSV reporting rate
  real<lower=0> b;                // base transmission rate of RSV
  real<lower=0> r;                // relative transmission rate of HMPV vs RSV
  real<lower=-1> c;               // effect of RSV forcing
  real<lower=0, upper=1> a;       // amplitude of seasonal forcing
  real<lower=0, upper=1> p;       // phase of seasonal forcing    
  real<lower=0> phi_mpv;
  real<lower=0> phi_rsv;
}

transformed parameters {
  vector<lower=0>[N] S_mpv;
  vector<lower=0>[N] I_mpv;
  vector<lower=0>[N] Ipred_mpv;
  vector<lower=0>[N] R_mpv;
  vector<lower=0>[N] S_rsv;
  vector<lower=0>[N] I_rsv;
  vector<lower=0>[N] Ipred_rsv;
  vector<lower=0>[N] R_rsv;
  vector<lower=0>[N] beta_mpv;
  vector<lower=0>[N] beta_rsv;
  real<lower=0> foi_mpv;
  real<lower=0> foi_rsv;
  real<lower=0> Sdeath_mpv;
  real<lower=0> Ideath_mpv;
  real<lower=0> Rdeath_mpv;
  real<lower=0> Sout_mpv;
  real<lower=0> Iout_mpv;
  real<lower=0> Rout_mpv;
  real<lower=0> Sdeath_rsv;
  real<lower=0> Ideath_rsv;
  real<lower=0> Rdeath_rsv;
  real<lower=0> Sout_rsv;
  real<lower=0> Iout_rsv;
  real<lower=0> Rout_rsv;
  real<lower=0> birth;
  real<lower=0> max_RSV_inc;
  
  S_mpv[1] = S0_mpv * pop;
  I_mpv[1] = I0_mpv * pop;
  Ipred_mpv[1] = I0_mpv * pop * rho_mpv * tests_mpv[1];
  R_mpv[1] = (1-S0_mpv-I0_mpv) * pop;
  
  S_rsv[1] = S0_rsv * pop;
  I_rsv[1] = I0_rsv * pop;
  Ipred_rsv[1] = I0_rsv * pop * rho_rsv * tests_rsv[1];
  R_rsv[1] = (1-S0_rsv-I0_rsv) * pop;
  
  beta_rsv[1] = b;
  beta_mpv[1] = b*r;
  
  birth = (1-exp(-mu))*pop;

  for(i in 2:N){
    // use +40 so the center is not near likely seasonal peak
    beta_rsv[i] = b * (a * cos (2 * pi() * ((i - 40) / 52.0 - p) ) + 1); 
    foi_rsv = beta_rsv[i] * I_rsv[i-1]/pop;
    
    Sout_rsv = (1-exp(-foi_rsv - mu))*S_rsv[i-1]*(foi_rsv/(mu + foi_rsv));
    Iout_rsv = (1-exp(-gamma - mu))*I_rsv[i-1]*(gamma/(mu + gamma));
    Rout_rsv = (1-exp(-omega - mu))*R_rsv[i-1]*(omega/(mu + omega));
    
    Sdeath_rsv = (1-exp(-foi_rsv - mu))*S_rsv[i-1]*(mu/(mu + foi_rsv));
    Ideath_rsv = (1-exp(-gamma - mu))*I_rsv[i-1]*(mu/(mu + gamma));
    Rdeath_rsv = (1-exp(-omega - mu))*R_rsv[i-1]*(mu/(mu + omega));
    
    S_rsv[i] = S_rsv[i-1] + birth - Sout_rsv + Rout_rsv - Sdeath_rsv;
    I_rsv[i] = I_rsv[i-1] + Sout_rsv - Iout_rsv - Ideath_rsv;
    Ipred_rsv[i] = Sout_rsv * rho_rsv * tests_rsv[i]; // use incidence (Eout)
    R_rsv[i] = R_rsv[i-1] + Iout_rsv - Rout_rsv - Rdeath_rsv;
    
  }
  
  max_RSV_inc = max(I_rsv);
  
  for (i in 2:N) {
    
    beta_mpv[i] = r * b * (1 + c * I_rsv[i]/max_RSV_inc) * (a * cos (2 * pi() * ((i - 40) / 52.0 - p) ) + 1); 
    foi_mpv = beta_mpv[i] * I_mpv[i-1]/pop;
    
    Sout_mpv = (1-exp(-foi_mpv - mu))*S_mpv[i-1]*(foi_mpv/(mu + foi_mpv));
    Iout_mpv = (1-exp(-gamma - mu))*I_mpv[i-1]*(gamma/(mu + gamma));
    Rout_mpv = (1-exp(-omega - mu))*R_mpv[i-1]*(omega/(mu + omega));
    
    Sdeath_mpv = (1-exp(-foi_mpv - mu))*S_mpv[i-1]*(mu/(mu + foi_mpv));
    Ideath_mpv = (1-exp(-gamma - mu))*I_mpv[i-1]*(mu/(mu + gamma));
    Rdeath_mpv = (1-exp(-omega - mu))*R_mpv[i-1]*(mu/(mu + omega));
    
    S_mpv[i] = S_mpv[i-1] + birth - Sout_mpv + Rout_mpv - Sdeath_mpv;
    I_mpv[i] = I_mpv[i-1] + Sout_mpv - Iout_mpv - Ideath_mpv;
    Ipred_mpv[i] = Sout_mpv * rho_mpv * tests_mpv[i]; // use incidence (Eout)
    R_mpv[i] = R_mpv[i-1] + Iout_mpv - Rout_mpv - Rdeath_mpv;
    
  }
}

model {
  S0_mpv ~ beta(2, 98);
  I0_mpv ~ beta(2, 98);
  S0_rsv ~ beta(2, 98);
  I0_rsv ~ beta(2, 98);
  rho_mpv ~ beta(1, 99);
  rho_rsv ~ beta(1, 99);
  b ~ normal(3, 1);
  r ~ normal(1, 0.2);
  a ~ normal(0.5, 0.1);
  c ~ normal(0, 0.2);
  p ~ normal(0.5, 0.1);
  phi_mpv ~ normal(0, 10);
  phi_rsv ~ normal(0, 10);
  
  cases_mpv ~ neg_binomial_2(Ipred_mpv, phi_mpv);
  cases_rsv ~ neg_binomial_2(Ipred_rsv, phi_rsv);
}
