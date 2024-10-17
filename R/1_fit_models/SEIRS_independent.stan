data {
  int N;                     // number of weeks
  int cases[N];              // detections 
  real tests[N];             // moving average of tests scaled [0,1] 
  int pop;                   // population size
  real gamma;                // assumed recovery rate
  real sigma;                // assumed 1/duration of latent period
  real mu;                   // assumed birth/death rate
  real omega;                // waning rate (assumed same for both)
}

parameters {
  real<lower=0, upper=1> S0; // initial conditions
  real<lower=0, upper=1> E0;
  real<lower=0, upper=1> I0;
  real<lower=0, upper=1> rho; // reporting rate
  real<lower=0> b;            // base transmission rate
  real<lower=0, upper=1> a;   // amplitude of seasonal forcing
  real<lower=0, upper=1> p;   // phase of seasonal forcing    
  real<lower=0> phi;
}

transformed parameters {
  vector<lower=0>[N] S;
  vector<lower=0>[N] E;
  vector<lower=0>[N] I;
  vector<lower=0>[N] Ipred;
  vector<lower=0>[N] R;
  vector<lower=0>[N] beta;
  real<lower=0> foi;
  real<lower=0> Sdeath;
  real<lower=0> Edeath;
  real<lower=0> Ideath;
  real<lower=0> Rdeath;
  real<lower=0> Sout;
  real<lower=0> Eout;
  real<lower=0> Iout;
  real<lower=0> Rout;
  real<lower=0> birth;
  
  S[1] = S0 * pop;
  E[1] = E0 * pop;
  I[1] = I0 * pop;
  Ipred[1] = I0 * pop * rho * tests[1];
  R[1] = (1-S0-E0-I0) * pop;
  
  beta[1] = b;
  
  birth = (1-exp(-mu))*pop;

  for(i in 2:N){
    // use +40 so the center is not near likely seasonal peak
    beta[i] = b * (a * cos (2 * pi() * ((i - 40) / 52.0 - p) ) + 1); 
    foi = beta[i] * I[i-1]/pop;
    
    Sout = (1-exp(-foi - mu))*S[i-1]*(foi/(mu + foi));
    Eout = (1-exp(-sigma - mu))*E[i-1]*(sigma/(mu + sigma));
    Iout = (1-exp(-gamma - mu))*I[i-1]*(gamma/(mu + gamma));
    Rout = (1-exp(-omega - mu))*R[i-1]*(omega/(mu + omega));
    
    Sdeath = (1-exp(-foi - mu))*S[i-1]*(mu/(mu + foi));
    Edeath = (1-exp(-sigma - mu))*E[i-1]*(mu/(mu + sigma));
    Ideath = (1-exp(-gamma - mu))*I[i-1]*(mu/(mu + gamma));
    Rdeath = (1-exp(-omega - mu))*R[i-1]*(mu/(mu + omega));
    
    S[i] = S[i-1] + birth - Sout + Rout - Sdeath;
    E[i] = E[i-1] + Sout - Eout - Edeath;
    I[i] = I[i-1] + Eout - Iout - Ideath;
    Ipred[i] = Eout * rho * tests[i]; // use incidence (Eout)
    R[i] = R[i-1] + Iout - Rout - Rdeath;
    
  }
}

model {
  S0 ~ beta(2, 98);
  E0 ~ beta(2, 98);
  I0 ~ beta(2, 98);
  rho ~ beta(1, 99);
  b ~ normal(3, 1);
  a ~ normal(0.5, 0.1);
  p ~ normal(0.5, 0.1);
  phi ~ normal(0, 10);
  
  cases ~ neg_binomial_2(Ipred, phi);
}
