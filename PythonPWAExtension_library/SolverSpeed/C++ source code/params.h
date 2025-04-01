class params{
public:
    
    double pT; //terminal pressure
    double p0; //reference pressure
    
    double dx;
    double dt;
    
    double rho;
    double mu;
    
    //parameters assoctited with inflow condition
    double q0, T, tau, qB;
    
    //parameters associated with elasticity
    double k1, k2, k3;
    
    bool terminate;

    //parameters associated with elastance model 
    double tm;      // time to the oneset constant elastance 

    double a, b;    // parameters in Phi function 
    double Emin, Emax; 

    double V0;      
    double VlvInit, plvInit;    
    double pa;      
    double Vb, Vbup;      
    double R;       
    double Llv; 
    double Time;
    double xF0; 

    double Qprev, Qlaprev, Vprev;
    int isQZero;
    int isQgret; 
    double Vfin, VfinUp; 
    
    double value;
    double valuePrev;
    bool aorticValve;
    int phase; 
    double pla; 
    double constPhase3; 
    double Lla; 
    double Rla;
    double div;
    double Ela; 
    double V0a; 
    double Qparam1, Qparam2; 
    double l;   // length of the jump

    int iterator; 

    // parameters related with the viscoelastic model of the artery 
    // (for brachial cuff purposes)
    double lambda_s, lambda_c; 
    double E_s, E_c; 
    double P_b; 


    // parameters connected with opening/closing valve 
    double Kvo; 
    double Kvc; 
    double plao, plac; // pressure close, open
    double Zprev; 
    double A; 
    double Mst, Mrg; 

    double model;
    bool message;

    // parameters connected with autoregulation 
    double r_d0;
    double r_dRA;
    double r_dLA;
    double r_dRM;
    double r_dLM;
    double r_dRP;
    double r_dLP;
    double h_d0;
    double h_dRA;
    double h_dLA;
    double h_dRM;
    double h_dLM;
    double h_dRP;
    double h_dLP;
    double sigma_e0;
    double sigma_eRA;
    double sigma_eLA;
    double sigma_eRM;
    double sigma_eLM;
    double sigma_eRP;
    double sigma_eLP;
    double P_dRA; 
    double P_dLM; 
    double P_dRM; 
    double P_dLP; 
    double P_dRP; 
    double P_dLA; 
    double drRA_dt; 
    double drLA_dt; 
    double drRM_dt; 
    double drLM_dt; 
    double drRP_dt; 
    double drLP_dt;
    double P_cRA; 
    double P_cLA; 
    double P_cRM; 
    double P_cLM; 
    double P_cRP; 
    double P_cLP;
    double x_RA; 
    double x_LA; 
    double x_RM; 
    double x_LM; 
    double x_RP; 
    double x_LP; 
    double K_sigma; 
    double sigma_coll; 
    double T_0; 
    double r_m; 
    double r_t; 
    double n_m; 
    double eta; 
    double K_v; 
    double K_gM; 
    double K_gA; 
    double K_gP; 
    double K_e; 
    double C_m; 
    double t_CA; 
    double G_CA; 
    double G_pv; 
    double G_f; 
    double G_cAA; 
    double G_cPP; 
    double G_cAM; 
    double G_cMP; 
    double P_an; 
    double G_o; 
    double P_s; 
    double Q_nM; 
    double Q_nA; 
    double Q_nP;
    double P_icn;
    double P_ic;

    // interstitial space
    double Q_ISRA;
    double Q_ISLA;
    double Q_ISRM;
    double Q_ISLM;
    double Q_ISRP;
    double Q_ISLP;
    double R_bfRA; 
    double R_bfLA; 
    double R_bfRM; 
    double R_bfLM; 
    double R_bfRP; 
    double R_bfLP; 
    double R_bbbRA; 
    double R_bbbLA; 
    double R_bbbRM; 
    double R_bbbLM; 
    double R_bbbRP; 
    double R_bbbLP; 
    double C_brRA; 
    double C_brLA; 
    double C_brRM; 
    double C_brLM; 
    double C_brRP; 
    double C_brLP;
  
    params(); //default constructor
};

params::params(){}; //default constructor
