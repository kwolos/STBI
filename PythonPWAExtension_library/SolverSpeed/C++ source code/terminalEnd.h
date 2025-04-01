struct rparamsT //structure with the parameters
{ 
    double dt; 
    double dxdt; 
    double dtdx; 
    double Anm; 
    double Anm_1; 
    double An1m_1; 
    double Qnm; 
    double Qn1m_1; 
    double A0m; 
    double fm; 
    double Anhm_1; 
    double p0; 
    double pC; 
    double pT; 
    double R1; 
    double R2; 
    double C; 
};


double rhs_terminal(double x, void *pars)
{
    //getting parameters
    //rewritning params
    double dt = ((struct rparamsT *) pars)->dt;  
    double dxdt = ((struct rparamsT *) pars)->dxdt;  
    double dtdx = ((struct rparamsT *) pars)->dtdx;  
    double Anm = ((struct rparamsT *) pars)->Anm;  
    double Anm_1 = ((struct rparamsT *) pars)->Anm_1;  
    double An1m_1 = ((struct rparamsT *) pars)->An1m_1;  
    double Qnm = ((struct rparamsT *) pars)->Qnm;  
    double Qn1m_1 = ((struct rparamsT *) pars)->Qn1m_1;  
    double A0m = ((struct rparamsT *) pars)->A0m;  
    double fm = ((struct rparamsT *) pars)->fm; 
    double Anhm_1 = ((struct rparamsT *) pars)->Anhm_1;  
    double p0 = ((struct rparamsT *) pars)->p0;  
    double pC = ((struct rparamsT *) pars)->pC;  
    double pT = ((struct rparamsT *) pars)->pT;  
    double R1 = ((struct rparamsT *) pars)->R1;  
    double R2 = ((struct rparamsT *) pars)->R2;  
    double C = ((struct rparamsT *) pars)->C;

    // double Pnm   = p0 + fm * (sqrt(Anm/A0m)-1);
    double Pnm  = p0 + fm * (1 - sqrt(A0m/Anm));
    
    double An1m = x;   // A^(n+1)_m <-- searched A
    double Pn1m = p0 + fm * (1 - sqrt(A0m/An1m));  // P^(n+1)_m

    const double f0 = ((Pn1m - Pnm)/dt 
                        - (Qn1m_1 + dxdt*((Anm_1 + Anm - An1m_1 - An1m)/2.0) - Qnm)*R1/dt)*C 
                            + (Pnm - pT - Qnm*R1)/R2 - Qnm; 
    
    return f0;
}

double rhs_terminal_deriv(double x, void *pars)
{
    //getting parameters
    //rewritning params
    double dt = ((struct rparamsT *) pars)->dt;  
    double dxdt = ((struct rparamsT *) pars)->dxdt;  
    double dtdx = ((struct rparamsT *) pars)->dtdx;  
    double Anm = ((struct rparamsT *) pars)->Anm;  
    double Anm_1 = ((struct rparamsT *) pars)->Anm_1;  
    double An1m_1 = ((struct rparamsT *) pars)->An1m_1;  
    double Qnm = ((struct rparamsT *) pars)->Qnm;  
    double Qn1m_1 = ((struct rparamsT *) pars)->Qn1m_1;  
    double A0m = ((struct rparamsT *) pars)->A0m;  
    double fm = ((struct rparamsT *) pars)->fm;  
    double p0 = ((struct rparamsT *) pars)->p0;  
    double pC = ((struct rparamsT *) pars)->pC;  
    double pT = ((struct rparamsT *) pars)->pT;  
    double R1 = ((struct rparamsT *) pars)->R1;  
    double R2 = ((struct rparamsT *) pars)->R2;  
    double C = ((struct rparamsT *) pars)->C;

    double An1m = x; 
    
    double dPn1m = fm * sqrt(A0m)/(2 * pow(sqrt(An1m), 3)); 
    double df_dAn1_m = dPn1m/dt * C + dxdt/2./dt * R1 * C; 

    return df_dAn1_m;
}

bool updateTerminalEnd(vessel* P, params par, double t) {
    //this function updates the inflow sondition
    
    int Np = P->N;
    
    if (P->Q[Np-1] != 0 ) {//if there is any flow

        double dt     = P->dt;
        double dxdt   = P->dx/P->dt; 
        double dtdx   = P->dt/P->dx;

        double Anm    = P->A[Np];       // A^n_m         <-- not computed yet, in memory is stored Anm in place An1m
        double Anm_1  = P->Anm_1;       // A^n_(m-1)     <-- I need store it in the memory
        double An1m_1 = P->A[Np-1];     // A^(n+1)_(m-1) <-- already computed in 'updateInterior'
        double Qnm    = P->Q[Np];       // Q^n_m         <-- not computed yet, in memory is stored Qnm in place Qn1m
        double Qn1m_1 = P->Q[Np-1];     // Q^(n+1)_(m-1) <-- already computed in 'updateInterior'
        double A0m    = P->A0[Np+1];    // A0_m          <-- invariant in time 
        double fm     = P->f[Np+1];     // f_m           <-- invariant in time
        double Anhm_1 = P->Anh[Np-1];   //

        double p0  = P->p0; 
        double pC  = P->pC;
        double pT  = par.pT*pow(P->rc,4)/pow(P->qc,2)/P->rho;
        double R1  = P->R1;
        double R2  = P->R2;
        double C   = P->C;
        
        struct rparamsT p = {dt, dxdt, dtdx, Anm, Anm_1, An1m_1, Qnm, Qn1m_1, A0m, fm, Anhm_1, p0, pC, pT, R1, R2, C};
        
        int status;
        size_t iter = 0;
        
        double xt, x = P->A[Np];
        
        
        double f = rhs_terminal(x,&p);
        double th = 1.;
        double df = rhs_terminal_deriv(x,&p);
        
        do
        {
            iter++;
            xt = x;
            x -= th*f/df;
            
            f = rhs_terminal(x,&p);
            df = rhs_terminal_deriv(x,&p);
            
            if (f == 0) {
                x = xt;
                th/=2.;
                
                f = rhs_terminal(x,&p);
                df = rhs_terminal_deriv(x,&p);
            } else {
                
                if (std::abs(f) < 1e-7) {
                    status = GSL_SUCCESS;
                    break;
                }
            }
            
        }
        while (iter < 20);
        
        if(status!=GSL_SUCCESS){
            // printf("Warning: lack of convergence at terminal end, ID = %d, time = %f, err = %f.\n", P->ID,t,std::abs(f));
            //return true;
        }
        
        //assigning the output
        const double An1m = x;
        const double Qn1m = Qn1m_1 + dxdt/2. * (Anm_1 + Anm - An1m_1 - An1m);

        P->A[Np]  = An1m;
        P->Q[Np]  = Qn1m;
        
    }
    
    return false;
}