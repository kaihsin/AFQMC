#include <iostream>
#include <string>
#include <cstdlib>
#include <cmath> 
#include <deque>
#include "cytnx.hpp"
#include "dqmc.hpp"
#include "Parser.hpp"
#include "rng.hpp"

using namespace std;
typedef cytnx::Accessor ac;

vector<cytnx::Tensor> Stable_Contract(const deque<int> &Rot, const vector<cytnx::Tensor> &Bs, const unsigned int &grpl){
    unsigned int timeL = Bs.size();
    unsigned int r_start, r_end;
    int NGrp = (timeL + grpl -1)/grpl;
    vector<cytnx::Tensor> rBs(NGrp);

    // contract grpl Matrices to NGrp:
    for(int rl=0;rl<NGrp;rl++){
        r_start = rl * grpl;
        r_end = ((r_start + grpl) > timeL) ? timeL:(r_start + grpl);
        
        rBs[rl] = Bs[Rot[r_start]];
        for(int l=r_start+1; l < r_end;l++){
            rBs[rl] = cytnx::linalg::Dot(rBs[rl],Bs[Rot[l]]);
        }
    }

    // Stablization by QR:
    auto udv = cytnx::linalg::Qdr(rBs[NGrp-1]);
    auto V = udv[2];
    for(int rl=NGrp-2;rl >=0;rl--){
        rBs[rl] = cytnx::linalg::Dot(rBs[rl],cytnx::linalg::Dot(udv[0],cytnx::linalg::Diag(udv[1])));
        udv = cytnx::linalg::Qdr(rBs[rl]);
        V = cytnx::linalg::Dot(udv[2],V);
    }
    udv[2] = V;
    return udv;

}

cytnx::Tensor makeB(const cytnx::Tensor &aluxf, const double &lambda, const double &dt, const cytnx::Tensor &eK, const int &sigma){
    double lt = lambda*dt;
    return cytnx::linalg::Dot(cytnx::linalg::Diag(cytnx::linalg::Exp(-sigma*lt*aluxf)),eK);   
     
}
cytnx::Tensor makeG(const deque<int> &Rot, const vector<cytnx::Tensor> &Bs, unsigned int &grpl){
    auto UDV = Stable_Contract(Rot,Bs,grpl);

    cytnx::linalg::InvM_(UDV[0]);
    cytnx::linalg::InvM_(UDV[2]);

    auto TP1 = cytnx::linalg::Dot(UDV[0],UDV[2]) + cytnx::linalg::Diag(UDV[1]);
    auto udv = cytnx::linalg::Qdr(TP1);
    cytnx::linalg::InvM_(udv[0]);
    cytnx::linalg::Inv_(udv[1],1.0e-15);
    cytnx::linalg::InvM_(udv[2]);

    UDV[2] = cytnx::linalg::Dot(UDV[2],udv[2]);
    udv[0] = cytnx::linalg::Dot(udv[0],UDV[0]);
    return cytnx::linalg::Dot(UDV[2], cytnx::linalg::Dot(cytnx::linalg::Diag(udv[1]),udv[0]));

}
cytnx::Tensor evolveG(const cytnx::Tensor &Bsgm_curr, const cytnx::Tensor &Gsgm_bef){
    auto invB = cytnx::linalg::InvM(Bsgm_curr);
    return cytnx::linalg::Dot(Bsgm_curr,cytnx::linalg::Dot(Gsgm_bef,invB));
}
cytnx::Tensor updateG(const double &delta, const double &Prob, const int &site_i, const cytnx::Tensor &Gsgm){
    //auto out = Gsgm.clone();
    
    cytnx::Tensor tmp = Gsgm(site_i);
    tmp(site_i) = tmp(site_i) + 1;

    auto out = cytnx::Tensor(Gsgm(0)) - (Gsgm(0,site_i).item<double>()*delta/Prob)*tmp;
    out.reshape_(1,-1);
    for(int mu = 1; mu < Gsgm.shape()[1];mu++){
        out.append(cytnx::Tensor(Gsgm(mu)) - (Gsgm(mu,site_i).item<double>()*delta/Prob)*tmp);
    }
    return out;

}
int main(int argc, char* argv[]){

    if(argc != 6){
        cout << "./dqmc <rc file> <t> <mu> <U> <beta>" << endl;
        exit(99);
    }

    /// global vars:
    int S_UP   =  1;
    int S_DOWN = -1;    
    

    /// simulation parameter:
    unsigned int BINNUM,BINSZ,EQUIN; //MC statistic tune.
    unsigned int nX; // linear system size
    //unsigned int nY; 
    unsigned int N    ; // total system size 
    unsigned int seed ; // seed for RNG
    unsigned int TrotM; // Trotter decomp. steps
    unsigned int CRTM ; // citical number m for sequential matrix multiplication.
    string ID; //unique ID for this simulation. 


    ///model parameter:
    double t,mu,U,beta;
    
   

    /// Use parser to read the rc file:
    Parser pars; 
    pars.Bind("ID",ID);
    pars.Bind("BINNUM",BINNUM);
    pars.Bind("BINSZ",BINSZ  );
    pars.Bind("EQUIN",EQUIN  );
    pars.Bind("nX",nX  );
    pars.Bind("seed",seed);
    pars.Bind("TrotM",TrotM);
    pars.Bind("CRTM",CRTM);
    pars.Parse(argv[1]);
    pars.Check_All();
    
    t = atof(argv[2]);
    mu = atof(argv[3]);
    U = atof(argv[4]);
    beta = atof(argv[5]);

    /// additional materials:
    N = nX; //nX*nY for 2D
    double dt = beta/TrotM; 
    double lambda = acos(exp(-dt*U*0.5))/dt;
    RNG rng(seed);

    pars.PrintVars();
    cout << "t  :" << t << endl;
    cout << "mu :" << mu << endl;
    cout << "U  :" << t << endl;
    cout << "Beta :" << beta << endl;

    //===================================================
    // Prepare materials:

    // 1. create K matrix (kinetic):
    auto K = cytnx::zeros({nX,nX},cytnx::Type.Double);
    if(nX==2){
        K(0,0) = K(1,1) = -mu;
        K(0,1) = K(1,0) = -t;
    }else{
        for(int i=0;i<nX;i++){
            K(i,(i + nX- 1) % nX) += -t;
            K(i,i) += -mu;
        }
    }
    auto eK = cytnx::linalg::ExpH(K,-dt);
    
    // 2. create auxiliary fields:
    vector<cytnx::Tensor> AulxS(TrotM,cytnx::Tensor({N},cytnx::Type.Double));
    for(int i=0;i<TrotM;i++){
        for(int j=0;j<N;j++)
            AulxS[i](j) = rng() < 0.5? 1:-1;
    }

    // 3. Construct B mats vector & rotator:
    deque<int> Rotator;
    vector<cytnx::Tensor> Bup(TrotM);
    vector<cytnx::Tensor> Bdown = cytnx::vec_clone(Bup);
    
    for(int it=0;it<TrotM;it++){
        Rotator.push_back((TrotM-it)%TrotM);
        Bup[Rotator.back()]   = makeB(AulxS[Rotator.back()],lambda,dt,eK,S_UP);
        Bdown[Rotator.back()] = makeB(AulxS[Rotator.back()],lambda,dt,eK,S_DOWN);
    }

    // 4. Greens fucntion series:
    vector<cytnx::Tensor> Gup(TrotM);
    vector<cytnx::Tensor> Gdown(TrotM);

    // 5. temp variables:
    double delUp,delDwn,probUp,probDwn,TransProb;
    int curr_sign;


    //======================================================
    // Main algos:

    // a: equillibrum stage:
    for(int eq=0;eq<EQUIN;eq++){
        //cout << eq << endl;    
        //Time evolution:
        for(int it=0;it<TrotM;it++){
            //cout << it << endl;
            if(it%CRTM==0){
                Gup[it] = makeG(Rotator,Bup,CRTM);
                Gdown[it] = makeG(Rotator,Bdown,CRTM);
            }else{
                Gup[it]   = evolveG(Bup[Rotator.front()]  ,Gup[it-1]);
                Gdown[it] = evolveG(Bdown[Rotator.front()],Gdown[it-1]);
            }
            
        }
        
        // each sites in one time update.
        for( int s=0;s<N;s++){
            delUp  = exp(2.*dt*lambda*AulxS[t](s).item<double>()) -1;
            delDwn = exp(-2.*dt*lambda*AulxS[t](s).item<double>()) -1;
            probUp  = (1 + (1 - Gup[t](s, s).item<double>() )*delUp );
            probDwn = (1 + (1 - Gdown[t](s,s).item<double>())*delDwn);
            TransProb = probUp*probDwn;
        
            if(rng() < abs(TransProb)){
                AulxS[t](s)*=-1;//update AuxField.
                Gup[t] = updateG(delUp,probUp,s,Gup[t]);
                Gdown[t] = updateG(delDwn,probDwn,s,Gdown[t]);
            }

        }


    }//eq





    

    return 0;
}
