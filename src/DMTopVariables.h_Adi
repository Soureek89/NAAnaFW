#ifndef _Top_Utilities_h_
#define _Top_Utilities_h_

/**
 *\Function EquationSolver:
 *
 * Reconstruct semileptonic tops observables.
 *
 * \Author A. Orso M. Iorio
 * 
 *
 *\version  $Id: 
 *
 *
*/


#include<cmath>
#include<string>
#include<iostream>
#include<vector>
#include<complex>

#include "DataFormats/Math/interface/LorentzVector.h"
#include "DataFormats/Math/interface/deltaPhi.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "TLorentzVector.h"
#include <Math/VectorUtil.h>
#include "./EquationSolver.h"

using namespace std;
//using namespace math;
class TopUtilities {
  
public:
  TopUtilities(){;}
  ~TopUtilities(){;}
  
  vector<math::PtEtaPhiELorentzVector> top4MomentaPTEtaPhiE(vector<TLorentzVector> , vector<TLorentzVector>, double, double );
  vector<TLorentzVector> top4Momenta(vector<TLorentzVector> , vector<TLorentzVector>, double, double );

  math::PtEtaPhiELorentzVector top4Momentum(TLorentzVector , TLorentzVector, double, double );
  math::PtEtaPhiELorentzVector top4Momentum(math::PtEtaPhiELorentzVector , math::PtEtaPhiELorentzVector, double, double );
  math::PtEtaPhiELorentzVector top4Momentum(float leptonPx, float leptonPy, float leptonPz, float leptonE, float jetPx, float jetPy, float jetPz, float jetE, float metPx, float metPy);
  math::XYZTLorentzVector NuMomentum(float leptonPx, float leptonPy, float leptonPz, float leptonPt, float leptonE, float metPx, float metPy );

  double topMtw(TLorentzVector, TLorentzVector, float metPx, float metPy);  
  double topMtw(math::PtEtaPhiELorentzVector lepton, math::PtEtaPhiELorentzVector jet, float metPx, float metPy);  
  double mlj(math::PtEtaPhiELorentzVector lepton, math::PtEtaPhiELorentzVector jet);
  double mlj(TLorentzVector , TLorentzVector);
  double mindeltaphi(double metphi,  vector<math::PtEtaPhiELorentzVector> jets);
  double mindeltaphi(double metphi,  vector<TLorentzVector> jets);
  double costhetapol(TLorentzVector, TLorentzVector, TLorentzVector);
  double costhetael(TLorentzVector,TLorentzVector, float metPx, float metPy);
    

};

vector <math::PtEtaPhiELorentzVector> TopUtilities::top4MomentaPTEtaPhiE(vector<TLorentzVector> leptons, vector<TLorentzVector> bjets, double metPx, double metPy){
  vector <math::PtEtaPhiELorentzVector> topCandidates;
  for (size_t i=0; i < (size_t)leptons.size();++i ){
    math::PtEtaPhiELorentzVector lep(leptons.at(i).Pt(),leptons.at(i).Eta(),leptons.at(i).Phi(),leptons.at(i).Energy());
    
    for (size_t b=0; b < (size_t)bjets.size();++b ){
      math::PtEtaPhiELorentzVector bjet(bjets.at(b).Pt(),bjets.at(b).Eta(),bjets.at(b).Phi(),bjets.at(b).Energy());
      topCandidates.push_back(top4Momentum(lep,bjet,metPx,metPy));
    }  
  }
  return topCandidates;
}

vector <TLorentzVector> TopUtilities::top4Momenta(vector<TLorentzVector> leptons, vector<TLorentzVector> bjets, double metPx, double metPy){
  vector <TLorentzVector> topCandidates;
  for (size_t i=0; i < (size_t)leptons.size();++i ){
    math::PtEtaPhiELorentzVector lep(leptons.at(i).Pt(),leptons.at(i).Eta(),leptons.at(i).Phi(),leptons.at(i).Energy());
    
    for (size_t b=0; b < (size_t)bjets.size();++b ){
      math::PtEtaPhiELorentzVector bjet(bjets.at(b).Pt(),bjets.at(b).Eta(),bjets.at(b).Phi(),bjets.at(b).Energy());
      TLorentzVector top;
      math::PtEtaPhiELorentzVector top_tmp = top4Momentum(lep,bjet,metPx,metPy);
      top.SetPtEtaPhiE(top_tmp.Pt(), top_tmp.Eta(), top_tmp.Phi(), top_tmp.E());
      topCandidates.push_back(top);

    }  
  }
  return topCandidates;
}

math::PtEtaPhiELorentzVector TopUtilities::top4Momentum(TLorentzVector lepton, TLorentzVector jet, double metPx, double metPy){
  math::PtEtaPhiELorentzVector lep(lepton.Pt(),lepton.Eta(),lepton.Phi(),lepton.Energy());
  math::PtEtaPhiELorentzVector bjet(jet.Pt(),jet.Eta(),jet.Phi(),jet.Energy());
  return top4Momentum(lep,bjet,metPx,metPy);
}

//top quark 4-momentum given lepton, met and b-jet
math::PtEtaPhiELorentzVector TopUtilities::top4Momentum(math::PtEtaPhiELorentzVector lepton, math::PtEtaPhiELorentzVector jet, double metPx, double metPy)
{
    return top4Momentum(lepton.px(), lepton.py(), lepton.pz(), lepton.energy(), jet.px(), jet.py(), jet.pz(), jet.energy(), metPx, metPy);
}

//top quark 4-momentum original function given the necessary parameters
math::PtEtaPhiELorentzVector TopUtilities::top4Momentum(float leptonPx, float leptonPy, float leptonPz, float leptonE, float jetPx, float jetPy, float jetPz, float jetE, float metPx, float metPy)
{
    float lepton_Pt = sqrt( (leptonPx * leptonPx) +  (leptonPy * leptonPy) );

    math::XYZTLorentzVector neutrino = NuMomentum(leptonPx, leptonPy, leptonPz, lepton_Pt, leptonE, metPx, metPy); //.at(0);;

    math::XYZTLorentzVector lep(leptonPx, leptonPy, leptonPz, leptonE);
    math::XYZTLorentzVector jet(jetPx, jetPy, jetPz, jetE);

    math::XYZTLorentzVector top = lep + jet + neutrino;
    return math::PtEtaPhiELorentzVector(top.pt(), top.eta(), top.phi(), top.E());
}


/////What it does:
//w boson mass put to pdg value
//obtained neutrino pz from kinematics
//We get a second order equation
/////In case of two positive Delta solutions:
//we choose solution with minimum |pz|
/////In case of two negative Delta solutions:
//in such case: mtw > mw
//To solve this: put mtw = mw
//Solve the equations
//In this way we must
//drop the constraints px_Nu = MET_x and py_Nu = MET_y
//Solve this by chosing the px_Nu and py_Nu that
//minimize the distance from the MET in the px-py plane
//Such minimization can be done analytically with derivatives
//and much patience. Here we exploit such analytical minimization
/////
//More detailed inline description: work in progress!
math::XYZTLorentzVector TopUtilities::NuMomentum(float leptonPx, float leptonPy, float leptonPz, float leptonPt, float leptonE, float metPx, float metPy )
{

    double  mW = 80.399;

    math::XYZTLorentzVector result;

    //  double Wmt = sqrt(pow(Lepton.et()+MET.pt(),2) - pow(Lepton.px()+metPx,2) - pow(leptonPy+metPy,2) );

    double MisET2 = (metPx * metPx + metPy * metPy);
    double mu = (mW * mW) / 2 + metPx * leptonPx + metPy * leptonPy;
    double a  = (mu * leptonPz) / (leptonE * leptonE - leptonPz * leptonPz);
    double a2 = TMath::Power(a, 2);
    double b  = (TMath::Power(leptonE, 2.) * (MisET2) - TMath::Power(mu, 2.)) / (TMath::Power(leptonE, 2) - TMath::Power(leptonPz, 2));
    double pz1(0), pz2(0), pznu(0);
    int nNuSol(0);
    if(nNuSol);

    math::XYZTLorentzVector p4nu_rec;
    math::XYZTLorentzVector p4W_rec;
    math::XYZTLorentzVector p4b_rec;
    math::XYZTLorentzVector p4Top_rec;
    math::XYZTLorentzVector p4lep_rec;

    p4lep_rec.SetPxPyPzE(leptonPx, leptonPy, leptonPz, leptonE);

    math::XYZTLorentzVector p40_rec(0, 0, 0, 0);

    if (a2 - b > 0 )
    {
        //if(!usePositiveDeltaSolutions_)
        //  {
        //  result.push_back(p40_rec);
        //  return result;
        //  }
        double root = sqrt(a2 - b);
        pz1 = a + root;
        pz2 = a - root;
        nNuSol = 2;

        //    if(usePzPlusSolutions_)pznu = pz1;
        //    if(usePzMinusSolutions_)pznu = pz2;
        //if(usePzAbsValMinimumSolutions_){
        pznu = pz1;
        if (fabs(pz1) > fabs(pz2)) pznu = pz2;
        //}


        double Enu = sqrt(MisET2 + pznu * pznu);

        p4nu_rec.SetPxPyPzE(metPx, metPy, pznu, Enu);

        //    result =.push_back(p4nu_rec);
        result = p4nu_rec;

    }
    else
    {

        // if(!useNegativeDeltaSolutions_){
        //result.push_back(p40_rec);
        //  return result;
        //    }
        //    double xprime = sqrt(mW;


        double ptlep = leptonPt, pxlep = leptonPx, pylep = leptonPy, metpx = metPx, metpy = metPy;

        double EquationA = 1;
        double EquationB = -3 * pylep * mW / (ptlep);
        double EquationC = mW * mW * (2 * pylep * pylep) / (ptlep * ptlep) + mW * mW - 4 * pxlep * pxlep * pxlep * metpx / (ptlep * ptlep) - 4 * pxlep * pxlep * pylep * metpy / (ptlep * ptlep);
        double EquationD = 4 * pxlep * pxlep * mW * metpy / (ptlep) - pylep * mW * mW * mW / ptlep;

        std::vector<long double> solutions = EquationSolve<long double>((long double)EquationA, (long double)EquationB, (long double)EquationC, (long double)EquationD);

        std::vector<long double> solutions2 = EquationSolve<long double>((long double)EquationA, -(long double)EquationB, (long double)EquationC, -(long double)EquationD);


        double deltaMin = 14000 * 14000;
        double zeroValue = -mW * mW / (4 * pxlep);
        double minPx = 0;
        double minPy = 0;

        //    std::cout<<"a "<<EquationA << " b " << EquationB  <<" c "<< EquationC <<" d "<< EquationD << std::endl;

        //  if(usePxMinusSolutions_){
        for ( int i = 0; i < (int)solutions.size(); ++i)
        {
            if (solutions[i] < 0 ) continue;
            double p_x = (solutions[i] * solutions[i] - mW * mW) / (4 * pxlep);
            double p_y = ( mW * mW * pylep + 2 * pxlep * pylep * p_x - mW * ptlep * solutions[i]) / (2 * pxlep * pxlep);
            double Delta2 = (p_x - metpx) * (p_x - metpx) + (p_y - metpy) * (p_y - metpy);

            //      std:://cout<<"intermediate solution1 met x "<<metpx << " min px " << p_x  <<" met y "<<metpy <<" min py "<< p_y << std::endl;

            if (Delta2 < deltaMin && Delta2 > 0)
            {
                deltaMin = Delta2;
                minPx = p_x;
                minPy = p_y;
            }
            //     std:://cout<<"solution1 met x "<<metpx << " min px " << minPx  <<" met y "<<metpy <<" min py "<< minPy << std::endl;
        }

        //    }

        //if(usePxPlusSolutions_){
        for ( int i = 0; i < (int)solutions2.size(); ++i)
        {
            if (solutions2[i] < 0 ) continue;
            double p_x = (solutions2[i] * solutions2[i] - mW * mW) / (4 * pxlep);
            double p_y = ( mW * mW * pylep + 2 * pxlep * pylep * p_x + mW * ptlep * solutions2[i]) / (2 * pxlep * pxlep);
            double Delta2 = (p_x - metpx) * (p_x - metpx) + (p_y - metpy) * (p_y - metpy);
            //  std:://cout<<"intermediate solution2 met x "<<metpx << " min px " << minPx  <<" met y "<<metpy <<" min py "<< minPy << std::endl;
            if (Delta2 < deltaMin && Delta2 > 0)
            {
                deltaMin = Delta2;
                minPx = p_x;
                minPy = p_y;
            }
            //  std:://cout<<"solution2 met x "<<metpx << " min px " << minPx  <<" met y "<<metpy <<" min py "<< minPy << std::endl;
        }
        //}

        double pyZeroValue = ( mW * mW * pxlep + 2 * pxlep * pylep * zeroValue);
        double delta2ZeroValue = (zeroValue - metpx) * (zeroValue - metpx) + (pyZeroValue - metpy) * (pyZeroValue - metpy);

        if (deltaMin == 14000 * 14000)return result;
        //    else std:://cout << " test " << std::endl;

        if (delta2ZeroValue < deltaMin)
        {
            deltaMin = delta2ZeroValue;
            minPx = zeroValue;
            minPy = pyZeroValue;
        }

        //    std:://cout<<" MtW2 from min py and min px "<< sqrt((minPy*minPy+minPx*minPx))*ptlep*2 -2*(pxlep*minPx + pylep*minPy)  <<std::endl;
        ///    ////Y part

        double mu_Minimum = (mW * mW) / 2 + minPx * pxlep + minPy * pylep;
        double a_Minimum  = (mu_Minimum * leptonPz) / (leptonE * leptonE - leptonPz * leptonPz);
        pznu = a_Minimum;

        //if(!useMetForNegativeSolutions_){
        double Enu = sqrt(minPx * minPx + minPy * minPy + pznu * pznu);
        p4nu_rec.SetPxPyPzE(minPx, minPy, pznu , Enu);

        //    }
        //    else{
        //      pznu = a;
        //      double Enu = sqrt(metpx*metpx+metpy*metpy + pznu*pznu);
        //      p4nu_rec.SetPxPyPzE(metpx, metpy, pznu , Enu);
        //    }

        //      result.push_back(p4nu_rec);
        result = p4nu_rec;
    }
    return result;
}

double  TopUtilities::topMtw(TLorentzVector lepton, TLorentzVector jet, float metPx, float metPy){
  math::PtEtaPhiELorentzVector lep(lepton.Pt(),lepton.Eta(),lepton.Phi(),lepton.Energy());
  math::PtEtaPhiELorentzVector bjet(jet.Pt(),jet.Eta(),jet.Phi(),jet.Energy());
  return topMtw(lep,bjet,metPx,metPy);
}

double  TopUtilities::mlj(TLorentzVector lepton, TLorentzVector jet){
  math::PtEtaPhiELorentzVector lep(lepton.Pt(),lepton.Eta(),lepton.Phi(),lepton.Energy());
  math::PtEtaPhiELorentzVector bjet(jet.Pt(),jet.Eta(),jet.Phi(),jet.Energy());
  return mlj(lep,bjet);
}

double  TopUtilities::mlj(math::PtEtaPhiELorentzVector lepton, math::PtEtaPhiELorentzVector jet){
  return (lepton + jet).mass();
}


double TopUtilities::mindeltaphi(double metphi,  vector<math::PtEtaPhiELorentzVector> jets){
  double mindphi = 100.0;
  for (size_t s =0; s<jets.size();++s){
    double dphi = deltaPhi(jets.at(s).phi(), metphi);
    if(fabs(dphi)<mindphi)mindphi = dphi;
  }
  return mindphi;
}

double TopUtilities::mindeltaphi(double metphi,  vector<TLorentzVector> jets){
  double mindphi = 100.0;
  for (size_t s =0; s<jets.size();++s){
    double dphi = deltaPhi(jets.at(s).Phi(), metphi);
    if(fabs(dphi)<mindphi)mindphi = dphi;
  }
  return mindphi;
}

double TopUtilities::topMtw(math::PtEtaPhiELorentzVector lepton, math::PtEtaPhiELorentzVector jet, float metPx, float metPy)
{
    math::PtEtaPhiELorentzVector lb = lepton + jet;
    double mlb2 = lb.mass() * lb.mass();
    double etlb = sqrt(mlb2 + lb.pt() * lb.pt());
    double metPT = sqrt(metPx * metPx + metPy * metPy);

    return sqrt( mlb2 + 2 * ( etlb * metPT - lb.px() * metPx - lb.py() * metPy ) );
}

double TopUtilities::costhetael(TLorentzVector top, TLorentzVector lepton, float metPx, float metPy)
{
  math::XYZTLorentzVector Numomentum = NuMomentum(lepton.Px(), lepton.Py(), lepton.Pz(), sqrt(lepton.Px()*lepton.Px()+lepton.Py()*lepton.Py()), lepton.E(), metPx, metPy);
  TLorentzVector Nu;
  Nu.SetPxPyPzE(Numomentum.Px(),Numomentum.Py(),Numomentum.Pz(),Numomentum.E());
  TLorentzVector W = lepton + Nu;
  TLorentzVector W_1;
  W_1.SetPxPyPzE(-W.Px(),-W.Py(),-W.Pz(),W.E());
  lepton.Boost(W_1.BoostVector());
  TLorentzVector top_1;
  top_1.SetPxPyPzE(-top.Px(),-top.Py(),-top.Pz(),top.E());
  W.Boost(top_1.BoostVector());
  lepton.Boost(top_1.BoostVector());
  double costheta=(W.Vect()*lepton.Vect())/((W.Vect()).Mag()*(lepton.Vect()).Mag());
  return costheta;
}

double TopUtilities::costhetapol(TLorentzVector lepton, TLorentzVector jet, TLorentzVector top)
{
  TLorentzVector top_1;
  top_1.SetPxPyPzE(-top.Px(),-top.Py(),-top.Pz(),top.E());
  jet.Boost(top.BoostVector());
  lepton.Boost(top_1.BoostVector());
  double costheta=(jet.Vect()*lepton.Vect())/((jet.Vect()).Mag()*(lepton.Vect()).Mag());
  return costheta;
}
  
#endif
