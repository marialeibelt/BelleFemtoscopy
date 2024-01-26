#include <iostream>
#include "TFile.h"
#include "TLorentzVector.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TTree.h"
#include "TDatabasePDG.h"
#include <deque>
#include <vector>
#define _USE_MATH_DEFINES
#include <cmath>
#include <math.h>

auto dbPDG = new TDatabasePDG();

float RelMom(TLorentzVector PartOne, TLorentzVector PartTwo){
    auto trackSum = PartOne + PartTwo;                                                     
    auto boost = trackSum.BoostVector();
    PartOne.Boost(-boost);
    PartTwo.Boost(-boost);
    return 0.5*(PartOne-PartTwo).P();
}    

void fillOneRChoice(TString inputfile, TH1D* SE, TH1D* ME, TH1D* pHist, TH1D* pTHist, TH2D* PhTh, TH2D* PhThME, float sphericityCutMin, float sphericityCutMax, float Nmin, float Nmax){
    TFile *input = new TFile(inputfile, "read");
    if(!input){ return;}

    TTree *tree = (TTree*)input->Get("tree");
    if(!tree){ return;}

    Int_t nentries = (Int_t)tree->GetEntries();
    int NRun = nentries; 

    auto nameCandidate = "__candidate__";
    auto nameCandidates = "__ncandidates__";
    auto nameEvent = "__event__";
    auto namePx = "px";
    auto namePy = "py";
    auto namePz = "pz";
    auto nameprobPoverK = "probPoverK";
    auto nameprobPoverPi = "probPoverPi";
    auto namenSVDHits = "nSVDHits";
    auto namenCDCHits = "nCDCHits";
    auto namedr = "dr";
    auto namedz = "dz";
    auto namefoxWolframR2 = "foxWolframR2";
    auto namecharge = "charge";
    auto namenTracks = "nTracks";
    
    int candidate = -9999, ncandidates = -9999, Event = -9999, event = -2837, Nextevent = -9898;
    double px = -9999., py = -9999., pz = -9999., probPoverK = -9999., probPoverPi = -9999., nSVDHits = -9999., nCDCHits = -9999., dr = -9999., dz = -9999., foxWolframR2 = -9999., charge = -9999., nTracks = -9999;
    
    tree->SetBranchAddress(nameCandidate, &candidate);
    tree->SetBranchAddress(nameCandidates, &ncandidates);
    tree->SetBranchAddress(nameEvent, &event);
    tree->SetBranchAddress(namePx, &px);
    tree->SetBranchAddress(namePy, &py);
    tree->SetBranchAddress(namePz, &pz);
    tree->SetBranchAddress(nameprobPoverK, &probPoverK);
    tree->SetBranchAddress(nameprobPoverPi, &probPoverPi);
    tree->SetBranchAddress(namenSVDHits, &nSVDHits);
    tree->SetBranchAddress(namenCDCHits, &nCDCHits);
    tree->SetBranchAddress(namedr, &dr);
    tree->SetBranchAddress(namedz, &dz);
    tree->SetBranchAddress(namefoxWolframR2, &foxWolframR2);
    tree->SetBranchAddress(namecharge, &charge);
    tree->SetBranchAddress(namenTracks, &nTracks);  
    
    TLorentzVector PartOne;
    TLorentzVector PartTwo;
    int pdg1 = 2212; //proton
    int pdg2 = 2212;
    int pdg3 = 11; //electron
    int NLayers = 5;   
    float pT = -999;
    float eventOfPreviousOuterLoop = -9999;

    std::vector<TLorentzVector>LayerCon; 
    std::deque<std::vector<TLorentzVector>>Con; 

    for(int i=0;i<NRun;i++){
        tree->GetEntry(i);
        Event=event;
        pT = sqrt(px*px+py*py);
        if(charge<0&& abs(pT)>0.1 && nSVDHits>3 && nCDCHits>20 && dr<0.1 && abs(dz)<0.5 && foxWolframR2>sphericityCutMin && foxWolframR2<sphericityCutMax && nTracks>=Nmin && nTracks<=Nmax){
            float numberOfSelectedAntiprotons = 1;
            for(int j=i+1; j<NRun; j++){
                tree->GetEntry(j);
                if(Event!=event){break;}
                pT = sqrt(px*px+py*py);
                if(charge<0){
                    if(abs(pT)>0.1 && nSVDHits>3 && nCDCHits>20 && dr<0.1 && abs(dz)<0.5 && foxWolframR2>sphericityCutMin && foxWolframR2<sphericityCutMax &&  nTracks>=Nmin && nTracks<=Nmax){numberOfSelectedAntiprotons++;}
                }
            }
            if(numberOfSelectedAntiprotons<2){continue;}

            //start analysis
            tree->GetEntry(i);
            double charge1 = charge;
            float p = sqrt(px*px+py*py+pz*pz);
            float pT1 = sqrt(px*px+py*py);
            pHist->Fill(p);
            pTHist->Fill(pT1);
            double svd1 = nSVDHits;
            double cdc1 = nCDCHits;
            double DR1 = dr;
            double DZ1 =dz;
            double R21 =foxWolframR2;
            int NPinLayer = numberOfSelectedAntiprotons; 
            Event=event;
            auto Px1=px;
            auto Py1=py;
            auto Pz1=pz;
            PartOne.SetXYZM(Px1,Py1,Pz1,dbPDG->GetParticle(pdg1)->Mass());
            LayerCon.push_back(PartOne);

            for(int j=i+1; j<NRun; j++){
                tree->GetEntry(j);
                double charge2 = charge;
                if(Event!=event){break;}
                if(charge2<0){
                    auto Px2=px;
                    auto Py2=py;
                    auto Pz2=pz;
                    float pT2 = sqrt(px*px+py*py);
                    double svd2 = nSVDHits;
                    double cdc2 = nCDCHits;
                    double DR2 = dr;
                    double DZ2 =dz;
                    double R22 =foxWolframR2;
                    if(abs(pT2)>0.1 && svd2>3 && cdc2>20 && DR2<0.1&&abs(DZ2)<0.5  && R22>sphericityCutMin && R22<sphericityCutMax && nTracks>=Nmin && nTracks<=Nmax){ 
                        PartTwo.SetXYZM(Px2,Py2,Pz2,dbPDG->GetParticle(pdg2)->Mass());
                        float relmomSE = RelMom(PartOne, PartTwo);
                        SE->Fill(relmomSE);
                        double DeltaPhi = PartOne.Phi()-PartTwo.Phi();
                        double DeltaTheta = PartOne.Theta()-PartTwo.Theta();
                        if(DeltaPhi>-9 && DeltaPhi<9 && DeltaTheta>-9 && DeltaTheta<9){
                            PhTh->Fill(DeltaPhi,DeltaTheta);
                        }
                        LayerCon.push_back(PartTwo);
                    }
                }
            } 
            if(eventOfPreviousOuterLoop != Event && LayerCon.size() == NPinLayer){
                Con.push_back(LayerCon);
            }  else if (eventOfPreviousOuterLoop != Event && LayerCon.size() != NPinLayer){
                std::cout<<"SOMETHING IS WRONG"<<std::endl;
            }

            //Start Event Mixing
            if(Con.size() == NLayers){
                std::vector<TLorentzVector>Layer = Con[0];
                for(int ll=0; ll < Layer.size(); ll++){
                    TLorentzVector proton = Layer[ll];
                    for(int j=1; j<NLayers; j++){
                        std::vector<TLorentzVector>Layer2 = Con[j];
                        int NPinLayer2 = Layer2.size();
                        for(int k=0; k < NPinLayer2; k++){
                            TLorentzVector mixproton = Layer2[k];
                            float relmomME = RelMom(proton, mixproton);
                            ME->Fill(relmomME);
                            double DeltaPhiME = proton.Phi()-mixproton.Phi();
                            double DeltaThetaME = proton.Theta()-mixproton.Theta();
                            if(DeltaPhiME>-9 && DeltaPhiME<9 && DeltaThetaME>-9 && DeltaThetaME<9){
                                PhThME->Fill(DeltaPhiME,DeltaThetaME);
                            }  
                        }
                    }
               }
               Con.pop_front();
            }
            LayerCon.clear();
        }
        eventOfPreviousOuterLoop = Event;
    }
    input->Close();
}

//Store Information for Histograms
int main(int argc, char** argv){                                                                     
    auto inputfile=argv[1];                                                             
    auto outputfile=argv[2];   
    const int nHists =  3; 
    const int nHistsforN = 11;
    TString pNamesNBins[nHistsforN] = {"0_2","3_4","5_6","7_8","9_10","11_12","13_15","16_18","19_21", "22_24","25_100"};  
    float nBinsMom =  1000;
    std::vector<TH1D*> pVec;
    std::vector<TH1D*> pTVec;
    TString pNames[3] = {"p","pY","pJet"};
    TString pTNames[3] = {"pT","pTY","pTJet"}; 
    float limitsMom[2] = {0,10};

    for(int i = 0; i<nHists; i++){
        for(int j = 0; j<nHistsforN;j++){
            pVec.push_back(new TH1D(pNames[i]+pNamesNBins[j],pNames[i]+pNamesNBins[j], nBinsMom, limitsMom[0], limitsMom[1]));
            pVec[i]->SetTitle("Total momentum;p (GeV/c);Entries");
            pTVec.push_back(new TH1D(pTNames[i]+pNamesNBins[j],pTNames[i]+pNamesNBins[j], nBinsMom, limitsMom[0], limitsMom[1]));
            pTVec[i]->SetTitle("Transverse momentum;p_{T} (GeV/c);Entries");
        }
    }
    float nBins =  4000;
    std::vector<TH1D*> SEVec;
    std::vector<TH1D*> MEVec;
    TString SENames[3] = {"SE","SEY","SEJet"};
    TString MENames[3] = {"ME","MEY","MEJet"};
    float limitsSEME[2] = {0,3};

    for(int i = 0; i<nHists; i++){
        for(int j = 0; j<nHistsforN;j++){
            SEVec.push_back(new TH1D(SENames[i]+pNamesNBins[j],SENames[i]+pNamesNBins[j], nBins, limitsSEME[0], limitsSEME[1]));
            MEVec.push_back(new TH1D(MENames[i]+pNamesNBins[j],MENames[i]+pNamesNBins[j], nBins, limitsSEME[0], limitsSEME[1]));
            //std::cout<<SENames[i]+pNamesNBins[j]<<std::endl;
        }
    }
    float nBinsPhiDelta =  1000;
    std::vector<TH2D*> PhThSEVec;
    std::vector<TH2D*> PhThMEVec;
    TString PhThSENames[3] = {"PhTh","PhThY","PhThJet"};
    TString PhThMENames[3] = {"PhThME","PhThMEY","PhThMEJet"};
    float limitsPhiTheta[2] = {7,4};

    for(int i = 0; i<nHists; i++){
        for(int j = 0; j<nHistsforN;j++){
            PhThSEVec.push_back(new TH2D(PhThSENames[i]+pNamesNBins[j],PhThSENames[i]+pNamesNBins[j], nBinsPhiDelta, -limitsPhiTheta[0], limitsPhiTheta[0], nBinsPhiDelta, -limitsPhiTheta[1], limitsPhiTheta[1]));
            PhThMEVec.push_back(new TH2D(PhThMENames[i]+pNamesNBins[j],PhThMENames[i]+pNamesNBins[j], nBinsPhiDelta, -limitsPhiTheta[0], limitsPhiTheta[0], nBinsPhiDelta, -limitsPhiTheta[1], limitsPhiTheta[1]));
        }
    }
    
    // Run analysis for different R2 values
    float sphericityCutMin[nHists] = {0, 0, 0.7};
    float sphericityCutMax[nHists] = {1,0.2,1};
    float Nmin[nHistsforN] = {0,3,5,7,9,11,13,16,19,22,25};
    float Nmax[nHistsforN] = {2,4,6,8,10,12,15,18,21,24,100};
    
    for(int i = 0; i<nHists; i++){
        for(int j = 0; j<nHistsforN; j++){
            fillOneRChoice(inputfile, SEVec[i*nHistsforN+j], MEVec[i*nHistsforN+j], pVec[i*nHistsforN+j], pTVec[i*nHistsforN+j], PhThSEVec[i*nHistsforN+j], PhThMEVec[i*nHistsforN+j], sphericityCutMin[i], sphericityCutMax[i],Nmin[j],Nmax[j]);
            //std::cout<<i*nHistsforN+j<<" Nmin[j]: "<<Nmin[j]<<"     Nmax[j]: "<<Nmax[j]<<std::endl;
        }
    }
    std::cout<<"writing"<<std::endl; 
    
    TFile *output = new TFile(outputfile, "RECREATE");
    
    for(int i = 0; i<nHists*nHistsforN; i++){
        SEVec[i]->Write();
        MEVec[i]->Write();
        PhThSEVec[i]->Write();
        PhThMEVec[i]->Write();
        pVec[i]->Write();
        pTVec[i]->Write();
    }
    output->Close();
    return 0;
}