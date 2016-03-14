#include "Analysis/EDSHyFT/plugins/EDSHyFTSelector.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include <iostream>
#include <algorithm>
#include <cxxabi.h>
using namespace std;

bool EDSHyFTSelector::filter( edm::Event & event, const edm::EventSetup& eventSetup )
{
    bool passed = edm::FilterWrapper<SHyFTSelector>::filter( event, eventSetup );

    std::vector<reco::ShallowClonePtrCandidate> const & ijets = filter_->cleanedJets();
    reco::ShallowClonePtrCandidate const & imet = filter_->selectedMET();
    std::vector<reco::ShallowClonePtrCandidate> const & imuons = filter_->selectedMuons();
    std::vector<reco::ShallowClonePtrCandidate> const & ielectrons = filter_->selectedElectrons();
    std::vector<edm::Ptr<pat::Tau> > const & itaus = filter_->selectedTaus();

    std::auto_ptr< std::vector<pat::Jet> > jets ( new std::vector<pat::Jet> );
    std::auto_ptr< std::vector<pat::MET> > mets ( new std::vector<pat::MET> );
    std::auto_ptr< std::vector<pat::Muon> > muons ( new std::vector<pat::Muon> );
    std::auto_ptr< std::vector<pat::Electron> > electrons ( new std::vector<pat::Electron> );
    std::auto_ptr< std::vector<pat::Tau> > taus ( new std::vector<pat::Tau> );

    pat::MET const * patmet = dynamic_cast<pat::MET const *>( imet.masterClonePtr().get() ); 
    if ( patmet != 0 ){  
        mets->push_back( *patmet );
        mets->back().setP4( imet.p4() );//set back the P4 to the clonned met
    }

    typedef std::vector<reco::ShallowClonePtrCandidate>::const_iterator clone_iter;
    typedef std::vector<edm::Ptr<pat::Tau> >::const_iterator tau_iter;
    for ( clone_iter ibegin = ijets.begin(), iend = ijets.end(), i = ibegin;i != iend; ++i ) {
        pat::Jet const * ijet = dynamic_cast<pat::Jet const *>( i->masterClonePtr().get() );
        jets->push_back( *ijet );
        float deltaRMin = 999.0;
        pat::Tau const * bestTau = NULL;
        if ( ijet != 0 ) {
            // match jet object with a (potential) Tau matching it
            for ( tau_iter jbegin = itaus.begin(), jend = itaus.end(), j = jbegin; j != jend; ++j ) {
                pat::Tau const * jtau = dynamic_cast<pat::Tau const *>( j->get() );
                if (jtau == 0) {
                    continue;
                }
                float deltaR = reco::deltaR(jtau->eta(), jtau->phi(),
                                            ijet->eta(), ijet->phi());
                if ((deltaRMin > deltaR) && (deltaR < 0.3)) {
                    bestTau = jtau;
                    deltaRMin = deltaR;
                }
            }
            if (bestTau != NULL) {
                jets->back().addUserFloat("tauJetPt", bestTau->pt());
                jets->back().addUserFloat("tauJetPhi", bestTau->phi());
                jets->back().addUserFloat("tauJetEta", bestTau->eta());
                jets->back().addUserFloat("tauJetMass", bestTau->mass());
                jets->back().addUserFloat("byLooseCombinedIsolationDeltaBetaCorr3Hits", bestTau->tauID("byLooseCombinedIsolationDeltaBetaCorr3Hits") );
                jets->back().addUserInt("byMediumCombinedIsolationDeltaBetaCorr3Hits", bestTau->tauID("byMediumCombinedIsolationDeltaBetaCorr3Hits"));
                jets->back().addUserInt("byTightCombinedIsolationDeltaBetaCorr3Hits", bestTau->tauID("byTightCombinedIsolationDeltaBetaCorr3Hits"));
            }
        }
    }

    edm::Handle<std::vector<reco::Vertex> > primVtxHandle;                 
    event.getByLabel("offlinePrimaryVertices", primVtxHandle);
    for ( clone_iter jbegin = imuons.begin(), jend = imuons.end(), j = jbegin; j != jend; ++j ) {
        pat::Muon const * jmuon = dynamic_cast<pat::Muon const *>( j->masterClonePtr().get() );
        if ( jmuon == NULL )
            continue;
        muons->push_back( *jmuon );
        muons->back().addUserInt("isTight", ( (jmuon->dB() < 2) && 
                                              (fabs( jmuon->muonBestTrack()->dz(primVtxHandle->at(0).position()) ) < 0.5 )));
    }

    for ( clone_iter jbegin = ielectrons.begin(), jend = ielectrons.end(), j = jbegin; j != jend; ++j ) {
        pat::Electron const * jelectron = dynamic_cast<pat::Electron const *>( j->masterClonePtr().get() );
        if ( jelectron != 0 )
            electrons->push_back( *jelectron );
    }

    for ( tau_iter jbegin = itaus.begin(), jend = itaus.end(), j = jbegin; j != jend; ++j ) {
        pat::Tau const * jtau = dynamic_cast<pat::Tau const *>( j->get() );
        if ( jtau != 0)
            taus->push_back( *jtau );
    }
    std::auto_ptr<float> puWeightPointer(new float);
    (*puWeightPointer) = filter_->puWeight();
    std::auto_ptr<int> genPVPointer(new int);
    (*genPVPointer) = filter_->genPV();
    std::auto_ptr<int> passTrigPointer(new int);
    (*passTrigPointer) = filter_->passTrig();
    event.put( jets, "jets");
    event.put( mets, "MET");
    event.put( muons, "muons");
    event.put( electrons, "electrons");
    event.put( taus, "taus" );
    event.put( puWeightPointer , "pileUp" );
    event.put( genPVPointer , "genpv" );
    event.put( passTrigPointer, "passTrig" );
    return passed; 
}


typedef edm::FilterWrapper<SHyFTSelector> EDSHyFTSelectorBase;
DEFINE_FWK_MODULE(EDSHyFTSelectorBase);
DEFINE_FWK_MODULE(EDSHyFTSelector);
