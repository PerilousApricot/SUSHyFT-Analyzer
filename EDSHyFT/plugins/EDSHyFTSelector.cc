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
        if (ijet->hasTagInfo("secondaryVertexTagInfos") && ijet->tagInfoSecondaryVertex("secondaryVertexTagInfos")->nVertices() > 0) {
            jets->back().addUserFloat("secvtxMass",ijet->tagInfoSecondaryVertex("secondaryVertexTagInfos")->secondaryVertex(0).p4().mass());
        }
        if ( ijet != 0 ) {
            // match jet object with a (potential) Tau matching it
            for ( tau_iter jbegin = itaus.begin(), jend = itaus.end(), j = jbegin; j != jend; ++j ) {
                pat::Tau const * jtau = dynamic_cast<pat::Tau const *>( j->get() );
                if (jtau == 0) {
                    continue;
                }
                if (reco::deltaR( jtau->eta(), jtau->phi(),
                                  ijet->eta(), ijet->phi() ) < 0.3) {
                    jets->back().addUserFloat("tauJetPt", jtau->pt());
                    jets->back().addUserFloat("tauJetPhi", jtau->phi());
                    jets->back().addUserFloat("tauJetEta", jtau->eta());
                    jets->back().addUserFloat("tauJetMass", jtau->mass());
                    jets->back().addUserFloat("byLooseCombinedIsolationDeltaBetaCorr3Hits", jtau->tauID("byLooseCombinedIsolationDeltaBetaCorr3Hits") );
                    jets->back().addUserInt("byMediumCombinedIsolationDeltaBetaCorr3Hits", jtau->tauID("byMediumCombinedIsolationDeltaBetaCorr3Hits"));
                    jets->back().addUserInt("byTightCombinedIsolationDeltaBetaCorr3Hits", jtau->tauID("byTightCombinedIsolationDeltaBetaCorr3Hits"));
                }
            }
        }
    }

    for ( clone_iter jbegin = imuons.begin(), jend = imuons.end(), j = jbegin; j != jend; ++j ) {
        pat::Muon const * jmuon = dynamic_cast<pat::Muon const *>( j->masterClonePtr().get() );
        if ( jmuon != 0 )
            muons->push_back( *jmuon );
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
    event.put( jets, "jets");
    event.put( mets, "MET");
    event.put( muons, "muons");
    event.put( electrons, "electrons");
    event.put( taus, "taus" );
    event.put( puWeightPointer , "pileUp" );
    return passed; 
}


typedef edm::FilterWrapper<SHyFTSelector> EDSHyFTSelectorBase;
DEFINE_FWK_MODULE(EDSHyFTSelectorBase);
DEFINE_FWK_MODULE(EDSHyFTSelector);
