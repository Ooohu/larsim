////////////////////////////////////////////////////////////////////////
// Class:       MergeSimSources
// Module Type: producer
// File:        MergeSimSources_module.cc
//
// Generated at Tue Feb 17 12:16:35 2015 by Wesley Ketchum using artmod
// from cetpkgsupport v1_08_02.
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"

#include "canvas/Persistency/Common/Assns.h"
#include "canvas/Persistency/Common/FindOneP.h"
#include "lardata/Utilities/AssociationUtil.h"

#include <memory>

#include "nusimdata/SimulationBase/MCTruth.h"
#include "larsim/Simulation/LArG4Parameters.h"
#include "MergeSimSources.h"

namespace sim {
  class MergeSimSources;
}

class sim::MergeSimSources : public art::EDProducer {
public:
  explicit MergeSimSources(fhicl::ParameterSet const & p);
  // The destructor generated by the compiler is fine for classes
  // without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  MergeSimSources(MergeSimSources const &) = delete;
  MergeSimSources(MergeSimSources &&) = delete;
  MergeSimSources & operator = (MergeSimSources const &) = delete;
  MergeSimSources & operator = (MergeSimSources &&) = delete;

  // Required functions.
  void produce(art::Event & e) override;

private:

  std::vector<std::string> fInputSourcesLabels;
  std::vector<int>         fTrackIDOffsets;
  bool                     fUseLitePhotons;

  MergeSimSourcesUtility   fMergeUtility;

};


sim::MergeSimSources::MergeSimSources(fhicl::ParameterSet const & p)
  : EDProducer{p}
  , fInputSourcesLabels(p.get< std::vector<std::string> >("InputSourcesLabels"))
  , fTrackIDOffsets(p.get< std::vector<int> >("TrackIDOffsets"))
  , fMergeUtility(fTrackIDOffsets)
{

  if(fInputSourcesLabels.size() != fTrackIDOffsets.size())
    throw std::runtime_error("ERROR in MergeSimSources: Unequal input vector sizes.");

  art::ServiceHandle<sim::LArG4Parameters const> lgp;
  fUseLitePhotons = lgp->UseLitePhotons();
  
  if(!fUseLitePhotons) produces< std::vector<sim::SimPhotons>     >();
  else                 produces< std::vector<sim::SimPhotonsLite> >();
  
  produces< std::vector<simb::MCParticle> >();
  produces< std::vector<sim::SimChannel>  >();
  produces< std::vector<sim::AuxDetSimChannel> >();
  produces< art::Assns<simb::MCTruth, simb::MCParticle> >();
}

void sim::MergeSimSources::produce(art::Event & e)
{
  art::ServiceHandle<sim::LArG4Parameters const> lgp;

  std::unique_ptr< std::vector<simb::MCParticle> > partCol  (new std::vector<simb::MCParticle  >);
  std::unique_ptr< std::vector<sim::SimChannel>  > scCol    (new std::vector<sim::SimChannel>);
  std::unique_ptr< std::vector<sim::SimPhotons>  > PhotonCol(new std::vector<sim::SimPhotons>);
  std::unique_ptr< std::vector<sim::SimPhotonsLite>  > LitePhotonCol(new std::vector<sim::SimPhotonsLite>);
  std::unique_ptr< art::Assns<simb::MCTruth, simb::MCParticle> > tpassn(new art::Assns<simb::MCTruth, simb::MCParticle>);
  std::unique_ptr< std::vector< sim::AuxDetSimChannel > > adCol (new  std::vector<sim::AuxDetSimChannel> );

  fMergeUtility.Reset();
  
  for(size_t i_source=0; i_source<fInputSourcesLabels.size(); i_source++){

    std::string const& input_label = fInputSourcesLabels[i_source];

    art::Handle< std::vector<simb::MCParticle> > input_partCol;
    e.getByLabel(input_label,input_partCol);
    std::vector<simb::MCParticle> const& input_partColVector(*input_partCol);
    fMergeUtility.MergeMCParticles(*partCol,input_partColVector,i_source);

    const std::vector< std::vector<size_t> >& assocVectorPrimitive(fMergeUtility.GetMCParticleListMap());
    art::FindOneP<simb::MCTruth> mctAssn(input_partCol,e,input_label);
    for(size_t i_p=0; i_p<mctAssn.size(); i_p++)
      util::CreateAssn(*this,e,
		       *(partCol.get()),mctAssn.at(i_p),*(tpassn.get()),
		       assocVectorPrimitive[i_source][i_p]);
    
        
    art::Handle< std::vector<sim::SimChannel> > input_scCol;
    e.getByLabel(input_label,input_scCol);
    fMergeUtility.MergeSimChannels(*scCol,*input_scCol,i_source);

    art::Handle< std::vector<sim::AuxDetSimChannel> > input_adCol;
    e.getByLabel(input_label,input_adCol);
    fMergeUtility.MergeAuxDetSimChannels(*adCol,*input_adCol,i_source);

    if(!fUseLitePhotons){
      art::Handle< std::vector<sim::SimPhotons> > input_PhotonCol;
      e.getByLabel(input_label,input_PhotonCol);
      fMergeUtility.MergeSimPhotons(*PhotonCol,*input_PhotonCol);
    }
    else{
      art::Handle< std::vector<sim::SimPhotonsLite> > input_LitePhotonCol;
      e.getByLabel(input_label,input_LitePhotonCol);
      fMergeUtility.MergeSimPhotonsLite(*LitePhotonCol,*input_LitePhotonCol);
    }

    //truth-->particle assoc stuff here

  }

  e.put(std::move(partCol));
  e.put(std::move(scCol));
  e.put(std::move(adCol));
  if(!fUseLitePhotons) e.put(std::move(PhotonCol));
  else                 e.put(std::move(LitePhotonCol));
  e.put(std::move(tpassn));
  
}

DEFINE_ART_MODULE(sim::MergeSimSources)
