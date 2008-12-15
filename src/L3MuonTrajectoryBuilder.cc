/**
 *  Class: L3MuonTrajectoryBuilder
 *
 *  Description:
 *   Reconstruct muons starting
 *   from a muon track reconstructed
 *   in the standalone muon system (with DT, CSC and RPC
 *   information).
 *   It tries to reconstruct the corresponding
 *   track in the tracker and performs
 *   matching between the reconstructed tracks
 *   in the muon system and the tracker.
 *
 *
 *  $Date: 2008/12/12 21:20:00 $
 *  $Revision: 1.13 $
 *
 *  Authors :
 *  N. Neumeister            Purdue University
 *  C. Liu                   Purdue University
 *  A. Everett               Purdue University
 *  with contributions from: S. Lacaprara, J. Mumford, P. Traczyk
 *
 **/

#include "RecoMuon/L3TrackFinder/interface/L3MuonTrajectoryBuilder.h"


//---------------
// C++ Headers --
//---------------

#include <iostream>
#include <iomanip>
#include <algorithm>

//-------------------------------
// Collaborating Class Headers --
//-------------------------------

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/Event.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "DataFormats/TrackReco/interface/Track.h"

#include "RecoMuon/TrackingTools/interface/MuonCandidate.h"
#include "RecoMuon/TrackingTools/interface/MuonServiceProxy.h"
#include "RecoMuon/GlobalTrackingTools/interface/GlobalMuonTrackMatcher.h"

#include "FWCore/ServiceRegistry/interface/Service.h"

#include "RecoTracker/Record/interface/CkfComponentsRecord.h"
#include "RecoTracker/Record/interface/NavigationSchoolRecord.h"

#include "RecoTracker/TkTrackingRegions/interface/RectangularEtaPhiTrackingRegion.h"
#include "TrackingTools/TrajectoryCleaning/interface/TrajectoryCleanerBySharedHits.h"

#include "TrackingTools/DetLayers/interface/NavigationSetter.h"

using namespace std;
using namespace edm;

//----------------
// Constructors --
//----------------

L3MuonTrajectoryBuilder::L3MuonTrajectoryBuilder(const edm::ParameterSet& par,
							 const MuonServiceProxy* service) : GlobalTrajectoryBuilderBase(par, service) {

  theFirstEvent = true;
    
  theTkBuilderName = par.getParameter<std::string>("TkTrackBuilder");

  theTrajectoryCleaner = new TrajectoryCleanerBySharedHits();    

  theTkCollName = par.getParameter<edm::InputTag>("tkTrajLabel");

}


//--------------
// Destructor --
//--------------

L3MuonTrajectoryBuilder::~L3MuonTrajectoryBuilder() {
  if (theTrajectoryCleaner) delete theTrajectoryCleaner;
}

//
// get information from event
//
void L3MuonTrajectoryBuilder::setEvent(const edm::Event& event) {
  
  const std::string category = "Muon|RecoMuon|L3MuonTrajectoryBuilder|setEvent";
    
  GlobalTrajectoryBuilderBase::setEvent(event);
    
    
  if (theFirstEvent) {
      
    theFirstEvent = false;
    LogInfo(category) << "Constructing a Tk Trajectory Builder";
    GlobalTrajectoryBuilderBase::service()->eventSetup().get<CkfComponentsRecord>().get(theTkBuilderName,theTkBuilder);  
  }
    
  theTkBuilder->setEvent(event);
    
  // get tracker TrackCollection from Event
  edm::Handle<std::vector<Trajectory> > handleTrackerTrajs;
  event.getByLabel(theTkCollName,allTrackerTracks);
  LogInfo(category) 
      << "Found " << allTrackerTracks->size() 
      << " tracker Tracks with label "<< theTkCollName;  
  if (event.getByLabel(theTkCollName,handleTrackerTrajs) && event.getByLabel(theTkCollName,tkAssoMap)) {
    theTkTrajsAvailableFlag = true;
    allTrackerTrajs = &*handleTrackerTrajs;  
  }

  theTrajsAvailable = event.getByLabel(theTkCollName,theTkTrajCollection);
  LogDebug(category)<<"theTrajsAvailableFlag " << theTrajsAvailable ;
  theTkCandsAvailable = event.getByLabel(theTkCollName,theTkTrackCandCollection);
  
}

//
// reconstruct trajectories
//
MuonCandidate::CandidateContainer L3MuonTrajectoryBuilder::trajectories(const TrackCand& staCandIn) {

  const std::string category = "Muon|RecoMuon|L3MuonTrajectoryBuilder|trajectories";

  // cut on muons with low momenta
  if ( (staCandIn).second->pt() < thePtCut || (staCandIn).second->innerMomentum().Rho() < thePtCut || (staCandIn).second->innerMomentum().R() < 2.5 ) return CandidateContainer();

  // convert the STA track into a Trajectory if Trajectory not already present
  TrackCand staCand(staCandIn);
  addTraj(staCand);

  vector<TrackCand> trackerTracks;

  vector<TrackCand> regionalTkTracks = makeTkCandCollection(staCand);
  LogInfo(category) << "Found " << regionalTkTracks.size() << " tracks within region of interest";  
  
  // match tracker tracks to muon track
  trackerTracks = trackMatcher()->match(staCand, regionalTkTracks);
  
  LogInfo(category) << "Found " << trackerTracks.size() << " matching tracker tracks within region of interest";
  if ( trackerTracks.empty() ) return CandidateContainer();
  
  // build a combined tracker-muon MuonCandidate
  //
  // turn tkMatchedTracks into MuonCandidates
  //
  LogInfo(category) << "turn tkMatchedTracks into MuonCandidates";
  CandidateContainer tkTrajs;
  for (vector<TrackCand>::const_iterator tkt = trackerTracks.begin(); tkt != trackerTracks.end(); tkt++) {
    if ((*tkt).first != 0 && (*tkt).first->isValid()) {

      MuonCandidate* muonCand = new MuonCandidate( 0 ,staCand.second,(*tkt).second, new Trajectory(*(*tkt).first));
      tkTrajs.push_back(muonCand);
      LogTrace(category) << "tpush";

    }
  }
  
  if ( tkTrajs.empty() )  {
    LogInfo(category) << "tkTrajs empty";
    return CandidateContainer();
  }
  
  CandidateContainer result = build(staCand, tkTrajs);  
  LogInfo(category) << "Found "<< result.size() << " L3Muons from one L2Cand";

  // free memory
  if ( staCandIn.first == 0) delete staCand.first;

  for( CandidateContainer::const_iterator it = tkTrajs.begin(); it != tkTrajs.end(); ++it) {
    if ( (*it)->trajectory() ) delete (*it)->trajectory();
    if ( (*it)->trackerTrajectory() ) delete (*it)->trackerTrajectory();
    if ( *it ) delete (*it);
  }
  tkTrajs.clear();  

  for ( vector<TrackCand>::const_iterator is = regionalTkTracks.begin(); is != regionalTkTracks.end(); ++is) {
    delete (*is).first;   
  }
  
  return result;
  
}

//
// make a TrackCand collection using tracker Track, Trajectory information
//
vector<L3MuonTrajectoryBuilder::TrackCand> L3MuonTrajectoryBuilder::makeTkCandCollection(const TrackCand& staCand) {

  const std::string category = "Muon|RecoMuon|L3MuonTrajectoryBuilder|makeTkCandCollection";

  vector<TrackCand> tkCandColl;
  
  vector<TrackCand> tkTrackCands;
  
  if ( theTkTrajsAvailableFlag ) {
    for(TrajTrackAssociationCollection::const_iterator it = tkAssoMap->begin(); it != tkAssoMap->end(); ++it){	
      const Ref<vector<Trajectory> > traj = it->key;
      const reco::TrackRef tk = it->val;
      TrackCand tkCand = TrackCand(0,tk);
      //      if( traj->isValid() ) tkCand.first = &*traj ;
      if( traj->isValid() ) tkCand.first = new Trajectory(*traj) ;
      tkTrackCands.push_back(tkCand);
      //
      tkCandColl.push_back(tkCand);
    }
  } else {
    for ( unsigned int position = 0; position != allTrackerTracks->size(); ++position ) {
      LogDebug(category);
      reco::TrackRef tkTrackRef(allTrackerTracks,position);
      TrackCand tkCand = TrackCand(0,tkTrackRef);
      tkTrackCands.push_back(tkCand); 
      //
      tkCandColl.push_back(tkCand);
    }
  }
  
  //tkCandColl = chooseRegionalTrackerTracks(staCand,tkTrackCands);
  
  return tkCandColl;
  
}

