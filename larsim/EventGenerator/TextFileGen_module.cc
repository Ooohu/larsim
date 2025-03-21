/**
 * @file TextFileGen_module.cc
 * @brief Producer generating Monte Carlo truth record in LArSoft format from a text file
 * @date Mon Apr  8 09:20:02 2013
 * @author Brian Rebel
 */
/**
 * @class evgen::TextFileGen
 *
 *  This module assumes that the input file has the hepevt format for
 *  each event to be simulated.  See
 *
 *  http://cepa.fnal.gov/psm/simulation/mcgen/lund/pythia_manual/pythia6.3/pythia6301/node39.html
 *
 *  for details on the format.  In brief each event contains at least two
 *  lines.  The first line contains two entries, the event number (which is
 *  ignored in ART/LArSoft) and the number of particles in the event.  Each
 *  following line containes 15 entries to describe each particle.  The entries
 *  are:
 *
 *  1.  status code (should be set to 1 for any particle to be tracked, others
 *      won't be tracked)
 *  2.  the pdg code for the particle
 *  3.  the entry of the first mother for this particle in the event,
 *      0 means no mother
 *  4.  the entry of the second mother for this particle in the event,
 *      0 means no mother
 *  5. the entry of the first daughter for this particle in the event,
 *      0 means no daughter
 *  6. the entry of the second daughter for this particle in the event,
 *      0 means no daughter
 *  7. x component of the particle momentum
 *  8. y component of the particle momentum
 *  9. z component of the particle momentum
 *  10. energy of the particle
 *  11. mass of the particle
 *  12. x position of the particle initial position
 *  13. y position of the particle initial position
 *  14. z position of the particle initial position
 *  15. time of the particle production
 *
 *  For example, if you want to simulate a single muon with a 5 GeV energy
 *  moving only in the z direction, the entry would be
 *  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 *  0 1
 *  1 13 0 0 0 0 0. 0. 1.0 5.0011 0.105 1.0 1.0 1.0 0.0
 *  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 *  There are some assumptions that go into using this format that may not
 *  be obvious.  The first is that only particles with status code = 1
 *  are tracked in the LArSoft/Geant4 combination making the mother daughter
 *  relations somewhat irrelevant.  That also means that you should let
 *  Geant4 handle any decays.
 *
 *  The units in LArSoft are cm for distances and ns for time.
 *  The use of `TLorentzVector` below does not imply space and time have the same units
 *   (do not use `TLorentzVector::Boost()`).
 *
 *===========================================================================
 *
 * File and server input
 * ---------------------
 *
 * Orginally TextFileGen was written to read events from a hepevt format text file.
 * TextFileGen has now been modified to be able to read hepevt events from a server
 * using a specially formatted url.  The server should return one hepevt event per 
 * request.
 *
 * The choice between file and server input is controlled by fcl parameters
 * "InputFileName" and "InputURL."  One of these should be spceified, and the
 * other left blank or left out of the fcl configuraiton.
 *
 * This module does not make any assumption about how the url should be constructed,
 * except that it will read one hepevt event using using a single http(s) GET request.
 * Any parameters that need to be passed to the server should be embedded in the url
 * sent to the server via fcl parameter "InputURL."  This module does not modify
 * the specified url in any way.
 *
 * If the server is protected behind Fermilab SSO, set fcl parameter UseSSOAuth
 * to be true.  The remaining fcl parameters (Cert, Key, CertType, KeyType, and
 * KeyPasswd) give access to most of the relevant curl options for dealing with
 * certificates and keys.  When running interactively, these remaining fcl parameters
 * can be left as their default values if a certificate is obtained using command "kx509."
 * To directly use a p12 certficate obtained from cilogon.org, set the name of the 
 * certificate file using fcl parameter "Cert" or by setting an environment variable,
 * set "CertType" to "P12," and set "KeyPasswd" to the password you entered when you 
 * obtained the certificate from cilogon.org.  Note that grid proxies obtained using
 * "voms_proxy_init," including managed proxies, and proxies supplied to batch jobs,
 * do not (currently) work for accessing SSO-protected web pages.
 *
 * Http request statuses 5xx are treated as retriable indefinitely, up to some
 * maximum cumulative timeout specified by fcl parameter "Timeout."  Other http 
 * error statuses are treated as fatal errors, and will result in an exception being
 * thrown.
 *
 * The details of how to set up a server, and how to construct a server url, are
 * beyond the scope of this comment.  However, here are some hints.  An example
 * cgi server script can be found in larsim/scripts.  This script was installed on
 * the MicroBooNE web server for testing.  Here is a typical server url using this
 * installation.
 *
 * https://microboone-exp.fnal.gov/cgi-bin/hepevt.py?file=HEPevents.txt
 *
 * H. Greenlee, 17-May-2022
 *
 *===========================================================================
 *
 * FCL parameters.
 *
 * InputFileName - Name of hepevt input file (no default).
 * Offset        - Number of events to skip (for file input, default 0).
 * InputURL      - Server url (no default).
 * Timeout       - Maximum server cumulative timeout (seconds, default 7200, 0=none).
 * MoveY         - Propagate particles to y-plane (default don't propagate).
 * UseSSOAuth    - Use SSO Auth (certificat/proxy based, default false).
 * Cert          - Certificate file (defaults: $X509_USER_CERT, $X509_USER_PROXY, /tmp/x509up_u<uid>
 * Key           - Private key file (defaults: $X509_USER_KEY, $X509_USER_PROXY, /tmp/x509up_u<uid>
 * CertType      - Type of certificate file (default: libcurl decides).
 * KeyType       - Type of key file (default: libcurl decides).
 * KeyPasswd     - Key file password (default: none).
 *
 *===========================================================================
 */
#include <fstream>
#include <memory>
#include <string>
#include <utility>

#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "cetlib_except/exception.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "TLorentzVector.h"

#include "larcore/Geometry/Geometry.h"
#include "larcoreobj/SummaryData/RunData.h"
#include "nusimdata/SimulationBase/MCParticle.h"
#include "nusimdata/SimulationBase/MCTruth.h"

namespace evgen {
  class TextFileGen;
}

class evgen::TextFileGen : public art::EDProducer {
public:
  explicit TextFileGen(fhicl::ParameterSet const& p);

  void produce(art::Event& e) override;
  void beginJob() override;
  void beginRun(art::Run& run) override;

private:
  std::pair<unsigned, unsigned> readEventInfo(std::istream& is);
  simb::MCTruth readNextHepEvt();
  unsigned long int fOffset;
  std::ifstream* fInputFile;
  std::string fInputFileName; ///< Name of text file containing events to simulate
  double fMoveY;              ///< Project particles to a new y plane.
};

//------------------------------------------------------------------------------
evgen::TextFileGen::TextFileGen(fhicl::ParameterSet const& p)
  : EDProducer{p}
  , fOffset{p.get<unsigned long int>("Offset")}
  , fInputFile(0)
  , fInputFileName{p.get<std::string>("InputFileName")}
  , fMoveY{p.get<double>("MoveY", -1e9)}
{
  if (fMoveY > -1e8) {
    mf::LogWarning("TextFileGen") << "Particles will be moved to a new plane y = " << fMoveY
                                  << " cm.\n";
  }

  produces<std::vector<simb::MCTruth>>();
  produces<sumdata::RunData, art::InRun>();
}

//------------------------------------------------------------------------------
void evgen::TextFileGen::beginJob()
{
  if(!fInputFileName.empty()) {
    fInputFile = new std::ifstream(fInputFileName.c_str());

  // check that the file is a good one
  if (!fInputFile->good())
    throw cet::exception("TextFileGen")
      << "input text file " << fInputFileName << " cannot be read.\n";

  for (unsigned i = 0; i != fOffset; ++i) {
    auto const [eventNo, nparticles] = readEventInfo(*fInputFile);
    for (unsigned p = 0; p != nparticles; ++p) {
      constexpr auto all_chars_until = std::numeric_limits<unsigned>::max();
      fInputFile->ignore(all_chars_until, '\n');
    }
  }
}

//------------------------------------------------------------------------------
void evgen::TextFileGen::beginRun(art::Run& run)
{
  art::ServiceHandle<geo::Geometry const> geo;
  run.put(std::make_unique<sumdata::RunData>(geo->DetectorName()), art::fullRun());
}

//------------------------------------------------------------------------------
void evgen::TextFileGen::produce(art::Event& e)
{
  // check that the file is still good
  if (!fInputFile->good())
    throw cet::exception("TextFileGen")
      << "input text file " << fInputFileName << " cannot be read in produce().\n";

  //Now, read the Event to be used.

  // check that the file is still good
  if (!fInputFile->good())
    throw cet::exception("TextFileGen")
      << "input text file " << fInputFileName << " cannot be read in produce().\n";
  auto truthcol = std::make_unique<std::vector<simb::MCTruth>>();
  truthcol->push_back(readNextHepEvt());

  e.put(std::move(truthcol));
}

simb::MCTruth evgen::TextFileGen::readNextHepEvt()
{

  // declare the variables for reading in the event record
  int status = 0;
  int pdg = 0;
  int firstMother = 0;
  int secondMother = 0;
  int firstDaughter = 0;
  int secondDaughter = 0;
  double xMomentum = 0.;
  double yMomentum = 0.;
  double zMomentum = 0.;
  double energy = 0.;
  double mass = 0.;
  double xPosition = 0.;
  double yPosition = 0.;
  double zPosition = 0.;
  double time = 0.;



  // read in line to get event number and number of particles
  std::string oneLine;
  std::istringstream inputLine;
  simb::MCTruth nextEvent;
  auto const [eventNo, nParticles] = readEventInfo(*fInputFile);

  // now read in all the lines for the particles
  // in this interaction. only particles with
  // status = 1 get tracked in Geant4.
  for (unsigned short i = 0; i < nParticles; ++i) {
    std::getline(*fInputFile, oneLine);
    inputLine.clear();
    inputLine.str(oneLine);

    inputLine >> status >> pdg >> firstMother >> secondMother >> firstDaughter >> secondDaughter >>
      xMomentum >> yMomentum >> zMomentum >> energy >> mass >> xPosition >> yPosition >>
      zPosition >> time;



    //Project the particle to a new y plane
    if (fMoveY > -1e8) {
      double totmom = sqrt(pow(xMomentum, 2) + pow(yMomentum, 2) + pow(zMomentum, 2));
      double kx = xMomentum / totmom;
      double ky = yMomentum / totmom;
      double kz = zMomentum / totmom;
      if (ky) {
        double l = (fMoveY - yPosition) / ky;
        xPosition += kx * l;
        yPosition += ky * l;
        zPosition += kz * l;
      }
    }

    TLorentzVector pos(xPosition, yPosition, zPosition, time);
    TLorentzVector mom(xMomentum, yMomentum, zMomentum, energy);

    simb::MCParticle part(i, pdg, "primary", firstMother, mass, status);
    part.AddTrajectoryPoint(pos, mom);

    nextEvent.Add(part);

  } //  end loop on particles.

  return nextEvent;
}

std::pair<unsigned, unsigned> evgen::TextFileGen::readEventInfo(std::istream& iss)
{
  std::string line;
  getline(iss, line);
  std::istringstream buffer{line};

  // Parse read line for the event number and particles per event
  unsigned event, nparticles;
  buffer >> event >> nparticles;
  return {event, nparticles};
}


size_t evgen::TextFileGen::curl_callback(char* p, size_t size, size_t nmemb, void* userdata)
{
  // Add data to stringstram pointed to by userdata.
  
  std::stringstream* ss = reinterpret_cast<std::stringstream*>(userdata);
  ss->write(p, nmemb);

  // Done.

  return nmemb;
}






DEFINE_ART_MODULE(evgen::TextFileGen)
