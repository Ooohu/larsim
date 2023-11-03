////////////////////////////////////////////////////////////////////////
/// \file  TwoBodyDecayGen_module.cc
/// \brief Generator for hypothetical particles that decay into two particles
///
/// Module designed to produce two random particles decay from one mother particle
/// Compared to SingleGen,  void SampleMany(simb::MCTruth &mct); is heavily modified.
///
/// \author keng.lin@rutgers.edu
////////////////////////////////////////////////////////////////////////
#ifndef EVGEN_TWOBODYGEN
#define EVGEN_TWOBODYGEN

// C++ includes.
#include <iostream>
#include <sstream>
#include <string>
#include <cmath>
#include <memory>
#include <iterator>
#include <vector>
#include <map>
#include <initializer_list>
#include <cctype> // std::tolower()


// Framework includes
#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Principal/Event.h"
#include "fhiclcpp/ParameterSet.h"
#include "fhiclcpp/types/Name.h"
#include "fhiclcpp/types/Comment.h"
#include "fhiclcpp/types/Atom.h"
#include "fhiclcpp/types/OptionalAtom.h"
#include "fhiclcpp/types/Sequence.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "art/Framework/Services/Optional/TFileDirectory.h"
#include "art/Framework/Services/Optional/RandomNumberGenerator.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "cetlib_except/exception.h"
#include "cetlib/exempt_ptr.h"
#include "cetlib/filesystem.h"
#include "cetlib/search_path.h"

// art extensions
#include "nutools/RandomUtils/NuRandomService.h"

// nutools includes
#include "nusimdata/SimulationBase/MCTruth.h"
#include "nusimdata/SimulationBase/MCParticle.h"
#include "nutools/EventGeneratorBase/evgenbase.h"

// lar includes
//#include "larcore/Geometry/Geometry.h"
//#include "larcoreobj/SummaryData/RunData.h"

#include "TVector3.h"
#include "TDatabasePDG.h"
#include "TMath.h"
#include "TFile.h"
#include "TH1.h"
#include "TH2.h"

#include "CLHEP/Random/RandFlat.h"
#include "CLHEP/Random/RandGaussQ.h"

namespace simb { class MCTruth; }

namespace evgen {

  /// module to produce single or multiple specified particles in the detector
  class TwoBodyDecayGen : public art::EDProducer {

  public:
    
    struct Config {
      using Name = fhicl::Name;
      using Comment = fhicl::Comment;
      
      fhicl::Atom<std::string> ParticleSelectionMode{
        Name("ParticleSelectionMode"),
        Comment("generate one particle, or one particle per PDG ID: " + presentOptions(ParticleSelectionModeNames))
      };
      
      fhicl::Atom<bool> PadOutVectors{
        Name("PadOutVectors"),
        Comment("if true, all per-PDG-ID quantities must contain only one value, which is then used for all PDG IDs")
      };
      
      fhicl::Sequence<int> PDG{
        Name("PDG"),
        Comment("PDG ID of the particles to be generated; this is the key for the other options marked as \"per PDG ID\"")
      };
      
      fhicl::Atom<std::string> PDist{
        Name("PDist"),
        Comment("momentum distribution type: " + presentOptions(DistributionNames)),
        optionName(kHIST, DistributionNames)
      };
      
      fhicl::Sequence<double> P0{
        Name("P0"),
        Comment("central momentum (GeV/c) to generate"),
//        [this](){ return !fromHistogram(PDist()); }//remove this helps removing "Unsupported parameters" error
      };
	  //CHECK
	  fhicl::Atom<std::string> MotherMassDist{
		  Name("MotherMassDist"),
			  Comment("mother mass distribution type: " + presentOptions(DistributionNames)),
			  optionName(kHIST, DistributionNames)
		};

	  fhicl::Sequence<double> MotherMass{Name("MotherMass"),
		  Comment("mother mass (GeV/c^2) to generate"),
//		  [this]() { return !fromHistogram(PDist()); }
		};

	  fhicl::Sequence<double> SigmaMotherMass{Name("SigmaMotherMass"),
		  Comment("variation in mother mass (GeV/c^2)"),
//		  [this]() { return !fromHistogram(PDist()); }
		};

	  fhicl::Atom<std::string> AngEXTDist{
		  Name("AngEXTDist"),
			  Comment("angular distribution type: " + presentOptions(DistributionNames)),
			  optionName(kHIST, DistributionNames)
		};

	  fhicl::Sequence<double> Theta0XZEXT{Name("Theta0XZEXT"),
		  Comment("Theta0XZEXT angle respected to mother particle momentum (degrees)")
		};

	  fhicl::Sequence<double> SigmaTheta0XZEXT{Name("SigmaTheta0XZEXT"),
		  Comment("variation in beta angle (degrees)")
		};

	  fhicl::Sequence<double> Theta0YZEXT{Name("Theta0YZEXT"),
		  Comment("Theta0YZEXT angle respected to mother particle momentum (degrees)")
		};

	  fhicl::Sequence<double> SigmaTheta0YZEXT{Name("SigmaTheta0YZEXT"),
		  Comment("variation in beta angle (degrees)")
		};
 //CHECK


      fhicl::Sequence<double> SigmaP{
        Name("SigmaP"),
        Comment("variation in momenta (GeV/c)"),
//        [this](){ return !fromHistogram(PDist()); }
      };
      
      fhicl::Sequence<double> X0{
        Name("X0"),
        Comment("central x position (cm) in world coordinates [per PDG ID]")
      };
      
      fhicl::Sequence<double> Y0{
        Name("Y0"),
        Comment("central y position (cm) in world coordinates [per PDG ID]")
      };
      
      fhicl::Sequence<double> Z0{
        Name("Z0"),
        Comment("central z position (cm) in world coordinates [per PDG ID]")
      };
      
      fhicl::Sequence<double> T0{
        Name("T0"),
        Comment("central time (s) [per PDG ID]")
      };
      
      fhicl::Sequence<double> SigmaX{
        Name("SigmaX"),
        Comment("variation (radius or RMS) in x position (cm) [per PDG ID]")
      };
      
      fhicl::Sequence<double> SigmaY{
        Name("SigmaY"),
        Comment("variation (radius or RMS) in y position (cm) [per PDG ID]")
      };
      
      fhicl::Sequence<double> SigmaZ{
        Name("SigmaZ"),
        Comment("variation (radius or RMS) in z position (cm) [per PDG ID]")
      };
      
      fhicl::Sequence<double> SigmaT{
        Name("SigmaT"),
        Comment("variation (semi-interval or RMS) in time (s) [per PDG ID]")
      };
      
      fhicl::Atom<std::string> PosDist{
        Name("PosDist"),
        Comment("distribution of starting position: " + presentOptions(DistributionNames, true, { kHIST }))
      };
      
      fhicl::Atom<std::string> TDist{
        Name("TDist"),
        Comment("time distribution type: " + presentOptions(DistributionNames, true, { kHIST }))
      };
      
      fhicl::Atom<bool> SingleVertex{
        Name("SingleVertex"),
        Comment("if true, all particles are produced at the same location"),
        false
      };
      
      fhicl::Sequence<double> Theta0XZ{
        Name("Theta0XZ"),
        Comment("angle from Z axis on world X-Z plane (degrees)")
      };
      
      fhicl::Sequence<double> Theta0YZ{
        Name("Theta0YZ"),
        Comment("angle from Z axis on world Y-Z plane (degrees)")
      };
      
      fhicl::Sequence<double> SigmaThetaXZ{
        Name("SigmaThetaXZ"),
        Comment("variation in angle in X-Z plane (degrees)")
      };
      
      fhicl::Sequence<double> SigmaThetaYZ{
        Name("SigmaThetaYZ"),
        Comment("variation in angle in Y-Z plane (degrees)")
      };
      
      fhicl::Atom<std::string> AngleDist{
        Name("AngleDist"),
        Comment("angular distribution type: " + presentOptions(DistributionNames)),
        optionName(kHIST, DistributionNames)
      };
      
      fhicl::Atom<std::string> HistogramFile{
        Name("HistogramFile"),
        Comment("ROOT file containing the required distributions for the generation"),
//        [this](){ return fromHistogram(AngleDist()) || fromHistogram(PDist()); }
      };
      
      fhicl::Sequence<std::string> PHist{
        Name("PHist"),
        Comment("name of the histograms of momentum distributions"),
//        [this](){ return fromHistogram(PDist()); }
      };
      
      /*
      fhicl::Sequence<std::string> ThetaPhiHist{
        Name("ThetaPhiHist"),
        Comment("name of the histograms of angular (theta/phi) distribution"),
        [this](){ return fromHistogram(AngleDist()); }
      };
      */
	  /*
      fhicl::Sequence<std::string> ThetaXzYzHist{
        Name("ThetaXzYzHist"),
        Comment("name of the histograms of angular (X-Z and Y-Z) distribution"),
        [this](){ return fromHistogram(AngleDist()); }
      };
	  */
      
      fhicl::OptionalAtom<rndm::NuRandomService::seed_t> Seed{
        Name("Seed"),
        Comment("override the random number generator seed")
      };
      
      
        private:
      
      /// Returns whether the specified mode is an histogram distribution.
      bool fromHistogram(std::string const& key) const;
      
    }; // struct Config
    
    
    using Parameters = art::EDProducer::Table<Config>;
    
    
    explicit TwoBodyDecayGen(Parameters const& config);

    // This is called for each event.
    void produce(art::Event& evt);
    void beginRun(art::Run& run);

  private:

    /// Names of all particle selection modes.
    static const std::map<int, std::string> ParticleSelectionModeNames;
    /// Names of all distribution modes.
    static const std::map<int, std::string> DistributionNames;
    
    void SampleOne(unsigned int   i, 
		   simb::MCTruth &mct);        
    void SampleMany(simb::MCTruth &mct);        
    void Sample(simb::MCTruth &mct);        
    void printVecs(std::vector<std::string> const& list);
    bool PadVector(std::vector<double> &vec);      
    double SelectFromHist(const TH1& h);
    void SelectFromHist(const TH2& h, double &x, double &y);
    
    /// @{
    /// @name Constants for particle type extraction mode (`ParticleSelectionMode` parameter).
    
    static constexpr int kSelectAllParts    = 0; ///< One particle per entry is generated
    static constexpr int kSelectOneRandPart = 1; ///< One particle is generated, extracted from the provided options.
    /// @}
    
    /// @{
    /// @name Constants for kinematic distribution options.
    
    static constexpr int kUNIF = 0;    ///< Uniform distribution.
    static constexpr int kGAUS = 1;    ///< Gaussian distribution.
    static constexpr int kHIST = 2;    ///< Distribution from histograms.
    /// @}

    int                 fMode;           ///< Particle Selection Mode 
                                         ///< 0--generate a list of all particles, 
                                         ///< 1--generate a single particle selected randomly from the list
    bool                fPadOutVectors;  ///< Select to pad out configuration vectors if they are not of 
					 ///< of the same length as PDG  
                                         ///< false: don't pad out - all values need to specified
                                         ///< true: pad out - default values assumed and printed out
    std::vector<int>    fPDG;            ///< PDG code of particles to generate    
    std::vector<double> fP0;             ///< Central momentum (GeV/c) to generate    
    std::vector<double> fSigmaP;         ///< Variation in momenta (GeV/c)    
    int                 fPDist;          ///< How to distribute momenta (gaus or uniform)    
    std::vector<double> fX0;             ///< Central x position (cm) in world coordinates 
    std::vector<double> fY0;             ///< Central y position (cm) in world coordinates
    std::vector<double> fZ0;             ///< Central z position (cm) in world coordinates
    std::vector<double> fT0;             ///< Central t position (s) in world coordinates
    std::vector<double> fSigmaX;         ///< Variation in x position (cm)    
    std::vector<double> fSigmaY;         ///< Variation in y position (cm)    
    std::vector<double> fSigmaZ;         ///< Variation in z position (cm)    
    std::vector<double> fSigmaT;         ///< Variation in t position (s)    
    int                 fPosDist;        ///< How to distribute xyz (gaus, or uniform)        
    int                 fTDist;          ///< How to distribute t  (gaus, or uniform)        
    bool                fSingleVertex;   ///< if true - all particles produced at the same location        
    std::vector<double> fTheta0XZ;       ///< Angle in XZ plane (degrees)    
    std::vector<double> fTheta0YZ;       ///< Angle in YZ plane (degrees)    
    std::vector<double> fSigmaThetaXZ;   ///< Variation in angle in XZ plane    
    std::vector<double> fSigmaThetaYZ;   ///< Variation in angle in YZ plane    
    int                 fAngleDist;      ///< How to distribute angles (gaus, uniform)

//CHECK
	std::vector<double> fMotherMass;
	std::vector<double> fSigmaMotherMass;
	int                 fMotherMassDist;
	std::vector<double> fTheta0XZEXT;
	std::vector<double> fSigmaTheta0XZEXT;
	int                 fAngEXTDist;
	std::vector<double> fTheta0YZEXT;
	std::vector<double> fSigmaTheta0YZEXT;
//CHECK

    std::string fHistFileName;               ///< Filename containing histogram of momenta
    std::vector<std::string> fPHist;     ///< name of histogram of momenta
//    std::vector<std::string> fThetaPhiHist; ///< name of histogram for theta/phi distribution
    std::vector<std::string> fThetaXzYzHist;   ///< name of histogram for thetaxz/thetayz distribution

    std::vector<std::unique_ptr<TH1>> hPHist ;     /// actual TH1 for momentum distributions
//    std::vector<TH2*> hThetaPhiHist ;  /// actual TH1 for theta distributions - Theta on x axis
    std::vector<std::unique_ptr<TH2>> hThetaXzYzHist ; /// actual TH2 for angle distributions - Xz on x axis . 
    // FYI - thetaxz and thetayz are related to standard polar angles as follows:
    // thetaxz = atan2(math.sin(theta) * cos(phi), cos(theta))
    // thetayz = asin(sin(theta) * sin(phi));
    
    cet::exempt_ptr<CLHEP::HepRandomEngine> fEngine; // FIXME: This should be a reference.

    
    /// Returns a vector with the name of particle selection mode keywords.
    static std::map<int, std::string> makeParticleSelectionModeNames();
    
    /// Returns a vector with the name of distribution keywords.
    static std::map<int, std::string> makeDistributionNames();
    
    
    /// Performs checks and initialization based on the current configuration.
    void setup();
    
    /**
     * @brief Parses an option string and returns the corresponding option number.
     * @tparam OptionList type of list of options (e.g. `std::map<int, std::string>`)
     * @param Option the string of the option to be parsed
     * @param allowedOptions list of valid options, as key/name pairs
     * @return the key of the `Option` string from `allowedOptions`
     * @throws std::runtime_error if `Option` is not in the option list
     * 
     * The option string `Option` represent a single one among the supported
     * options as defined in `allowedOptions`. The option string can be either
     * one of the option names (the matching is not case-sensitive) or the
     * number of the option itself.
     * 
     * `OptionList` requirements
     * --------------------------
     * 
     * `OptionList` must behave like a sequence with forward iterators.
     * Each element must behave as a pair, whose first element is the option key
     * and the second element is the option name, equivalent to a string in that
     * it must be forward-iterable and its elements can be converted by
     * `std::tolower()`. The key type has no requirements beside being copiable.
     */
    template <typename OptionList>
    static auto selectOption
      (std::string Option, OptionList const& allowedOptions) -> decltype(auto);
    
    /**
     * @brief Returns a string describing all options in the list
     * @tparam OptionList type of list of options (e.g. `std::map<int, std::string>`)
     * @param allowedOptions the list of allowed options
     * @param printKey whether to print the key of the option beside its name
     * @param excludeKeys list of keys to be ignored (none by default)
     * @return a string with all options in a line
     * 
     * The result string is a list of option names, separated by commas, like in
     * `"'apple', 'orange', 'banana'"`. If `printKey` is `true`, the key of each
     * option is also written in parentheses, like in
     * `"'apple' (1), 'orange' (7), 'banana' (2)"`.
     */
    template <typename OptionList>
    static std::string presentOptions(
      OptionList const& allowedOptions, bool printKey,
      std::initializer_list<typename OptionList::value_type::first_type> exclude
      );
    
    template <typename OptionList>
    static std::string presentOptions
      (OptionList const& allowedOptions, bool printKey = true)
      { return presentOptions(allowedOptions, printKey, {}); }
    
    
    /// Returns the name of the specified option key, or `defName` if not known.
    template <typename OptionList>
    static std::string optionName(
      typename OptionList::value_type::first_type optionKey,
      OptionList const& allowedOptions,
      std::string defName = "<unknown>"
      );
    
  }; // class TwoBodyDecayGen
}

namespace evgen{

  std::map<int, std::string> TwoBodyDecayGen::makeParticleSelectionModeNames() {
    std::map<int, std::string> names;
    names[int(kSelectAllParts   )] = "all";
    names[int(kSelectOneRandPart)] = "singleRandom";
    return names;
  } // TwoBodyDecayGen::makeParticleSelectionModeNames()
  
  std::map<int, std::string> TwoBodyDecayGen::makeDistributionNames() {
    std::map<int, std::string> names;
    names[int(kUNIF)] = "uniform";
    names[int(kGAUS)] = "Gaussian";
    names[int(kHIST)] = "histograms";
    return names;
  } // TwoBodyDecayGen::makeDistributionNames()
  
  const std::map<int, std::string> TwoBodyDecayGen::ParticleSelectionModeNames
    = TwoBodyDecayGen::makeParticleSelectionModeNames();
  const std::map<int, std::string> TwoBodyDecayGen::DistributionNames
    = TwoBodyDecayGen::makeDistributionNames();
  
  
  template <typename OptionList>
  auto TwoBodyDecayGen::selectOption
    (std::string Option, OptionList const& allowedOptions) -> decltype(auto)
  {
    using key_type = typename OptionList::value_type::first_type;
    using tolower_type = int(*)(int);
    auto toLower = [](auto const& S)
      {
        std::string s;
        s.reserve(S.size());
        std::transform(S.cbegin(), S.cend(), std::back_inserter(s),
          (tolower_type) &std::tolower);
        return s;
      };
    auto option = toLower(Option);
    for (auto const& candidate: allowedOptions) {
      if (toLower(candidate.second) == option) return candidate.first;
    }
    try {
      std::size_t end;
      key_type num = std::stoi(Option, &end);
      if (allowedOptions.count(num) && (end == Option.length())) return num;
    }
    catch (std::invalid_argument const&) {}
    throw std::runtime_error("Option '" + Option + "' not supported.");
  } // TwoBodyDecayGen::selectOption()
  
  
  template <typename OptionList>
  std::string TwoBodyDecayGen::presentOptions(
    OptionList const& allowedOptions, bool printKey /* = true */,
    std::initializer_list<typename OptionList::value_type::first_type> exclude /* = {} */
  ) {
    std::string msg;
    
    unsigned int n = 0;
    for (auto const& option: allowedOptions) {
      auto const& key = option.first;
      if (std::find(exclude.begin(), exclude.end(), key) != exclude.end())
        continue;
      if (n++ > 0) msg += ", ";
      msg += '\"' + std::string(option.second) + '\"';
      if (printKey)
        msg += " (" + std::to_string(key) + ")";
    } // for
    return msg;
  } // TwoBodyDecayGen::presentOptions()
  
  
  template <typename OptionList>
  std::string TwoBodyDecayGen::optionName(
    typename OptionList::value_type::first_type optionKey,
    OptionList const& allowedOptions,
    std::string defName /* = "<unknown>" */
  ) {
    auto iOption = allowedOptions.find(optionKey);
    return (iOption != allowedOptions.end())? iOption->second: defName;
  } // TwoBodyDecayGen::optionName()
  
  
  //____________________________________________________________________________
  bool TwoBodyDecayGen::Config::fromHistogram(std::string const& key) const {
    return selectOption(PDist(), DistributionNames) == kHIST;
  } // TwoBodyDecayGen::Config::fromHistogram()
  
  //____________________________________________________________________________
  TwoBodyDecayGen::TwoBodyDecayGen(Parameters const& config)
    : EDProducer{config}
    , fMode         (selectOption(config().ParticleSelectionMode(), ParticleSelectionModeNames))
    , fPadOutVectors(config().PadOutVectors())
    , fPDG          (config().PDG())
    , fP0           (config().P0())
    , fSigmaP       (config().SigmaP())
    , fPDist        (selectOption(config().PDist(), DistributionNames))
    , fX0           (config().X0())
    , fY0           (config().Y0())
    , fZ0           (config().Z0())
    , fT0           (config().T0())
    , fSigmaX       (config().SigmaX())
    , fSigmaY       (config().SigmaY())
    , fSigmaZ       (config().SigmaZ())
    , fSigmaT       (config().SigmaT())
    , fPosDist      (selectOption(config().PosDist(), DistributionNames))
    , fTDist        (selectOption(config().TDist(), DistributionNames))
    , fSingleVertex (config().SingleVertex())
    , fTheta0XZ     (config().Theta0XZ())
    , fTheta0YZ     (config().Theta0YZ())
    , fSigmaThetaXZ (config().SigmaThetaXZ())
    , fSigmaThetaYZ (config().SigmaThetaYZ())
    , fAngleDist    (selectOption(config().AngleDist(), DistributionNames))
//CHECK
    , fMotherMass(config().MotherMass())
    , fSigmaMotherMass(config().SigmaMotherMass())
    , fMotherMassDist(selectOption(config().MotherMassDist(), DistributionNames))
    , fTheta0XZEXT(config().Theta0XZEXT())
    , fSigmaTheta0XZEXT(config().SigmaTheta0XZEXT())
    , fAngEXTDist(selectOption(config().AngEXTDist(), DistributionNames))
    , fTheta0YZEXT(config().Theta0YZEXT())
    , fSigmaTheta0YZEXT(config().SigmaTheta0YZEXT())
//CHECK
    , fHistFileName (config().HistogramFile())
    , fPHist        (config().PHist())
//    , fThetaPhiHist (config().ThetaPhiHist())
//    , fThetaXzYzHist(config().ThetaXzYzHist())
  {
    setup();
    
    // create a default random engine; obtain the random seed from NuRandomService,
    // unless overridden in configuration with key "Seed"
    (void)art::ServiceHandle<rndm::NuRandomService>()->createEngine(*this);
    art::ServiceHandle<art::RandomNumberGenerator> rng;
    auto& engine = rng->getEngine(art::ScheduleID::first(),
                                  config.get_PSet().get<std::string>("module_label"));
    fEngine = cet::make_exempt_ptr(&engine);
    rndm::NuRandomService::seed_t seed;
    if (config().Seed(seed)) {
      fEngine->setSeed(seed, 0 /* dummy? */);
    }

    produces< std::vector<simb::MCTruth> >();
    //    produces< sumdata::RunData, art::InRun >();

  }
  
  
  //____________________________________________________________________________
  void TwoBodyDecayGen::setup()
  {
    // do not put seed in reconfigure because we don't want to reset 
    // the seed midstream
    std::vector<std::string> vlist(21);
    vlist[0]  = "PDG";
    vlist[1]  = "P0";
    vlist[2]  = "SigmaP";
    vlist[3]  = "X0";
    vlist[4]  = "Y0";
    vlist[5]  = "Z0";
    vlist[6]  = "SigmaX";
    vlist[7]  = "SigmaY";
    vlist[8]  = "SigmaZ";
    vlist[9]  = "Theta0XZ";
    vlist[10] = "Theta0YZ";
    vlist[11] = "SigmaThetaXZ";
    vlist[12] = "SigmaThetaYZ";
    vlist[13] = "T0";
    vlist[14] = "SigmaT";

	//CHECK
	vlist[15] = "MotherMass";             
	vlist[16] = "SigmaMotherMass";        
	vlist[17] = "Theta0XZEXT";           
	vlist[18] = "SigmaTheta0XZEXT";      
	vlist[19] = "Theta0YZEXT";     
	vlist[20] = "SigmaTheta0YZEXT";


//    vlist[15] = "PHist";
//    vlist[16] = "ThetaHist";
//    vlist[17] = "PhiHist";
    
    // begin tests for multiple particle error possibilities  
    std::string list;
    if (fPDist != kHIST) {
      if( !this->PadVector(fP0            ) ){ list.append(vlist[1].append(", \n")); }
      if( !this->PadVector(fSigmaP        ) ){ list.append(vlist[2].append(", \n")); }
    }
    if( !this->PadVector(fX0              ) ){ list.append(vlist[3].append(", \n")); }
    if( !this->PadVector(fY0              ) ){ list.append(vlist[4].append(", \n")); }
    if( !this->PadVector(fZ0              ) ){ list.append(vlist[5].append(", \n")); }
    if( !this->PadVector(fSigmaX          ) ){ list.append(vlist[6].append(", \n")); }
    if( !this->PadVector(fSigmaY          ) ){ list.append(vlist[7].append(", \n")); }
    if( !this->PadVector(fSigmaZ          ) ){ list.append(vlist[8].append(", \n")); }
    if( !this->PadVector(fTheta0XZ        ) ){ list.append(vlist[9].append(", \n")); }
    if( !this->PadVector(fTheta0YZ        ) ){ list.append(vlist[10].append(", \n")); }
    if( !this->PadVector(fSigmaThetaXZ    ) ){ list.append(vlist[11].append(", \n")); }
    if( !this->PadVector(fSigmaThetaYZ    ) ){ list.append(vlist[12].append("  \n")); }
    if( !this->PadVector(fT0              ) ){ list.append(vlist[13].append(", \n")); }
    if( !this->PadVector(fSigmaT          ) ){ list.append(vlist[14].append(", \n")); }
//CHECK
    if( !this->PadVector(fMotherMass              ) ){ list.append(vlist[15].append(", \n")); }
    if( !this->PadVector(fSigmaMotherMass         ) ){ list.append(vlist[16].append(", \n")); }
    if( !this->PadVector(fTheta0XZEXT            ) ){ list.append(vlist[17].append(", \n")); }
    if( !this->PadVector(fSigmaTheta0XZEXT       ) ){ list.append(vlist[18].append(", \n")); }
    if( !this->PadVector(fTheta0YZEXT      ) ){ list.append(vlist[19].append(", \n")); }
    if( !this->PadVector(fSigmaTheta0YZEXT ) ){ list.append(vlist[20].append(", \n")); }
    

    if(list.size() > 0)
      throw cet::exception("TwoBodyDecayGen") << "The "<< list 
					<< "\n vector(s) defined in the fhicl files has/have "
					<< "a different size than the PDG vector "
					<< "\n and it has (they have) more than one value, "
					<< "\n disallowing sensible padding "
					<< " and/or you have set fPadOutVectors to false. \n";
    
    if(fPDG.size() > 1 && fPadOutVectors) this->printVecs(vlist);

    // If needed, get histograms for momentum and angle distributions
    TFile* histFile = nullptr;
    if (!fHistFileName.empty()) {
      if (fHistFileName[0] == '/') {
        // We have an absolute path, use given name exactly.
        if (cet::file_exists(fHistFileName)) {
          histFile = new TFile(fHistFileName.c_str());
          if (!histFile || histFile->IsZombie() || !histFile->IsOpen()) {
            delete histFile;
            histFile = nullptr;
            throw art::Exception(art::errors::NotFound) << "Cannot open ROOT file specified in parameter HistogramFile: \"" << fHistFileName << "\"";
          }
        }
        else {
          throw art::Exception(art::errors::NotFound) << "ROOT file specified in parameter HistogramFile: \"" << fHistFileName << "\" does not exist!";
        }
      }
      else {
        // We have a relative path, search starting from current directory.
        std::string relative_filename{"./"};
        relative_filename += fHistFileName;
        if (cet::file_exists(relative_filename)) {
          histFile = new TFile(relative_filename.c_str());
          if (!histFile || histFile->IsZombie() || !histFile->IsOpen()) {
            delete histFile;
            histFile = nullptr;
            throw art::Exception(art::errors::NotFound) << "Cannot open ROOT file found using relative path and originally specified in parameter HistogramFile: \"" << relative_filename << '"';
          }
        }
        else {
          cet::search_path sp{"FW_SEARCH_PATH"};
          std::string found_filename;
          auto found = sp.find_file(fHistFileName, found_filename);
          if (!found) {
            throw art::Exception(art::errors::NotFound) << "Cannot find ROOT file in current directory nor on FW_SEARCH_PATH specified in parameter HistogramFile: \"" << fHistFileName << '"';
          }
          histFile = new TFile(found_filename.c_str());
          if (!histFile || histFile->IsZombie() || !histFile->IsOpen()) {
            delete histFile;
            histFile = nullptr;
            throw art::Exception(art::errors::NotFound) << "Cannot open ROOT file found on FW_SEARCH_PATH and originally specified in parameter HistogramFile: \"" << found_filename << '"';
          }
        }
      }
    }
    
    //
    // deal with position distribution
    //
    switch (fPosDist) {
      case kGAUS: case kUNIF: break; // supported, no further action needed
      default:
        throw art::Exception(art::errors::Configuration)
          << "Position distribution of type '"
          << optionName(fPosDist, DistributionNames)
          << "' (" << std::to_string(fPosDist) << ") is not supported.";
    } // switch(fPosDist)
    
    //
    // deal with time distribution
    //
    switch (fTDist) {
      case kGAUS: case kUNIF: break; // supported, no further action needed
      default:
        throw art::Exception(art::errors::Configuration)
          << "Time distribution of type '"
          << optionName(fTDist, DistributionNames)
          << "' (" << std::to_string(fTDist) << ") is not supported.";
    } // switch(fTDist)
    
    //
    // deal with momentum distribution
    //
    switch (fPDist) {
      case kHIST:
        //if (fPHist.size() != fPDG.size()) {
        //  throw art::Exception(art::errors::Configuration)
        //    << fPHist.size() << " momentum histograms to describe " << fPDG.size() << " particle types...";
        //}
        hPHist.reserve(fPHist.size());
        for (auto const& histName: fPHist) {
          TH1* pHist = dynamic_cast<TH1*>(histFile->Get(histName.c_str()));
          if (!pHist) {
            throw art::Exception(art::errors::NotFound)
             << "Failed to read momentum histogram '" << histName << "' from '" << histFile->GetPath() << "\'";
          }
          pHist->SetDirectory(nullptr); // make it independent of the input file
          hPHist.emplace_back(pHist);
        } // for
        break;
      default: // supported, no further action needed
        break;
    } // switch(fPDist)
    
    switch (fAngleDist) {
      case kHIST:
        if (fThetaXzYzHist.size() != fPDG.size()) {
          throw art::Exception(art::errors::Configuration)
            << fThetaXzYzHist.size() << " direction histograms to describe " << fPDG.size() << " particle types...";
        }
        hThetaXzYzHist.reserve(fThetaXzYzHist.size());
        for (auto const& histName: fThetaXzYzHist) {
          TH2* pHist = dynamic_cast<TH2*>(histFile->Get(histName.c_str()));
          if (!pHist) {
            throw art::Exception(art::errors::NotFound)
             << "Failed to read direction histogram '" << histName << "' from '" << histFile->GetPath() << "\'";
          }
          pHist->SetDirectory(nullptr); // make it independent of the input file
          hThetaXzYzHist.emplace_back(pHist);
        } // for
      default: // supported, no further action needed
        break;
    } // switch(fAngleDist)
    
    delete histFile;
    
  }

  //____________________________________________________________________________
  bool TwoBodyDecayGen::PadVector(std::vector<double> &vec)
  {
    // check if the vec has the same size as fPDG
    if( vec.size() != fPDG.size() ){
      // if not padding out the vectors always cause an 
      // exception to be thrown if the vector in question
      // is not the same size as the fPDG vector
      // the exception is thrown in the reconfigure method
      // that calls this one
      if     (!fPadOutVectors) return false;
      else if( fPadOutVectors){
	// if padding of vectors is desired but the vector in
	// question has more than one entry it isn't clear
	// what the padded values should be so cause
	// an exception
	if(vec.size() != 1) return false;

	// pad it out
	vec.resize(fPDG.size(), vec[0]);

      }// end if padding out vectors
    }// end if the vector size is not the same as fPDG

    return true;
  }

  //____________________________________________________________________________
  void TwoBodyDecayGen::beginRun(art::Run& run)
  {

    // grab the geometry object to see what geometry we are using
    //    art::ServiceHandle<geo::Geometry> geo;
    // std::unique_ptr<sumdata::RunData> runcol(new sumdata::RunData(geo->DetectorName()));

    //run.put(std::move(runcol));

    return;
  }

  //____________________________________________________________________________
  void TwoBodyDecayGen::produce(art::Event& evt)
  {

    ///unique_ptr allows ownership to be transferred to the art::Event after the put statement
    std::unique_ptr< std::vector<simb::MCTruth> > truthcol(new std::vector<simb::MCTruth>);

    simb::MCTruth truth;
    truth.SetOrigin(simb::kSingleParticle);
    Sample(truth);

    MF_LOG_DEBUG("TwoBodyDecayGen") << truth;

    truthcol->push_back(truth);

    evt.put(std::move(truthcol));

    return;
  }

  //____________________________________________________________________________
  // Draw the type, momentum and position of a single particle from the 
  // FCIHL description
  void TwoBodyDecayGen::SampleOne(unsigned int i, simb::MCTruth &mct){

    CLHEP::RandFlat   flat(*fEngine);
    CLHEP::RandGaussQ gauss(*fEngine);

    // Choose momentum
    double p = 0.0;
    double m = 0.0;
    if (fPDist == kGAUS) {
      p = gauss.fire(fP0[i], fSigmaP[i]);
    }
    else if (fPDist == kHIST){
      p = SelectFromHist(*(hPHist[i]));
    }
    else if (fPDist == kUNIF){
      p = fP0[i] + fSigmaP[i]*(2.0*flat.fire()-1.0);
    }
//    else {std::cout << "do not understand the value of PDist!";}

    static TDatabasePDG  pdgt;
    TParticlePDG* pdgp = pdgt.GetParticle(fPDG[i]);
    if (pdgp) m = pdgp->Mass();
    
    // Choose position
    TVector3 x;
    if (fPosDist == kGAUS) {
      x[0] = gauss.fire(fX0[i], fSigmaX[i]);;
      x[1] = gauss.fire(fY0[i], fSigmaY[i]);
      x[2] = gauss.fire(fZ0[i], fSigmaZ[i]);
    }
    else {
      x[0] = fX0[i] + fSigmaX[i]*(2.0*flat.fire()-1.0);
      x[1] = fY0[i] + fSigmaY[i]*(2.0*flat.fire()-1.0);
      x[2] = fZ0[i] + fSigmaZ[i]*(2.0*flat.fire()-1.0);
    }

    double t = 0.;
    if(fTDist==kGAUS){
      t = gauss.fire(fT0[i], fSigmaT[i]);
    }
    else{
      t = fT0[i] + fSigmaT[i]*(2.0*flat.fire()-1.0);
    }

    TLorentzVector pos(x[0], x[1], x[2], t);
    
    // Choose angles
    double thxz = 0;
    double thyz = 0;
    
    double thyzrads = 0; 
    double thyzradsplussigma = 0;
    double thyzradsminussigma = 0;

    if (fAngleDist == kGAUS) {
      thxz = gauss.fire(fTheta0XZ[i], fSigmaThetaXZ[i]);
      thyz = gauss.fire(fTheta0YZ[i], fSigmaThetaYZ[i]);
    }
    else if (fAngleDist == kHIST){ // Select thetaxz and thetayz from histogram
      double thetaxz = 0;
      double thetayz = 0;
      SelectFromHist(*(hThetaXzYzHist[i]), thetaxz, thetayz);
      thxz = (180./M_PI)*thetaxz;
      thyz = (180./M_PI)*thetayz;
    }
    else {
      
      // Choose angles flat in phase space, which is flat in theta_xz 
      // and flat in sin(theta_yz).
   
      thxz = fTheta0XZ[i] + fSigmaThetaXZ[i]*(2.0*flat.fire()-1.0);
     
      thyzrads = std::asin(std::sin((M_PI/180.)*(fTheta0YZ[i]))); //Taking asin of sin gives value between -Pi/2 and Pi/2 regardless of user input
      thyzradsplussigma = TMath::Min((thyzrads + ((M_PI/180.)*fabs(fSigmaThetaYZ[i]))), M_PI/2.);
      thyzradsminussigma = TMath::Max((thyzrads - ((M_PI/180.)*fabs(fSigmaThetaYZ[i]))), -M_PI/2.);
         
      //uncomment line to print angular variation info
      //std::cout << "Central angle: " << (180./M_PI)*thyzrads << " Max angle: " << (180./M_PI)*thyzradsplussigma << " Min angle: " << (180./M_PI)*thyzradsminussigma << std::endl; 

      double sinthyzmin = std::sin(thyzradsminussigma);
      double sinthyzmax = std::sin(thyzradsplussigma);
      double sinthyz = sinthyzmin + flat.fire() * (sinthyzmax - sinthyzmin);
      thyz = (180. / M_PI) * std::asin(sinthyz);
    }
    
    double thxzrad=thxz*M_PI/180.0;	
    double thyzrad=thyz*M_PI/180.0;

    TLorentzVector pvec(p*std::cos(thyzrad)*std::sin(thxzrad),
			p*std::sin(thyzrad),
			p*std::cos(thxzrad)*std::cos(thyzrad),
			std::sqrt(p*p+m*m));
 
    // set track id to -i as these are all primary particles and have id <= 0
    int trackid = -1*(i+1);
    std::string primary("primary");

    simb::MCParticle part(trackid, fPDG[i], primary);
    part.AddTrajectoryPoint(pos, pvec);

    //std::cout << "Px: " <<  pvec.Px() << " Py: " << pvec.Py() << " Pz: " << pvec.Pz() << std::endl;
    //std::cout << "x: " <<  pos.X() << " y: " << pos.Y() << " z: " << pos.Z() << " time: " << pos.T() << std::endl;
    //std::cout << "YZ Angle: " << (thyzrad * (180./M_PI)) << " XZ Angle: " << (thxzrad * (180./M_PI)) << std::endl; 
     
    mct.Add(part);
  }

  //____________________________________________________________________________
  // Draw the type, momentum and position for all particles from the 
  // FCIHL description.  Start positions will all match but momenta and angles drawn from
  // distributions defined in the fhicls
	void TwoBodyDecayGen::SampleMany(simb::MCTruth &mct){

		bool fverbose = true;

		//CHECK, modify this one! 
		CLHEP::RandFlat   flat(*fEngine);
		CLHEP::RandGaussQ gauss(*fEngine);

		//Choose position, all particles come from the same position
		TVector3 x;
		if (fPosDist == kGAUS) {
			x[0] = gauss.fire(fX0[0], fSigmaX[0]);;
			x[1] = gauss.fire(fY0[0], fSigmaY[0]);
			x[2] = gauss.fire(fZ0[0], fSigmaZ[0]);
		}
		else {
			x[0] = fX0[0] + fSigmaX[0]*(2.0*flat.fire()-1.0);
			x[1] = fY0[0] + fSigmaY[0]*(2.0*flat.fire()-1.0);
			x[2] = fZ0[0] + fSigmaZ[0]*(2.0*flat.fire()-1.0);
		}

		double t = 0.;
		if(fTDist==kGAUS){
			t = gauss.fire(fT0[0], fSigmaT[0]);
		}
		else{
			t = fT0[0] + fSigmaT[0]*(2.0*flat.fire()-1.0);
		}

		TLorentzVector pos(x[0], x[1], x[2], t);


		//STEP2 Prepare the mother particle, we only need the momentum, p0vec
		// Choose momentum
		double p = 0.0;

		if (fPDist == kGAUS) {
			p = abs(gauss.fire(fP0[0], fSigmaP[0]));// CHECK, use absolute value
		}
		else if (fPDist == kHIST){
			p = SelectFromHist(*(hPHist[0]));
		}
		else {
			p = fP0[0] + fSigmaP[0]*(2.0*flat.fire()-1.0);
		}

		// Choose mass
		double m = 0.0;
		if (fMotherMassDist == kGAUS) {//p1vec[last_index] = p0vec upon initialization;
			m = gauss.fire(fMotherMass[0], fSigmaMotherMass[0]);
		}
		else {
			m = fMotherMass[0] + fSigmaMotherMass[0]*(2.0*flat.fire()-1.0);
		}


		// Choose angles
		double thxz = 0;//0,360degrees
		double thyz = 0;//-90degrees,90degrees

		if (fAngleDist == kGAUS) {
			thxz = gauss.fire(fTheta0XZ[0], fSigmaThetaXZ[0]);
			thyz = gauss.fire(fTheta0YZ[0], fSigmaThetaYZ[0]);
		}
		else if (fAngleDist == kHIST){
			double thetaxz = 0;
			double thetayz = 0;
			SelectFromHist(*(hThetaXzYzHist[0]), thetaxz, thetayz);
			thxz = (180./M_PI)*thetaxz;
			thyz = (180./M_PI)*thetayz;
		} 
		else { // Choose angles flat in phase space, which is flat in theta_xz 
			// and flat in sin(theta_yz).

			thxz = fTheta0XZ[0] + fSigmaThetaXZ[0]*(2.0*flat.fire()-1.0);

			double thyzrads = std::asin(std::sin((M_PI/180.)*(fTheta0YZ[0]))); //Taking asin of sin gives value between -Pi/2 and Pi/2 regardless of user input
			double thyzradsplussigma = TMath::Min((thyzrads + ((M_PI/180.)*fabs(fSigmaThetaYZ[0]))), M_PI/2.);
			double thyzradsminussigma = TMath::Max((thyzrads - ((M_PI/180.)*fabs(fSigmaThetaYZ[0]))), -M_PI/2.);

			//std::cout << "Central angle: " << (180./M_PI)*thyzrads << " Max angle: " << (180./M_PI)*thyzradsplussigma << " Min angle: " << (180./M_PI)*thyzradsminussigma << std::endl; 

			double sinthyzmin = std::sin(thyzradsminussigma);
			double sinthyzmax = std::sin(thyzradsplussigma);
			double sinthyz = sinthyzmin + flat.fire() * (sinthyzmax - sinthyzmin);
			thyz = (180. / M_PI) * std::asin(sinthyz);

			if(fverbose) std::cout<<"Mother particle angles: XZ "<<thxz<<" YZ: "<<thyz<<" Momentum: "<<p<<" mass:" <<m<<std::endl;
		}

		TLorentzVector p0vec(p*std::cos(thyz*M_PI/180.0)*std::sin(thxz*M_PI/180.0),
				p*std::sin(thyz*M_PI/180.0),
				p*std::cos(thxz*M_PI/180.0)*std::cos(thyz*M_PI/180.0),
				std::sqrt(p*p+m*m));



		//STEP2, daughter particles at rest frame

		// Choose daughter 1 angles
		double Theta0XZEXT = 0;//0,360degrees
		double Theta0YZEXT = 0;//-90degrees,90degrees

		if (fAngEXTDist == kGAUS) {
			Theta0XZEXT = gauss.fire(fTheta0XZEXT[0], fSigmaTheta0XZEXT[0]);
			Theta0YZEXT = gauss.fire(fTheta0YZEXT[0], fSigmaTheta0YZEXT[0]);
		}
		else { // Choose angles flat in phase space, which is flat in theta_xz 
			// and flat in sin(theta_yz).

			Theta0XZEXT = fTheta0XZEXT[0] + fSigmaTheta0XZEXT[0]*(2.0*flat.fire()-1.0);

			double thyzrads = std::asin(std::sin((M_PI/180.)*(fTheta0YZEXT[0]))); //Taking asin of sin gives value between -Pi/2 and Pi/2 regardless of user input
			double thyzradsplussigma = TMath::Min((thyzrads + ((M_PI/180.)*fabs(fSigmaTheta0YZEXT[0]))), M_PI/2.);
			double thyzradsminussigma = TMath::Max((thyzrads - ((M_PI/180.)*fabs(fSigmaTheta0YZEXT[0]))), -M_PI/2.);

//			std::cout << "Central angle: " << (180./M_PI)*thyzrads << " Max angle: " << (180./M_PI)*thyzradsplussigma << " Min angle: " << (180./M_PI)*thyzradsminussigma << std::endl; 

			double sinthyzmin = std::sin(thyzradsminussigma);
			double sinthyzmax = std::sin(thyzradsplussigma);
			double sinthyz = sinthyzmin + flat.fire() * (sinthyzmax - sinthyzmin);
			Theta0YZEXT = (180. / M_PI) * std::asin(sinthyz);

			if(fverbose) std::cout<<"At mother's rest frame daughter 1 angles XZ "<<Theta0XZEXT<<" YZ: "<<Theta0YZEXT<<std::endl;
		}


		// Choose daughter mass (0 for photon)
		static TDatabasePDG  pdgt;
		TParticlePDG* pdg1 = pdgt.GetParticle(fPDG[0]);
		TParticlePDG* pdg2 = pdgt.GetParticle(fPDG[1]);
		double dm1 = (pdg1)? pdg1->Mass(): 0 ;//daugther 1 mass
		double dm2 = (pdg2)? pdg2->Mass(): 0 ;

		// daughter momentum is set after daughter masses are set;
		double p1 = TMath::Sqrt( 
				pow( (pow(m,2) - pow(dm1,2) + pow(dm2,2) )/(2*m) , 2) - pow(dm2,2)
				);//conserve p & E in rest frame
		// daughter 4-momentum rest frame 
		TLorentzVector p1vec(
				p1*std::cos(Theta0YZEXT*M_PI/180.0)*std::sin(Theta0XZEXT*M_PI/180.0),
				p1*std::sin(Theta0YZEXT*M_PI/180.0),
				p1*std::cos(Theta0XZEXT*M_PI/180.0)*std::cos(Theta0YZEXT*M_PI/180.0),
				std::sqrt(p1*p1+dm1*dm1));

//		std::cout<<" Daughter stat : momentum (rest) "<< p1<<" px "<<p1vec.X()<<std::endl;

		// boost daughter 1 to the lab frame beta^2=1/((m/p^2+1))
		p1vec.Boost( p0vec.X()/p0vec.E(),
			p0vec.Y()/p0vec.E(),
			p0vec.Z()/p0vec.E()
			);
//		std::cout<<" Daughter stat : mother mass "<< m <<" boost x "<<-p0vec.X()/p0vec.E()<<std::endl;
//		std::cout<<" Daughter stat : mother mass "<< m <<" boost y "<<-p0vec.Y()/p0vec.E()<<std::endl;
//		std::cout<<" Daughter stat : mother mass "<< m <<" boost z "<<-p0vec.Z()/p0vec.E()<<std::endl;
//		std::cout<<" Daughter stat : mass 1 (lab) "<< dm1<<"  m2 "<<dm2<<std::endl;

		TLorentzVector p2vec(
				-p1*std::cos(Theta0YZEXT*M_PI/180.0)*std::sin(Theta0XZEXT*M_PI/180.0),
				-p1*std::sin(Theta0YZEXT*M_PI/180.0),
				-p1*std::cos(Theta0XZEXT*M_PI/180.0)*std::cos(Theta0YZEXT*M_PI/180.0),
				std::sqrt(p1*p1+dm2*dm2));

		p2vec.Boost( p0vec.X()/p0vec.E(),
			p0vec.Y()/p0vec.E(),
			p0vec.Z()/p0vec.E()
			);


		simb::MCParticle part1(-1, fPDG[0], "primary");
		part1.AddTrajectoryPoint(pos, p1vec);//use mother particle's location
		mct.Add(part1);

		simb::MCParticle part2(-1, fPDG[1], "primary");
		part2.AddTrajectoryPoint(pos, p2vec);
		mct.Add(part2);


		//Print daughter's information
		if(fverbose){
			std::cout<<"------------- Summary of Particles Kinematics --------------"<<std::endl;
			std::cout<<std::setw(12)<<"Part ";
			std::cout<<std::setw(13)<<"px ";
			std::cout<<std::setw(13)<<"py ";
			std::cout<<std::setw(13)<<"pz ";
			std::cout<<std::setw(11)<<"E ";
			std::cout<<std::setw(11)<<"x ";
			std::cout<<std::setw(11)<<"y ";
			std::cout<<std::setw(11)<<"z ";
			std::cout<<std::setw(11)<<"T ";
			std::cout<<std::endl;

			std::cout<<std::setw(12)<<"Mother";
			std::cout<<std::setw(13)<<p0vec.Px();
			std::cout<<std::setw(13)<<p0vec.Py();
			std::cout<<std::setw(13)<<p0vec.Pz();
			std::cout<<std::setw(11)<<p0vec.E();
			std::cout<<std::setw(11)<<pos.X();
			std::cout<<std::setw(11)<<pos.Y();
			std::cout<<std::setw(11)<<pos.Z();
			std::cout<<std::setw(11)<<pos.T();
			std::cout<<std::endl;
		
			std::cout<<std::setw(12)<<"Daug. 1";
			std::cout<<std::setw(13)<<p1vec.Px();
			std::cout<<std::setw(13)<<p1vec.Py();
			std::cout<<std::setw(13)<<p1vec.Pz();
			std::cout<<std::setw(11)<<p1vec.E();
			std::cout<<std::endl;

			std::cout<<std::setw(12)<<"Daug. 2";
			std::cout<<std::setw(13)<<p2vec.Px();
			std::cout<<std::setw(13)<<p2vec.Py();
			std::cout<<std::setw(13)<<p2vec.Pz();
			std::cout<<std::setw(11)<<p2vec.E();
			std::cout<<"\n"<<std::endl;
		}
	}


  //____________________________________________________________________________
  void TwoBodyDecayGen::Sample(simb::MCTruth &mct) 
  {

    switch (fMode) {
    case 0: // List generation mode: every event will have one of each
	    // particle species in the fPDG array
        if (fSingleVertex){
          SampleMany(mct);
        }
        else{
          for (unsigned int i=0; i<fPDG.size(); ++i) {
            SampleOne(i,mct);
          }//end loop over particles
        }
        break;
    case 1: // Random selection mode: every event will exactly one particle
            // selected randomly from the fPDG array
      {
        CLHEP::RandFlat flat(*fEngine);

	unsigned int i=flat.fireInt(fPDG.size());
	SampleOne(i,mct);
      }
      break;
    default:
      mf::LogWarning("UnrecognizeOption") << "TwoBodyDecayGen does not recognize ParticleSelectionMode "
					  << fMode;
      break;
    } // switch on fMode

    return;
  }

  //____________________________________________________________________________
  void TwoBodyDecayGen::printVecs(std::vector<std::string> const& list)
  {
 
    mf::LogInfo("TwoBodyDecayGen") << " You are using vector values for TwoBodyDecayGen configuration.\n   " 
			     << " Some of the configuration vectors may have been padded out ,"
			     << " because they (weren't) as long as the pdg vector"
			     << " in your configuration. \n"
			     << " The new input particle configuration is:\n" ;

    std::string values;
    for(size_t i = 0; i <=1; ++i){// list.size(); ++i){

      values.append(list[i]);
      values.append(": [ ");      
      
      for(size_t e = 0; e < fPDG.size(); ++e){
        std::stringstream buf;
        buf.width(10);
	if(i == 0 ) buf << fPDG[e]          << ", ";
	buf.precision(5);
	if(i == 1 ) buf << fP0[e]           << ", ";
	if(i == 2 ) buf << fSigmaP[e] 	    << ", ";
	if(i == 3 ) buf << fX0[e]     	    << ", ";
	if(i == 4 ) buf << fY0[e]     	    << ", ";
	if(i == 5 ) buf << fZ0[e]	    << ", ";
	if(i == 6 ) buf << fSigmaX[e] 	    << ", ";
	if(i == 7 ) buf << fSigmaY[e] 	    << ", ";
	if(i == 8 ) buf << fSigmaZ[e] 	    << ", ";
	if(i == 9 ) buf << fTheta0XZ[e]     << ", ";
	if(i == 10) buf << fTheta0YZ[e]     << ", ";
	if(i == 11) buf << fSigmaThetaXZ[e] << ", ";
	if(i == 12) buf << fSigmaThetaYZ[e] << ", ";
	if(i == 13) buf << fT0[e]           << ", ";
	if(i == 14) buf << fSigmaT[e]       << ", ";
        values.append(buf.str());
      }

      values.erase(values.find_last_of(","));
      values.append(" ] \n");

    }// end loop over vector names in list

    mf::LogInfo("TwoBodyDecayGen") << values;

    return;
  }
  
  
  //____________________________________________________________________________
  double TwoBodyDecayGen::SelectFromHist(const TH1& h) // select from a 1D histogram
  {
    CLHEP::RandFlat   flat(*fEngine);
    
    double throw_value = h.Integral() * flat.fire();
    double cum_value(0);
    for (int i(0); i < h.GetNbinsX()+1; ++i){
      cum_value += h.GetBinContent(i);
      if (throw_value < cum_value){
        return flat.fire()*h.GetBinWidth(i) + h.GetBinLowEdge(i);
      }
    }
    return throw_value; // for some reason we've gone through all bins and failed?
  }
  //____________________________________________________________________________
  void TwoBodyDecayGen::SelectFromHist(const TH2& h, double &x, double &y) // select from a 2D histogram
  {
    CLHEP::RandFlat   flat(*fEngine);
    
    double throw_value = h.Integral() * flat.fire();
    double cum_value(0);
    for (int i(0); i < h.GetNbinsX()+1; ++i){
      for (int j(0); j < h.GetNbinsY()+1; ++j){
        cum_value += h.GetBinContent(i, j);
        if (throw_value < cum_value){
          x = flat.fire()*h.GetXaxis()->GetBinWidth(i) + h.GetXaxis()->GetBinLowEdge(i);
          y = flat.fire()*h.GetYaxis()->GetBinWidth(j) + h.GetYaxis()->GetBinLowEdge(j);
          return;
        }
      }
    }
    return; // for some reason we've gone through all bins and failed?
  }
  //____________________________________________________________________________


}//end namespace evgen

namespace evgen{

  DEFINE_ART_MODULE(TwoBodyDecayGen)

}//end namespace evgen

#endif
////////////////////////////////////////////////////////////////////////
