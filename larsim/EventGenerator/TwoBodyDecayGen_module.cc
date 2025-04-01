////////////////////////////////////////////////////////////////////////
/// \file  TwoBodyDecayGen_module.cc
/// \brief Generator for hypothetical particles that decay into two particles
/// 
/// 
/// Update in Mar. 2025
/// ----
/// Use a set of 2D histogram to generate particle kinematics.
/// Add timing information into the simulation. 
///
/// First version
/// ----
/// Module designed to produce two random particles decay from one mother particle.
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
#include "nusimdata/SimulationBase/MCFlux.h"
#include "nusimdata/SimulationBase/MCParticle.h"
#include "nutools/EventGeneratorBase/evgenbase.h"

// lar includes
//#include "larcore/Geometry/Geometry.h"
//#include "larcoreobj/SummaryData/RunData.h"

#include "lardata/Utilities/AssociationUtil.h"

#include "TVector3.h"
#include "TDatabasePDG.h"
#include "TMath.h"
#include "TFile.h"
#include "TH1.h"
#include "TH2.h"

#include "CLHEP/Random/RandFlat.h"
#include "CLHEP/Random/RandGaussQ.h"


//-----------------
//---- Overview----
//beginRun() kicks off the generation in a run
//produce(), then Sample().
//art::EDProducer define the class, declare variables
//setup() loads the FHiCL files
//SampleMany() is the actual calculation for the twobodydecay.
//produce() store particles
//-----------------

namespace simb { class MCTruth; }

namespace evgen {

  /// module to produce single or multiple specified particles in the detector
  class TwoBodyDecayGen : public art::EDProducer {

    public:

        //--- Validation of configurations
      struct Config {
        using Name = fhicl::Name;
        using Comment = fhicl::Comment;
//Set up FHiCL Parameters fhicl::Atom<type> X( <NAME>, <Note>, <Value>)
//See https://mu2ewiki.fnal.gov/wiki/FclIntro
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
            [this](){ return fromHistogram(AngleDist()) || fromHistogram(PDist()); }
        };

//      fhicl::Sequence<std::string> PHist{
//        Name("PHist"),
//        Comment("name of the histograms of momentum distributions"),
//        [this](){ return fromHistogram(PDist()); }
//      };

//        fhicl::Sequence<std::string> ThetaXzYzHist{
//          Name("ThetaXzYzHist"),
//            Comment("name of the histograms of angular (X-Z and Y-Z) distribution"),
//            [this](){ return fromHistogram(AngleDist()); }
//        };

        fhicl::OptionalAtom<rndm::NuRandomService::seed_t> Seed{
          Name("Seed"),
            Comment("override the random number generator seed")
        };


        //-----
        //Begin Exclusive Parameters
        //-----
        //Mother Particle Mass
        fhicl::Atom<std::string> MotherMassDist{
          Name("MotherMassDist"),
            Comment("mother mass distribution type: " + presentOptions(DistributionNames)),
            optionName(kHIST, DistributionNames)
        };

        fhicl::Sequence<double> MotherMass{Name("MotherMass"),
          Comment("mother mass (GeV/c^2) to generate"),
          //          [this]() { return !fromHistogram(PDist()); }
        };

        fhicl::Sequence<double> SigmaMotherMass{Name("SigmaMotherMass"),
          Comment("variation in mother mass (GeV/c^2)"),
          //          [this]() { return !fromHistogram(PDist()); }
        };

         //Mother Particle Momentum
        fhicl::Sequence<std::string> ZPHist{
          Name("ZPHist"),
            Comment("name of the histograms of creation Z [cm] - momentum [GeV] distributions"),
            [this](){ return fromHistogram(PDist()); }
        };

        fhicl::Sequence<std::string> ZTHist{
          Name("ZTHist"),
            Comment("name of the histograms of creation Z [cm] - time [ns] distributions"),
            [this](){ return fromHistogram(PDist()); }
        };

        fhicl::Sequence<std::string> ZTheta0XZHist{
          Name("ZTheta0XZHist"),
            Comment("name of the histograms of creation Z [cm] - XZ Angle [degree] distributions"),
            [this](){ return fromHistogram(AngleDist()); }
        };


        fhicl::Sequence<std::string> ZTheta0YZHist{
          Name("ZTheta0YZHist"),
            Comment("name of the histograms of creation Z [cm] - YZ Angle [degree] distributions"),
            [this](){ return fromHistogram(AngleDist()); }
        };

        //Mother particle travel distance before entering the detector, dist1; dist2 is the exiting distance, ignore it for now;
        fhicl::Atom<std::string> Traveldist1Dist{
          Name("Traveldist1Dist"),
            Comment("Travel distance estimator type: " + presentOptions(DistributionNames)),
            optionName(kHIST, DistributionNames)
        };

        fhicl::Sequence<std::string> Zdist1Hist{
          Name("Zdist1Hist"),
            Comment("name of the histograms of creation Z [cm] - Travel distance [cm] distributions"),
            [this](){ return fromHistogram(AngleDist()); }
        };


        //First Daughter Outgoing Angle in the rest frame
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



        //-----
        //End Exclusive Parameters
        //-----


        private:

        /// Returns whether the specified mode is an histogram distribution.
        bool fromHistogram(std::string const& key) const;

      }; // struct Config


      using Parameters = art::EDProducer::Table<Config>;//Using the description service


      explicit TwoBodyDecayGen(Parameters const& config);//apply the validation

      // This is called for each event.
      void produce(art::Event& evt);
      void beginRun(art::Run& run);

      //Begin exclusive function
      //Extra piece information to store where the particle is created.
      TLorentzVector GetCeationPoint(std::vector<double> boxDim, TLorentzVector dkvertex, TLorentzVector dir, double dist);
      //End exclusive function

    private:

      /// Names of all particle selection modes.
      static const std::map<int, std::string> ParticleSelectionModeNames;
      /// Names of all distribution modes.
      static const std::map<int, std::string> DistributionNames;

      void SampleOne(unsigned int   i, 
          simb::MCTruth &mct);        
      void SampleMany(simb::MCTruth &mct, simb::MCFlux &flux);
      void Sample(simb::MCTruth &mct, simb::MCFlux &flux);
      void printVecs(std::vector<std::string> const& list);
      bool PadVector(std::vector<double> &vec);      
      double SelectFromHist(const TH1& h);
      void SelectFromHist(const TH2& h, double &x, double &y);
      void SelectFromHistFixX(const TH2& h, double &x, double &y);

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

      //Begin Exclusive Parameters
      std::vector<double> fMotherMass;       ///Mother particle can be massive
      std::vector<double> fSigmaMotherMass;  ///The variation on the mass
      int                 fMotherMassDist;   ///How to distribute the mass
      std::vector<double> fTheta0XZEXT;      ///
      std::vector<double> fSigmaTheta0XZEXT;
      int                 fAngEXTDist;
      std::vector<double> fTheta0YZEXT;
      std::vector<double> fSigmaTheta0YZEXT;
      
      //Histogram Objects with prefex of h
      std::vector<std::unique_ptr<TH2>> hZPHist ;     /// actual TH2 for Creation Z (CZ) - momentum(GeV) distributions
      std::vector<std::unique_ptr<TH2>> hZTHist ;     /// actual TH2 for Creation Z (CZ) - time(ns) distributions
      std::vector<std::unique_ptr<TH2>> hZTheta0XZHist ; /// actual TH2 for angle distributions - Xz on x axis with fixed CZ
      std::vector<std::unique_ptr<TH2>> hZTheta0YZHist ; /// actual TH2 for angle distributions - Yz on x axis with fixed CZ
      std::vector<std::unique_ptr<TH2>> hZdist1Hist ;     /// actual TH2 for Creation Z (CZ) - travel distance distributions

      std::vector<std::string> fZPHist;     ///< name of histogram of momenta
      std::vector<std::string> fZTHist;     ///< name of histogram of time
      std::vector<std::string> fZTheta0XZHist;     ///< name of histogram of xz angle
      std::vector<std::string> fZTheta0YZHist;     ///< name of histogram of yz angle
      std::vector<std::string> fZdist1Hist;     ///< name of histogram of travel distance

      int                        fTraveldist1Dist; //distribution types for travel distance
      //End Exclusive Parameters

      std::string fHistFileName;               ///< Filename containing histogram of momenta
      std::vector<std::string> fPHist;     ///< name of histogram of momenta
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

    //Begin exclusive function
    TLorentzVector TwoBodyDecayGen::GetCeationPoint(std::vector<double> boxDim, TLorentzVector dkvertex, TLorentzVector dir, double dist){
            
        TVector3 tmp_pos(dkvertex.X(), dkvertex.Y(), dkvertex.Z());//decay vertex
        TVector3 tmp_dir(dir.X(), dir.Y(), dir.Z());//3-momentum
        tmp_dir = tmp_dir.Unit();

        double t_min = -1, t_max = 1e20;

        for( int index = 0; index < 3; index++){//0,1,2 = x,y,z
            double pmin = boxDim[2*index];//get the corner of the boxDim
            double pmax = boxDim[2*index+1];

            if (dir[index] != 0) {  // Avoid division by zero
                double t1 = (pmin - dkvertex[index]) / dir[index];
                double t2 = (pmax - dkvertex[index]) / dir[index];

                if (t1 > t2) std::swap(t1, t2);
                t_min = std::max(t_min, t1);//get the closest edge in,x, y, or z opposite to the dir
                t_max = std::min(t_max, t2);//get the closest edge in,x, y, or z along the dir
            } else if (dkvertex[index] < pmin || dkvertex[index] > pmax) {
                std::cout<<"Error, particles are simulated outside the defined region. Return the decay vertex."<<std::endl;
                return dkvertex;
            }
        }//t_min is used to projected the intersecting point for particle entering.

        if (t_min > t_max || t_max < 0) {
            std::cout<<"Error, intersection calculation does not make sense. Return the decay vertex."<<std::endl;
            return dkvertex;  // No valid intersection
        }

        TVector3 entry = tmp_pos - t_min * tmp_dir;  // Compute entry point

        // Find point at a given dist along the path
        double t_length = dist;//because tmp_dir is a unit vector;
        TVector3 creation_pos = entry - t_length * tmp_dir;

        TLorentzVector mpcp(creation_pos, dkvertex.T());//mother particle creation point, to be calculated

        return mpcp;

    }
    //end exclusive function

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

  //_______Note that the order needs to follow the order in the private class___________________________________
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
    //Begin Exclusive Parameters
    , fMotherMass(config().MotherMass())
    , fSigmaMotherMass(config().SigmaMotherMass())
    , fMotherMassDist(selectOption(config().MotherMassDist(), DistributionNames))
    , fTheta0XZEXT(config().Theta0XZEXT())
    , fSigmaTheta0XZEXT(config().SigmaTheta0XZEXT())
    , fAngEXTDist(selectOption(config().AngEXTDist(), DistributionNames))
    , fTheta0YZEXT(config().Theta0YZEXT())
    , fSigmaTheta0YZEXT(config().SigmaTheta0YZEXT())
    , fZPHist        (config().ZPHist())
    , fZTHist        (config().ZTHist())
    , fZTheta0XZHist  (config().ZTheta0XZHist())
    , fZTheta0YZHist  (config().ZTheta0YZHist())
    , fZdist1Hist     (config().Zdist1Hist())
    , fTraveldist1Dist(selectOption(config().Traveldist1Dist(), DistributionNames))
    //End Exclusive Parameters
    , fHistFileName (config().HistogramFile())
//    , fPHist        (config().PHist())
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

    //Take care of parameters for numerical values
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

    // Begin Exclusive Parameters
    vlist[15] = "MotherMass";             
    vlist[16] = "SigmaMotherMass";        
    vlist[17] = "Theta0XZEXT";           
    vlist[18] = "SigmaTheta0XZEXT";      
    vlist[19] = "Theta0YZEXT";     
    vlist[20] = "SigmaTheta0YZEXT";
    // End Exclusive Parameters
    //
    // Beacuse I want MCFlux
    produces< std::vector<simb::MCFlux>  >();
    produces< art::Assns<simb::MCTruth, simb::MCFlux> >();

    //    vlist[15] = "ZPHist";
    //    vlist[16] = "ThetaHist";
    //    vlist[17] = "PhiHist";

    // begin tests for multiple particle error possibilities  
    std::string list;
    if (fPDist != kHIST) {
      if( !this->PadVector(fP0            ) ){ list.append(vlist[1].append(", \n")); }
      if( !this->PadVector(fSigmaP        ) ){ list.append(vlist[2].append(", \n")); }
    }

    if(fAngleDist != kHIST){ 
        if( !this->PadVector(fTheta0XZ        ) ){ list.append(vlist[9].append(", \n")); }
        if( !this->PadVector(fTheta0YZ        ) ){ list.append(vlist[10].append(", \n")); }
        if( !this->PadVector(fSigmaThetaXZ    ) ){ list.append(vlist[11].append(", \n")); }
        if( !this->PadVector(fSigmaThetaYZ    ) ){ list.append(vlist[12].append("  \n")); }
    }

    if( !this->PadVector(fX0              ) ){ list.append(vlist[3].append(", \n")); }
    if( !this->PadVector(fY0              ) ){ list.append(vlist[4].append(", \n")); }
    if( !this->PadVector(fZ0              ) ){ list.append(vlist[5].append(", \n")); }
    if( !this->PadVector(fSigmaX          ) ){ list.append(vlist[6].append(", \n")); }
    if( !this->PadVector(fSigmaY          ) ){ list.append(vlist[7].append(", \n")); }
    if( !this->PadVector(fSigmaZ          ) ){ list.append(vlist[8].append(", \n")); }
    if( !this->PadVector(fT0              ) ){ list.append(vlist[13].append(", \n")); }
    if( !this->PadVector(fSigmaT          ) ){ list.append(vlist[14].append(", \n")); }

    // Begin Exclusive Parameters
    if( !this->PadVector(fMotherMass              ) ){ list.append(vlist[15].append(", \n")); }
    if( !this->PadVector(fSigmaMotherMass         ) ){ list.append(vlist[16].append(", \n")); }
    if( !this->PadVector(fTheta0XZEXT            ) ){ list.append(vlist[17].append(", \n")); }
    if( !this->PadVector(fSigmaTheta0XZEXT       ) ){ list.append(vlist[18].append(", \n")); }
    if( !this->PadVector(fTheta0YZEXT      ) ){ list.append(vlist[19].append(", \n")); }
    if( !this->PadVector(fSigmaTheta0YZEXT ) ){ list.append(vlist[20].append(", \n")); }
    // End Exclusive Parameters


    if(list.size() > 0)
      throw cet::exception("TwoBodyDecayGen") << "The "<< list 
        << "\n vector(s) defined in the fhicl files has/have "
        << "a different size than the PDG vector "
        << "\n and it has (they have) more than one value, "
        << "\n disallowing sensible padding "
        << " and/or you have set fPadOutVectors to false. \n";

    if(fPDG.size() > 1 && fPadOutVectors) this->printVecs(vlist);

    //Prepare Histogram File if it needed
    // If needed, get histograms for momentum, angle, and Time distributions
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

    // Files are ready, then Fill in histograms 
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
        case kHIST:
            hZTHist.reserve(fZTHist.size());
            for (auto const& histName: fZTHist) {
                TH2* temHist = dynamic_cast<TH2*>(histFile->Get(histName.c_str()));
                if (!temHist) {
                    throw art::Exception(art::errors::NotFound)
                        << "Failed to read direction histogram '" << histName << "' from '" << histFile->GetPath() << "\'";
                }
                temHist->SetDirectory(nullptr); // make it independent of the input file
                hZTHist.emplace_back(temHist);
            } // for
            break;

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
       //Keng disables the reading of 1d spectrum of PHist
       // hPHist.reserve(fPHist.size());
       // for (auto const& histName: fPHist) {
       //   TH1* temHist = dynamic_cast<TH1*>(histFile->Get(histName.c_str()));
       //   if (!temHist) {
       //     throw art::Exception(art::errors::NotFound)
       //       << "Failed to read momentum histogram '" << histName << "' from '" << histFile->GetPath() << "\'";
       //   }
       //   temHist->SetDirectory(nullptr); // make it independent of the input file
       //   hPHist.emplace_back(temHist);
       // } // for

        hZPHist.reserve(fZPHist.size());
        for (auto const& histName: fZPHist) {
          TH2* temHist = dynamic_cast<TH2*>(histFile->Get(histName.c_str()));
          if (!temHist) {
            throw art::Exception(art::errors::NotFound)
              << "Failed to read XZ direction histogram '" << histName << "' from '" << histFile->GetPath() << "\'";
          }
          temHist->SetDirectory(nullptr); // make it independent of the input file
          hZPHist.emplace_back(temHist);
        } // for ZP histograms


        break;
      default: // supported, no further action needed
        break;
    } // switch(fPDist)

    switch (fAngleDist) {
      case kHIST:
        // if (fThetaXzYzHist.size() != fPDG.size()) {
        //   throw art::Exception(art::errors::Configuration)
        //     << fThetaXzYzHist.size() << " direction histograms to describe " << fPDG.size() << " particle types...";
        // }
        //Keng disables the drawing of both angles in 1 2dhist.
        //hThetaXzYzHist.reserve(fThetaXzYzHist.size());
        //for (auto const& histName: fThetaXzYzHist) {
        //  TH2* temHist = dynamic_cast<TH2*>(histFile->Get(histName.c_str()));
        //  if (!temHist) {
        //    throw art::Exception(art::errors::NotFound)
        //      << "Failed to read direction histogram '" << histName << "' from '" << histFile->GetPath() << "\'";
        //  }
        //  temHist->SetDirectory(nullptr); // make it independent of the input file
        //  hThetaXzYzHist.emplace_back(temHist);
        //} // Two more histograms for angles
        //First Xz
        hZTheta0XZHist.reserve(fZTheta0XZHist.size());
        for (auto const& histName: fZTheta0XZHist) {
          TH2* temHist = dynamic_cast<TH2*>(histFile->Get(histName.c_str()));
          if (!temHist) {
            throw art::Exception(art::errors::NotFound)
              << "Failed to read XZ direction histogram '" << histName << "' from '" << histFile->GetPath() << "\'";
          }
          temHist->SetDirectory(nullptr); // make it independent of the input file
          hZTheta0XZHist.emplace_back(temHist);
        } // Two more histograms for angles

        //Then Yz
        hZTheta0YZHist.reserve(fZTheta0YZHist.size());
        for (auto const& histName: fZTheta0YZHist) {
          TH2* temHist = dynamic_cast<TH2*>(histFile->Get(histName.c_str()));
          if (!temHist) {
            throw art::Exception(art::errors::NotFound)
              << "Failed to read direction histogram '" << histName << "' from '" << histFile->GetPath() << "\'";
          }
          temHist->SetDirectory(nullptr); // make it independent of the input file
          hZTheta0YZHist.emplace_back(temHist);
        } // Two more histograms for angles


        break;
      default: // supported, no further action needed
        break;
    } // switch(fAngleDist)

    switch(fTraveldist1Dist){//Mother particle travel distance1
        case(kHIST):
            hZdist1Hist.reserve(fZdist1Hist.size());
            for (auto const& histName: fZdist1Hist) {
                TH2* temHist = dynamic_cast<TH2*>(histFile->Get(histName.c_str()));
                if (!temHist) {
                    throw art::Exception(art::errors::NotFound)
                        << "Failed to read direction histogram '" << histName << "' from '" << histFile->GetPath() << "\'";
                }
                temHist->SetDirectory(nullptr); // make it independent of the input file
                hZdist1Hist.emplace_back(temHist);
            } // Two more histograms for angles

            break;
        default:
            break;
    }

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
      std::unique_ptr< std::vector<simb::MCFlux>  > fluxcol   (new std::vector<simb::MCFlux >);
      std::unique_ptr< art::Assns<simb::MCTruth, simb::MCFlux> > tfassn(new art::Assns<simb::MCTruth, simb::MCFlux>);

      simb::MCTruth truth;
      simb::MCFlux flux;
      truth.SetOrigin(simb::kSingleParticle);
      Sample(truth, flux);

      MF_LOG_DEBUG("TwoBodyDecayGen") << truth;

      truthcol->push_back(truth);
      fluxcol->push_back(flux);

      util::CreateAssn(*this, evt, *truthcol, *fluxcol, *tfassn, fluxcol->size()-1, fluxcol->size());

      evt.put(std::move(truthcol));
      evt.put(std::move(fluxcol));
	  evt.put(std::move(tfassn));

      return;
  }

  //____________________________________________________________________________
  // Draw the type, momentum and position of a single particle from the 
  // FCIHL description
  void TwoBodyDecayGen::SampleOne(unsigned int i, simb::MCTruth &mct){
  }

  //____________________________________________________________________________
  // Draw the type, momentum and position for all particles from the 
  // FCIHL description.  
  // Use TH2D for P, Time, XZ, and YZ angles.
  void TwoBodyDecayGen::SampleMany(simb::MCTruth &mct, simb::MCFlux &flux){

    bool fverbose = true;

    CLHEP::RandFlat   flat(*fEngine);
    CLHEP::RandGaussQ gauss(*fEngine);

    //STEP1 Choose the production vertex, all particles come from the same position
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


    //STEP2 Prepare the mother particle kinematics, 
    // define the momentum & Creation Z --> Arrival time, XZ, and YZ angle
    // Choose momentum
    double p = 0.0;

    double CZ = 0.0;//This will be the key variable connecting all relative variables. 
    //This is special, and it is measured in the NuMI coordinate. Do not use CZ to measure distance without conversion

    if (fPDist == kGAUS) {
      p = abs(gauss.fire(fP0[0], fSigmaP[0]));
    }
    else if (fPDist == kHIST){

        // p = SelectFromHist(*(hPHist[0]));
        double tmpCZ = 0;
        double tmpp = 0;
        SelectFromHist(*(hZPHist[0]), tmpCZ, tmpp);
        CZ = tmpCZ;
        p = tmpp;

    }
    else {
      std::cout<<"Error, momentum has to be drawn from the histogram"<<std::endl;
      p = fP0[0] + fSigmaP[0]*(2.0*flat.fire()-1.0);
    }

    //Choose Time
    double t = 0.;
    if(fTDist==kGAUS){
      t = gauss.fire(fT0[0], fSigmaT[0]);
    }
    else if (fTDist == kHIST){

        double tmpt;
        SelectFromHistFixX(*(hZTHist[0]), CZ, tmpt);
        t = tmpt;
    }
    else{
      t = fT0[0] + fSigmaT[0]*(2.0*flat.fire()-1.0);
    }

    TLorentzVector pos(x[0], x[1], x[2], t);


    // Determine mother particle mass
    double m = 0.0;
    if (fMotherMassDist == kGAUS) {//p1vec[last_index] = p0vec upon initialization;
      m = gauss.fire(fMotherMass[0], fSigmaMotherMass[0]);
    }
    else {
      m = fMotherMass[0] + fSigmaMotherMass[0]*(2.0*flat.fire()-1.0);
    }

    //Need to correct p, which was taken as energy
    if(p > m ){//Physical, p was actually energy
        p = sqrt(p*p - m*m);  
    } else {
        return;
    }


    // Choose Mother angles
    double thxz = 0;//0,360degrees
    double thyz = 0;//-90degrees,90degrees

    if (fAngleDist == kGAUS) {
      thxz = gauss.fire(fTheta0XZ[0], fSigmaThetaXZ[0]);
      thyz = gauss.fire(fTheta0YZ[0], fSigmaThetaYZ[0]);
    }
    else if (fAngleDist == kHIST){
      double thetaxz = 0;
      double thetayz = 0;
      SelectFromHistFixX(*(hZTheta0XZHist[0]), CZ, thetaxz);
      SelectFromHistFixX(*(hZTheta0YZHist[0]), CZ, thetayz);
      thxz = thetaxz;
      thyz = thetayz;
      //thxz = (180./M_PI)*thetaxz;
      //thyz = (180./M_PI)*thetayz;
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

    }
    if(fverbose) std::cout<<"Mother particle angles: XZ "<<thxz<<" YZ: "<<thyz<<" Momentum: "<<p<<" mass:" <<m<<std::endl;

    //Mother 4-Momentum
    TLorentzVector p0vec(p*std::cos(thyz*M_PI/180.0)*std::sin(thxz*M_PI/180.0),
        p*std::sin(thyz*M_PI/180.0),
        p*std::cos(thxz*M_PI/180.0)*std::cos(thyz*M_PI/180.0),
        std::sqrt(p*p+m*m));

    //We also want to know the mother travel distance before entering the detector
    double dist1 = 0;
    if(fTraveldist1Dist == kHIST){
        double tmpdist = 0;
        SelectFromHistFixX( *(hZdist1Hist[0]), CZ, tmpdist);
        dist1= tmpdist;
    }


    //CHECK, we want to save t and dist1 as an object of MCFlux
    //Option 1: translate dist1 usnig pos, boundary & p0vec.
    std::vector<double> boxDim = {fX0[0] - fSigmaX[0], fX0[0] + fSigmaX[0],
        fY0[0] - fSigmaY[0], fY0[0] + fSigmaY[0],
        fZ0[0] - fSigmaZ[0], fZ0[0] + fSigmaZ[0]};

    
    TLorentzVector CreationP = GetCeationPoint(boxDim, pos, p0vec, dist1);
    std::cout<<"Projected creation point at "<<CreationP.X()<<","<<CreationP.Y()<<","<<CreationP.Z()<<" at T="<<CreationP.T()<<std::endl;
    //Now CreationP is the MCFlux particle
	//What these variables hould be? See https://internal.dunescience.org/doxygen/MCFlux_8h_source.html
    flux.Reset();
    flux.fFluxType = simb::kNtuple;
    
    flux.fvx = pos.X();
    flux.fvy = pos.Y();
    flux.fvz = pos.Z();

    flux.fgenx = CreationP.X();//mother particle creation point
    flux.fgeny = CreationP.Y();
    flux.fgenz = CreationP.Z();
	flux.fdk2gen = 0;// distance from decay to ray origin .. 0?
	flux.fgen2vtx = (CreationP.Vect() - pos.Vect()).Mag(); //distance from ray origin to event vtx.
    //Finish adding MCFlux;


    //STEP2, daughter particles at rest frame

    // Choose daughter 1 angles
    double Theta0XZEXT = 0;//0,360 degrees
    double Theta0YZEXT = 0;//-90 degrees,90degrees

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

      //std::cout << "Central angle: " << (180./M_PI)*thyzrads << " Max angle: " << (180./M_PI)*thyzradsplussigma << " Min angle: " << (180./M_PI)*thyzradsminussigma << std::endl; 

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

    //        std::cout<<" Daughter stat : momentum (rest) "<< p1<<" px "<<p1vec.X()<<std::endl;

    // boost daughter 1 to the lab frame beta^2=1/((m/p^2+1))
    p1vec.Boost( p0vec.X()/p0vec.E(),
        p0vec.Y()/p0vec.E(),
        p0vec.Z()/p0vec.E()
        );
    //        std::cout<<" Daughter stat : mother mass "<< m <<" boost x "<<-p0vec.X()/p0vec.E()<<std::endl;
    //        std::cout<<" Daughter stat : mother mass "<< m <<" boost y "<<-p0vec.Y()/p0vec.E()<<std::endl;
    //        std::cout<<" Daughter stat : mother mass "<< m <<" boost z "<<-p0vec.Z()/p0vec.E()<<std::endl;
    //        std::cout<<" Daughter stat : mass 1 (lab) "<< dm1<<"  m2 "<<dm2<<std::endl;

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
      std::cout<<std::setw(11)<<"decay_x ";
      std::cout<<std::setw(11)<<"dk_y ";
      std::cout<<std::setw(11)<<"dk_z ";
      std::cout<<std::setw(11)<<"T ";
      std::cout<<std::setw(11)<<"CZ ";
      std::cout<<std::setw(11)<<"Dist1 ";
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
      std::cout<<std::setw(11)<<CZ;
      std::cout<<std::setw(11)<<dist1;
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
  void TwoBodyDecayGen::Sample(simb::MCTruth &mct, simb::MCFlux &flux) 
  {

    switch (fMode) {
      case 0: // List generation mode: every event will have one of each
        // particle species in the fPDG array
        if (fSingleVertex){
          SampleMany(mct, flux);
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
  //Want to see what have been generated
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
        if(i == 2 ) buf << fSigmaP[e]         << ", ";
        if(i == 3 ) buf << fX0[e]             << ", ";
        if(i == 4 ) buf << fY0[e]             << ", ";
        if(i == 5 ) buf << fZ0[e]        << ", ";
        if(i == 6 ) buf << fSigmaX[e]         << ", ";
        if(i == 7 ) buf << fSigmaY[e]         << ", ";
        if(i == 8 ) buf << fSigmaZ[e]         << ", ";
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
    //For all histogram types: nbins, xlow, xup
    //bin = 0;       underflow bin
    //bin = 1;       first bin with low-edge xlow INCLUDED
    //bin = nbins;   last bin with upper-edge xup EXCLUDED
    //bin = nbins+1; overflow bin


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
    // select from a 2D histogram with a fixed X
    //____________________________________________________________________________
    void TwoBodyDecayGen::SelectFromHistFixX(const TH2& h, double &x, double &y) 
    {
      CLHEP::RandFlat   flat(*fEngine);

      double cum_value(0);
      for (int i(0); i < h.GetNbinsX()+1; ++i){
          if( !(h.GetXaxis()->GetBinLowEdge(i) < x &&  h.GetXaxis()->GetBinLowEdge(i+1) > x)) continue; 
//          std::cout<<"Search X value "<<x<<" btw bins "<<h.GetXaxis()->GetBinLowEdge(i) <<" and "<<h.GetXaxis()->GetBinLowEdge(i+1)<<std::endl;

          double throw_value = h.Integral(i,i,1, h.GetNbinsY()) * flat.fire();
          for (int j(0); j < h.GetNbinsY()+1; ++j){
              cum_value += h.GetBinContent(i, j);
              if (throw_value < cum_value){
                  y = flat.fire()*h.GetYaxis()->GetBinWidth(j) + h.GetYaxis()->GetBinLowEdge(j);
//                  std::cout<<"Got value: "<<y<<" at ("<<i<<","<<j<<") bin\n"<<std::endl;
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
