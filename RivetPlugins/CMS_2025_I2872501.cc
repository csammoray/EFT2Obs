// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/VisibleFinalState.hh"
#include "Rivet/Tools/RivetYODA.hh"

constexpr double ZMASS = 91.1876;
constexpr double MIN_MZ1 = 40.0;
constexpr double MAX_MZ1 = 120.0;
constexpr double MIN_MZ2 = 12.0;
constexpr double MAX_MZ2 = 120.0;

struct ZMassResult {
    bool passFidSel;
    std::vector<int> z_leps_idx;
};

namespace Rivet {

  int getTrueMotherPID(const Particle& part, int originalPID = 0) {
    if (originalPID == 0) originalPID = part.pid();

    for (const Particle& parent : part.parents()) {
      int pid = parent.pid();
      if (pid == 22) continue;
      if (pid != originalPID and pid != 22) return pid;

      int result = getTrueMotherPID(parent, originalPID);
      if (result != 0) return result;
    }

    return 0;
  }

  /// @brief H->ZZ analysis at 13.6 TeV with 2022 dataset
  class CMS_2025_I2872501 : public Analysis {
  public:
    /// Constructor
    // DEFAULT_RIVET_ANALYSIS_CTOR(CMS_2025_I2872501);
    RIVET_DEFAULT_ANALYSIS_CTOR(CMS_2025_I2872501);
    void init() {
      sumW_ = 0.;
      _nEventsTotal = 0;
      _nEventsFinal = 0;
      //---All final state particles
      FinalState fs;
      declare(fs, "FS");

      //---Visible final state (no neutrinos)
      VisibleFinalState vfs(fs);
      declare(vfs, "VFS");

      //---Photons
      FinalState fs_photons(Cuts::abspid == PID::PHOTON);
      declare(fs_photons, "FS_PHOTONS");

      //---Leptons
      FinalState fs_leptons(Cuts::abspid == PID::ELECTRON || Cuts::abspid == PID::MUON || Cuts::abspid == PID::TAU);
      declare(fs_leptons, "FS_LEPTONS");

      FinalState fs_electrons(Cuts::abspid == PID::ELECTRON);
      declare(fs_electrons, "FS_ELECTRONS");

      FinalState fs_muons(Cuts::abspid == PID::MUON);
      declare(fs_muons, "FS_MUONS");

      book(_histo, "hist", 1, 0, 100000);
      book(_h_ZZ_pth, "pt_h", {0,10,20,30,45,60,80,120,200,10000});
      book(_h_ZZ_pth_incl, "pt_h_incl", 1, 0, 10000);
      book(_h_ZZ_mz2, "m_z2", {12,20,24,28,32,40,55,65,120});
      // TODO: Define these two
      // book(_h_ZZ_deta, "deta_jj", {0.0,1.6,3.0,1000});
      // book(_h_ZZ_deltaphijj, "deltaphijj", {-M_PI, -M_PI/2, 0, M_PI/2, M_PI});
    }

    void analyze(const Event& event) {
      sumW_ += event.weights()[0];
      _histo->fill(sumW_);

      ++_nEventsTotal;

      Particles fsr_photons;

      Particles vfs = apply<FinalState>(event, "VFS").particlesByPt();

      Particles leptons = apply<FinalState>(event, "FS_LEPTONS").particlesByPt();
      Particles photons = apply<FinalState>(event, "FS_PHOTONS").particlesByPt();

      std::vector<int> fsr_idx;
      std::vector<Particle> dressed_leptons;
      std::vector<double> dressed_leptons_iso;
      int nFidDressedLeps = 0;

      for(const Particle& lep : vfs) {
        if(!((lep.abspid() == PID::MUON) || (lep.abspid() == PID::ELECTRON) || (lep.abspid() == PID::TAU))) continue;
        if(!(lep.genParticle()->status() == 1 || lep.abspid() == PID::TAU)) continue;
        const Particles& parents = lep.parents();
        if(parents.empty()) continue;

        int mom_id = getTrueMotherPID(lep,0);
        if (!(mom_id == 25 || mom_id == 23 || mom_id == 443 || mom_id == 553 || abs(mom_id) == 24)) continue;        

        FourMomentum lep_dressed = lep.momentum();

        int _fsr_idx = -1;
        for(const Particle& fsr_part : vfs) {
          _fsr_idx += 1;
          if(fsr_part.genParticle()->status() != 1) continue;
          if(fsr_part.pid() != 22) continue;
          const Particles& fsr_parents = fsr_part.parents();
          bool idMatch = false;
          for (const Particle& fsr_parent : fsr_parents) {
            if (fsr_parent.abspid() == lep.abspid()) {
              idMatch = true;
              break;
            }
          }
          if (!idMatch) continue;
          if(deltaR(lep, fsr_part) < 0.3){
            fsr_photons.push_back(fsr_part);
            fsr_idx.push_back(_fsr_idx);
            lep_dressed += fsr_part.momentum();
          }
        }

        nFidDressedLeps += 1;
        int _iso_idx = -1;
        double genIso = 0.0;
        for(const Particle& fs_part : vfs) {
          _iso_idx += 1;
          if(!(fs_part.genParticle()->status() == 1)) continue;
          if((fs_part.abspid() == PID::ELECTRON) || (fs_part.abspid() == PID::MUON)) continue;
          if((!fsr_idx.empty()) && (std::find(fsr_idx.begin(), fsr_idx.end(), _iso_idx) != fsr_idx.end())) continue;
          double dRvL = deltaR(lep_dressed, fs_part);
          if(dRvL < 0.3){
            genIso += fs_part.pT();
          }
        }
        genIso = genIso/lep_dressed.pT();
        Particle dressed_lep = lep;
        dressed_lep.setMomentum(lep_dressed);
        dressed_leptons.push_back(dressed_lep);
        dressed_leptons_iso.push_back(genIso);
      }

      int nFidLeps = 0;
      int nFidPtLead = 0;
      int nFidPtSubLead = 0;

      for (long unsigned int i = 0; i < dressed_leptons.size(); ++i) {
        const Particle& lep = dressed_leptons[i];
        double iso = dressed_leptons_iso[i];

        double pt = lep.pT();
        double eta = lep.eta();
        int pdgid = lep.abspid();

        bool passMu = (pdgid == 13 && pt > 5.0 && std::abs(eta) < 2.4);
        bool passEle = (pdgid == 11 && pt > 7.0 && std::abs(eta) < 2.5);

        if ((passMu || passEle) && iso < 0.35) {
          ++nFidLeps;
          if(pt > 20.0) ++nFidPtLead;
          if(pt > 10.0) ++nFidPtSubLead;
        }
      }

      if(!(nFidLeps>=4 && nFidPtLead>=1 && nFidPtSubLead>=2)) vetoEvent;

      ZMassResult zmass_result = buildZMasses(dressed_leptons, dressed_leptons_iso, true);
      if (!zmass_result.passFidSel) vetoEvent;

      auto [passMassOS, passElMuDeltaR, passDeltaR] = checkEventTopology(dressed_leptons, zmass_result.z_leps_idx);

      if (!(passMassOS && passElMuDeltaR && passDeltaR)) zmass_result.passFidSel = false;
      if (!zmass_result.passFidSel) vetoEvent;

      FourMomentum ZZsystem;
      FourMomentum Z2cand;
      for (int idx : zmass_result.z_leps_idx) {
        ZZsystem += dressed_leptons[idx].momentum();
  if (idx > 1) {
    std::cout << idx << std::endl;
          Z2cand += dressed_leptons[idx].momentum();    
        }
      }

      double m4l = ZZsystem.mass();
      double pT4l = ZZsystem.pT();
      double mz2 = Z2cand.mass();

      if (m4l < 105.0 || m4l > 160.0) vetoEvent;

      ++_nEventsFinal;
      _h_ZZ_pth->fill(pT4l / GeV);
      _h_ZZ_pth_incl->fill(pT4l / GeV);
      _h_ZZ_mz2->fill(mz2 / GeV);
    }

    void finalize(){
      MSG_INFO("Events: " << _nEventsTotal);
      MSG_INFO("Selected: " << _nEventsFinal);
      scale(_h_ZZ_pth, crossSection() / femtobarn * BR / sumOfWeights());
      scale(_h_ZZ_pth_incl, crossSection() / femtobarn * BR / sumOfWeights());
      scale(_h_ZZ_mz2, crossSection() / femtobarn * BR /sumOfWeights());
    }

  private:
    ZMassResult buildZMasses(const std::vector<Particle>& leptons,
                             const std::vector<double>& iso,
                             bool makeCuts);

    std::tuple<double, bool, int, int> buildZ1Mass(const std::vector<Particle>& leptons,
                                                   const std::vector<double>& iso,
                                                   bool makeCuts);

    std::tuple<bool, int, int> buildZ2Mass(const std::vector<Particle>& leptons,
                                           const std::vector<double>& iso,
                                           int idx1, int idx2,
                                           bool makeCuts);

    std::tuple<FourMomentum, FourMomentum> buildLLPair(const Particle& lep1, const Particle& lep2);

    bool checkCuts(const std::vector<Particle>& leptons,
                   const std::vector<double>& iso,
                   int idx1, int idx2);

    std::tuple<bool, bool, bool> checkEventTopology(const std::vector<Particle>& leptons,
                                                    const std::vector<int>& z_leps_idx);

    double sumW_;
    size_t _nEventsTotal;
    size_t _nEventsFinal;
    Histo1DPtr _histo;
    Histo1DPtr _h_ZZ_pth;
    Histo1DPtr _h_ZZ_pth_incl;
    Histo1DPtr _h_ZZ_mz2;
    // H > 4l BR
    const double BR = 0.000128;
    // TODO: Define histos as
    // map<string, Histo1DPtr> _histo;
    // In init: _histo["var"]
  };
  RIVET_DECLARE_PLUGIN(CMS_2025_I2872501);

  std::tuple<FourMomentum, FourMomentum> CMS_2025_I2872501::buildLLPair(const Particle& lep1, const Particle& lep2) {
    return std::make_tuple(lep1.momentum(), lep2.momentum());
  }

  bool CMS_2025_I2872501::checkCuts(const std::vector<Particle>& leptons,
                                    const std::vector<double>& iso,
                                    int idx1, int idx2) {
    const Particle& l1 = leptons[idx1];
    const Particle& l2 = leptons[idx2];
    double iso1 = iso[idx1];
    double iso2 = iso[idx2];

    int id1 = l1.abspid();
    int id2 = l2.abspid();

    bool pass1 = (id1 == PID::MUON && l1.pT() > 5.0 && std::abs(l1.eta()) < 2.4) ||
                 (id1 == PID::ELECTRON && l1.pT() > 7.0 && std::abs(l1.eta()) < 2.5);
    bool pass2 = (id2 == PID::MUON && l2.pT() > 5.0 && std::abs(l2.eta()) < 2.4) ||
                 (id2 == PID::ELECTRON && l2.pT() > 7.0 && std::abs(l2.eta()) < 2.5);

    return pass1 && pass2 && iso1 < 0.35 && iso2 < 0.35;
  }

  ZMassResult CMS_2025_I2872501::buildZMasses(const std::vector<Particle>& leptons,
                                              const std::vector<double>& iso,
                                              bool makeCuts) {
    ZMassResult result;
    result.passFidSel = false;

    auto [offshell, findZ1, idx1, idx2] = buildZ1Mass(leptons, iso, makeCuts);

    if (!findZ1) return result;

    auto [l1, l2] = buildLLPair(leptons[idx1], leptons[idx2]);
    double mZ1 = (l1 + l2).mass();

    bool passZ1 = (!makeCuts) || (mZ1 > MIN_MZ1 && mZ1 < MAX_MZ1);

    auto [findZ2, idx3, idx4] = buildZ2Mass(leptons, iso, idx1, idx2, makeCuts);

    if (passZ1 && findZ2) {
      result.passFidSel = true;
      result.z_leps_idx = {idx1, idx2, idx3, idx4};
    }

    return result;
  }


  std::tuple<double, bool, int, int> CMS_2025_I2872501::buildZ1Mass(const std::vector<Particle>& leptons,
                                 const std::vector<double>& iso,
                                 bool makeCuts) {
    double offshell = 999.0;
    bool findZ1 = false;
    int idx_l1 = -1, idx_l2 = -1;

    for (size_t i = 0; i < leptons.size(); ++i) {
      for (size_t j = i + 1; j < leptons.size(); ++j) {
        if (leptons[i].pid() + leptons[j].pid() != 0) continue;

        if (makeCuts && !checkCuts(leptons, iso, i, j)) continue;

        auto [l_i, l_j] = buildLLPair(leptons[i], leptons[j]);
        double mll = (l_i + l_j).mass();

        if (std::abs(mll - ZMASS) <= offshell) {
          offshell = std::abs(mll - ZMASS);
          idx_l1 = i;
          idx_l2 = j;
          findZ1 = true;
        }
      }
    }

    return {offshell, findZ1, idx_l1, idx_l2};
  }

  std::tuple<bool, int, int> CMS_2025_I2872501::buildZ2Mass(const std::vector<Particle>& leptons,
                                 const std::vector<double>& iso,
                                 int idx1, int idx2,
                                 bool makeCuts) {
    bool findZ2 = false;
    int idx_l3 = -1, idx_l4 = -1;
    double maxPtSum = 0.0;

    for (size_t i = 0; i < leptons.size(); ++i) {
      if (int(i) == idx1 || int(i) == idx2) continue;

      for (size_t j = i + 1; j < leptons.size(); ++j) {
        if (int(j) == idx1 || int(j) == idx2) continue;
        if (leptons[i].pid() + leptons[j].pid() != 0) continue;

        if (makeCuts && !checkCuts(leptons, iso, i, j)) continue;

        auto [l_i, l_j] = buildLLPair(leptons[i], leptons[j]);
        FourMomentum Z2 = l_i + l_j;
        double ptSum = l_i.pT() + l_j.pT();
        double mass_Z2 = Z2.mass();

        if (ptSum >= maxPtSum) {
          if ((mass_Z2 >= MIN_MZ2 && mass_Z2 <= MAX_MZ2) || !makeCuts) {
            idx_l3 = i;
            idx_l4 = j;
            findZ2 = true;
            maxPtSum = ptSum;
          } else if (!findZ2) {
            idx_l3 = i;
            idx_l4 = j;
          }
        }
      }
    }

    return {findZ2, idx_l3, idx_l4};
  }

  std::tuple<bool, bool, bool> CMS_2025_I2872501::checkEventTopology(const std::vector<Particle>& leptons,
                                        const std::vector<int>& z_leps_idx) {

    bool passedMassOS = true;
    bool passedElMuDeltaR = true;
    bool passedDeltaR = true;

    for (size_t i = 0; i < leptons.size(); ++i) {
      if (std::find(z_leps_idx.begin(), z_leps_idx.end(), i) == z_leps_idx.end()) continue;

      for (size_t j = i + 1; j < leptons.size(); ++j) {
        if (std::find(z_leps_idx.begin(), z_leps_idx.end(), j) == z_leps_idx.end()) continue;

        const Particle& l1 = leptons[i];
        const Particle& l2 = leptons[j];

        FourMomentum l1mom = l1.momentum();
        FourMomentum l2mom = l2.momentum();
        FourMomentum mll = l1mom + l2mom;

        if ((l1.pid() * l2.pid() < 0) && (mll.mass() <= 4.0)) {
          passedMassOS = false;
          break;
        }

        double dRll = deltaR(l1, l2);

        if (std::abs(l1.abspid()) != std::abs(l2.abspid())) {
          if (dRll <= 0.02) {
            passedElMuDeltaR = false;
            break;
          }
        }

        if (dRll <= 0.02) {
          passedDeltaR = false;
          break;
        }
      }

      if (!passedMassOS || !passedElMuDeltaR || !passedDeltaR) break;
    }

    return std::make_tuple(passedMassOS, passedElMuDeltaR, passedDeltaR);
  }

}