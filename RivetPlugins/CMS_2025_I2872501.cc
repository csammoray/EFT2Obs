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
      // sumW_ = 0.;
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

      //---Jets
      FastJets fs_jets(fs, FastJets::ANTIKT, 0.4);
      declare(fs_jets, "JETS");

      book(_h_sigma, "h_sigma", 1, 0, 2); // Histogram to get the cross-section without fiducial cuts
      book(_h_ZZ_incl, "incl", 1, 0, 10000);
      book(_h_ZZ_pth, "pt_h", {0,10,16,22,28,36,46,60,80,106,146,10000});
      book(_h_ZZ_mz2, "m_z2", {12,22,26,28,32,34,40,50,65});
      book(_h_ZZ_mz2_incl, "m_z2_incl", 1, 0, 65);
      book(_h_ZZ_dphi_jj, "dphi_jj", {-100, -M_PI, -M_PI/2, 0, M_PI/2, M_PI});
      book(_h_ZZ_deta_jj, "deta_jj", {-100,0,1.1,2.9,4.4,10});
      book(_h_ZZ_phi, "phi", {-M_PI, -3*M_PI/4, -M_PI/2, -M_PI/4, 0, M_PI/4, M_PI/2, 3*M_PI/4, M_PI});
      book(_h_ZZ_y4l_pT4l, "yh_pt_h", 12, 0, 12);
      book(_h_ZZ_mz1_mz2, "mz1_mz2", 7, 0, 7);
    }

    void analyze(const Event& event) {
      _h_sigma->fill(1.0);

      auto jets = apply<JetAlg>(event, "JETS").jetsByPt(Cuts::abseta < 4.7 && Cuts::pT > 30*GeV);

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
      for (int idx : zmass_result.z_leps_idx) {
        ZZsystem += dressed_leptons[idx].momentum();
      }

      FourMomentum Z1cand;
      int idx1 = zmass_result.z_leps_idx[0];
      int idx2 = zmass_result.z_leps_idx[1];
      Z1cand += dressed_leptons[idx1].momentum();
      Z1cand += dressed_leptons[idx2].momentum();

      FourMomentum Z2cand;
      int idx3 = zmass_result.z_leps_idx[2];
      int idx4 = zmass_result.z_leps_idx[3];
      Z2cand += dressed_leptons[idx3].momentum();
      Z2cand += dressed_leptons[idx4].momentum();
      
      double m4l = ZZsystem.mass();
      double pT4l = ZZsystem.pT();
      double y4l = fabs(ZZsystem.rapidity());
      double mz1 = Z1cand.mass();
      double mz2 = Z2cand.mass();

      if (mz2 < MIN_MZ2) vetoEvent;
      if (m4l < 105.0 || m4l > 160.0) vetoEvent;

      ++_nEventsFinal;
      _h_ZZ_incl->fill(pT4l / GeV);
      _h_ZZ_pth->fill(pT4l / GeV);
      _h_ZZ_mz2->fill(mz2 / GeV);
      _h_ZZ_mz2_incl->fill(mz2 / GeV);
      _h_ZZ_phi->fill(computePhi(dressed_leptons, zmass_result.z_leps_idx));

      if(jets.size() > 1) {
        _h_ZZ_dphi_jj->fill(deltaphi_jj(jets[0].momentum(), jets[1].momentum()));
        _h_ZZ_deta_jj->fill(fabs(deltaEta(jets[0], jets[1])));
      }
      else {
        _h_ZZ_dphi_jj->fill((-4.0)); // Underflow bin for events with 0 or 1 jet
        _h_ZZ_deta_jj->fill((-0.5)); // Underflow bin for events with 0 or 1 jet
      }

      // Fill 2D histogram y4l vs pT4l
      if (pT4l < 50.0) {
        if      (y4l < 0.2) _h_ZZ_y4l_pT4l->fill(0.5);
        else if (y4l < 0.4) _h_ZZ_y4l_pT4l->fill(1.5);
        else if (y4l < 0.65) _h_ZZ_y4l_pT4l->fill(2.5);
        else if (y4l < 0.9) _h_ZZ_y4l_pT4l->fill(3.5);
        else if (y4l < 1.2) _h_ZZ_y4l_pT4l->fill(4.5);
        else if (y4l <= 2.5) _h_ZZ_y4l_pT4l->fill(5.5);
      }
      else if (pT4l < 105.0) {
        if      (y4l < 0.5) _h_ZZ_y4l_pT4l->fill(6.5);
        else if (y4l < 1.15) _h_ZZ_y4l_pT4l->fill(7.5);
        else if (y4l <= 2.5) _h_ZZ_y4l_pT4l->fill(8.5);
      }
      else if (pT4l < 10000.0) {
        if      (y4l < 0.45) _h_ZZ_y4l_pT4l->fill(9.5);
        else if (y4l < 1.0) _h_ZZ_y4l_pT4l->fill(10.5);
        else if (y4l <= 2.5) _h_ZZ_y4l_pT4l->fill(11.5);
      }

      // Fill 2D histogram mz1 vs mz2
      if (mz1 >= 40.0 && mz1 < 88.0) {
          if (mz2 >= 12.0 && mz2 < 28.0) _h_ZZ_mz1_mz2->fill(0.5);
          else if (mz2 < 34.0) _h_ZZ_mz1_mz2->fill(1.5);
          else if (mz2 < 40.0) _h_ZZ_mz1_mz2->fill(2.5);
          else if (mz2 < 65.0) _h_ZZ_mz1_mz2->fill(3.5);
      }
      else if (mz1 >= 88.0 && mz1 < 120.0) {
          if (mz2 >= 12.0 && mz2 < 25.0) _h_ZZ_mz1_mz2->fill(4.5);
          else if (mz2 < 28.0) _h_ZZ_mz1_mz2->fill(5.5);
          else if (mz2 < 65.0) _h_ZZ_mz1_mz2->fill(6.5);
      
      }
    }    
    void finalize(){
      MSG_INFO("Events: " << _nEventsTotal);
      MSG_INFO("Selected: " << _nEventsFinal);
      // scale(_histo, crossSection() / femtobarn * BR / sumOfWeights());
      scale(_h_sigma, crossSection() / femtobarn * BR / sumOfWeights());
      scale(_h_ZZ_incl, crossSection() / femtobarn * BR / sumOfWeights());
      scale(_h_ZZ_pth, crossSection() / femtobarn * BR / sumOfWeights());
      scale(_h_ZZ_mz2, crossSection() / femtobarn * BR / sumOfWeights());
      scale(_h_ZZ_mz2_incl, crossSection() / femtobarn * BR / sumOfWeights());
      scale(_h_ZZ_dphi_jj, crossSection() / femtobarn * BR / sumOfWeights());
      scale(_h_ZZ_deta_jj, crossSection() / femtobarn * BR / sumOfWeights());
      scale(_h_ZZ_phi, crossSection() / femtobarn * BR / sumOfWeights());
      scale(_h_ZZ_y4l_pT4l, crossSection() / femtobarn * BR / sumOfWeights());
      scale(_h_ZZ_mz1_mz2, crossSection() / femtobarn * BR / sumOfWeights());
    }

  private:
    static void constrainedRemovePairMass(FourMomentum& p1, FourMomentum& p2, double m1=0.0, double m2=0.0) {
      const FourMomentum nullp(0, 0, 0, 0);
      if (p1 == nullp || p2 == nullp) return;

      const FourMomentum p1old = p1;
      const FourMomentum p2old = p2;

      const FourMomentum p12 = p1old + p2old;
      const FourMomentum diff_p2p1 = p2old - p1old;

      const double p1sq = p1old.mass2();
      const double p2sq = p2old.mass2();
      const double p1p2 = p1old.dot(p2old);
      const double m1sq = m1 * fabs(m1);
      const double m2sq = m2 * fabs(m2);
      const double p12sq = p12.mass2();

      FourMomentum avec = p2old;
      avec *= p1sq;
      FourMomentum tmp = p1old;
      tmp *= p2sq;
      avec -= tmp;
      tmp = diff_p2p1;
      tmp *= p1p2;
      avec += tmp;

      const double a = avec.mass2();
      const double b = (p12sq + m2sq - m1sq) * (p1p2*p1p2 - p1sq*p2sq);
      const double c = 0.25 * std::pow(p12sq + m2sq - m1sq, 2) * p1sq - std::pow(p1sq + p1p2, 2) * m2sq;

      const double eta = (-b - std::sqrt(fabs(b*b - 4.0 * a * c))) / (2.0 * a);
      const double xi = (p12sq + m2sq - m1sq - 2.0 * eta * (p2sq + p1p2)) / (2.0 * (p1sq + p1p2));

      FourMomentum p1hat = p1old;
      p1hat *= (1. - xi);
      tmp = p2old;
      tmp *= (1. - eta);
      p1hat += tmp;

      FourMomentum p2hat = p1old;
      p2hat *= xi;
      tmp = p2old;
      tmp *= eta;
      p2hat += tmp;

      p1 = p1hat;
      p2 = p2hat;
    }

    double computePhi(const std::vector<Particle>& leptons, const std::vector<int>& z_leps_idx) {
      if (z_leps_idx.size() != 4) return -999;

      const Particle& a1 = leptons[z_leps_idx[0]];
      const Particle& a2 = leptons[z_leps_idx[1]];
      const Particle& b1 = leptons[z_leps_idx[2]];
      const Particle& b2 = leptons[z_leps_idx[3]];

      // Sorting the leptons
      FourMomentum Z1_neg_lep = (a1.pid() > 0) ? a1.momentum() : a2.momentum();
      FourMomentum Z1_pos_lep = (a1.pid() > 0) ? a2.momentum() : a1.momentum();

      FourMomentum Z2_neg_lep = (b1.pid() > 0) ? b1.momentum() : b2.momentum();
      FourMomentum Z2_pos_lep = (b1.pid() > 0) ? b2.momentum() : b1.momentum();

      constrainedRemovePairMass(Z1_neg_lep, Z1_pos_lep, 0.0, 0.0);
      constrainedRemovePairMass(Z2_neg_lep, Z2_pos_lep, 0.0, 0.0);


      const FourMomentum Z1 = Z1_neg_lep + Z1_pos_lep;
      const FourMomentum Z2 = Z2_neg_lep + Z2_pos_lep;
      const FourMomentum H  = Z1 + Z2;

      const LorentzTransform boost_H_restframe = LorentzTransform::mkFrameTransformFromBeta(H.betaVec());

      const FourMomentum Z1_neg_lep_boosted = boost_H_restframe.transform(Z1_neg_lep);
      const FourMomentum Z1_pos_lep_boosted = boost_H_restframe.transform(Z1_pos_lep);
      const FourMomentum Z2_neg_lep_boosted = boost_H_restframe.transform(Z2_neg_lep);
      const FourMomentum Z2_pos_lep_boosted = boost_H_restframe.transform(Z2_pos_lep);

      const FourMomentum Z1_boosted = Z1_neg_lep_boosted + Z1_pos_lep_boosted;

      const Vector3 Z1_p3 = Z1_boosted.vector3();
      const Vector3 n1 = Z1_neg_lep_boosted.vector3().cross(Z1_pos_lep_boosted.vector3());
      const Vector3 n2 = Z2_neg_lep_boosted.vector3().cross(Z2_pos_lep_boosted.vector3());

      if (Z1_p3.mod() == 0.0 || n1.mod() == 0.0 || n2.mod() == 0.0)
        return -999;

      const Vector3 Z1_p3_norm = Z1_p3.unit();
      const Vector3 n1_norm = n1.unit();
      const Vector3 n2_norm = n2.unit();

      const double cosPhi = std::max(-1.0, std::min(1.0, -n1_norm.dot(n2_norm)));
      const double sinPhi = Z1_p3_norm.dot(n1_norm.cross(n2_norm));

      return std::atan2(sinPhi, cosPhi); 
    }

    double deltaphi_jj(const FourMomentum& h1, const FourMomentum& h2) {
      //Direction of the two jets - vectors in the lab frame
      Vector3 j1dir(h1.x(), h1.y(), h1.z());
      Vector3 j2dir(h2.x(), h2.y(), h2.z());
      //Transverse component in the xy plane
      Vector3 jt1(h1.x(), h1.y(), 0);
      Vector3 jt2(h2.x(), h2.y(), 0);
      //Unit vectors of the transverse components
      Vector3 jt1_norm   = jt1 * (1/jt1.mod());
      Vector3 jt2_norm   = jt2 * (1/jt2.mod());
      //Unit vector of the z axis
      Vector3 z(0,0,1);
      //Cross product between transverse components
      double cross      = jt1_norm.cross(jt2_norm).dot(z);
      double cross_norm = cross * (1 / abs(cross));
      //Dot product between transverse components
      double dot         = jt1_norm.dot(jt2_norm);
      //Difference between the direction of the two jets
      double diff       = (j1dir - j2dir).dot(z);
      double diff_norm  = diff * (1 / abs(diff));
      return acos(dot) * diff_norm * cross_norm;
    }

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

    // double sumW_;
    size_t _nEventsTotal;
    size_t _nEventsFinal;
    // Histo1DPtr _histo;
    Histo1DPtr _h_sigma;
    Histo1DPtr _h_ZZ_incl;
    Histo1DPtr _h_ZZ_pth;
    Histo1DPtr _h_ZZ_mz2;
    Histo1DPtr _h_ZZ_mz2_incl;
    Histo1DPtr _h_ZZ_dphi_jj;
    Histo1DPtr _h_ZZ_deta_jj;
    Histo1DPtr _h_ZZ_phi;
    Histo1DPtr _h_ZZ_y4l_pT4l;
    Histo1DPtr _h_ZZ_mz1_mz2;
    // H > 4l BR
    const double BR = 0.000128;
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
