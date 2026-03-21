//
// APFEL++ 2017
//
// Authors: Valerio Bertone: valerio.bertone@cern.ch
//

#include <LHAPDF/LHAPDF.h>
#include <apfel/apfelxx.h>

// b* prescription
double bstar(double const& b)
{
  const double b0 = 2 * exp( - apfel::emc);
  return b / sqrt( 1 + pow(b / b0, 2) );
}

// mu* prescription
double mustar(double const& mu)
{
  const double b0 = 2 * exp( - apfel::emc);
  return b0 / bstar(b0 / mu);
}

// Non-perturbative function
double fNP(double const&, double const& b, double const& zetaf)
{
  const double Q02 = 3.2;
  const double g1  = 0;//0.20;
  const double g2  = 0;//0.09;
  return exp( - ( g1 + g2 * log(zetaf / Q02) / 2 ) * pow(b, 2) / 2 );
}

// Main program
int main()
{
  // Open LHAPDF set
  LHAPDF::pathsPrepend(".");
  const std::string ffset = "DSS14_NLO_Pip";
  LHAPDF::PDF* distff = LHAPDF::mkPDF(ffset, 0);

  // Heavy-quark thresholds
  const std::vector<double> Thresholds{0, 0, 0, 0, 0};
  //for (auto const& v : distff->flavors())
  //  if (v > 0 && v < 7)
  //    Thresholds.push_back(distff->quarkThreshold(v));

  // x-space grid (to setup APFEL++ computation)
  const apfel::Grid g{{{100, 1e-2, 3}, {80, 6e-1, 3}, {50, 8e-1, 3}}};

  // Rotate FF set into the QCD evolution basis
  const auto RotFFs = [=] (double const& x, double const& mu) -> std::map<int,double>{ return apfel::PhysToQCDEv(distff->xfxQ(x,mu)); };

  // Construct set of distributions as a function of the scale to be
  // tabulated
  const auto EvolvedFFs = [=,&g] (double const& mu) -> apfel::Set<apfel::Distribution>
    {
      return apfel::Set<apfel::Distribution>{apfel::EvolutionBasisQCD{apfel::NF(mu, Thresholds)}, DistributionMap(g, RotFFs, mu)};
    };

  // Tabulate collinear FFs
  const apfel::TabulateObject<apfel::Set<apfel::Distribution>> TabFFs{EvolvedFFs, 100, 1, distff->qMax(), 3, Thresholds};
  const auto CollFFs = [&] (double const& mu) -> apfel::Set<apfel::Distribution> { return TabFFs.Evaluate(mustar(mu)); };

  // Initialize TMD objects
  const auto TmdObj = apfel::InitializeTmdObjects(g, Thresholds);

  // Alpha_em
  const double aref = 0.007555310522369057;
  const double Qref = 91.1876;

  // Perturbative order
  const int PerturbativeOrder = 3;

  // Kinematics
  const double yb = 1;
  const double Qb = 91.1876;
  const double z1 = 0.7;
  const double z2 = 0.7;
  const int nqT = 100;
  const double qTmin = 0.001;
  const double qTmax = 10;
  const double qTstp = ( qTmax - qTmin ) / ( nqT - 1 );

  // Double exponential quadrature
  apfel::DoubleExponentialQuadrature DEObj{};

  // Alpha_s (from FFs)
  const auto Alphas = [&] (double const& mu) -> double{ return distff->alphasQ(mu); };

  // Build evolved TMD FFs
  const auto EvTMDFFs = BuildTmdFFs(TmdObj, CollFFs, Alphas, PerturbativeOrder, 1);

  // Get hard-factor function
  const auto Hf = HardFactor("DY", TmdObj, Alphas, PerturbativeOrder, 1);

  // Renormalisation scales
  const double muf   = Qb;
  const double zetaf = Qb * Qb;
	  
  // Number of active flavours at 'muf'
  const int nf = apfel::NF(muf, Thresholds);

  // EW charges
  const std::vector<double> Bq = apfel::ElectroWeakCharges(Qb, true);

  // Electromagnetic coupling squared
  const double aem2 = pow(aref, 2);
	  
  // Compute the hard factor
  const double hcs = Hf(muf);

  // Construct the TMD luminosities in b space to be fed to be
  // trasformed in qT space.
  // Opposite sign
  const std::function<double(double const&)> TMDLumibOS = [=] (double const& b) -> double
    {
      // Get Evolved TMD FFs and rotate them into the physical
      // basis
      const std::map<int,apfel::Distribution> xF = QCDEvToPhys(EvTMDFFs(bstar(b), muf, zetaf).GetObjects());

      // Combine TMDs through the EW charges
      double lumi = 0;
      for (int i = 1; i <= nf; i++)
	lumi += Bq[i-1] * ( xF.at(i).Evaluate(z1) * xF.at(i).Evaluate(z2) + xF.at(-i).Evaluate(z1) * xF.at(-i).Evaluate(z2) );

      // Combine all pieces and return
      return b * lumi * fNP(z1, b, zetaf) * fNP(z2, b, zetaf);
    };
  // Same sign
  const std::function<double(double const&)> TMDLumibSS = [=] (double const& b) -> double
    {
      // Get Evolved TMD FFs and rotate them into the physical
      // basis
      const std::map<int,apfel::Distribution> xF = QCDEvToPhys(EvTMDFFs(bstar(b), muf, zetaf).GetObjects());

      // Combine TMDs through the EW charges
      double lumi = 0;
      for (int i = 1; i <= nf; i++)
	lumi += Bq[i-1] * ( xF.at(i).Evaluate(z1) * xF.at(-i).Evaluate(z2) + xF.at(-i).Evaluate(z1) * xF.at(i).Evaluate(z2) );

      // Combine all pieces and return
      return b * lumi * fNP(z1, b, zetaf) * fNP(z2, b, zetaf);
    };
  std::cout << "#    qT            OS            SS" << std::endl;
  for (double qT = qTmin; qT <= 1.000001 * qTmax; qT += qTstp)
    {
      // Perform Fourier transform and obtain cross section
      const double prefactor = apfel::ConvFact * qT * 12 * M_PI * aem2 * hcs / pow(Qb, 2) * ( 1 + pow(1 - yb, 2) ) / 2 / z1 / z2;

      // Alternative integration
      const double deOS = prefactor * DEObj.transform(TMDLumibOS, qT);
      const double deSS = prefactor * DEObj.transform(TMDLumibSS, qT);

      // Print results
      std::cout << std::scientific << qT << "  " << deOS << "  " << deSS << std::endl;
    }
  std::cout << std::endl;

  const int nz = 10;
  const double zmin = 0.01;
  const double zmax = 0.9;
  const double zstp = exp(log(zmax / zmin) / ( nz - 1 ));
  std::cout << "#    z          d(u)/d(ubar) " << std::endl;
  for (double z = zmin; z <= 1.000001 * zmax; z *= zstp)
    std::cout << z << "\t" << distff->xfxQ(2, z, Qb) / distff->xfxQ(-2, z, Qb) << std::endl;
  std::cout << std::endl;

  delete distff;
  return 0;
}
