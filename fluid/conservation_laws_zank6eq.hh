/*!
\file conservation_laws_zank6eq.hh
\brief Declares rules for Zank's 6-equation turbulence model
\author Vladimir Florinski

This file is part of the SPECTRUM suite of scientific numerical simulation codes. SPECTRUM stands for Space Plasma and Energetic Charged particle TRansport on Unstructured Meshes. The code simulates plasma or neutral particle flows using MHD equations on a grid, transport of cosmic rays using stochastic or grid based methods. The "unstructured" part refers to the use of a geodesic mesh providing a uniform coverage of the surface of a sphere.
*/

#ifndef SPECTRUM_CONSERVATION_LAWS_ZANK6EQ_HH
#define SPECTRUM_CONSERVATION_LAWS_ZANK6EQ_HH

namespace Spectrum {

#define CL_ZANK6EQ_NVARS 6

struct PrimitiveStateZank6eq;
struct ConservedStateZank6eq;
struct FluxFunctionZank6eq;
struct SourceFunctionZank6eq;

//----------------------------------------------------------------------------------------------------------------------------------------------------
// PrimitiveStateZank6eq class declaration
//----------------------------------------------------------------------------------------------------------------------------------------------------

/*!
\brief Codifies physics for Zank et al. 6-equation turbulence model
\author Keyvan Ghanbari
\author Vladimir Florinski
*/
struct PrimitiveStateZank6eq : public SimpleArray<double, CL_ZANK6EQ_NVARS>
{
   using SimpleArray<double, CL_ZANK6EQ_NVARS>::data;

//! A trait to be used in template specializations
   static constexpr bool is_cons_law_active = false;

//! Alias for total turbulent energy (RO)
   const double& tte(void) const {return data[0];};

//! Alias for total turbulent energy (RW)
   double& tte(void) {return data[0];};

//! Alias for cross helicity (RO)
   const double& che(void) const {return data[1];};

//! Alias for cross helicity (RW)
   double& che(void) {return data[1];};

//! Alias for spatial integral of forward energy L+ (RO)
   const double& sif(void) const {return data[2];};

//! Alias for spatial integral of forward energy L+ (RW)
   double& sif(void) {return data[2];};

//! Alias for spatial integral of backward energy L- (RO)
   const double& sib(void) const {return data[3];};

//! Alias for spatial integral of backward energy L- (RW)
   double& sib(void) {return data[3];};

//! Alias for residual energy (RO)
   const double& ren(void) const {return data[4];};

//! Alias for residual energy (RW)
   double& ren(void) {return data[4];};

//! Alias for spatial integral of residual energy L_d (RO)
   const double& sir(void) const {return data[5];};

//! Alias for spatial integral of residual energy L_d (RO)
   double& sir(void) {return data[5];};

//! Return the number of variables
   SPECTRUM_DEVICE_FUNC static constexpr int Nvars(void) {return CL_PASSIVE_NVARS;};

//! Default constructor
   SPECTRUM_DEVICE_FUNC PrimitiveStateZank6eq(void) = default;

//! Constructor from components
   SPECTRUM_DEVICE_FUNC PrimitiveStateZank6eq(double tte_in, double che_in, double sif_in, double sib_in, double ren_in, double sir_in);

//! Constructor from the base class
   SPECTRUM_DEVICE_FUNC PrimitiveStateZank6eq(const SimpleArray<double, CL_ZANK6EQ_NVARS>& other);

//! Calculate conserved state
   template <typename cl_prim, std::enable_if_t<cl_prim::is_cons_law_active, bool> = true>
   SPECTRUM_DEVICE_FUNC PrimitiveStateZank6eq ToConserved(cl_prim active_prim) const;

//! Calculate flux function
   template <typename cl_prim, std::enable_if_t<cl_prim::is_cons_law_active, bool> = true>
   SPECTRUM_DEVICE_FUNC FluxFunctionMHDturb ToFlux(cl_prim active_prim) const;
};

/*!
\brief Codifies physics for Zank et al. 6-equation turbulence model
\author Keyvan Ghanbari
\author Vladimir Florinski
*/
struct ConservedStateMHDturb :  public SimpleArray<double, CL_ZANK6EQ_NVARS>
{
   using SimpleArray<double, CL_ZANK6EQ_NVARS>::data;

//! A trait to be used in template specializations
   static constexpr bool is_cons_law_active = false;

//! Alias for conserved total turbulent energy (RO)
   const double& ctte(void) const {return data[0];};

//! Alias for conserved total turbulent energy (RW)
   double& ctte(void) {return data[0];};

//! Alias for conserved cross helicity (RO)
   const double& cche(void) const {return data[1];};

//! Alias for conserved cross helicity (RO)
   double& cche(void) {return data[1];};

//! Alias for conserved spatial integral of forward energy (L+)
   const double& csif(void) const {return data[11];};
         double& csif(void)       {return data[11];};

//! Alias for conserved spatial integral of backward energy (L-)
   const double& csib(void) const {return data[12];};
         double& csib(void)       {return data[12];};

//! Alias for conserved residual energy (E_d)
   const double& cren(void) const {return data[13];};
         double& cren(void)       {return data[13];};

//! Alias for conserved spatial integral of residual energy (L_d)
   const double& csir(void) const {return data[14];};
         double& csir(void)       {return data[14];};

//! Calculate primitive state
   SPECTRUM_DEVICE_FUNC PrimitiveStateMHDturb ToPrimitive(void) const;
};

struct FluxFunctionMHDturb : public CLStateBase<CL_MHD_TURB_TOTAL> {

//! Alias for density flux
   const double& den(void) const {return data[0];};
         double& den(void)       {return data[0];};

//! Alias for momentum flux
   const GeoVector& mom(void) const {return (const GeoVector&)data[1];};
         GeoVector& mom(void)       {return (GeoVector&)data[1];};

//! Alias for energy flux
   const double& enr(void) const {return data[4];};
         double& enr(void)       {return data[4];};

//! Alias for magnetic flux
   const GeoVector& mag(void) const {return (const GeoVector&)data[5];};
         GeoVector& mag(void)       {return (GeoVector&)data[5];};

#if CL_MHD_GLM != 0
//! Alias for Largange multiplier
   const double& glm(void) const {return data[8];};
         double& glm(void)       {return data[8];};
#endif

//! Alias for total turbulent energy flux
   const double& ttef(void) const {return data[9];};
         double& ttef(void)       {return data[9];};
         
//! Alias for cross helicity flux
   const double& chef(void) const {return data[10];};
         double& chef(void)       {return data[10];};
         
//! Alias for spatial integral of forward energy flux
   const double& siff(void) const {return data[11];};
         double& siff(void)       {return data[11];};
         
//! Alias for spatial integral of backward energy flux
   const double& sibf(void) const {return data[12];};
         double& sibf(void)       {return data[12];};
         
//! Alias for residual energy flux
   const double& renf(void) const {return data[13];};
         double& renf(void)       {return data[13];};
         
//! Alias for  spatial integral of residual energy flux
   const double& sirf(void) const {return data[14];};
         double& sirf(void)       {return data[14];};

};

/*!
\brief Codifies physics for Zank et al. 6-equation turbulence model
\author Keyvan Ghanbari
\author Vladimir Florinski
*/
struct SourceFunctionMHDturb : public CLStateBase<CL_TURB_NVARS> {

//! Alias for total turbulent energy source
   const double& ttes(void) const {return data[0];};
         double& ttes(void)       {return data[0];};
         
//! Alias for cross helicity source
   const double& ches(void) const {return data[1];};
         double& ches(void)       {return data[1];};
         
//! Alias for spatial integral of forward energy source
   const double& sifs(void) const {return data[2];};
         double& sifs(void)       {return data[2];};
         
//! Alias for spatial integral of backward energy source
   const double& sibs(void) const {return data[3];};
         double& sibs(void)       {return data[3];};
         
//! Alias for residual energy source
   const double& rens(void) const {return data[4];};
         double& rens(void)       {return data[4];};
         
//! Alias for  spatial integral of residual energy source
   const double& sirs(void) const {return data[5];};
         double& sirs(void)       {return data[5];};

};
//----------------------------------------------------------------------------------------------------------------------------------------------------
// MHD inline methods
//----------------------------------------------------------------------------------------------------------------------------------------------------

/*!
\author Vladimir Florinski
\date 03/06/2024
\return Fastest speed (fast magnetosonic)
*/
SPECTRUM_DEVICE_FUNC inline double PrimitiveStateMHDturb::FastestSpeed(void) const
{
   return sqrt(Alfven2(den(), mag().Norm2()) + Sound2(den(), pre(), CL_MHD_FLUID));
};

/*!
\author Vladimir Florinski
\date 03/06/2024
\return Vector of conserved variables
*/
SPECTRUM_DEVICE_FUNC inline ConservedStateMHDturb PrimitiveStateMHDturb::ToConserved(void) const
{
   ConservedStateMHDturb cons;
   cons.den() = den();
   cons.mom() = vel() * den();
   cons.enr() = Energy(den(), vel().Norm2(), mag().Norm2(), pre(), CL_MHD_FLUID);
   cons.mag() = mag();

#if CL_MHD_GLM != 0
   cons.glm() = glm();
#endif

   cons.ctte() = den() * tte();
   cons.cche() = den() * che();
   cons.csif() = den() * sif();
   cons.csib() = den() * sib();
   cons.cren() = den() * ren();
   cons.csir() = den() * sir();

// Indicator variables
   for(auto i = CL_MHD_NVARS + CL_TURB_NVARS + CL_MHD_GLM; i < CL_MHD_NVARS + CL_TURB_NVARS + CL_MHD_GLM + CL_MHD_IND; i++) {
      cons.data[i] = den() * data[i];
   };
   return cons;
};

/*!
\author Vladimir Florinski
\date 03/06/2024
\return Flux vector
*/
SPECTRUM_DEVICE_FUNC inline FluxFunctionMHDturb PrimitiveStateMHDturb::ToFlux(void) const
{
   double enr = Energy(den(), vel().Norm2(), mag().Norm2(), pre(), CL_MHD_FLUID);
   double pre_tot = pre() + mag().Norm2() / eightpi;
   FluxFunctionMHDturb flux;
   flux.den() = den() * vel()[0];
   flux.mom() = den() * vel()[0] * vel() + pre_tot * gv_nx - (mag()[0] / fourpi) * mag();
   flux.enr() = (enr + pre_tot) * vel()[0] - (vel() * mag()) * mag()[0] / fourpi;
   flux.mag() = vel()[0] * mag() - mag()[0] * vel();

#if CL_MHD_GLM != 0
// The GLM flux is proportional to the square of the fastest speed. The safety factor "CL_MHD_GLMSAFETY" is used to make sure it is slightly outside of the "normal" Riemann fan.
   flux.mag() += glm() * gv_nx;
   flux.glm() = Sqr(CL_MHD_GLMSAFETY * FastestSpeed()) * mag()[0];
#endif
   
   flux.ttef() = tte() * den() * vel()[0];
   flux.chef() = che() * den() * vel()[0];
   flux.siff() = sif() * den() * vel()[0];
   flux.sibf() = sib() * den() * vel()[0];
   flux.renf() = ren() * den() * vel()[0];
   flux.sirf() = sir() * den() * vel()[0];

// Indicator variables
   for(auto i = CL_MHD_NVARS + CL_TURB_NVARS + CL_MHD_GLM; i < CL_MHD_NVARS + CL_TURB_NVARS + CL_MHD_GLM + CL_MHD_IND; i++) {
      flux.data[i] = den() * data[i] * vel()[0];
   };
   return flux;
};
/*!
\author Vladimir Florinski - Keyvan Ghanbari
\date 03/06/2024 - 03/13/2024
\return source terms ONY for turbulence conservation laws
*/
SPECTRUM_DEVICE_FUNC inline FluxFunctionMHDturb PrimitiveStateMHDturb::ToSource(void) const
{
   double divu, mixing, nlt_p, nlt_m, karman_taylor_constant;
   PrimitiveStateMHDturb prim;
   // Warning::: where and how to define/compute the following three parameters? ????????????????????
   karman_taylor_constant = 0.05; //  From CIR simulation
   divu = 0.0;                    // (∇ · u)
   mixing = 0.0;                  // M = (∇ · u)/2 − b̂ · (b̂ · ∇)u

   nlt_p = tte() + che();
   nlt_m = tte() - che();
   SourceFunctionMHDturb source;
   source.ttes() = -0.5 * divu * den() * tte() - mixing * den() * res() - karman_taylor_constant * den() * (Sqr(nlt_p)*sqrt(nlt_m)/sif() + Sqr(nlt_m)*sqrt(nlt_p)/sib() );
   source.ches() = -0.5 * divu * den() * che()                          - karman_taylor_constant * den() * (Sqr(nlt_p)*sqrt(nlt_m)/sif() - Sqr(nlt_m)*sqrt(nlt_p)/sib() );
   source.sifs() = -0.5 * divu * den() * sif() - 0.5 * mixing * den() * sir();
   source.sibs() = -0.5 * divu * den() * sib() - 0.5 * mixing * den() * sir();
   source.rens() = -0.5 * divu * den() * ren() - mixing * den() * tte() - karman_taylor_constant * den() * ren() * (nlt_p * sqrt(nlt_p)/sif() + nlt_m * sqrt(nlt_p)/sib() ); 
   source.sirs() = -0.5 * divu * den() * sir() - mixing * den() * (sif() + sib());

   return source;
};
/*!
\author Vladimir Florinski
\date 03/06/2024
\return Vector of primitive variables
*/
SPECTRUM_DEVICE_FUNC inline PrimitiveStateMHDturb ConservedStateMHDturb::ToPrimitive(void) const
{
   PrimitiveStateMHDturb prim;
   prim.den() = den();
   prim.vel() = mom() / den();
   prim.pre() = Pressure(den(), prim.vel().Norm2(), mag().Norm2(), enr(), CL_MHD_FLUID);
   prim.mag() = mag();

#if CL_MHD_LAGR != 0
   prim.glm() = glm();
#endif

   prim.tte() = ctte() / den();
   prim.che() = cche() / den();
   prim.sif() = csif() / den();
   prim.sib() = csib() / den();
   prim.ren() = cren() / den();
   prim.sir() = csir() / den();

// Indicator variables
   for(auto i = CL_MHD_NVARS + CL_TURB_NVARS + CL_MHD_GLM; i < CL_MHD_NVARS + CL_TURB_NVARS + CL_MHD_GLM + CL_MHD_IND; i++) {
      prim.data[i] = data[i] / den();
   };
   return prim;
};

};

#endif
