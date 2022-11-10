/*--------------------------------------------------------------------
REEF3D
Copyright 2021 Hans Bihs
Copyright 2020 SINTEF Ocean

This file is part of REEF3D.

REEF3D is free software; you can redistribute it and/or modify it
under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License
for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, see <http://www.gnu.org/licenses/>.
--------------------------------------------------------------------
Author: Csaba Pakozdi
--------------------------------------------------------------------*/

#include"wave_lib_ssgw.h"
#include <iomanip>

void wave_lib_ssgw::setWave(double k, double d,double H)
{
    defineInput2Iteration(k,d,H);
    PetviashviliIteration();
    getPhysicsParameters();
    surfaceCalculated = computeSurfaceVariables();
}

void wave_lib_ssgw::defineInput2Iteration(double kk, double dd, double HH)
{
    // Normalized water depth
    kd = kk*dd;
    // Wave steepness
    kH2 = 0.5*kk*HH;

    
    if((1.0-std::tanh(kd)) < tol) // Deep water case.
    {
        scaleVel = sqrt(GRAV/kk);
        scaleLength = 1/kk;
        // Depth
        d =  std::numeric_limits<double>::infinity();
        // Wavenumber
        k = 1.0;
        // Characteristic wavelength lambda
        lambda = 1/k;
        infiniteDepth = true;
    }
    else
    {
        scaleVel = sqrt(GRAV*dd);
        scaleLength = dd;
        // Depth
        d = 1.0;
        // Wavenumber
        k = kd/d;
        // Characteristic wavelength lambda
        lambda = std::tanh(kd)/k;
        infiniteDepth = false;
    }
    // Acceleration due to gravity
    g      = 1.0;
    // Linear phase velocity squared
    c02    = g*lambda;
    // Total wave height
    H      = 2.0*kH2/k;
    // Half-length of the computational domain (with \f$ c_r=c_e \f$)
    L      = M_PI/k;
    // \f$ \Delta \alpha\f$
    dalpha = L/N;
    // \f$ \Delta k\f$
    dk     = M_PI/L;
    
    // Vector of abscissas in the conformal space.
    vectorAlpha.setLinSpaced(2*N,0.0,static_cast<double>(2*N-1)*dalpha);
    
    // Vector of wavenumbers.
    vectorK.head(N) = Eigen::VectorXd::LinSpaced(N,0.0,static_cast<double>(N-1)*dk);
    vectorK.tail(N) = Eigen::VectorXd::LinSpaced(N,-static_cast<double>(N)*dk,-dk);
    
    // Iterations counter
    iter = 0;
    // Enforce loop entry
    err  = std::numeric_limits<double>::infinity();
    //fft.SetFlag(fft.HalfSpectrum);
}

void wave_lib_ssgw::PetviashviliIteration()
{
    // Iteration variables
    double mean_ys, sigma;
    double Bg=1.0;
    
    // Initial guess for the solution:
    
    // Airy solution for Upsilon.
    Upsilon = 0.5*H*(1.0+(k*vectorAlpha.array()).array().cos());
    
    // Parameter sigma. EQ. 5.6
    sigma = 1.0;
    
    while (err > tol)
    {
        // << Upsilon >>
        double mean_Upsilon = Upsilon.mean();
        // Y_s. Under Section 5.4 first
        Ys = Upsilon.array() - mean_Upsilon;
        fft.fwd(Ys_hat, Ys);
        
        // Compute sigma and delta
        
        // Deep water
        if(infiniteDepth)
        {
            sigma = 1.0; // sigma
            
            // C{ Y_s }
            dummy_hat =(vectorK.array().abs()*Ys_hat.array());
            fft.inv(dummy,dummy_hat);
            CYs = dummy.array().real();
            
            // % << y_s >>. EQ. 5.7
            mean_ys = (Ys.adjoint()*CYs);
            mean_ys = -0.5*mean_ys/static_cast<double>(N);
        } // if d
        // Finite depth.
        else
        {
            // Operator C in Fourier space. EQ. 3.17
            C_hat = vectorK.array()/(sigma*d*vectorK.array()).tanh();
            C_hat(0) = 1/(sigma*d);
            
            // Operator S^2 in Fourier space
            S2_hat = (vectorK.array()/(sigma*d*vectorK.array()).sinh()).pow(2);
            S2_hat(0) = 1/(sigma*sigma*d*d);
            
            // Equation for sigma. EQ. 5.8
            dummy_hat = C_hat.array()*Ys_hat.array();
            fft.inv(dummy,dummy_hat);
            double E = (Ys.array()*dummy.array().real()).mean() + (sigma-1.0)*d;

            // / Its derivative. EQ. 5.9
            dummy_hat = S2_hat.array()*Ys_hat.array();
            fft.inv(dummy,dummy_hat);
            double dE = d - d*(Ys.array()*dummy.array().real()).mean();
            
            // Newton new sigma. EQ. 5.9
            sigma -= E/dE;
            
            // << y_s >>. EQ. 5.8
            mean_ys = (sigma - 1.0)*d;
        } // else of if d
        
        // Parameter delta. Section 5.4 last
        double delta = mean_ys - mean_Upsilon;
        
        // Operator C in Fourier space. EQ. 3.17
        C_hat = vectorK.array()/(sigma*d*vectorK.array()).tanh();
        C_hat(0) = 1/(sigma*d);
        //std::cout << std::endl<< C_hat << std::endl;
        
        // Compute Bernoulli constant B.
        Upsilon2 = Upsilon.array().pow(2); // Upsilon^2
        //double mean_Upsilon2 = Upsilon2.mean(); // << Upsilon^2 >>
        
        // Fourier transform of Upsilon
        fft.fwd(Upsilon_hat, Upsilon);
        
        // Fourier transform of Upsilon^2
        fft.fwd(Upsilon2_hat, Upsilon2);
        
        // C{ Upsilon }.
        dummy_hat = C_hat.array()*Upsilon_hat.array();
        fft.inv(dummy,dummy_hat);
        CUpsilon = dummy.array().real();
        
        //C{ Upsilon^2 }
        fft.fwd(dummy, Upsilon2);
        dummy_hat = C_hat.array()*dummy.array();
        fft.inv(dummy,dummy_hat);
        CUpsilon2 = dummy.array().real();
        
        //C{ Upsilon }_trough - C{ Upsilon }_crest.
        double DCU = CUpsilon(N) - CUpsilon(0);
        
        //C{ Upsilon^2 }_trough - C{ Upsilon^2 }_crest
        double DCU2 = CUpsilon2(N) - CUpsilon2(0);
        
        // B/g. EQ. 4.7
        Bg    = 2.0*delta - H/sigma*(1.0+delta/d+sigma*CUpsilon(0))/DCU + 0.5*DCU2/DCU;
        
        //Define linear operators in Fourier space.
        //Operator C_inf.  EQ.3.17
        Cinf_hat = vectorK.array().abs();
        Cinf_hat(0) = 0.0;
        
        //Operator C_inf o C^{-1}. Text under EQ. 3.17
        CIC_hat  = ((sigma*d)*vectorK.array().abs()).tanh();
        // Regularisation.
        if(infiniteDepth)
        {
            CIC_hat(0) = 1.0;
        }
        // Operator L. EQ. 5.2
        L_hat = (Bg-2.0*delta)*Cinf_hat.array() - ((1.0+delta/d)/sigma)*CIC_hat.array();
        
        // Operator L^-1
        IL_hat = 1.0/L_hat.array();
        IL_hat(0) = 1.0;
        
        // Petviashvili's iteration.
        CUps_hat = C_hat.array()*Upsilon_hat.array();
        
        // L{Upsilon}. EQ. 5.2
        dummy_hat = L_hat.array()*Upsilon_hat.array();
        fft.inv(dummy,dummy_hat);
        LUps=dummy.array().real();
        
        // Nonlinear term in Fourier space. EQ. 5.3
        fft.inv(dummy,CUps_hat);
        dummy_real = Upsilon.array()*dummy.array().real();
        fft.fwd(Nupsilon_hat,dummy_real);
        Nupsilon_hat = CIC_hat.array()*Nupsilon_hat.array() + 0.5*Cinf_hat.array()*Upsilon2_hat.array();
        
        // N{ Upsilon }. EQ.5.3
        fft.inv(dummy,Nupsilon_hat);
        Nupsilon=dummy.array().real();
        
        //Weight. EQ.5.4 second part
        double Snom   = (Upsilon.adjoint()*LUps);
        double Sdenom = (Upsilon.adjoint()*Nupsilon);
        double S = Snom/Sdenom;
        // New Upsilon. EQ.5.4 first part
        dummy_hat = Nupsilon_hat.array()*IL_hat.array();
        fft.inv(dummy,dummy_hat);
        U = S*S*dummy.array().real();
        
        // Enforce mean value. EQ. 5.5
        U = H*(U.array() - U(N)) / (U(0)-U(N));
        
        // Update values.
        err = (U-Upsilon).lpNorm<Eigen::Infinity>();
        Upsilon = U;
        iter++;
    } // while loop
    //std::cout << "err= " << err << " iter= " << iter << std::endl;
    // Post processing
    
    // Inverse Hilbert transform. EQ.3.11
    std::complex<double> i1(0.0,1.0);
    IH_hat = -i1/(sigma*d*vectorK.array()).tanh();
    IH_hat(0) = 0;
    //  Under Section 5.4 first
    Ys = Upsilon.array() - Upsilon.array().mean();
    fft.fwd(Ys_hat, Ys);
    
    dummy_hat = C_hat.array()*Ys_hat.array();
    fft.inv(dummy,dummy_hat);
    CYs = dummy.array().real();
    
    // EQ. 3.11
    dummy_hat = IH_hat.array()*Ys_hat.array();
    fft.inv(dummy,dummy_hat);
    Xs = dummy.array().real();

    // << y_s >>. EQ. 5.7
    mean_ys = (Ys.adjoint()*CYs);
    mean_ys = -0.5*mean_ys/static_cast<double>(N);
    
    // Section 2, p. 8, first parag
    Zs = Xs.array() + i1*Ys.array();

    // tilde(z)(alpha) Appendix D. under EQ: D.2. Need for W(z) and F(z)
    zs = vectorAlpha.array() + i1*mean_ys + Zs.array();
    
    // FFT Zs Eigen accept only rael time series
    dummy_real = Zs.array().real();
    fft.fwd(dummy,dummy_real);
    dummy_real = Zs.array().imag();
    fft.fwd(dummy_hat,dummy_real);
    dummy = dummy.array()+i1*dummy_hat.array();
    
    dummy_hat = i1*vectorK.array()*dummy.array();
    fft.inv(dZs,dummy_hat);
    
    // Appendix D. under EQ.D.2. Need for W(z) and F(z)
    dzs = 1.0 + dZs.array();
    B = g*Bg;
    ce = ((1.0 + CYs.array())/(dzs.array().abs().pow(2))).sum()/2.0/static_cast<double>(N);
    ce = std::sqrt(B/ce);
    cs = sigma*ce;
    ws = -ce/dzs.array();
    a  = zs.imag().maxCoeff();
    b  = -zs.imag().minCoeff();
    
    double Bce2d;
    if (infiniteDepth)
    {
        Bce2d = 0.0;
        // Inverse C-operator.
        IC = 1/vectorK.array().abs();
        IC(0) = 0.0;
    }
    else
    {
        Bce2d = (B-ce*ce)*d;
        // Inverse C-operator.
        IC = (sigma*d*vectorK.array()).tanh()/vectorK.array();
        IC(0) = sigma*d;
    }
    
    // Appendix B
    ydx = dzs.array().real()*zs.array().imag();
    // Impulse
    intImpulse = -ce*mean_ys;
    // Potential energy
    intPotentialEnergy = 0.5*g*(ydx.array()*zs.array().imag()).mean();
    // Kinetic energy
    intKineticEnergy = 0.5*intImpulse*ce;
    // Radiation stress.
    intRadiationStress = 2.0*ce*intImpulse - 2.0*intPotentialEnergy + Bce2d;
    // Momentum flux
    intMomentumFlux = intRadiationStress - intPotentialEnergy + 0.5*g*d*d;
    // Energy flux
    intEnergyFlux = 0.5*Bce2d*ce + 0.5*(B+ce*ce)*intImpulse + (intKineticEnergy-2.0*intPotentialEnergy)*ce;
    // Group velocity.
    cg = intEnergyFlux/(intKineticEnergy+intPotentialEnergy);
    
    //Constant K1 EQ 4.3
    double K1 = zs.array().imag().matrix().adjoint()*(0.5*zs.array().imag()-Bg).matrix();
    K1 = 0.5*K1/static_cast<double>(N);
    
    // Residual
    fft.fwd(dummy,ydx);
    dummy_hat = dummy.array()*IC.array();
    fft.inv(dummy,dummy_hat);
    ICydx = dummy.array().real();
    errfun = (zs.array().imag() - Bg + (Bg*Bg+2.0*K1-2.0*ICydx.array()).sqrt()).matrix().lpNorm<Eigen::Infinity>();
}

bool wave_lib_ssgw::computeSurfaceVariables()
{
    
    Eigen::Map<Eigen::VectorXd> Xs(&xs[0],2*N);
    Xs.head(N) = scaleLength*zs.tail(N).array().real() - ParameterValue.waveLength;
    Xs.tail(N) = scaleLength*zs.head(N).array().real();
    
    Eigen::Map<Eigen::VectorXd> Ys(&ys[0],2*N);
    Ys.head(N) = scaleLength*zs.tail(N).array().imag();
    Ys.tail(N) = scaleLength*zs.head(N).array().imag();
        
    Eigen::Map<Eigen::VectorXd> Us(&us[0],2*N);
    Us.head(N) = scaleVel*(ws.tail(N).array().real() + ce);
    Us.tail(N) = scaleVel*(ws.head(N).array().real() + ce);
    
    Eigen::Map<Eigen::VectorXd> Vs(&vs[0],2*N);
    Vs.head(N) = -scaleVel*ws.tail(N).array().imag();
    Vs.tail(N) = -scaleVel*ws.head(N).array().imag();
    
    Eigen::Map<Eigen::VectorXd> Phis(&phis[0],2*N);
    dummy_real = ce*scaleVel*scaleLength*(zs.real().array()-vectorAlpha.array());
    Phis.head(N) = dummy_real.tail(N);
    Phis.tail(N) = dummy_real.head(N);
    
    return true;
}

void wave_lib_ssgw::writeResult(const std::string folderName)
{
    if (!surfaceCalculated)
    {
        surfaceCalculated = computeSurfaceVariables();
    }

    std::string fileName = folderName + "/ssgw.csv";
    std::ofstream myFile(fileName);
    Eigen::Map<Eigen::VectorXd> Ys(&ys[0],2*N);

    fft.fwd(dummy_hat,Ys);
    dummy_real = dummy_hat.array().abs();
    myFile << std::setw(16) << "xs," << std::setw(16) << "eta," << std::setw(16) << "fi," << std::setw(7) << "i," << std::setw(15) << "S(k)" << std::endl;
    for (int i = 0;i<xs.size();i++)
    {
        myFile << std::scientific;
        myFile <<std::setw(15) << xs[i] << "," << std::setw(15) << ys[i] << "," << std::setw(15) << phis[i] << "," << std::setw(6) << i << "," << std::setw(15) << dummy_real(i) << std::endl;
    }
    myFile.close();
}

void wave_lib_ssgw::computePotentialField(std::vector<double>& x, std::vector<double>& y, std::vector<double>& phi)
{
    if (!surfaceCalculated)
    {
        getPhysicsParameters();
        surfaceCalculated = computeSurfaceVariables();
    }
    
    double dal = 0.5*ParameterValue.waveLength/static_cast<double>(N);
    std::complex<double> i1(0.0,1.0);
    std::complex<double> integralScaler = 0.5*i1*ParameterValue.phaseVelocity*dal/M_PI;
    std::complex<double> kd2i = 2.0*i1*ParameterValue.waterDepth*ParameterValue.waveNumber;
    std::complex<double> kdi  = i1*ParameterValue.waterDepth*ParameterValue.waveNumber;

    Eigen::VectorXcd kzs       = ParameterValue.waveNumber*zs.array()*scaleLength;
    Eigen::VectorXcd kzsconj   = kzs.array().conjugate();
    Eigen::RowVectorXcd BB     = dzs.array()-1.0;
    Eigen::RowVectorXcd BBconj = BB.array().conjugate();
    Eigen::VectorXcd DDnominator = (0.5*(kzs.array()+kdi)).sin();
    Eigen::VectorXcd EEnominator = (0.5*(kzsconj.array()-kdi)).sin();
    
    for (int j=0;j<phi.size();j++)
    {
        std::complex<double> kz   = ParameterValue.waveNumber*(x[j] + i1*y[j]);
        Eigen::VectorXcd DD = (DDnominator.array()/(0.5*(kzs.array()-kz)).sin()).log();
        Eigen::VectorXcd EE_ = (EEnominator.array()/(0.5*(kzsconj.array()-kd2i-kz)).sin()).log();
                
        std::complex<double> f  = BB*DD;
        std::complex<double> fm = BBconj*EE_;
        std::complex<double> F = f - fm;// Because of periodity this sum is trapez integration !!!!!
        F *= integralScaler;
        phi[j] = F.real();
    }
}

void wave_lib_ssgw::computeVelocityField(std::vector<double>& x, std::vector<double>& y, std::vector<double>& u, std::vector<double>& v)
{
    if (!surfaceCalculated)
    {
        getPhysicsParameters();
        surfaceCalculated = computeSurfaceVariables();
    }
    
    double dal = 0.5*ParameterValue.waveLength/static_cast<double>(N);
    std::complex<double> i1(0.0,1.0);
    std::complex<double> integralScaler = 0.25*i1*ParameterValue.waveNumber*ParameterValue.phaseVelocity*dal/M_PI;
    std::complex<double> kd2i = 2.0*i1*ParameterValue.waterDepth*ParameterValue.waveNumber;
    
    Eigen::VectorXcd kzs       = ParameterValue.waveNumber*zs.array()*scaleLength;
    Eigen::VectorXcd kzsconj   = kzs.array().conjugate();
    Eigen::RowVectorXcd BB     = dzs.array()-1.0;
    Eigen::RowVectorXcd BBconj = BB.array().conjugate();
    
    for (int j=0;j<u.size();j++)
    {
        std::complex<double> kz = ParameterValue.waveNumber*(x[j] + i1*y[j]);
        Eigen::VectorXcd AA = 1.0/(0.5*(kzs.array()-kz)).tan();
        Eigen::VectorXcd CC = 1.0/(0.5*(kzsconj.array()-kd2i-kz)).tan();
                
        std::complex<double> w  = BB*AA;
        std::complex<double> wm = BBconj*CC;
        std::complex<double> W = w - wm;// Because of periodity this sum is trapez integration !!!!!
        W *= integralScaler;
        u[j] = W.real();
        v[j] = -W.imag();
    }
}

void wave_lib_ssgw::getPhysicsParameters()
{
    ParameterValue.N   = N;
    ParameterValue.tol = tol;
    
    if(infiniteDepth)
    {
        ParameterValue.waterDepth   = kd*scaleLength;
        ParameterValue.waveNumber   = k/scaleLength;
    }
    else
    {
        ParameterValue.waterDepth   = scaleLength;
        ParameterValue.waveNumber   = kd/scaleLength;
    }
    
    ParameterValue.waveHeight   = H*scaleLength;
    ParameterValue.waveLength   = 2.0*L*scaleLength;
    ParameterValue.gravityAcc   = g*scaleVel*scaleVel/scaleLength;
    
    ParameterValue.phaseVelocity     = ce*scaleVel;
    ParameterValue.celerityCs        = cs*scaleVel;
    ParameterValue.creastHeight      = a*scaleLength;
    ParameterValue.troughHeight      = b*scaleLength;
    ParameterValue.groupVelocity     = cg*scaleVel;
    ParameterValue.BernoulliConstant = B*scaleVel*scaleVel;
    
    ParameterValue.Impulse         = intImpulse*scaleVel*scaleLength;
    ParameterValue.PotentialEnergy = intPotentialEnergy*scaleVel*scaleVel*scaleLength;
    ParameterValue.KineticEnergy   = intKineticEnergy*scaleVel*scaleVel*scaleLength;
    ParameterValue.RadiationStress = intRadiationStress*scaleVel*scaleVel*scaleLength;
    ParameterValue.MomentumFlux    = intMomentumFlux*scaleVel*scaleVel*scaleLength;
    ParameterValue.EnergyFlux      = intEnergyFlux*scaleVel*scaleVel*scaleVel*scaleLength;
}

bool wave_lib_ssgw::resizing()
{
    //complex
    zs.resize(2*N);
    dzs.resize(2*N);
    //Real
    vectorAlpha.resize(2*N);
    vectorK.resize(2*N);
    
    // Variable for the iteration
    //Real
    dummy_real.resize(2*N);
    Upsilon.resize(2*N);
    Ys.resize(2*N);
    CYs.resize(2*N);
    C_hat.resize(2*N);
    S2_hat.resize(2*N);
    Upsilon2.resize(2*N);
    CUpsilon.resize(2*N);
    CUpsilon2.resize(2*N);
    Cinf_hat.resize(2*N);
    CIC_hat.resize(2*N);
    L_hat.resize(2*N);
    
    // Variable for the post processing
    IL_hat.resize(2*N);
    LUps.resize(2*N);
    Nupsilon.resize(2*N);
    IC.resize(2*N);
    ydx.resize(2*N);
    ICydx.resize(2*N);
    
    U.resize(2*N);
    Xs.resize(2*N);
    phis.resize(2*N);
    xs.resize(2*N);
    ys.resize(2*N);
    us.resize(2*N);
    vs.resize(2*N);

    // Complex
    dummy_hat.resize(2*N);
    dummy.resize(2*N);
    Ys_hat.resize(2*N);
    Upsilon_hat.resize(2*N);
    Upsilon2_hat.resize(2*N);
    CUps_hat.resize(2*N);
    Nupsilon_hat.resize(2*N);
    IH_hat.resize(2*N);
    Zs.resize(2*N);
    dZs.resize(2*N);
    ws.resize(2*N);
    return true;
}

double wave_lib_ssgw::modulo(double a, double b)
{
    double m = fmod(a, b);
    return m + b * (m < 0.f);
}

