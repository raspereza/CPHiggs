#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>
#include <math.h>
#include <map>
#include <vector>
#include <string> 

using namespace std;
namespace py = pybind11;

map<string, py::array_t<double>> fastmtt(unsigned int N,
					 py::array_t<double> pt_1_vec,
					 py::array_t<double> eta_1_vec,
					 py::array_t<double> phi_1_vec,
					 py::array_t<double> mass_1_vec,
					 py::array_t<int> decay_type_1_vec,
					 py::array_t<double> pt_2_vec,
					 py::array_t<double> eta_2_vec,
					 py::array_t<double> phi_2_vec,
					 py::array_t<double> mass_2_vec,
					 py::array_t<int> decay_type_2_vec,
					 py::array_t<double> met_pt_vec,
					 py::array_t<double> met_phi_vec,
					 py::array_t<double> metcov_xx_vec,
					 py::array_t<double> metcov_xy_vec,
					 py::array_t<double> metcov_yy_vec,
					 bool verbosity,
					 double delta,
					 double reg_order,
					 double mX,
					 double GammaX) {
  //  bool verbosity = true;
  //  double delta=1/1.15;
  //  double reg_order=6;
  //  double mX = 125.10; // resonance mass (default 125.10 GeV)
  //  double GammaX = 0.004; // resonance width (default 4 MeV)

  double m_ele = 0.51100e-3;
  double m_muon = 0.10566;
  double m_tau = 1.77685;
  double m_pion = 0.13957;

  // Casting of variables
  auto pt_1 = pt_1_vec.unchecked<1>();
  auto eta_1 = eta_1_vec.unchecked<1>();
  auto phi_1 = phi_1_vec.unchecked<1>();
  auto mass_1 = mass_1_vec.unchecked<1>();
  auto decay_type_1 = decay_type_1_vec.unchecked<1>();
  
  auto pt_2 = pt_2_vec.unchecked<1>();
  auto eta_2 = eta_2_vec.unchecked<1>();
  auto phi_2 = phi_2_vec.unchecked<1>();
  auto mass_2 = mass_2_vec.unchecked<1>();
  auto decay_type_2 = decay_type_2_vec.unchecked<1>();

  auto met_pt  = met_pt_vec.unchecked<1>();
  auto met_phi = met_phi_vec.unchecked<1>();
  auto metcov_xx = metcov_xx_vec.unchecked<1>();
  auto metcov_xy = metcov_xy_vec.unchecked<1>();
  auto metcov_yy = metcov_yy_vec.unchecked<1>();

  py::buffer_info buffer = pt_1_vec.request();
  auto x_1_vec = py::array_t<double>(buffer.size);
  auto x_2_vec = py::array_t<double>(buffer.size);
  auto x_1_BW_vec = py::array_t<double>(buffer.size);
  auto x_2_BW_vec = py::array_t<double>(buffer.size);
  auto mass_vec = py::array_t<double>(buffer.size);
  auto mass_BW_vec = py::array_t<double>(buffer.size);
  
  py::buffer_info x_1_buffer = x_1_vec.request();
  py::buffer_info x_2_buffer = x_2_vec.request();
  py::buffer_info x_1_BW_buffer = x_1_BW_vec.request();
  py::buffer_info x_2_BW_buffer = x_2_BW_vec.request();
  py::buffer_info mass_buffer = mass_vec.request();
  py::buffer_info mass_BW_buffer = mass_BW_vec.request();
  
  double * x_1 = static_cast<double *>(x_1_buffer.ptr);
  double * x_2 = static_cast<double *>(x_2_buffer.ptr);
  double * x_1_BW = static_cast<double *>(x_1_BW_buffer.ptr);
  double * x_2_BW = static_cast<double *>(x_2_BW_buffer.ptr);
  double * mass = static_cast<double *>(mass_buffer.ptr);
  double * mass_BW = static_cast<double *>(mass_BW_buffer.ptr); 
  
  unsigned int counter = 0;
  for (unsigned int i=0; i<N; ++i) {  

    
    // grab the correct masses based on tau decay type
    // tau decay_type: 0 ==> leptonic to electron, 
    // 1 ==> leptonic to muon, 
    // 2 ==> leptonic to hadronic
    double m1 = mass_1(i);
    if (decay_type_1(i)==0)
      m1 = m_ele;
    else if (decay_type_1(i)==1)
      m1 = m_muon;
    
    double m2 = mass_2(i);
    if (decay_type_2(i)==0)
      m2 = m_ele;
    else if (decay_type_2(i)==1)
      m2 = m_muon;
    
    // store visible masses
    double m_vis_1 = m1;
    double m_vis_2 = m2;
    
    // determine minimum and maximum possible masses
    double m_vis_min_1 = 0; double m_vis_max_1 = 0;
    double m_vis_min_2 = 0; double m_vis_max_2 = 0; 
    if (decay_type_1(i) == 0) { m_vis_min_1 = m_ele ; m_vis_max_1 = m_ele; }
    if (decay_type_1(i) == 1) { m_vis_min_1 = m_muon; m_vis_max_1 = m_muon;}
    if (decay_type_1(i) == 2) { m_vis_min_1 = m_pion; m_vis_max_1 = 1.5;   }
    if (decay_type_2(i) == 0) { m_vis_min_2 = m_ele ; m_vis_max_2 = m_ele; }
    if (decay_type_2(i) == 1) { m_vis_min_2 = m_muon; m_vis_max_2 = m_muon;}
    if (decay_type_2(i) == 2) { m_vis_min_2 = m_pion; m_vis_max_2 = 1.5;   }
    
    if (m_vis_1 < m_vis_min_1) { m_vis_1 = m_vis_min_1;}
    if (m_vis_1 > m_vis_max_1) { m_vis_1 = m_vis_max_1;}
    if (m_vis_2 < m_vis_min_2) { m_vis_2 = m_vis_min_2;}
    if (m_vis_2 > m_vis_max_2) { m_vis_2 = m_vis_max_2;}
    
    // 4-vectors of taus
    double px_1 = pt_1(i)*cos(phi_1(i));
    double py_1 = pt_1(i)*sin(phi_1(i));
    double pz_1 = pt_1(i)*sinh(eta_1(i));
    double E_1 = sqrt(px_1*px_1+py_1*py_1+pz_1*pz_1+m_vis_1*m_vis_1);

    double px_2 = pt_2(i)*cos(phi_2(i));
    double py_2 = pt_2(i)*sin(phi_2(i));
    double pz_2 = pt_2(i)*sinh(eta_2(i));
    double E_2 = sqrt(px_2*px_2+py_2*py_2+pz_2*pz_2+m_vis_2*m_vis_2);
    
    // avoiding lorentzvectors from ROOT and akward
    double px_vis = px_1 + px_2;
    double py_vis = py_1 + py_2;
    double pz_vis = pz_1 + pz_2;
    double E_vis  = E_1 + E_2;

    double met_x = met_pt(i)*cos(met_phi(i));
    double met_y = met_pt(i)*sin(met_phi(i));
      
    double p_vis = sqrt(px_vis*px_vis+
			py_vis*py_vis+
			pz_vis*pz_vis);
    double m_vis = sqrt(E_vis*E_vis-p_vis*p_vis);
  
    // correct initial visible masses
    if (decay_type_1(i) == 2 && m_vis_1 > 1.5)
      m_vis_1 = 0.3;
    if (decay_type_2(i) == 2 && m_vis_2 > 1.5)
      m_vis_2 = 0.3;
  
    // invert met covariance matrix, calculate determinant
    double metcovinv_xx = metcov_yy(i);
    double metcovinv_yy = metcov_xx(i);
    double metcovinv_xy = -metcov_xy(i);
    double metcovinv_yx = -metcov_xy(i);
    double metcovinv_det = (metcovinv_xx*metcovinv_yy -
			    metcovinv_yx*metcovinv_xy);
    if (metcovinv_det<1e-10) { 
      printf("Warning! Ill-conditioned MET covariance at event index");
    }
    
    // perform likelihood scan 
    // see http://cms.cern.ch/iCMS/jsp/openfile.jsp?tp=draft&files=AN2019_032_v3.pdf
    double met_const = 1/sqrt(metcovinv_det);
    double min_likelihood = 0.;
    double x1_opt = 0.5;
    double x2_opt = 0.5;
    double min_likelihood_BW = 0.;
    double x1_opt_BW = 0.02;
    double x2_opt_BW = 0.02;
    
    // scan over weights for each ditau four-vector
    bool initialise = true;
    for (unsigned int i1 = 1; i1 < 100; ++i1) {
      double x1 = 0.01*double(i1);
      for (unsigned int i2 = 1; i2 < 100; ++i2) {
	double x2 = 0.01*double(i2);
	double test_mass = m_vis/sqrt(x1*x2);
	double x1_min = min(double(1.0), (m_vis_1/m_tau)*(m_vis_1/m_tau));
	double x2_min = min(double(1.0), (m_vis_2/m_tau)*(m_vis_2/m_tau));

	double zeroFactor = 1.0;
	if ((x1 < x1_min) || (x2 < x2_min)) 
	  zeroFactor = 0.0;
	
	// test weighted four-vectors
	double nu_test_px = px_1/x1 + px_2/x2 - px_vis;
	double nu_test_py = py_1/x1 + py_2/x2 - py_vis;
	
	// calculate mass likelihood integral 
	double m_shift = test_mass * delta;

	if (m_shift < m_vis)
	  zeroFactor = 0.0;
	
	x1_min = min(double(1.0), (m_vis_1/m_tau)*(m_vis_1/m_tau));
	x2_min = max((m_vis_2/m_tau)*(m_vis_2/m_tau), 
		     (m_vis/m_shift)*(m_vis/m_shift));
	double x2_max = min(double(1.0), (m_vis/m_shift)*(m_vis/m_shift)/x1_min);

	if (x2_max < x2_min)
	  zeroFactor = 0.0;

	double J = 2*m_vis*m_vis * pow(m_shift, -reg_order);
	double I_x2 = log(x2_max) - log(x2_min);
	double I_tot = I_x2;
	if (decay_type_1(i) != 2) {
	  double I_m_nunu_1 = (m_vis/m_shift) * (m_vis/m_shift) * (1.0/x2_max - 1.0/x2_min);
	  I_tot += I_m_nunu_1;
	}
	if (decay_type_2(i) != 2) {
	  double I_m_nunu_2 = (m_vis/m_shift) * (m_vis/m_shift) * I_x2 - (x2_max - x2_min);
	  I_tot += I_m_nunu_2;
	}
	double mass_likelihood = 1.0E9 * J * I_tot;
	
	// calculate MET transfer function 
	double residual_x = met_x - nu_test_px;
	double residual_y = met_y - nu_test_py;
	double pull2 = (residual_x*(metcovinv_xx*residual_x + 
				    metcovinv_xy*residual_y) +
			residual_y*(metcovinv_yx*residual_x +
				    metcovinv_yy*residual_y));
	if (metcovinv_det<1e-10)
	  zeroFactor = 0.0;
	
	pull2 /= metcovinv_det;
	double met_transfer = met_const*exp(-0.5*pull2);
	double deltaM = test_mass*test_mass-mX*mX;
	double mG = mX*GammaX;
	double BreitWigner_likelihood = mG*mG/(deltaM*deltaM + mG*mG); 
	double likelihood = -met_transfer * mass_likelihood * zeroFactor;
	double likelihood_BW = likelihood*BreitWigner_likelihood;
	if (initialise) {
	  min_likelihood = likelihood;
	  min_likelihood_BW = likelihood_BW;
	  x1_opt = x1;
	  x2_opt = x2;
	  x1_opt_BW = x1;
	  x2_opt_BW = x2;
	  initialise = false;
	  //	  printf("init : x1 = %5.3f  x2 = %5.3f\n",x1,x2);
	}
	else {
	  if (likelihood < min_likelihood) {
	    min_likelihood = likelihood;
	    x1_opt = x1;
	    x2_opt = x2;
	    //	    printf("min : x1 = %5.3f  x2 = %5.3f\n",x1,x2);
	  }
	  if (likelihood_BW < min_likelihood_BW) {
	    min_likelihood_BW = likelihood_BW;
	    x1_opt_BW = x1;
	    x2_opt_BW = x2;
	  }
	}
      }
    }
    x_1[i] = x1_opt;
    x_2[i] = x2_opt;
    x_1_BW[i] = x1_opt_BW;
    x_2_BW[i] = x2_opt_BW;
    mass[i] = m_vis/sqrt(x1_opt*x2_opt);
    mass_BW[i] = m_vis/sqrt(x1_opt_BW*x2_opt_BW);

    if (counter%1000==0 && verbosity) {
      printf("Processed %1i events out of %1i\n",counter,N);
      printf("m(vis) = %5.1f   m(tautau) = %5.1f   m(cons) = %5.1f\n",m_vis,mass[i],mass_BW[i]);
    }
    counter++;
  }

  map<string, py::array_t<double> > results = {
    {"x1",x_1_vec},
    {"x2",x_2_vec},
    {"x1_cons",x_1_BW_vec},
    {"x2_cons",x_2_BW_vec},
    {"mass",mass_vec},
    {"mass_cons",mass_BW_vec}
  };
  return results;
    
}

PYBIND11_MODULE(fastmtt, m) {
  m.def("fastmtt", fastmtt, "FastMTT");
}
