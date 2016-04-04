#ifndef FD_Sn_Transport_hh
#define FD_Sn_Transport_hh

#include <string>
#include <vector>

using std::string;
using std::vector;

/*
  Create and solve finite difference Sn equations for
  one energy group and isotropic scattering
  
  See FD_Sn_Transport.cc for full documentation
*/

class FD_Sn_Transport
{
private:
    string filename_in_;
    string filename_out_;
    
    int max_iterations_;
    int tolerance_;
    int p_norm_;
    int number_of_boundaries_;
    int number_of_cells_;
    int number_of_edges_;
    int number_of_ordinates_;
    int number_of_regions_;

    vector<double> ordinates_;
    vector<double> weights_;
    vector<double> boundary_source_;
    vector<double> alpha_;
    vector<double> cell_length_;
    vector<double> sigma_t_;
    vector<double> sigma_s_;
    vector<double> internal_source_;

    bool converged_;
    int total_iterations_;
    double time_;
    vector<double> spectral_radius_;
    
    vector<double> phi_;
    vector<double> phi_old_;
    vector<double> psi_average_;
    vector<double> psi_edge_;
    vector<double> source_;
    
    void parse_xml();
    void write_xml();
    void check_class_invariants();
    void update_flux();
    void update_source();
    void calculate_spectral_radius();
    bool check_convergence();
    double lp_norm(vector<double> phi);
    void sweep();
    void solve();
    
public:

    FD_Sn_Transport(string filename_in);
};

#endif
