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

    // Input and output files
    string filename_in_;
    string filename_out_;

    // Problem size data
    int number_of_boundaries_;
    int number_of_cells_;
    int number_of_edges_;
    int number_of_ordinates_;
    int number_of_regions_;

    // Constant problem data
    vector<double> ordinates_;
    vector<double> weights_;
    vector<double> boundary_source_;
    vector<double> alpha_;
    vector<double> cell_length_;
    vector<double> sigma_t_;
    vector<double> sigma_s_;
    vector<double> internal_source_;

    // Convergence data
    int p_norm_;
    int max_iterations_;
    double tolerance_;
    bool converged_;
    int total_iterations_;
    double time_;
    vector<double> spectral_radius_;

    // Solution-dependent problem data
    vector<double> phi_;
    vector<double> phi_old_;
    vector<double> psi_average_;
    vector<double> psi_edge_;
    vector<double> source_;

    // Functions
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

    // Creator
    FD_Sn_Transport(string filename_in);
};

#endif
