#include "FD_Sn_Transport.hh"

#include <cmath>
#include <iostream>
#include <iterator>
#include <sstream>
#include <string>
#include <vector>

#include "Check.hh"
#include "Gauss_Legendre.hh"
#include "XML_Child_Value.hh"
#include "Timer.hh"

using namespace std;

/*
  Create and solve finite difference Sn equations for
  one energy group and isotropic scattering
*/
FD_Sn_Transport::
FD_Sn_Transport(string filename_in):
    filename_in_(filename_in)
{
    // Parse problem
    parse_xml();

    // Check array sizes
    check_class_invariants();

    // Solve problem
    solve();

    // Make sure array sizes haven't changed
    check_class_invariants();

    // Write output data
    write_xml();
}

/*
  Check array sizes and other class invariants
*/
void FD_Sn_Transport::
check_class_invariants()
{
    Assert(number_of_edges_ == number_of_cells_ + 1);
    Assert(ordinates_.size() == number_of_ordinates_);
    Assert(weights_.size() == number_of_ordinates_);
    Assert(boundary_source_.size() == number_of_boundaries_);
    Assert(alpha_.size() == number_of_boundaries_);
    Assert(cell_length_.size() == number_of_cells_);
    Assert(sigma_t_.size() == number_of_cells_);
    Assert(internal_source_.size() == number_of_cells_);
    Assert(spectral_radius_.size() == total_iterations_);
    Assert(phi_.size() == number_of_cells_);
    Assert(psi_average_.size() == number_of_cells_ * number_of_ordinates_);
    Assert(psi_edge_.size() == number_of_edges_ * number_of_ordinates_);
    Assert(source_.size() == number_of_cells_);
}

/*
  Parse problem from xml input file
*/
void FD_Sn_Transport::
parse_xml()
{
    // The output filename will be filename-out.xml
    
    filename_out_ = filename_in_.substr(0, filename_in_.find_last_of(".")) + "-out.xml";

    // Open input document and parse
    
    pugi::xml_document doc;
    
    if (!doc.load_file(filename_in_.c_str()))
    {
        AssertMsg(false, "could not open document");
    }

    // Create additional xml nodes
    
    pugi::xml_node input = doc.child("input");
    pugi::xml_node problem = input.child("problem");
    pugi::xml_node boundary = input.child("boundary");
    pugi::xml_node regions = input.child("regions");

    // Parse global problem data
    
    max_iterations_ = child_value<int>(problem, "max_iterations");
    tolerance_ = child_value<double>(problem, "tolerance");
    delta_tolerance_ = tolerance_ * tolerance_;
    p_norm_ = 2;
    number_of_boundaries_ = 2;
    number_of_cells_ = child_value<int>(problem, "number_of_cells");
    number_of_edges_ = number_of_cells_ + 1;
    number_of_ordinates_ = child_value<int>(problem, "number_of_ordinates");
    number_of_regions_ = child_value<int>(regions, "number_of_regions");

    // Create ordinate system

    ordinates_.resize(number_of_ordinates_);
    weights_.resize(number_of_ordinates_);
    
    gauss_legendre_vec(number_of_ordinates_, ordinates_, weights_);
    
    // Parse boundary data
    
    boundary_source_ = child_vector<double>(boundary, "source", number_of_boundaries_);
    alpha_ = child_vector<double>(boundary, "alpha", number_of_boundaries_);
    
    // Parse region (cellwise) data
    
    int cell_sum = 0;
    int region_sum = 0;
    for (pugi::xml_node region = regions.child("region"); region; region = region.next_sibling("region"))
    {
        int region_cells = child_value<int>(region, "number_of_cells");
        double cell_length = child_value<double>(region, "cell_length");
        double sigma_t = child_value<double>(region, "sigma_t");
        double sigma_s = child_value<double>(region, "sigma_s");
        double internal_source = child_value<double>(region, "internal_source");
        
        for (int i = cell_sum; i < cell_sum + region_cells; ++i)
        {
            cell_length_.push_back(cell_length);
            sigma_t_.push_back(sigma_t);
            sigma_s_.push_back(sigma_s);
            internal_source_.push_back(internal_source);
        }

        cell_sum += region_cells;
        region_sum += 1;
    }
    Check(cell_sum == number_of_cells_);
    Check(region_sum == number_of_regions_);

    // Set initial values of solution data

    total_iterations_ = 0;
    time_ = 0;
    spectral_radius_.resize(0);
    
    phi_.assign(number_of_cells_, 0);
    phi_old_.assign(number_of_cells_, 0);
    phi_older_.assign(number_of_cells_, 0);
    psi_average_.assign(number_of_cells_ * number_of_ordinates_, 0);
    psi_edge_.assign(number_of_edges_ * number_of_ordinates_ , 0);
    source_.assign(number_of_cells_, 0);
}

/*
  Output data to xml file
*/
void FD_Sn_Transport::
write_xml()
{
    pugi::xml_document doc;

    pugi::xml_node output = doc.append_child("output");
    pugi::xml_node data = output.append_child("data");
    pugi::xml_node solution = output.append_child("solution");

    append_child(data, filename_in_, "input_filename");
    append_child(data, filename_out_, "output_filename");
    
    append_child(data, number_of_boundaries_, "number_of_boundaries");
    append_child(data, number_of_cells_, "number_of_cells");
    append_child(data, number_of_edges_, "number_of_edges");
    append_child(data, number_of_ordinates_, "number_of_ordinates");
    append_child(data, number_of_regions_, "number_of_regions");

    append_child(data, ordinates_, "ordinates");
    append_child(data, weights_, "weights");
    append_child(data, boundary_source_, "boundary_source");
    append_child(data, alpha_, "alpha");
    append_child(data, cell_length_, "cell_length");
    append_child(data, sigma_t_, "sigma_t");
    append_child(data, sigma_s_, "sigma_s");
    append_child(data, internal_source_, "internal_source");
    
    append_child(data, p_norm_, "p_norm");
    append_child(data, max_iterations_, "max_iterations");
    append_child(data, tolerance_, "tolerance");
    append_child(data, delta_tolerance_, "delta_tolerance");
    append_child(solution, converged_, "converged");
    append_child(solution, total_iterations_, "total_iterations");
    append_child(solution, time_, "time");
    append_child(solution, spectral_radius_, "spectral_radius");
    
    append_child(solution, phi_, "phi");
    append_child(solution, psi_average_, "psi", "ordinate-cell");
    append_child(solution, psi_edge_, "psi_edge", "ordinate-edge");
    append_child(solution, source_, "source");
    
    doc.save_file(filename_out_.c_str());
}

/* 
   Update scalar flux to most recent value
*/
void FD_Sn_Transport::
update_flux()
{
    phi_older_ = phi_old_;
    phi_old_ = phi_;

    for (int i = 0; i < number_of_cells_; ++i)
    {
        double sum = 0;
        
        for (unsigned o = 0; o < number_of_ordinates_; ++o)
        {
            int k = o + number_of_ordinates_ * i;

            sum += weights_[o] * psi_average_[k];
        }

        phi_[i] = sum;
    }
}

/*
  Update source
*/
void FD_Sn_Transport::
update_source()
{
    for (int i = 0; i < number_of_cells_; ++i)
    {
        source_[i] = 0.5 * sigma_s_[i] * phi_[i] + 0.5 * internal_source_[i];
    }
}

/*
  Calculate spectral radius
*/
void FD_Sn_Transport::
calculate_spectral_radius()
{
    double num = 0; // numerator
    double den = 0; // denomintor
    for (int i = 0; i < number_of_cells_; ++i)
    {
        num += pow(phi_[i] - phi_old_[i], 2) * cell_length_[i];
        den += pow(phi_old_[i] - phi_older_[i], 2) * cell_length_[i];
    }
    num = sqrt(num);
    den = sqrt(den);
    
    spectral_radius_.push_back(num / den);
}

/*
  Check whether relative error is less than tolerance
*/
bool FD_Sn_Transport::
check_convergence()
{
    converged_ = false;
    
    for (unsigned i = 0; i < number_of_cells_; ++i)
    {
        if (abs(phi_[i] - phi_old_[i]) / (abs(phi_old_[i]) + delta_tolerance_) > tolerance_)
        {
            return false;
        }
    }

    converged_ = true;
    
    return true;
}

/* 
   Perform a single transport sweep
*/
void FD_Sn_Transport::
sweep()
{
    // Assign left boundary condition
    {
        int i = 0;

        for (int o = number_of_ordinates_ / 2; o < number_of_ordinates_; ++o)
        {
            int k1 = (number_of_ordinates_ - 1 - o) + number_of_ordinates_ * i;
            int k2 = o + number_of_ordinates_ * i;
            
            psi_edge_[k2] = boundary_source_[0] + psi_edge_[k1] * alpha_[0];
        }
    }

    // Sweep right over cells

    for (int i = 0; i < number_of_cells_; ++i)
    {
        for (int o = number_of_ordinates_ / 2; o < number_of_ordinates_; ++o)
        {
            int k1 = o + number_of_ordinates_ * i;
            int k2 = o + number_of_ordinates_ * (i + 1);

            double t1 = (ordinates_[o] - 0.5 * sigma_t_[i] * cell_length_[i]) * psi_edge_[k1];
            double t2 = source_[i] * cell_length_[i];
            double t3 = ordinates_[o] + 0.5 * sigma_t_[i] * cell_length_[i];
            
            psi_edge_[k2] = (t1 + t2) / t3;
            psi_average_[k1] = 0.5 * (psi_edge_[k1] + psi_edge_[k2]);
        }
    }
    
    // Sweep final cell
    {
        int i = number_of_cells_ - 1;

        for (int o = 0; o < number_of_ordinates_ / 2; ++o)
        {
            int k1 = (number_of_ordinates_ - 1 - o) + number_of_ordinates_ * (i + 1);
            int k2 = o + number_of_ordinates_ * (i + 1);
            
            psi_edge_[k2] = boundary_source_[1] + psi_edge_[k1] * alpha_[1];
        }
    }
    
    // Sweep left over cells
    for (int i = 0; i < number_of_cells_; ++i)
    {
        for (int o = 0; o < number_of_ordinates_ / 2; ++o)
        {
            int k1 = o + number_of_ordinates_ * i;
            int k2 = o + number_of_ordinates_ * (i + 1);

            double t1 = (-ordinates_[o] - 0.5 * sigma_t_[i] * cell_length_[i]) * psi_edge_[k2];
            double t2 = source_[i] * cell_length_[i];
            double t3 = -ordinates_[o] + 0.5 * sigma_t_[i] * cell_length_[i];
            
            psi_edge_[k1] = (t1 + t2) / t3;
            psi_average_[k1] = 0.5 * (psi_edge_[k1] + psi_edge_[k2]);
        }
    }
}

/* 
   Perform transport sweeps until convergence
*/
void FD_Sn_Transport::
solve()
{
    // Ensure iteration data is initialized correctly
    total_iterations_ = 0;
    time_ = 0;
    spectral_radius_.resize(0);

    // Get flux and source for first iteration
    update_flux();
    update_source();

    // Start iteration timer
    Timer timer;
    timer.start();
    
    // Perform source iteration
    for (int i = 0; i < max_iterations_; ++i)
    {
        // Calculate new angular flux
        sweep();

        // Calculate new scalar flux
        update_flux();

        // Calculate new source
        update_source();
        
        // Calculate spectral radius
        total_iterations_ += 1;
        calculate_spectral_radius();

        // Check for convergence
        if(check_convergence())
        {
            break;
        }
    }

    timer.stop();
    time_ = timer.time();
}
