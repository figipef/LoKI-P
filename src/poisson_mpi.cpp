#include "poisson_mpi.h"
#include <iostream>
#include <cmath>
#include <stdexcept>
#include <algorithm> // for std::min

// ----------------------------------------------------
// Constructor: initialize MPI info, local domain
PoissonMPI::PoissonMPI(const Grid& grid,
                       const std::vector<double>& eps,
                       bool is_axial)
    : grid_(grid), eps_(eps), is_axial_(is_axial)
{
    MPI_Comm_rank(MPI_COMM_WORLD, &rank_);
    MPI_Comm_size(MPI_COMM_WORLD, &size_);

    int N = grid_.size();
    int base = N / size_;
    int rem = N % size_;

    local_start_ = rank_ * base + std::min(rank_, rem);
    local_N_ = base + (rank_ < rem ? 1 : 0);
    local_end_ = local_start_ + local_N_;

    phiw_.resize(local_N_, 0.0);
    phie_.resize(local_N_, 0.0);
    phic_.resize(local_N_, 0.0);

    rhs_surface_charges.resize(local_N_, 0.0);
    rhs_boundary_potential.resize(local_N_, 0.0);

    boundary_conditions_.resize(2, 0);
    boundary_values_.resize(2, 0);
    surface_charges_.resize(grid_.size() + 1, 0.0);

    initialize_coefficients();

    if(rank_==0)
        std::cout << "Initialized PoissonMPI operator on " << size_ << " MPI ranks.\n";
}

// ----------------------------------------------------
// Initialize local coefficients
void PoissonMPI::initialize_coefficients()
{
    const auto& centers = grid_.cell_centers();
    const auto& boundaries = grid_.boundaries();
    int N = grid_.size();

    for (int local = 0; local < local_N_; ++local) {
        int i = local_start_ + local; // global index

        double eps_i = eps_[i];
        double eps_ip1 = (i + 1 < N) ? eps_[i + 1] : eps_i;

        double dx_i   = (i + 1 < N) ? (centers[i+1] - boundaries[i+1]) : 0.0;
        double dx_ip1 = (i + 1 < N) ? (boundaries[i+1] - centers[i]) : 0.0;

        double east_coef = (i + 1 < N) ? eps_i * eps_ip1 / (eps_ip1*dx_ip1 + eps_i*dx_i) : 0.0;

        if (!is_axial_) east_coef *= 2 * M_PI * boundaries[i+1];
        std::cout << "qpnfpnqpnfqp\n";
        double west_coef = (local == 0 && i > 0) ? -1* eps_i * eps_[i-1] / (eps_i * (boundaries[i] - centers[i-1]) + eps_[i-1] * (centers[i] - boundaries[i])) : 0;
        std::cout <<west_coef<<"\n";
        phiw_[local] = (local > 0 ) ?  phie_[local-1] : west_coef;
        phie_[local] = (-1) * east_coef;
        phic_[local] = -(phiw_[local] + phie_[local]);
    }

    // Adjust for ghost cells at boundaries
    if (local_start_ == 0) {
        double dx_ip1 = (centers[0] - boundaries[0]);

        double west_coef = eps_[0] / dx_ip1;
        phiw_[0] = (-1) *west_coef; //phie_[0];         // left ghost
        phic_[0] -= phiw_[0];
    }
    if (local_end_ == N) {

        double dx_ip1 = (boundaries[N] - centers[N-1]);
        std::cout <<(boundaries[local_N_] - centers[local_N_-1])<<" help\n";
        double east_coef = eps_[local_N_-1] / dx_ip1;
        phie_[local_N_-1] = (-1) *east_coef;//phiw_[local_N_-1]; // right ghost
        phic_[local_N_-1] -= phie_[local_N_-1];
    }

    // Debug print
    for(int local = 0; local < local_N_; ++local){
        int i = local_start_ + local;
        std::cout << "rank " << rank_ << " i_local=" << local 
                  << " phiw=" << phiw_[local] 
                  << " phic=" << phic_[local] 
                  << " phie=" << phie_[local] << "\n";
    }
}


// ----------------------------------------------------
// Update boundary conditions
void PoissonMPI::update_boundary_values(const std::vector<int>& boundary_conditions,
                                        const std::vector<double>& boundary_values)
{
    boundary_conditions_ = boundary_conditions;
    boundary_values_ = boundary_values;

    if (rank_ == 0) {
        if (boundary_conditions[0])
            phic_[0] += phiw_[0];
        else
            rhs_boundary_potential[0] = -boundary_values[0] * phiw_[0];
    }

    if (rank_ == size_ - 1) {
        if (boundary_conditions[1])
            phic_[local_N_-1] += phie_[local_N_-1];
        else
            rhs_boundary_potential[local_N_-1] = -boundary_values[1] * phie_[local_N_-1];
            std::cout <<"açlamkafpae"<< rhs_boundary_potential[local_N_-1]<<"\n";

    }
}

// ----------------------------------------------------
// Update surface charges
void PoissonMPI::update_surface_charges(const std::vector<int>& boundary_idxs,
                                        const std::vector<double>& surface_charges)
{
    const auto& centers = grid_.cell_centers();
    const auto& boundaries = grid_.boundaries();
    const auto& permitivity = grid_.permitivity();

    int N = boundary_idxs.size();
    for(int i=0;i<N;i++){
        surface_charges_[boundary_idxs[i]] = surface_charges[i];
        if(boundary_idxs[i]>=local_start_ && boundary_idxs[i]<local_end_)
            rhs_surface_charges[boundary_idxs[i]-local_start_] = 0.0;
    }

    for(int i=0;i<N;i++){
        if(boundary_idxs[i]>=local_start_ && boundary_idxs[i]<local_end_)
            calculate_rhs_charges(boundary_idxs[i], surface_charges[i], centers, boundaries, permitivity);
    }
}

// ----------------------------------------------------
// Calculate RHS charges
void PoissonMPI::calculate_rhs_charges(int boundary_idx, double surface_charge,
                                       const std::vector<double>& centers,
                                       const std::vector<double>& boundaries,
                                       const std::vector<double>& permitivity)
{
    int i = boundary_idx - local_start_;
    if(i<0 || i>=local_N_) return;

    double nume = centers[boundary_idx] - boundaries[boundary_idx];
    double denom = permitivity[boundary_idx-1]*(centers[boundary_idx]-boundaries[boundary_idx])
                 + permitivity[boundary_idx]*(boundaries[boundary_idx]-centers[boundary_idx-1]);

    rhs_surface_charges[i] += surface_charge * permitivity[boundary_idx]*(-1)*nume/denom;
    if(i>0) rhs_surface_charges[i-1] += surface_charge * permitivity[boundary_idx-1]*nume/denom;
}


// ----------------------------------------------------
// Solve (MPI-safe): gather global coefficients & RHS on rank 0,
// form classical tridiagonal a,b,c,d and solve there, scatter solution back.
std::vector<double> PoissonMPI::solve(const std::vector<double>& charge_density)
{
    int N = grid_.size();
    std::vector<double> local_rhs(local_N_, 0.0);
    const auto& boundaries = grid_.boundaries();

    // Build local RHS (include geometry factor consistently)
    for (int gi = local_start_; gi < local_end_; ++gi) {
        int li = gi - local_start_;
        double geom_factor = (!is_axial_) ?
            (boundaries[gi+1]*boundaries[gi+1] - boundaries[gi]*boundaries[gi]) * M_PI :
            (boundaries[gi+1] - boundaries[gi]);
        local_rhs[li] = geom_factor * charge_density[gi]
                      + rhs_surface_charges[li]
                      + rhs_boundary_potential[li];
    }

    // Prepare counts/displacements
    std::vector<int> counts(size_), displs(size_);
    int base = N / size_, rem = N % size_;
    for (int r = 0; r < size_; ++r) {
        counts[r] = base + (r < rem ? 1 : 0);
        displs[r] = r * base + std::min(r, rem);
    }

    // Gather full RHS on rank 0
    std::vector<double> full_rhs;
    if (rank_ == 0) full_rhs.assign(N, 0.0);
    MPI_Gatherv(local_rhs.data(), local_N_, MPI_DOUBLE,
                full_rhs.data(), counts.data(), displs.data(), MPI_DOUBLE,
                0, MPI_COMM_WORLD);

    // Gather global coefficient arrays on rank 0
    std::vector<double> phiw_global, phic_global, phie_global;
    if (rank_ == 0) {
        phiw_global.assign(N, 0.0);
        phic_global.assign(N, 0.0);
        phie_global.assign(N, 0.0);
    }

    MPI_Gatherv(phiw_.data(), local_N_, MPI_DOUBLE,
                phiw_global.data(), counts.data(), displs.data(), MPI_DOUBLE,
                0, MPI_COMM_WORLD);
    MPI_Gatherv(phic_.data(), local_N_, MPI_DOUBLE,
                phic_global.data(), counts.data(), displs.data(), MPI_DOUBLE,
                0, MPI_COMM_WORLD);
    MPI_Gatherv(phie_.data(), local_N_, MPI_DOUBLE,
                phie_global.data(), counts.data(), displs.data(), MPI_DOUBLE,
                0, MPI_COMM_WORLD);

    // Rank 0: form standard tridiagonal a,b,c,d and solve
    std::vector<double> full_solution;
    if (rank_ == 0) {
        for (int i= 0; i< N; i++){
            std::cout << phie_global[i]<<"\n";
        }
        // sanity
        if ((int)phiw_global.size() != N || (int)phic_global.size() != N || (int)phie_global.size() != N || (int)full_rhs.size() != N) {
            throw std::runtime_error("PoissonMPI::solve: incorrect global sizes on rank 0");
        }

        full_solution = solve_thomasalg(phiw_global,phic_global,phie_global,full_rhs);
    }

    // Scatter solution to ranks
    std::vector<double> local_solution(local_N_, 0.0);
    MPI_Scatterv(full_solution.data(), counts.data(), displs.data(), MPI_DOUBLE,
                 local_solution.data(), local_N_, MPI_DOUBLE,
                 0, MPI_COMM_WORLD);

    // Return full solution on rank0 (so main prints it), local slice on other ranks
    if (rank_ == 0) return full_solution;
    return local_solution;
}



std::vector<double> PoissonMPI::solve_thomasalg(
    const std::vector<double>& gphiw_,
    const std::vector<double>& gphic_,
    const std::vector<double>& gphie_,
    const std::vector<double>& grhs)
{
    const int N = gphic_.size();
    const double epsilon = 1e-50;

    // Trim diagonals to match tridiagonal structure
    std::vector<double> gphiw_trimmed(gphiw_.begin() + 1, gphiw_.end()); // remove first
    std::vector<double> gphie_trimmed(gphie_.begin(), gphie_.end() - 1); // remove last

    std::vector<double> c_prime(N - 1, 0.0);  // Modified upper diagonal
    std::vector<double> d_prime(N, 0.0);      // Modified RHS
    std::vector<double> solution(N, 0.0);     // Final solution

    // Forward sweep
    c_prime[0] = gphie_trimmed[0] / (gphic_[0] + epsilon);
    d_prime[0] = grhs[0] / (gphic_[0] + epsilon);

    for (int i = 1; i < N - 1; ++i) {
        double denom = gphic_[i] - gphiw_trimmed[i - 1] * c_prime[i - 1];
        denom = (std::abs(denom) < epsilon) ? epsilon : denom;
        c_prime[i] = gphie_trimmed[i] / denom;
    }

    for (int i = 1; i < N; ++i) {
        double denom = gphic_[i] - gphiw_trimmed[i - 1] * c_prime[i - 1];
        denom = (std::abs(denom) < epsilon) ? epsilon : denom;
        d_prime[i] = (grhs[i] - gphiw_trimmed[i - 1] * d_prime[i - 1]) / denom;
    }

    // Back substitution
    solution[N - 1] = d_prime[N - 1];
    for (int i = N - 2; i >= 0; --i) {
        solution[i] = d_prime[i] - c_prime[i] * solution[i + 1];
    }

    return solution;
}



// ----------------------------------------------------
// Print summary
void PoissonMPI::print_summary() const
{
    if(rank_==0){
        std::cout << "\n=== PoissonMPI Operator Summary ===\n";
        std::cout << "Grid size: " << grid_.size() << "\n";
        std::cout << "Total MPI ranks: " << size_ << "\n";
        std::cout << "Boundary values: Left=" << boundary_values_[0]
                  << ", Right=" << boundary_values_[1] << "\n";
    }
}
