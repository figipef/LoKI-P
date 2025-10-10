#include "grid.h"
#include <iomanip>
#include <numeric>

Grid::Grid(int grid_size,
           int plasma_start,
           int plasma_end,
           const std::vector<int>& region_cells,
           const std::vector<double>& region_lengths)
    : grid_size_(grid_size),
      plasma_start_(plasma_start),
      plasma_end_(plasma_end),
      region_cells_(region_cells),
      region_lengths_(region_lengths)
{
    if (region_cells_.size() != region_lengths_.size())
        throw std::runtime_error("region_cells and region_lengths must have the same size.");

    int total_cells = std::accumulate(region_cells_.begin(), region_cells_.end(), 0);
    if (total_cells != grid_size_)
        throw std::runtime_error("Sum of region_cells must equal GRID_SIZE.");

    total_length_ = std::accumulate(region_lengths_.begin(), region_lengths_.end(), 0.0);

    if (plasma_start_ < 0 || plasma_end_ >= grid_size_ || plasma_start_ >= plasma_end_)
        throw std::runtime_error("Invalid plasma region indices.");

    build_grid();
}

void Grid::build_grid()
{
    boundaries_.resize(grid_size_ + 1);
    cell_centers_.resize(grid_size_);
    cell_lengths_.resize(grid_size_);

    double x = 0.0;
    int current_cell = 0;

    for (size_t region = 0; region < region_cells_.size(); ++region) {
        int n = region_cells_[region];
        double L = region_lengths_[region];
        double dx = L / n;

        for (int i = 0; i < n; ++i) {
            boundaries_[current_cell] = x;
            x += dx;
            cell_lengths_[current_cell] = dx;
            cell_centers_[current_cell] = x - 0.5 * dx;
            ++current_cell;
        }
    }

    boundaries_[grid_size_] = x; // last boundary
}

void Grid::print_summary() const
{
    std::cout << std::fixed << std::setprecision(6);
    std::cout << "Grid Summary:\n";
    std::cout << "  Total cells:   " << grid_size_ << "\n";
    std::cout << "  Total length:  " << total_length_ << " m\n";
    std::cout << "  Regions:       " << region_cells_.size() << "\n";

    for (size_t i = 0; i < region_cells_.size(); ++i) {
        std::cout << "    Region " << i+1
                  << ": " << region_cells_[i] << " cells, "
                  << region_lengths_[i] << " m, "
                  << "dx = " << region_lengths_[i] / region_cells_[i] << " m\n";
    }

    std::cout << "  Plasma region:\n";
    std::cout << "    Cell indices: [" << plasma_start_ << ", " << plasma_end_ << "]\n";
    std::cout << "    Cell centers: [" 
              << cell_centers_[plasma_start_] << " m, "
              << cell_centers_[plasma_end_] << " m]\n";
    std::cout << "  Boundaries stored: " << boundaries_.size()
              << " (N+1 = " << grid_size_ + 1 << ")\n\n";
}
