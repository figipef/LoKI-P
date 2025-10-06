#pragma once
#include <vector>
#include <string>
#include <stdexcept>
#include <iostream>

class Grid {
public:
    Grid(int grid_size,
         int plasma_start,
         int plasma_end,
         const std::vector<int>& region_cells,
         const std::vector<double>& region_lengths);

    // --- Accessors ---
    int size() const { return grid_size_; }

    int plasma_start_index() const { return plasma_start_; }
    int plasma_end_index() const { return plasma_end_; }

    double total_length() const { return total_length_; }

    const std::vector<double>& cell_centers() const { return cell_centers_; }
    const std::vector<double>& cell_lengths() const { return cell_lengths_; }
    const std::vector<double>& boundaries() const { return boundaries_; }

    double x_center(int i) const { return cell_centers_.at(i); }
    double dx(int i) const { return cell_lengths_.at(i); }
    double x_boundary(int i) const { return boundaries_.at(i); }

    void print_summary() const;

private:
    // --- Grid configuration ---
    int grid_size_;
    int plasma_start_;
    int plasma_end_;

    std::vector<int> region_cells_;
    std::vector<double> region_lengths_;
    double total_length_;

    // --- Computed geometry ---
    std::vector<double> boundaries_;    // N+1 boundaries
    std::vector<double> cell_centers_;  // N centers
    std::vector<double> cell_lengths_;  // N cell widths

    void build_grid();
};
