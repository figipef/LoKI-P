#ifndef SPECIE_H
#define SPECIE_H

#include <string>
#include <vector>
#include <sstream>
#include <stdexcept>
#include <Eigen/Dense>

constexpr double E_CHARGE = 1.60217662e-19; // Elementary charge in Coulombs

struct SpecieData {
    int id;
    std::string name;
    int charge;
    double mass;
    double qm_ratio;
    std::vector<int> react_net;
    Eigen::MatrixXd n;
    Eigen::MatrixXd mu;
    Eigen::MatrixXd De;
};

class Specie {
private:
    int id;
    std::string name;
    int charge;
    double mass;
    double qm_ratio;
    std::vector<int> react_net;

    Eigen::MatrixXd n;
    Eigen::MatrixXd mu;
    Eigen::MatrixXd De;

public:
    // Constructors
    Specie();
    Specie(int _id, const std::string& _name, int _charge, double _mass,
           const std::vector<int>& _react_net, const Eigen::MatrixXd& _n);

    // Getters
    int get_id() const;
    std::string get_name() const;
    int get_charge() const;
    double get_mass() const;
    double get_qm_ratio() const;
    const std::vector<int>& get_react_net() const;
    Eigen::MatrixXd get_density() const;
    Eigen::MatrixXd get_mob_coef() const;
    Eigen::MatrixXd get_dif_coef() const;

    // Setters / Updaters
    void update_density(const Eigen::MatrixXd& new_n);
    void update_mob_coef(const Eigen::MatrixXd& new_mu);
    void update_dif_coef(const Eigen::MatrixXd& new_De);

    // Snapshot
    SpecieData snapshot() const;

    // Serialization
    std::string to_json() const;

    // Comparison
    bool operator==(const Specie& other) const;
};

#endif // SPECIE_H
