#include "specie.h"

// Default constructor
Specie::Specie()
    : id(0), name(""), charge(0), mass(0.0), qm_ratio(0.0), react_net({}) {}

// Parameterized constructor
Specie::Specie(int _id, const std::string& _name, int _charge, double _mass,
               const std::vector<int>& _react_net, const Eigen::MatrixXd& _n)
    : id(_id), name(_name), charge(_charge), mass(_mass),
      react_net(_react_net), n(_n) {
    qm_ratio = charge / mass * E_CHARGE;
}

// Getters
int Specie::get_id() const { return id; }
std::string Specie::get_name() const { return name; }
int Specie::get_charge() const { return charge; }
double Specie::get_mass() const { return mass; }
double Specie::get_qm_ratio() const { return qm_ratio; }
const std::vector<int>& Specie::get_react_net() const { return react_net; }
Eigen::MatrixXd Specie::get_density() const { return n; }
Eigen::MatrixXd Specie::get_mob_coef() const { return mu; }
Eigen::MatrixXd Specie::get_dif_coef() const { return De; }

// Updaters with validation
void Specie::update_density(const Eigen::MatrixXd& new_n) {
    if (new_n.rows() != n.rows() || new_n.cols() != n.cols()) {
        throw std::invalid_argument("Density matrix size mismatch");
    }
    n = new_n;
}

void Specie::update_mob_coef(const Eigen::MatrixXd& new_mu) {
    if (new_mu.rows() != n.rows() || new_mu.cols() != n.cols()) {
        throw std::invalid_argument("Mobility matrix size mismatch");
    }
    mu = new_mu;
}

void Specie::update_dif_coef(const Eigen::MatrixXd& new_De) {
    if (new_De.rows() != n.rows() || new_De.cols() != n.cols()) {
        throw std::invalid_argument("Diffusion matrix size mismatch");
    }
    De = new_De;
}

// Snapshot
SpecieData Specie::snapshot() const {
    return SpecieData{id, name, charge, mass, qm_ratio, react_net, n, mu, De};
}

// Serialization
std::string Specie::to_json() const {
    std::ostringstream oss;
    oss << "{ \"id\": " << id
        << ", \"name\": \"" << name << "\""
        << ", \"charge\": " << charge
        << ", \"mass\": " << mass
        << ", \"qm_ratio\": " << qm_ratio
        << ", \"react_net_size\": " << react_net.size()
        << " }";
    return oss.str();
}

// Comparison
bool Specie::operator==(const Specie& other) const {
    return id == other.id &&
           name == other.name &&
           charge == other.charge &&
           mass == other.mass;
}
