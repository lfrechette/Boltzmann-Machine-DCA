/* model.h Created by Layne Frechette on Feb. 17, 2021
 *
 * Implements DCA model class.
 */

#ifndef MODEL_H
#define MODEL_H

//Headers
#include <armadillo>

//Class for model
class model {
  public:

    int N; //Length of sequences
    int q; //Number of possible states (residue types)
    double lambda; //l2 regularization constant
    bool symmetrize_on;
    double avg_ene; //Average energy estimated by most recent MC run
    int num_seqs;
    std::vector<std::string> seqs; //Sample of sequences obtained by most recent MC run
    int nwell; //Number of wells (1 or 2)
    double Tmix; //Mixing temperature (only used for nwell=2)
    arma::mat h1, h2; //Single-site propensities
    arma::cube J1, J2; //Coupling parameters
    arma::mat mom1, mom1_ene1, mom1_ene2; //First moment
    arma::cube mom2, mom2_ene1, mom2_ene2; //Second moment
    double mom1_err; //Std error in first moment via block averaging
    double mom2_err; //Std error in second moment via block averaging    
    std::string mc_init;

    //Constructors
    model();
    model(int N1, int q1, int nwell1, double lam, bool symon, double Tm, std::string mc_init1);
    //Destructor
    ~model();

    //Other Methods

    /** Get the energy of a sequence.
     * @param seq, Protein sequence. */
    double get_energy(std::vector<int> &seq);

    double get_delta_energy_single(std::vector<int> &seq, int well, int m, int r);
    double get_energy_single(std::vector<int> &seq, int well);

    /** Get the energy of a selected site within a sequence.
     * @param seq, Protein sequence.
     * @param m, Index of site. */
    double get_single_site_energy(std::vector<int> &seq, int m);

    /** Compute partition function. */
    double get_Z();

    /** Compute energies of each sequence in an MSA.
     *
    void calc_energy_dists(const char *name);
    */

    /** Convert to zero-sum gauge. */
    void convert_to_zero_sum();

    void dump_params(int iter);
};

#endif
