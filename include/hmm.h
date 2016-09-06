#ifndef HMM_H
#define HMM_H

#include <map>

class InferenceManager;
struct InferenceBundle;

class HMM
{
    friend class InferenceManager;

    public:
    HMM(const int hmm_num, const Matrix<int> &obs, const InferenceBundle *ib, 
        std::vector<double*> rho_vals, const std::vector<int> stitch_to_block, 
        const std::map<double, Matrix<adouble>> *transition_map);
    void Estep(bool);
    double loglik(void);
    Vector<adouble> Q(void);

    private:
    HMM(HMM const&) = delete;
    HMM& operator=(HMM const&) = delete;
    // Methods
    void domain_error(double);
    inline block_key ob_key(int i) { return block_key(obs.row(i).transpose().tail(obs.cols() - 1)); }
    double * map_to_rho(int i);

    // Instance variables
    const int hmm_num;
    const Matrix<int> obs;
    const InferenceBundle *ib;
    const int M, L;
    double ll;
    Matrix<double> alpha_hat, xisum, gamma;
    Vector<double> c, gamma0;
    std::map<block_key, Vector<double> > gamma_sums;

    // Extra stuff
    std::vector<double*> rho_vals;
    const std::vector<int> stitch_to_block;
    const std::map<double, Matrix<adouble>> *transition_map;
};

#endif
