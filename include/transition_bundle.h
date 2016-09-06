#ifndef TRANSITION_BUNDLE_H
#define TRANSITION_BUNDLE_H

#include "common.h"
#include "block_key.h"

struct eigensystem
{
    Eigen::MatrixXcd P, Pinv;
    Eigen::VectorXcd d;
    Matrix<double> P_r, Pinv_r;
    double scale;
    Eigen::VectorXcd d_scaled;
    Vector<double> d_r;
    Vector<double> d_r_scaled;
    bool cplx;
    eigensystem(const Eigen::EigenSolver<Matrix<double> > &es) :
        P(es.eigenvectors()), Pinv(P.inverse()), d(es.eigenvalues()),
        P_r(P.real()), Pinv_r(Pinv.real()), 
        scale(d.cwiseAbs().maxCoeff()), d_scaled(d / scale),
        d_r(d.real()), d_r_scaled(d_r / scale),
        cplx(d.imag().cwiseAbs().maxCoeff() > 0)
    {
        DEBUG << "max imag: " << d.imag().cwiseAbs().maxCoeff();
    }
};

class TransitionBundle
{
    public:
    TransitionBundle(
            const std::set<std::tuple<int, block_key, double* const> > &targets_s,
            const std::map<block_key, Vector<adouble> >* emission_probs) : 
        targets(targets_s.begin(), targets_s.end()),
        emission_probs(emission_probs) {}

    void update(const std::map<double, Matrix<adouble>> &new_map);
    //void update(const Matrix<adouble> &new_T);
    //Matrix<adouble> T;
    //Matrix<double> Td;
    Eigen::VectorXcd d;
    Eigen::MatrixXcd P, Pinv;
    std::map<std::tuple<int, block_key, double*>, Matrix<std::complex<double> > > span_Qs;
    std::map<std::pair<block_key, double>, eigensystem> eigensystems;
    std::map<double, Matrix<adouble>> transitions;

    private:
    const std::vector<std::tuple<int, block_key, double* const> > targets;
    const std::map<block_key, Vector<adouble> >* emission_probs;
};

#endif

