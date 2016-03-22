#ifndef TRANSITION_BUNDLE_H
#define TRANSITION_BUNDLE_H

struct eigensystem
{
    Eigen::MatrixXcd P, Pinv;
    Eigen::VectorXcd d;
    Matrix<double> P_r, Pinv_r;
    Vector<double> d_r;
    eigensystem(const Eigen::EigenSolver<Matrix<double> > &es) :
        P(es.eigenvectors()), Pinv(P.inverse()), d(es.eigenvalues()),
        P_r(P.real()), Pinv_r(Pinv.real()), d_r(d.real())
    {
        if (Pinv.imag().cwiseAbs().maxCoeff() > 1e-10 or P.imag().cwiseAbs().maxCoeff() > 1e-10)
        {
            std::cout << "imag parts: Pinv=" << Pinv.imag().cwiseAbs().maxCoeff() << " P=" << P.imag().cwiseAbs().maxCoeff() << std::endl;
            // throw std::runtime_error("Non-negligible imaginary part of eigendecomposition");
        }
    }
};

class TransitionBundle
{
    public:
    TransitionBundle(const std::set<std::pair<int, block_key> > &targets_s,
            const std::map<block_key, Vector<adouble> >* emission_probs) : 
        targets(targets_s.begin(), targets_s.end()),
        emission_probs(emission_probs) {}

    void update(const Matrix<adouble> &new_T)
    {
        T = new_T;
        Td = T.template cast<double>();
        const int M = T.rows();
        eigensystems.clear();
        span_Qs.clear();
#pragma omp parallel for
        for (auto it = targets.begin(); it < targets.end(); ++it)
        {
            Matrix<double> tmp;
            int span = it->first;
            block_key key = it->second;
#pragma omp critical(checkEigensystem)
            {
                if (eigensystems.count(key) == 0)
                {
                    tmp = emission_probs->at(key).template cast<double>().asDiagonal() * Td.transpose();
                    Eigen::EigenSolver<Matrix<double> > es(tmp);
                    eigensystems.emplace(key, es);
                }
            }
            eigensystem eig = eigensystems.at(key);
            Eigen::MatrixXcd ctmp = Eigen::MatrixXcd::Zero(M, M);
            for (int a = 0; a < M; ++a)
                for (int b = 0; b < M; ++b)
                    ctmp(a, b) = (a == b) ? (double)span * std::pow(eig.d(a), span - 1) : 
                        (std::pow(eig.d(a), span) - std::pow(eig.d(b), span)) / (eig.d(a) - eig.d(b));
#pragma omp critical(spanQinsert)
            span_Qs.emplace(*it, ctmp);
        }
    }
    Matrix<adouble> T;
    Matrix<double> Td;
    Eigen::VectorXcd d;
    Eigen::MatrixXcd P, Pinv;
    std::map<std::pair<int, block_key>, Matrix<std::complex<double> > > span_Qs;
    std::map<block_key, eigensystem> eigensystems;

    private:
    const std::vector<std::pair<int, block_key> > targets;
    const std::map<block_key, Vector<adouble> >* emission_probs;
};

#endif

