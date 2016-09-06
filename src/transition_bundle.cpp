#include "transition_bundle.h"

// void TransitionBundle::update(const Matrix<adouble> &new_T)
// {
//     T = new_T;
//     Td = T.template cast<double>();
//     const int M = T.rows();
//     eigensystems.clear();
//     span_Qs.clear();
//     Matrix<double> tmp;
//     for (auto it = targets.begin(); it < targets.end(); ++it)
//     {
//         block_key key = it->second;
//         if (this->eigensystems.count(key) == 0)
//         {
//             tmp = this->emission_probs->at(key).template cast<double>().asDiagonal() * this->Td.transpose();
//             Eigen::EigenSolver<Matrix<double> > es(tmp);
//             this->eigensystems.emplace(key, es);
//         }
//     }
// #pragma omp parallel for
//     for (auto it = targets.begin(); it < targets.end(); ++it)
//     {
//         int span = it->first;
//         block_key key = it->second;
//         eigensystem eig = this->eigensystems.at(key);
//         Matrix<std::complex<double> > Q(M, M);
//         for (int a = 0; a < M; ++a)
//             for (int b = 0; b < M; ++b)
//                 if (a == b)
//                     Q(a, b) = std::pow(eig.d_scaled(a), span - 1) * (double)span;
//                 else
//                     Q(a, b) = (std::pow(eig.d_scaled(a), span) - std::pow(eig.d_scaled(b), span)) / 
//                         (eig.d_scaled(a) - eig.d_scaled(b));
// #pragma omp critical(emplace_Q)
//         {
//             this->span_Qs.emplace(*it, Q);
//         }
//     }
// }

void TransitionBundle::update(const std::map<double, Matrix<adouble>> &new_map)
{
    const int M = new_map.begin()->second.rows();
    transitions = new_map;
    eigensystems.clear();
    span_Qs.clear();
    Matrix<double> tmp;
    for (auto it = targets.begin(); it < targets.end(); ++it)
    {
        block_key bk = std::get<1>(*it);
        double rho = *std::get<2>(*it);
        std::pair<block_key, double> key = std::make_pair(bk, rho);
        if (this->eigensystems.count(key) == 0)
        {
            Matrix<double> Td = transitions.at(rho).template cast<double>();
            tmp = this->emission_probs->at(bk).template cast<double>().asDiagonal() * Td.transpose();
            Eigen::EigenSolver<Matrix<double> > es(tmp);
            this->eigensystems.emplace(key, es);
        }
    }
#pragma omp parallel for
    for (auto it = targets.begin(); it < targets.end(); ++it)
    {
        int span = std::get<0>(*it);
        block_key bk = std::get<1>(*it);
        double rho = *std::get<2>(*it);
        eigensystem eig = this->eigensystems.at(std::make_pair(bk, rho));
        Matrix<std::complex<double> > Q(M, M);
        for (int a = 0; a < M; ++a)
            for (int b = 0; b < M; ++b)
                if (a == b)
                    Q(a, b) = std::pow(eig.d_scaled(a), span - 1) * (double)span;
                else
                    Q(a, b) = (std::pow(eig.d_scaled(a), span) - std::pow(eig.d_scaled(b), span)) / 
                        (eig.d_scaled(a) - eig.d_scaled(b));
#pragma omp critical(emplace_Q)
        {
            this->span_Qs.emplace(*it, Q);
        }
    }
}
