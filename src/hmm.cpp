#include "common.h"
#include "inference_manager.h"
#include "inference_bundle.h"
#include "hmm.h"

HMM::HMM(const int hmm_num, const Matrix<int> &obs, const InferenceBundle* ib, 
        std::vector<double*> rho_vals, const std::vector<int> stitch_to_block, 
        const std::map<double, Matrix<adouble>> *transition_map) : 
    hmm_num(hmm_num), obs(obs), ib(ib), M(ib->pi->rows()), L(obs.rows()), ll(0.),
    alpha_hat(M, L), xisum(M, M), c(L), rho_vals(rho_vals), stitch_to_block(stitch_to_block), transition_map(transition_map)
{}

double HMM::loglik()
{
    return ll;
}

void HMM::domain_error(double ret)
{
    if (std::isinf(ret) or std::isnan(ret))
    {
        std::cout << ib->pi->template cast<double>() << std::endl << std::endl;
        //std::cout << ib->tb->Td << std::endl << std::endl;
        throw std::domain_error("badness encountered");
    }
}

double * HMM::map_to_rho(int i)
{
    for(size_t bk = 0; bk < stitch_to_block.size() - 1; bk++)
    {
        if(stitch_to_block[bk] <= i && stitch_to_block[bk+1] > i){
            return rho_vals[bk];
        }
    }
    throw std::runtime_error("Should map to one of the rho values");
}

void HMM::Estep(bool fbOnly)
{
    TransitionBundle *tb = ib->tb;
    alpha_hat = Matrix<double>::Zero(M, L);
    if (*(ib->saveGamma))
        gamma = Matrix<double>::Zero(M, L);
    gamma_sums.clear();
    const Vector<double> z = Vector<double>::Zero(M);
    gamma_sums.emplace(ob_key(0), z);
    alpha_hat.col(0) = ib->pi->template cast<double>().cwiseProduct(ib->emission_probs->at(ob_key(0)).template cast<double>());
    c(0) = alpha_hat.col(0).sum();
    alpha_hat.col(0) /= c(0);
    Eigen::DiagonalMatrix<double, Eigen::Dynamic, Eigen::Dynamic> B;
    DEBUG << "forward algorithm (HMM #" << hmm_num << ")";
    int prog = (int)((double)L * 0.1);
    ll = 0.;
    for (int ell = 1; ell < L; ++ell)
    {
        if (ell == prog)
        {
            DEBUG << "hmm " << hmm_num << ": " << (int)(100. * (double)ell / (double)L) << "%";
            prog += (int)((double)L * 0.1);
        }
        block_key bk = ob_key(ell);
        gamma_sums.emplace(bk, z);
        B = ib->emission_probs->at(bk).template cast<double>().asDiagonal();
        int span = obs(ell, 0);
        //grab rho here to get the correct eigensystem
        double * rho = map_to_rho(ell);
        std::pair<block_key, double> key = std::make_pair(bk, *rho);
        if (span > 1 and tb->eigensystems.count(key) > 0)
        {
            eigensystem es = tb->eigensystems.at(key);
            // alpha_hat.col(ell) = (es.P * (es.d.array().pow(span).matrix().asDiagonal() * 
            //             (es.Pinv * alpha_hat.col(ell - 1).template cast<std::complex<double> >()))).real();
            if (es.cplx)
                alpha_hat.col(ell) = (es.P * (es.d_scaled.array().pow(span).matrix().asDiagonal() * 
                            (es.Pinv * alpha_hat.col(ell - 1).template cast<std::complex<double> >()))).real();
            else
                alpha_hat.col(ell) = (es.P_r * (es.d_r_scaled.array().pow(span).matrix().asDiagonal() * 
                            (es.Pinv_r * alpha_hat.col(ell - 1))));

            // c(ell) = c_true(ell) * scale**(-span)
            ll += span * log(es.scale);
        }
        else
        {
            Matrix<double> T = (*transition_map).at(*rho).template cast<double>();
            if (span != 1)
                throw std::runtime_error("span != 1");
            Matrix<double> M = (B * T.transpose()).pow(span);
            alpha_hat.col(ell) = M * alpha_hat.col(ell - 1);
        }
        CHECK_NAN(alpha_hat.col(ell));
        c(ell) = alpha_hat.col(ell).sum();
        alpha_hat.col(ell) /= c(ell);
        ll += log(c(ell));
    }
    Vector<double> beta = Vector<double>::Ones(M), v(M), alpha(M);
    xisum.setZero();
    Matrix<std::complex<double> > Q(M, M);
    Matrix<double> Q_r(M, M);
    DEBUG << "backward algorithm (HMM #" << hmm_num << ")";
    for (int ell = L - 1; ell > 0; --ell)
    {
        v.setZero();
        int span = obs(ell, 0);
        block_key bk = ob_key(ell);
        //TODO: grab rho again and change the name of key
        B = ib->emission_probs->at(bk).template cast<double>().asDiagonal();

        double * rho = map_to_rho(ell);
        std::pair<block_key, double> key = std::make_pair(bk, *rho);
        if (span > 1 and tb->eigensystems.count(key) > 0)
        {
            eigensystem es = tb->eigensystems.at(key);
            if (es.cplx)
            {
                Q = es.Pinv * 
                    (alpha_hat.col(ell - 1) * beta.transpose()).template cast<std::complex<double> >() * 
                    es.P;
                Q = Q.cwiseProduct(tb->span_Qs.at(std::make_tuple(span, key.first, &key.second)));
                v = (es.P * es.d_scaled.asDiagonal() * Q * es.Pinv).diagonal().real() / c(ell);
                xisum += ((es.P * Q * es.Pinv).real() * B) / (c(ell) * es.scale);
                beta = (es.Pinv.transpose() * (es.d_scaled.array().pow(span).matrix().asDiagonal() * 
                            (es.P.transpose() * beta.template cast<std::complex<double> >()))).real();
            }
            else
            {
                Q_r = es.Pinv_r * 
                    (alpha_hat.col(ell - 1) * beta.transpose()) * 
                    es.P_r;
                Q_r = Q_r.cwiseProduct(tb->span_Qs.at(std::make_tuple(span, key.first, &key.second)).real());
                v = (es.P_r * es.d_r_scaled.asDiagonal() * Q_r * es.Pinv_r).diagonal() / c(ell);
                xisum += ((es.P_r * Q_r * es.Pinv_r) * B) / (c(ell) * es.scale);
                beta = (es.Pinv_r.transpose() * (es.d_r_scaled.array().pow(span).matrix().asDiagonal() * 
                            (es.P_r.transpose() * beta)));
            }
        }
        else
        {
            if (span != 1)
                throw std::runtime_error("span");
            v = alpha_hat.col(ell).cwiseProduct(beta);
            xisum += alpha_hat.col(ell - 1) * beta.transpose() * B / c(ell);
            Matrix<double> T = (*transition_map).at(*rho).template cast<double>();
            beta = T * (B * beta);
        }
        beta /= c(ell);
        CHECK_NAN(xisum);
        CHECK_NAN(v);
        CHECK_NAN(beta);
        gamma_sums.at(bk) += v;
        if (*(ib->saveGamma))
            gamma.col(ell) = v;
    }
    gamma0 = alpha_hat.col(0).cwiseProduct(beta);
    gamma_sums.at(ob_key(0)) += gamma0;
    if (*(ib->saveGamma))
        gamma.col(0) = gamma0;
    double * rho = map_to_rho(0);
    Matrix<double> T = (*transition_map).at(*rho).template cast<double>();
    xisum = xisum.cwiseProduct(T);
}

// void HMM::Estep(bool fbOnly)
// {
//     TransitionBundle *tb = ib->tb;
//     alpha_hat = Matrix<double>::Zero(M, L);
//     if (*(ib->saveGamma))
//         gamma = Matrix<double>::Zero(M, L);
//     Matrix<double> T = tb->Td;
//     gamma_sums.clear();
//     const Vector<double> z = Vector<double>::Zero(M);
//     gamma_sums.emplace(ob_key(0), z);
//     alpha_hat.col(0) = ib->pi->template cast<double>().cwiseProduct(ib->emission_probs->at(ob_key(0)).template cast<double>());
// 	c(0) = alpha_hat.col(0).sum();
//     alpha_hat.col(0) /= c(0);
//     Eigen::DiagonalMatrix<double, Eigen::Dynamic, Eigen::Dynamic> B;
//     DEBUG << "forward algorithm (HMM #" << hmm_num << ")";
//     int prog = (int)((double)L * 0.1);
//     ll = 0.;
//     for (int ell = 1; ell < L; ++ell)
//     {
//         if (ell == prog)
//         {
//             DEBUG << "hmm " << hmm_num << ": " << (int)(100. * (double)ell / (double)L) << "%";
//             prog += (int)((double)L * 0.1);
//         }
//         block_key key = ob_key(ell);
//         gamma_sums.emplace(key, z);
//         B = ib->emission_probs->at(key).template cast<double>().asDiagonal();
//         int span = obs(ell, 0);
//         if (span > 1 and tb->eigensystems.count(key) > 0)
//         {
//             eigensystem es = tb->eigensystems.at(key);
//             // alpha_hat.col(ell) = (es.P * (es.d.array().pow(span).matrix().asDiagonal() * 
//             //             (es.Pinv * alpha_hat.col(ell - 1).template cast<std::complex<double> >()))).real();
//             if (es.cplx)
//                 alpha_hat.col(ell) = (es.P * (es.d_scaled.array().pow(span).matrix().asDiagonal() * 
//                             (es.Pinv * alpha_hat.col(ell - 1).template cast<std::complex<double> >()))).real();
//             else
//                 alpha_hat.col(ell) = (es.P_r * (es.d_r_scaled.array().pow(span).matrix().asDiagonal() * 
//                             (es.Pinv_r * alpha_hat.col(ell - 1))));

//             // c(ell) = c_true(ell) * scale**(-span)
//             ll += span * log(es.scale);
//         }
//         else
//         {
//             if (span != 1)
//                 throw std::runtime_error("span != 1");
//             Matrix<double> M = (B * T.transpose()).pow(span);
//             alpha_hat.col(ell) = M * alpha_hat.col(ell - 1);
//         }
//         CHECK_NAN(alpha_hat.col(ell));
//         c(ell) = alpha_hat.col(ell).sum();
//         alpha_hat.col(ell) /= c(ell);
//         ll += log(c(ell));
//     }
//     Vector<double> beta = Vector<double>::Ones(M), v(M), alpha(M);
//     xisum.setZero();
//     Matrix<std::complex<double> > Q(M, M);
//     Matrix<double> Q_r(M, M);
//     DEBUG << "backward algorithm (HMM #" << hmm_num << ")";
//     for (int ell = L - 1; ell > 0; --ell)
//     {
//         v.setZero();
//         int span = obs(ell, 0);
//         block_key key = ob_key(ell);
//         B = ib->emission_probs->at(key).template cast<double>().asDiagonal();
//         if (span > 1 and tb->eigensystems.count(key) > 0)
//         {
//             eigensystem es = tb->eigensystems.at(key);
//             if (es.cplx)
//             {
//                 Q = es.Pinv * 
//                     (alpha_hat.col(ell - 1) * beta.transpose()).template cast<std::complex<double> >() * 
//                     es.P;
//                 Q = Q.cwiseProduct(tb->span_Qs.at({span, key}));
//                 v = (es.P * es.d_scaled.asDiagonal() * Q * es.Pinv).diagonal().real() / c(ell);
//                 xisum += ((es.P * Q * es.Pinv).real() * B) / (c(ell) * es.scale);
//                 beta = (es.Pinv.transpose() * (es.d_scaled.array().pow(span).matrix().asDiagonal() * 
//                             (es.P.transpose() * beta.template cast<std::complex<double> >()))).real();
//             }
//             else
//             {
//                 Q_r = es.Pinv_r * 
//                     (alpha_hat.col(ell - 1) * beta.transpose()) * 
//                     es.P_r;
//                 Q_r = Q_r.cwiseProduct(tb->span_Qs.at({span, key}).real());
//                 v = (es.P_r * es.d_r_scaled.asDiagonal() * Q_r * es.Pinv_r).diagonal() / c(ell);
//                 xisum += ((es.P_r * Q_r * es.Pinv_r) * B) / (c(ell) * es.scale);
//                 beta = (es.Pinv_r.transpose() * (es.d_r_scaled.array().pow(span).matrix().asDiagonal() * 
//                             (es.P_r.transpose() * beta)));
//             }
//         }
//         else
//         {
//             if (span != 1)
//                 throw std::runtime_error("span");
//             v = alpha_hat.col(ell).cwiseProduct(beta);
//             xisum += alpha_hat.col(ell - 1) * beta.transpose() * B / c(ell);
//             beta = T * (B * beta);
//         }
//         beta /= c(ell);
//         CHECK_NAN(xisum);
//         CHECK_NAN(v);
//         CHECK_NAN(beta);
//         gamma_sums.at(key) += v;
//         if (*(ib->saveGamma))
//             gamma.col(ell) = v;
//     }
//     gamma0 = alpha_hat.col(0).cwiseProduct(beta);
//     gamma_sums.at(ob_key(0)) += gamma0;
//     if (*(ib->saveGamma))
//         gamma.col(0) = gamma0;
//     xisum = xisum.cwiseProduct(T);
// }

Vector<adouble> HMM::Q(void)
{
    DEBUG << "HMM::Q";
    Vector<adouble> ret(3);
    DEBUG << "Shouldn't be calling Q";
    throw std::runtime_error("Q isn't properly implemented for stitchpoints");
    // ret.setZero();
    // Vector<adouble> pi = *(ib->pi);
    // ret(0) = (pi.array().log() * gamma0.array().template cast<adouble_base_type>()).sum();

    // std::vector<adouble> gs;
    // for (auto &p : gamma_sums)
    // {
    //     Vector<adouble> ep = ib->emission_probs->at(p.first);
    //     Vector<adouble> c = ep.array().log().matrix().cwiseProduct(p.second);
    //     if (ep.minCoeff() <= 0.0)
    //     {
    //         WARNING << "zeros detected in emission probability, key=" << p.first;
    //         gs.clear();
    //         gs[0] = adouble(-INFINITY);
    //         break;
    //     }
    //     gs.insert(std::end(gs), c.data(), c.data() + M);
    // }
    // ret(1) = doubly_compensated_summation(gs);

    // Matrix<adouble> T = ib->tb->T;
    // Matrix<adouble> log_T = T.array().log().matrix();
    // CHECK_NAN(log_T);
    // Matrix<adouble> prod = log_T.cwiseProduct(xisum.template cast<adouble_base_type>());
    // std::vector<adouble> es(prod.data(), prod.data() + M * M);
    // ret(2) = doubly_compensated_summation(es);

    // CHECK_NAN(ret);
    return ret;
}
