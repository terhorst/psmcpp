#include <map>

#include "inference_manager.h"
#include "transition.h"
#include "bin_key.h"
#include "marginalize_key.h"
#include "tensorslice.h"
#include "jcsfs.h"

PiecewiseConstantRateFunction<adouble>* defaultEta(const std::vector<double> &hidden_states)
{
    std::vector<std::vector<adouble> > params;
    std::vector<adouble> p{adouble(1.0)};
    params.push_back(p);
    params.push_back(p);
    return new PiecewiseConstantRateFunction<adouble>(params, hidden_states);
}

InferenceManager::InferenceManager(
        const int npop,
        const int sfs_dim,
        const std::vector<int> obs_lengths,
        const std::vector<int*> observations,
        const std::vector<double> hidden_states,
        double* const rho_vals,
        const std::vector<int> stitchpoints, // (0, stitch_1, stitch_2 ...stitch_n = L) --- stitchpoints are indexed by base position 
        ConditionedSFS<adouble> *csfs) :
    Stitchable(rho_vals, stitchpoints),
    saveGamma(false), folded(false),
    hidden_states(hidden_states),
    npop(npop),
    sfs_dim(sfs_dim),
    M(hidden_states.size() - 1),
    obs(map_obs(observations, obs_lengths)),
    stitchpoints(stitchpoints),
    rho_vals(rho_vals),
    stitch_to_block(process_stitchpoints()),
    csfs(csfs),
    hmms(obs.size()),
    pi(M),
    targets(fill_targets()),
    tb(targets, &emission_probs),
    ib{&pi, &tb, &emission_probs, &saveGamma},
    dirty({true, true, true}),
    eta(defaultEta(hidden_states))
{
    if (*std::min_element(hidden_states.begin(), hidden_states.end()) != 0.)
        throw std::runtime_error("first hidden interval should be [0, <something>)");
    // pi = Vector<adouble>::Zero(M);
    recompute_initial_distribution();
    transition = Matrix<adouble>::Zero(M, M);
    transition.setZero();
    InferenceBundle *ibp = &ib;
    std::map<double, Matrix<adouble>> *tm_ptr = &transition_map;
#pragma omp parallel for
    for (unsigned int i = 0; i < obs.size(); ++i)
    {
        DEBUG << "creating HMM";
        hmms[i].reset(new HMM(i, this->obs[i], ibp, rho_vals, stitch_to_block, tm_ptr));
    }

    // Collect all the block keys for recomputation later
    populate_emission_probs();
}

void InferenceManager::recompute_initial_distribution()
{
    for (int m = 0; m < M - 1; ++m)
    {
        pi(m) = exp(-(eta->R(hidden_states[m]))) - exp(-(eta->R(hidden_states[m + 1])));
        assert(pi(m) >= 0.0);
        assert(pi(m) <= 1.0);
    }
    pi(M - 1) = exp(-(eta->R(hidden_states[M - 1])));
    adouble small = eta->zero() + 1e-20;
    pi = pi.unaryExpr([small] (const adouble &x) { if (x < 1e-20) return small; return x; });
    CHECK_NAN(pi);
}

void InferenceManager::setRho(const double rho)
{
    this->rho = rho;
    dirty.rho = true;
}

void InferenceManager::setTheta(const double theta)
{
    this->theta = theta;
    dirty.theta = true;
}

void InferenceManager::parallel_do(std::function<void(hmmptr&)> lambda)
{
#pragma omp parallel for
    for (auto it = hmms.begin(); it < hmms.end(); ++it)
        lambda(*it);
}

template <typename T>
std::vector<T> InferenceManager::parallel_select(std::function<T(hmmptr &)> lambda)
{
    std::vector<T> ret(hmms.size());
#pragma omp parallel for
    for (unsigned int i = 0; i < hmms.size(); ++i)
        ret[i] = lambda(hmms[i]);
    return ret;
}
template std::vector<double> InferenceManager::parallel_select(std::function<double(hmmptr &)>);
template std::vector<adouble> InferenceManager::parallel_select(std::function<adouble(hmmptr &)>);

void InferenceManager::Estep(bool fbonly)
{
    DEBUG << "E step";
    do_dirty_work();
    parallel_do([fbonly] (hmmptr &hmm) { hmm->Estep(fbonly); });
}

std::vector<adouble> InferenceManager::Q(void)
{
    DEBUG << "InferenceManager::Q";
    do_dirty_work();
    std::vector<Vector<adouble> > ps = parallel_select<Vector<adouble> >([] (hmmptr &hmm) { return hmm->Q(); });
    adouble q1 = 0, q2 = 0, q3 = 0;
    for (unsigned int i = 0; i < ps.size(); ++i)
    {
        q1 += ps[i][0];
        q2 += ps[i][1];
        q3 += ps[i][2];
    }
    DEBUG << "\nq1:" << q1.value() << " [" << q1.derivatives().transpose() << "]\nq2:"
        << q2.value() << " [" << q2.derivatives().transpose() << "]\nq3:" << q3.value()
        << " [" << q3.derivatives().transpose() << "]\n";
    return {q1, q2, q3};
}

std::vector<std::map<block_key, Vector<double> >*> InferenceManager::getGammaSums()
{
    std::vector<std::map<block_key, Vector<double> >*> ret;
    for (auto &hmm : hmms)
        ret.push_back(&hmm->gamma_sums);
    return ret;
}

std::vector<Matrix<double>*> InferenceManager::getGammas()
{
    std::vector<Matrix<double>*> ret;
    for (auto &hmm : hmms)
        ret.push_back(&hmm->gamma);
    return ret;
}

std::vector<Matrix<double>*> InferenceManager::getXisums()
{
    std::vector<Matrix<double>*> ret;
    for (auto &hmm : hmms)
        ret.push_back(&hmm->xisum);
    return ret;
}

Matrix<adouble>& InferenceManager::getPi(void)
{
    static Matrix<adouble> mat;
    mat = pi;
    return mat;
}

Matrix<adouble>& InferenceManager::getTransition(void)
{
    return transition;
}

Matrix<adouble>& InferenceManager::getEmission(void)
{
    return emission;
}

std::map<block_key, Vector<adouble> >& InferenceManager::getEmissionProbs()
{
    return emission_probs;
}

std::vector<double> InferenceManager::loglik(void)
{
    return parallel_select<double>([] (hmmptr &hmm) { return hmm->loglik(); });
}

// Begin stuff for NPop inference manager
std::vector<Eigen::Matrix<int, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> > 
    InferenceManager::map_obs(const std::vector<int*> &observations, const std::vector<int> &obs_lengths)
{
    std::vector<Eigen::Matrix<int, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> > ret;
    for (unsigned int i = 0; i < observations.size(); ++i)
        ret.push_back(Eigen::Matrix<int, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>::Map(
                    observations[i], obs_lengths[i], 1 + 3 * npop));
    return ret;
}

void InferenceManager::populate_emission_probs()
{
    Vector<adouble> tmp;
    for (auto ob : obs)
    {
        const int q = ob.cols() - 1;
        for (int i = 0; i < ob.rows(); ++i)
        {
            block_key key(ob.row(i).tail(q).transpose());
            if (emission_probs.count(key) == 0)
            {
                emission_probs.insert({key, tmp});
                bpm_keys.push_back(key);
            }
        }
    }
}

void InferenceManager::do_dirty_work()
{
    // Figure out what changed and recompute accordingly.
    if (dirty.eta)
    {
        recompute_initial_distribution();
        sfss = csfs->compute(*eta);
    }
    if (dirty.theta or dirty.eta)
        recompute_emission_probs();
    if (dirty.eta or dirty.rho)
        recompute_transitions();
        //transition = compute_transition(*eta, rho);
    if (dirty.theta or dirty.eta or dirty.rho)
        tb.update(transition_map);
        //tb.update(transition);
    // restore pristine status
    dirty = {false, false, false};
}

void InferenceManager::recompute_transitions()
{
    transition_map.clear();
    for (size_t i = 0; i < stitchpoints.size() - 1; ++i)
    {
        double rho = *(rho_vals + i);
        if (transition_map.count(rho) == 0)
            transition_map.emplace(rho, compute_transition(*eta, rho));
    }
}

// Finds block index of each stitchpoint
std::vector<int> InferenceManager::process_stitchpoints()
{
    //TOFIX: pretty wasteful and should really change obs to break up spans in which stitchpoints occur
    Eigen::Matrix<int, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> ob;
    ob = obs.at(0);
    int current = 0;
    int prev = 0;
    std::vector<int>::const_iterator it = stitchpoints.begin();
    std::vector<int> ret;
    for (int i = 0; i < ob.rows(); ++i)
    {
        current += ob(i, 0);
        if(current > *it && prev <= *it)
        {
            ret.push_back(i);
            ++it;
        }
        prev = current;
    }
    ret.push_back(ob.rows());
    return ret;
}

// Returns targets over span, block_key, rho addresses
std::set<std::tuple<int, block_key, double* const> > InferenceManager::fill_targets()
{
    std::set<std::tuple<int, block_key, double* const> > ret;
    for (auto ob : obs)
    {
        const int q = ob.cols() - 1;
        for (int i = 0; i < ob.rows(); ++i)
        {
            if (ob(i, 0) > 1)
            {
                ret.emplace(ob(i, 0), block_key(ob.row(i).tail(q).transpose()), const_cast<double* const>(map_to_rho(i)));
            }
        }
    }
    return ret;
}

void InferenceManager::setParams(const ParameterVector &params)
{
    eta.reset(new PiecewiseConstantRateFunction<adouble>(params, hidden_states));
    dirty.eta = true;
}

template <size_t P>
FixedVector<int, 2 * P> NPopInferenceManager<P>::make_tensordims()
{
    FixedVector<int, 2 * P> ret;
    for (int p = 0; p < P; ++p)
    {
        ret(2 * p) = na(p) + 1;
        ret(2 * p + 1) = n(p) + 1;
    }
    return ret;
}

template <size_t P>
block_key NPopInferenceManager<P>::bk_to_map_key(const block_key &bk)
{
    Vector<int> ret(2 * P);
    for (int p = 0; p < P; ++p)
    {
        ret(2 * p) = bk(3 * p);
        ret(2 * p + 1) = bk(3 * p + 1);
    }
    return block_key(ret);
}

template <size_t P>
std::map<block_key, std::map<block_key, double> > NPopInferenceManager<P>::construct_bins(const bool binning)
{
    std::map<block_key, std::map<block_key, double> > ret;
    for (auto ob : obs)
    {
        const int q = ob.cols() - 1;
        for (int i = 0; i < ob.rows(); ++i)
        {
            block_key bk(ob.row(i).tail(q).transpose());
            if (ret.count(bk) == 0)
            {
                std::map<block_key, double> m;
                for (const block_key &kbin : bin_key<P>::run(bk, na, binning))
                    for (const auto &p : marginalize_key<P>::run(kbin.vals, n, na))
                        m[bk_to_map_key(p.first)] += p.second;
                ret[bk] = m;
            }
        }
    }
    return ret;
}

template <size_t P>
void NPopInferenceManager<P>::recompute_emission_probs()
{
    // Initialize emission matrix
    // Due to lack of good support for tensors, we store the emission
    // tensor in "flattened" matrix form. Note that this is actually
    // 1 larger along each axis than the true number, because the SFS
    // ranges in {0, 1, ..., n_pop_k}.
    emission = Matrix<adouble>::Zero(M, (na(0) + 1) * sfs_dim);
    emission.setZero();

    Eigen::Matrix<adouble, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> em_tmp(na(0) + 1, sfs_dim);
    std::vector<Matrix<adouble> > new_sfss = incorporate_theta(sfss, theta);
    for (int m = 0; m < M; ++m)
    {
        CHECK_NAN(new_sfss[m]);
        em_tmp = new_sfss[m];
        emission.row(m) = Matrix<adouble>::Map(em_tmp.data(), 1, (na(0) + 1) * sfs_dim);
    }
    
    DEBUG << "recompute B";
    Matrix<adouble> e2 = Matrix<adouble>::Zero(M, 2);
    std::vector<adouble> avg_ct = eta->average_coal_times();
    adouble small = eta->zero() + 1e-20;
    for (int m = 0; m < M; ++m)
    {
        if (std::isnan(avg_ct[m].value()))
        {
            // if the two lineages are separated by a split, their
            // average coalescence time within each interval before
            // (more recently than) the split is undefined. in this
            // case, assign very low probabilities to all such
            // observations.
            e2(m, 0) = small;
            e2(m, 1) = small;
        }
        else
        {
            e2(m, 1) = 2. * theta * avg_ct[m];
            e2(m, 0) = 1. - e2(m, 1);
        }
        // CHECK_NAN(e2(m, 1));
    }
    const adouble zero = eta->zero();
    const adouble one = zero + 1.;
#pragma omp parallel for
    for (auto it = bpm_keys.begin(); it < bpm_keys.end(); ++it)
    {
        // std::set<block_key> keys;
        // keys.insert(key);
        /*
        if (this->folded)
        {
            Vector<int> new_key(it->size());
            for (size_t p = 0; p < P; ++p)
            {
                int a = key(3 * p);
                int b = key(3 * p + 1);
                int nb = key(3 * p + 2);
                new_key(1 + 2 * p) = nb - b;
                new_key(2 + 2 * p) = nb;
            }
            keys.emplace(new_key);
        }
        */
        const block_key k = *it;
        std::array<std::set<FixedVector<int, 3> >, P> s;
        Vector<adouble> tmp(M);
        tmp.fill(zero);
        bool reduced = true;
        FixedVector<int, P> a, b, nb;
        for (unsigned int p = 0; p < P; ++p)
        {
            a(p) = k(3 * p);
            b(p) = k(1 + 3 * p);
            nb(p) = k(2 + 3 * p);
            reduced &= nb(p) == 0;
        }
        if (reduced and (a.isConstant(-1) or (a.minCoeff() >= 0)))
        {
            if (a.isConstant(-1))
                tmp.fill(one);
            else // if (a.minCoeff() >= 0)
                tmp = e2.col(a.sum() % 2);
        }
        else
            for (const auto &p : bins.at(k))
                tmp += p.second * tensorRef(p.first);
        if (tmp.maxCoeff() > 1.0 or tmp.minCoeff() <= 0.0)
        {
            std::cout << k << std::endl;
            std::cout << tmp.template cast<double>().transpose() << std::endl;
            std::cout << tmp.maxCoeff() << std::endl;
            throw std::runtime_error("probability vector not in [0, 1]");
        }
        CHECK_NAN(tmp);
        this->emission_probs.at(k) = tmp;
    }
    DEBUG << "recompute done";
}

Vector<adouble> TwoPopInferenceManager::tensorRef(const block_key &key)
{
    Vector<int> vals = key.vals;
    // Support the case of (0, 2) by flipping the key
    if (na(0) == 0 and na(1) == 2)
    {
        vals.head(2) = key.vals.tail(2);
        vals.tail(2) = key.vals.head(2);
    }
    return tensorSlice<4>::run(emission, key.vals, tensordims);
}

template <size_t P>
Vector<adouble> NPopInferenceManager<P>::tensorRef(const block_key &key)
{
    return tensorSlice<2 * P>::run(emission, key.vals, tensordims);
}


Matrix<adouble> sfs_cython(const int n, const ParameterVector p, 
        const double t1, const double t2, bool below_only)
{
    std::vector<double> hs{t1, t2};
    OnePopConditionedSFS<adouble> csfs(n);
    std::vector<Matrix<adouble> > v;
    PiecewiseConstantRateFunction<adouble> eta(p, hs);
    if (below_only)
        v = csfs.compute_below(eta);
    else
        v = csfs.compute(eta);
    return v[0];
}

OnePopInferenceManager::OnePopInferenceManager(
            const int n,
            const std::vector<int> obs_lengths,
            const std::vector<int*> observations,
            const std::vector<double> hidden_states,
            double* const rho_vals,
            const std::vector<int> stitchpoints, // (0, stitch_1, stitch_2 ...stitch_n = L) --- stitchpoints are indexed by base position 
            const bool binning) :
        NPopInferenceManager(
                FixedVector<int, 1>::Constant(n),
                FixedVector<int, 1>::Constant(2),
                obs_lengths, observations, hidden_states, 
                rho_vals, stitchpoints,
                new OnePopConditionedSFS<adouble>(n),
                binning) {}

JointCSFS<adouble>* create_jcsfs(int n1, int n2, int a1, int a2, const std::vector<double> &hidden_states)
{
    if (a1 == 0 and a2 == 2)
    {
        std::swap(n1, n2);
        std::swap(a1, a2);
    }
    return new JointCSFS<adouble>(n1, n2, a1, a2, hidden_states);

}

TwoPopInferenceManager::TwoPopInferenceManager(
            const int n1, const int n2,
            const int a1, const int a2,
            const std::vector<int> obs_lengths,
            const std::vector<int*> observations,
            const std::vector<double> hidden_states,
            double* const rho_vals,
            const std::vector<int> stitchpoints, // (0, stitch_1, stitch_2 ...stitch_n = L) --- stitchpoints are indexed by base position 
            const bool binning) :
        NPopInferenceManager(
                (FixedVector<int, 2>() << n1, n2).finished(),
                (FixedVector<int, 2>() << a1, a2).finished(),
                obs_lengths, observations, hidden_states, 
                rho_vals, stitchpoints,
                create_jcsfs(n1, n2, a1, a2, hidden_states),
                binning), a1(a1), a2(a2)
{
    if (a1 + a2 != 2)
        throw std::runtime_error("configuration not supported");
}

void TwoPopInferenceManager::setParams(
        const ParameterVector &distinguished_params,
        const ParameterVector &params1,
        const ParameterVector &params2,
        const double split)
{
    InferenceManager::setParams(distinguished_params);
    dynamic_cast<JointCSFS<adouble>*>(csfs.get())->pre_compute(params1, params2, split);
}

template class NPopInferenceManager<1>;
template class NPopInferenceManager<2>;
