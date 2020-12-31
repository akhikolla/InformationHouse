#include <RcppArmadillo.h>

const double pi = 3.141592653589793238463;

using namespace Rcpp;
using namespace arma;

// [[Rcpp::depends(RcppArmadillo)]]

//' Implement Kalman smoothing
//'
//' Estimate the hidden state and expected log-likelihood given the observations, exogeneous input and system parameters. This is an internal function and should not be called directly.
//'
//' @param y Observation matrix (may need to be normalized and centered before hand) (q rows, T columns)
//' @param u Input matrix for the state equation (m_u rows, T columns)
//' @param v Input matrix for the output equation (m_v rows, T columns)
//' @param theta A list of system parameters (A, B, C, D, Q, R)'
//' @param stdlik Boolean, whether the likelihood is divided by the number of observations. Standardizing the likelihood this way may speed up convergence in the case of long time series.
//' @return A list of fitted elements (X, Y, V, J, and lik)
//' @section Note: This code only works on one dimensional state and output at the moment. Therefore, transposing is skipped, and matrix inversion is treated as /, and log(det(Sigma)) is treated as log(Sigma).
// [[Rcpp::export]]
List Kalman_smoother(arma::mat y, arma::mat u, arma::mat v, List theta, bool stdlik = true) {

    // Model parameters
    mat A = theta["A"];
    mat B = theta["B"];
    mat C = theta["C"];
    mat D = theta["D"];
    mat Q = theta["Q"];
    mat R = theta["R"];
    mat x1 = theta["mu1"];
    mat V1 = theta["V1"];

    // Declare variables
    int T = y.n_cols; // Time series length
    mat Xp, Vp, Yp, Xu, Vu, Xs, Vs, Ys, Cov, I, J, K, Sigma;

    I.eye(1,1);

    // FORWARD PASS ==================

    // Prior t|t-1, p stands for prior
    Xp.zeros(1, T);
    Vp.zeros(1, T);
    Yp.zeros(1, T);

    // Remember C++ starts from 0
    Xp.col(0) = x1;
    Vp.col(0) = V1;
    if (v.n_cols == 1) {
        Yp.col(0) = C * Xp.col(0);
    } else {
        Yp.col(0) = C * Xp.col(0) + D * v.col(0);
    }

    // Posterior t|t, u stands for updated
    Xu.zeros(1, T);
    Vu.zeros(1, T);

    // First time step
    if (NumericMatrix::is_na(y(0,0))) {
        Xu.col(0) = Xp.col(0);
        Vu.col(0) = Vp.col(0);
    } else {
        K = Vp.col(0) * C * inv(C * Vp.col(0) * C + R);
        Xu.col(0) = Xp.col(0) + K * (y.col(0) - Yp.col(0));
        Vu.col(0) = (I - K * C) * Vp.col(0);
    }
    // Subsequent time steps
    for (int t=1; t<T; t++) {
        if (u.n_cols == 1) {
            Xp.col(t) = A * Xu.col(t-1);
        } else {
            Xp.col(t) = A * Xu.col(t-1) + B * u.col(t-1);
        }
        Vp.col(t) = A * Vu.col(t-1) * A + Q;
        if (v.n_cols == 1) {
            Yp.col(t) = C * Xp.col(t);
        } else {
            Yp.col(t) = C * Xp.col(t) + D * v.col(t);
        }
        if (NumericMatrix::is_na(y(0,t))) {
            Xu.col(t) = Xp.col(t);
            Vu.col(t) = Vp.col(t);
        } else {
            K = Vp.col(t) * C * inv(C * Vp.col(t) * C + R);
            Xu.col(t) = Xp.col(t) + K * (y.col(t) - Yp.col(t));
            Vu.col(t) = (I - K * C) * Vp.col(t);
        }
    }

    // BACKWARD PASS =========================
    // Smoothing, s stands for smoothed
    Xs = Xu;
    Vs = Vu;
    // Cov = Vs; // Cov(X_t+1, X_t)
    J.zeros(1, T);
    J.col(T-1) = Vu.col(T-1) * A * inv(A * Vu.col(T-1) * A + Q);
    for (int t=T-2; t>=0; t--) {
        J.col(t) = Vu.col(t) * A * inv(Vp.col(t+1));
        Xs.col(t) = Xu.col(t) + J.col(t) * (Xs.col(t+1) - Xp.col(t+1));
        Vs.col(t) = Vu.col(t) + J.col(t) * (Vs.col(t+1) - Vp.col(t+1)) * trans(J.col(t));
        // Cov.col(t) = Vs.col(t+1)*J;
    }
    // Final prediction
    if (v.n_cols == 1) {
        Ys = C*Xs;
    } else {
        Ys = C*Xs + D*v;
    }

    // Likelihood using the incomplete form, see equation (A.9) of Cheng and Sabes (2006)
    uvec obs = find_finite(y);  // Find indices where there are observations
    int n_obs = obs.size();     // Number of observations
    mat delta = y.cols(obs) - Yp.cols(obs);  // Innovations

    Sigma = Vp.cols(obs);
    for (int i=0; i < n_obs; i++) {
        Sigma.col(i) = C * Sigma.col(i) * C + R;
    }

    double lik = -0.5 * n_obs * log(2*pi) - 0.5 * accu(delta / Sigma % delta + log(Sigma));

    if (stdlik) lik = lik / n_obs;

    return List::create(_["X"] = Xs,
                        _["Y"] = Ys,
                        _["V"] = Vs,
                        _["J"] = J,
                        _["lik"] = lik);
}

//' Maximizing expected likelihood using analytical solution
//'
//' @inheritParams Kalman_smoother
//' @param fit result of [Kalman_smoother]
//' @return An object of class `theta`: a list of
// [[Rcpp::export]]
List Mstep(arma::mat y, arma::mat u, arma::mat v, List fit) {

    mat X = fit["X"];
    mat V = fit["V"];
    mat J = fit["J"];
    int T = y.n_cols;
    int p = u.n_rows;
    int q = v.n_rows;
    uvec obs = find_finite(y); // indices of observations

    // C, D, and R ================================
    // Sum for observed part
    mat Syx = y.cols(obs) * trans(X.cols(obs));
    mat Sxx = X.cols(obs) * trans(X.cols(obs)) + accu(V.cols(obs));
    mat C(1, 1);
    mat D; D.zeros(1, q);
    mat y_hat(1, T);

    if (v.n_cols > 1) {
        mat Syv = y.cols(obs) * trans(v.cols(obs));
        mat Sxv = X.cols(obs) * trans(v.cols(obs));
        mat Svx = trans(Sxv);
        mat Svv = v.cols(obs) * trans(v.cols(obs));

        // Solve linear system
        mat P1 = join_horiz(Syx, Syv);
        mat P2 = join_vert(join_horiz(Sxx, Sxv), join_horiz(Svx, Svv));
        mat CD  = P1 * inv(P2);

        C = CD.col(0);
        D = CD.cols(1, q);
        y_hat = C * X.cols(obs) + D * v.cols(obs);
    } else {
        C = Syx * inv(Sxx); // D = 0 as we didn't touch it
        y_hat = C * X.cols(obs);
    }

    // mat R = (delta_y * trans(delta_y) + C*accu(V.cols(obs))*C) / obs.size();
    mat R = ((y.cols(obs) - y_hat) * trans(y.cols(obs))) / obs.size();

    // A, B, and Q ================================
    mat Tx1x = X.cols(1, T-1) * trans(X.cols(0, T-2)) + V.cols(1, T-1) * trans(J.cols(0, T-2));
    mat Txx  = X.cols(0, T-2) * trans(X.cols(0, T-2)) + accu(V.cols(0, T-2));
    mat Txx1 = Tx1x.t();
    mat Tx1x1 = X.cols(1, T-1) * trans(X.cols(1, T-1)) + accu(V.cols(1, T-1));

    mat A(1, 1);
    mat B; B.zeros(1, p);
    mat Q(1, 1);

    if (u.n_cols > 1) {
        mat Tx1u = X.cols(1, T-1) * trans(u.cols(0, T-2));
        mat Tux  = u.cols(0, T-2) * trans(X.cols(0, T-2));
        mat Txu  = Tux.t();
        mat Tuu  = u.cols(0, T-2) * trans(u.cols(0, T-2));
        mat Tux1 = trans(Tx1u);

        mat P3 = join_horiz(Tx1x, Tx1u);
        mat P4 = join_vert(join_horiz(Txx, Txu), join_horiz(Tux, Tuu));
        mat AB  = P3 * inv(P4);

        A = AB.col(0);
        B = AB.cols(1, p);

        // Note on updating Q ===================================
        // Txx1 = Tx1x is not true for the first few time steps, but using P(t+1,t) is better
        // because this is the quantity derived from the backward recursion.
        // P(t, t+1) is not a standard quantity that is derived and therefore we shouldn't use
        // V_t J_{t+1} to calculate it.
        // Using the symmetric form matches the alternate formula for Q
        // mat Q = (Tx1x1 - A*Txx*A - Tx1u*trans(B) - B*Tux1 + B*Tuu*trans(B)) / (T-1);
        Q = (Tx1x1 - A*Txx1 - B*Tux1) / (T-1);
    } else {
        A = Tx1x * inv(Txx); // B = 0 as we didn't touch it
        Q = (Tx1x1 - A*Txx1) / (T-1);
    }

    // Initial state ====================================
    // Explicitly assign to matrices, otherwise "Not a matrix" error.
    mat mu1 = X.col(0);
    mat V1 = V.col(0);

    return List::create(_["A"] = A,
                        _["B"] = B,
                        _["C"] = C,
                        _["D"] = D,
                        _["Q"] = Q,
                        _["R"] = R,
                        _["mu1"] = mu1,
                        _["V1"] = V1);
}

//' Learn LDS model
//'
//' Estimate the hidden state and model parameters given observations and exogenous inputs using the EM algorithm. This is the key backend routine of this package.
//'
//' @inheritParams Kalman_smoother
//' @param theta0 A vector of initial values for the parameters
//' @param niter Maximum number of iterations, default 1000
//' @param tol Tolerance for likelihood convergence, default 1e-5. Note that the log-likelihood is normalized
//' @return A list of model results
//' * theta: model parameters (A, B, C, D, Q, R, mu1, V1) resulted from Mstep
//' * fit: results of Estep
//' * liks : vector of loglikelihood over the iteration steps
//' @section Note: This code only works on one dimensional state and output at the moment. Therefore, transposing is skipped, and matrix inversion is treated as /, and log(det(Sigma)) is treated as log(Sigma).
// [[Rcpp::export]]
List LDS_EM(arma::mat y, arma::mat u, arma::mat v, List theta0, int niter = 1000, double tol = 1e-5) {

    vec lik(niter);
    List theta = theta0;
    // The exit condition relies on two consecutive iterations so we need to manually do the first two.
    // i = 0
    List fit = Kalman_smoother(y, u, v, theta);
    lik[0] = fit["lik"];
    theta = Mstep(y, u, v, fit);
    // i = 1
    fit = Kalman_smoother(y, u, v, theta);
    lik[1] = fit["lik"];
    // subsequent iterations
    int lastIter = 2;
    for (int i=2; i<niter; i++) {
        // Check user interruption every 100 iterations; otherwise R can crash upon interruption.
        if (i % 100 == 0)
            Rcpp::checkUserInterrupt();
        // Iterations
        theta = Mstep(y, u, v, fit);
        fit = Kalman_smoother(y, u, v, theta);
        lik[i] = fit["lik"];
        lastIter++;
        // Check for convergence: terminates when the change in likelihood is less than tol
        //     for two consecutive iterations.
        // Use std::abs, otherwise the compiler may understand abs as
        //     int abs(int x) and returns 0, which stops the iterations immediately.
        if (std::abs(lik[i] - lik[i-1]) < tol && std::abs(lik[i-1] - lik[i-2]) < tol) {
            break;
        }
    }
    return List::create(_["theta"] = theta,
                        _["fit"] = fit,
                        _["liks"] = lik.head(lastIter),
                        _["lik"] = fit["lik"]); // return the final likelihood for convenient access
}

//' State propagation
//'
//' This function propagates the state trajectory based on the exogenous inputs only
//' (without measurement update), and calculates the corresponding log-likelihood
//'
//' @param theta A list of system parameters (A, B, C, D, Q, R)'
//' @param u Input matrix for the state equation (m_u rows, T columns)
//' @param v Input matrix for the output equation (m_v rows, T columns)
//' @param y Observations
//' @param stdlik Boolean, whether the likelihood is divided by the number of observations. Standardizing the likelihood this way may speed up convergence in the case of long time series.
//' @section Note: This code only works on one dimensional state and output at the moment. Therefore, transposing is skipped, and matrix inversion is treated as /, and log(det(Sigma)) is treated as log(Sigma).
//' @return A list of predictions and log-likelihood (X, Y, V, lik)
// [[Rcpp::export]]
List propagate(List theta, arma::mat u, arma::mat v, arma::mat y, bool stdlik = true) {

    // Model parameters
    mat A = theta["A"];
    mat B = theta["B"];
    mat C = theta["C"];
    mat D = theta["D"];
    mat Q = theta["Q"];
    mat R = theta["R"];
    mat x1 = theta["mu1"];
    mat V1 = theta["V1"];

    // Declare variables
    int T = u.n_cols; // Time series length
    mat Xp, Vp, Yp, Sigma, delta;

    // FORWARD PASS ==================

    // Prior t|t-1, p stands for prior
    Xp.zeros(1, T);
    Vp.zeros(1, T);

    // First time step, remember C++ starts from 0
    Xp.col(0) = x1;
    Vp.col(0) = V1;

    // Next time steps
    for (int t=1; t<T; t++) {
        if (u.n_cols == 1) {
            Xp.col(t) = A*Xp.col(t-1);
        } else {
            Xp.col(t) = A*Xp.col(t-1) + B*u.col(t-1);
        }
        Vp.col(t) = A*Vp.col(t-1)*A + Q;
    }

    // Prediction =====================
    if (v.n_cols == 1) {
        Yp = C*Xp;
    } else {
        Yp = C*Xp + D*v;
    }

    // Likelihood ======================
    uvec obs = find_finite(y);  // Find indices where there are observations
    int n_obs = obs.size();     // Number of observations
    delta = y.cols(obs) - Yp.cols(obs);  // Innovations

    Sigma = Vp.cols(obs);
    for (int i=0; i < n_obs; i++) {
        Sigma.col(i) = C * Sigma.col(i) * C + R;
    }

    double lik = -0.5*n_obs*log(2*pi) - 0.5*accu(delta / Sigma % delta + log(Sigma));

    if (stdlik) lik = lik / n_obs;

    return List::create(_["X"] = Xp,
                        _["Y"] = Yp,
                        _["V"] = Vp,
                        _["lik"] = lik);
}
